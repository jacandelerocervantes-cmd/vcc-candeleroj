#!/usr/bin/env python
"""PASO 5a PARTE 2, etapa de DATOS (pesada). Jiang (6 lineas) y VCC 2025 (test).

CONSTRUCCION, DECLARADA ANTES DE MEDIR
  JIANG
  * Poblacion real: `comun.arma_linea` (= `suelo_arnes.build_fold_adatas` sin
    modificar), `seed_real = 20260001 + fold`, verificada por sha256 contra la
    huella guardada en `venn_ed.json`. Control = TODO el pool NT de la linea.
  * R(g,L) CRUDO: `compute_de` + `prep_de_side` del preset vcc2026 (parametros
    leidos de EvalConfig), `p_adj < 0.05`, y se quita la fila del propio gen
    (feature == target). Se COMPRUEBA gen a gen contra `venn_ed.json['ed']`
    (que salio del CLI `cell-eval2 run --write-degenes`).
  * Estadisticas de control por gen (34,025): CPM medio, tasa de deteccion,
    Fano crudo; y por linea: n de control, UMI mediana por celula de control.
  * R(g,L) IGUALADO EN POTENCIA: (i) control submuestreado sin reemplazo a
    n_min = min_L n_control(L) (semilla 20281201+fold); (ii) TODAS las celulas
    de la linea (control y perturbadas) adelgazadas binomialmente con
    p_L = D_min / D_L, D = UMI mediana por celula del control de la linea
    (semilla 20281211+fold). Se re-corre la misma puerta. Las estadisticas de
    control de la version igualada salen del control igualado.
  VCC 2025
  * Poblacion real: `v25.selecciona` (semillas 20280201 / 20280202, REUSADAS a
    proposito, la misma poblacion de `vcc2025_suelo`), leida en memoria (no se
    escribe a disco) y verificada por sha256 contra `real_2025_test_meta.json`.
  * R(g) con la misma puerta. Filtrables = CPM medio del control (18,400) > 5.
"""
import json
import os
import sys
import time

import numpy as np
import polars as pl
import scipy.sparse as sp

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/ntc_nulo")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C  # noqa: E402
from puerta import CFG, params  # noqa: E402
import semillas as S  # noqa: E402

OUT = "/datos/modelo/out/seleccion/p2_datos.json"
VENN = "/home/ubuntu/vcc/modelo/out/venn_ed.json"
META25 = "/datos/modelo/out/noche/b4/real_2025_test_meta.json"
UMBRAL = float(CFG.filter.filter_gene_min_cpm_cell)
T = float(CFG.de.p_adj_threshold)
CTRL = CFG.control  # "non-targeting"


def log(m):
    sys.stderr.write(f"[{time.strftime('%H:%M:%S')}] {m}\n"); sys.stderr.flush()


def stats(X):
    Xw = X.astype(np.float64).tocsr()
    n = Xw.shape[0]
    L = np.asarray(Xw.sum(axis=1)).ravel()
    Lc = L.copy(); Lc[Lc <= 0] = 1.0
    cpm = np.asarray(Xw.T @ ((1e6 / Lc) / n)).ravel()
    det = np.diff(Xw.tocsc().indptr) / n
    mu = np.asarray(Xw.mean(axis=0)).ravel()
    X2 = Xw.copy(); X2.data **= 2
    var = (np.asarray(X2.mean(axis=0)).ravel() - mu ** 2) * n / (n - 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        fano = np.where(mu > 0, var / mu, np.nan)
    return {"cpm": cpm, "det": det, "fano": fano, "n": int(n),
            "L_mediana": float(np.median(L))}


def puerta_sets(X, labels, var_names, perts):
    """La puerta del preset sobre una poblacion con control `CTRL`. Devuelve
    ({pert: sorted genes}, n_probados_por_pert)."""
    import anndata as ad
    import pandas as pd
    from cell_eval2.de_compute import compute_de
    from cell_eval2.de import prep_de_side
    c = CFG
    a = ad.AnnData(X=sp.csr_matrix(X).astype(np.float32),
                   obs=pd.DataFrame({c.pert_col: list(labels)},
                                    index=[f"c{i}" for i in range(X.shape[0])]),
                   var=pd.DataFrame(index=np.asarray(var_names, dtype=str)))
    df = compute_de(a, backend=c.de.backend, groupby=c.pert_col, reference=CTRL,
                    mean_calc=c.de.mean_calc, epsilon=c.de.epsilon, input_type=c.input_type,
                    target_sum=c.target_sum, clip_value=c.de.clip_value,
                    filter_gene_min_cpm_cell=c.filter.filter_gene_min_cpm_cell,
                    fdr_scope=c.de.fdr_scope, threads=c.num_threads,
                    replicate_col=c.de.replicate_col, device="cpu")
    tab, _ = prep_de_side(df, name="real", sort_by=c.de.sort_by,
                          nan_lfc_policy=c.de.nan_lfc_policy, min_abs_log2fc=c.de.min_abs_log2fc)
    tcol = "target" if "target" in tab.columns else c.pert_col
    nprob = {str(k): int(v) for k, v in
             tab.group_by(tcol).len().iter_rows()}
    sig = (tab.filter(pl.col("p_adj") < T)
              .filter(pl.col("feature").cast(pl.Utf8) != pl.col(tcol).cast(pl.Utf8))
              .select(tcol, "feature").unique())
    out = {p: [] for p in perts}
    for t, f in sig.iter_rows():
        out.setdefault(str(t), []).append(str(f))
    return {k: sorted(v) for k, v in out.items()}, nprob, tab.columns


def adelgaza(X, p, rng):
    X = sp.csr_matrix(X).astype(np.int64)
    if p >= 1.0:
        return X
    Y = X.copy()
    Y.data = rng.binomial(Y.data, p)
    Y.eliminate_zeros()
    return Y


def jiang(res):
    from suelo_arnes import load_global_obs, load_var_names, IFNB_PATH
    venn = json.load(open(VENN))
    panel = C.carga_panel()
    assert panel == venn["panel"]
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = np.asarray(load_var_names(IFNB_PATH), dtype=str)
    n_genes = len(var_names)
    pobl = {}
    for i, L in enumerate(C.LINES):
        t = time.time()
        real_ad, _pred, ctrl_pool, meta = C.arma_linea(L, i, panel, list(var_names),
                                                       cell_type, gene, n_genes)
        # la huella de venn_ed se tomo DESPUES de renombrar NT -> non-targeting
        real_ad.obs["target_gene"] = real_ad.obs["target_gene"].replace({C.CONTROL_LABEL: CTRL})
        h = C.huella_adata(real_ad)
        lab = real_ad.obs["target_gene"].to_numpy()
        X = real_ad.X.tocsr()
        del real_ad, _pred
        st = stats(ctrl_pool)
        ed, nprob, cols = puerta_sets(X, lab, var_names, panel)
        ref = venn["ed"][L]
        dif = {p: [len(set(ed.get(p, [])) - set(ref.get(p, []))),
                   len(set(ref.get(p, [])) - set(ed.get(p, [])))] for p in panel}
        pobl[L] = (X, lab)
        res["jiang"]["crudo"][L] = {
            "huella_real": h, "huella_venn": venn["verif"][L]["huella_real"],
            "huella_coincide": h == venn["verif"][L]["huella_real"],
            "n_control": st["n"], "L_mediana": st["L_mediana"],
            "cpm": st["cpm"].tolist(), "det": st["det"].tolist(),
            "fano": [None if not np.isfinite(x) else float(x) for x in st["fano"]],
            "n_pasa": int((st["cpm"] > UMBRAL).sum()),
            "ed": ed, "n_probados": nprob, "cols_tabla": list(cols),
            "perts_identicas_venn": int(sum(1 for p in panel if dif[p] == [0, 0])),
            "dif_venn": dif, "dt_s": time.time() - t}
        log(f"J crudo {L}: huella={h == venn['verif'][L]['huella_real']} "
            f"idem_venn={res['jiang']['crudo'][L]['perts_identicas_venn']}/33 "
            f"nctrl={st['n']} Lmed={st['L_mediana']:.0f} dt={time.time()-t:.0f}s")
        escribe(res)
    n_min = min(res["jiang"]["crudo"][L]["n_control"] for L in C.LINES)
    D_min = min(res["jiang"]["crudo"][L]["L_mediana"] for L in C.LINES)
    res["jiang"]["igualacion"] = {"n_min": n_min, "D_min": D_min}
    for i, L in enumerate(C.LINES):
        t = time.time()
        X, lab = pobl.pop(L)
        D_L = res["jiang"]["crudo"][L]["L_mediana"]
        p = D_min / D_L
        rc = np.random.default_rng(S.SEMILLAS["p2_ig_control"][i])
        ra = np.random.default_rng(S.SEMILLAS["p2_ig_adelgaza"][i])
        ic = np.where(lab == CTRL)[0]
        ip = np.where(lab != CTRL)[0]
        keep_c = np.sort(rc.choice(ic, n_min, replace=False))
        idx = np.concatenate([ip, keep_c])
        Xe = adelgaza(X[idx], p, ra)
        labe = lab[idx]
        st = stats(Xe[labe == CTRL])
        ed, nprob, _ = puerta_sets(Xe, labe, var_names, panel)
        res["jiang"]["igualado"][L] = {
            "p_adelgaza": p, "n_control": st["n"], "L_mediana": st["L_mediana"],
            "cpm": st["cpm"].tolist(), "det": st["det"].tolist(),
            "n_pasa": int((st["cpm"] > UMBRAL).sum()), "ed": ed, "n_probados": nprob,
            "semilla_control": S.SEMILLAS["p2_ig_control"][i],
            "semilla_adelgaza": S.SEMILLAS["p2_ig_adelgaza"][i], "dt_s": time.time() - t}
        log(f"J igualado {L}: p={p:.3f} nctrl={st['n']} Lmed={st['L_mediana']:.0f} "
            f"|ED|={sum(len(v) for v in ed.values())} dt={time.time()-t:.0f}s")
        escribe(res)
        del X, Xe
    res["jiang"]["panel"] = panel
    res["jiang"]["genes"] = var_names.tolist()


def vcc2025(res):
    import v25
    from suelo_arnes import read_rows_csr
    import anndata as ad
    import pandas as pd
    t = time.time()
    old = json.load(open(META25))
    tg, _indptr, shape, var_names, *_ = v25.lee_obs(v25.SPLIT_SUP)
    s = v25.selecciona(tg)
    X = read_rows_csr(v25.ruta(v25.SPLIT_SUP), s["idx"], shape[1]).astype(np.int64)
    lab = np.asarray(tg[s["idx"]], dtype=object)
    a = ad.AnnData(X=X, obs=pd.DataFrame({"target_gene": list(lab)}),
                   var=pd.DataFrame(index=var_names))
    a.obs_names_make_unique()
    h = C.huella_adata(a)
    del a
    var_names = np.asarray(var_names, dtype=str)
    lab = lab.astype(str)
    st = stats(X[lab == CTRL])
    ed, nprob, cols = puerta_sets(X, lab, var_names, s["perts"])
    res["v2025"] = {"huella": h, "huella_meta": old["huella"], "huella_coincide": h == old["huella"],
                    "n_control": st["n"], "L_mediana": st["L_mediana"],
                    "cpm": st["cpm"].tolist(), "det": st["det"].tolist(),
                    "n_pasa": int((st["cpm"] > UMBRAL).sum()), "perts": s["perts"],
                    "por_pert": s["por_pert"], "ed": ed, "n_probados": nprob,
                    "genes": var_names.tolist(), "dt_s": time.time() - t}
    log(f"2025: huella={h == old['huella']} pasa={res['v2025']['n_pasa']} "
        f"|ED| suma={sum(len(v) for v in ed.values())} dt={time.time()-t:.0f}s")
    escribe(res)


def escribe(res):
    tmp = OUT + ".tmp"
    json.dump(res, open(tmp, "w"))
    os.replace(tmp, OUT)


def main():
    t0 = time.time()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    sem = S.comprueba()
    assert sem["sin_choques"], sem
    res = {"preset": params(), "umbral": UMBRAL, "p_adj": T, "semillas": sem,
           "jiang": {"crudo": {}, "igualado": {}}, "t_inicio": time.strftime("%F %T")}
    que = sys.argv[1:] or ["jiang", "2025"]
    if os.path.exists(OUT):
        prev = json.load(open(OUT))
        if "jiang" not in que:
            res["jiang"] = prev.get("jiang", res["jiang"])
        if "2025" not in que and "v2025" in prev:
            res["v2025"] = prev["v2025"]
    if "jiang" in que:
        jiang(res)
    if "2025" in que:
        vcc2025(res)
    res["dt_total_s"] = time.time() - t0
    escribe(res)


if __name__ == "__main__":
    main()
