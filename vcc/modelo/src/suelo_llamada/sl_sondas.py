#!/usr/bin/env python
"""SONDAS DE DE por mitades, con la maquinaria del marcador (`cell-eval2 run
--write-degenes`). Uso: sl_sondas.py jiang | v2025

CONSTRUCCION, DECLARADA ANTES DE MEDIR
  * Poblacion: Jiang = `comun.arma_linea` (seed_real = 20260001 + indice de
    linea, la de venn_ed, verificada por sha256); 2025 = `v25.selecciona`
    (20280201/20280202, la de vcc2025_suelo, verificada por sha256).
  * Particion POR CELULA, una rng por (linea, j) con semilla 20282000+100*fold+j:
    primero el pool de control, despues cada perturbacion en orden de panel.
    Cada grupo se permuta y se parte en dos mitades disjuntas. 🔴 El CONTROL
    tambien se parte: mitad A lleva controles cA, mitad B lleva cB, disjuntos.
    Desempate de impares ALTERNADO: la celula sobrante va a A si j es par, a B
    si j es impar.
  * Una corrida por (linea, j): `-ar mitadA -ap mitadB --preset vcc2026
    --profile vcc2026 --set control_source=pred --write-degenes`. 🔴 Con el
    preset puro (`control_source=real`) el lado pred se llamaria contra el
    control del real (run.py, `_pred_reference`), compartiendo su ruido; con
    `control_source=pred` cada mitad usa su propio medio control. Se asserta en
    `run_params.yaml`.
  * VERIFICACION en j = 0 de cada linea/superficie: corrida con los papeles
    cambiados (`-ar mitadB -ap mitadA`); el conjunto de cada mitad por el lado
    REAL (preset sin tocar en ese lado) tiene que ser identico gen a gen, y con
    el mismo p_adj, al obtenido por el lado pred.
  * R = `prepare_de` + `p_adj < 0.05` + `_exclude_own_gene`, igual que venn_ed.
"""
import json
import os
import sys
import time

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C  # noqa: E402
import sl_semillas as S  # noqa: E402

sys.path.insert(0, C.SRC)
WORK = "/datos/modelo/out/suelo_llamada"
SETS = os.path.join(WORK, "sets")
CTRL = C.CONTROL_COMP


def log(m):
    sys.stderr.write(f"[{time.strftime('%H:%M:%S')}] {m}\n"); sys.stderr.flush()


def parte(n, rng, j):
    perm = rng.permutation(n)
    nA = n // 2 + (n % 2 if j % 2 == 0 else 0)
    return np.sort(perm[:nA]), np.sort(perm[nA:])


def mitades(X, lab, grupos, seed, j):
    """grupos: [CTRL] + perts en orden. Devuelve (idxA, idxB, cortes)."""
    rng = np.random.default_rng(seed)
    iA, iB, cortes = [], [], {}
    for g in grupos:
        idx = np.where(lab == g)[0]
        a, b = parte(idx.size, rng, j)
        iA.append(idx[a]); iB.append(idx[b])
        cortes[g] = [int(a.size), int(b.size)]
    iA, iB = np.concatenate(iA), np.concatenate(iB)
    assert np.intersect1d(iA, iB).size == 0
    return iA, iB, cortes


def escribe_h5(X, lab, var_names, path, dtype):
    import anndata as ad
    import pandas as pd
    a = ad.AnnData(X=sp.csr_matrix(X).astype(dtype),
                   obs=pd.DataFrame({"target_gene": list(lab)}),
                   var=pd.DataFrame(index=list(var_names)))
    a.obs_names_make_unique()
    a.write_h5ad(path)


def conjuntos(rd, perts):
    import polars as pl
    from cell_eval2.de import prepare_de
    from cell_eval2.metrics.de import _exclude_own_gene
    import techo_arnes as TA
    pr = prepare_de(os.path.join(rd, "de_pred.parquet"), os.path.join(rd, "de_real.parquet"),
                    control=CTRL, sort_by=TA.DE_SORT_BY, p_adj_threshold=TA.DE_P_ADJ_THRESHOLD,
                    nan_lfc_policy=TA.DE_NAN_LFC_POLICY, min_abs_log2fc=TA.DE_MIN_ABS_LOG2FC)
    T = pr.p_adj_threshold
    out = {}
    for side, df in (("real", pr.real_df), ("pred", pr.pred_df)):
        sig = _exclude_own_gene(df.filter(pl.col("p_adj") < T).select("target", "feature", "p_adj").unique(),
                                pr, who="suelo_llamada", side=side)
        d = {p: {} for p in perts}
        for t, f, q in sig.iter_rows():
            d.setdefault(str(t), {})[str(f)] = float(q)
        out[side] = d
    return out, float(T)


def corre_cli(pa, pb, rd, logp):
    cmd = [C.CELLEVAL2_PY, "run", "-ar", pa, "-ap", pb, "--preset", "vcc2026", "--profile", "vcc2026",
           "--set", "control_source=pred", "--pert-col", "target_gene", "--write-degenes", "-o", rd]
    p, dt = C.corre(cmd, logp)
    if p.returncode != 0:
        raise RuntimeError(f"cli fallo {logp}: {p.stderr[-1500:]}")
    rp = open(os.path.join(rd, "run_params.yaml")).read()
    assert "control_source: pred" in rp, "control_source=pred no aplicado"
    return dt


def sonda(tag, j, X, lab, var_names, perts, dtype, verifica):
    out_json = os.path.join(SETS, f"{tag}_j{j}.json")
    if os.path.exists(out_json):
        log(f"{tag} j={j}: ya hecho"); return
    C.guarda_disco()
    seed = S.semilla(tag, j)
    t = time.time()
    iA, iB, cortes = mitades(X, lab, [CTRL] + perts, seed, j)
    tmp = os.path.join(WORK, "tmp", f"{tag}_j{j}"); os.makedirs(tmp, exist_ok=True)
    pa, pb = os.path.join(tmp, "A.h5ad"), os.path.join(tmp, "B.h5ad")
    escribe_h5(X[iA], lab[iA], var_names, pa, dtype)
    escribe_h5(X[iB], lab[iB], var_names, pb, dtype)
    rd = os.path.join(tmp, "run")
    dt = corre_cli(pa, pb, rd, os.path.join(WORK, "logs", f"{tag}_j{j}.log"))
    sets, T = conjuntos(rd, perts)
    rec = {"tag": tag, "j": j, "semilla": seed, "cortes": cortes, "p_adj_threshold": T,
           "A": sets["real"], "B": sets["pred"], "dt_cli_s": dt}
    if verifica:
        rd2 = os.path.join(tmp, "run_swap")
        dt2 = corre_cli(pb, pa, rd2, os.path.join(WORK, "logs", f"{tag}_j{j}_swap.log"))
        s2, _ = conjuntos(rd2, perts)
        def igual(x, y):
            return all(x.get(p, {}) == y.get(p, {}) for p in perts)
        def ngen(x, y):
            return sum(len(set(x.get(p, {})) ^ set(y.get(p, {}))) for p in perts)
        rec["verif"] = {"B_real_igual_B_pred": igual(s2["real"], sets["pred"]),
                        "A_real_igual_A_pred": igual(sets["real"], s2["pred"]),
                        "genes_distintos_B": ngen(s2["real"], sets["pred"]),
                        "genes_distintos_A": ngen(sets["real"], s2["pred"]), "dt_cli_s": dt2}
        log(f"{tag} j={j} VERIF: {rec['verif']}")
    C.borra(tmp)
    rec["dt_total_s"] = time.time() - t
    C.escribe_json(out_json, rec)
    log(f"{tag} j={j}: |A|={sum(len(v) for v in rec['A'].values())} |B|={sum(len(v) for v in rec['B'].values())} "
        f"cli={dt:.0f}s total={rec['dt_total_s']:.0f}s")


def jiang(solo=None):
    """`solo` = una linea; si se da, solo se sondea esa (el runner por linea)."""
    from suelo_arnes import load_global_obs, load_var_names, IFNB_PATH
    venn = json.load(open(os.path.join(C.OUT, "venn_ed.json")))
    panel = C.carga_panel()
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    meta = {}
    for i, L in enumerate(S.LINES):
        if solo is not None and L != solo:
            continue
        if all(os.path.exists(os.path.join(SETS, f"{L}_j{j}.json")) for j in range(S.J)):
            continue
        real_ad, _p, _c, m = C.arma_linea(L, i, panel, var_names, cell_type, gene, len(var_names))
        real_ad.obs["target_gene"] = real_ad.obs["target_gene"].replace({C.CONTROL_LABEL: CTRL})
        h = C.huella_adata(real_ad)
        ok = h == venn["verif"][L]["huella_real"]
        meta[L] = {"huella": h, "coincide_venn": ok, "n_control": m["n_control_cells"]}
        log(f"{L}: huella coincide con venn_ed = {ok}")
        assert ok
        X = real_ad.X.tocsr(); lab = real_ad.obs["target_gene"].to_numpy().astype(str)
        dtype = X.dtype
        del real_ad, _p
        for j in range(S.J):
            sonda(L, j, X, lab, var_names, panel, dtype, verifica=(j == 0))
        C.escribe_json(os.path.join(SETS, f"{L}_meta.json"), meta[L])


def v2025():
    import v25
    import anndata as ad
    import pandas as pd
    from suelo_arnes import read_rows_csr
    old = json.load(open("/datos/modelo/out/noche/b4/real_2025_test_meta.json"))
    tg, _ip, shape, var_names, *_ = v25.lee_obs(v25.SPLIT_SUP)
    s = v25.selecciona(tg)
    X = read_rows_csr(v25.ruta(v25.SPLIT_SUP), s["idx"], shape[1]).astype(np.int64)
    lab = np.asarray(tg[s["idx"]], dtype=object)
    a = ad.AnnData(X=X, obs=pd.DataFrame({"target_gene": list(lab)}), var=pd.DataFrame(index=var_names))
    a.obs_names_make_unique()
    h = C.huella_adata(a); del a
    ok = h == old["huella"]
    log(f"2025: huella coincide = {ok}")
    assert ok
    C.escribe_json(os.path.join(SETS, "v2025_meta.json"), {"huella": h, "coincide": ok, "n_control": s["n_control_usadas"]})
    lab = lab.astype(str)
    X = X.astype(np.float32)          # el X de 2025 es float32 de origen
    for j in range(S.J):
        sonda("v2025", j, X, lab, var_names, s["perts"], np.float32, verifica=(j == 0))


if __name__ == "__main__":
    os.makedirs(SETS, exist_ok=True); os.makedirs(os.path.join(WORK, "logs"), exist_ok=True)
    sem = S.comprueba(); assert sem["sin_choques"], sem
    C.escribe_json(os.path.join(WORK, "semillas.json"), sem)
    modo = sys.argv[1]
    assert modo in ("jiang", "v2025"), modo
    if modo == "jiang":
        solo = sys.argv[2] if len(sys.argv) > 2 else None
        assert solo is None or solo in S.LINES, solo
        jiang(solo)
    else:
        v2025()
