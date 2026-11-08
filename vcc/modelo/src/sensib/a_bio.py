#!/usr/bin/env python
"""BRAZO A — LA BIOLOGIA: que estructura tiene el delta real.

🔴 Esto se mide en superficies que NO son la de entrega: Jiang (6 contextos,
somera, panel de 33) y VCC 2025 (1 contexto, profunda, 100 perturbaciones).
Ninguna cifra de aqui es biologia de D/E/F, que no existe en disco.

LA CONSTRUCCION, DECLARADA ANTES DE MEDIR
  A1/A3  La puerta es la del preset vcc2026, con los parametros leidos de
     `EvalConfig` (`seleccion/puerta.py`), igual que el PASO 5a: `compute_de` +
     `prep_de_side`, `p_adj < 0.05`, y fuera la fila del propio gen. Se
     COMPRUEBA gen a gen contra `p2_datos.json` (que ya se cruzo con el CLI).
     Poblaciones IDENTICAS a las ya verificadas por sha256: Jiang
     `comun.arma_linea` con `seed_real = 20260001 + fold`; 2025 `v25.selecciona`
     con 20280201/20280202.
     De la tabla se guarda `log2_fold_change`, que `p2_datos` no guardo.
  A2  Dos repartos DISTINTOS, declarados aparte y nunca sumados:
     (a) ENTRE CONTEXTOS, mismo blanco: `V_sh/V_tot` de magnitud_centro.md §4,
         ANCLADA en su control de blancos permutados. Se LEE de
         `magnitud_centro/calc.json` y `ctrl.json` (no se recalcula). Solo
         existe en Jiang: 2025 es UN contexto y ahi no se puede medir.
     (b) ENTRE PERTURBACIONES, dentro de un contexto: que fraccion del delta es
         un EJE COMUN a las perturbaciones de ese contexto.
             V_com = <mean_g h1_g, mean_g h2_g>/G     (producto cruzado de mitades)
             V_tot = (1/n) Σ_g <h1_g, h2_g>/G
         🔴 Esta es la que el brazo C puede inyectar: la inyeccion ocurre
         DENTRO de un contexto, asi que lo unico que se puede repartir es lo
         comun ENTRE PSEUDO-PERTURBACIONES. Se declara que (a) y (b) no son la
         misma cantidad.
         Su CONTROL NEGATIVO: signos ±1 por perturbacion, el MISMO en las dos
         mitades (20283900+i, 20 sorteos). Deja `V_tot` intacta y destruye el
         eje comun, asi que mide el suelo de no-comparticion del estimador
         (~1/n). La fraccion ANCLADA es (obs − control)/(1 − control).
  A4  TODO como RANGO entre superficies y entre perturbaciones. Jiang se parte
     en NUCLEO-5 y OTRAS-28: el nucleo son las 5 con dJ mediana > 0.10 en
     suelo_llamada.md §2.3 (IFNAR1, IFNAR2, JAK1, STAT1, TYK2; la siguiente es
     0.0836), que ese reporte midio como regimen distinto. El criterio es el
     corte, no la lista.
"""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/seleccion")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/ntc_nulo")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/magnitud_centro")
sys.path.insert(0, "/home/ubuntu/vcc/modelo/src")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C  # noqa: E402
import sem as S  # noqa: E402

WORK = "/datos/modelo/out/sensib"
OUT_DE = os.path.join(WORK, "a_de.json")
OUT_REP = os.path.join(WORK, "a_reparto.json")
OUT_RES = os.path.join(WORK, "a_resumen.json")
MC = "/datos/modelo/out/magnitud_centro"
MCO = "/home/ubuntu/vcc/modelo/out/magnitud_centro"
P2 = "/datos/modelo/out/seleccion/p2_datos.json"
NUCLEO = ["IFNAR1", "IFNAR2", "JAK1", "STAT1", "TYK2"]
CORTES = [0.1, 0.25, 0.5, 1.0, 2.0]
PCT = [10, 25, 50, 75, 90, 99]


def log(m):
    C.log(m)


# ───────────────────────────── A1 + A3 ─────────────────────────────
def puerta_tab(X, labels, var_names):
    """La puerta del preset, devolviendo la TABLA con log2FC y p_adj."""
    import anndata as ad
    import pandas as pd
    import polars as pl
    import scipy.sparse as sp
    from cell_eval2.de_compute import compute_de
    from cell_eval2.de import prep_de_side
    from puerta import CFG
    c = CFG
    a = ad.AnnData(X=sp.csr_matrix(X).astype(np.float32),
                   obs=pd.DataFrame({c.pert_col: list(labels)},
                                    index=[f"c{i}" for i in range(X.shape[0])]),
                   var=pd.DataFrame(index=np.asarray(var_names, dtype=str)))
    df = compute_de(a, backend=c.de.backend, groupby=c.pert_col, reference=c.control,
                    mean_calc=c.de.mean_calc, epsilon=c.de.epsilon, input_type=c.input_type,
                    target_sum=c.target_sum, clip_value=c.de.clip_value,
                    filter_gene_min_cpm_cell=c.filter.filter_gene_min_cpm_cell,
                    fdr_scope=c.de.fdr_scope, threads=c.num_threads,
                    replicate_col=c.de.replicate_col, device="cpu")
    tab, _ = prep_de_side(df, name="real", sort_by=c.de.sort_by,
                          nan_lfc_policy=c.de.nan_lfc_policy,
                          min_abs_log2fc=c.de.min_abs_log2fc)
    tcol = "target" if "target" in tab.columns else c.pert_col
    tab = tab.filter(pl.col("feature").cast(pl.Utf8) != pl.col(tcol).cast(pl.Utf8))
    return tab, tcol, float(c.de.p_adj_threshold)


def por_pert(tab, tcol, T, perts):
    import polars as pl
    out = {}
    for p in perts:
        sub = tab.filter(pl.col(tcol).cast(pl.Utf8) == p)
        n_prob = sub.height
        if n_prob == 0:
            out[p] = {"n_probados": 0, "n_ed": 0}
            continue
        lf = sub["log2_fold_change"].to_numpy().astype(float)
        pa = sub["p_adj"].to_numpy().astype(float)
        sig = (pa < T) & np.isfinite(lf)
        al = np.abs(lf)
        als = al[sig]
        r = {"n_probados": int(n_prob), "n_ed": int(sig.sum()),
             "frac_ed": float(sig.sum() / n_prob),
             "lfc_abs_pct_sig": {str(q): (float(np.percentile(als, q)) if als.size else None)
                                 for q in PCT},
             "lfc_abs_max_sig": float(als.max()) if als.size else None,
             "lfc_abs_mediana_todos": float(np.median(al[np.isfinite(al)])),
             "frac_sobre_corte": {str(c_): float(((als >= c_).sum()) / n_prob) for c_ in CORTES},
             "frac_sig_sobre_corte": {str(c_): (float((als >= c_).mean()) if als.size else None)
                                      for c_ in CORTES},
             "signo_neg_sig": float((lf[sig] < 0).mean()) if als.size else None}
        out[p] = r
    return out


def de_jiang(res):
    from suelo_arnes import load_global_obs, load_var_names, IFNB_PATH
    p2 = json.load(open(P2))
    panel = p2["jiang"]["panel"]
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = np.asarray(load_var_names(IFNB_PATH), dtype=str)
    for i, L in enumerate(C.LINES):
        if L in res["jiang"]:
            continue
        t = time.time()
        real_ad, _p, _c, _m = C.arma_linea(L, i, panel, list(var_names), cell_type, gene,
                                           len(var_names))
        real_ad.obs["target_gene"] = real_ad.obs["target_gene"].replace(
            {C.CONTROL_LABEL: C.CONTROL_COMP})
        h = C.huella_adata(real_ad)
        lab = real_ad.obs["target_gene"].to_numpy().astype(str)
        X = real_ad.X.tocsr()
        del real_ad, _p
        tab, tcol, T = puerta_tab(X, lab, var_names)
        pp = por_pert(tab, tcol, T, panel)
        ref = p2["jiang"]["crudo"][L]["ed"]
        coincide = sum(1 for g in panel if pp[g]["n_ed"] == len(ref.get(g, [])))
        res["jiang"][L] = {"huella": h, "huella_p2": p2["jiang"]["crudo"][L]["huella_real"],
                           "huella_coincide": h == p2["jiang"]["crudo"][L]["huella_real"],
                           "p_adj_threshold": T, "por_pert": pp,
                           "perts_con_|ED|_igual_a_p2": coincide, "n_perts": len(panel),
                           "dt_s": time.time() - t}
        log(f"A1 {L}: |ED| coincide con 5a en {coincide}/33, dt={time.time()-t:.0f}s")
        C.escribe_json(OUT_DE, res)


def de_2025(res):
    import v25
    from suelo_arnes import read_rows_csr
    if "v2025" in res:
        return
    t = time.time()
    p2 = json.load(open(P2))
    tg, _ip, shape, var_names, *_ = v25.lee_obs(v25.SPLIT_SUP)
    s = v25.selecciona(tg)
    X = read_rows_csr(v25.ruta(v25.SPLIT_SUP), s["idx"], shape[1]).astype(np.int64)
    lab = np.asarray(tg[s["idx"]], dtype=object).astype(str)
    var_names = np.asarray(var_names, dtype=str)
    tab, tcol, T = puerta_tab(X, lab, var_names)
    pp = por_pert(tab, tcol, T, s["perts"])
    ref = p2["v2025"]["ed"]
    coincide = sum(1 for g in s["perts"] if pp[g]["n_ed"] == len(ref.get(g, [])))
    res["v2025"] = {"p_adj_threshold": T, "por_pert": pp, "perts": s["perts"],
                    "perts_con_|ED|_igual_a_p2": coincide, "n_perts": len(s["perts"]),
                    "truncadas": [d["pert"] for d in s["por_pert"] if d["truncada"]],
                    "dt_s": time.time() - t}
    log(f"A1 2025: |ED| coincide con 5a en {coincide}/{len(s['perts'])}, dt={time.time()-t:.0f}s")
    C.escribe_json(OUT_DE, res)


# ───────────────────────────── A2 ─────────────────────────────
def reparto_entre_contextos():
    """(a) Se LEE de magnitud_centro. No se recalcula nada."""
    calc = json.load(open(os.path.join(MCO, "calc.json")))
    ctrl = json.load(open(os.path.join(MCO, "ctrl.json")))
    return {"fuente": [os.path.join(MCO, "calc.json"), os.path.join(MCO, "ctrl.json")],
            "claves_calc": list(calc)[:40], "claves_ctrl": list(ctrl)[:40],
            "calc": {k: v for k, v in calc.items() if not isinstance(v, (list,))
                     or len(json.dumps(v)) < 200000},
            "ctrl": ctrl}


def reparto_entre_perts_jiang():
    """(b) del npz de deltas: (8 particiones, 33 perts, 2 mitades, 33,992 genes)."""
    out = {}
    for L in C.LINES:
        z = np.load(os.path.join(MC, f"deltas_{L}.npz"))
        D = z["D"].astype(np.float64)
        panel = [str(x) for x in z["panel"]]
        nj, ng, _, G = D.shape
        filas = []
        for j in range(nj):
            h1, h2 = D[j, :, 0, :], D[j, :, 1, :]
            v_tot = float(np.mean(np.einsum("gi,gi->g", h1, h2) / G))
            m1, m2 = h1.mean(0), h2.mean(0)
            v_com = float(m1 @ m2 / G)
            ctr = []
            for i in range(20):
                rng = np.random.default_rng(S.s_boot(i))
                sg = rng.choice([-1.0, 1.0], size=ng)[:, None]
                a1, a2 = h1 * sg, h2 * sg
                ctr.append(float(a1.mean(0) @ a2.mean(0) / G) / v_tot if v_tot else None)
            c = [x for x in ctr if x is not None]
            obs = v_com / v_tot if v_tot else None
            anc = ((obs - np.mean(c)) / (1 - np.mean(c))) if (obs is not None and c) else None
            filas.append({"j": j, "V_tot": v_tot, "V_com": v_com, "V_esp": v_tot - v_com,
                          "frac_obs": obs, "control_medio": float(np.mean(c)) if c else None,
                          "control_sd": float(np.std(c)) if c else None,
                          "suelo_analitico_1_sobre_n": 1.0 / ng,
                          "frac_anclada": float(anc) if anc is not None else None})
        out[L] = {"n_perts": ng, "n_genes_espacio": int(G), "por_particion": filas,
                  "frac_anclada_mediana": float(np.median([f["frac_anclada"] for f in filas])),
                  "frac_anclada_rango": [float(min(f["frac_anclada"] for f in filas)),
                                         float(max(f["frac_anclada"] for f in filas))],
                  "frac_obs_mediana": float(np.median([f["frac_obs"] for f in filas]))}
        log(f"A2b {L}: frac comun anclada {out[L]['frac_anclada_mediana']:+.4f} "
            f"rango {out[L]['frac_anclada_rango']}")
    return out


def reparto_entre_perts_2025():
    """(b) en 2025: deltas propios, mitades de 200 y control PARTIDO, con el
    mismo espacio (`bulk_lognorm`, target_sum 50000, panel fuera, centrado)."""
    import anndata as ad
    import pandas as pd
    import v25
    from suelo_arnes import read_rows_csr
    from cell_eval2.prep import pseudobulk_bulk_lognorm
    from cell_eval2.distances import resolve_panel_columns, panel_reduced
    tg, _ip, shape, var_names, *_ = v25.lee_obs(v25.SPLIT_SUP)
    s = v25.selecciona(tg)
    X = read_rows_csr(v25.ruta(v25.SPLIT_SUP), s["idx"], shape[1]).astype(np.int64)
    lab = np.asarray(tg[s["idx"]], dtype=object).astype(str)
    var_names = np.asarray(var_names, dtype=str)
    perts = [p for p in s["perts"]]
    cols = resolve_panel_columns(np.asarray(perts), var_names, target_gene_map=None)
    G = panel_reduced(var_names[None, :], cols)[0].size
    idx_c = np.flatnonzero(lab == "non-targeting")
    idx_p = {p: np.flatnonzero(lab == p) for p in perts}
    filas = []
    for j in range(4):
        rng = np.random.default_rng(S.s_part2025(j))
        grp = np.empty(X.shape[0], dtype=object)
        pc = rng.permutation(idx_c)
        n1 = len(pc) // 2 + (len(pc) % 2 if j % 2 == 0 else 0)
        grp[pc[:n1]] = "NT|1"; grp[pc[n1:]] = "NT|2"
        for p in perts:
            pp = rng.permutation(idx_p[p])
            k = pp.size // 2 + (pp.size % 2 if j % 2 == 0 else 0)
            grp[pp[:k]] = f"{p}|A"; grp[pp[k:]] = f"{p}|B"
        tmp = ad.AnnData(X=X, obs=pd.DataFrame({"grp": grp.astype(str)}),
                         var=pd.DataFrame(index=var_names))
        gs, means = pseudobulk_bulk_lognorm(tmp, "grp", bulk_target_sum=50000.0)
        del tmp
        gs = np.asarray(gs).astype(str)
        red = panel_reduced(means, cols)
        pos = {g: i for i, g in enumerate(gs)}
        h1 = np.stack([red[pos[f"{p}|A"]] - red[pos["NT|1"]] for p in perts])
        h2 = np.stack([red[pos[f"{p}|B"]] - red[pos["NT|2"]] for p in perts])
        h1 = h1 - h1.mean(axis=1, keepdims=True)
        h2 = h2 - h2.mean(axis=1, keepdims=True)
        v_tot = float(np.mean(np.einsum("gi,gi->g", h1, h2) / G))
        v_com = float(h1.mean(0) @ h2.mean(0) / G)
        ctr = []
        for i in range(20):
            r2 = np.random.default_rng(S.s_boot(i))
            sg = r2.choice([-1.0, 1.0], size=len(perts))[:, None]
            ctr.append(float((h1 * sg).mean(0) @ (h2 * sg).mean(0) / G) / v_tot)
        obs = v_com / v_tot
        filas.append({"j": j, "semilla": S.s_part2025(j), "V_tot": v_tot, "V_com": v_com,
                      "frac_obs": obs, "control_medio": float(np.mean(ctr)),
                      "control_sd": float(np.std(ctr)),
                      "suelo_analitico_1_sobre_n": 1.0 / len(perts),
                      "frac_anclada": float((obs - np.mean(ctr)) / (1 - np.mean(ctr)))})
        log(f"A2b 2025 j={j}: obs {obs:+.4f} control {np.mean(ctr):+.4f} "
            f"anclada {filas[-1]['frac_anclada']:+.4f}")
    return {"n_perts": len(perts), "n_genes_espacio": int(G), "por_particion": filas,
            "frac_anclada_mediana": float(np.median([f["frac_anclada"] for f in filas])),
            "frac_anclada_rango": [float(min(f["frac_anclada"] for f in filas)),
                                   float(max(f["frac_anclada"] for f in filas))]}


# ───────────────────── A4: rangos y el puente al brazo C ─────────────────────
def rango(vals):
    v = [x for x in vals if x is not None and np.isfinite(x)]
    if not v:
        return None
    return {"n": len(v), "min": float(min(v)), "p25": float(np.percentile(v, 25)),
            "mediana": float(np.median(v)), "p75": float(np.percentile(v, 75)),
            "max": float(max(v))}


def resumen():
    de = json.load(open(OUT_DE))
    rep = json.load(open(OUT_REP))
    out = {"grupos": {}, "nucleo": NUCLEO, "criterio_nucleo":
           "dJ mediana > 0.10 en suelo_llamada.md §2.3 (las 5; la siguiente es 0.0836)"}
    grupos = {}
    for L in C.LINES:
        if L not in de.get("jiang", {}):
            continue
        pp = de["jiang"][L]["por_pert"]
        grupos[f"jiang_{L}_nucleo5"] = [pp[g] for g in NUCLEO if g in pp]
        grupos[f"jiang_{L}_otras28"] = [v for g, v in pp.items() if g not in NUCLEO]
    if "v2025" in de:
        grupos["v2025_100"] = list(de["v2025"]["por_pert"].values())
    for k, lst in grupos.items():
        out["grupos"][k] = {
            "n_perts": len(lst),
            "frac_ed": rango([x.get("frac_ed") for x in lst]),
            "n_ed": rango([x.get("n_ed") for x in lst]),
            "lfc_mediana_sig": rango([(x.get("lfc_abs_pct_sig") or {}).get("50") for x in lst]),
            "lfc_p90_sig": rango([(x.get("lfc_abs_pct_sig") or {}).get("90") for x in lst]),
            **{f"frac_sobre_{c_}": rango([(x.get("frac_sobre_corte") or {}).get(str(c_))
                                          for x in lst]) for c_ in CORTES}}
    com = {}
    for L, v in rep.get("entre_perts_jiang", {}).items():
        com[f"jiang_{L}"] = v["frac_anclada_mediana"]
    if "entre_perts_2025" in rep:
        com["v2025"] = rep["entre_perts_2025"]["frac_anclada_mediana"]
    out["frac_comun_anclada_por_superficie"] = com
    out["frac_comun_anclada_rango"] = rango(list(com.values()))

    # ── el puente al brazo C: tres puntos de estructura MEDIDA, declarados
    amp = rango([v["lfc_mediana_sig"]["mediana"] for v in out["grupos"].values()
                 if v.get("lfc_mediana_sig")])
    fr = rango([v["frac_ed"]["mediana"] for v in out["grupos"].values() if v.get("frac_ed")])
    cm = out["frac_comun_anclada_rango"]
    def pt(nombre, a, f, c):
        return {"nombre": nombre, "amplitud": round(float(a), 4),
                "frac_genes": round(float(f), 5), "frac_comun": round(float(max(0.0, min(1.0, c))), 4)}
    out["puntos_C"] = [pt("baja", amp["min"], fr["min"], cm["min"]),
                       pt("mediana", amp["mediana"], fr["mediana"], cm["mediana"]),
                       pt("alta", amp["max"], fr["max"], cm["max"])]
    out["puente_declarado"] = (
        "amplitud <- mediana de |log2FC| de los genes significativos de una "
        "perturbacion; frac_genes <- |ED|/n_probados; frac_comun <- la fraccion "
        "ANCLADA del eje comun ENTRE PERTURBACIONES (A2b), que es la unica "
        "inyectable dentro de un contexto. Los tres puntos son el MIN, la "
        "MEDIANA y el MAX de esos rangos sobre los 13 grupos (6 lineas x "
        "nucleo/otras + 2025), no una eleccion.")
    C.escribe_json(OUT_RES, out)
    log(f"A resumen: puntos C = {json.dumps(out['puntos_C'])}")
    return out


def main():
    os.makedirs(WORK, exist_ok=True)
    modo = sys.argv[1]
    if modo == "de":
        res = json.load(open(OUT_DE)) if os.path.exists(OUT_DE) else {"jiang": {}}
        de_jiang(res)
        de_2025(res)
    elif modo == "reparto":
        res = json.load(open(OUT_REP)) if os.path.exists(OUT_REP) else {}
        if "entre_contextos" not in res:
            res["entre_contextos"] = reparto_entre_contextos()
            C.escribe_json(OUT_REP, res)
        if "entre_perts_jiang" not in res:
            res["entre_perts_jiang"] = reparto_entre_perts_jiang()
            C.escribe_json(OUT_REP, res)
        if "entre_perts_2025" not in res:
            res["entre_perts_2025"] = reparto_entre_perts_2025()
            C.escribe_json(OUT_REP, res)
    elif modo == "resumen":
        resumen()


if __name__ == "__main__":
    main()
