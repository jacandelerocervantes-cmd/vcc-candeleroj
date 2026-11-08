#!/usr/bin/env python
"""PARTE 2c: cuanto del nulo lo explica la IDENTIDAD DE LA GUIA.

Dos mediciones, las dos sin partir nunca una guia:

(1) PARTICION QUE AGRUPA GUIAS contra PARTICIONES ALEATORIAS DE GUIAS.
    La particion aleatoria de la Parte 1 reparte las 46 guias al azar, asi que
    la diferencia entre mitades solo recoge el ruido de muestreo de guias. La
    particion que AGRUPA las ordena por su propio perfil y corta por la mitad:
    las 23 guias de un extremo contra las 23 del otro. Es el peor caso que la
    estructura de guia permite SIN partir ninguna guia.
      * perfil de guia = media aritmetica del CPM por celula sobre sus 400
        celulas, restringido al panel probable (CPM medio del pool > 5.0),
        log1p, centrado por gen.
      * se ordena por la PRIMERA componente principal de esa matriz 46 x G.
      * mitad REFERENCIA = las 23 de PC1 mas baja; BLANCO = las 23 mas alta.
    Determinista: no lleva semilla. Se corre en los mismos tres brazos de `n`.

(2) DESCOMPOSICION DE VARIANZA POR GEN (ANOVA de un factor, 46 guias x 400).
    En el espacio log1p(CPM) -- el mismo que ve el Wilcoxon antes de rangos:
        MSB = SSB/45,  MSW = SSW/(18400-46),  F = MSB/MSW
        ICC = (MSB - MSW) / (MSB + 399*MSW)
    `F` mide si la identidad de la guia mueve la media del gen mas de lo que
    mueve el muestreo de 400 celulas; `ICC` mide QUE FRACCION de la varianza
    por celula se la lleva la guia. Se reporta la distribucion sobre el panel
    probable y cuantos genes pasan BH al 0.05 del preset.
"""
import json, os, sys
import numpy as np
import scipy.sparse as sp
import h5py
from scipy.stats import f as fdist
from scipy.stats import false_discovery_control

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from puerta import corre_puerta, CFG, params

RETO = "/datos/reto"
CTX = [("A", "val", 0), ("B", "val", 1), ("C", "val", 2),
       ("D", "final", 3), ("E", "final", 4), ("F", "final", 5)]
OUT = "/home/ubuntu/vcc/modelo/out/ntc_nulo/p2c.json"
N_ENTREGA = 400
UMBRAL_CPM = CFG.filter.filter_gene_min_cpm_cell
ALPHA = CFG.de.p_adj_threshold
SEM_SUB = 20274500          # solo para submuestrear a 400 en esta particion


def carga(letra, part):
    p = os.path.join(RETO, part, f"context_{letra}.h5ad")
    with h5py.File(p, "r") as f:
        shape = tuple(int(x) for x in f["X"].attrs["shape"])
        X = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]),
                          shape=shape)
        vn = np.array([x.decode() if isinstance(x, bytes) else x
                       for x in f["var/_index/values"][:]])
        cats = np.array([x.decode() if isinstance(x, bytes) else x
                         for x in f["obs/ntc_id/categories"][:]])
        cods = f["obs/ntc_id/codes"][:].astype(np.int64)
    return X, vn, cats, cods


def cpm_log1p(X):
    """CPM por celula (x/L*1e6) y luego log1p, en disperso."""
    Xw = X.astype(np.float64).tocsr()
    L = np.asarray(Xw.sum(axis=1)).ravel(); L[L <= 0] = 1.0
    Y = sp.diags(1e6 / L) @ Xw
    Y = Y.tocsr()
    Y.data = np.log1p(Y.data)
    return Y


def cpm_medio_pool(X):
    Xw = X.astype(np.float64)
    L = np.asarray(Xw.sum(axis=1)).ravel(); L[L <= 0] = 1.0
    return np.asarray(Xw.T @ ((1e6 / L) / X.shape[0])).ravel()


res = {"puerta": params(), "umbral_cpm": UMBRAL_CPM, "alpha": ALPHA,
       "contextos": []}

for letra, part, k in CTX:
    sys.stderr.write(f"\n===== 2c contexto {letra} =====\n")
    X, vn, cats, cods = carga(letra, part)
    N, G = X.shape
    pool_cpm = cpm_medio_pool(X)
    panel = pool_cpm > UMBRAL_CPM
    sys.stderr.write(f"  panel probable: {int(panel.sum())} de {G}\n")

    Y = cpm_log1p(X)
    # indicador 46 x N de guia
    Gi = sp.csr_matrix((np.ones(N), (cods, np.arange(N))), shape=(46, N))
    S = np.asarray((Gi @ Y).todense())                 # 46 x G sumas
    Y2 = Y.copy(); Y2.data = Y2.data ** 2
    S2 = np.asarray((Gi @ Y2).todense())               # 46 x G sumas de cuadrados
    del Y2
    n_g = np.asarray(Gi.sum(axis=1)).ravel()           # 400 cada una
    tot = S.sum(axis=0); tot2 = S2.sum(axis=0)
    SSB = (S ** 2 / n_g[:, None]).sum(axis=0) - tot ** 2 / N
    SSW = tot2 - (S ** 2 / n_g[:, None]).sum(axis=0)
    df_b, df_w = 45, N - 46
    MSB, MSW = SSB / df_b, SSW / df_w
    with np.errstate(divide="ignore", invalid="ignore"):
        F = MSB / MSW
        ICC = (MSB - MSW) / (MSB + (400 - 1) * MSW)
    pF = fdist.sf(F, df_b, df_w)
    ok = panel & np.isfinite(pF) & (MSW > 0)
    padj = np.full(G, np.nan)
    padj[ok] = false_discovery_control(pF[ok], method="bh")
    n_sig_guia = int((padj[ok] < ALPHA).sum())
    iccp = ICC[ok]

    # ---- particion que AGRUPA guias: PC1 de los perfiles de guia
    M = np.log1p(S[:, panel] / n_g[:, None] * 0.0 + (S[:, panel] / n_g[:, None]))
    # S/n ya es la media de log1p(CPM); no se vuelve a log1p
    M = S[:, panel] / n_g[:, None]
    Mc = M - M.mean(axis=0, keepdims=True)
    U, sv, Vt = np.linalg.svd(Mc, full_matrices=False)
    pc1 = U[:, 0] * sv[0]
    ordpc = np.argsort(pc1)
    g1 = set(ordpc[:23].tolist()); g2 = set(ordpc[23:].tolist())
    i1 = np.where(np.isin(cods, list(g1)))[0]
    i2 = np.where(np.isin(cods, list(g2)))[0]
    assert len(g1 & g2) == 0 and len(g1 | g2) == 46
    assert i1.size == 9200 and i2.size == 9200
    rng = np.random.default_rng(SEM_SUB + k)
    s1 = np.sort(rng.choice(i1, N_ENTREGA, replace=False))
    s2 = np.sort(rng.choice(i2, N_ENTREGA, replace=False))

    corridas = []
    for etq, ir, it in (("completo", i1, i2), ("n400_400", s1, s2),
                        ("n400_ref_completa", i1, s2)):
        ns, npr, ex, _ = corre_puerta(X[ir], X[it], vn)
        corridas.append({"brazo": etq, "n_ref": int(ir.size), "n_tgt": int(it.size),
                         "n_sig": ns, "n_probados": npr, "p_min": ex["p_min"],
                         "abs_lfc_max": ex["abs_lfc_max"],
                         "abs_lfc_p99": ex["abs_lfc_p99"],
                         "genes_sig": ex["genes_sig"]})
        sys.stderr.write(f"  [AGRUPA {etq:18s}] n_sig={ns} probados={npr} "
                         f"p_min={ex['p_min']:.3e}\n")

    res["contextos"].append({
        "contexto": letra, "particion": "agrupa_guias_PC1",
        "n_panel_probable": int(panel.sum()), "n_genes": int(G),
        "pc1_varianza_explicada": float(sv[0] ** 2 / (sv ** 2).sum()),
        "pc1_rango": [float(pc1.min()), float(pc1.max())],
        "guias_ref_pc1_baja": sorted(cats[sorted(g1)].tolist()),
        "guias_tgt_pc1_alta": sorted(cats[sorted(g2)].tolist()),
        "guia_partida": False, "disjuntas": True, "union_46": True,
        "corridas": corridas,
        "anova": {
            "n_genes_evaluados": int(ok.sum()),
            "n_genes_guia_significativa_BH": n_sig_guia,
            "frac_genes_guia_significativa": round(float(n_sig_guia / ok.sum()), 5),
            "F_mediana": float(np.median(F[ok])),
            "F_p99": float(np.percentile(F[ok], 99)),
            "F_max": float(np.nanmax(F[ok])),
            "ICC_mediana": float(np.median(iccp)),
            "ICC_p75": float(np.percentile(iccp, 75)),
            "ICC_p99": float(np.percentile(iccp, 99)),
            "ICC_max": float(np.nanmax(iccp)),
            "ICC_frac_negativo": float(np.mean(iccp < 0)),
            "df_entre": df_b, "df_dentro": int(df_w),
        },
    })
    del X, Y, S, S2
    json.dump(res, open(OUT, "w"), indent=1, sort_keys=True, default=str)

json.dump(res, open(OUT, "w"), indent=1, sort_keys=True, default=str)
sys.stderr.write(f"\nlisto: {OUT}\n")
