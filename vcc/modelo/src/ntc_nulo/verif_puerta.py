#!/usr/bin/env python
"""Verificacion del arnes de la puerta, contra la DEFINICION del preset:
 1. el conjunto de genes probados == {genes con CPM medio de la REFERENCIA > 5.0}
 2. `p_adj` == Benjamini-Hochberg de `p_value` DENTRO de ese conjunto
 3. el conjunto DE == filas con `p_adj < 0.05`
 4. el log2FC == log2((media_CPM_blanco + eps)/(media_CPM_ref + eps))
Si los cuatro cierran, el arnes no esta reimplementando la puerta: la esta corriendo.
"""
import sys, os
import numpy as np, scipy.sparse as sp, h5py, polars as pl
from scipy.stats import false_discovery_control
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from puerta import corre_puerta, CFG

P = "/datos/reto/final/context_D.h5ad"
with h5py.File(P, "r") as f:
    shape = tuple(int(x) for x in f["X"].attrs["shape"])
    X = sp.csr_matrix((f["X/data"][:], f["X/indices"][:], f["X/indptr"][:]), shape=shape)
    vn = np.array([x.decode() if isinstance(x, bytes) else x
                   for x in f["var/_index/values"][:]])
    cods = f["obs/ntc_id/codes"][:].astype(np.int64)
    cats = np.array([x.decode() if isinstance(x, bytes) else x
                     for x in f["obs/ntc_id/categories"][:]])

orden = np.argsort(cats)
rng = np.random.default_rng(20274300)
perm = rng.permutation(46)
g1 = set(orden[perm[:23]].tolist()); g2 = set(orden[perm[23:]].tolist())
i1 = np.where(np.isin(cods, list(g1)))[0]
i2 = np.where(np.isin(cods, list(g2)))[0]
s1 = np.sort(rng.choice(i1, 400, replace=False))
s2 = np.sort(rng.choice(i2, 400, replace=False))

ns, npr, ex, tabla = corre_puerta(X[s1], X[s2], vn, devuelve_tabla=True)
print(f"corrida: n_sig={ns} m={npr}")


def cpm_medio(Xb):
    Xw = Xb.astype(np.float64)
    L = np.asarray(Xw.sum(axis=1)).ravel(); L[L <= 0] = 1.0
    return np.asarray(Xw.T @ ((1e6 / L) / Xb.shape[0])).ravel()


ref = cpm_medio(X[s1]); tgt = cpm_medio(X[s2])
esperados = set(vn[ref > CFG.filter.filter_gene_min_cpm_cell].tolist())
obtenidos = set(tabla["feature"].to_list())
print(f"1. conjunto probado: esperado {len(esperados)}, obtenido {len(obtenidos)}, "
      f"IGUAL = {esperados == obtenidos}")
if esperados != obtenidos:
    print("   sobran:", sorted(obtenidos - esperados)[:5],
          " faltan:", sorted(esperados - obtenidos)[:5])

p = tabla["p_value"].to_numpy().astype(float)
padj_obt = tabla["p_adj"].to_numpy().astype(float)
val = ~np.isnan(p)
padj_esp = np.full(p.shape, np.nan)
padj_esp[val] = false_discovery_control(p[val], method="bh")
d = np.nanmax(np.abs(padj_esp - padj_obt))
print(f"2. p_adj == BH(p_value) dentro del conjunto: max|dif| = {d:.3e} -> "
      f"IGUAL = {d < 1e-12}")

n_sig_esp = int((padj_esp[val] < CFG.de.p_adj_threshold).sum())
print(f"3. conjunto DE: n_sig del arnes = {ns}, recontado de BH = {n_sig_esp} -> "
      f"IGUAL = {ns == n_sig_esp}")

pos = {g: i for i, g in enumerate(vn)}
gi = np.array([pos[g] for g in tabla["feature"].to_list()])
lfc_esp = np.log2((tgt[gi] + CFG.de.epsilon) / (ref[gi] + CFG.de.epsilon))
lfc_obt = tabla["log2_fold_change"].to_numpy().astype(float)
dl = np.nanmax(np.abs(lfc_esp - lfc_obt))
print(f"4. log2FC == log2((tgt+eps)/(ref+eps)): max|dif| = {dl:.3e} -> "
      f"IGUAL = {dl < 1e-6}")
print(f"   (eps = {CFG.de.epsilon}, mean_calc = {CFG.de.mean_calc}, "
      f"target_sum = {CFG.target_sum})")
