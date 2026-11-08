#!/usr/bin/env python
"""LA PUERTA DE SIGNIFICANCIA DEL PRESET vcc2026, reusada SIN reimplementarla.

Todos los parametros se LEEN de `EvalConfig.from_preset("vcc2026")` -- ninguna
constante numerica se escribe aqui. El camino es el mismo que recorre el
marcador en su lado de CPU:

  1. `cell_eval2.de_compute.compute_de`
       `_to_linear`   -> normalize_total(target_sum)   (CPM)
       medias aritmeticas por grupo -> log2FC con `epsilon`
       `_log1p_view`  -> MWU del backend (`auto` -> pdex sin CUDA)
       `_apply_cpm_filter` -> conserva los genes cuyo CPM MEDIO EN EL GRUPO DE
          REFERENCIA supera `filter_gene_min_cpm_cell`, y RECALCULA BH por
          blanco sobre las filas que quedan (`_bh_per_target`)
  2. `cell_eval2.de.prep_de_side`
       `normalize_de_schema` -> `apply_nan_policy` ('mask': p_adj=1 si log2FC es
          NaN) -> `apply_lfc_floor` (no-op con min_abs_log2fc=0.0)
  3. conjunto DE = filas con `p_adj < p_adj_threshold`, el mismo filtro de
     `de._rank_matrix`.
"""
import numpy as np
import scipy.sparse as sp
import polars as pl

from cell_eval2.config import EvalConfig
from cell_eval2.de_compute import compute_de
from cell_eval2.de import prep_de_side

CFG = EvalConfig.from_preset("vcc2026")


def params():
    c = CFG
    return {
        "preset": "vcc2026",
        "pert_col": c.pert_col, "control": c.control, "version": c.version,
        "input_type": c.input_type, "target_sum": c.target_sum,
        "filter_gene_min_cpm_cell": c.filter.filter_gene_min_cpm_cell,
        "de.method": c.de.method, "de.backend": c.de.backend,
        "de.mean_calc": c.de.mean_calc, "de.epsilon": c.de.epsilon,
        "de.clip_value": c.de.clip_value, "de.fdr_scope": c.de.fdr_scope,
        "de.p_adj_threshold": c.de.p_adj_threshold,
        "de.min_abs_log2fc": c.de.min_abs_log2fc,
        "de.nan_lfc_policy": c.de.nan_lfc_policy, "de.sort_by": c.de.sort_by,
        "num_threads": c.num_threads,
    }


def corre_puerta(X_ref, X_tgt, var_names, tgt_label="mitad_B", devuelve_tabla=False):
    """Corre la puerta entre dos bloques de cuentas CRUDAS. `X_ref` es el grupo
    de REFERENCIA (el que decide el conjunto de genes de la puerta CPM, #351).
    Devuelve (n_sig, n_probados, extra)."""
    import anndata as ad
    import pandas as pd

    c = CFG
    ref_label = c.control
    X = sp.vstack([sp.csr_matrix(X_ref), sp.csr_matrix(X_tgt)], format="csr")
    if X.dtype != np.float32:
        X = X.astype(np.float32)
    obs = pd.DataFrame(
        {c.pert_col: [ref_label] * X_ref.shape[0] + [tgt_label] * X_tgt.shape[0]},
        index=[f"c{i}" for i in range(X.shape[0])],
    )
    adata = ad.AnnData(X=X, obs=obs,
                       var=pd.DataFrame(index=np.asarray(var_names, dtype=str)))

    df = compute_de(
        adata, backend=c.de.backend, groupby=c.pert_col, reference=ref_label,
        mean_calc=c.de.mean_calc, epsilon=c.de.epsilon, input_type=c.input_type,
        target_sum=c.target_sum, clip_value=c.de.clip_value,
        filter_gene_min_cpm_cell=c.filter.filter_gene_min_cpm_cell,
        fdr_scope=c.de.fdr_scope, threads=c.num_threads,
        replicate_col=c.de.replicate_col, device="cpu",
    )
    tabla, _perts = prep_de_side(df, name="nulo", sort_by=c.de.sort_by,
                                 nan_lfc_policy=c.de.nan_lfc_policy,
                                 min_abs_log2fc=c.de.min_abs_log2fc)
    sig = tabla.filter(pl.col("p_adj") < c.de.p_adj_threshold)
    lfc = tabla["log2_fold_change"].to_numpy()
    p = tabla["p_value"].to_numpy()
    extra = {
        "n_probados": int(tabla.height),
        "n_sig": int(sig.height),
        "p_min": float(np.nanmin(p)) if tabla.height else None,
        "abs_lfc_max": float(np.nanmax(np.abs(lfc))) if tabla.height else None,
        "abs_lfc_p99": float(np.nanpercentile(np.abs(lfc), 99)) if tabla.height else None,
        "abs_lfc_mediana": float(np.nanmedian(np.abs(lfc))) if tabla.height else None,
        "n_lfc_nan": int(np.isnan(lfc).sum()),
        "genes_sig": sig["feature"].to_list() if sig.height <= 200 else None,
        "umbral_bh_rango1": float(c.de.p_adj_threshold / max(tabla.height, 1)),
    }
    return (int(sig.height), int(tabla.height), extra,
            tabla if devuelve_tabla else None)
