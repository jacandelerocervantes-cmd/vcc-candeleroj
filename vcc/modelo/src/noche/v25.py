#!/usr/bin/env python
"""Comun a B4a/B4b/B4c: la superficie VCC 2025.

CONSTRUCCION DE LA POBLACION REAL 2025, DECLARADA ANTES DE MEDIR
  * split `test` (100 perturbaciones CON VERDAD; es el que mas se parece en
    forma al reto 2026, que lleva 300)
  * celulas por perturbacion: `min(400, disponibles)`. El reto 2026 fija 400
    (`manifest.json: cells_per_pert`). 🔴 15 de las 100 perturbaciones del split
    tienen menos de 400 celulas (minimo 67): se toman TODAS y se declara cual
    es cual. ESTO DIFIERE del arnes de Jiang, que exige >= 400 y por eso
    descarta perturbaciones.
  * control: las celulas `non-targeting`, submuestreadas a 18,400 -- el mismo
    numero que el reto 2026 lleva por contexto (`manifest.json: control_cells`)
  * etiqueta de control: `non-targeting`, la que exige la regla de competencia
  * semillas: 20280201 (celulas por perturbacion) y 20280202 (submuestreo de
    control), familia NUEVA
"""
import json
import os
import sys

import numpy as np
import scipy.sparse as sp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C

sys.path.insert(0, C.SRC)

RAIZ = "/datos/vcc2025"
SPLITS = {
    "train": ("train/adata_Training.h5ad", "train/pert_counts_Training.csv"),
    "validation": ("validation/adata_Validation.h5ad",
                   "validation/pert_counts_Validation.csv"),
    "test": ("test/adata_Test.h5ad", "test/pert_counts_Test.csv"),
}
GENES_2025 = os.path.join(RAIZ, "gene_names.csv")
RETO_FINAL = "/datos/reto/final"
RETO_VAL = "/datos/reto/val"
CONTROL = "non-targeting"
SPLIT_SUP = "test"
N_POR_PERT = 400
N_CONTROL = 18400
SEED_PERT = C.SEMILLAS_NUEVAS["b4_real_pert"][0]
SEED_CTRL = C.SEMILLAS_NUEVAS["b4_control_sub"][0]

WORK = os.path.join(C.WORK_BASE, "b4")
REAL_H5 = os.path.join(WORK, f"real_2025_{SPLIT_SUP}.h5ad")
REAL_META = os.path.join(WORK, f"real_2025_{SPLIT_SUP}_meta.json")


def ruta(split):
    return os.path.join(RAIZ, SPLITS[split][0])


def ruta_counts(split):
    return os.path.join(RAIZ, SPLITS[split][1])


def lee_obs(split):
    """Devuelve (target_gene como array de str, indptr, shape, var_names, cols)."""
    import h5py
    from suelo_arnes import decode_categorical
    with h5py.File(ruta(split), "r") as f:
        shape = [int(x) for x in f["X"].attrs["shape"]]
        tg = np.array(decode_categorical(f["obs"]["target_gene"]), dtype=object)
        cols = list(f["obs"].keys())
        otros = {}
        for c in ("batch", "guide_id"):
            if c in f["obs"]:
                v = f["obs"][c]
                otros[c] = (np.array(decode_categorical(v), dtype=object)
                            if hasattr(v, "keys") else None)
        indptr = f["X"]["indptr"][:].astype(np.int64)
        vnode = f["var"]["_index"]
        var_names = [x.decode() if isinstance(x, bytes) else x for x in vnode[:]]
        var_cols = list(f["var"].keys())
    return tg, indptr, shape, var_names, cols, var_cols, otros


def selecciona(tg, seed_pert=SEED_PERT, seed_ctrl=SEED_CTRL,
               n_por_pert=N_POR_PERT, n_control=N_CONTROL):
    """La seleccion de filas, reproducible y declarada."""
    perts = sorted(set(tg.tolist()) - {CONTROL})
    rng = np.random.default_rng(seed_pert)
    sel, por_pert = [], []
    for p in perts:
        idx = np.where(tg == p)[0]
        k = min(n_por_pert, len(idx))
        ch = rng.choice(idx, size=k, replace=False)
        sel.append(ch)
        por_pert.append({"pert": p, "disponibles": int(len(idx)), "tomadas": int(k),
                         "truncada": bool(k < n_por_pert)})
    sel = np.concatenate(sel)
    ctrl = np.where(tg == CONTROL)[0]
    k = min(n_control, len(ctrl))
    csub = np.random.default_rng(seed_ctrl).choice(ctrl, size=k, replace=False)
    idx_all = np.sort(np.concatenate([sel, csub]))
    return {"perts": perts, "por_pert": por_pert, "idx": idx_all,
            "n_pert_cells": int(len(sel)), "n_control_pool": int(len(ctrl)),
            "n_control_usadas": int(k), "n_filas": int(len(idx_all))}


def arma_real(log=True):
    """Construye (o reusa) el h5ad real 2025. Devuelve (ruta, meta)."""
    import anndata as ad
    import pandas as pd
    from suelo_arnes import read_rows_csr
    os.makedirs(WORK, exist_ok=True)
    if os.path.exists(REAL_H5) and os.path.exists(REAL_META):
        return REAL_H5, json.load(open(REAL_META))
    tg, indptr, shape, var_names, _c, _vc, _o = lee_obs(SPLIT_SUP)
    s = selecciona(tg)
    if log:
        C.log(f"B4 real 2025: leyendo {s['n_filas']} filas de {shape[0]}")
    X = read_rows_csr(ruta(SPLIT_SUP), s["idx"], shape[1])
    lab = tg[s["idx"]]
    adata = ad.AnnData(X=X.astype(np.int64),
                       obs=pd.DataFrame({"target_gene": list(lab)}),
                       var=pd.DataFrame(index=var_names))
    adata.obs_names_make_unique()
    meta = {
        "split": SPLIT_SUP, "shape": list(adata.shape),
        "n_perturbaciones": len(s["perts"]),
        "n_celulas_perturbadas": s["n_pert_cells"],
        "n_control_pool": s["n_control_pool"],
        "n_control_usadas": s["n_control_usadas"],
        "n_por_pert_objetivo": N_POR_PERT,
        "n_truncadas": sum(1 for x in s["por_pert"] if x["truncada"]),
        "por_pert": s["por_pert"],
        "seed_pert": SEED_PERT, "seed_control": SEED_CTRL,
        "etiqueta_control": CONTROL,
        "nnz": int(X.nnz), "dtype_origen": "float32 -> int64",
        "n_genes": int(shape[1]),
        "huella": C.huella_adata(adata),
    }
    adata.write_h5ad(REAL_H5)
    del adata
    meta["bytes_h5ad"] = os.path.getsize(REAL_H5)
    C.escribe_json(REAL_META, meta)
    if log:
        C.log(f"B4 real 2025 escrito: {meta['bytes_h5ad']/2**30:.2f} GiB, "
              f"nnz {meta['nnz']/1e6:.1f} M")
    return REAL_H5, meta


def tabla_construccion(meta):
    return C.tabla(["# ", "cantidad", "valor"], [
        [1, "split", f"`{meta['split']}` (CON verdad)"],
        [2, "perturbaciones", meta["n_perturbaciones"]],
        [3, "celulas por perturbacion", f"`min({meta['n_por_pert_objetivo']}, disponibles)`"],
        [4, "🔴 perturbaciones truncadas (< 400 celulas)", f"**{meta['n_truncadas']}** de {meta['n_perturbaciones']}"],
        [5, "celulas perturbadas", f"{meta['n_celulas_perturbadas']:,}"],
        [6, "pool de control `non-targeting`", f"{meta['n_control_pool']:,}"],
        [7, "control usado", f"{meta['n_control_usadas']:,} (el mismo numero que el reto 2026 lleva por contexto)"],
        [8, "forma del h5ad real", f"{meta['shape'][0]:,} x {meta['shape'][1]:,}"],
        [9, "`nnz`", f"{meta['nnz']:,} ({meta['nnz']/1e6:.1f} M)"],
        [10, "semillas", f"`{meta['seed_pert']}` (celulas/pert), `{meta['seed_control']}` (control) -- familia NUEVA 20280000+"],
        [11, "huella sha256 del real", meta["huella"][:32] + "..."],
        [12, "dtype", "el X de 2025 es `float32` con valores enteros; se convierte a `int64`, como el de Jiang"],
    ])
