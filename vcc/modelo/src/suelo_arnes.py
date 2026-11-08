#!/usr/bin/env python
"""Arnes reutilizable: dado (funcion de prediccion, linea excluida, condicion)
devuelve las seis metricas vcc2026 y el avg_score, puntuando con cell_eval2
0.20.0 (~/celleval2_venv). El suelo (resample con reemplazo de los controles)
es UNA funcion de prediccion entre otras que este arnes acepta -- no un caso
especial.

Contrato de la funcion de prediccion:
    predict_fn(control_pool_X, labels, rng) -> dict[label] = csr_matrix (400, n_genes)
donde control_pool_X es SOLO el perfil de control de la linea (nunca su
identidad ni su nombre), labels es la lista de etiquetas a predecir
(perturbaciones del panel + "NT"), y rng es un numpy Generator ya fijado con
semilla documentada.
"""
import json
import os
import subprocess
import sys
import time
from collections import Counter

import h5py
import numpy as np
import scipy.sparse as sp
import anndata as ad
import pandas as pd

IFNB_PATH = "/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad"
CONTROL_LABEL = "NT"
N_PER_LABEL = 400
CELLEVAL2_PY = os.path.expanduser("~/celleval2_venv/bin/cell-eval2")

VCC2026_METRICS = [
    "pds_cosine",
    "expr_mse_unbiased_capped_norm",
    "de_wilcoxon_lfc_nmae",
    "de_wilcoxon_direction_fidelity_yield_raw",
    "de_wilcoxon_direction_reach_raw",
    "de_wilcoxon_sig_jaccard",
]


def decode_categorical(node):
    codes = node["codes"][:]
    cats_raw = node["categories"][:]
    cats = [c.decode() if isinstance(c, bytes) else c for c in cats_raw]
    return [cats[c] if c >= 0 else None for c in codes]


def load_global_obs(h5path=IFNB_PATH):
    with h5py.File(h5path, "r") as f:
        n_obs = int(f["X"].attrs["shape"][0])
        cell_type = decode_categorical(f["obs"]["cell_type"])
        gene = decode_categorical(f["obs"]["gene"])
    cell_type = np.array(cell_type, dtype=object)
    gene = np.array(gene, dtype=object)
    assert len(cell_type) == n_obs == len(gene)
    return cell_type, gene


def _group_runs(sorted_idx):
    """sorted_idx: array 1D ordenado ascendente, unico. Devuelve lista de
    (start, stop) de bloques contiguos (stop exclusivo)."""
    runs = []
    start = sorted_idx[0]
    prev = sorted_idx[0]
    for v in sorted_idx[1:]:
        if v == prev + 1:
            prev = v
            continue
        runs.append((int(start), int(prev) + 1))
        start = v
        prev = v
    runs.append((int(start), int(prev) + 1))
    return runs


def read_rows_csr(h5path, row_idx_sorted_unique, n_genes):
    """Lee filas arbitrarias (ordenadas, unicas) de un X CSR en disco sin
    densificar nunca la matriz completa. Agrupa en bloques contiguos para
    minimizar llamadas a h5py."""
    runs = _group_runs(row_idx_sorted_unique)
    with h5py.File(h5path, "r") as f:
        indptr_full = f["X"]["indptr"]
        data_ds = f["X"]["data"]
        indices_ds = f["X"]["indices"]
        all_indptr = indptr_full[:]  # int32, (n_obs+1,) -- liviano, nunca se toca /X con esto
        data_chunks = []
        indices_chunks = []
        new_indptr = [0]
        for start, stop in runs:
            byte_start = int(all_indptr[start])
            byte_stop = int(all_indptr[stop])
            d = data_ds[byte_start:byte_stop]
            ix = indices_ds[byte_start:byte_stop]
            data_chunks.append(d)
            indices_chunks.append(ix)
            row_lens = all_indptr[start + 1:stop + 1] - all_indptr[start:stop]
            new_indptr.extend((np.cumsum(row_lens) + new_indptr[-1]).tolist())
        data_all = np.concatenate(data_chunks) if data_chunks else np.array([], dtype=np.int64)
        indices_all = np.concatenate(indices_chunks) if indices_chunks else np.array([], dtype=np.int32)
        indptr_all = np.array(new_indptr, dtype=np.int64)
    n_rows = len(row_idx_sorted_unique)
    return sp.csr_matrix((data_all, indices_all, indptr_all), shape=(n_rows, n_genes))


def load_var_names(h5path=IFNB_PATH):
    with h5py.File(h5path, "r") as f:
        node = f["var"]["_index"]
        if hasattr(node, "keys"):
            attrs = dict(node.attrs)
            enc = attrs.get("encoding-type")
            if enc == "categorical":
                codes = node["codes"][:]
                cats_raw = node["categories"][:]
                cats = [c.decode() if isinstance(c, bytes) else c for c in cats_raw]
                return [cats[c] if c >= 0 else None for c in codes]
            elif enc == "nullable-string-array":
                vals = node["values"][:]
                mask = node["mask"][:] if "mask" in node.keys() else None
                decoded = [v.decode() if isinstance(v, bytes) else v for v in vals]
                if mask is not None:
                    decoded = [None if m else v for v, m in zip(decoded, mask)]
                return decoded
            raise ValueError(f"encoding desconocido para var/_index: {enc}")
        vals = node[:]
        return [v.decode() if isinstance(v, bytes) else v for v in vals]


def load_line_block(line, cell_type, gene, n_genes, h5path=IFNB_PATH):
    """Para una linea: lee (sin densificar) TODAS sus filas (control +
    perturbadas) de IFNB, agrupando en los bloques contiguos que ya forma esa
    linea en el archivo. Devuelve (csr_local, gene_local)."""
    idx = np.where(cell_type == line)[0]
    csr_local = read_rows_csr(h5path, idx, n_genes)
    gene_local = gene[idx]
    return csr_local, gene_local


def build_fold_adatas(line, panel, predict_fn, var_names, seed_real, seed_pred,
                       cell_type, gene, n_genes, h5path=IFNB_PATH):
    """Construye VERDAD y SUELO (o cualquier prediccion) para una linea
    excluida, condicion IFNB, panel comun ya fijado. Devuelve (real_ad,
    pred_ad, meta_construccion)."""
    csr_local, gene_local = load_line_block(line, cell_type, gene, n_genes, h5path)

    control_mask = gene_local == CONTROL_LABEL
    control_pool = csr_local[control_mask]
    n_control_cells = control_pool.shape[0]

    rng_real = np.random.default_rng(seed_real)
    real_blocks = []
    real_labels = []
    per_pert_pool_sizes = {}
    for pert in panel:
        pert_mask = gene_local == pert
        pool_idx_local = np.where(pert_mask)[0]
        per_pert_pool_sizes[pert] = int(len(pool_idx_local))
        assert len(pool_idx_local) >= N_PER_LABEL, (
            f"{line}/{pert}: pool={len(pool_idx_local)} < {N_PER_LABEL} -- el panel comun "
            "deberia garantizar esto, revisar suelo_panel.json"
        )
        chosen = rng_real.choice(pool_idx_local, size=N_PER_LABEL, replace=False)
        real_blocks.append(csr_local[chosen])
        real_labels.extend([pert] * N_PER_LABEL)
    # control real: TODO el pool de control de la linea, sin submuestrear ni modificar
    real_blocks.append(control_pool)
    real_labels.extend([CONTROL_LABEL] * n_control_cells)

    real_X = sp.vstack(real_blocks, format="csr")
    real_obs = pd.DataFrame({"target_gene": real_labels})
    real_ad = ad.AnnData(X=real_X, obs=real_obs, var=pd.DataFrame(index=var_names))
    real_ad.obs_names_make_unique()

    rng_pred = np.random.default_rng(seed_pred)
    pred_dict = predict_fn(control_pool, panel + [CONTROL_LABEL], rng_pred)
    pred_blocks = []
    pred_labels = []
    for label in panel + [CONTROL_LABEL]:
        block = pred_dict[label]
        assert block.shape[0] == N_PER_LABEL, (label, block.shape)
        pred_blocks.append(block)
        pred_labels.extend([label] * N_PER_LABEL)
    pred_X = sp.vstack(pred_blocks, format="csr")
    pred_obs = pd.DataFrame({"target_gene": pred_labels})
    pred_ad = ad.AnnData(X=pred_X, obs=pred_obs, var=pd.DataFrame(index=var_names))
    pred_ad.obs_names_make_unique()

    meta = {
        "line": line,
        "n_control_cells": int(n_control_cells),
        "per_pert_pool_sizes": per_pert_pool_sizes,
        "seed_real": seed_real,
        "seed_pred": seed_pred,
        "real_shape": list(real_ad.shape),
        "pred_shape": list(pred_ad.shape),
    }
    return real_ad, pred_ad, meta


def run_cell_eval2(pred_path, real_path, outdir, pert_col="target_gene", control=CONTROL_LABEL):
    os.makedirs(outdir, exist_ok=True)
    cmd = [
        CELLEVAL2_PY, "run",
        "-ap", pred_path, "-ar", real_path,
        "--preset", "vcc2026", "--profile", "vcc2026",
        "--pert-col", pert_col, "--control", control,
        "-o", outdir,
    ]
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    log_path = os.path.join(outdir, "cell_eval2_run.log")
    with open(log_path, "w") as fh:
        fh.write("CMD: " + " ".join(cmd) + "\n\n")
        fh.write("STDOUT:\n" + proc.stdout + "\n\nSTDERR:\n" + proc.stderr + "\n")
        fh.write(f"\nEXITCODE:{proc.returncode}\nDURATION_S:{dt:.2f}\n")
    if proc.returncode != 0:
        raise RuntimeError(
            f"cell-eval2 run fallo (exit={proc.returncode}) para {outdir}. Ver {log_path}"
        )
    agg_path = os.path.join(outdir, "agg_results.csv")
    meta_path = os.path.join(outdir, "run_meta.json")
    assert os.path.exists(agg_path), f"no se escribio {agg_path}"
    assert os.path.exists(meta_path), f"no se escribio {meta_path}"
    return agg_path, meta_path, dt


def score_vs_self(agg_path, meta_path, outdir):
    out_csv = os.path.join(outdir, "scored.csv")
    cmd = [
        CELLEVAL2_PY, "score",
        "--user-agg", agg_path, "--baseline-agg", agg_path,
        "--baseline-meta", meta_path, "--user-meta", meta_path,
        "--allow-config-mismatch",
        "-o", out_csv,
    ]
    log_path = os.path.join(outdir, "cell_eval2_score.log")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    with open(log_path, "w") as fh:
        fh.write("CMD: " + " ".join(cmd) + "\n\n")
        fh.write("STDOUT:\n" + proc.stdout + "\n\nSTDERR:\n" + proc.stderr + "\n")
        fh.write(f"\nEXITCODE:{proc.returncode}\n")
    if proc.returncode != 0:
        raise RuntimeError(f"cell-eval2 score fallo (exit={proc.returncode}). Ver {log_path}")
    return out_csv


def read_agg_metrics(agg_path):
    """agg_results.csv es ancho: una fila por estadistico (count, mean, std,
    min, max, median, null_count), una columna por metrica. El valor
    reportado es la fila 'mean' -- la misma que usa `score`."""
    df = pd.read_csv(agg_path)
    row = df.loc[df["statistic"] == "mean"].iloc[0]
    return {m: float(row[m]) for m in VCC2026_METRICS if m in df.columns}


def read_avg_score(scored_path):
    """scored.csv de `score` es largo: columnas metric,from_baseline, con una
    fila final metric=='avg_score'."""
    df = pd.read_csv(scored_path)
    row = df.loc[df["metric"] == "avg_score"]
    if row.empty:
        raise KeyError(f"no se encontro fila avg_score en {scored_path}: {df}")
    return float(row["from_baseline"].iloc[0])
