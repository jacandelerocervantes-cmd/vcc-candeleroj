#!/usr/bin/env python
"""Paso 1, Parte 2 — construye la CARGA TRIVIAL (contratos A/B/C, partición val)
y la escribe en /datos/modelo/out/ (artefactos grandes, fuera de syncout.service
y del volumen raíz — ver /datos/modelo/out/LEEME.txt). No envía nada, no toca
benchmark/, no modifica /datos/reto ni ninguna otra ruta de origen.

Carga trivial: para cada una de las 300 perturbaciones de val/pert_counts.csv y
cada uno de los contextos A, B, C, se muestrean 400 células CON REEMPLAZO de las
18,400 células de control de ESE contexto (val/context_{A,B,C}.h5ad), copiando su
X tal cual (sin modificar ningún conteo). obs.target_gene = símbolo de la
perturbación; obs.context = la letra, reutilizada exactamente como aparece en el
archivo de control de origen.

Antes de escribir nada, calcula y reporta las entradas-no-cero (nnz) previstas,
por archivo real de indptr (sin cargar data), contra el tope del contrato
(4,750,000,000 entradas almacenadas en todo el .h5ad) y el presupuesto por
célula (~13,200 = 4.75e9 / 360,000). Si no cabe, se detiene ANTES de tocar X.
"""
import sys
import json
import argparse
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata
import csv

VAL_DIR = "/datos/reto/val"
CONTEXTS = ["A", "B", "C"]
GENE_NAMES_CSV = f"{VAL_DIR}/gene_names.csv"
PERT_COUNTS_CSV = f"{VAL_DIR}/pert_counts.csv"
CELLS_PER_PERT = 400
CAP_NNZ = 4_750_000_000
DEFAULT_OUT_H5AD = "/datos/modelo/out/prueba.h5ad"

SEEDS = {"A": 1029, "B": 1030, "C": 1031}  # documentado, reproducible


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--out",
        default=DEFAULT_OUT_H5AD,
        help=f"ruta de salida del .h5ad (default: {DEFAULT_OUT_H5AD} — "
        "artefactos grandes van a /datos, fuera de syncout.service)",
    )
    return p.parse_args()


def read_csv_col(path):
    with open(path, newline="") as fh:
        rows = list(csv.reader(fh))
    return [r[0] for r in rows[1:]]


def main():
    args = parse_args()
    out_h5ad = args.out
    gene_names = read_csv_col(GENE_NAMES_CSV)
    perts = read_csv_col(PERT_COUNTS_CSV)
    n_genes = len(gene_names)
    n_perts = len(perts)
    print(f"[info] n_genes={n_genes} n_perts={n_perts}", file=sys.stderr)
    assert n_perts == 300, f"se esperaban 300 perturbaciones, hay {n_perts}"

    # --- paso previo: nnz_per_cell por contexto, leído de indptr, sin tocar data ---
    nnz_per_cell = {}
    n_obs_ctx = {}
    for c in CONTEXTS:
        path = f"{VAL_DIR}/context_{c}.h5ad"
        with h5py.File(path, "r") as f:
            indptr = f["X"]["indptr"][:]
            nnz_per_cell[c] = np.diff(indptr).astype(np.int64)
            n_obs_ctx[c] = len(nnz_per_cell[c])
            assert f["X"]["indptr"].shape[0] - 1 == 18400

    # --- muestreo CON REEMPLAZO, determinista, documentado ---
    sampled_idx = {}
    for c in CONTEXTS:
        rng = np.random.default_rng(SEEDS[c])
        n_pool = n_obs_ctx[c]
        blocks = [rng.integers(0, n_pool, size=CELLS_PER_PERT) for _ in range(n_perts)]
        sampled_idx[c] = np.concatenate(blocks)
        assert sampled_idx[c].shape[0] == n_perts * CELLS_PER_PERT

    # --- entradas previstas, calculadas ANTES de tocar X ---
    predicted_nnz_by_ctx = {}
    for c in CONTEXTS:
        predicted_nnz_by_ctx[c] = int(nnz_per_cell[c][sampled_idx[c]].sum())
    total_predicted_nnz = sum(predicted_nnz_by_ctx.values())
    n_rows_total = n_perts * CELLS_PER_PERT * len(CONTEXTS)
    per_cell_avg = total_predicted_nnz / n_rows_total
    per_cell_budget = CAP_NNZ / n_rows_total

    report = {
        "n_rows_total": n_rows_total,
        "n_genes": n_genes,
        "predicted_nnz_by_context": predicted_nnz_by_ctx,
        "total_predicted_nnz": total_predicted_nnz,
        "cap_nnz": CAP_NNZ,
        "fits_under_cap": total_predicted_nnz <= CAP_NNZ,
        "per_cell_avg_predicted": per_cell_avg,
        "per_cell_budget_from_cap": per_cell_budget,
        "margin_factor": CAP_NNZ / total_predicted_nnz,
    }
    print(json.dumps(report, indent=2))
    sys.stdout.flush()

    if not report["fits_under_cap"]:
        print("[STOP] la carga trivial excede el tope de dispersión del contrato "
              "— no se escribe nada.", file=sys.stderr)
        sys.exit(2)

    print(f"[ok] {total_predicted_nnz:,} nnz previstas <= {CAP_NNZ:,} "
          f"(margen x{report['margin_factor']:.2f}); "
          f"{per_cell_avg:.1f} nnz/celula promedio vs presupuesto "
          f"{per_cell_budget:.1f} nnz/celula. Se procede a construir y escribir.",
          file=sys.stderr)

    # --- construir X real, por contexto, con fancy-indexing de filas ---
    X_blocks = []
    obs_target_gene = []
    obs_context = []
    obs_index = []
    umi_max_seen = 0

    for c in CONTEXTS:
        path = f"{VAL_DIR}/context_{c}.h5ad"
        adata_c = anndata.read_h5ad(path)  # en memoria, ~220MB, cabe de sobra
        assert list(adata_c.var_names) == gene_names, f"orden de genes distinto en {c}"
        X_full = adata_c.X.tocsr()
        X_sub = X_full[sampled_idx[c], :].astype(np.int32)
        row_sums = np.asarray(X_sub.sum(axis=1)).ravel()
        umi_max_seen = max(umi_max_seen, int(row_sums.max()))
        X_blocks.append(X_sub)

        for i, pert in enumerate(perts):
            obs_target_gene.extend([pert] * CELLS_PER_PERT)
            obs_context.extend([c] * CELLS_PER_PERT)
            obs_index.extend([f"{c}_{pert}_{j:04d}" for j in range(CELLS_PER_PERT)])
        del adata_c, X_full

    X_total = sp.vstack(X_blocks, format="csr")
    X_total = X_total.astype(np.int32)

    print(f"[check] X_total.shape={X_total.shape} nnz_real={X_total.nnz:,} "
          f"umi_max_observado_en_la_muestra={umi_max_seen:,}", file=sys.stderr)
    assert X_total.shape == (n_rows_total, n_genes)
    assert X_total.nnz == total_predicted_nnz, (
        f"nnz real ({X_total.nnz}) no coincide con lo previsto ({total_predicted_nnz})"
    )
    assert int(X_total.data.min()) >= 0
    assert umi_max_seen < 1_000_000

    obs = pd.DataFrame(
        {
            "target_gene": pd.Categorical(obs_target_gene),
            "context": pd.Categorical(obs_context, categories=CONTEXTS),
        },
        index=pd.Index(obs_index, name=None),
    )
    assert (obs["target_gene"] == "non-targeting").sum() == 0
    assert obs["target_gene"].nunique() == 300
    for c in CONTEXTS:
        vc = obs.loc[obs["context"] == c, "target_gene"].value_counts()
        assert len(vc) == 300 and (vc == 400).all(), f"conteo roto en contexto {c}"

    var = pd.DataFrame(index=pd.Index(gene_names, name=None))

    adata_out = anndata.AnnData(X=X_total, obs=obs, var=var)
    print(f"[check] adata_out.shape={adata_out.shape}", file=sys.stderr)

    adata_out.write_h5ad(out_h5ad, compression="gzip", compression_opts=4)
    print(f"[done] escrito {out_h5ad}", file=sys.stderr)


if __name__ == "__main__":
    main()
