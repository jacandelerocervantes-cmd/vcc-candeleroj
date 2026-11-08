#!/usr/bin/env python
"""Parte 2 — construye, para la condicion IFNB (elegida en la parte 1-e), los
6 pliegues (uno por linea excluida). Para cada linea escribe:

  /datos/modelo/out/bucle_suelo/IFNB/<LINEA>/real.h5ad
      controles NT de la linea (todos, sin modificar) + hasta 400 celulas
      reales por perturbacion (todas si hay menos de 400 -- documentado).
  /datos/modelo/out/bucle_suelo/IFNB/<LINEA>/pred.h5ad
      los MISMOS controles NT (sin modificar) + 400 celulas por perturbacion,
      muestreadas CON REEMPLAZO de esos controles, sin tocar ningun conteo.

Solo lee /datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad. No lee
nada de benchmark/. No construye ningun modelo -- el suelo es un resample.
"""
import sys
import json
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

SRC = "/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad"
OUT_ROOT = "/datos/modelo/out/bucle_suelo/IFNB"
CONTROL_LABEL = "NT"
CELLS_PER_PERT = 400
BASE_SEED = 20261029  # fecha de este encargo, documentado

ADOPTION_REPORT_PATH = "/home/ubuntu/vcc/modelo/out/bucle_suelo_build_counts.json"


def main():
    sys.stderr.write(f"[cargando] {SRC} en memoria (sin backed; 11.8 GiB en disco)\n")
    adata = ad.read_h5ad(SRC)
    sys.stderr.write(f"[cargado] shape={adata.shape}\n")

    cell_type = adata.obs["cell_type"].astype(str).to_numpy()
    gene = adata.obs["gene"].astype(str).to_numpy()
    lines = sorted(set(cell_type))
    assert len(lines) == 6, lines

    X = adata.X.tocsr()
    var = adata.var.copy()

    report = {"lines": {}}

    for line_i, line in enumerate(lines):
        line_mask = cell_type == line
        control_mask = line_mask & (gene == CONTROL_LABEL)
        control_idx = np.flatnonzero(control_mask)
        n_control = len(control_idx)

        perts = sorted(set(gene[line_mask & (gene != CONTROL_LABEL)]))
        n_perts = len(perts)

        rng = np.random.default_rng(BASE_SEED + line_i)

        pred_blocks = [X[control_idx, :]]
        pred_labels = [CONTROL_LABEL] * n_control
        real_blocks = [X[control_idx, :]]
        real_labels = [CONTROL_LABEL] * n_control

        per_pert_detail = {}
        for pert in perts:
            real_idx_full = np.flatnonzero(line_mask & (gene == pert))
            n_avail = len(real_idx_full)

            if n_avail >= CELLS_PER_PERT:
                real_sample_idx = rng.choice(real_idx_full, size=CELLS_PER_PERT, replace=False)
                shortfall = False
            else:
                real_sample_idx = real_idx_full  # todas las disponibles -- documentado
                shortfall = True

            pred_sample_idx = rng.choice(control_idx, size=CELLS_PER_PERT, replace=True)

            real_blocks.append(X[real_sample_idx, :])
            real_labels += [pert] * len(real_sample_idx)
            pred_blocks.append(X[pred_sample_idx, :])
            pred_labels += [pert] * CELLS_PER_PERT

            per_pert_detail[pert] = {
                "n_avail_real": int(n_avail),
                "n_usado_real": int(len(real_sample_idx)),
                "shortfall_bajo_400": shortfall,
                "n_pred": CELLS_PER_PERT,
            }

        X_real = sp.vstack(real_blocks, format="csr")
        X_pred = sp.vstack(pred_blocks, format="csr")

        obs_real = pd.DataFrame({"gene": pd.Categorical(real_labels)})
        obs_pred = pd.DataFrame({"gene": pd.Categorical(pred_labels)})
        obs_real.index = [f"{line}_real_{i:06d}" for i in range(len(real_labels))]
        obs_pred.index = [f"{line}_pred_{i:06d}" for i in range(len(pred_labels))]

        adata_real = ad.AnnData(X=X_real, obs=obs_real, var=var)
        adata_pred = ad.AnnData(X=X_pred, obs=obs_pred, var=var)

        outdir = f"{OUT_ROOT}/{line}"
        import os
        os.makedirs(outdir, exist_ok=True)
        real_path = f"{outdir}/real.h5ad"
        pred_path = f"{outdir}/pred.h5ad"
        adata_real.write_h5ad(real_path, compression="gzip", compression_opts=4)
        adata_pred.write_h5ad(pred_path, compression="gzip", compression_opts=4)

        n_real_total = X_real.shape[0]
        n_pred_total = X_pred.shape[0]
        sys.stderr.write(
            f"[{line}] n_control={n_control} n_perts={n_perts} "
            f"n_real_total={n_real_total} n_pred_total={n_pred_total} "
            f"shortfalls={sum(1 for v in per_pert_detail.values() if v['shortfall_bajo_400'])} "
            f"-> {real_path} / {pred_path}\n"
        )

        report["lines"][line] = {
            "n_control": int(n_control),
            "n_perts": int(n_perts),
            "n_real_total": int(n_real_total),
            "n_pred_total": int(n_pred_total),
            "n_shortfalls": sum(1 for v in per_pert_detail.values() if v["shortfall_bajo_400"]),
            "real_path": real_path,
            "pred_path": pred_path,
            "per_pert_detail": per_pert_detail,
            "seed": BASE_SEED + line_i,
        }

    with open(ADOPTION_REPORT_PATH, "w") as fh:
        json.dump(report, fh, indent=2)
    sys.stderr.write(f"[done] reporte en {ADOPTION_REPORT_PATH}\n")


if __name__ == "__main__":
    main()
