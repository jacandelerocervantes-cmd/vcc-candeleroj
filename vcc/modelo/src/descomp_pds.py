#!/usr/bin/env python
"""Parte 2 de la descomposicion: mide la dispersion entre perfiles de
perturbacion en las TRES emisiones (nuestro suelo, el b oficial, el real), para
decidir si el escalado de baseline.py:344 destruye discriminacion.

Criterio declarado ANTES de medir: las SEIS lineas, sin seleccion.

`pds_cosine` rankea el DELTA (efecto) predicho contra los deltas reales
(`metrics/discrimination.py:50-53`), y el comparador de vcc2026 es
`bulk_lognorm` con `bulk_target_sum=50000` (`configs/vcc2026.yaml`). Asi que la
dispersion se mide sobre exactamente eso: `prep.pseudobulk_bulk_lognorm` y
`prep.delta`, las funciones del propio marcador.

No modifica nada. Solo lee y cuenta.
"""
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (
    IFNB_PATH, CONTROL_LABEL, CELLEVAL2_PY,
    build_fold_adatas, load_global_obs, load_var_names,
)
from suelo_floor import predict_floor

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT_JSON = "/home/ubuntu/vcc/modelo/out/descomp_pds.json"
WORK = "/datos/modelo/out/descomp"
CONTROL_COMP = "non-targeting"
SEED_BASE_REAL = 20260001
SEED_BASE_PRED = 20261001
BULK_TARGET_SUM = 50000.0


def renombra(adata):
    adata.obs["target_gene"] = adata.obs["target_gene"].replace(
        {CONTROL_LABEL: CONTROL_COMP}
    )
    return adata


def dispersion(adata, etiqueta):
    """Dispersion entre los deltas de perturbacion, con las funciones del
    marcador. Devuelve:
      norma_media      : |delta| medio -- cuanto se aparta del control cada pert
      cos_par_media    : distancia cosine media entre PARES de deltas -- cuanto
                         se separan las perturbaciones ENTRE SI. 0 = todas en la
                         misma direccion (sin discriminacion posible).
      frac_norma_cero  : fraccion de deltas con norma ~0 (un delta de norma cero
                         empata toda su fila a distancia 1.0 bajo cosine, ver
                         discrimination.py:66-70)
    """
    from cell_eval2.prep import pseudobulk_bulk_lognorm, delta
    perts, means = pseudobulk_bulk_lognorm(
        adata, "target_gene", bulk_target_sum=BULK_TARGET_SUM
    )
    perts_nc, eff = delta(means, perts, CONTROL_COMP)
    eff = np.asarray(eff, dtype=np.float64)
    normas = np.linalg.norm(eff, axis=1)
    nz = normas > 0
    E = eff[nz] / normas[nz][:, None]
    G = E @ E.T
    n = G.shape[0]
    if n > 1:
        iu = np.triu_indices(n, k=1)
        cos_d = 1.0 - G[iu]
        cos_par_media = float(cos_d.mean())
        cos_par_min = float(cos_d.min())
        cos_par_max = float(cos_d.max())
    else:
        cos_par_media = cos_par_min = cos_par_max = float("nan")
    return {
        "etiqueta": etiqueta,
        "n_perts": int(eff.shape[0]),
        "norma_media": float(normas.mean()),
        "norma_min": float(normas.min()),
        "norma_max": float(normas.max()),
        "frac_norma_cero": float((~nz).sum() / normas.size),
        "cos_par_media": cos_par_media,
        "cos_par_min": cos_par_min,
        "cos_par_max": cos_par_max,
    }


def main():
    os.makedirs(WORK, exist_ok=True)
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    sys.stderr.write(f"panel={len(panel)} n_genes={n_genes}\n")

    salida = {}
    if os.path.exists(OUT_JSON):
        salida = json.load(open(OUT_JSON))

    for fold_idx, line in enumerate(LINES):
        if line in salida:
            sys.stderr.write(f"{line}: ya hecho, se salta\n")
            continue
        t0 = time.time()
        sys.stderr.write(f"\n=== {line} ===\n")
        real_ad, pred_ad, _ = build_fold_adatas(
            line, panel, predict_floor, var_names,
            SEED_BASE_REAL + fold_idx, SEED_BASE_PRED + fold_idx,
            cell_type, gene, n_genes, IFNB_PATH,
        )
        renombra(real_ad)
        renombra(pred_ad)

        real_path = os.path.join(WORK, f"real_{line}.h5ad")
        real_ad.write_h5ad(real_path)

        res = {}
        sys.stderr.write("  dispersion del REAL... ")
        res["real"] = dispersion(real_ad, "real (la verdad)")
        sys.stderr.write("ok\n")
        sys.stderr.write("  dispersion de NUESTRO SUELO... ")
        res["suelo"] = dispersion(pred_ad, "nuestro suelo (control sin escalar)")
        sys.stderr.write("ok\n")
        del real_ad, pred_ad

        # La emision del b oficial: hay que pedirsela a `cell-eval2 baseline`
        # con --save-pred, porque es el unico que construye el perfil y la escala.
        base_dir = os.path.join(WORK, f"BASE_{line}")
        base_pred = os.path.join(WORK, f"basepred_{line}.h5ad")
        os.makedirs(base_dir, exist_ok=True)
        cmd = [
            CELLEVAL2_PY, "baseline", "-ar", real_path,
            "--preset", "vcc2026", "--profile", "vcc2026",
            "--pert-col", "target_gene",
            "--emit", "dispersed", "--seed", "0",
            "--save-pred", base_pred, "-o", base_dir,
        ]
        sys.stderr.write("  cell-eval2 baseline --save-pred... ")
        sys.stderr.flush()
        p = subprocess.run(cmd, capture_output=True, text=True)
        with open(os.path.join(WORK, f"baseline_{line}.log"), "w") as fh:
            fh.write("CMD: " + " ".join(cmd) + "\n\nSTDOUT:\n" + p.stdout
                     + "\n\nSTDERR:\n" + p.stderr + f"\nEXIT:{p.returncode}\n")
        if p.returncode != 0:
            raise RuntimeError(f"baseline fallo para {line}; ver el log")
        sys.stderr.write("ok\n")

        import anndata as ad
        b_ad = ad.read_h5ad(base_pred)
        sys.stderr.write("  dispersion del B OFICIAL... ")
        res["b_oficial"] = dispersion(b_ad, "b oficial (control escalado al perfil)")
        sys.stderr.write("ok\n")
        del b_ad

        for q in (real_path, base_pred):
            if os.path.exists(q):
                os.remove(q)

        res["duracion_s"] = round(time.time() - t0, 2)
        salida[line] = res
        with open(OUT_JSON, "w") as fh:
            json.dump(salida, fh, indent=2, sort_keys=True)
        sys.stderr.write(
            f"  [{line}] cos_par_media -- real={res['real']['cos_par_media']:.6f} "
            f"suelo={res['suelo']['cos_par_media']:.6f} "
            f"b={res['b_oficial']['cos_par_media']:.6f}\n"
        )

    sys.stderr.write(f"\nlisto: {OUT_JSON}\n")


if __name__ == "__main__":
    main()
