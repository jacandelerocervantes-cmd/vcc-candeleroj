#!/usr/bin/env python
"""PASO 2 / parte 1: los deltas por (linea, particion, blanco, mitad).

LA CONSTRUCCION, declarada antes de medir
-----------------------------------------
* **Datos**: Jiang IFNB, las seis lineas, el panel COMUN de 33 perturbaciones
  (100 % de solape), de `suelo_panel.json`.
* **Poblacion real fija por linea**: `build_fold_adatas` con
  `seed_real = 20260001 + fold`. 400 celulas por perturbacion (sin reemplazo) y
  TODO el pool de control de la linea. Se guarda la huella sha256 del objeto.
* **El espacio**: `bulk_lognorm` con `bulk_target_sum = 50000.0`, leido de
  `run_params.yaml` de las corridas ya hechas, via
  `cell_eval2.prep.pseudobulk_bulk_lognorm`. Los 33 genes blanco del panel
  SALEN del espacio (`resolve_panel_columns` + `panel_reduced`, que es lo que
  hace `exclusion_scope="panel"`). Quedan 33,992 genes de 34,025.
  🔴 Si el espacio cambia, el reparto cambia.
* **El delta**: `delta_g,L = pb(g) - pb(control)` en ese espacio, y luego
  **centrado sobre los genes** (como `b_contra_a.md` §5), asi que `<h1,h2>/G`
  es una covarianza.

LAS MITADES (TRAMPA 1)
----------------------
* **Unidad: la celula.** Las 400 celulas perturbadas de `g` en `L` se permutan
  con la semilla de particion y se parten **200 / 200**. 400 es par, asi que
  aqui NO hay desempate de impares.
* 🔴 **El control va PARTIDO en dos mitades disjuntas**, no compartido. Es el
  confundidor que `b_contra_a.md` §4 destapo citando `anchor.py`: con un control
  comun el ruido de muestreo del control entra identico en las dos cantidades
  cuyo acuerdo se mide y lo infla (midio 0.54 -> 0.74 de fiabilidad). Asi que
  `delta_g,L,A = pb(mitad_A) - pb(c1)` y `delta_g,L,B = pb(mitad_B) - pb(c2)`.
* **Desempate de impares del control, declarado y ALTERNADO**: si el pool de
  control tiene `n` impar, la celula sobrante va a `c1` cuando `j` es par y a
  `c2` cuando `j` es impar.
* **Semilla de particion**: familia NUEVA `20278000 + 100*fold + j`, j = 0..7.
  Comprobado sin choque en `mc_semillas.py` (`semillas.json`).

NO pasa por la puerta de DE: esto es sobre la MAGNITUD del delta.
NO se leyo nada de ~/vcc/benchmark/ ni de /datos/benchmark_VEDADO_*.
"""
import hashlib, json, os, sys
import numpy as np
import scipy.sparse as sp

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src")
from suelo_arnes import (IFNB_PATH, CONTROL_LABEL, build_fold_adatas,
                         load_global_obs, load_var_names)
import gwps_via as GV

LINES = GV.LINES
PANEL_PATH = GV.PANEL_PATH
BULK_TARGET_SUM = 50000.0
SEED_BASE_REAL = 20260001
SEED_BASE_PART = 20278000
N_PART = 8
N_PER_LABEL = 400
WORK = "/datos/modelo/out/magnitud_centro"
OUTDIR = "/home/ubuntu/vcc/modelo/out/magnitud_centro"


def predict_trivial():
    def fn(control_pool_X, labels, rng):
        n = control_pool_X.shape[0]
        return {l: control_pool_X[rng.integers(0, n, size=N_PER_LABEL)].tocsr()
                for l in labels}
    return fn


def main():
    from cell_eval2.prep import pseudobulk_bulk_lognorm
    from cell_eval2.distances import resolve_panel_columns, panel_reduced

    panel = json.load(open(PANEL_PATH))["panel_comun"]
    panel = sorted(map(str, panel))
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = np.asarray(load_var_names(IFNB_PATH)).astype(str)
    n_genes = len(var_names)

    # el espacio reducido es el MISMO en las seis lineas: mismo panel, mismo var
    cols = resolve_panel_columns(np.asarray(panel), var_names, target_gene_map=None)
    var_r = panel_reduced(var_names[None, :], cols)[0]
    G = var_r.size
    sys.stderr.write(f"espacio: {G} genes de {n_genes} "
                     f"({n_genes - G} del panel fuera); panel {len(panel)} perts\n")

    meta_all = {"espacio": {"n_genes_total": int(n_genes), "n_genes_espacio": int(G),
                            "n_panel_fuera": int(n_genes - G),
                            "bulk_target_sum": BULK_TARGET_SUM,
                            "comparator": "bulk_lognorm",
                            "centrado_sobre_genes": True},
                "panel": panel, "lineas": {}}

    for line in (sys.argv[1:] or LINES):
        fold = LINES.index(line)
        sys.stderr.write(f"\n===== {line} (fold {fold}) =====\n")
        real_ad, _, meta = build_fold_adatas(
            line, panel, predict_trivial(), var_names.tolist(),
            SEED_BASE_REAL + fold, 999_700 + fold,
            cell_type, gene, n_genes, IFNB_PATH)
        GV.renombra(real_ad)
        huella = GV.huella(real_ad)
        lab = real_ad.obs["target_gene"].to_numpy().astype(str)
        X = real_ad.X.tocsr()
        n_ctrl = int((lab == GV.CONTROL_COMP).sum())
        idx_ctrl = np.flatnonzero(lab == GV.CONTROL_COMP)
        idx_pert = {g: np.flatnonzero(lab == g) for g in panel}
        for g in panel:
            assert idx_pert[g].size == N_PER_LABEL, (g, idx_pert[g].size)
        sys.stderr.write(f"  real {X.shape}, control {n_ctrl} celulas, "
                         f"huella {huella[:16]}...\n")
        del real_ad

        D = np.zeros((N_PART, len(panel), 2, G), dtype=np.float32)
        info_part = []
        for j in range(N_PART):
            sem = SEED_BASE_PART + 100 * fold + j
            rng = np.random.default_rng(sem)
            grp = np.empty(X.shape[0], dtype=object)

            # --- control partido, desempate de impares ALTERNADO por paridad de j
            pc = rng.permutation(idx_ctrl)
            if n_ctrl % 2 == 0:
                n1 = n_ctrl // 2
                desempate = "par: 50/50 exacto"
            elif j % 2 == 0:
                n1 = n_ctrl // 2 + 1
                desempate = "impar, j par: la sobrante a c1"
            else:
                n1 = n_ctrl // 2
                desempate = "impar, j impar: la sobrante a c2"
            c1, c2 = pc[:n1], pc[n1:]
            grp[c1] = "NT|1"; grp[c2] = "NT|2"

            # --- perturbadas: 200 / 200, sin desempate (400 es par)
            for g in panel:
                pp = rng.permutation(idx_pert[g])
                grp[pp[:200]] = f"{g}|A"
                grp[pp[200:]] = f"{g}|B"

            import pandas as pd
            import anndata as ad
            tmp = ad.AnnData(X=X, obs=pd.DataFrame({"grp": grp.astype(str)}),
                             var=pd.DataFrame(index=var_names))
            gs, means = pseudobulk_bulk_lognorm(tmp, "grp",
                                                bulk_target_sum=BULK_TARGET_SUM)
            del tmp
            gs = np.asarray(gs).astype(str)
            red = panel_reduced(means, cols)
            pos = {g: i for i, g in enumerate(gs)}
            ctrl1 = red[pos["NT|1"]]
            ctrl2 = red[pos["NT|2"]]
            for gi, g in enumerate(panel):
                for hi, (suf, cr) in enumerate((("A", ctrl1), ("B", ctrl2))):
                    d = red[pos[f"{g}|{suf}"]] - cr
                    d = d - d.mean()                 # centrado sobre genes
                    D[j, gi, hi] = d.astype(np.float32)
            info_part.append({
                "j": j, "semilla": sem,
                "n_control_c1": int(c1.size), "n_control_c2": int(c2.size),
                "desempate_impares": desempate,
                "n_pert_mitad_A": 200, "n_pert_mitad_B": 200,
                "n_grupos_pseudobulk": int(gs.size),
            })
            sys.stderr.write(f"  j={j} sem={sem} c1={c1.size} c2={c2.size} "
                             f"grupos={gs.size} ({desempate})\n")

        np.savez_compressed(os.path.join(WORK, f"deltas_{line}.npz"),
                            D=D, panel=np.asarray(panel), var_r=var_r)
        meta_all["lineas"][line] = {
            "fold": fold, "seed_real": SEED_BASE_REAL + fold,
            "huella_real_sha256": huella,
            "n_control_cells": n_ctrl,
            "n_celulas_real": int(X.shape[0]),
            "per_pert_pool_sizes": meta["per_pert_pool_sizes"],
            "particiones": info_part,
        }
        del X, D
        json.dump(meta_all, open(os.path.join(OUTDIR, "deltas_meta.json"), "w"),
                  indent=1, sort_keys=True, default=str)

    json.dump(meta_all, open(os.path.join(OUTDIR, "deltas_meta.json"), "w"),
              indent=1, sort_keys=True, default=str)
    sys.stderr.write(f"\nlisto: {WORK}/deltas_*.npz + {OUTDIR}/deltas_meta.json\n")


if __name__ == "__main__":
    main()
