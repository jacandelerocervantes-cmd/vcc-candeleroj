#!/usr/bin/env python
"""¿DE QUE ESTA HECHO el delta que lee `pds_cosine`? Descomposicion del delta
real en COLA (genes del ED real) y BULTO (el resto), en el MISMO espacio que
usa la metrica.

🔴 TRANSDUCTIVO Y DIAGNOSTICO: lee el lado real perturbado. No es un candidato
y ninguna de sus cifras puede entrar en una submission.

El espacio es el de `pds_cosine` en una corrida v2 counts/counts:
`bulk_lognorm` con `bulk_target_sum = 50000.0` (leido de `run_params.yaml` de
las corridas ya hechas), y con los 33 genes blanco del panel FUERA del espacio
de rasgos, que es lo que hace `exclusion_scope="panel"`.

Como `pds_cosine` usa distancia COSINE, la cantidad que decide es la direccion
del vector: se reporta el reparto de la NORMA AL CUADRADO entre cola y bulto, y
el coseno entre el delta completo y cada una de sus dos partes.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (IFNB_PATH, CONTROL_LABEL, build_fold_adatas,
                         load_global_obs, load_var_names)
import gwps_via as GV

LINES = GV.LINES
PANEL_PATH = GV.PANEL_PATH
BULK_TARGET_SUM = 50000.0          # run_params.yaml: bulk_target_sum
ED_CACHE = "/datos/modelo/out/pds_obj/ed_real.json"
OUT = "/home/ubuntu/vcc/modelo/out/pds_composicion.json"
SEED_BASE_REAL = 20260001


def predict_trivial():
    def fn(control_pool_X, labels, rng):
        n = control_pool_X.shape[0]
        return {l: control_pool_X[rng.integers(0, n, size=400)].tocsr()
                for l in labels}
    return fn


def main():
    from cell_eval2.prep import pseudobulk_bulk_lognorm
    from cell_eval2.distances import resolve_panel_columns, panel_reduced

    panel = json.load(open(PANEL_PATH))["panel_comun"]
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    ed = json.load(open(ED_CACHE))

    # 🔴 FUSIONA en vez de sobrescribir: cada linea se mide cuando su cache de
    # ED real existe, asi que el JSON se completa en varias pasadas.
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for line in (sys.argv[1:] or LINES):
        fold = LINES.index(line)
        if line not in ed:
            sys.stderr.write(f"{line}: sin ED real en cache, se omite\n")
            continue
        real_ad, _, _ = build_fold_adatas(
            line, panel, predict_trivial(), var_names,
            SEED_BASE_REAL + fold, 999_800 + fold,
            cell_type, gene, n_genes, IFNB_PATH)
        GV.renombra(real_ad)
        huella = GV.huella(real_ad)
        perts, means = pseudobulk_bulk_lognorm(
            real_ad, "target_gene", bulk_target_sum=BULK_TARGET_SUM)
        perts = np.asarray(perts).astype(str)
        del real_ad

        # los 33 genes blanco del panel salen del espacio, como hace la metrica
        no_ctrl = perts[perts != GV.CONTROL_COMP]
        cols = resolve_panel_columns(no_ctrl, np.asarray(var_names).astype(str),
                                     target_gene_map=None)
        ictrl = int(np.flatnonzero(perts == GV.CONTROL_COMP)[0])
        ctrl = means[ictrl]
        reducido = panel_reduced(means, cols)
        ctrl_r = reducido[ictrl]
        var_r = panel_reduced(np.asarray(var_names).astype(str)[None, :], cols)[0]
        pos = {g: i for i, g in enumerate(var_r)}

        filas = []
        for h in no_ctrl:
            i = int(np.flatnonzero(perts == h)[0])
            d = reducido[i] - ctrl_r
            en_ed = np.zeros(d.size, dtype=bool)
            gs = [pos[g] for g in ed[line].get(str(h), ()) if g in pos]
            en_ed[gs] = True
            n2 = float(d @ d)
            cola = d * en_ed
            bulto = d * ~en_ed
            n2c, n2b = float(cola @ cola), float(bulto @ bulto)
            nc, nb, nt = np.linalg.norm(cola), np.linalg.norm(bulto), np.linalg.norm(d)
            filas.append({
                "pert": str(h),
                "n_ed": int(len(ed[line].get(str(h), ()))),
                "n_ed_en_espacio": int(en_ed.sum()),
                "n_genes_espacio": int(d.size),
                "norma2_total": n2, "norma2_cola": n2c, "norma2_bulto": n2b,
                "frac_norma2_cola": (n2c / n2 if n2 else None),
                "cos_completo_cola": (float(d @ cola) / (nt * nc)
                                      if nt and nc else None),
                "cos_completo_bulto": (float(d @ bulto) / (nt * nb)
                                       if nt and nb else None),
            })
        fr = [f["frac_norma2_cola"] for f in filas if f["frac_norma2_cola"] is not None]
        out[line] = {
            "huella_real": huella,
            "n_genes_espacio": filas[0]["n_genes_espacio"],
            "n_panel_fuera": int(n_genes - filas[0]["n_genes_espacio"]),
            "frac_norma2_cola_media": float(np.mean(fr)),
            "frac_norma2_cola_min": float(np.min(fr)),
            "frac_norma2_cola_max": float(np.max(fr)),
            "n_perts_ed_vacio": int(sum(1 for f in filas if f["n_ed"] == 0)),
            "por_pert": filas,
        }
        sys.stderr.write(
            f"{line}: espacio {out[line]['n_genes_espacio']} genes "
            f"({out[line]['n_panel_fuera']} del panel fuera), "
            f"frac norma^2 en la cola media {out[line]['frac_norma2_cola_media']:.4f} "
            f"[{out[line]['frac_norma2_cola_min']:.4f}, "
            f"{out[line]['frac_norma2_cola_max']:.4f}], "
            f"{out[line]['n_perts_ed_vacio']} perts con ED vacio\n")
        with open(OUT, "w") as fh:
            json.dump(out, fh)
    sys.stderr.write(f"listo: {OUT}\n")


if __name__ == "__main__":
    main()
