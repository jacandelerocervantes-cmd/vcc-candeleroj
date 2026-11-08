#!/usr/bin/env python
"""Verificacion de la premisa central de la descomposicion de `jac`:
¿`jaccard == 1.0` por perturbacion equivale a `|union| == 0`?

El conteo de la descomposicion usa `count(jac == 1.0)` como E. Eso es una
INFERENCIA: un jaccard de 1.0 tambien podria venir de un acuerdo perfecto
(tp == union > 0). Aqui se mide, para una linea, con las cuentas de union del
propio marcador.

Reusa `techo_arnes.contar_de_por_perturbacion` y
`techo_arnes.run_cell_eval2_techo` SIN modificarlos (el segundo agrega
`--write-degenes`, que es lo que hace falta para tener las tablas de DE).

Emision: `suelo_floor.predict_floor`, la misma del suelo. Poblacion: la de
siempre, seed_real = 20260001 (A549).
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (
    IFNB_PATH, CONTROL_LABEL, build_fold_adatas, load_global_obs, load_var_names,
)
from suelo_floor import predict_floor
from techo_arnes import run_cell_eval2_techo, contar_de_por_perturbacion

LINE = "A549"
FOLD = 0
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT_JSON = "/home/ubuntu/vcc/modelo/out/descomp_verif.json"
WORK = "/datos/modelo/out/descomp_verif"


def main():
    os.makedirs(WORK, exist_ok=True)
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)

    # La etiqueta de control se deja en `NT`, porque `run_cell_eval2_techo`
    # pasa `--control NT` (es el arnes del techo, sin modificar). Eso no afecta
    # a lo que se verifica aqui -- la equivalencia es entre dos cantidades de la
    # MISMA corrida -- y el paso 0 de anclas_nativas_v2.md ya probo que la
    # etiqueta no mueve ningun numero.
    real_ad, pred_ad, _ = build_fold_adatas(
        LINE, panel, predict_floor, var_names, 20260001 + FOLD, 20261001 + FOLD,
        cell_type, gene, n_genes, IFNB_PATH,
    )
    real_path = os.path.join(WORK, "real.h5ad")
    pred_path = os.path.join(WORK, "pred.h5ad")
    real_ad.write_h5ad(real_path)
    pred_ad.write_h5ad(pred_path)
    del real_ad, pred_ad

    outdir = os.path.join(WORK, "run")
    sys.stderr.write("cell-eval2 run --write-degenes... ")
    sys.stderr.flush()
    run_cell_eval2_techo(pred_path, real_path, outdir, "real")
    sys.stderr.write("ok\n")

    cuentas = contar_de_por_perturbacion(outdir)

    import pandas as pd
    res = pd.read_csv(os.path.join(outdir, "results.csv"))
    jac = {str(r.perturbation): float(r.value)
           for r in res[res.metric == "de_wilcoxon_sig_jaccard"].itertuples()}

    tabla = []
    for p, c in cuentas.items():
        tabla.append({
            "pert": p,
            "n_DE_real": c["n_DE_real"],
            "n_DE_pred": c["n_DE_pred"],
            "n_union": c["n_DE_union"],
            "union_vacia": c["union_vacia"],
            "jac_marcador": jac.get(p),
        })
    tabla.sort(key=lambda d: (not d["union_vacia"], d["pert"]))

    # La tabla de contingencia 2x2 que decide la premisa.
    uno_y_vacia = sum(1 for t in tabla if t["jac_marcador"] == 1.0 and t["union_vacia"])
    uno_no_vacia = sum(1 for t in tabla if t["jac_marcador"] == 1.0 and not t["union_vacia"])
    no_uno_vacia = sum(1 for t in tabla if t["jac_marcador"] != 1.0 and t["union_vacia"])
    no_uno_no_vacia = sum(1 for t in tabla
                          if t["jac_marcador"] != 1.0 and not t["union_vacia"])
    equivale = (uno_no_vacia == 0 and no_uno_vacia == 0)

    salida = {
        "linea": LINE,
        "emision": "suelo_floor.predict_floor",
        "contingencia": {
            "jac==1.0 y union vacia": uno_y_vacia,
            "jac==1.0 y union NO vacia": uno_no_vacia,
            "jac!=1.0 y union vacia": no_uno_vacia,
            "jac!=1.0 y union NO vacia": no_uno_no_vacia,
        },
        "equivalencia_verificada": equivale,
        "n_perts": len(tabla),
        "tabla": tabla,
    }
    with open(OUT_JSON, "w") as fh:
        json.dump(salida, fh, indent=2, sort_keys=True)

    for q in (real_path, pred_path):
        if os.path.exists(q):
            os.remove(q)
    for nombre in ("de_real.parquet", "de_pred.parquet"):
        q = os.path.join(outdir, nombre)
        if os.path.exists(q):
            os.remove(q)

    sys.stderr.write(
        f"\ncontingencia: {salida['contingencia']}\n"
        f"equivalencia jac==1.0 <=> union vacia: {equivale}\n"
        f"escrito {OUT_JSON}\n"
    )


if __name__ == "__main__":
    main()
