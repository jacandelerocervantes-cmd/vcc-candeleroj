#!/usr/bin/env python
"""Arma modelo/out/bucle_suelo.md a partir de modelo/out/bucle_suelo.json
(JSON-lines, un objeto por pliegue) y modelo/out/suelo_panel.json."""
import json

JSONL_PATH = "/home/ubuntu/vcc/modelo/out/bucle_suelo.json"
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT_MD = "/home/ubuntu/vcc/modelo/out/bucle_suelo.md"

METRICS = [
    "pds_cosine",
    "expr_mse_unbiased_capped_norm",
    "de_wilcoxon_lfc_nmae",
    "de_wilcoxon_direction_fidelity_yield_raw",
    "de_wilcoxon_direction_reach_raw",
    "de_wilcoxon_sig_jaccard",
]


def load_records():
    records = []
    with open(JSONL_PATH) as fh:
        for ln in fh:
            ln = ln.strip()
            if ln:
                records.append(json.loads(ln))
    return records


def main():
    records = load_records()
    panel = json.load(open(PANEL_PATH))

    lines_md = []
    lines_md.append("# Bucle del suelo sobre 6 pliegues -- condicion IFNB\n")
    lines_md.append(
        f"Panel comun (perturbaciones con >=400 celulas en las SEIS lineas de IFNB): "
        f"**{panel['n_panel_comun']}** de {panel['n_perts_union_ifnb']} vistas en total.\n"
    )
    lines_md.append("## Descartes por linea (perturbaciones con <400 celulas EN ESA linea)\n")
    lines_md.append("| linea | n descartadas | perturbaciones |")
    lines_md.append("|---|---|---|")
    for line, d in panel["descartes_por_linea"].items():
        lines_md.append(f"| {line} | {d['n_descartadas']} | {', '.join(d['perts'])} |")

    lines_md.append("\n## Panel comun -- 33 perturbaciones, celulas por linea (intermedio)\n")
    lines_md.append("| perturbacion | " + " | ".join(panel["lineas"]) + " |")
    lines_md.append("|---" * (1 + len(panel["lineas"])) + "|")
    table = panel["per_pert_counts_table"]
    by_pert = {}
    for row in table:
        by_pert.setdefault(row["pert"], {})[row["linea"]] = row["n_cells"]
    for pert in panel["panel_comun"]:
        counts = by_pert[pert]
        lines_md.append(
            f"| {pert} | " + " | ".join(str(counts[l]) for l in panel["lineas"]) + " |"
        )

    lines_md.append("\n## La tabla -- 6 pliegues x metricas (el entregable)\n")
    header = (
        "| linea excluida | n control cells | n control guides | " + " | ".join(METRICS)
        + " | avg_score |"
    )
    lines_md.append(header)
    lines_md.append("|---" * (4 + len(METRICS)) + "|")
    for r in sorted(records, key=lambda x: x.get("fold_idx", 0)):
        if "error" in r:
            lines_md.append(f"| {r['linea_excluida']} | ERROR: {r['error']} | | | | | | | |")
            continue
        vals = [f"{r.get(m):.6f}" if r.get(m) is not None else "NaN" for m in METRICS]
        lines_md.append(
            f"| {r['linea_excluida']} | {r['n_control_cells']} | {r['n_control_guides']} | "
            + " | ".join(vals) + f" | {r['avg_score']:.6f} |"
        )

    lines_md.append("\n## Intermedios de tiempo y limpieza, por pliegue\n")
    lines_md.append("| linea | fold_duration_s | run_duration_s | free_gib_before | real_shape | pred_shape | seed_real | seed_pred |")
    lines_md.append("|---|---|---|---|---|---|---|---|")
    for r in sorted(records, key=lambda x: x.get("fold_idx", 0)):
        if "error" in r:
            continue
        lines_md.append(
            f"| {r['linea_excluida']} | {r['fold_duration_s']} | {r['run_duration_s']} | "
            f"{r['free_gib_before']} | {r['real_shape']} | {r['pred_shape']} | "
            f"{r['seed_real']} | {r['seed_pred']} |"
        )

    import numpy as np
    ctrl = np.array([r["n_control_cells"] for r in records], dtype=float)
    lines_md.append(
        "\n## Correlacion con el grosor del control (n=6, solo descriptiva)\n"
    )
    lines_md.append("| metrica | pearson r vs n_control_cells |")
    lines_md.append("|---|---|")
    for m in METRICS:
        v = np.array([r[m] for r in records], dtype=float)
        if v.std() == 0:
            lines_md.append(f"| {m} | constante en 0, sin correlacion definida |")
        else:
            r = float(np.corrcoef(ctrl, v)[0, 1])
            lines_md.append(f"| {m} | {r:.3f} |")
    lines_md.append(
        "\n`expr_mse_unbiased_capped_norm` correlaciona fuerte y negativamente "
        "(r=-0.93): con un control mas fino (menos celulas), el estimador "
        "insesgado de MSE normalizado sale mas alto (peor) -- consistente con "
        "mas ruido de muestreo en el perfil de control. `pds_cosine` "
        "correlaciona positivo (r=0.79): mas celulas de control, mejor "
        "separabilidad. Las columnas de DE (`lfc_nmae`, `sig_jaccard`) "
        "correlacionan negativo, mas debil. Con solo 6 puntos esto es "
        "descriptivo, no una prueba estadistica.\n"
    )

    lines_md.append("## avg_score = 0.0 en los seis pliegues -- no es un error\n")
    lines_md.append(
        "El suelo se puntua CONTRA SI MISMO (`--user-agg` y `--baseline-agg` son "
        "el mismo `agg_results.csv`), por el camino diagnostico de `score` que "
        "reporta `from_baseline` sin exigir ancla (`--allow-config-mismatch`, "
        "sin `--anchor` ni `--real-bundle`). En la escala vcc2026, 0 es por "
        "definicion el suelo de respuesta media -- asi que un suelo puntuado "
        "contra si mismo cae en 0.0 por construccion, en los seis pliegues. Lo "
        "informativo no es esta columna (siempre 0 aqui): son las seis columnas "
        "de metricas crudas de arriba, que SI varian por pliegue y que quedan "
        "fijadas como el `agg_results.csv` de referencia (`--baseline-agg`) para "
        "puntuar futuros candidatos contra el suelo.\n"
    )

    lines_md.append("## Decisiones de construccion, no triviales\n")
    lines_md.append(
        "- `cell_eval2.io.validate_pair` exige que PRED y REAL compartan "
        "exactamente el mismo conjunto de etiquetas de perturbacion, y que la "
        "etiqueta de control este presente en REAL. El encargo solo describe "
        "la VERDAD como \"las celulas perturbadas\" -- se le agrego el bloque "
        "de control de la linea (TODO el pool real de esa linea, sin "
        "submuestrear, sin modificar) para cumplir ese contrato de la "
        "herramienta. El bloque de control de SUELO se construye con la MISMA "
        "regla que cualquier otra etiqueta: 400 celulas remuestreadas con "
        "reemplazo del propio pool de control (semilla distinta a la de la "
        "VERDAD), no una copia de las celulas reales -- asi la funcion de "
        "prediccion no necesita un caso especial para el control.\n"
        "- El arnes (`modelo/src/suelo_arnes.py`) nunca pasa el nombre ni la "
        "identidad de la linea a la funcion de prediccion: `predict_fn(control_pool_X, "
        "labels, rng)` solo recibe la matriz de control y la lista de "
        "etiquetas a predecir. `suelo_floor.predict_floor` (modelo/src/suelo_floor.py) "
        "es una implementacion de esa firma entre las que acepta "
        "`suelo_arnes.build_fold_adatas` -- un candidato real se conecta sin "
        "tocar el bucle (`modelo/src/bucle_suelo.py`).\n"
        "- `--preset vcc2026 --profile vcc2026 --pert-col target_gene --control NT` "
        "en `cell-eval2 run`, usando `~/celleval2_venv` (0.20.0, instalado con "
        "`pdex` para DE en CPU -- sin GPU en esta VM, confirmado por "
        "`de.backend='auto'` cayendo a `pdex` en el log de cada pliegue).\n"
        "- Lectura de `/X` sin densificar nunca la matriz completa: para cada "
        "linea se agrupan sus filas (dispersas en ~32 bloques contiguos dentro "
        "del h5ad de 328,542x34,025) y se leen esos bloques con `h5py`; el "
        "submuestreo (VERDAD sin reemplazo, SUELO con reemplazo) se hace en "
        "memoria sobre ese bloque ya leido, no celda por celda contra disco.\n"
    )

    lines_md.append("## Criterio de adopcion (recordatorio, no aplica aqui)\n")
    lines_md.append(
        "Esta corrida establece el suelo, no lo compara contra nada: no hay "
        "ningun candidato todavia. El criterio declarado para cuando exista uno "
        "sigue siendo el mismo -- un candidato sustituye al suelo solo si gana "
        "en los SEIS pliegues; si gana en 5, se reporta en cual perdio; si "
        "ninguno domina en los seis, no hay metodo.\n"
    )

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(lines_md) + "\n")
    print(f"escrito {OUT_MD} con {len(records)} pliegues")


if __name__ == "__main__":
    main()
