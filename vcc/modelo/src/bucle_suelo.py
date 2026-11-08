#!/usr/bin/env python
"""Bucle de 6 pliegues (leave-one-line-out) sobre IFNB, puntuando el SUELO
con cell_eval2 0.20.0. Incremental, con las guardas de disco/tiempo/limpieza
declaradas en el encargo. Escribe cada pliegue apenas se puntua, en
modelo/out/bucle_suelo.json (JSON-lines, un objeto por linea -- asi un fallo
cuesta un pliegue, no todos).
"""
import json
import os
import shutil
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import numpy as np

from suelo_arnes import (
    IFNB_PATH, build_fold_adatas, load_global_obs, load_var_names,
    run_cell_eval2, score_vs_self, read_agg_metrics, read_avg_score,
    VCC2026_METRICS,
)
from suelo_floor import predict_floor

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT_JSONL = "/home/ubuntu/vcc/modelo/out/bucle_suelo.json"
DATOS_OUT = "/datos/modelo/out"
MIN_FREE_GIB = 10.0
TIME_BUDGET_S = 2 * 3600

SEED_BASE_REAL = 20260001
SEED_BASE_PRED = 20261001


def free_gib(path):
    total, used, free = shutil.disk_usage(path)
    return free / (1024 ** 3)


def append_record(record):
    with open(OUT_JSONL, "a") as fh:
        fh.write(json.dumps(record) + "\n")


def run_one_fold(fold_idx, line, panel, cell_type, gene, var_names, n_genes):
    t_fold0 = time.time()
    free_before = free_gib(DATOS_OUT)
    if free_before < MIN_FREE_GIB:
        raise SystemExit(
            f"PARAR: solo {free_before:.1f} GiB libres en /datos (<{MIN_FREE_GIB}) antes "
            f"del pliegue {line}."
        )

    seed_real = SEED_BASE_REAL + fold_idx
    seed_pred = SEED_BASE_PRED + fold_idx

    real_ad, pred_ad, build_meta = build_fold_adatas(
        line, panel, predict_floor, var_names, seed_real, seed_pred,
        cell_type, gene, n_genes, IFNB_PATH,
    )

    real_path = os.path.join(DATOS_OUT, f"fold_{line}_real.h5ad")
    pred_path = os.path.join(DATOS_OUT, f"fold_{line}_pred.h5ad")
    real_ad.write_h5ad(real_path)
    pred_ad.write_h5ad(pred_path)

    run_outdir = os.path.join(DATOS_OUT, f"run_{line}")
    agg_path, meta_path, run_dt = run_cell_eval2(pred_path, real_path, run_outdir)
    scored_path = score_vs_self(agg_path, meta_path, run_outdir)

    metrics = read_agg_metrics(agg_path)
    avg_score = read_avg_score(scored_path)

    # limpieza: borrar los .h5ad grandes de este pliegue (guarda b)
    for p in (real_path, pred_path):
        if os.path.exists(p):
            os.remove(p)

    t_fold_total = time.time() - t_fold0

    record = {
        "fold_idx": fold_idx,
        "linea_excluida": line,
        "n_control_cells": build_meta["n_control_cells"],
        "n_control_guides": None,  # se completa desde suelo_panel.json en el orquestador
        "seed_real": seed_real,
        "seed_pred": seed_pred,
        "real_shape": build_meta["real_shape"],
        "pred_shape": build_meta["pred_shape"],
        "run_duration_s": round(run_dt, 2),
        "fold_duration_s": round(t_fold_total, 2),
        "free_gib_before": round(free_before, 2),
        **{m: metrics.get(m) for m in VCC2026_METRICS},
        "avg_score": avg_score,
    }
    return record


def main():
    os.makedirs(DATOS_OUT, exist_ok=True)
    panel_data = json.load(open(PANEL_PATH))
    panel = panel_data["panel_comun"]
    per_line_control = panel_data["per_line_control"]
    sys.stderr.write(f"panel comun: {len(panel)} perturbaciones\n")

    sys.stderr.write("leyendo obs global de IFNB (cell_type, gene)...\n")
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    sys.stderr.write(f"n_obs={len(cell_type)} n_genes={n_genes}\n")

    # reanudar: no repetir pliegues ya escritos
    done_lines = set()
    if os.path.exists(OUT_JSONL):
        with open(OUT_JSONL) as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln:
                    continue
                done_lines.add(json.loads(ln)["linea_excluida"])
    if done_lines:
        sys.stderr.write(f"pliegues ya escritos, se saltan: {sorted(done_lines)}\n")

    fold1_duration = None
    for fold_idx, line in enumerate(LINES):
        if line in done_lines:
            continue
        sys.stderr.write(f"\n=== pliegue {fold_idx+1}/6: excluyendo {line} ===\n")
        t0 = time.time()
        try:
            record = run_one_fold(fold_idx, line, panel, cell_type, gene, var_names, n_genes)
        except Exception as e:
            sys.stderr.write(f"FALLO en pliegue {line}: {e}\n")
            append_record({
                "fold_idx": fold_idx, "linea_excluida": line, "error": str(e),
            })
            raise
        record["n_control_guides"] = per_line_control[line]["n_control_guides"]
        dt = time.time() - t0
        sys.stderr.write(
            f"pliegue {line} OK en {dt:.1f}s -- avg_score={record['avg_score']:.6f}\n"
        )
        append_record(record)

        if fold_idx == 0:
            fold1_duration = dt
            proyeccion_total = dt * len(LINES)
            sys.stderr.write(
                f"pliegue 1 tardo {dt:.1f}s -- proyeccion de los 6: {proyeccion_total:.1f}s "
                f"({proyeccion_total/3600:.2f} h)\n"
            )
            if proyeccion_total > TIME_BUDGET_S:
                sys.stderr.write(
                    "PARAR: la proyeccion supera las 2 horas. Deteniendo el bucle tras el "
                    "pliegue 1, como pide la guarda de tiempo.\n"
                )
                sys.exit(3)

    sys.stderr.write("\nbucle completo: 6 pliegues escritos en " + OUT_JSONL + "\n")


if __name__ == "__main__":
    main()
