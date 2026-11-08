#!/usr/bin/env python
"""Bucle del TECHO: 6 pliegues (leave-one-line-out) sobre IFNB x 2 variantes de
control (S compartido / I partido) x 2 direcciones (A como prediccion / B como
prediccion) = 24 corridas de cell-eval2.

Incremental: cada corrida se escribe en modelo/out/bucle_techo.json (JSON-lines)
apenas se puntua, asi un fallo cuesta una corrida y no las 24. Reanudable.
"""
import json
import os
import shutil
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (
    IFNB_PATH, VCC2026_METRICS,
    load_global_obs, load_var_names, score_vs_self,
    read_agg_metrics, read_avg_score,
)
from techo_arnes import (
    SEED_BASE_SPLIT, SEED_BASE_CTRLS,
    partir_linea, build_techo_pair, run_cell_eval2_techo,
    contar_de_por_perturbacion, leer_results_por_perturbacion,
)

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
VARIANTES = ["S", "I"]
DIRECCIONES = ["A", "B"]

PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT_JSONL = "/home/ubuntu/vcc/modelo/out/bucle_techo.json"
DATOS_OUT = "/datos/modelo/out/techo"
MIN_FREE_GIB = 10.0
TIME_BUDGET_S = 4 * 3600


def free_gib(path):
    total, used, free = shutil.disk_usage(path)
    return free / (1024 ** 3)


def append_record(record):
    with open(OUT_JSONL, "a") as fh:
        fh.write(json.dumps(record) + "\n")


def run_one(line, variante, direccion, mitades, control_pool, panel, var_names,
            seed_ctrl_boot, split_meta):
    t0 = time.time()
    libre = free_gib(DATOS_OUT)
    if libre < MIN_FREE_GIB:
        raise SystemExit(
            f"PARAR: solo {libre:.1f} GiB libres en /datos (<{MIN_FREE_GIB}) antes de "
            f"{line}/{variante}/{direccion}."
        )

    real_ad, pred_ad, pair_meta = build_techo_pair(
        mitades, control_pool, panel, var_names, variante, direccion, seed_ctrl_boot
    )

    tag = f"{line}_{variante}_{direccion}"
    real_path = os.path.join(DATOS_OUT, f"techo_{tag}_real.h5ad")
    pred_path = os.path.join(DATOS_OUT, f"techo_{tag}_pred.h5ad")
    real_ad.write_h5ad(real_path)
    pred_ad.write_h5ad(pred_path)
    del real_ad, pred_ad

    outdir = os.path.join(DATOS_OUT, f"run_{tag}")
    agg_path, meta_path, params_path, run_dt = run_cell_eval2_techo(
        pred_path, real_path, outdir, pair_meta["control_source"]
    )
    # `score_vs_self` solo produce avg_score, que vale 0.0 por construccion al
    # puntuar contra si mismo. En el techo puede fallar legitimamente: cuando las
    # dos mitades concuerdan, `expr_mse_unbiased_capped_norm` (estimador
    # INSESGADO, con signo) cae a ~0 y puede salir negativo, y `score` se niega a
    # tomar como base un denominador que no sea positivo finito. Eso no invalida
    # las seis metricas, que ya estan en agg_results.csv -- asi que se registra
    # el motivo y se sigue, en vez de perder la corrida por una columna que no
    # informa nada.
    try:
        scored_path = score_vs_self(agg_path, meta_path, outdir)
        avg_score = read_avg_score(scored_path)
        score_error = None
    except Exception as e:
        avg_score = None
        score_error = str(e).splitlines()[0][:300]

    metrics = read_agg_metrics(agg_path)
    de_counts = contar_de_por_perturbacion(outdir)
    per_pert = leer_results_por_perturbacion(outdir)

    # Cruce: el jaccard recomputado desde las tablas DE debe reproducir el del
    # marcador. Si no, las columnas del punto 5 no describen lo que se puntuo.
    max_dif = 0.0
    n_cruzadas = 0
    for pert, fila in de_counts.items():
        vm = per_pert.get(pert, {}).get("de_wilcoxon_sig_jaccard")
        if vm is None:
            continue
        n_cruzadas += 1
        max_dif = max(max_dif, abs(vm - fila["jaccard_recomputado"]))

    for p in (real_path, pred_path):
        if os.path.exists(p):
            os.remove(p)
    for nombre in ("de_real.parquet", "de_pred.parquet"):
        p = os.path.join(outdir, nombre)
        if os.path.exists(p):
            os.remove(p)

    record = {
        "linea_excluida": line,
        "variante": variante,
        "direccion": direccion,
        "control_source": pair_meta["control_source"],
        "seed_split": split_meta["seed_split"],
        "seed_ctrl_boot": seed_ctrl_boot if variante == "S" else None,
        "n_control_cells": split_meta["n_control_cells"],
        "n_pred_ctrl": pair_meta["n_pred_ctrl"],
        "n_real_ctrl": pair_meta["n_real_ctrl"],
        "real_shape": pair_meta["real_shape"],
        "pred_shape": pair_meta["pred_shape"],
        "run_duration_s": round(run_dt, 2),
        "total_duration_s": round(time.time() - t0, 2),
        "free_gib_before": round(libre, 2),
        "jaccard_cruce_max_dif": max_dif,
        "jaccard_cruce_n": n_cruzadas,
        **{m: metrics.get(m) for m in VCC2026_METRICS},
        "avg_score": avg_score,
        "score_error": score_error,
        "n_pred_por_pert": pair_meta["n_pred_por_pert"],
        "n_real_por_pert": pair_meta["n_real_por_pert"],
        "de_counts": de_counts,
        "per_pert_metrics": per_pert,
    }
    return record


def main():
    os.makedirs(DATOS_OUT, exist_ok=True)
    panel_data = json.load(open(PANEL_PATH))
    panel = panel_data["panel_comun"]
    sys.stderr.write(f"panel comun: {len(panel)} perturbaciones\n")

    sys.stderr.write("leyendo obs global de IFNB (cell_type, gene)...\n")
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    sys.stderr.write(f"n_obs={len(cell_type)} n_genes={n_genes}\n")

    hechas = set()
    if os.path.exists(OUT_JSONL):
        with open(OUT_JSONL) as fh:
            for ln in fh:
                ln = ln.strip()
                if not ln:
                    continue
                r = json.loads(ln)
                if "error" not in r:
                    hechas.add((r["linea_excluida"], r["variante"], r["direccion"]))
    if hechas:
        sys.stderr.write(f"corridas ya escritas, se saltan: {len(hechas)}\n")

    t_inicio = time.time()
    n_hechas = 0
    for fold_idx, line in enumerate(LINES):
        pendientes = [
            (v, d) for v in VARIANTES for d in DIRECCIONES
            if (line, v, d) not in hechas
        ]
        if not pendientes:
            continue

        seed_split = SEED_BASE_SPLIT + fold_idx
        seed_ctrl_boot = SEED_BASE_CTRLS + fold_idx
        sys.stderr.write(
            f"\n=== pliegue {fold_idx+1}/6: {line} -- particion seed={seed_split} ===\n"
        )
        mitades, control_pool, split_meta = partir_linea(
            line, panel, cell_type, gene, n_genes, seed_split, IFNB_PATH
        )
        sys.stderr.write(
            f"{line}: control={split_meta['n_control_cells']} celulas "
            f"(A={split_meta['n_mitadA_por_pert']['NT']}, "
            f"B={split_meta['n_mitadB_por_pert']['NT']})\n"
        )

        for variante, direccion in pendientes:
            sys.stderr.write(f"  -> {line} variante {variante} direccion {direccion}... ")
            sys.stderr.flush()
            try:
                record = run_one(
                    line, variante, direccion, mitades, control_pool, panel,
                    var_names, seed_ctrl_boot, split_meta,
                )
            except Exception as e:
                sys.stderr.write(f"FALLO: {e}\n")
                append_record({
                    "linea_excluida": line, "variante": variante,
                    "direccion": direccion, "error": str(e),
                })
                raise
            record["n_celulas_por_pert"] = split_meta["n_celulas_por_pert"]
            record["n_mitadA_por_pert"] = split_meta["n_mitadA_por_pert"]
            record["n_mitadB_por_pert"] = split_meta["n_mitadB_por_pert"]
            append_record(record)
            n_hechas += 1
            sys.stderr.write(
                f"OK {record['total_duration_s']:.0f}s "
                f"jac={record['de_wilcoxon_sig_jaccard']:.4f} "
                f"pds={record['pds_cosine']:.4f} "
                f"cruce_dif={record['jaccard_cruce_max_dif']:.2e}\n"
            )

            transcurrido = time.time() - t_inicio
            if n_hechas == 1:
                total_pend = sum(
                    1 for L in LINES for v in VARIANTES for d in DIRECCIONES
                    if (L, v, d) not in hechas
                )
                proy = transcurrido * total_pend
                sys.stderr.write(
                    f"  primera corrida: {transcurrido:.0f}s -- proyeccion de "
                    f"{total_pend}: {proy/60:.1f} min\n"
                )
                if proy > TIME_BUDGET_S:
                    sys.stderr.write("PARAR: la proyeccion supera el presupuesto.\n")
                    sys.exit(3)

        del mitades, control_pool

    sys.stderr.write(f"\nbucle completo: {n_hechas} corridas nuevas en {OUT_JSONL}\n")


if __name__ == "__main__":
    main()
