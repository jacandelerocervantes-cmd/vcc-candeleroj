#!/usr/bin/env python
"""Huecos 2 y 3: la cohorte del lado real, y el piso de azar de `reach`.

QUE MIDE, y por que hace falta correr algo en vez de leer artefactos:

  * HUECO 2 -- la cohorte. `results.csv` y `baseline_results.csv` dicen el VALOR
    por perturbacion, no los conteos DE que lo producen. Para saber si la
    cohorte "la referencia tiene senal DE" es la misma entre objetos hay que ver
    el lado real de cada uno, y eso solo sale de `de_real.parquet`
    (`--write-degenes`). Se corre el MISMO par (real, prediccion) de cada objeto
    sobre la MISMA poblacion de siempre.

  * HUECO 3 -- `n_conf`. Es el denominador de `reach`
    (`metrics/direction.py:1054`, `k_star/n_conf`) y lo construye
    `_reference_stats` desde `prepared.real_df` SOLO (direction.py:293-301,
    343-350). No se persiste en ningun artefacto. Se obtiene llamando a
    `_reference_stats` -- la funcion del marcador, sin reimplementarla -- sobre
    el `prepared` de las tablas DE que escribe `--write-degenes`.

LOS TRES OBJETOS, por linea:
  * nuestro suelo: `predict_floor` sin modificar, semillas de siempre
  * el `b` oficial: se REGENERA su h5ad de prediccion con `cell-eval2 baseline
    --seed 0 --save-pred` (el original se borro al terminar `anclas_nativas`) y
    se verifica que `baseline_results.csv` reproduce el de `BASE_<linea>` antes
    de usarlo para nada
  * cada particion del anchor: la cubre `descomp_rgen.py`, que ya guarda
    `n_DE_real` por perturbacion y particion

CRITERIO DE SELECCION, declarado antes de medir: las SEIS lineas, los 33 de
`suelo_panel.json`, sin excluir ninguna celda por su valor ni por su signo.

No modifica ninguna metrica ni ningun arnes: `build_fold_adatas`,
`predict_floor`, `contar_de_por_perturbacion` y `_reference_stats` se usan tal
cual.
"""
import json
import os
import subprocess
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (
    IFNB_PATH, CONTROL_LABEL, CELLEVAL2_PY,
    build_fold_adatas, load_global_obs, load_var_names,
)
from suelo_floor import predict_floor
from techo_arnes import contar_de_por_perturbacion

FOLD_IDX = {"A549": 0, "BXPC3": 1, "HAP1": 2, "HT29": 3, "K562": 4, "MCF7": 5}
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
BASE_OFICIAL = "/datos/modelo/out/nativas_comp/BASE_{L}/baseline_results.csv"
AGG_OFICIAL = "/datos/modelo/out/nativas_comp/BASE_{L}/baseline_agg.csv"
OUT_JSON = "/home/ubuntu/vcc/modelo/out/huecos_cohorte.json"
WORK = "/datos/modelo/out/huecos"
CONTROL_COMP = "non-targeting"
SEED_BASE_REAL = 20260001
SEED_BASE_PRED = 20261001

REACH = "de_wilcoxon_direction_reach_raw"
JAC = "de_wilcoxon_sig_jaccard"
PDS = "pds_cosine"
FID = "de_wilcoxon_direction_fidelity_yield_raw"
LFC = "de_wilcoxon_lfc_nmae"


def corre(cmd, log):
    t0 = time.time()
    p = subprocess.run(cmd, capture_output=True, text=True)
    with open(log, "w") as fh:
        fh.write("CMD: " + " ".join(cmd) + "\n\nSTDOUT:\n" + p.stdout
                 + "\n\nSTDERR:\n" + p.stderr
                 + f"\nEXIT:{p.returncode}\nDUR:{time.time()-t0:.1f}\n")
    if p.returncode != 0:
        raise RuntimeError(f"fallo: {' '.join(cmd[:4])}... ver {log}")
    return p


def refstats(outdir):
    """`n_conf`, `q` y `d` por perturbacion, con la funcion del marcador.

    `prepare_de` se invoca con los MISMOS parametros que
    `contar_de_por_perturbacion` (los del preset vcc2026), y `_reference_stats`
    se llama tal cual: no se reimplementa el conteo.
    """
    from cell_eval2.de import prepare_de
    from cell_eval2.metrics.direction import _reference_stats
    from techo_arnes import (DE_SORT_BY, DE_P_ADJ_THRESHOLD,
                             DE_NAN_LFC_POLICY, DE_MIN_ABS_LOG2FC)

    prepared = prepare_de(
        os.path.join(outdir, "de_pred.parquet"),
        os.path.join(outdir, "de_real.parquet"),
        control=CONTROL_COMP, sort_by=DE_SORT_BY,
        p_adj_threshold=DE_P_ADJ_THRESHOLD, nan_lfc_policy=DE_NAN_LFC_POLICY,
        min_abs_log2fc=DE_MIN_ABS_LOG2FC,
    )
    st = _reference_stats(prepared).to_dicts()
    return {
        str(r["target"]): {
            "n_conf": int(r["n_conf"]),
            "q": None if r["q"] is None else float(r["q"]),
            "d": None if r["d"] is None else float(r["d"]),
        }
        for r in st
    }


def por_pert(path):
    """`results.csv` / `baseline_results.csv` -> dict[pert][metric] = valor.
    NaN se conserva como None para poder distinguir "no entro" de 0.0."""
    df = pd.read_csv(path)
    out = {}
    for pert, metric, value in zip(df["perturbation"], df["metric"], df["value"]):
        v = float(value)
        out.setdefault(str(pert), {})[str(metric)] = None if np.isnan(v) else v
    return out


def un_objeto(nombre, real_path, pred_path, linea, extra_set=None):
    outdir = os.path.join(WORK, f"{linea}_{nombre}_run")
    cmd = [
        CELLEVAL2_PY, "run",
        "-ap", pred_path, "-ar", real_path,
        "--preset", "vcc2026", "--profile", "vcc2026",
        "--pert-col", "target_gene",
        "--write-degenes",
        "-o", outdir,
    ]
    if extra_set:
        cmd += ["--set", extra_set]
    t0 = time.time()
    sys.stderr.write(f"  {linea}/{nombre}: cell-eval2 run --write-degenes... ")
    sys.stderr.flush()
    corre(cmd, os.path.join(WORK, f"{linea}_{nombre}_run.log"))
    cuentas = contar_de_por_perturbacion(outdir, control=CONTROL_COMP)
    stats = refstats(outdir)
    valores = por_pert(os.path.join(outdir, "results.csv"))
    agg = pd.read_csv(os.path.join(outdir, "agg_results.csv"))
    fila = agg.loc[agg["statistic"] == "mean"].iloc[0]
    sys.stderr.write(f"ok ({time.time()-t0:.0f}s)\n")
    for q in ("de_real.parquet", "de_pred.parquet"):
        p = os.path.join(outdir, q)
        if os.path.exists(p):
            os.remove(p)
    return {
        "outdir": outdir,
        "agg_mean": {m: (None if pd.isna(fila[m]) else float(fila[m]))
                     for m in (PDS, "expr_mse_unbiased_capped_norm", LFC, FID, REACH, JAC)
                     if m in agg.columns},
        "por_pert": {
            p: {
                "n_DE_real": cuentas[p]["n_DE_real"],
                "n_DE_pred": cuentas[p]["n_DE_pred"],
                "n_DE_union": cuentas[p]["n_DE_union"],
                "n_conf": stats.get(p, {}).get("n_conf"),
                "q": stats.get(p, {}).get("q"),
                "jac": valores.get(p, {}).get(JAC),
                "reach": valores.get(p, {}).get(REACH),
                "fid": valores.get(p, {}).get(FID),
                "lfc_nmae": valores.get(p, {}).get(LFC),
                "pds": valores.get(p, {}).get(PDS),
            }
            for p in sorted(cuentas)
        },
    }


def main():
    os.makedirs(WORK, exist_ok=True)
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)

    lineas = sys.argv[1:] or list(FOLD_IDX)
    salida = json.load(open(OUT_JSON)) if os.path.exists(OUT_JSON) else {}

    for linea in lineas:
        if linea in salida:
            sys.stderr.write(f"{linea}: ya hecho\n")
            continue
        fold_idx = FOLD_IDX[linea]
        sys.stderr.write(f"\n=== {linea} (fold {fold_idx}) ===\n")
        real_ad, pred_ad, meta = build_fold_adatas(
            linea, panel, predict_floor, var_names,
            SEED_BASE_REAL + fold_idx, SEED_BASE_PRED + fold_idx,
            cell_type, gene, n_genes, IFNB_PATH,
        )
        for a in (real_ad, pred_ad):
            a.obs["target_gene"] = a.obs["target_gene"].replace(
                {CONTROL_LABEL: CONTROL_COMP}
            )
        real_path = os.path.join(WORK, f"{linea}_real.h5ad")
        suelo_path = os.path.join(WORK, f"{linea}_suelo.h5ad")
        real_ad.write_h5ad(real_path)
        pred_ad.write_h5ad(suelo_path)
        del real_ad, pred_ad

        # --- el `b` oficial: se regenera su prediccion y se VERIFICA ---
        base_dir = os.path.join(WORK, f"{linea}_BASE2")
        base_pred = os.path.join(WORK, f"{linea}_basepred.h5ad")
        sys.stderr.write(f"  {linea}: cell-eval2 baseline --seed 0 --save-pred... ")
        sys.stderr.flush()
        t0 = time.time()
        corre([
            CELLEVAL2_PY, "baseline",
            "-ar", real_path,
            "--preset", "vcc2026", "--profile", "vcc2026",
            "--pert-col", "target_gene",
            "--emit", "dispersed", "--seed", "0",
            "--save-pred", base_pred,
            "-o", base_dir,
        ], os.path.join(WORK, f"{linea}_baseline.log"))
        sys.stderr.write(f"ok ({time.time()-t0:.0f}s)\n")

        of = por_pert(BASE_OFICIAL.format(L=linea))
        nu = por_pert(os.path.join(base_dir, "baseline_results.csv"))
        difs = []
        for p in sorted(set(of) | set(nu)):
            for m in sorted(set(of.get(p, {})) | set(nu.get(p, {}))):
                a, b = of.get(p, {}).get(m), nu.get(p, {}).get(m)
                if a is None and b is None:
                    continue
                if a is None or b is None or a != b:
                    difs.append({"pert": p, "metric": m, "oficial": a, "regen": b})
        verif = {
            "n_celdas_comparadas": sum(len(v) for v in of.values()),
            "n_difieren": len(difs),
            "identico": len(difs) == 0,
            "difs": difs[:40],
        }
        sys.stderr.write(
            f"  {linea}: baseline regenerado vs oficial -> "
            f"{'IDENTICO' if verif['identico'] else str(len(difs)) + ' DIFIEREN'}\n"
        )

        objetos = {
            "suelo": un_objeto("suelo", real_path, suelo_path, linea),
            # `allow_fractional_counts=true` NO es una desviacion nuestra: es lo
            # que `cell-eval2 baseline` se pone a si mismo (`baseline.py:1320-1321`,
            # "the profile is a MEAN, hence fractional in any space") y esta en
            # `DIGEST_EXEMPT_FIELDS` (`baseline.py:692-703`: "a validation
            # ALLOWANCE, not a scoring semantic"). Sin el, `run` rechaza la
            # emision del baseline con "declared input_type='counts' but values
            # are fractional" (`norm.py:482`). Aplica al lado PRED solamente
            # (`run.py:397-401`), asi que el lado real sigue validado estricto.
            "b_oficial": un_objeto("b", real_path, base_pred, linea,
                                   extra_set="allow_fractional_counts=true"),
        }

        salida[linea] = {
            "fold_idx": fold_idx,
            "seed_real": SEED_BASE_REAL + fold_idx,
            "seed_pred": SEED_BASE_PRED + fold_idx,
            "real_shape": meta["real_shape"],
            "n_control_cells": meta["n_control_cells"],
            "verif_baseline_regenerado": verif,
            "objetos": objetos,
        }
        with open(OUT_JSON, "w") as fh:
            json.dump(salida, fh, indent=2, sort_keys=True)

        for p in (real_path, suelo_path, base_pred):
            if os.path.exists(p):
                os.remove(p)

    sys.stderr.write(f"\nlisto: {OUT_JSON}\n")


if __name__ == "__main__":
    main()
