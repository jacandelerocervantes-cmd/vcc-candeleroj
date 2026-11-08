#!/usr/bin/env python
"""Construye los DOS extremos de la escala vcc2026 por el camino NATIVO del
marcador, una vez por linea celular de IFNB, con la ETIQUETA DE CONTROL DE
COMPETENCIA.

Nada se calcula a mano. Este modulo solo:
  (a) arma el h5ad real de cada linea reusando `suelo_arnes.build_fold_adatas`
      SIN modificarlo -- con las MISMAS semillas del suelo, para que la
      poblacion sea exactamente la que puntuaron bucle_suelo.md y el lado real
      del techo -- y renombra la etiqueta de control de `NT` a
      `non-targeting`;
  (b) invoca `cell-eval2 baseline` (produce el `b`) y
      `cell-eval2 prep-real-bundle` (produce el `r` y empaqueta los dos
      extremos con la identidad que los liga);
  (c) captura las guardas y lee los artefactos.

POR QUE EL RENOMBRADO: `competition.py:112` declara
`_RULE_EXCLUDED = ("device", "pert_chunk", "pert_col")`. `pert_col` esta
excluido del hash de la regla, pero `control` NO, asi que `--control NT`
desviaba el `config_hash` y los seis bundles salian DIAGNOSTICO. La regla exige
`control: non-targeting`. `--pert-col target_gene` se conserva porque esta
excluido. El modo `paso0` comprueba que el renombrado no mueve ningun numero
ANTES de gastar el computo completo.

ADVERTENCIA TRANSDUCTIVA (baseline.py:6-10): el perfil del baseline se promedia
de las perturbaciones reales evaluadas, asi que es informacion de ORACULO. Este
artefacto es de PUNTUACION unicamente. Ninguna cifra suya -- ni el perfil, ni la
escala r=profile/ctrl, ni `b` -- puede entrar nunca a un modelo.
"""
import json
import os
import shutil
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (  # reusados sin modificar
    IFNB_PATH, CONTROL_LABEL, CELLEVAL2_PY, VCC2026_METRICS,
    build_fold_adatas, load_global_obs, load_var_names,
)
from suelo_floor import predict_floor

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT_JSONL = "/home/ubuntu/vcc/modelo/out/anclas_nativas_comp.json"
PASO0_JSON = "/home/ubuntu/vcc/modelo/out/anclas_paso0.json"
WORK = "/datos/modelo/out/nativas_comp"
MIN_FREE_GIB = 10.0

# La etiqueta que exige la regla de competencia (configs/vcc2026.yaml).
CONTROL_COMP = "non-targeting"

# Las MISMAS semillas del suelo, a proposito: reproducen su poblacion celda por
# celda. No son semillas nuevas.
SEED_BASE_REAL = 20260001
SEED_BASE_PRED = 20261001

# Convencion del marcador, no nuestra (competition.py:38,41).
ANCHOR_BASE_SEED = 0
ANCHOR_SPLITS = 5


def free_gib(path):
    return shutil.disk_usage(path)[2] / (1024 ** 3)


def append_record(rec, path=OUT_JSONL):
    with open(path, "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def corre(cmd, log_path):
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    with open(log_path, "w") as fh:
        fh.write("CMD: " + " ".join(cmd) + "\n\n")
        fh.write("STDOUT:\n" + proc.stdout + "\n\nSTDERR:\n" + proc.stderr + "\n")
        fh.write(f"\nEXITCODE:{proc.returncode}\nDURATION_S:{dt:.2f}\n")
    return proc, dt


def cmd_baseline(real_path, outdir, pred_path, control=None):
    """`--save-pred` NO es opcional: `prep-real-bundle --baseline` recibe el h5ad
    de la PREDICCION del baseline, no el directorio
    (`real_bundle.build_real_bundle(real, baseline_pred, ...)`, real_bundle.py:262).

    `control=None` OMITE `--control`, de modo que el valor queda el del preset
    (`non-targeting`) y el `config_hash` no se desvia."""
    cmd = [
        CELLEVAL2_PY, "baseline",
        "-ar", real_path,
        "--preset", "vcc2026", "--profile", "vcc2026",
        "--pert-col", "target_gene",
        "--emit", "dispersed", "--seed", "0",
        "--save-pred", pred_path,
        "-o", outdir,
    ]
    if control is not None:
        cmd += ["--control", control]
    return cmd


def cmd_bundle(real_path, base_pred_path, bundle_dir, control=None):
    cmd = [
        CELLEVAL2_PY, "prep-real-bundle",
        "--real", real_path,
        "--baseline", base_pred_path,
        "--preset", "vcc2026", "--profile", "vcc2026",
        "--pert-col", "target_gene",
        "--anchor-base-seed", str(ANCHOR_BASE_SEED),
        "--anchor-splits", str(ANCHOR_SPLITS),
        "-o", bundle_dir,
    ]
    if control is not None:
        cmd += ["--control", control]
    return cmd


def leer_b(dirpath, nombre="baseline_agg.csv"):
    """El `b` por metrica: la fila `mean` -- la misma que lee `score`
    (score.py:731)."""
    import pandas as pd
    df = pd.read_csv(os.path.join(dirpath, nombre))
    row = df.loc[df["statistic"] == "mean"].iloc[0]
    return {m: (float(row[m]) if m in df.columns else None) for m in VCC2026_METRICS}


def leer_r(bundle_dir):
    import polars as pl
    df = pl.read_parquet(os.path.join(bundle_dir, "anchor_agg.parquet"))
    out = {}
    for row in df.iter_rows(named=True):
        out[row["metric"]] = {
            "replicate": row["replicate"],
            "replicate_sd": row["replicate_sd"],
            "replicate_min": row["replicate_min"],
            "replicate_max": row["replicate_max"],
            "n_perturbations_min": row["n_perturbations_min"],
            "n_perturbations_max": row["n_perturbations_max"],
            "estimator": row["estimator"],
        }
    return out


def leer_splits(bundle_dir):
    import polars as pl
    df = pl.read_parquet(os.path.join(bundle_dir, "anchor_splits.parquet"))
    out = {}
    for row in df.iter_rows(named=True):
        out.setdefault(row["metric"], []).append({
            "split_index": row["split_index"],
            "seed": row["seed"],
            "value": row["value"],
            "n_perturbations": row["n_perturbations"],
        })
    for m in out:
        out[m].sort(key=lambda d: d["split_index"])
    return out


def leer_meta(path, claves):
    with open(path) as fh:
        m = json.load(fh)
    return {k: m.get(k) for k in claves}


def hash_de_la_regla(bundle_dir):
    """Recalcula `rule_config_hash` sobre el config.yaml que el bundle guardo, y
    lo compara contra el de referencia. Sirve para reportar el hash OBTENIDO si
    la guarda vuelve a saltar, en vez de solo el mensaje."""
    from cell_eval2.config import EvalConfig
    from cell_eval2 import competition
    cfg_path = os.path.join(bundle_dir, "config.yaml")
    if not os.path.exists(cfg_path):
        return {"obtenido": None, "referencia": None, "coincide": None}
    cfg = EvalConfig.from_yaml(cfg_path)
    obtenido = competition.rule_config_hash(cfg)
    referencia = competition.competition_payload()["config_hash"]
    return {
        "obtenido": obtenido,
        "referencia": referencia,
        "coincide": obtenido == referencia,
        "control_en_config": cfg.control,
        "pert_col_en_config": cfg.pert_col,
    }


def arma_real(line, fold_idx, panel, cell_type, gene, var_names, n_genes,
              control_label):
    """Arma el h5ad real de la linea reusando build_fold_adatas SIN modificarlo,
    y renombra la etiqueta de control si se pide.

    build_fold_adatas exige una funcion de emision por contrato; se le pasa
    predict_floor y se DESCARTA su salida. Lo que se usa es real_ad, que asi sale
    identico celda por celda al que puntuo bucle_suelo.md."""
    seed_real = SEED_BASE_REAL + fold_idx
    seed_pred = SEED_BASE_PRED + fold_idx
    real_ad, descartado, build_meta = build_fold_adatas(
        line, panel, predict_floor, var_names, seed_real, seed_pred,
        cell_type, gene, n_genes, IFNB_PATH,
    )
    del descartado
    n_antes = int((real_ad.obs["target_gene"] == CONTROL_LABEL).sum())
    if control_label != CONTROL_LABEL:
        real_ad.obs["target_gene"] = real_ad.obs["target_gene"].replace(
            {CONTROL_LABEL: control_label}
        )
    n_despues = int((real_ad.obs["target_gene"] == control_label).sum())
    assert n_antes == n_despues, (
        f"el renombrado movio filas: {n_antes} -> {n_despues}"
    )
    build_meta["n_filas_control_renombradas"] = n_despues
    build_meta["etiqueta_control"] = control_label
    build_meta["seed_real"] = seed_real
    return real_ad, build_meta


# --------------------------------------------------------------------------
# PASO 0: la invariancia del renombrado, solo A549, antes de gastar el computo
# --------------------------------------------------------------------------
def paso0(panel, cell_type, gene, var_names, n_genes):
    line = "A549"
    fold_idx = 0
    os.makedirs(WORK, exist_ok=True)
    sys.stderr.write(f"\n=== PASO 0: invariancia del renombrado, {line} ===\n")

    brazos = [
        # (nombre, etiqueta en el obs, valor para --control)
        ("NT_con_flag", CONTROL_LABEL, CONTROL_LABEL),
        ("nontargeting_sin_flag", CONTROL_COMP, None),
    ]
    res = {}
    for nombre, etiqueta, flag in brazos:
        real_ad, bmeta = arma_real(line, fold_idx, panel, cell_type, gene,
                                   var_names, n_genes, etiqueta)
        real_path = os.path.join(WORK, f"paso0_{nombre}_real.h5ad")
        pred_path = os.path.join(WORK, f"paso0_{nombre}_pred.h5ad")
        outdir = os.path.join(WORK, f"PASO0_{nombre}")
        real_ad.write_h5ad(real_path)
        shape = list(real_ad.shape)
        del real_ad
        os.makedirs(outdir, exist_ok=True)
        sys.stderr.write(f"  [{nombre}] cell-eval2 baseline... ")
        sys.stderr.flush()
        p, dt = corre(cmd_baseline(real_path, outdir, pred_path, flag),
                      os.path.join(WORK, f"paso0_{nombre}.log"))
        if p.returncode != 0:
            raise RuntimeError(
                f"paso0/{nombre} fallo (exit={p.returncode}). Ver "
                f"{WORK}/paso0_{nombre}.log"
            )
        sys.stderr.write(f"OK {dt:.0f}s\n")
        res[nombre] = {
            "etiqueta_control": etiqueta,
            "flag_control": flag,
            "real_shape": shape,
            "n_filas_control": bmeta["n_filas_control_renombradas"],
            "duracion_s": round(dt, 2),
            "b": leer_b(outdir),
            "meta": leer_meta(
                os.path.join(outdir, "baseline_meta.json"),
                ["source_fingerprint_strict", "source_fingerprint",
                 "config_digest", "comparator"],
            ),
        }
        for pth in (real_path, pred_path):
            if os.path.exists(pth):
                os.remove(pth)

    a, b = res["NT_con_flag"]["b"], res["nontargeting_sin_flag"]["b"]
    difs = {}
    identicas = True
    for m in VCC2026_METRICS:
        va, vb = a.get(m), b.get(m)
        if va is None or vb is None:
            difs[m] = None
            identicas = False
            continue
        d = abs(va - vb)
        difs[m] = d
        # BIT A BIT: igualdad exacta de float, no una tolerancia.
        if va != vb:
            identicas = False

    salida = {
        "linea": line,
        "brazos": res,
        "difs_absolutas": difs,
        "identicas_bit_a_bit": identicas,
        "n_filas_control_iguales": (
            res["NT_con_flag"]["n_filas_control"]
            == res["nontargeting_sin_flag"]["n_filas_control"]
        ),
    }
    with open(PASO0_JSON, "w") as fh:
        json.dump(salida, fh, indent=2, sort_keys=True)

    sys.stderr.write("\n  metrica                                   NT        non-targeting   dif\n")
    for m in VCC2026_METRICS:
        sys.stderr.write(
            f"  {m:40s} {a.get(m)!r:>12} {b.get(m)!r:>14} "
            f"{difs[m]!r}\n"
        )
    sys.stderr.write(
        f"\n  IDENTICAS BIT A BIT: {identicas}\n"
        f"  escrito {PASO0_JSON}\n"
    )
    return identicas


# --------------------------------------------------------------------------
# PASO 1: las 6 lineas con la etiqueta de competencia
# --------------------------------------------------------------------------
def una_linea(fold_idx, line, panel, cell_type, gene, var_names, n_genes):
    t0 = time.time()
    libre = free_gib(WORK)
    if libre < MIN_FREE_GIB:
        raise SystemExit(f"PARAR: {libre:.1f} GiB libres en /datos (<{MIN_FREE_GIB}).")

    real_ad, build_meta = arma_real(line, fold_idx, panel, cell_type, gene,
                                    var_names, n_genes, CONTROL_COMP)
    real_path = os.path.join(WORK, f"real_{line}.h5ad")
    real_ad.write_h5ad(real_path)
    real_shape = list(real_ad.shape)
    del real_ad

    base_dir = os.path.join(WORK, f"BASE_{line}")
    base_pred_path = os.path.join(WORK, f"basepred_{line}.h5ad")
    bundle_dir = os.path.join(WORK, f"BUNDLE_{line}")
    os.makedirs(base_dir, exist_ok=True)

    sys.stderr.write(f"  [{line}] cell-eval2 baseline (control del preset)... ")
    sys.stderr.flush()
    p_base, dt_base = corre(cmd_baseline(real_path, base_dir, base_pred_path, None),
                            os.path.join(WORK, f"baseline_{line}.log"))
    if p_base.returncode != 0:
        raise RuntimeError(
            f"baseline fallo para {line} (exit={p_base.returncode}). Ver "
            f"{WORK}/baseline_{line}.log"
        )
    sys.stderr.write(f"OK {dt_base:.0f}s\n")

    base_meta = leer_meta(
        os.path.join(base_dir, "baseline_meta.json"),
        ["source_fingerprint_strict", "source_fingerprint", "cell_eval2_version",
         "comparator", "config_digest"],
    )

    sys.stderr.write(f"  [{line}] cell-eval2 prep-real-bundle (k={ANCHOR_SPLITS})... ")
    sys.stderr.flush()
    p_bun, dt_bun = corre(cmd_bundle(real_path, base_pred_path, bundle_dir, None),
                          os.path.join(WORK, f"bundle_{line}.log"))
    if p_bun.returncode != 0:
        raise RuntimeError(
            f"prep-real-bundle fallo para {line} (exit={p_bun.returncode}). Ver "
            f"{WORK}/bundle_{line}.log"
        )
    sys.stderr.write(f"OK {dt_bun:.0f}s\n")

    # GUARDA (cli.py:788-793): la linea literal que distingue competencia de
    # DIAGNOSTICO. Se captura tal cual, no se interpreta.
    lineas = p_bun.stdout.strip().splitlines()
    linea_guarda = lineas[-1].strip() if lineas else ""
    es_diagnostico = "DIAGNOSTIC" in linea_guarda

    manifest = leer_meta(
        os.path.join(bundle_dir, "manifest.json"),
        ["manifest_version", "real_bundle_id", "rule_digest", "rule_mismatches",
         "source_fingerprint_strict", "source_fingerprint", "cell_eval2_version",
         "base_seed", "n_splits", "derived_seeds", "seed_derivation",
         "estimators", "members", "profile", "control_source_effective",
         "anchor_config_hash", "anchor_digest", "bulk_target_sum",
         "resolved_de_backend", "resolved_device", "real_fingerprint"],
    )

    rec = {
        "fold_idx": fold_idx,
        "linea": line,
        "etiqueta_control": CONTROL_COMP,
        "flag_control_pasado": None,
        "seed_real": build_meta["seed_real"],
        "real_shape": real_shape,
        "n_control_cells": build_meta["n_control_cells"],
        "n_filas_control": build_meta["n_filas_control_renombradas"],
        "per_pert_pool_sizes": build_meta["per_pert_pool_sizes"],
        "anchor_base_seed": ANCHOR_BASE_SEED,
        "anchor_splits": ANCHOR_SPLITS,
        "baseline_duration_s": round(dt_base, 2),
        "bundle_duration_s": round(dt_bun, 2),
        "total_duration_s": round(time.time() - t0, 2),
        "free_gib_before": round(libre, 2),
        "guarda_linea_literal": linea_guarda,
        "bundle_diagnostico": es_diagnostico,
        "hash_regla": hash_de_la_regla(bundle_dir),
        "baseline_meta": base_meta,
        "manifest": manifest,
        "b": leer_b(bundle_dir),
        "b_corrida_suelta": leer_b(base_dir),
    }

    if es_diagnostico:
        rec["r"] = None
        rec["splits"] = None
        rec["detenido_por_guarda"] = True
    else:
        rec["r"] = leer_r(bundle_dir)
        rec["splits"] = leer_splits(bundle_dir)
        rec["detenido_por_guarda"] = False

    for p in (real_path, base_pred_path):
        if os.path.exists(p):
            os.remove(p)
    return rec


def main():
    modo = sys.argv[1] if len(sys.argv) > 1 else "paso1"
    os.makedirs(WORK, exist_ok=True)
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    sys.stderr.write(f"panel comun: {len(panel)} perturbaciones\n")
    sys.stderr.write("leyendo obs global de IFNB...\n")
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    sys.stderr.write(f"n_obs={len(cell_type)} n_genes={n_genes}\n")

    if modo == "paso0":
        ok = paso0(panel, cell_type, gene, var_names, n_genes)
        if not ok:
            sys.stderr.write(
                "\n*** PARAR: las seis metricas NO son identicas bit a bit. El "
                "renombrado no es un renombrado, es un cambio. No se corre el "
                "paso 1.\n"
            )
            sys.exit(4)
        sys.stderr.write(
            "\n*** OK: el renombrado no mueve ningun numero. Todo lo ya medido "
            "sigue comparable. Se puede correr el paso 1.\n"
        )
        return

    hechas = set()
    if os.path.exists(OUT_JSONL):
        with open(OUT_JSONL) as fh:
            for ln in fh:
                ln = ln.strip()
                if ln:
                    r = json.loads(ln)
                    if "error" not in r:
                        hechas.add(r["linea"])
    if hechas:
        sys.stderr.write(f"lineas ya hechas, se saltan: {sorted(hechas)}\n")

    for fold_idx, line in enumerate(LINES):
        if line in hechas:
            continue
        sys.stderr.write(
            f"\n=== {fold_idx+1}/6: {line} (seed_real={SEED_BASE_REAL+fold_idx}, "
            f"control={CONTROL_COMP!r}) ===\n"
        )
        try:
            rec = una_linea(fold_idx, line, panel, cell_type, gene, var_names, n_genes)
        except Exception as e:
            sys.stderr.write(f"FALLO en {line}: {e}\n")
            append_record({"fold_idx": fold_idx, "linea": line, "error": str(e)})
            raise
        append_record(rec)
        marca = "DIAGNOSTICO" if rec["bundle_diagnostico"] else "competencia"
        sys.stderr.write(
            f"  [{line}] bundle {marca} -- {rec['total_duration_s']:.0f}s\n"
            f"  [{line}] guarda: {rec['guarda_linea_literal']}\n"
        )
        if rec["bundle_diagnostico"]:
            hr = rec["hash_regla"]
            sys.stderr.write(
                f"  [{line}] DETENIDO POR GUARDA: no se leyo el anchor.\n"
                f"  [{line}] motivo estampado: {rec['manifest'].get('rule_mismatches')}\n"
                f"  [{line}] hash obtenido  : {hr.get('obtenido')}\n"
                f"  [{line}] hash referencia: {hr.get('referencia')}\n"
                f"  [{line}] control en config: {hr.get('control_en_config')!r}\n"
            )
            sys.stderr.write(
                "\n*** PARAR: un bundle salio DIAGNOSTICO. No se siguen las "
                "lineas restantes ni se lee ningun anchor.\n"
            )
            sys.exit(5)

    sys.stderr.write(f"\nlisto: {OUT_JSONL}\n")


if __name__ == "__main__":
    main()
