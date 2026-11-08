#!/usr/bin/env python
"""Prueba del cable: estrena `run` -> `score --real-bundle` con una emision de la
que ya sabemos la respuesta. NO es un candidato: la emision es desechable.

Emision: `suelo_floor.predict_floor` SIN modificar -- el control remuestreado con
reemplazo, sin escalar. El `b` oficial es el control ESCALADO al perfil de
respuesta promedio (`baseline.py:344`), asi que esta emision deberia quedar POR
DEBAJO de `b` en casi todo y su `s` salir ~0 o negativo. Si sale muy por encima
de 0, el problema esta en la cadena, no en el suelo.

Poblacion: la MISMA de siempre -- `build_fold_adatas` sin modificar,
`seed_real = 20260001 + fold_idx`, etiqueta de control `non-targeting` (la que
exige la regla de competencia; ver anclas_nativas_v2.md paso 0 para la prueba de
que el renombrado no mueve ningun numero).

No repara nada: si `score` se niega, se registra el mensaje literal y se para.
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
from suelo_floor import predict_floor  # la emision, sin modificar

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT_JSONL = "/home/ubuntu/vcc/modelo/out/cable_suelo.json"
WORK = "/datos/modelo/out/cable"
BUNDLES = "/datos/modelo/out/nativas_comp"
MIN_FREE_GIB = 10.0

CONTROL_COMP = "non-targeting"
SEED_BASE_REAL = 20260001
SEED_BASE_PRED = 20261001


def free_gib(path):
    return shutil.disk_usage(path)[2] / (1024 ** 3)


def append_record(rec):
    with open(OUT_JSONL, "a") as fh:
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


def cmd_run(pred_path, real_path, outdir):
    """Sin `--control`: el valor queda el del preset (`non-targeting`) y el
    config_hash no se desvia. `--pert-col` SI se pasa porque esta en
    `_RULE_EXCLUDED` (competition.py:112)."""
    return [
        CELLEVAL2_PY, "run",
        "-ap", pred_path, "-ar", real_path,
        "--preset", "vcc2026", "--profile", "vcc2026",
        "--pert-col", "target_gene",
        "-o", outdir,
    ]


def cmd_score(run_dir, bundle_dir, out_csv):
    return [
        CELLEVAL2_PY, "score",
        "--user-agg", os.path.join(run_dir, "agg_results.csv"),
        "--user-meta", os.path.join(run_dir, "run_meta.json"),
        "--real-bundle", bundle_dir,
        "-o", out_csv,
    ]


def hash_del_run(run_dir):
    """Recalcula `rule_config_hash` sobre el run_params.yaml que escribio `run`,
    y lo compara contra el de referencia. `run` no imprime una linea de
    competencia como `prep-real-bundle`, asi que esta es la unica forma de saber
    si su config califica."""
    from cell_eval2.config import EvalConfig
    from cell_eval2 import competition
    p = os.path.join(run_dir, "run_params.yaml")
    if not os.path.exists(p):
        return {"obtenido": None, "referencia": None, "coincide": None}
    cfg = EvalConfig.from_yaml(p)
    obtenido = competition.rule_config_hash(cfg)
    referencia = competition.competition_payload()["config_hash"]
    return {
        "obtenido": obtenido,
        "referencia": referencia,
        "coincide": obtenido == referencia,
        "control_en_config": cfg.control,
        "pert_col_en_config": cfg.pert_col,
        "cache_strict": cfg.cache_strict,
    }


def leer_scored(path):
    """scored.csv es largo. Con un bundle de competencia deberia traer
    `from_replicate` ademas de `from_baseline`, y el `--help` de `score` dice que
    en ese caso `from_replicate` toma `avg_score` y el de `from_baseline` queda
    null. Se lee TODO sin suponer que columnas hay."""
    import pandas as pd
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    cols = [c for c in df.columns if c != "metric"]
    filas = {}
    for _, row in df.iterrows():
        filas[str(row["metric"])] = {c: (None if pd.isna(row[c]) else
                                        (float(row[c]) if isinstance(row[c], (int, float))
                                         else row[c]))
                                     for c in cols}
    return {"columnas": cols, "filas": filas}


def renombra(adata):
    n_antes = int((adata.obs["target_gene"] == CONTROL_LABEL).sum())
    adata.obs["target_gene"] = adata.obs["target_gene"].replace(
        {CONTROL_LABEL: CONTROL_COMP}
    )
    n_despues = int((adata.obs["target_gene"] == CONTROL_COMP).sum())
    assert n_antes == n_despues, f"el renombrado movio filas: {n_antes} -> {n_despues}"
    return n_despues


def una_linea(fold_idx, line, panel, cell_type, gene, var_names, n_genes):
    t0 = time.time()
    libre = free_gib(WORK)
    if libre < MIN_FREE_GIB:
        raise SystemExit(f"PARAR: {libre:.1f} GiB libres en /datos (<{MIN_FREE_GIB}).")

    bundle_dir = os.path.join(BUNDLES, f"BUNDLE_{line}")
    if not os.path.isdir(bundle_dir):
        raise SystemExit(f"PARAR: no existe el bundle {bundle_dir}")

    seed_real = SEED_BASE_REAL + fold_idx
    seed_pred = SEED_BASE_PRED + fold_idx

    # build_fold_adatas sin modificar: devuelve la VERDAD y la emision de
    # predict_floor en un solo paso, con las semillas del suelo.
    real_ad, pred_ad, build_meta = build_fold_adatas(
        line, panel, predict_floor, var_names, seed_real, seed_pred,
        cell_type, gene, n_genes, IFNB_PATH,
    )
    # La etiqueta se renombra en LOS DOS lados: `io.validate_pair` exige que los
    # conjuntos de etiquetas coincidan entre PRED y REAL.
    n_ctrl_real = renombra(real_ad)
    n_ctrl_pred = renombra(pred_ad)

    real_path = os.path.join(WORK, f"real_{line}.h5ad")
    pred_path = os.path.join(WORK, f"pred_{line}.h5ad")
    real_ad.write_h5ad(real_path)
    pred_ad.write_h5ad(pred_path)
    real_shape, pred_shape = list(real_ad.shape), list(pred_ad.shape)
    del real_ad, pred_ad

    run_dir = os.path.join(WORK, f"run_{line}")
    sys.stderr.write(f"  [{line}] cell-eval2 run... ")
    sys.stderr.flush()
    p_run, dt_run = corre(cmd_run(pred_path, real_path, run_dir),
                          os.path.join(WORK, f"run_{line}.log"))
    run_ok = p_run.returncode == 0
    sys.stderr.write(f"{'OK' if run_ok else 'FALLO'} {dt_run:.0f}s\n")

    rec = {
        "fold_idx": fold_idx,
        "linea": line,
        "emision": "suelo_floor.predict_floor (control remuestreado, sin escalar)",
        "seed_real": seed_real,
        "seed_pred": seed_pred,
        "etiqueta_control": CONTROL_COMP,
        "real_shape": real_shape,
        "pred_shape": pred_shape,
        "n_control_real": n_ctrl_real,
        "n_control_pred": n_ctrl_pred,
        "n_control_cells_pool": build_meta["n_control_cells"],
        "run_returncode": p_run.returncode,
        "run_duration_s": round(dt_run, 2),
        "run_stderr_ultimas": p_run.stderr.strip().splitlines()[-6:] if p_run.stderr else [],
        "free_gib_before": round(libre, 2),
    }

    if not run_ok:
        rec["fallo_en"] = "run"
        rec["mensaje_literal"] = p_run.stderr.strip()[-2000:]
        for p in (real_path, pred_path):
            if os.path.exists(p):
                os.remove(p)
        return rec

    rec["hash_run"] = hash_del_run(run_dir)
    rec["run_meta"] = {}
    rm_path = os.path.join(run_dir, "run_meta.json")
    if os.path.exists(rm_path):
        with open(rm_path) as fh:
            rm = json.load(fh)
        rec["run_meta"] = {
            k: rm.get(k) for k in
            ("source_fingerprint_strict", "source_fingerprint", "cell_eval2_version",
             "config_digest", "comparator", "anchor_semantic_identity",
             "anchor_metric_names", "de_pred_fingerprint")
        }

    # Las metricas crudas del run, para poder leer el `s` con su `u` al lado.
    import pandas as pd
    agg = os.path.join(run_dir, "agg_results.csv")
    if os.path.exists(agg):
        df = pd.read_csv(agg)
        row = df.loc[df["statistic"] == "mean"].iloc[0]
        rec["u"] = {m: (float(row[m]) if m in df.columns else None)
                    for m in VCC2026_METRICS}

    out_csv = os.path.join(WORK, f"scored_{line}.csv")
    sys.stderr.write(f"  [{line}] cell-eval2 score --real-bundle... ")
    sys.stderr.flush()
    p_sc, dt_sc = corre(cmd_score(run_dir, bundle_dir, out_csv),
                        os.path.join(WORK, f"score_{line}.log"))
    score_ok = p_sc.returncode == 0
    sys.stderr.write(f"{'OK' if score_ok else 'SE NEGO'} {dt_sc:.0f}s\n")

    rec["score_returncode"] = p_sc.returncode
    rec["score_duration_s"] = round(dt_sc, 2)
    rec["score_stdout"] = p_sc.stdout.strip()[-2000:]
    rec["score_stderr_literal"] = p_sc.stderr.strip()[-2000:]
    rec["score_ok"] = score_ok
    rec["scored"] = leer_scored(out_csv) if score_ok else None
    if not score_ok:
        rec["fallo_en"] = "score"
        rec["mensaje_literal"] = p_sc.stderr.strip()[-2000:]

    for p in (real_path, pred_path):
        if os.path.exists(p):
            os.remove(p)
    rec["total_duration_s"] = round(time.time() - t0, 2)
    return rec


def main():
    os.makedirs(WORK, exist_ok=True)
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    sys.stderr.write(f"panel comun: {len(panel)} perturbaciones\n")
    sys.stderr.write("leyendo obs global de IFNB...\n")
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    sys.stderr.write(f"n_obs={len(cell_type)} n_genes={n_genes}\n")

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
        sys.stderr.write(f"\n=== {fold_idx+1}/6: {line} ===\n")
        try:
            rec = una_linea(fold_idx, line, panel, cell_type, gene, var_names, n_genes)
        except Exception as e:
            sys.stderr.write(f"EXCEPCION en {line}: {e}\n")
            append_record({"fold_idx": fold_idx, "linea": line, "error": str(e)})
            raise
        append_record(rec)

        if rec.get("fallo_en"):
            sys.stderr.write(
                f"\n*** PARAR: fallo en `{rec['fallo_en']}` para {line}. No se "
                "repara nada y no se siguen las lineas restantes.\n"
                f"*** mensaje literal:\n{rec['mensaje_literal']}\n"
            )
            sys.exit(6)

        sc = rec.get("scored") or {}
        avg = (sc.get("filas") or {}).get("avg_score", {})
        sys.stderr.write(
            f"  [{line}] columnas de scored.csv: {sc.get('columnas')}\n"
            f"  [{line}] avg_score: {avg}\n"
            f"  [{line}] config_hash de competencia: "
            f"{(rec.get('hash_run') or {}).get('coincide')}\n"
        )

    sys.stderr.write(f"\nlisto: {OUT_JSONL}\n")


if __name__ == "__main__":
    main()
