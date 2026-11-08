#!/usr/bin/env python
"""Medicion del brazo PERMUTADO: 6 lineas x 8 desarreglos, puntuadas contra los
6 bundles de competencia. Mismo arnes, misma escala, misma poblacion real que
`gwps_via.py`; lo unico que cambia es que cada etiqueta recibe el perfil de OTRO
blanco (`gwps_permutado.py`).
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (
    IFNB_PATH, CONTROL_LABEL, build_fold_adatas, load_global_obs,
    load_var_names, VCC2026_METRICS,
)
from gwps_permutado import hace_predict_fn_perm
import gwps_via as GV

LINES = GV.LINES
PANEL_PATH = GV.PANEL_PATH
BUNDLES = GV.BUNDLES
OUT_JSONL = os.environ.get("PERM_OUT",
                           "/home/ubuntu/vcc/modelo/out/control_permutado.json")
WORK = os.environ.get("PERM_WORK", "/datos/modelo/out/gwps_perm")
SEED_BASE_REAL = 20260001
SEED_BASE_EMISION = 20272000          # familia nueva
N_SEMILLAS = int(os.environ.get("PERM_N", "8"))


def append_record(rec):
    with open(OUT_JSONL, "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def una_corrida(line, fold_idx, j, panel, cell_type, gene, var_names, n_genes,
                real_path, ref):
    t0 = time.time()
    if GV.free_gib(WORK) < GV.MIN_FREE_GIB:
        raise SystemExit("PARAR: poco espacio en /datos.")
    seed_real = SEED_BASE_REAL + fold_idx
    seed_em = SEED_BASE_EMISION + 100 * fold_idx + j
    registro = {}
    predict_fn = hace_predict_fn_perm(var_names, registro, j)

    real_ad, pred_ad, meta = build_fold_adatas(
        line, panel, predict_fn, var_names, seed_real, seed_em,
        cell_type, gene, n_genes, IFNB_PATH)
    GV.renombra(real_ad); GV.renombra(pred_ad)
    h_real = GV.huella(real_ad); h_pred = GV.huella(pred_ad)
    if ref is None:
        real_ad.write_h5ad(real_path)
    else:
        assert h_real == ref, f"{line}/j={j}: la POBLACION REAL se movio"
        if not os.path.exists(real_path):
            real_ad.write_h5ad(real_path)
    del real_ad

    datos = pred_ad.X.tocsr().data
    contrato = {"nnz": int(pred_ad.X.nnz), "dtype": str(pred_ad.X.dtype),
                "todo_entero": bool(np.all(datos == np.round(datos))),
                "min_valor": float(datos.min()) if datos.size else None,
                "max_valor": float(datos.max()) if datos.size else None,
                "hay_negativos": bool((datos < 0).any())}
    assert contrato["todo_entero"] and not contrato["hay_negativos"], contrato

    pred_path = os.path.join(WORK, f"pred_{line}_{j}.h5ad")
    pred_ad.write_h5ad(pred_path)
    pred_shape = list(pred_ad.shape)
    del pred_ad

    run_dir = os.path.join(WORK, f"run_{line}_{j}")
    sys.stderr.write(f"  [{line} j={j} em={seed_em} des={20272500+j}] run... ")
    sys.stderr.flush()
    p_run, dt_run = GV.corre(GV.cmd_run(pred_path, real_path, run_dir),
                             os.path.join(WORK, f"run_{line}_{j}.log"))
    if p_run.returncode != 0:
        sys.stderr.write("FALLO\n")
        os.remove(pred_path)
        return {"linea": line, "j": j, "fallo_en": "run",
                "mensaje_literal": p_run.stderr.strip()[-2000:]}

    agg = pd.read_csv(os.path.join(run_dir, "agg_results.csv"))
    fila = agg.loc[agg["statistic"] == "mean"].iloc[0]
    u = {m: (None if (m not in agg.columns or pd.isna(fila[m])) else float(fila[m]))
         for m in VCC2026_METRICS}
    coh = GV.cohortes(run_dir)
    desg = GV.desglose_replegados(run_dir)

    out_csv = os.path.join(WORK, f"scored_{line}_{j}.csv")
    p_sc, dt_sc = GV.corre(
        GV.cmd_score(run_dir, os.path.join(BUNDLES, f"BUNDLE_{line}"), out_csv),
        os.path.join(WORK, f"score_{line}_{j}.log"))
    score_ok = p_sc.returncode == 0
    sc = GV.leer_scored(out_csv) if score_ok else None
    sys.stderr.write(f"run {dt_run:.0f}s score {'OK' if score_ok else 'SE NEGO'}")

    rec = {"linea": line, "fold_idx": fold_idx, "j": j, "brazo": "PERMUTADO",
           "seed_real": seed_real, "seed_emision": seed_em,
           "seed_desarreglo": 20272500 + j,
           "huella_real": h_real, "huella_pred": h_pred,
           "poblacion_real_identica": (ref is None or h_real == ref),
           "pred_shape": pred_shape, "contrato_emision": contrato,
           "n_control_cells_pool": meta["n_control_cells"],
           "hash_run": GV.hash_del_run(run_dir),
           "u": u, "cohortes": coh, "desglose_replegados": desg,
           "desarreglo": registro.get("desarreglo"),
           "score_ok": score_ok, "scored": sc,
           "run_duration_s": round(dt_run, 2),
           "total_duration_s": round(time.time() - t0, 2)}
    if not score_ok:
        rec["fallo_en"] = "score"
        rec["mensaje_literal"] = p_sc.stderr.strip()[-2000:]

    for q in (pred_path, os.path.join(run_dir, "de_real.parquet"),
              os.path.join(run_dir, "de_pred.parquet")):
        if os.path.exists(q):
            os.remove(q)

    fr = ((sc or {}).get("filas") or {}).get("avg_score", {})
    sys.stderr.write(f"  avg={fr.get('from_replicate')}  "
                     f"fid={coh['n_fid']} union={coh['n_union_no_vacia']} "
                     f"DEpred={coh['suma_n_DE_pred']}\n")
    return rec


def main():
    os.makedirs(WORK, exist_ok=True)
    panel = json.load(open(PANEL_PATH))["panel_comun"]
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)
    sys.stderr.write(f"panel={len(panel)} n_genes={n_genes} "
                     f"n_desarreglos={N_SEMILLAS} out={OUT_JSONL}\n")

    hechas = {}
    if os.path.exists(OUT_JSONL):
        for ln in open(OUT_JSONL):
            if ln.strip():
                r = json.loads(ln)
                if not r.get("fallo_en"):
                    hechas[(r["linea"], r["j"])] = r.get("huella_real")

    for line in (sys.argv[1:] or LINES):
        fold_idx = LINES.index(line)
        real_path = os.path.join(WORK, f"real_{line}.h5ad")
        ref = None
        sys.stderr.write(f"\n=== {line} (fold {fold_idx}) PERMUTADO ===\n")
        for j in range(N_SEMILLAS):
            if (line, j) in hechas:
                ref = hechas[(line, j)] or ref
                sys.stderr.write(f"  [{line} j={j}] ya hecha\n")
                continue
            rec = una_corrida(line, fold_idx, j, panel, cell_type, gene,
                              var_names, n_genes, real_path, ref)
            append_record(rec)
            if rec.get("fallo_en"):
                sys.stderr.write(f"\n*** PARAR: {rec['fallo_en']} en {line} j={j}\n"
                                 f"{rec['mensaje_literal']}\n")
                sys.exit(7)
            ref = ref or rec["huella_real"]
        if os.path.exists(real_path):
            os.remove(real_path)
    sys.stderr.write(f"\nlisto: {OUT_JSONL}\n")


if __name__ == "__main__":
    main()
