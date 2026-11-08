#!/usr/bin/env python
"""CANDIDATO 01: medicion. 6 lineas x 8 semillas de emision, puntuadas contra
los 6 bundles de competencia.

LA EMISION es `cand01_covar.hace_predict_fn` -- respuesta lineal a lo largo de
la covarianza del propio pool de control. Ese modulo documenta su estructura y
las dos versiones de magnitud que se midieron y se descartaron ANTES de puntuar
nada.

🔴 LA REFUTACION YA ESTA MEDIDA, y esta medicion se corre a pesar de ella, no
para rescatarla. Sobre los 6934 genes expresados de A549, con Spearman del delta
de pseudobulk:
    A) mitad-contra-mitad del REAL (techo alcanzable) ... mediana +0.1575
    B) real de un blanco vs real de OTRO (generica) ..... mediana +0.2004
    C) esta prediccion vs el real del MISMO blanco ...... mediana +0.0014
    D) esta prediccion vs el real de OTRO blanco ........ mediana -0.0043
A > 0 demuestra que el test tiene potencia; C ~ 0 es la refutacion. No se itera
ni se ajusta nada: se mide una vez y se reporta.

QUE VARIA Y QUE NO:
  * VARIA: la semilla de emision. Familia NUEVA `20266000 + 100*fold + j`,
    j = 0..7. Ninguna colision con 2026000x (verdad), 2026100x (suelo),
    2026200x/2026300x (techo), 2026400x (dispersion de u), 2026500x
    (dispersion de b) ni la 0 del baseline/anchor.
  * NO VARIA: la poblacion real. `build_fold_adatas` SIN MODIFICAR con
    `seed_real = 20260001 + fold_idx`, etiqueta `non-targeting`. Se comprueba
    con sha256: una huella del real por linea, ocho del predicho.
  * NO VARIA: la ESTRUCTURA del candidato. `construye_respuesta` no usa rng --
    los ratios por blanco son deterministicos dado el pool de control -- asi que
    `sd_s` mide solo el muestreo de la emision, que es lo que se quiere.

NO TRANSDUCTIVO, verificado en el codigo del arnes: `build_fold_adatas` llama
`predict_fn(control_pool, panel + [CONTROL_LABEL], rng_pred)`
(`suelo_arnes.py:176`). El pool de control sale de `gene_local == CONTROL_LABEL`
(147-148), asi que la funcion no recibe ni la identidad de la linea ni una sola
celula perturbada. Tampoco entra ninguna cifra del `b` ni del anchor.
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (  # reusados sin modificar
    IFNB_PATH, CONTROL_LABEL, CELLEVAL2_PY, VCC2026_METRICS,
    build_fold_adatas, load_global_obs, load_var_names,
)
from cand01_covar import hace_predict_fn
from techo_arnes import contar_de_por_perturbacion   # sin modificar

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
OUT_JSONL = os.environ.get("CAND01_OUT",
                           "/home/ubuntu/vcc/modelo/out/candidato_01.json")
WORK = os.environ.get("CAND01_WORK", "/datos/modelo/out/cand01")
BUNDLES = "/datos/modelo/out/nativas_comp"
MIN_FREE_GIB = 10.0

CONTROL_COMP = "non-targeting"
SEED_BASE_REAL = 20260001
SEED_BASE_EMISION = 20266000     # familia nueva: 20266000 + 100*fold + j
N_SEMILLAS = 8

FID = "de_wilcoxon_direction_fidelity_yield_raw"
REACH = "de_wilcoxon_direction_reach_raw"
LFC = "de_wilcoxon_lfc_nmae"
JAC = "de_wilcoxon_sig_jaccard"


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
        fh.write("CMD: " + " ".join(cmd) + "\n\nSTDOUT:\n" + proc.stdout
                 + "\n\nSTDERR:\n" + proc.stderr
                 + f"\nEXITCODE:{proc.returncode}\nDURATION_S:{dt:.2f}\n")
    return proc, dt


def huella(adata):
    X = adata.X.tocsr()
    m = hashlib.sha256()
    for arr in (X.data, X.indices, X.indptr):
        m.update(np.ascontiguousarray(arr).tobytes())
    m.update("|".join(map(str, adata.obs["target_gene"].tolist())).encode())
    return m.hexdigest()


def renombra(adata):
    n = int((adata.obs["target_gene"] == CONTROL_LABEL).sum())
    adata.obs["target_gene"] = adata.obs["target_gene"].replace(
        {CONTROL_LABEL: CONTROL_COMP}
    )
    assert int((adata.obs["target_gene"] == CONTROL_COMP).sum()) == n, \
        "el renombrado movio filas"
    return n


def cmd_run(pred_path, real_path, outdir):
    """Sin `--control`: el valor lo pone el preset (`non-targeting`) y el
    `config_hash` no se desvia. `--write-degenes` esta verificado que no lo
    desvia (`dispersion_u.md` §1) y es lo que permite MEDIR las dos cohortes
    moviles en vez de inferirlas."""
    return [
        CELLEVAL2_PY, "run",
        "-ap", pred_path, "-ar", real_path,
        "--preset", "vcc2026", "--profile", "vcc2026",
        "--pert-col", "target_gene",
        "--write-degenes",
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
    from cell_eval2.config import EvalConfig
    from cell_eval2 import competition
    p = os.path.join(run_dir, "run_params.yaml")
    if not os.path.exists(p):
        return {"coincide": None}
    cfg = EvalConfig.from_yaml(p)
    obtenido = competition.rule_config_hash(cfg)
    referencia = competition.competition_payload()["config_hash"]
    return {"obtenido": obtenido, "referencia": referencia,
            "coincide": obtenido == referencia}


def leer_scored(path):
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path)
    cols = [c for c in df.columns if c != "metric"]
    filas = {}
    for _, row in df.iterrows():
        filas[str(row["metric"])] = {
            c: (None if pd.isna(row[c]) else
                (float(row[c]) if isinstance(row[c], (int, float, np.floating))
                 else row[c]))
            for c in cols
        }
    return {"columnas": cols, "filas": filas}


def cohortes(run_dir):
    res = pd.read_csv(os.path.join(run_dir, "results.csv"))

    def no_nan(m):
        d = res[res.metric == m]["value"].astype(float)
        return int((~d.isna()).sum()), int(len(d))

    n_fid, t_fid = no_nan(FID)
    n_reach, t_reach = no_nan(REACH)
    n_lfc, t_lfc = no_nan(LFC)
    n_jac, t_jac = no_nan(JAC)

    cuentas = contar_de_por_perturbacion(run_dir, control=CONTROL_COMP)
    vacias = sorted(p for p, c in cuentas.items() if c["union_vacia"])
    con_pred = sorted(p for p, c in cuentas.items() if c["n_DE_pred"] > 0)
    return {
        "n_fid": n_fid, "filas_fid": t_fid,
        "n_reach": n_reach, "filas_reach": t_reach,
        "n_lfc": n_lfc, "filas_lfc": t_lfc,
        "n_jac": n_jac, "filas_jac": t_jac,
        "n_union_no_vacia": len(cuentas) - len(vacias),
        "E_union_vacia": len(vacias),
        "perts_union_vacia": vacias,
        "n_perts_con_pred_sig": len(con_pred),
        "perts_con_pred_sig": con_pred,
        "suma_n_DE_pred": int(sum(c["n_DE_pred"] for c in cuentas.values())),
        "suma_n_DE_real": int(sum(c["n_DE_real"] for c in cuentas.values())),
        "suma_n_DE_interseccion": int(sum(c["n_DE_interseccion"]
                                          for c in cuentas.values())),
        "n_perts_total": len(cuentas),
        "por_pert": {p: {k: c[k] for k in ("n_DE_real", "n_DE_pred",
                                           "n_DE_interseccion", "n_DE_union")}
                     for p, c in sorted(cuentas.items())},
    }


def una_corrida(line, fold_idx, j, panel, cell_type, gene, var_names, n_genes,
                real_path, real_huella_ref, predict_fn, registro):
    t0 = time.time()
    libre = free_gib(WORK)
    if libre < MIN_FREE_GIB:
        raise SystemExit(f"PARAR: {libre:.1f} GiB libres en /datos.")

    seed_real = SEED_BASE_REAL + fold_idx
    seed_em = SEED_BASE_EMISION + 100 * fold_idx + j

    real_ad, pred_ad, meta = build_fold_adatas(
        line, panel, predict_fn, var_names, seed_real, seed_em,
        cell_type, gene, n_genes, IFNB_PATH,
    )
    renombra(real_ad)
    renombra(pred_ad)
    h_real = huella(real_ad)
    h_pred = huella(pred_ad)
    if real_huella_ref is None:
        real_ad.write_h5ad(real_path)
    else:
        assert h_real == real_huella_ref, (
            f"{line}/j={j}: la POBLACION REAL se movio -- "
            f"{h_real[:16]} != {real_huella_ref[:16]}"
        )
        if not os.path.exists(real_path):
            real_ad.write_h5ad(real_path)
    del real_ad

    # el contrato de emision, comprobado en cada corrida y no asumido
    datos = pred_ad.X.tocsr().data
    contrato = {
        "nnz": int(pred_ad.X.nnz),
        "dtype": str(pred_ad.X.dtype),
        "todo_entero": bool(np.all(datos == np.round(datos))),
        "min_valor": float(datos.min()) if datos.size else None,
        "max_valor": float(datos.max()) if datos.size else None,
        "hay_negativos": bool((datos < 0).any()),
    }
    assert contrato["todo_entero"] and not contrato["hay_negativos"], contrato

    pred_path = os.path.join(WORK, f"pred_{line}_{j}.h5ad")
    pred_ad.write_h5ad(pred_path)
    pred_shape = list(pred_ad.shape)
    del pred_ad

    run_dir = os.path.join(WORK, f"run_{line}_{j}")
    sys.stderr.write(f"  [{line} j={j} seed={seed_em}] run... ")
    sys.stderr.flush()
    p_run, dt_run = corre(cmd_run(pred_path, real_path, run_dir),
                          os.path.join(WORK, f"run_{line}_{j}.log"))
    if p_run.returncode != 0:
        sys.stderr.write("FALLO\n")
        os.remove(pred_path)
        return {"linea": line, "fold_idx": fold_idx, "j": j,
                "seed_emision": seed_em, "fallo_en": "run",
                "mensaje_literal": p_run.stderr.strip()[-2000:]}

    agg = pd.read_csv(os.path.join(run_dir, "agg_results.csv"))
    fila = agg.loc[agg["statistic"] == "mean"].iloc[0]
    u = {m: (None if (m not in agg.columns or pd.isna(fila[m])) else float(fila[m]))
         for m in VCC2026_METRICS}

    coh = cohortes(run_dir)

    out_csv = os.path.join(WORK, f"scored_{line}_{j}.csv")
    p_sc, dt_sc = corre(cmd_score(run_dir, os.path.join(BUNDLES, f"BUNDLE_{line}"),
                                  out_csv),
                        os.path.join(WORK, f"score_{line}_{j}.log"))
    score_ok = p_sc.returncode == 0
    sc = leer_scored(out_csv) if score_ok else None
    sys.stderr.write(f"run {dt_run:.0f}s score {'OK' if score_ok else 'SE NEGO'}")

    rec = {
        "linea": line, "fold_idx": fold_idx, "j": j,
        "seed_real": seed_real, "seed_emision": seed_em,
        "huella_real": h_real, "huella_pred": h_pred,
        "poblacion_real_identica": (real_huella_ref is None
                                    or h_real == real_huella_ref),
        "pred_shape": pred_shape,
        "contrato_emision": contrato,
        "n_control_cells_pool": meta["n_control_cells"],
        "hash_run": hash_del_run(run_dir),
        "u": u,
        "cohortes": coh,
        "score_ok": score_ok,
        "scored": sc,
        "score_stderr_literal": p_sc.stderr.strip()[-1500:] if not score_ok else "",
        "run_duration_s": round(dt_run, 2),
        "score_duration_s": round(dt_sc, 2),
        "total_duration_s": round(time.time() - t0, 2),
    }
    if j == 0:
        rec["diag_estructura"] = registro.get("diag")
    if not score_ok:
        rec["fallo_en"] = "score"
        rec["mensaje_literal"] = p_sc.stderr.strip()[-2000:]

    for q in (pred_path,
              os.path.join(run_dir, "de_real.parquet"),
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
    sys.stderr.write(f"panel={len(panel)} n_genes={n_genes}\n")

    hechas = {}
    if os.path.exists(OUT_JSONL):
        with open(OUT_JSONL) as fh:
            for ln in fh:
                if ln.strip():
                    r = json.loads(ln)
                    if not r.get("fallo_en"):
                        hechas[(r["linea"], r["j"])] = r.get("huella_real")
    if hechas:
        sys.stderr.write(f"ya hechas: {len(hechas)} corridas\n")

    lineas = sys.argv[1:] or LINES
    for line in lineas:
        fold_idx = LINES.index(line)
        real_path = os.path.join(WORK, f"real_{line}.h5ad")
        ref = None
        # un predict_fn por linea: su cache guarda los ratios del pool de esa
        # linea, que son deterministicos y se reusan en las 8 semillas
        registro = {}
        predict_fn = hace_predict_fn(var_names, registro)
        sys.stderr.write(f"\n=== {line} (fold {fold_idx}) ===\n")
        for j in range(N_SEMILLAS):
            if (line, j) in hechas:
                ref = hechas[(line, j)] or ref
                sys.stderr.write(f"  [{line} j={j}] ya hecha\n")
                continue
            rec = una_corrida(line, fold_idx, j, panel, cell_type, gene,
                              var_names, n_genes, real_path, ref,
                              predict_fn, registro)
            append_record(rec)
            if rec.get("fallo_en"):
                sys.stderr.write(
                    f"\n*** PARAR: fallo en `{rec['fallo_en']}` en {line} j={j}. "
                    "No se repara nada.\n*** mensaje literal:\n"
                    f"{rec['mensaje_literal']}\n"
                )
                sys.exit(7)
            ref = ref or rec["huella_real"]
        if os.path.exists(real_path):
            os.remove(real_path)

    sys.stderr.write(f"\nlisto: {OUT_JSONL}\n")


if __name__ == "__main__":
    main()
