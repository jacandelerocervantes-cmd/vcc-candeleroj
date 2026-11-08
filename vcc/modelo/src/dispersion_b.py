#!/usr/bin/env python
"""La dispersion de `b`: el tercer termino del margen, que hasta hoy faltaba.

POR QUE: el margen de `dispersion_u.md` §2 es
    margen = sqrt( sd_s^2 + (|s|*SE_r/|r-b|)^2 )
y OMITE el termino de `b`, porque `b` es UNA sola realizacion con `--seed 0`.
Mientras lo sea, todo margen es cota inferior. Aqui se estima esa cantidad.

QUE VARIA Y QUE NO -- es lo unico que hace valida la medicion:
  * VARIA: la semilla del baseline. Familia NUEVA, `20265000 + 100*fold + j`,
    j = 0..7. Ninguna colision con 2026000x (verdad), 2026100x (suelo),
    2026200x/2026300x (techo), 2026400x (dispersion de u) ni la 0 del
    baseline/anchor.
  * NO VARIA: la poblacion real. `build_fold_adatas` sin modificar con
    `seed_real = 20260001 + fold_idx`, etiqueta `non-targeting`. Se COMPRUEBA
    igual que en `dispersion_u.md`: sha256 de la matriz real (data, indices,
    indptr) y de sus etiquetas, una sola huella por linea, y ocho del h5ad que
    el baseline EMITE (`--save-pred`).
  * NO VARIA: preset, perfil, `--emit dispersed`, ausencia de `--control`.

EL `b` DE SEMILLA 0 que ya esta en los bundles entra como una observacion mas y
se marca como tal. No se sustituye ni se promedia dentro del bundle: los 6
bundles de competencia se quedan exactamente como estan.

LA `E` Y LAS COHORTES, sin una segunda pasada del marcador. `baseline_results.csv`
trae el valor por perturbacion de las seis metricas, asi que:
  * la cohorte de `fid`, `reach` y `lfc_nmae` es el conteo de filas no-NaN;
  * la `E` (perturbaciones con `|R u P| = 0`, que cobran `jac = 1.0` por
    convencion, D-023) se lee como `jac == 1.0`.
El segundo punto es un PROXY, y queda VALIDADO antes de usarse: contra los
conteos de union de `huecos_cohorte.json` -- que salieron de `de_real.parquet` y
`de_pred.parquet` con `--write-degenes` -- el conjunto `{jac == 1.0}` coincide
exactamente con `{|union| = 0}` en las seis lineas con semilla 0 (15/0/14/3/7/0).
Ademas cada semilla lleva su propia GUARDA, no una suposicion: si la union es
vacia entonces `n_real = n_pred = 0`, luego `max(n_pred, n_conf) = 0` y `fid`
sale NaN (`direction.py:776` con `_safe_div`, 689-699). Asi que toda
perturbacion con `jac == 1.0` tiene que tener `fid` NaN. Si alguna no la tiene,
el proxy fallo en esa semilla y se reporta en vez de taparse.

No modifica ningun arnes: `build_fold_adatas` se usa tal cual.
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
from suelo_floor import predict_floor   # solo para satisfacer la firma; su
                                        # emision NO se usa aqui

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
PANEL_PATH = "/home/ubuntu/vcc/modelo/out/suelo_panel.json"
# Se puede partir en flujos paralelos, uno por fichero, para no entrelazar
# escrituras en un fichero que sostiene cifras. Se fusionan al final.
OUT_JSONL = os.environ.get("DISPB_OUT",
                           "/home/ubuntu/vcc/modelo/out/dispersion_b.json")
WORK = "/datos/modelo/out/dispb"
BASE_OFICIAL = "/datos/modelo/out/nativas_comp/BASE_{L}"
MIN_FREE_GIB = 10.0

CONTROL_COMP = "non-targeting"
SEED_BASE_REAL = 20260001
SEED_BASE_BASELINE = 20265000     # familia nueva: 20265000 + 100*fold + j
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


def huella_ad(adata):
    """sha256 de la matriz y de las etiquetas."""
    X = adata.X.tocsr() if hasattr(adata.X, "tocsr") else adata.X
    m = hashlib.sha256()
    if hasattr(X, "data"):
        for arr in (X.data, X.indices, X.indptr):
            m.update(np.ascontiguousarray(arr).tobytes())
    else:
        m.update(np.ascontiguousarray(X).tobytes())
    m.update("|".join(map(str, adata.obs["target_gene"].tolist())).encode())
    return m.hexdigest()


def huella_fichero(path):
    m = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            m.update(chunk)
    return m.hexdigest()


def renombra(adata):
    n = int((adata.obs["target_gene"] == CONTROL_LABEL).sum())
    adata.obs["target_gene"] = adata.obs["target_gene"].replace(
        {CONTROL_LABEL: CONTROL_COMP}
    )
    assert int((adata.obs["target_gene"] == CONTROL_COMP).sum()) == n, \
        "el renombrado movio filas"
    return n


def cmd_baseline(real_path, outdir, pred_path, seed):
    """Sin `--control`: el valor lo pone el preset (`non-targeting`) y el
    `config_hash` no se desvia (`anclas_nativas.md` aislo que `'control'` NO
    esta en `_RULE_EXCLUDED`, `competition.py:112`)."""
    return [
        CELLEVAL2_PY, "baseline",
        "-ar", real_path,
        "--preset", "vcc2026", "--profile", "vcc2026",
        "--pert-col", "target_gene",
        "--emit", "dispersed", "--seed", str(seed),
        "--save-pred", pred_path,
        "-o", outdir,
    ]


def hash_del_dir(d):
    from cell_eval2.config import EvalConfig
    from cell_eval2 import competition
    p = os.path.join(d, "run_params.yaml")
    if not os.path.exists(p):
        return {"coincide": None}
    cfg = EvalConfig.from_yaml(p)
    obtenido = competition.rule_config_hash(cfg)
    referencia = competition.competition_payload()["config_hash"]
    return {"obtenido": obtenido, "referencia": referencia,
            "coincide": obtenido == referencia}


def lee_agg(outdir):
    agg = pd.read_csv(os.path.join(outdir, "baseline_agg.csv"))
    if "statistic" in agg.columns:
        fila = agg.loc[agg["statistic"] == "mean"].iloc[0]
        return {m: (None if (m not in agg.columns or pd.isna(fila[m]))
                    else float(fila[m])) for m in VCC2026_METRICS}
    # formato largo
    d = {}
    for _, r in agg.iterrows():
        d[str(r[agg.columns[0]])] = (None if pd.isna(r[agg.columns[-1]])
                                    else float(r[agg.columns[-1]]))
    return {m: d.get(m) for m in VCC2026_METRICS}


def lee_por_pert(outdir):
    """`baseline_results.csv` es largo: perturbation, metric, value."""
    df = pd.read_csv(os.path.join(outdir, "baseline_results.csv"))
    piv = df.pivot_table(index="perturbation", columns="metric",
                         values="value", aggfunc="first")
    return piv


def cohortes_y_E(piv):
    """Las cohortes son conteos de no-NaN. La `E` es `jac == 1.0`, con su guarda."""
    def no_nan(m):
        if m not in piv.columns:
            return None, None
        col = piv[m].astype(float)
        return int((~col.isna()).sum()), int(len(col))

    n_fid, t_fid = no_nan(FID)
    n_reach, t_reach = no_nan(REACH)
    n_lfc, t_lfc = no_nan(LFC)
    n_jac, t_jac = no_nan(JAC)

    jac = piv[JAC].astype(float)
    e_perts = sorted(jac.index[jac == 1.0].tolist())
    # GUARDA del proxy: union vacia => fid NaN. Si alguna E tiene fid no-NaN,
    # el proxy fallo en esta semilla.
    fid = piv[FID].astype(float) if FID in piv.columns else None
    viola = []
    if fid is not None:
        for p in e_perts:
            if not np.isnan(float(fid.loc[p])):
                viola.append(p)
    con_fid = sorted(piv.index[~piv[FID].astype(float).isna()].tolist()) \
        if FID in piv.columns else []
    return {
        "n_fid": n_fid, "filas_fid": t_fid,
        "n_reach": n_reach, "filas_reach": t_reach,
        "n_lfc": n_lfc, "filas_lfc": t_lfc,
        "n_jac": n_jac, "filas_jac": t_jac,
        "E_union_vacia": len(e_perts),
        "perts_union_vacia": e_perts,
        "n_union_no_vacia": int(len(jac)) - len(e_perts),
        "proxy_violaciones": viola,
        "proxy_ok": len(viola) == 0,
        "perts_cohorte_fid": con_fid,
        "n_perts_total": int(len(piv)),
    }


def una_corrida(line, fold_idx, j, panel, cell_type, gene, var_names, n_genes,
                real_path, real_huella_ref):
    t0 = time.time()
    libre = free_gib(WORK)
    if libre < MIN_FREE_GIB:
        raise SystemExit(f"PARAR: {libre:.1f} GiB libres en /datos.")

    seed_real = SEED_BASE_REAL + fold_idx
    seed_b = SEED_BASE_BASELINE + 100 * fold_idx + j

    # La poblacion real: identica en las 8, y se comprueba.
    real_ad, _pred_ad, meta = build_fold_adatas(
        line, panel, predict_floor, var_names, seed_real, seed_real,
        cell_type, gene, n_genes, IFNB_PATH,
    )
    del _pred_ad
    renombra(real_ad)
    h_real = huella_ad(real_ad)
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

    outdir = os.path.join(WORK, f"base_{line}_{j}")
    pred_path = os.path.join(WORK, f"bpred_{line}_{j}.h5ad")
    sys.stderr.write(f"  [{line} j={j} seed={seed_b}] baseline... ")
    sys.stderr.flush()
    p, dt = corre(cmd_baseline(real_path, outdir, pred_path, seed_b),
                  os.path.join(WORK, f"base_{line}_{j}.log"))
    if p.returncode != 0:
        sys.stderr.write("FALLO\n")
        return {"linea": line, "fold_idx": fold_idx, "j": j,
                "seed_baseline": seed_b, "fallo_en": "baseline",
                "mensaje_literal": p.stderr.strip()[-2000:]}

    b = lee_agg(outdir)
    piv = lee_por_pert(outdir)
    coh = cohortes_y_E(piv)
    h_pred = huella_fichero(pred_path) if os.path.exists(pred_path) else None

    # el valor por perturbacion, que es el intermedio que sostiene la E
    por_pert = {}
    for pname in piv.index:
        fila = {}
        for m, corto in ((JAC, "jac"), (FID, "fid"), (REACH, "reach"),
                         (LFC, "lfc_nmae"), ("pds_cosine", "pds")):
            if m in piv.columns:
                v = float(piv.loc[pname, m])
                fila[corto] = None if np.isnan(v) else v
        por_pert[str(pname)] = fila

    meta_json = {}
    mp = os.path.join(outdir, "baseline_meta.json")
    if os.path.exists(mp):
        meta_json = json.load(open(mp))

    rec = {
        "linea": line, "fold_idx": fold_idx, "j": j,
        "seed_real": seed_real, "seed_baseline": seed_b,
        "huella_real": h_real,
        "huella_pred_baseline": h_pred,
        "poblacion_real_identica": (real_huella_ref is None
                                    or h_real == real_huella_ref),
        "n_control_cells_pool": meta["n_control_cells"],
        "hash_regla": hash_del_dir(outdir),
        "b": b,
        "cohortes": coh,
        "por_pert": por_pert,
        "source_fingerprint_strict": meta_json.get("source_fingerprint_strict"),
        "baseline_duration_s": round(dt, 2),
        "total_duration_s": round(time.time() - t0, 2),
    }

    # limpieza: la emision del baseline no se conserva (ya tiene su huella)
    if os.path.exists(pred_path):
        os.remove(pred_path)

    sys.stderr.write(
        f"{dt:.0f}s  pds={b.get('pds_cosine')}  E={coh['E_union_vacia']}"
        f" fid={coh['n_fid']} proxy={'OK' if coh['proxy_ok'] else 'FALLO'}\n")
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
        sys.stderr.write(f"\n=== {line} (fold {fold_idx}) ===\n")
        for j in range(N_SEMILLAS):
            if (line, j) in hechas:
                ref = hechas[(line, j)] or ref
                sys.stderr.write(f"  [{line} j={j}] ya hecha\n")
                continue
            rec = una_corrida(line, fold_idx, j, panel, cell_type, gene,
                              var_names, n_genes, real_path, ref)
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
