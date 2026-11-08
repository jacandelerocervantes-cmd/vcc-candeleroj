#!/usr/bin/env python
"""Comun a los cuatro bloques de la noche. No mide nada: solo paths, la
familia de semillas NUEVA con su comprobacion de choque EN CODIGO, el tope de
disco, y los armadores de tabla markdown.

🔴 FAMILIA DE SEMILLAS DE ESTA NOCHE: 20280000 +. Las gastadas que se
comprueban estan declaradas en GASTADAS_RANGOS y GASTADAS_EXPLICITAS: los
rangos son los que el encargo nombra (2026000x-2026600x, 2027000x-2027900x) y
las explicitas son las que el codigo del proyecto realmente construyo.
"""
import hashlib
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np

SRC = "/home/ubuntu/vcc/modelo/src"
sys.path.insert(0, SRC)

OUT = "/home/ubuntu/vcc/modelo/out"
WORK_BASE = "/datos/modelo/out/noche"
MIN_FREE_GIB = 10.0
PY = "/home/ubuntu/mf/envs/vcc/bin/python"
CELLEVAL2_PY = os.path.expanduser("~/celleval2_venv/bin/cell-eval2")

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
PANEL_PATH = os.path.join(OUT, "suelo_panel.json")
BUNDLES = "/datos/modelo/out/nativas_comp"
CONTROL_LABEL = "NT"
CONTROL_COMP = "non-targeting"
N_PER_LABEL = 400
SEED_BASE_REAL = 20260001          # NO es nueva: se reusa a proposito

VCC2026_METRICS = [
    "pds_cosine",
    "expr_mse_unbiased_capped_norm",
    "de_wilcoxon_lfc_nmae",
    "de_wilcoxon_direction_fidelity_yield_raw",
    "de_wilcoxon_direction_reach_raw",
    "de_wilcoxon_sig_jaccard",
]
METRICA_CORTA = {
    "pds_cosine": "pds",
    "expr_mse_unbiased_capped_norm": "mse",
    "de_wilcoxon_lfc_nmae": "lfc",
    "de_wilcoxon_direction_fidelity_yield_raw": "fid",
    "de_wilcoxon_direction_reach_raw": "reach",
    "de_wilcoxon_sig_jaccard": "jac",
}

# ---------------------------------------------------------------- semillas
# Las semillas NUEVAS de esta noche, una por uso, declaradas antes de medir.
SEMILLAS_NUEVAS = {
    # B1: desempate del recorte por potencia
    "b1_desempate": [20280101],
    # B2: submuestreo de los 400 controles por linea, y el bootstrap sobre genes
    "b2_control": [20280401 + f for f in range(6)],
    "b2_bootstrap": [20280451 + f for f in range(6)],
    # B4: poblacion real 2025, submuestreo de control, y las 8 del suelo
    "b4_real_pert": [20280201],
    "b4_control_sub": [20280202],
    "b4_suelo": [20280300 + j for j in range(8)],
}

GASTADAS_RANGOS = [(20260000, 20266009), (20270000, 20279009)]
def _gastadas_explicitas():
    s = set()
    bases = [20260000, 20260001, 20261000, 20262000, 20263000, 20264000,
             20265000, 20266000, 20271000, 20272000, 20272500, 20273000,
             20273900, 20274000, 20275000, 20276000, 20277000]
    for b in bases:
        for f in range(6):
            s.add(b + f)
            for j in range(16):
                s.add(b + 100 * f + j)
        for j in range(16):
            s.add(b + j)
    return s


def comprueba_semillas():
    """🔴 Devuelve las semillas nuevas ENTERAS (no su conteo) y la prueba de que
    ninguna choca: ni entre si, ni contra las familias gastadas, por los dos
    criterios (rango nombrado por el encargo y enumeracion del codigo)."""
    usos = {}
    for fam, vals in SEMILLAS_NUEVAS.items():
        for v in vals:
            usos.setdefault(int(v), []).append(fam)
    internos = {k: v for k, v in usos.items() if len(v) > 1}
    expl = _gastadas_explicitas()
    en_rango = [k for k in usos if any(a <= k <= b for a, b in GASTADAS_RANGOS)]
    en_expl = sorted(set(usos) & expl)
    return {
        "familia_nueva": "20280000 +",
        "semillas_nuevas_enteras": {k: sorted(v) for k, v in SEMILLAS_NUEVAS.items()},
        "lista_plana_ordenada": sorted(usos),
        "n_semillas_nuevas": len(usos),
        "choques_internos": internos,
        "rangos_gastados_comprobados": GASTADAS_RANGOS,
        "n_gastadas_enumeradas_del_codigo": len(expl),
        "choques_con_rangos_gastados": sorted(en_rango),
        "choques_con_enumeracion_del_codigo": en_expl,
        "sin_choques": (not internos) and (not en_rango) and (not en_expl),
        "nota_seed_real": (f"{SEED_BASE_REAL} + fold NO entra: se reusa A PROPOSITO "
                           "para que la poblacion real sea la misma del suelo, gwps "
                           "REAL/PERMUTADO y los oraculos; se verifica por sha256."),
    }


# ---------------------------------------------------------------- utilidades
def free_gib(path="/datos"):
    return shutil.disk_usage(path)[2] / (1024 ** 3)


def guarda_disco(path="/datos"):
    libre = free_gib(path)
    if libre < MIN_FREE_GIB:
        raise SystemExit(f"PARAR: {libre:.2f} GiB libres en {path} (< {MIN_FREE_GIB}).")
    return libre


def log(msg):
    sys.stderr.write(f"[{time.strftime('%H:%M:%S')}] {msg}\n")
    sys.stderr.flush()


def corre(cmd, log_path):
    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    with open(log_path, "w") as fh:
        fh.write("CMD: " + " ".join(cmd) + "\n\nSTDOUT:\n" + proc.stdout +
                 "\n\nSTDERR:\n" + proc.stderr +
                 f"\nEXITCODE:{proc.returncode}\nDURATION_S:{dt:.2f}\n")
    return proc, dt


def huella_adata(adata):
    """sha256 de la matriz y de las etiquetas: la MISMA definicion que
    gwps_via.huella / dispersion_u.huella, para poder cruzar contra las huellas
    ya guardadas en los json de disco."""
    X = adata.X.tocsr()
    m = hashlib.sha256()
    for arr in (X.data, X.indices, X.indptr):
        m.update(np.ascontiguousarray(arr).tobytes())
    m.update("|".join(map(str, adata.obs["target_gene"].tolist())).encode())
    return m.hexdigest()


def sha256_archivo(path, bufsize=16 * 1024 * 1024):
    h = hashlib.sha256()
    n = 0
    with open(path, "rb") as fh:
        while True:
            b = fh.read(bufsize)
            if not b:
                break
            h.update(b)
            n += len(b)
    return h.hexdigest(), n


def borra(*paths):
    """Borra intermedios grandes. Nunca falla."""
    for p in paths:
        try:
            if p and os.path.isdir(p):
                shutil.rmtree(p, ignore_errors=True)
            elif p and os.path.exists(p):
                os.remove(p)
        except Exception:
            pass


def escribe_json(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=1, default=str)
    os.replace(tmp, path)


def escribe_md(path, lineas):
    tmp = path + ".tmp"
    with open(tmp, "w") as fh:
        fh.write("\n".join(lineas) + "\n")
    os.replace(tmp, path)


# ---------------------------------------------------------------- markdown
def tabla(encab, filas):
    out = ["| " + " | ".join(str(c) for c in encab) + " |",
           "|" + "|".join(["---"] * len(encab)) + "|"]
    for f in filas:
        out.append("| " + " | ".join("" if c is None else str(c) for c in f) + " |")
    return out


def f4(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/d"
    return f"{x:+.4f}"


def f4a(x):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/d"
    return f"{x:.4f}"


def fe(x, nd=3):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "n/d"
    return f"{x:.{nd}e}"


def rango(vals):
    v = [x for x in vals if x is not None and np.isfinite(x)]
    if not v:
        return "n/d"
    return f"[{min(v):+.4f}, {max(v):+.4f}]"


def rango_i(vals):
    v = [x for x in vals if x is not None]
    if not v:
        return "n/d"
    return f"[{min(v)}, {max(v)}]"


def dur(s):
    s = float(s)
    if s < 90:
        return f"{s:.1f} s"
    return f"{s/60:.1f} min ({s:.0f} s)"


# ---------------------------------------------------------------- arnes
def carga_panel():
    return json.load(open(PANEL_PATH))["panel_comun"]


def arma_linea(line, fold_idx, panel, var_names, cell_type, gene, n_genes):
    """Arma la poblacion real de la linea con `suelo_arnes.build_fold_adatas`
    SIN MODIFICARLO y DEVUELVE TAMBIEN el pool de control, capturandolo del
    propio contrato `predict_fn(control_pool_X, labels, rng)`. Asi una sola
    lectura del bloque de la linea sirve para todas las emisiones del bloque.

    La emision que se le pasa al arnes es el SUELO con la semilla de
    `dispersion_u` (20264000 + 100*fold), de modo que `pred_ad` es la emision
    de suelo j=0 ya puntuada, no una nueva."""
    from suelo_arnes import build_fold_adatas, IFNB_PATH
    caja = {}

    def fn_suelo_capturando(control_pool_X, labels, rng):
        caja["control_pool"] = control_pool_X
        n = control_pool_X.shape[0]
        return {l: control_pool_X[rng.integers(0, n, size=N_PER_LABEL)].tocsr()
                for l in labels}

    seed_real = SEED_BASE_REAL + fold_idx
    seed_pred = 20264000 + 100 * fold_idx          # suelo j=0 de dispersion_u
    real_ad, pred_suelo, meta = build_fold_adatas(
        line, panel, fn_suelo_capturando, var_names, seed_real, seed_pred,
        cell_type, gene, n_genes, IFNB_PATH)
    meta["seed_pred_suelo_j0"] = seed_pred
    return real_ad, pred_suelo, caja["control_pool"], meta


def huellas_guardadas():
    """Lee las `huella_real` ya guardadas en los json de disco, por linea, para
    poder verificar que la poblacion real de esta noche es la MISMA."""
    ref = {}
    for nombre in ("dispersion_u.json", "gwps_via.json", "pds_objetivo.json",
                   "control_permutado.json"):
        p = os.path.join(OUT, nombre)
        if not os.path.exists(p):
            continue
        for ln in open(p):
            ln = ln.strip()
            if not ln:
                continue
            try:
                r = json.loads(ln)
            except Exception:
                continue
            L, h = r.get("linea"), r.get("huella_real")
            if L and h:
                ref.setdefault(L, {}).setdefault(h, []).append(nombre)
    return ref
