#!/usr/bin/env python
"""Arnes del TECHO: el mismo marcador y la misma geometria que el suelo, con la
funcion de emision reemplazada. En vez de emitir una prediccion degenerada
(suelo_floor.predict_floor: 400 celulas del control remuestreadas con
reemplazo), `predict_techo` emite LA OTRA MITAD DE LAS CELULAS REALES de cada
perturbacion. Las seis metricas salen de `cell-eval2 run` igual que en el
suelo, asi que la tabla resta columna a columna contra bucle_suelo.md.

Lo que se reusa de suelo_arnes.py SIN modificarlo: load_global_obs,
load_var_names, load_line_block, read_rows_csr, score_vs_self,
read_agg_metrics, read_avg_score, VCC2026_METRICS, IFNB_PATH, CONTROL_LABEL,
N_PER_LABEL.

Lo que NO se pudo reusar, y por que (reportado en bucle_techo.md):
  * `build_fold_adatas`: su contrato de emision es
    `predict_fn(control_pool_X, labels, rng)` -- la funcion de emision recibe
    SOLO el pool de control, nunca las celulas perturbadas, asi que por
    construccion no puede emitir "la otra mitad de las celulas reales".
    Ademas fija el lado VERDAD el mismo (400 celulas sorteadas con seed_real),
    y el techo necesita gobernar LOS DOS lados para correr las dos
    direcciones. Aqui se escribe `build_techo_pair`, que ensambla los AnnData
    con las MISMAS reglas (columna obs `target_gene`, mismo var index,
    obs_names_make_unique) sin tocar el original.
  * `run_cell_eval2`: su argv es fijo y no admite banderas extra. El techo
    necesita dos: `--write-degenes` (para las columnas n_DE/|union| del punto
    5, que se cuentan con la definicion del propio marcador) y
    `--set control_source=pred` (la variante I). `run_cell_eval2_techo`
    reproduce ese argv y agrega solo esas banderas; `ARGV_SUELO_ESPERADO`
    deja la comparacion asertada.
"""
import json
import os
import subprocess
import sys
import time

import numpy as np
import scipy.sparse as sp
import anndata as ad
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (  # noqa: E402  -- reusados sin modificar
    IFNB_PATH, CONTROL_LABEL, N_PER_LABEL, CELLEVAL2_PY, VCC2026_METRICS,
    load_global_obs, load_var_names, load_line_block,
    score_vs_self, read_agg_metrics, read_avg_score,
)

# Semillas: la VERDAD del suelo usa 2026000x y la prediccion del suelo 2026100x.
# El techo no reusa ninguna de esas, y nunca 42.
SEED_BASE_SPLIT = 20262000   # particion A/B de las celulas reales, por pliegue
SEED_BASE_CTRLS = 20263000   # bootstrap del control de PRED en la variante S

# Umbrales de DE del perfil vcc2026 (configs/vcc2026.yaml), no inventados aqui.
DE_P_ADJ_THRESHOLD = 0.05
DE_MIN_ABS_LOG2FC = 0.0
DE_NAN_LFC_POLICY = "mask"
DE_SORT_BY = "abs_log2_fold_change"


def argv_suelo(pred_path, real_path, outdir, pert_col="target_gene", control=CONTROL_LABEL):
    """El argv EXACTO que suelo_arnes.run_cell_eval2 arma, reproducido aqui para
    poder asertar que el del techo es ese mas las banderas declaradas."""
    return [
        CELLEVAL2_PY, "run",
        "-ap", pred_path, "-ar", real_path,
        "--preset", "vcc2026", "--profile", "vcc2026",
        "--pert-col", pert_col, "--control", control,
        "-o", outdir,
    ]


def run_cell_eval2_techo(pred_path, real_path, outdir, control_source,
                         pert_col="target_gene", control=CONTROL_LABEL):
    """`cell-eval2 run` con el argv del suelo mas `--write-degenes` y, solo en la
    variante I, `--set control_source=pred`.

    control_source="real"  -> variante S: no se agrega ningun --set, asi que el
        argv es el del suelo mas --write-degenes. El preset vcc2026 ya trae
        control_source: real, de modo que la configuracion de puntuacion es
        identica a la del suelo.
    control_source="pred"  -> variante I: cada lado hace su DE contra su propio
        bloque de control.
    """
    os.makedirs(outdir, exist_ok=True)
    cmd = argv_suelo(pred_path, real_path, outdir, pert_col, control)
    cmd.append("--write-degenes")
    if control_source == "pred":
        cmd += ["--set", "control_source=pred"]
    elif control_source != "real":
        raise ValueError(f"control_source debe ser 'real' o 'pred', no {control_source!r}")

    t0 = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    dt = time.time() - t0
    log_path = os.path.join(outdir, "cell_eval2_run.log")
    with open(log_path, "w") as fh:
        fh.write("CMD: " + " ".join(cmd) + "\n\n")
        fh.write("STDOUT:\n" + proc.stdout + "\n\nSTDERR:\n" + proc.stderr + "\n")
        fh.write(f"\nEXITCODE:{proc.returncode}\nDURATION_S:{dt:.2f}\n")
    if proc.returncode != 0:
        raise RuntimeError(
            f"cell-eval2 run fallo (exit={proc.returncode}) para {outdir}. Ver {log_path}"
        )

    agg_path = os.path.join(outdir, "agg_results.csv")
    meta_path = os.path.join(outdir, "run_meta.json")
    params_path = os.path.join(outdir, "run_params.yaml")
    for p in (agg_path, meta_path, params_path):
        assert os.path.exists(p), f"no se escribio {p}"

    # Verificar que el --set surtio efecto: no asumirlo.
    with open(params_path) as fh:
        params_txt = fh.read()
    esperado = f"control_source: {control_source}"
    assert esperado in params_txt, (
        f"run_params.yaml de {outdir} no declara {esperado!r}; el --set no surtio efecto"
    )
    return agg_path, meta_path, params_path, dt


def split_mitades(idx_local, rng, ordinal_grupo):
    """Particion aleatoria SIN reemplazo de las celulas de un grupo en mitad A y
    mitad B, con desempate NO dirigido.

    Con n par las dos mitades tienen n/2. Con n impar sobra una celula: el
    desempate ALTERNA por grupo (ordinal par -> la extra va a A; ordinal impar
    -> va a B), de modo que ningun lado recibe sistematicamente mas potencia.
    Ninguna celula queda en las dos mitades y ninguna se descarta.
    """
    perm = rng.permutation(np.asarray(idx_local))
    n = perm.size
    h = n // 2
    if n % 2 == 0:
        return perm[:h], perm[h:]
    if ordinal_grupo % 2 == 0:
        return perm[:h + 1], perm[h + 1:]
    return perm[:h], perm[h:]


def predict_techo(mitad_emisora):
    """La funcion de emision del techo: devuelve LAS CELULAS REALES de la mitad
    emisora, sin modificarlas y sin remuestrearlas.

    Es la contraparte de suelo_floor.predict_floor, que para cada etiqueta
    devolvia 400 celulas del control remuestreadas con reemplazo. Aqui no hay
    remuestreo: emitir la mitad tal cual es lo que mide la reproducibilidad de
    la propia verdad. La diferencia de profundidad contra las 400 del suelo
    queda declarada, no corregida (ver bucle_techo.md).

    mitad_emisora: dict[label] -> csr_matrix. Devuelve un dict equivalente.
    """
    return {label: bloque for label, bloque in mitad_emisora.items()}


def _ensamblar_adata(bloques, var_names):
    """Mismas reglas de ensamblado que suelo_arnes.build_fold_adatas: columna obs
    `target_gene`, var index = var_names, obs_names_make_unique."""
    X = sp.vstack([b for _, b in bloques], format="csr")
    labels = []
    for label, b in bloques:
        labels.extend([label] * b.shape[0])
    obs = pd.DataFrame({"target_gene": labels})
    adata = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=var_names))
    adata.obs_names_make_unique()
    return adata


def partir_linea(line, panel, cell_type, gene, n_genes, seed_split, h5path=IFNB_PATH):
    """Lee el bloque de la linea y parte en mitades A/B las celulas reales de cada
    perturbacion del panel y del control. Devuelve (mitades, meta).

    mitades[label] = (csr_A, csr_B); mitades[CONTROL_LABEL] tambien, mas
    `control_pool` completo en meta para la variante S.
    """
    csr_local, gene_local = load_line_block(line, cell_type, gene, n_genes, h5path)
    rng = np.random.default_rng(seed_split)

    mitades = {}
    n_por_pert = {}
    for ordinal, pert in enumerate(panel):
        idx = np.where(gene_local == pert)[0]
        n_por_pert[pert] = int(idx.size)
        a_idx, b_idx = split_mitades(idx, rng, ordinal)
        mitades[pert] = (csr_local[a_idx], csr_local[b_idx])

    ctrl_idx = np.where(gene_local == CONTROL_LABEL)[0]
    n_por_pert[CONTROL_LABEL] = int(ctrl_idx.size)
    ca_idx, cb_idx = split_mitades(ctrl_idx, rng, len(panel))
    mitades[CONTROL_LABEL] = (csr_local[ca_idx], csr_local[cb_idx])
    control_pool = csr_local[ctrl_idx]

    meta = {
        "line": line,
        "seed_split": seed_split,
        "n_celulas_por_pert": n_por_pert,
        "n_control_cells": int(ctrl_idx.size),
        "n_mitadA_por_pert": {k: int(v[0].shape[0]) for k, v in mitades.items()},
        "n_mitadB_por_pert": {k: int(v[1].shape[0]) for k, v in mitades.items()},
    }
    return mitades, control_pool, meta


def build_techo_pair(mitades, control_pool, panel, var_names, variante, direccion,
                     seed_ctrl_boot):
    """Arma (real_ad, pred_ad) para una variante y una direccion.

    direccion="A": la mitad A es la PREDICCION y la mitad B es la VERDAD.
    direccion="B": al reves. No se promedian: se reportan las dos.

    variante="S" (control COMPARTIDO, techo optimista): el bloque de control de
        VERDAD es TODO el pool de control de la linea y el de PREDICCION son 400
        celulas de ese mismo pool remuestreadas con reemplazo -- la regla
        exacta del suelo. Con control_source=real (el default del preset
        vcc2026) los dos lados hacen su DE contra ese unico control, asi que sus
        errores quedan correlacionados.
    variante="I" (control PARTIDO, techo pesimista): el control tambien se parte
        y cada mitad usa el suyo, con control_source=pred.
    """
    if direccion not in ("A", "B"):
        raise ValueError(f"direccion debe ser 'A' o 'B', no {direccion!r}")
    i_pred, i_real = (0, 1) if direccion == "A" else (1, 0)

    emisor = {pert: mitades[pert][i_pred] for pert in panel}
    pred_perts = predict_techo(emisor)
    real_perts = {pert: mitades[pert][i_real] for pert in panel}

    if variante == "S":
        rng_boot = np.random.default_rng(seed_ctrl_boot)
        n_ctrl = control_pool.shape[0]
        elegidas = rng_boot.integers(0, n_ctrl, size=N_PER_LABEL)
        pred_ctrl = control_pool[elegidas]
        real_ctrl = control_pool
    elif variante == "I":
        pred_ctrl = mitades[CONTROL_LABEL][i_pred]
        real_ctrl = mitades[CONTROL_LABEL][i_real]
    else:
        raise ValueError(f"variante debe ser 'S' o 'I', no {variante!r}")

    pred_bloques = [(p, pred_perts[p]) for p in panel] + [(CONTROL_LABEL, pred_ctrl)]
    real_bloques = [(p, real_perts[p]) for p in panel] + [(CONTROL_LABEL, real_ctrl)]

    pred_ad = _ensamblar_adata(pred_bloques, var_names)
    real_ad = _ensamblar_adata(real_bloques, var_names)

    meta = {
        "variante": variante,
        "direccion": direccion,
        "control_source": "real" if variante == "S" else "pred",
        "n_pred_ctrl": int(pred_ctrl.shape[0]),
        "n_real_ctrl": int(real_ctrl.shape[0]),
        "pred_shape": list(pred_ad.shape),
        "real_shape": list(real_ad.shape),
        "n_pred_por_pert": {p: int(pred_perts[p].shape[0]) for p in panel},
        "n_real_por_pert": {p: int(real_perts[p].shape[0]) for p in panel},
    }
    return real_ad, pred_ad, meta


def contar_de_por_perturbacion(outdir, control=CONTROL_LABEL):
    """Columnas del punto 5: |R|, |P|, |R u P| y el jaccard recomputado, POR
    perturbacion, con la definicion del propio marcador.

    No reimplementa nada: llama a `cell_eval2.de.prepare_de` sobre las tablas
    que escribio `--write-degenes` y repite las dos lineas de filtrado de
    `cell_eval2.metrics.de.de_sig_jaccard` -- filtro `p_adj < T`, `.unique()`
    sobre (target, feature) y el anti-join que excluye el gen propio de cada
    target (issue #172) en LOS DOS lados.

    Razon de ser de estas columnas: `de_sig_jaccard` devuelve 1.0 cuando la
    union es vacia (convencion J(vacio, vacio)=1, D-023). Sin |R u P| a la
    vista, una perturbacion sin ninguna senal aparece con techo perfecto y es
    invisible.
    """
    import polars as pl
    from cell_eval2.de import prepare_de
    from cell_eval2.metrics.de import _exclude_own_gene

    de_real = os.path.join(outdir, "de_real.parquet")
    de_pred = os.path.join(outdir, "de_pred.parquet")
    for p in (de_real, de_pred):
        assert os.path.exists(p), f"falta {p} (--write-degenes no escribio la tabla)"

    prepared = prepare_de(
        de_pred, de_real, control=control, sort_by=DE_SORT_BY,
        p_adj_threshold=DE_P_ADJ_THRESHOLD, nan_lfc_policy=DE_NAN_LFC_POLICY,
        min_abs_log2fc=DE_MIN_ABS_LOG2FC,
    )
    T = prepared.p_adj_threshold

    real_sig = _exclude_own_gene(
        prepared.real_df.filter(pl.col("p_adj") < T).select("target", "feature").unique(),
        prepared, who="techo_arnes", side="real",
    )
    pred_sig = _exclude_own_gene(
        prepared.pred_df.filter(pl.col("p_adj") < T).select("target", "feature").unique(),
        prepared, who="techo_arnes", side="pred",
    )
    a_map = dict(zip(*real_sig.group_by("target").len().to_dict(as_series=False).values()))
    b_map = dict(zip(*pred_sig.group_by("target").len().to_dict(as_series=False).values()))
    inter = real_sig.join(pred_sig, on=["target", "feature"], how="inner").group_by("target").len()
    tp_map = dict(zip(*inter.to_dict(as_series=False).values()))

    filas = {}
    for p in prepared.perturbations:
        a = int(a_map.get(p, 0))
        b = int(b_map.get(p, 0))
        tp = int(tp_map.get(p, 0))
        union = a + b - tp
        filas[p] = {
            "n_DE_real": a,
            "n_DE_pred": b,
            "n_DE_interseccion": tp,
            "n_DE_union": union,
            "jaccard_recomputado": 1.0 if union == 0 else tp / union,
            "union_vacia": union == 0,
        }
    return filas


def leer_results_por_perturbacion(outdir):
    """results.csv es largo (perturbation, metric, value). Devuelve
    dict[pert][metric] = value, para cruzar el jaccard del marcador contra el
    recomputado y detectar cualquier divergencia."""
    path = os.path.join(outdir, "results.csv")
    assert os.path.exists(path), f"falta {path}"
    df = pd.read_csv(path)
    out = {}
    for pert, metric, value in zip(df["perturbation"], df["metric"], df["value"]):
        out.setdefault(str(pert), {})[str(metric)] = float(value)
    return out
