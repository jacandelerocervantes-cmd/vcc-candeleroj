#!/usr/bin/env python
"""CARACTERIZAR pds_cosine COMO OBJETIVO: los cuatro brazos diagnosticos.

🔴 ESTOS CUATRO BRAZOS NO SON CANDIDATOS. Los tres oraculos leen el LADO REAL
PERTURBADO de la MISMA poblacion que se puntua, asi que son TRANSDUCTIVOS y van
al mismo estante que el `b` oficial. SON UNA MEDICION DE LA METRICA, NUNCA UNA
EMISION QUE PUEDA ENVIARSE: ninguna de sus cifras puede entrar en una submission
ni en un rasgo de un candidato.

LOS BRAZOS
  1a. ORACULO COMPLETO POR COPIA: el bloque predicho de `h` ES el bloque real
      perturbado de `h`. Reproduce el delta real EXACTAMENTE, asi que es la cota
      superior y la comprobacion de cordura (pds_cosine debe dar 1.0).
  1b. ORACULO COMPLETO POR RAZON: `r_g = media_real_pert_g / media_control_g` en
      TODOS los genes, realizado con el mismo `_reescala_entero` que los brazos
      2 y 3. Es el DENOMINADOR comparable de la fraccion de cola: 1a y 1b se
      diferencian solo en la realizacion, y la distancia entre los dos mide lo
      que cuesta la realizacion por conteos enteros.
  2.  ORACULO DE COLA: la misma razon SOLO en los genes del conjunto ED real de
      cada perturbacion; `r = 1` en el resto.
  3.  ORACULO DE BULTO: la misma razon con los genes del ED real puestos a
      `r = 1`. El complemento exacto del brazo 2.
  4.  SIGNOS ALEATORIOS DENSOS: la misma densidad y el mismo multiconjunto de
      |s| que la via gwps, con el signo de cada gen sorteado. 8 semillas.

SEMILLAS (todas familias nuevas, declaradas antes de medir)
  * poblacion real:        20260001 + fold          (la misma de siempre)
  * 4 signos aleatorios:   emision 20273000 + 100*fold + j,  j = 0..7
  * 4 sorteo de signos:    20274000 + j   (los MISMOS en las seis lineas)
  * 1a copia:              20273900 + fold  (la emision no consume el rng)
  * 1b oraculo completo:   20275000 + 100*fold + j,  j = 0..3
  * 2  oraculo de cola:    20276000 + 100*fold + j,  j = 0..3
  * 3  oraculo de bulto:   20277000 + 100*fold + j,  j = 0..3

🔴 NINGUNA familia colisiona con otra ni con 2026000x-2026600x, 2027100x o
2027200x/2027250x. Se comprueba en `main` y se reporta.

TODO LO DEMAS SIN TOCAR: `suelo_arnes.build_fold_adatas`, el
`cand01_covar._reescala_entero`, el bundle de competencia de cada linea, el
preset/profile `vcc2026`.
"""
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from suelo_arnes import (
    IFNB_PATH, CONTROL_LABEL, N_PER_LABEL,
    build_fold_adatas, load_global_obs, load_var_names,
)
from cand01_covar import _reescala_entero
from gwps_emision import construye_ratios, SIN_GWPS
import gwps_via as GV
import techo_arnes as TA

LINES = GV.LINES
PANEL_PATH = GV.PANEL_PATH
BUNDLES = GV.BUNDLES
WORK = os.environ.get("PDS_WORK", "/datos/modelo/out/pds_obj")
OUT_JSONL = os.environ.get("PDS_OUT", "/home/ubuntu/vcc/modelo/out/pds_objetivo.json")
ED_CACHE = os.path.join(WORK, "ed_real.json")

SEED_BASE_REAL = 20260001
SEM = {            # brazo -> (base de emision, n_semillas)
    "copia":    (20273900, 1),
    "completo": (20275000, 4),
    "cola":     (20276000, 4),
    "bulto":    (20277000, 4),
    "signos":   (20273000, 8),
}
SEED_SIGNOS = 20274000      # sorteo de signos del brazo 4, igual en las 6 lineas


def seed_emision(brazo, fold_idx, j):
    base, _ = SEM[brazo]
    return base + fold_idx + 0 if brazo == "copia" else base + 100 * fold_idx + j


def comprueba_semillas():
    """🔴 Ninguna semilla de este estudio colisiona con otra ni con las familias
    ya usadas. Se devuelve el conteo para reportarlo."""
    usadas = {}
    for brazo, (base, n) in SEM.items():
        for fold in range(len(LINES)):
            for j in range(n):
                s = seed_emision(brazo, fold, j)
                usadas.setdefault(s, []).append(f"{brazo}/fold{fold}/j{j}")
    for j in range(SEM["signos"][1]):
        usadas.setdefault(SEED_SIGNOS + j, []).append(f"signos_sorteo/j{j}")
    # 🔴 `SEED_BASE_REAL + fold` NO entra en esta comprobacion: se reusa A
    # PROPOSITO, porque la poblacion real tiene que ser la MISMA que en el suelo,
    # gwps REAL y gwps PERMUTADO para que las tablas sean comparables. Se verifica
    # por sha256 en cada corrida (`poblacion_real_identica`), no por semilla.
    # familias ya gastadas en el proyecto
    previas = set()
    for base in (20260000, 20260001, 20261000, 20262000, 20263000, 20264000,
                 20265000, 20266000):
        for fold in range(6):
            for j in range(8):
                previas.add(base + 100 * fold + j); previas.add(base + fold)
    for base in (20271000, 20272000):
        for fold in range(6):
            for j in range(8):
                previas.add(base + 100 * fold + j)
    for j in range(8):
        previas.add(20272500 + j)
    choques_internos = {s: v for s, v in usadas.items() if len(v) > 1}
    choques_previos = sorted(set(usadas) & previas)
    return {"n_semillas": len(usadas),
            "choques_internos": choques_internos,
            "choques_con_familias_previas": choques_previos}


# ---------------------------------------------------------------- razones reales

def medias_reales(real_ad, panel):
    """`(media_real_pert[h], media_control)` en CUENTAS CRUDAS, del lado real
    perturbado de la poblacion que se puntua. 🔴 Esto es lo transductivo."""
    X = real_ad.X.tocsr()
    lab = np.asarray(real_ad.obs["target_gene"].tolist(), dtype=object)
    med = {}
    for h in panel:
        m = lab == h
        assert int(m.sum()) == N_PER_LABEL, (h, int(m.sum()))
        med[h] = np.asarray(X[m].mean(axis=0)).ravel()
    mc = lab == CONTROL_LABEL
    ctrl = np.asarray(X[mc].mean(axis=0)).ravel()
    return med, ctrl, int(mc.sum())


def razones_oraculo(real_ad, panel, ed_real, modo):
    """`r_g` por perturbacion. `modo`:
      'completo' -> todos los genes;  'cola' -> solo el ED real;
      'bulto'    -> todos MENOS el ED real.
    Los genes con media de control 0 quedan en `r = 1`: `_reescala_entero` solo
    toca entradas NO NULAS, asi que ahi la razon es inalcanzable por
    construccion. Se cuenta y se reporta, no se corrige."""
    med, ctrl, n_ctrl = medias_reales(real_ad, panel)
    var = np.asarray(real_ad.var.index.values, dtype=str)
    pos = {g: i for i, g in enumerate(var)}
    vivo = ctrl > 0
    ratios, diag = {}, {}
    for h in panel:
        r = np.ones(len(var))
        with np.errstate(divide="ignore", invalid="ignore"):
            cociente = np.where(vivo, med[h] / np.where(vivo, ctrl, 1.0), 1.0)
        if modo == "completo":
            sel = vivo.copy()
        else:
            en_ed = np.zeros(len(var), dtype=bool)
            cols = [pos[g] for g in ed_real.get(h, ()) if g in pos]
            en_ed[cols] = True
            sel = (vivo & en_ed) if modo == "cola" else (vivo & ~en_ed)
        r[sel] = cociente[sel]
        ratios[h] = r
        diag[h] = {
            "n_genes_tocados": int(sel.sum()),
            "n_ed_real": int(len(ed_real.get(h, ()))),
            "n_ed_real_mapeado": int(len([g for g in ed_real.get(h, ()) if g in pos])),
            "r_min": float(r.min()), "r_max": float(r.max()),
            "n_r_mayor_1": int((r > 1).sum()), "n_r_menor_1": int((r < 1).sum()),
        }
    diag["_global"] = {
        "n_var": int(len(var)),
        "n_control_mean_cero": int((~vivo).sum()),
        "n_control_cells": n_ctrl,
        "fraccion_delta_real_en_genes_sin_control": None,  # se llena abajo
    }
    # cuanto del |delta real| vive en genes que el control no puede alcanzar
    num = den = 0.0
    for h in panel:
        d = np.abs(med[h] - ctrl)
        den += float(d.sum()); num += float(d[~vivo].sum())
    diag["_global"]["fraccion_delta_real_en_genes_sin_control"] = (
        num / den if den else None)
    return ratios, diag


def razones_signos(var_names, panel, j):
    """Brazo 4: la densidad y el multiconjunto de |s| de la via gwps, con el
    signo de cada gen sorteado. Los 6 replegados siguen en `r = 1`."""
    base, _ = construye_ratios(var_names, panel)
    rng = np.random.default_rng(SEED_SIGNOS + j)
    ratios, diag = {}, {}
    for h in panel:
        r0 = base[h]
        s = r0 - 1.0
        signo = rng.integers(0, 2, size=len(s)) * 2 - 1
        s2 = np.abs(s) * signo
        r = np.clip(1.0 + s2, 0.0, 2.0)
        ratios[h] = r
        diag[h] = {
            "replegado": h in SIN_GWPS,
            "n_tocados_original": int((s != 0).sum()),
            "n_tocados_nuevo": int((s2 != 0).sum()),
            "suma_abs_s_original": float(np.abs(s).sum()),
            "suma_abs_s_nuevo": float(np.abs(s2).sum()),
            "n_sube": int((s2 > 0).sum()), "n_baja": int((s2 < 0).sum()),
        }
    diag["_global"] = {"semilla_signos": SEED_SIGNOS + j}
    return ratios, diag


# ---------------------------------------------------------------- predict_fn

def predict_copia(real_ad, panel):
    """1a: el bloque predicho de `h` ES el bloque real perturbado de `h`."""
    X = real_ad.X.tocsr()
    lab = np.asarray(real_ad.obs["target_gene"].tolist(), dtype=object)
    bloques = {h: X[lab == h] for h in panel}

    def fn(control_pool_X, labels, rng):
        out = {}
        for label in labels:
            if label == CONTROL_LABEL:
                idx = rng.integers(0, control_pool_X.shape[0], size=N_PER_LABEL)
                out[label] = control_pool_X[idx].tocsr()
            else:
                out[label] = bloques[label].tocsr()
        return out
    return fn


def predict_razon(ratios):
    """1b, 2, 3 y 4: remuestrea el control y lo reescala por `r`, con el
    `_reescala_entero` de `cand01_covar` SIN MODIFICAR."""
    def fn(control_pool_X, labels, rng):
        n = control_pool_X.shape[0]
        out = {}
        for label in labels:
            idx = rng.integers(0, n, size=N_PER_LABEL)
            bloque = control_pool_X[idx]
            if label == CONTROL_LABEL:
                out[label] = bloque.tocsr()
            else:
                out[label] = _reescala_entero(bloque, ratios[label], rng)
        return out
    return fn


def predict_suelo():
    def fn(control_pool_X, labels, rng):
        n = control_pool_X.shape[0]
        return {l: control_pool_X[rng.integers(0, n, size=N_PER_LABEL)].tocsr()
                for l in labels}
    return fn


# ---------------------------------------------------------------- ED real

def ed_real_de_sonda(line, fold_idx, panel, cell_type, gene, var_names, n_genes):
    """Una SONDA por linea: emite el suelo, corre con `--write-degenes` y saca
    el conjunto ED REAL por perturbacion con la misma maquinaria del marcador
    (`prepare_de` + filtro p_adj + `_exclude_own_gene`), igual que
    `techo_arnes.contar_de_por_perturbacion`. El lado real depende solo de
    `seed_real`, asi que una sonda por linea alcanza."""
    import polars as pl
    from cell_eval2.de import prepare_de
    from cell_eval2.metrics.de import _exclude_own_gene

    seed_real = SEED_BASE_REAL + fold_idx
    real_ad, pred_ad, _ = build_fold_adatas(
        line, panel, predict_suelo(), var_names, seed_real, 999_000 + fold_idx,
        cell_type, gene, n_genes, IFNB_PATH)
    GV.renombra(real_ad); GV.renombra(pred_ad)
    h_real = GV.huella(real_ad)
    rp = os.path.join(WORK, f"sonda_real_{line}.h5ad")
    pp = os.path.join(WORK, f"sonda_pred_{line}.h5ad")
    real_ad.write_h5ad(rp); pred_ad.write_h5ad(pp)
    del real_ad, pred_ad
    rd = os.path.join(WORK, f"sonda_run_{line}")
    p, dt = GV.corre(GV.cmd_run(pp, rp, rd), os.path.join(WORK, f"sonda_{line}.log"))
    if p.returncode != 0:
        raise SystemExit(f"sonda {line} fallo:\n{p.stderr[-2000:]}")
    prepared = prepare_de(
        os.path.join(rd, "de_pred.parquet"), os.path.join(rd, "de_real.parquet"),
        control=GV.CONTROL_COMP, sort_by=TA.DE_SORT_BY,
        p_adj_threshold=TA.DE_P_ADJ_THRESHOLD,
        nan_lfc_policy=TA.DE_NAN_LFC_POLICY, min_abs_log2fc=TA.DE_MIN_ABS_LOG2FC)
    T = prepared.p_adj_threshold
    real_sig = _exclude_own_gene(
        prepared.real_df.filter(pl.col("p_adj") < T)
        .select("target", "feature").unique(),
        prepared, who="pds_objetivo", side="real")
    d = real_sig.to_dict(as_series=False)
    ed = {}
    for t, f in zip(d["target"], d["feature"]):
        ed.setdefault(str(t), []).append(str(f))
    for h in panel:
        ed.setdefault(h, [])
    for q in (rp, pp, os.path.join(rd, "de_real.parquet"),
              os.path.join(rd, "de_pred.parquet")):
        if os.path.exists(q):
            os.remove(q)
    return {h: sorted(v) for h, v in ed.items()}, h_real, dt


# ---------------------------------------------------------------- una corrida

def append_record(rec):
    with open(OUT_JSONL, "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def una_corrida(brazo, line, fold_idx, j, panel, cell_type, gene, var_names,
                n_genes, real_path, ref, ed_real, cache_raz):
    t0 = time.time()
    if GV.free_gib(WORK) < GV.MIN_FREE_GIB:
        raise SystemExit("PARAR: poco espacio en /datos.")
    seed_real = SEED_BASE_REAL + fold_idx
    seed_em = seed_emision(brazo, fold_idx, j)

    # el lado real primero, para poder leerlo (brazos transductivos)
    real_ad, _, meta0 = build_fold_adatas(
        line, panel, predict_suelo(), var_names, seed_real, 999_900 + fold_idx,
        cell_type, gene, n_genes, IFNB_PATH)
    diag = None
    if brazo == "copia":
        pfn = predict_copia(real_ad, panel)
    elif brazo in ("completo", "cola", "bulto"):
        if brazo not in cache_raz:
            cache_raz[brazo] = razones_oraculo(real_ad, panel, ed_real, brazo)
        ratios, diag = cache_raz[brazo]
        pfn = predict_razon(ratios)
    else:
        ratios, diag = razones_signos(var_names, panel, j)
        pfn = predict_razon(ratios)
    del real_ad

    real_ad, pred_ad, meta = build_fold_adatas(
        line, panel, pfn, var_names, seed_real, seed_em,
        cell_type, gene, n_genes, IFNB_PATH)
    GV.renombra(real_ad); GV.renombra(pred_ad)
    h_real = GV.huella(real_ad); h_pred = GV.huella(pred_ad)
    if ref is None:
        real_ad.write_h5ad(real_path)
    else:
        assert h_real == ref, f"{line}: poblacion real se movio {h_real[:16]} != {ref[:16]}"
        if not os.path.exists(real_path):
            real_ad.write_h5ad(real_path)
    del real_ad

    datos = pred_ad.X.tocsr().data
    contrato = {
        "nnz": int(pred_ad.X.nnz), "dtype": str(pred_ad.X.dtype),
        "todo_entero": bool(np.all(datos == np.round(datos))),
        "min_valor": float(datos.min()) if datos.size else None,
        "max_valor": float(datos.max()) if datos.size else None,
        "hay_negativos": bool((datos < 0).any()),
    }
    assert contrato["todo_entero"] and not contrato["hay_negativos"], contrato

    pred_path = os.path.join(WORK, f"pred_{brazo}_{line}_{j}.h5ad")
    pred_ad.write_h5ad(pred_path)
    pred_shape = list(pred_ad.shape)
    del pred_ad

    rd = os.path.join(WORK, f"run_{brazo}_{line}_{j}")
    sys.stderr.write(f"  [{brazo} {line} j={j} em={seed_em}] run... ")
    sys.stderr.flush()
    p_run, dt_run = GV.corre(GV.cmd_run(pred_path, real_path, rd),
                             os.path.join(WORK, f"run_{brazo}_{line}_{j}.log"))
    if p_run.returncode != 0:
        sys.stderr.write("FALLO\n")
        os.remove(pred_path)
        return {"brazo": brazo, "linea": line, "fold_idx": fold_idx, "j": j,
                "seed_emision": seed_em, "fallo_en": "run",
                "mensaje_literal": p_run.stderr.strip()[-2000:]}

    agg = pd.read_csv(os.path.join(rd, "agg_results.csv"))
    fila = agg.loc[agg["statistic"] == "mean"].iloc[0]
    u = {m: (None if (m not in agg.columns or pd.isna(fila[m])) else float(fila[m]))
         for m in GV.VCC2026_METRICS}
    res = pd.read_csv(os.path.join(rd, "results.csv"))
    pds_pp = {str(r.perturbation): float(r.value)
              for r in res[res.metric == "pds_cosine"].itertuples()}
    coh = GV.cohortes(rd)

    out_csv = os.path.join(WORK, f"scored_{brazo}_{line}_{j}.csv")
    p_sc, _ = GV.corre(
        GV.cmd_score(rd, os.path.join(BUNDLES, f"BUNDLE_{line}"), out_csv),
        os.path.join(WORK, f"score_{brazo}_{line}_{j}.log"))
    score_ok = p_sc.returncode == 0
    sc = GV.leer_scored(out_csv) if score_ok else None
    sys.stderr.write(f"run {dt_run:.0f}s score {'OK' if score_ok else 'SE NEGO'}  "
                     f"pds_u={u['pds_cosine']:.4f}\n")

    rec = {"brazo": brazo, "linea": line, "fold_idx": fold_idx, "j": j,
           "seed_real": seed_real, "seed_emision": seed_em,
           "huella_real": h_real, "huella_pred": h_pred,
           "poblacion_real_identica": (ref is None or h_real == ref),
           "pred_shape": pred_shape, "contrato_emision": contrato,
           "n_control_cells_pool": meta["n_control_cells"],
           "hash_run": GV.hash_del_run(rd),
           "u": u, "pds_por_pert": pds_pp, "cohortes": coh,
           "diag_emision": (None if diag is None else
                            {k: v for k, v in diag.items() if k == "_global"}),
           "score_ok": score_ok, "scored": sc,
           "run_duration_s": round(dt_run, 2),
           "total_duration_s": round(time.time() - t0, 2)}
    if not score_ok:
        rec["fallo_en"] = "score"
        rec["mensaje_literal"] = p_sc.stderr.strip()[-2000:]
    for q in (pred_path, os.path.join(rd, "de_real.parquet"),
              os.path.join(rd, "de_pred.parquet")):
        if os.path.exists(q):
            os.remove(q)
    return rec


def main():
    os.makedirs(WORK, exist_ok=True)
    chk = comprueba_semillas()
    sys.stderr.write(f"semillas: {chk['n_semillas']} distintas, "
                     f"choques internos={len(chk['choques_internos'])}, "
                     f"con familias previas={len(chk['choques_con_familias_previas'])}\n")
    assert not chk["choques_internos"] and not chk["choques_con_familias_previas"], chk
    with open(os.path.join(WORK, "semillas.json"), "w") as fh:
        json.dump(chk, fh, indent=1)

    panel = json.load(open(PANEL_PATH))["panel_comun"]
    cell_type, gene = load_global_obs(IFNB_PATH)
    var_names = load_var_names(IFNB_PATH)
    n_genes = len(var_names)

    ed = json.load(open(ED_CACHE)) if os.path.exists(ED_CACHE) else {}
    brazos = (os.environ.get("PDS_BRAZOS")
              or "copia,completo,cola,bulto,signos").split(",")
    lineas = (os.environ.get("PDS_LINEAS") or ",".join(LINES)).split(",")

    hechas = {}
    if os.path.exists(OUT_JSONL):
        for ln in open(OUT_JSONL):
            if ln.strip():
                r = json.loads(ln)
                if not r.get("fallo_en"):
                    hechas[(r["brazo"], r["linea"], r["j"])] = r["huella_real"]

    for line in lineas:
        fold_idx = LINES.index(line)
        if line not in ed:
            sys.stderr.write(f"\n=== sonda ED real {line} ===\n")
            ed[line], h, dt = ed_real_de_sonda(
                line, fold_idx, panel, cell_type, gene, var_names, n_genes)
            ed[line + "__huella_real"] = h
            with open(ED_CACHE, "w") as fh:
                json.dump(ed, fh)
            sys.stderr.write(f"  ED real: {sum(len(v) for v in ed[line].values())} "
                             f"pares (target, gen) en {dt:.0f}s\n")
        real_path = os.path.join(WORK, f"real_{line}.h5ad")
        ref = ed.get(line + "__huella_real")
        cache_raz = {}
        for brazo in brazos:
            n = SEM[brazo][1]
            if os.environ.get("PDS_NMAX"):      # solo para la prueba de humo
                n = min(n, int(os.environ["PDS_NMAX"]))
            sys.stderr.write(f"\n=== {line} (fold {fold_idx}) brazo {brazo} ===\n")
            for j in range(n):
                if (brazo, line, j) in hechas:
                    sys.stderr.write(f"  [{brazo} {line} j={j}] ya hecha\n")
                    continue
                rec = una_corrida(brazo, line, fold_idx, j, panel, cell_type,
                                  gene, var_names, n_genes, real_path, ref,
                                  ed[line], cache_raz)
                append_record(rec)
                if rec.get("fallo_en"):
                    sys.stderr.write(f"\n*** PARAR: {rec['fallo_en']} "
                                     f"{brazo}/{line}/j={j}\n"
                                     f"{rec['mensaje_literal']}\n")
                    sys.exit(7)
                ref = ref or rec["huella_real"]
        if os.path.exists(real_path):
            os.remove(real_path)
    sys.stderr.write(f"\nlisto: {OUT_JSONL}\n")


if __name__ == "__main__":
    main()
