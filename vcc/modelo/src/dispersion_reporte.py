#!/usr/bin/env python
"""Escribe modelo/out/dispersion_u.md desde dispersion_u.json.

Lee y cuenta. No mide nada nuevo y no propone ningun margen ni umbral.
"""
import json
import os
import sys

import numpy as np

OUT_MD = "/home/ubuntu/vcc/modelo/out/dispersion_u.md"
DISP = "/home/ubuntu/vcc/modelo/out/dispersion_u.json"
NATIVAS = "/home/ubuntu/vcc/modelo/out/anclas_nativas_comp.json"

LINEAS = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
# nombre canonico -> etiqueta corta, en el orden de bucle_suelo.md
METRICAS = [
    ("pds_cosine", "pds_cosine"),
    ("expr_mse_unbiased_capped_norm", "mse_capped_norm"),
    ("de_wilcoxon_lfc_nmae", "lfc_nmae"),
    ("de_wilcoxon_direction_fidelity_yield_raw", "fid"),
    ("de_wilcoxon_direction_reach_raw", "reach"),
    ("de_wilcoxon_sig_jaccard", "jac"),
]
SENTIDO = {"pds_cosine": "alto", "expr_mse_unbiased_capped_norm": "bajo",
           "de_wilcoxon_lfc_nmae": "bajo",
           "de_wilcoxon_direction_fidelity_yield_raw": "alto",
           "de_wilcoxon_direction_reach_raw": "alto",
           "de_wilcoxon_sig_jaccard": "alto"}
N_SEM = 8


def leer_jsonl(path):
    out = []
    with open(path) as fh:
        for ln in fh:
            if ln.strip():
                out.append(json.loads(ln))
    return out


def f(x, n=6):
    if x is None:
        return "n/d"
    if isinstance(x, float) and np.isnan(x):
        return "NaN"
    return f"{x:.{n}f}"


def tabla(L, head, rows):
    L.append("| " + " | ".join(head) + " |")
    L.append("|" + "|".join(["---"] * len(head)) + "|")
    for r in rows:
        L.append("| " + " | ".join(str(c) for c in r) + " |")
    L.append("")


def main():
    recs = leer_jsonl(DISP)
    nat = {r["linea"]: r for r in leer_jsonl(NATIVAS)}

    # indice: (linea, j) -> registro
    idx = {(r["linea"], r["j"]): r for r in recs if not r.get("fallo_en")}
    fallos = [r for r in recs if r.get("fallo_en")]

    L = []
    L += [
        "# La dispersion de `u`",
        "",
        "Cuanto se mueve el `s` de una MISMA emision al cambiar **solo** la semilla de",
        "emision. 8 semillas x 6 lineas = 48 corridas `run` -> `score --real-bundle`",
        "contra los bundles de competencia de `/datos/modelo/out/nativas_comp/`.",
        "",
        "La emision es `suelo_floor.predict_floor` sin modificar -- la misma desechable de",
        "siempre, no un candidato. Las 8 corridas de una linea **no se promedian**: se",
        "puntuan por separado y lo que se reporta es su dispersion.",
        "",
    ]

    # ------------------------------------------------ 1. lo que vario y lo que no
    L += [
        "## 1. Lo que vario y lo que no, con la prueba",
        "",
        "**Vario solo la semilla de emision.** Familia nueva `20264000 + 100*fold + j`,",
        "`j = 0..7`. Ninguna colision con `2026000x` (verdad), `2026100x` (suelo),",
        "`2026200x`/`2026300x` (techo) ni la `0` del anchor.",
        "",
        "**No vario la poblacion real.** `build_fold_adatas` sin modificar con",
        "`seed_real = 20260001 + fold_idx`, etiqueta `non-targeting`, identica a",
        "`cable_suelo.md`. No se asume: de cada construccion se tomo un `sha256` de la",
        "matriz real (`data`, `indices`, `indptr`) y de sus etiquetas. Las 8 de cada linea",
        "tienen que coincidir, y el script aborta si no.",
        "",
    ]
    rows = []
    for ln in LINEAS:
        rr = [idx[(ln, j)] for j in range(N_SEM) if (ln, j) in idx]
        if not rr:
            rows.append([ln, "**faltan**"] + ["-"] * 5)
            continue
        hs = {r["huella_real"] for r in rr}
        hp = {r["huella_pred"] for r in rr}
        hashes_ok = all((r.get("hash_run") or {}).get("coincide") for r in rr)
        sc_ok = all(r.get("score_ok") for r in rr)
        rows.append([
            ln, len(rr),
            f"{20264000 + 100*rr[0]['fold_idx']}..{20264000 + 100*rr[0]['fold_idx'] + 7}",
            f"**{len(hs)}**", f"**{len(hp)}**",
            "**si**" if hashes_ok else "**NO**",
            "**si**" if sc_ok else "**NO**",
        ])
    tabla(L, ["linea", "corridas", "semillas de emision",
              "huellas distintas del REAL", "huellas distintas del PRED",
              "config_hash de competencia en todas", "score corrio en todas"], rows)
    L += [
        "**Una sola huella del real por linea y ocho del predicho**: la poblacion esta",
        "fija y la emision cambia en las ocho. Es la condicion que hace que la dispersion",
        "medida tenga una sola fuente.",
        "",
        "El `config_hash` recalculado desde el `run_params.yaml` de cada corrida coincide",
        "con el de referencia `a5519017556705ab...` en las 48. `--write-degenes` -- que se",
        "pasa para poder MEDIR las dos cohortes moviles en vez de inferirlas -- no lo",
        "desvia.",
        "",
    ]
    if fallos:
        L += [f"⚠️ **{len(fallos)} corridas fallaron.** Mensajes literales:", ""]
        for r in fallos:
            L += [f"- {r['linea']} j={r['j']} en `{r['fallo_en']}`:",
                  "  ```", "  " + (r.get("mensaje_literal") or "")[:600], "  ```"]
        L.append("")
    else:
        L += ["Ninguna corrida fallo: 48 de 48 `run` y 48 de 48 "
              "`score --real-bundle`.", ""]

    # ------------------------------------------------ 2. la formula del margen
    L += [
        "## 2. El margen combinado: la formula, y por que esa",
        "",
        "La escala es `s = (u - b)/(r - b)` (`score.py:395-399`). De las tres cantidades,",
        "esta medicion mueve **una**:",
        "",
        "| cantidad | que es | su ruido aqui |",
        "|---|---|---|",
        "| `u` | la emision del candidato | **lo que mide este encargo**: 8 semillas, `sd_u` |",
        "| `r` | el ancla replicada, media de 5 particiones | `SE_r = sd_r/sqrt(5)`, de `anclas_nativas_v2.md` |",
        "| `b` | el baseline generico, **una sola realizacion** con `seed 0` | no estimado: ver la reserva de abajo |",
        "",
        "Propagacion de primer orden. Con `b` fijo y `u`, `r` independientes:",
        "",
        "```",
        "ds/du = 1/(r - b)",
        "ds/dr = -(u - b)/(r - b)^2 = -s/(r - b)",
        "",
        "sd(s) de la emision   = sd_u / |r - b|          <- ya viene medida directa en s",
        "sd(s) del ancla       = |s| * SE_r / |r - b|",
        "",
        "margen_celda = sqrt( sd_s^2 + ( |s| * SE_r / |r - b| )^2 )",
        "```",
        "",
        "Dos detalles que importan y que no se aplican en silencio:",
        "",
        "1. **`sd_s` se mide sobre los `s`, no sobre los `u`.** Las 8 corridas se puntuaron",
        "   por separado contra el mismo bundle, asi que la division por `(r - b)` ya esta",
        "   dentro de cada `s`. Dividir otra vez seria contarla dos veces. El `sd_u` crudo",
        "   se reporta igual, en su columna, porque es la cantidad que no depende del",
        "   bundle.",
        "2. **El termino del ancla escala con `|s|`.** En estas 48 corridas los `s` son ~0 o",
        "   negativos y pequenos, asi que ese termino es pequeno **por el valor de `s`, no",
        "   porque `r` sea estable**. Para un candidato con `s` grande pesaria mas.",
        "",
        "**La reserva sobre `b`.** `b` es una medicion, no una constante: `cell-eval2",
        "baseline --seed 0` sobre la poblacion evaluada. Su propia dispersion no se estimo",
        "-- habria hecho falta correr el baseline con varias semillas, que no es este",
        "encargo -- asi que el margen de arriba **omite un termino**. Queda declarado: el",
        "margen es una cota INFERIOR del ruido total de `s`.",
        "",
        "**Y el caso que no necesita el termino del ancla.** Para comparar DOS candidatos en",
        "la MISMA celda del MISMO bundle, `r` y `b` son comunes y se cancelan:",
        "",
        "```",
        "s1 - s2 = (u1 - u2)/(r - b)",
        "```",
        "",
        "El ruido del ancla no entra en la diferencia: mueve las dos `s` juntas. Ahi la",
        "cantidad relevante es solo `sd_s`. Las dos columnas estan en la tabla para que se",
        "pueda usar la que corresponda; **no se propone ningun umbral**.",
        "",
        "**Convencion de `sd`, declarada.** `sd_r` se toma tal como lo estampa",
        "`anclas_nativas_v2.md`, que uso `ddof=0` (verificado: su 0.035008 para A549/jac se",
        "reproduce con `ddof=0` y no con `ddof=1`, que da 0.039140). Para `u` se reporta",
        "`ddof=1` como principal -- es el estimador no sesgado de la dispersion de UNA",
        "emision, que es la pregunta -- y `ddof=0` al lado. Con 8 draws la razon entre los",
        "dos es `sqrt(8/7) = 1.069`.",
        "",
    ]

    # ------------------------------------------------ 3. las 36 celdas
    def serie(ln, met):
        vals = []
        for j in range(N_SEM):
            r = idx.get((ln, j))
            if r is None:
                vals.append(None)
                continue
            sc = (r.get("scored") or {}).get("filas", {}).get(met, {})
            vals.append(sc.get("from_replicate"))
        return vals

    def serie_u(ln, met):
        return [(idx[(ln, j)]["u"].get(met) if (ln, j) in idx else None)
                for j in range(N_SEM)]

    def sd_r_de(ln, met):
        sp = nat[ln]["splits"].get(met)
        if not sp:
            return None, None
        v = [s["value"] for s in sp]
        return float(np.std(v, ddof=0)), v

    celdas = []
    for ln in LINEAS:
        for met, corto in METRICAS:
            s = [x for x in serie(ln, met) if x is not None]
            uu = [x for x in serie_u(ln, met) if x is not None]
            b = nat[ln]["b"].get(met)
            rr = nat[ln]["r"].get(met)
            r = (rr.get("replicate") if isinstance(rr, dict) else rr)
            sd_r, splits_r = sd_r_de(ln, met)
            se_r = (sd_r / np.sqrt(5)) if sd_r is not None else None
            span = (r - b) if (r is not None and b is not None) else None
            media = float(np.mean(s)) if s else None
            sd1 = float(np.std(s, ddof=1)) if len(s) > 1 else None
            sd0 = float(np.std(s, ddof=0)) if s else None
            sdu1 = float(np.std(uu, ddof=1)) if len(uu) > 1 else None
            term_r = (abs(media) * se_r / abs(span)
                      if (media is not None and se_r is not None and span) else None)
            margen = (float(np.hypot(sd1, term_r))
                      if (sd1 is not None and term_r is not None) else None)
            celdas.append({
                "linea": ln, "metrica": corto, "canon": met, "s": s, "u": uu,
                "media": media, "sd1": sd1, "sd0": sd0, "sdu1": sdu1,
                "rango": (max(s) - min(s)) if s else None,
                "min": min(s) if s else None, "max": max(s) if s else None,
                "b": b, "r": r, "span": span, "sd_r": sd_r, "se_r": se_r,
                "term_r": term_r, "margen": margen,
                "splits_r": splits_r,
                "u_media": float(np.mean(uu)) if uu else None,
            })

    L += [
        "## 3. Las 36 celdas",
        "",
        "Antes de desglosarlas, el agregado, porque es el numero que decide y sale gratis",
        "de las mismas 48 corridas. `avg_score` es el promedio plano de las seis `s`",
        "(`from_replicate`; `from_baseline` queda null con un bundle de competencia):",
        "",
    ]
    rows = []
    for ln in LINEAS:
        v = [((idx[(ln, j)].get("scored") or {}).get("filas", {})
              .get("avg_score", {}) or {}).get("from_replicate")
             for j in range(N_SEM) if (ln, j) in idx]
        v = [x for x in v if x is not None]
        if not v:
            continue
        rows.append([
            ln, " / ".join(f"{x:+.4f}" for x in v), f(float(np.mean(v))),
            f(float(np.std(v, ddof=1))), f(max(v) - min(v)),
            f(min(v), 4), f(max(v), 4),
        ])
    tabla(L, ["linea", "los 8 avg_score", "media", "sd (ddof=1)", "rango",
              "min", "max"], rows)
    if rows:
        rangos = [float(r[4]) for r in rows]
        sds = [float(r[3]) for r in rows]
        L += [
            f"**Solo cambiando la semilla de emision, el `avg_score` de una misma emision",
            f"se mueve en un rango de {min(rangos):.4f} a {max(rangos):.4f}** segun la",
            f"linea, con `sd` de {min(sds):.4f} a {max(sds):.4f}. Es la cifra que faltaba",
            "para poder separar una diferencia entre candidatos del sorteo.",
            "",
            "Y conviene decir de donde NO sale: las seis `s` que promedia no tienen la misma",
            "dispersion -- hay celdas con `sd = 0` exacta por razones distintas (§5) -- asi",
            "que todo el movimiento del agregado lo ponen las celdas que si se mueven,",
            "repartido a 1/6 cada una.",
            "",
            "⚠️ **El caso de A549 no es ruido continuo.** Su rango (0.2476) y su `sd` (0.0792)",
            "son los mayores de las seis, y 0.270 de ese rango los pone una sola celda en una",
            "sola semilla: `j=5` perdio el 1.0 de convencion de **STAT5A** en `jac` y cobro",
            "`s = -1.6197` en vez de 0, que a 1/6 son -0.270 de agregado. Las otras siete",
            "semillas caen entre -0.0714 y -0.1831. §5.2 lo desarma entero. Leer la `sd` de",
            "A549 como la anchura de una nube seria leer mal la celda.",
            "",
        ]
    L += [
        "Los 8 valores de `s` de cada celda, con su dispersion. `s` es la columna",
        "`from_replicate` de `scored.csv`, que con un bundle de competencia es la que toma",
        "`avg_score`.",
        "",
    ]
    rows = []
    for c in celdas:
        rows.append([
            c["linea"], f"`{c['metrica']}`",
            " / ".join(f"{x:+.4f}" for x in c["s"]),
            f(c["media"]), f(c["sd1"]), f(c["sd0"]), f(c["rango"]),
            f(c["min"], 4), f(c["max"], 4),
        ])
    tabla(L, ["linea", "metrica", "los 8 valores de s", "media", "sd (ddof=1)",
              "sd (ddof=0)", "rango", "min", "max"], rows)

    L += [
        "Y el margen combinado de cada celda, con sus piezas a la vista:",
        "",
    ]
    rows = []
    for c in celdas:
        rows.append([
            c["linea"], f"`{c['metrica']}`", SENTIDO[c["canon"]],
            f(c["u_media"]), f(c["sdu1"]), f(c["b"]), f(c["r"]),
            f(c["span"]), f(c["sd_r"]), f(c["se_r"]),
            f(c["media"]), f(c["sd1"]), f(c["term_r"]),
            f"**{f(c['margen'])}**",
        ])
    tabla(L, ["linea", "metrica", "sentido", "u medio", "sd_u (ddof=1)", "b", "r",
              "span r-b", "sd_r (5 splits)", "SE_r", "s medio", "sd_s",
              "termino de r", "margen combinado"], rows)

    # ------------------------------------------------ 4. sd acumulada
    L += [
        "## 4. La `sd` acumulada con k = 2..8",
        "",
        "La `sd` (ddof=1) de las primeras `k` semillas, en el orden `j = 0, 1, ... 7` --",
        "que es el orden de la familia de semillas, fijado antes de medir y no una",
        "reordenacion por valor. Responde a cuantas semillas hacen falta: si la columna se",
        "aplana en k = 4, pedir 8 es gasto.",
        "",
    ]
    rows = []
    vivos = []
    for c in celdas:
        s = c["s"]
        acc = []
        for k in range(2, N_SEM + 1):
            acc.append(float(np.std(s[:k], ddof=1)) if len(s) >= k else None)
        c["acum"] = acc
        est = None
        if acc[-1] and acc[-1] > 0:
            # k mas chico desde el cual la sd ya esta dentro del 20% de la de 8
            for k, v in zip(range(2, N_SEM + 1), acc):
                if v is not None and abs(v - acc[-1]) <= 0.20 * acc[-1]:
                    est = k
                    break
            vivos.append((c, est))
        rows.append([
            c["linea"], f"`{c['metrica']}`",
            *[("0" if (v == 0) else f(v)) if v is not None else "n/d" for v in acc],
            str(est) if est else ("--" if not acc[-1] else "n/d"),
        ])
    tabla(L, ["linea", "metrica"] + [f"k={k}" for k in range(2, N_SEM + 1)]
          + ["primer k dentro del 20% de k=8"], rows)
    if vivos:
        ks = [k for _, k in vivos if k]
        L += [
            f"De las **{len(vivos)} celdas con `sd > 0`**, el primer `k` cuya `sd` ya cae",
            f"dentro del 20% de la de `k = 8` es: "
            + ", ".join(f"{c['linea']}/`{c['metrica']}` k={k}" for c, k in vivos if k),
            ".",
            "",
            f"Mediana de ese `k`: **{int(np.median(ks))}**; maximo **{max(ks)}**. Es una",
            "descripcion de estas 48 corridas, no una recomendacion: el encargo mide la",
            "dispersion y no propone cuantas semillas pedir.",
            "",
            "⚠️ **Y una celda donde esa columna no significa lo que parece.** Donde la `sd`",
            "arranca en 0 y salta de golpe (A549/`jac` 0 hasta k=5 y 0.661 en k=6; HAP1/`fid`",
            "0 hasta k=5; HT29/`fid` 0 hasta k=5; K562/`fid` y K562/`jac` 0 hasta k=3) la",
            "dispersion no es una nube que se estreche con mas draws: es un evento discreto",
            "que ocurre o no (§5.2). Ahi `k` no dice cuando se estabiliza la `sd`, dice en que",
            "draw cayo el evento. La `sd` acumulada de una Bernoulli con 8 draws no converge",
            "en 4; estima una probabilidad con 8 observaciones.",
            "",
        ]

    # ------------------------------------------------ 5. las sd = 0
    L += [
        "## 5. Las `sd = 0`, separadas por causa",
        "",
        "Son las celdas donde la emision no mueve nada. Tienen causas distintas y una de",
        "ellas **no es estabilidad**.",
        "",
    ]
    ceros = [c for c in celdas if c["sd1"] == 0.0]
    novale = [c for c in celdas if c["sd1"] is not None and c["sd1"] > 0]
    L += [
        f"**{len(ceros)} de 36 celdas tienen `sd = 0` exacta**, y "
        f"{len(novale)} tienen dispersion.",
        "",
    ]
    L += [
        "**La particion es limpia y la decide una sola columna: `sd_u`.**",
        "",
        "- `sd_u = 0`: la emision no mueve el `u` crudo. La `sd_s = 0` es **estructural**",
        "  y es informacion real -- esa celda no tiene ruido de semilla.",
        "- `sd_u > 0` con `sd_s = 0`: el `u` crudo SI se mueve y el `s` no. Algo entre los",
        "  dos borra la variacion. Eso **no es estabilidad**.",
        "",
        "Antes de la tabla, el aviso que cambia como se lee: el encargo esperaba que `jac`",
        "diera `sd = 0` por razon estructural **en las seis lineas**, y solo la da en dos.",
        "Las otras cuatro estan entre las celdas con dispersion, por un mecanismo discreto",
        "que §5.2 nombra perturbacion por perturbacion.",
        "",
    ]
    rows = []
    for c in ceros:
        sdu = c["sdu1"]
        if sdu is not None and sdu > 0:
            if c["canon"] == "expr_mse_unbiased_capped_norm":
                causa = "**clamp a [0,1]** (`catalog.py:415-416`) -- ceguera, no estabilidad"
            else:
                causa = "**sin identificar** -- `sd_u > 0` y `sd_s = 0`"
        elif c["canon"] == "de_wilcoxon_sig_jaccard":
            causa = ("**estructural** -- `E = 0` en las 8: el lado real tiene senal en las 33, "
                     "no hay 1.0 de convencion que perder (§5.2)")
        elif c["canon"] == "de_wilcoxon_direction_fidelity_yield_raw":
            causa = ("**estructural** -- `u` constante en "
                     + (f"{f(c['u_media'],4)}" if c["u_media"] is not None else "n/d"))
        else:
            causa = "**estructural** -- `u` constante; mecanismo no identificado"
        rows.append([
            c["linea"], f"`{c['metrica']}`", f(c["u_media"]), f(c["sdu1"]),
            f(c["b"]), f(c["media"]), causa,
        ])
    tabla(L, ["linea", "metrica", "u medio", "**sd_u**", "b", "s", "causa"], rows)
    n_clamp = sum(1 for r in rows if "clamp" in r[6])
    n_estr = sum(1 for r in rows if "estructural" in r[6])
    n_ni = sum(1 for r in rows if "sin identificar" in r[6])
    L += [
        f"Reparto: **{n_estr} estructurales**, **{n_clamp} por el clamp** y "
        f"**{n_ni} sin identificar**.",
        "",
    ]

    # 5.1 mse
    mses = [c for c in celdas if c["canon"] == "expr_mse_unbiased_capped_norm"]
    todos_cero = all(c["sd1"] == 0.0 and c["media"] == 0.0 for c in mses)
    L += [
        "### 5.1 `mse_capped_norm`: ceguera, no estabilidad",
        "",
        ("**Las seis lineas dan `s = 0.000000` exacto en las 48 corridas.**"
         if todos_cero else
         "⚠️ **No todas las lineas dan 0.** Ver la tabla."),
        "",
    ]
    rows = []
    for c in mses:
        rows.append([
            c["linea"], f(c["u_media"]), f(c["b"]), f(c["r"]),
            f"si" if (c["u_media"] is not None and c["b"] is not None
                      and ((c["u_media"] > c["b"]) if SENTIDO[c["canon"]] == "bajo"
                           else (c["u_media"] < c["b"]))) else "no",
            " / ".join(f"{x:.6f}" for x in c["u"]),
            f(float(np.std(c["u"], ddof=1))) if len(c["u"]) > 1 else "n/d",
            f(c["media"]), f(c["sd1"]),
        ])
    tabla(L, ["linea", "u medio", "b", "r", "u peor que b",
              "los 8 valores de u (CRUDOS)", "sd de u", "s", "sd de s"], rows)
    L += [
        "La columna que lo dice todo es `sd de u` contra `sd de s`: **el `u` crudo SI se",
        "mueve de semilla a semilla**, y el `s` no se mueve nada. El clamp lo borra.",
        "",
        "`catalog.py:415-416` fija `[0, 1]` sobre esta metrica y solo sobre esta. El suelo",
        "esta por debajo de `b` en las seis lineas (columna `u peor que b`), asi que `s`",
        "sale negativo antes del clamp y el clamp lo devuelve exactamente 0. **Una celda",
        "que marca 0 con `sd = 0` aqui no informa de que un candidato este mejorando: dos",
        "candidatos distintos, los dos por debajo de `b`, son indistinguibles.** Es la",
        "propiedad que `descomposicion.md` ya habia declarado, ahora medida con 8 semillas:",
        "la dispersion del `u` crudo es la prueba de que hay senal que el clamp tira.",
        "",
        "### 5.2 `jac`: la prediccion del encargo sale REFUTADA en 4 de 6 lineas",
        "",
        "El encargo predijo que la `E` del suelo no depende de la semilla de emision,",
        "porque es el complemento de la cohorte del lado real (§3.15), y que por eso la",
        "`sd` de `jac` deberia ser 0 por razon estructural, **y pidio comprobarlo y, si NO",
        "era 0, decir que de la emision mueve la union.** No es 0:",
        "",
        "- `sd = 0` solo en **BXPC3 y MCF7**, y ahi es estructural por una razon mas fuerte",
        "  que la predicha: `E = 0` en las 8 semillas. El lado real tiene senal DE en las 33",
        "  perturbaciones, asi que no hay ningun 1.0 de convencion que perder.",
        "- En **A549, HAP1, HT29 y K562 la `E` SI se mueve con la semilla de emision**, y con",
        "  ella `jac`. `sd_u` = 0.010714, 0.015683, 0.010714, 0.010714. La emision mueve la",
        "  union.",
        "",
        "**Que la mueve, nombrado.** La emision predice algun gen significativo en "
        + str(sum(1 for r in recs if r["cohortes"]["perts_con_pred_sig"]))
        + " de las",
        "48 corridas. Eso por si solo no cambia `E`: solo la cambia cuando la perturbacion",
        "predicha es una cuyo lado real esta VACIO. Entonces esa perturbacion pierde su 1.0",
        "gratis y cobra un jaccard genuino de 0 (interseccion vacia con un predicho no",
        "vacio). Las ocurrencias, una por una:",
        "",
        "| linea | semilla | perturbacion que salio de la union vacia | E antes -> despues | mecanismo |",
        "|---|---|---|---|---|",
        "| A549 | j=5 | **STAT5A** | 15 -> 14 | predicho significativo sobre real vacio |",
        "| HAP1 | j=0 | **RAP1GAP** | 14 -> 13 | predicho significativo sobre real vacio |",
        "| HAP1 | j=5 | **BRD9** | 14 -> 13 | predicho significativo sobre real vacio |",
        "| HAP1 | j=7 | **SP100** | 14 -> 13 | predicho significativo sobre real vacio |",
        "| HT29 | j=1 | **RAP1GAP** | 6 -> 5 | predicho significativo sobre real vacio |",
        "| K562 | j=3 | **SOCS1** | 8 -> 7 | predicho significativo sobre real vacio |",
        "",
        "Son los mismos nombres que `huecos.md` §2.5 habia aislado para el `b` (BRD9, IRF7,",
        "RAP1GAP, ZBP1): el conjunto de perturbaciones donde esto puede pasar es el",
        "complemento de la cohorte del lado real, y es el lado real el que lo fija. Lo que",
        "cambia con la semilla es CUAL de ellas cae, no cuales son candidatas.",
        "",
        "**Y la amplificacion, que es el hallazgo caro.** Un cambio de `E` de uno mueve el",
        "`jac` crudo en exactamente `1/33 = 0.030303`. Dividido por el span de la celda, eso",
        "es el salto en `s`:",
        "",
        "| linea | span `r - b` | `(1/33)/span` | salto de `s` observado | coincide |",
        "|---|---|---|---|---|",
        "| A549 | 0.018709 | **1.619745** | 0.0000 -> -1.619745 | exacto |",
        "| HAP1 | 0.043914 | **0.690056** | 0.0000 -> -0.690056 | exacto |",
        "| HT29 | 0.181026 | **0.167396** | 0.501189 -> 0.333793 | exacto |",
        "| K562 | 0.277874 | **0.109053** | 0.109053 -> 0.0000 | exacto |",
        "",
        "En A549 el span entero (0.0187) es **menos que una unidad de convencion** (0.0303):",
        "una sola perturbacion de 33 que cruce el umbral de significancia mueve `s` en 1.62",
        "puntos. Repartido a 1/6 en el agregado son 0.270 de `avg_score` -- y el rango de",
        "`avg_score` de A549 es 0.2476, el mayor de las seis lineas. **Casi toda la",
        "dispersion del agregado de A549 es ese unico evento discreto**, no ruido continuo de",
        "la emision: `j=5` da -0.3190 y las otras siete caen entre -0.0714 y -0.1831.",
        "",
        "Consecuencia para un candidato, dicha sin proponer nada: en las celdas de `jac` con",
        "span pequeno la `sd` no describe una nube, describe una moneda. Promediar 8 semillas",
        "ahi no estrecha un intervalo gaussiano; estima la probabilidad de un salto.",
        "",
        "Conviene separar dos cosas que no son la misma: que `sd = 0` **no** quiere decir",
        "`s = 0`. En BXPC3 y MCF7, las dos con `sd = 0`, el `s` sale -0.004893 y -0.000399",
        "-- practicamente 0 porque `u = 0` y el `b` tambien es casi 0. En HT29 y K562 el `u`",
        "esta POR ENCIMA del `b` -- las cuatro celdas de `huecos.md` §2.5 -- y `s` sale",
        "positivo (+0.4803 y +0.0954 de media): un puntaje que una copia del control cobra",
        "por convencion, aunque ahi la semilla si lo mueva.",
        "",
    ]
    rows = []
    for ln in LINEAS:
        Es = [idx[(ln, j)]["cohortes"]["E_union_vacia"] for j in range(N_SEM)
              if (ln, j) in idx]
        uns = [idx[(ln, j)]["cohortes"]["n_union_no_vacia"] for j in range(N_SEM)
               if (ln, j) in idx]
        conj = [tuple(idx[(ln, j)]["cohortes"]["perts_union_vacia"])
                for j in range(N_SEM) if (ln, j) in idx]
        cj = [c for c in celdas if c["linea"] == ln
              and c["canon"] == "de_wilcoxon_sig_jaccard"][0]
        rows.append([
            ln, "/".join(map(str, Es)), len(set(Es)),
            "/".join(map(str, uns)), len(set(conj)),
            f(cj["u_media"]), f(cj["b"]), f(cj["sdu1"]), f(cj["sd1"]),
        ])
    tabla(L, ["linea", "E por semilla (j=0..7)", "valores distintos de E",
              "union no vacia por semilla", "CONJUNTOS distintos de union vacia",
              "u medio de jac", "b de jac", "sd de u", "sd de s"], rows)
    L += [
        "La columna `CONJUNTOS distintos` es la que hace la comprobacion: donde vale 1",
        "(BXPC3, MCF7) son exactamente las mismas perturbaciones en las 8 semillas; donde",
        "vale 2 o 4 (A549, HT29, K562 y HAP1) **no lo son**, y el conteo coincide con el",
        "numero de configuraciones distintas de la tabla de ocurrencias de arriba.",
        "",
        "La parte de la prediccion que SI se sostiene es la mas de fondo: el CONJUNTO de",
        "perturbaciones donde `E` puede moverse es el complemento de la cohorte del lado",
        "real, y eso no depende de la emision. `huecos.md` §2.5 midio para nuestro suelo la",
        "igualdad `E_u = 33 - (tamano de la cohorte del lado real)`; lo que aqui se corrige",
        "es que esa igualdad se cumple **por semilla pero no identicamente en todas**: el",
        "suelo deja de ser una copia exacta del `b` en 6 de las 48 corridas.",
        "",
    ]

    # ------------------------------------------------ 6. cohortes moviles
    L += [
        "## 6. Las dos cohortes moviles, por semilla",
        "",
        "`fid` es la unica de las seis cuyo gate puede moverlo la submission: su",
        "denominador es `max(n_pred, N_conf)` (`direction.py:776`) y `_safe_div` devuelve",
        "NaN cuando vale 0 (`direction.py:689-699`), asi que su cohorte es",
        "`max(n_pred, N_conf) > 0`. Las otras cinco dependen del lado real solo.",
        "",
    ]
    rows = []
    for ln in LINEAS:
        rr = [idx[(ln, j)] for j in range(N_SEM) if (ln, j) in idx]
        nf = [r["cohortes"]["n_fid"] for r in rr]
        nr = [r["cohortes"]["n_reach"] for r in rr]
        nl = [r["cohortes"]["n_lfc"] for r in rr]
        nu = [r["cohortes"]["n_union_no_vacia"] for r in rr]
        npred = [r["cohortes"]["n_perts_con_pred_sig"] for r in rr]
        spred = [r["cohortes"]["suma_n_DE_pred"] for r in rr]
        rows.append([
            ln, "/".join(map(str, nf)), f"{min(nf)}-{max(nf)}" if nf else "-",
            "/".join(map(str, nu)),
            "/".join(map(str, nr)), "/".join(map(str, nl)),
            "/".join(map(str, npred)), "/".join(map(str, spred)),
        ])
    tabla(L, ["linea", "cohorte de fid por semilla", "rango fid",
              "union no vacia (jac) por semilla", "cohorte de reach",
              "cohorte de lfc_nmae", "perts con algun gen predicho significativo",
              "suma de n_DE_pred"], rows)
    movidas = []
    for ln in LINEAS:
        rr = [idx[(ln, j)] for j in range(N_SEM) if (ln, j) in idx]
        nf = [r["cohortes"]["n_fid"] for r in rr]
        if len(set(nf)) > 1:
            movidas.append((ln, min(nf), max(nf)))
    L += [
        ("**La cohorte de `fid` SI se mueve sola, solo por la semilla de emision**, en "
         + ", ".join(f"{ln} ({a}-{b})" for ln, a, b in movidas) + "."
         if movidas else
         "**Ninguna cohorte se movio entre semillas en ninguna linea.**"),
        "",
        "Las de `reach` y `lfc_nmae` no se mueven en ninguna linea, como tiene que ser:",
        "son del lado real y el lado real esta fijo (§1).",
        "",
    ]
    if movidas:
        L += [
            "Esto es ruido del candidato, no del ancla: la misma emision con otra semilla",
            "mete o saca perturbaciones del promedio de `fid`. Hay que verlo en cada",
            "corrida, que es lo que §3.16 pide.",
            "",
        ]

    # ------------------------------------------------ 7. 20 filas trabajadas
    L += [
        "## 7. Las 20 filas trabajadas",
        "",
        "**Criterio de seleccion, declarado antes de mirar los valores**: las celdas de",
        "`pds_cosine`, `lfc_nmae` y `reach` en las seis lineas = 18 filas, mas las dos de",
        "`jac` de HT29 y K562 -- las dos lineas donde `huecos.md` §2.5 ya habia identificado,",
        "ANTES de este encargo, que la `E` del `b` no coincide con la de nuestro suelo",
        "(BRD9/IRF7/RAP1GAP y ZBP1) = **20**. Es un criterio por metrica y por hallazgo",
        "previo, no por el valor medido aqui.",
        "",
        "Cada fila trae los 8 `u` crudos, los 8 `s`, y la aritmetica `(u - b)/(r - b)`",
        "rehecha a mano sobre la media para que se vea que el `s` reportado es ese.",
        "",
    ]
    sel = [c for c in celdas if c["metrica"] in ("pds_cosine", "lfc_nmae", "reach")]
    sel += [c for c in celdas if c["metrica"] == "jac"
            and c["linea"] in ("HT29", "K562")]
    rows = []
    for c in sel:
        s_mano = ((c["u_media"] - c["b"]) / c["span"]) if c["span"] else None
        rows.append([
            c["linea"], f"`{c['metrica']}`",
            " / ".join(f"{x:.4f}" for x in c["u"]),
            f(c["u_media"]), f(c["sdu1"]),
            f(c["b"]), f(c["r"]), f(c["span"]),
            f(s_mano) if s_mano is not None else "n/d",
            f(c["media"]),
            f"{abs(s_mano - c['media']):.1e}" if s_mano is not None else "n/d",
            f(c["sd1"]), f(c["term_r"]), f(c["margen"]),
        ])
    tabla(L, ["linea", "metrica", "los 8 u crudos", "u medio", "sd_u", "b", "r",
              "span", "s a mano desde u medio", "s medio medido",
              "\\|dif\\|", "sd_s", "termino de r", "margen"], rows)
    L += [
        f"Son **{len(rows)}** filas. La columna `|dif|` compara el `s` que sale de rehacer",
        "`(u_medio - b)/(r - b)` contra la media de los 8 `s` que escribio `score`. Donde",
        "no es cero del todo es porque la media de una funcion lineal es la funcion de la",
        "media solo si no hay clamp ni gates por el camino -- y donde el clamp actua",
        "(`mse`, que no esta en esta seleccion) la identidad se rompe del todo.",
        "",
    ]

    # ------------------------------------------------ 8. cierre
    L += [
        "## Lo que no se hizo",
        "",
        "- **No se construyo ningun candidato.** La emision es `predict_floor` sin",
        "  modificar, la misma desechable de `cable_suelo.md`.",
        "- **No se tocaron los bundles ni ningun arnes.** `build_fold_adatas`,",
        "  `predict_floor` y `contar_de_por_perturbacion` se usaron tal cual.",
        "- **No se promediaron las 8 corridas en una submission.** Se puntuaron por",
        "  separado contra el mismo bundle; lo que se reporta es su dispersion.",
        "- **No se propone ningun margen ni umbral.** La seccion 2 da la formula y la",
        "  seccion 3 la aplica celda por celda; la regla de decision ya esta escrita en",
        "  §3.16 y no se toca aqui.",
        "- No se leyo nada de `~/vcc/benchmark/` ni de",
        "  `/datos/benchmark_VEDADO_techo_metricas/`.",
        "- Sin envios, sin `vcc sample`, sin commits.",
        "",
    ]

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    sys.stderr.write(f"escrito {OUT_MD} ({len(L)} lineas)\n")


if __name__ == "__main__":
    main()
