#!/usr/bin/env python
"""Escribe modelo/out/huecos.md: los cuatro huecos de §3.14.

Lee, no mide. Las mediciones estan en:
  * descomp_rgen.json    -- las 5 particiones del anchor recomputadas, 6 lineas
  * huecos_cohorte.json  -- lado real por perturbacion de los dos objetos
  * huecos_pds.json      -- distribucion de pds por perturbacion
  * anclas_nativas_comp.json -- `b` y `r` oficiales por metrica y linea
"""
import json
import os
import sys

import numpy as np
import pandas as pd

OUT_MD = "/home/ubuntu/vcc/modelo/out/huecos.md"
RGEN = "/home/ubuntu/vcc/modelo/out/descomp_rgen.json"
COH = "/home/ubuntu/vcc/modelo/out/huecos_cohorte.json"
PDSJ = "/home/ubuntu/vcc/modelo/out/huecos_pds.json"
NATIVAS = "/home/ubuntu/vcc/modelo/out/anclas_nativas_comp.json"
BUNDLES = "/datos/modelo/out/nativas_comp"

LINEAS = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
N_PANEL = 33
JAC = "de_wilcoxon_sig_jaccard"
REACH = "de_wilcoxon_direction_reach_raw"
FID = "de_wilcoxon_direction_fidelity_yield_raw"
LFC = "de_wilcoxon_lfc_nmae"
PDS = "pds_cosine"
C_FIT = 0.96          # direction.py:967 -- "the fitted ``c`` is ~0.96 for a coin-flip predictor"
MIN_GATE = 10         # metrics/de.py:601 -- min_gate_size


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


# ---------------------------------------------------------------- cabecera
def cabecera(L):
    L += [
        "# Los cuatro huecos de §3.14",
        "",
        "Cuatro preguntas independientes, medidas, no inferidas. Ningun criterio de",
        "adopcion se redacta aqui: no es de este encargo.",
        "",
        "**Fuentes de medicion** (todas de esta sesion, todas sobre la MISMA poblacion de",
        "siempre -- `build_fold_adatas` sin modificar, `seed_real = 20260001 + fold_idx`,",
        "panel comun de 33 de `suelo_panel.json`, etiqueta de control `non-targeting`):",
        "",
        "| fuente | que contiene | como se obtuvo |",
        "|---|---|---|",
        "| `descomp_rgen.json` | las 5 particiones del anchor recomputadas, **las 6 lineas** | `cell_eval2.ceiling._disjoint_halves` sin modificar + las semillas que estampo cada bundle |",
        "| `huecos_cohorte.json` | lado real por perturbacion (`n_DE_real`, `n_conf`, `q`) y valores por perturbacion de los 6 metricas, para nuestro suelo y para el `b` oficial | `cell-eval2 run --write-degenes` + `cell_eval2.metrics.direction._reference_stats` sin reimplementar |",
        "| `huecos_pds.json` | distribucion de `pds_cosine` por perturbacion | `results.csv` / `baseline_results.csv` ya en disco |",
        "| `anclas_nativas_comp.json` | `b` y `r` oficiales de las 6 metricas por linea | los 6 bundles de competencia de `anclas_nativas_v2.md` |",
        "",
        "**Criterios de seleccion, declarados antes de medir.** Las SEIS lineas. Los 33 del",
        "panel comun. Ninguna celda excluida por su valor ni por su signo. Donde una tabla",
        "muestra menos de 33 filas es porque la metrica no emitio la fila (gate del lado",
        "real), y el conteo de cuantas faltan ES el resultado, no un recorte.",
        "",
        "**Una desviacion instrumental, declarada.** Para obtener las tablas DE del `b`",
        "oficial hubo que volver a correr su emision por `cell-eval2 run`, y `run` la rechaza",
        "con `declared input_type='counts' but values are fractional` (`norm.py:482`): el",
        "perfil del baseline es una MEDIA, fraccionaria en cualquier espacio. Se paso",
        "`--set allow_fractional_counts=true`, que es **lo que `cell-eval2 baseline` se pone",
        "a si mismo** (`baseline.py:1320-1321`) y que esta en `DIGEST_EXEMPT_FIELDS`",
        "(`baseline.py:692-703`: *\"a validation ALLOWANCE, not a scoring semantic\"*). Aplica",
        "al lado PRED solamente (`run.py:397-401`), el lado real sigue validado estricto. No",
        "cambia ninguna semantica de puntuacion; la verificacion de abajo lo comprueba",
        "numero a numero.",
        "",
    ]


def resumen(L):
    L += [
        "## Los cuatro veredictos, en una linea cada uno",
        "",
        "| hueco | la pregunta | el veredicto |",
        "|---|---|---|",
        "| **1** | ¿la convencion resta del span, o eso pasa solo donde `E` es grande? | **No se generaliza.** Resta en 2 de 6 (A549, HAP1) y suma en 4. El signo lo decide cuanta senal DE tiene el panel, no el tamano de `E`. |",
        "| **2** | ¿la cohorte del lado real es la misma para los tres objetos? | **Si entre submissions** (medido celda por celda, 6 lineas), **no entre poblaciones**: las particiones del anchor corren sobre media poblacion y su cohorte de `fid` y `reach` se mueve. Y la cohorte que la componente genuina de `jac` necesita -- union no vacia -- NO es la del lado real y si la mueve la submission. |",
        "| **3** | ¿nuestro `b` de `reach` cae sobre nuestro `c/N_conf`? | **No, y por razones distintas en cada punta.** 4 lineas por encima (hasta 9.9x), K562 encima, A549 un 6% por debajo. El `b` esta hecho en un 53-88% de 4 a 7 perturbaciones con `N_conf <= 10` que sacan 1.0 exacto. |",
        "| **4** | ¿explica el empate en 0.5 el `b < u` de HAP1 y BXPC3? | **No. Hipotesis refutada.** Ningun delta tiene norma cero, el `b` tiene mas 0.5 exactos que el suelo, y 0.5 esta en la rejilla asi que ese conteo no mide empates. Lo que explica las seis diferencias es el `SE` propio de una media de 33 rangos casi uniformes, 0.0518: HAP1 queda en 0.97 sigma. |",
        "",
    ]


def verif_b(L, coh):
    """Que `cell-eval2 run` sobre la emision del baseline reproduzca
    `baseline_results.csv`. Si no lo hiciera, los numeros por perturbacion del
    `b` de este reporte no serian los del `b` oficial."""
    L += [
        "**Verificacion de que el `b` por perturbacion de este reporte ES el `b` oficial.**",
        "Dos comprobaciones encadenadas, las dos celda por celda:",
        "",
        "1. la emision del baseline regenerada con `--seed 0` reproduce el",
        "   `baseline_results.csv` del `BASE_<linea>` que construyo los bundles;",
        "2. `cell-eval2 run` sobre esa emision reproduce los mismos valores por",
        "   perturbacion que `cell-eval2 baseline`.",
        "",
    ]
    rows = []
    celdas = []
    for ln in LINEAS:
        if ln not in coh:
            continue
        v = coh[ln]["verif_baseline_regenerado"]
        for x in v["difs"]:
            celdas.append([ln, x["pert"], x["metric"], repr(x["oficial"]),
                           repr(x["regen"]),
                           f"{x['regen'] - x['oficial']:+.3e}"
                           if (x["oficial"] is not None and x["regen"] is not None)
                           else "n/d"])
        of = pd.read_csv(os.path.join(BUNDLES, f"BASE_{ln}", "baseline_results.csv"))
        of = {(str(a), str(b)): float(c)
              for a, b, c in zip(of["perturbation"], of["metric"], of["value"])}
        pp = coh[ln]["objetos"]["b_oficial"]["por_pert"]
        mapa = {"jac": JAC, "reach": REACH, "fid": FID, "lfc_nmae": LFC, "pds": PDS}
        n_cmp = n_dif = 0
        peor = 0.0
        for p, d in pp.items():
            for k, m in mapa.items():
                a, b = d[k], of.get((p, m))
                if a is None and (b is None or np.isnan(b)):
                    continue
                if a is None or b is None or np.isnan(b):
                    n_dif += 1
                    n_cmp += 1
                    continue
                n_cmp += 1
                if a != b:
                    n_dif += 1
                    peor = max(peor, abs(a - b))
        peor1 = max((abs(x["regen"] - x["oficial"]) for x in v["difs"]
                     if x["oficial"] is not None and x["regen"] is not None),
                    default=0.0)
        rows.append([
            ln, v["n_celdas_comparadas"], v["n_difieren"],
            f"{peor1:.2e}" if v["n_difieren"] else "0",
            n_cmp, n_dif, f"{peor:.2e}" if n_dif else "0",
        ])
    tabla(L, ["linea", "celdas (1) baseline regenerado vs oficial", "difieren",
              "max \\|dif\\|", "celdas (2) run vs baseline", "difieren",
              "max \\|dif\\|"], rows)
    peores = [float(r[3]) for r in rows] + [float(r[6]) for r in rows]
    tope = max(peores)
    L += [
        f"Ninguna diferencia pasa de **{tope:.2e}**. No son diferencias: son los 2-3 ULP",
        "de coma flotante que `anclas_nativas_v2.md` ya habia aislado en el PASO 0, y que",
        "vuelven a salir **en las mismas celdas**:",
        "",
    ]
    if celdas:
        tabla(L, ["linea", "perturbacion", "metrica", "oficial", "regenerado",
                  "dif"], celdas)
        L += [
            "Todas en `de_wilcoxon_lfc_nmae`. Que reaparezcan las mismas perturbaciones en",
            "una corrida nueva confirma lo que el PASO 0 habia concluido: la cadena no es",
            "bit-reproducible corrida a corrida, y el efecto esta localizado. **No depende",
            "de `allow_fractional_counts`**: la comprobacion (1) no lo usa en absoluto y",
            "produce la misma desviacion que la (2).",
            "",
        ]
    else:
        L += ["No hubo ninguna celda distinta: identidad bit a bit en las dos "
              "comprobaciones.", ""]
    L += [
        "El `b` por perturbacion de este reporte es, a 2-3 ULP, el `b` oficial, y",
        "`allow_fractional_counts` no movio nada.",
        "",
    ]


# ------------------------------------------------------- hueco 1
def hueco1(L, rgen, nat):
    L += [
        "---",
        "",
        "## Hueco 1 — la descomposicion en las cuatro lineas que faltaban",
        "",
        "`descomposicion.md` recomputo las 5 particiones del anchor solo en A549 y HAP1.",
        "Aqui estan BXPC3, HT29, K562 y MCF7 con el mismo metodo, y las seis juntas.",
        "",
        "### 1.1 Lo primero: el valor recomputado reproduce el del bundle",
        "",
        "Antes de descomponer nada. Si la reproduccion fallara, la descomposicion no",
        "describiria el `r` que puntua.",
        "",
    ]
    rows = []
    for ln in LINEAS:
        if ln not in rgen:
            rows.append([ln, "**falta**", "-", "-", "-"])
            continue
        d = rgen[ln]
        difs = [s["dif_vs_bundle"] for s in d["splits"]]
        rows.append([
            ln, len(d["splits"]), f(d["r_crudo_media"]),
            f"{max(difs):.2e}", "si" if max(difs) < 1e-12 else "**NO**",
        ])
    tabla(L, ["linea", "particiones", "r crudo (media)",
              "max \\|recomputado - bundle\\|", "reproduce"], rows)
    peor = max(max(s["dif_vs_bundle"] for s in rgen[ln]["splits"])
               for ln in LINEAS if ln in rgen)
    L += [
        f"La desviacion maxima en las 30 particiones es **{peor:.2e}** -- ruido de ultimo bit",
        "de coma flotante, no una diferencia. La descomposicion describe el `r` que puntua.",
        "",
        "### 1.2 Las 20 particiones nuevas, fila por fila",
        "",
        "Los intermedios crudos de los que sale todo lo demas. `E` = perturbaciones con",
        "`|union| = 0`, que cobran 1.0 cada una por la convencion D-023.",
        "",
    ]
    rows = []
    for ln in ["BXPC3", "HT29", "K562", "MCF7"]:
        if ln not in rgen:
            continue
        for s in rgen[ln]["splits"]:
            rows.append([
                ln, s["split_index"], s["seed"], s["n_half_a"], s["n_half_b"],
                s["E_union_vacia"], s["n_genuinos"], f(s["suma_genuinos"], 4),
                f(s["media_genuinos"]), f(s["componente_convencion"]),
                f(s["componente_genuina"]), f(s["valor_recomputado"]),
                f(s["valor_del_bundle"]), f"{s['dif_vs_bundle']:.1e}",
            ])
    tabla(L, ["linea", "split", "semilla", "n half_a", "n half_b", "E", "n genuinos",
              "suma genuinos", "media genuinos", "E/33", "suma/33", "jac recomputado",
              "jac del bundle", "dif"], rows)
    L += [
        f"Son **{len(rows)}** filas: 4 lineas x 5 particiones. La identidad",
        "`jac = (E + suma_genuinos)/33` se lee directamente en cada una: la columna `E/33`",
        "mas la columna `suma/33` da la columna `jac recomputado`.",
        "",
        "### 1.3 E y la media de los genuinos, las seis lineas",
        "",
    ]
    rows = []
    for ln in LINEAS:
        if ln not in rgen:
            continue
        d = rgen[ln]
        Es = [s["E_union_vacia"] for s in d["splits"]]
        gm = [s["media_genuinos"] for s in d["splits"]]
        rows.append([
            ln, "/".join(str(e) for e in Es), f"{np.mean(Es):.1f}",
            f"{min(Es)}-{max(Es)}",
            " / ".join(f(g, 4) for g in gm), f(float(np.mean(gm)), 6),
            f(d["r_crudo_media"]), f(d["r_crudo_sd"]),
            f(d["r_convencion_media"]), f(d["r_genuina_media"]),
        ])
    tabla(L, ["linea", "E por particion", "E medio", "rango E",
              "media genuinos por particion", "media genuinos",
              "r crudo", "sd r crudo", "componente convencion de r",
              "componente genuina de r"], rows)

    L += [
        "### 1.4 El aporte de la convencion al span, contra el aporte de lo genuino",
        "",
        "El span crudo es `r - b`. Las dos componentes son aditivas porque la convencion",
        "entra como `E/33` y lo genuino como `suma_genuinos/33`, y las dos suman el crudo.",
        "Asi que:",
        "",
        "```",
        "span_crudo      = r_crudo      - b_crudo",
        "aporte_conv     = E_r/33       - E_b/33",
        "aporte_genuino  = suma_gen_r/33 - suma_gen_b/33",
        "span_crudo      = aporte_conv  + aporte_genuino      (identidad, se verifica abajo)",
        "```",
        "",
        "`E_b` y `suma_gen_b` del `b` oficial salen de `huecos_cohorte.json`, que trae",
        "`|union|` por perturbacion MEDIDO (no inferido del `jac == 1.0`).",
        "",
    ]
    return rows


def hueco1_span(L, rgen, coh, nat):
    porb = {r["linea"]: r for r in nat}
    rows = []
    veredicto = []
    for ln in LINEAS:
        if ln not in rgen or ln not in coh:
            continue
        pp = coh[ln]["objetos"]["b_oficial"]["por_pert"]
        E_b = sum(1 for v in pp.values() if v["n_DE_union"] == 0)
        gen_b = sum(v["jac"] for v in pp.values()
                    if v["n_DE_union"] != 0 and v["jac"] is not None)
        b_crudo = porb[ln]["b"][JAC]
        conv_b, gb = E_b / N_PANEL, gen_b / N_PANEL
        d = rgen[ln]
        E_r = float(np.mean([s["E_union_vacia"] for s in d["splits"]]))
        gen_r = float(np.mean([s["suma_genuinos"] for s in d["splits"]]))
        conv_r, gr = E_r / N_PANEL, gen_r / N_PANEL
        span = d["r_crudo_media"] - b_crudo
        ap_c, ap_g = conv_r - conv_b, gr - gb
        rows.append([
            ln, E_b, f"{E_r:.1f}", f(conv_b), f(conv_r), f(gb), f(gr),
            f(b_crudo), f"{abs(b_crudo - (conv_b + gb)):.1e}",
            f(d["r_crudo_media"]), f(span),
            f"**{ap_c:+.6f}**", f"**{ap_g:+.6f}**",
            f"{ap_c + ap_g:.6f}", f"{abs(span - (ap_c + ap_g)):.1e}",
        ])
        veredicto.append((ln, E_b, E_r, ap_c, ap_g, span))
    tabla(L, ["linea", "E del b", "E del r (medio)", "conv b = E_b/33",
              "conv r = E_r/33", "genuina b", "genuina r", "b crudo",
              "\\|b - (conv+gen)\\|", "r crudo",
              "span crudo", "aporte CONVENCION", "aporte GENUINO", "suma de aportes",
              "\\|span - suma\\|"], rows)
    L += [
        "Dos verificaciones en la misma tabla. `|b - (conv+gen)|` comprueba que el `b`",
        "oficial del bundle se reconstruye desde los conteos por perturbacion que medimos",
        "aqui -- es la identidad `jac = (E + suma_genuinos)/33` aplicada al `b`. Y",
        "`|span - suma|` confirma que los dos aportes suman el span crudo en las seis",
        "lineas.",
        "",
        "#### Veredicto de §1.4: la afirmacion con las seis a la vista",
        "",
    ]
    restan = [v for v in veredicto if v[3] < 0]
    suman = [v for v in veredicto if v[3] > 0]
    L += [
        f"**La convencion resta del span en {len(restan)} de 6 lineas"
        f"{' y suma en ' + str(len(suman)) if suman else ''}.**",
        "",
    ]
    for ln, E_b, E_r, ap_c, ap_g, span in veredicto:
        signo = "resta" if ap_c < 0 else ("suma" if ap_c > 0 else "no mueve")
        L.append(
            f"- **{ln}**: E_b = {E_b}, E_r medio = {E_r:.1f} -> la convencion {signo} "
            f"{abs(ap_c):.6f} del span de {span:.6f} "
            f"({100*abs(ap_c)/abs(span) if span else 0:.1f}% de su magnitud); "
            f"lo genuino aporta {ap_g:+.6f}."
        )
    L.append("")
    # el mecanismo contra el tamano de E
    Ebs = np.array([v[1] for v in veredicto], dtype=float)
    apcs = np.array([v[3] for v in veredicto])
    if len(Ebs) > 2:
        rr = float(np.corrcoef(Ebs, apcs)[0, 1])
    else:
        rr = float("nan")
    L += [
        "**La afirmacion de `descomposicion.md` -- \"la convencion resta del span\" -- no se",
        "generaliza.** Estaba medida en A549 y HAP1, y en las otras cuatro lineas la",
        "convencion **suma**. Lo que hay no es un mecanismo con un signo fijo sino una",
        "resta de dos cantidades con origenes distintos, y con las seis a la vista se ve",
        "cual es cual:",
        "",
        "- **`E_b` es, salvo cuatro celdas, el complemento de la cohorte del lado real**:",
        "  `33 - |n_conf > 0|`. Son 15, 0, 14, 6, 8 y 0 (§2.5 lo mide y nombra las cuatro",
        "  excepciones). Mide **cuantas perturbaciones del panel no tienen senal DE en la",
        "  referencia**. No tiene nada que ver con el baseline: nuestro suelo da el mismo",
        "  numero.",
        "- **`E_r` lo decide la perdida de potencia de media poblacion**: 13.8, 2.0, 11.2,",
        "  6.2, 12.0 y 0.6. Es cuantas quedan sin senal **a los dos lados de la particion**.",
        "",
        "Y por eso el signo cambia: donde el panel tiene poca senal (A549 con 18 de 33 en",
        "cohorte, HAP1 con 19) `E_b` ya es alto y la particion no puede subirlo mucho mas,",
        "asi que `E_r < E_b` y la convencion resta. Donde el panel tiene mucha senal (BXPC3",
        "y MCF7, 33 de 33) `E_b = 0` por construccion y cualquier perdida de potencia lo",
        "sube, asi que la convencion suma. La correlacion entre `E_b` y el aporte de la",
        f"convencion sobre las seis lineas es **{rr:+.3f}**.",
        "",
        "**No hace falta que `E` sea grande para que reste; hace falta que el panel tenga",
        "poca senal.** Son la misma cosa en estas seis lineas, y por eso los dos casos donde",
        "resta son exactamente los dos de cohorte mas chica.",
        "",
        "Y la magnitud esta cuantizada: cada unidad de E mueve el crudo exactamente",
        f"`1/33 = {1/33:.6f}`. Por eso el `sd` entre particiones del crudo es grande (E",
        "cambia de particion a particion) mientras el de la componente genuina es chico.",
        "",
    ]
    return veredicto


# ------------------------------------------------------- hueco 2
def hueco2(L, coh, rgen):
    L += [
        "---",
        "",
        "## Hueco 2 — la cohorte comun de la componente genuina",
        "",
        "### 2.1 La hipotesis, y por que esa y no otra",
        "",
        "> La cohorte debe definirse por el LADO REAL unicamente: las perturbaciones donde",
        "> la referencia tiene senal DE. Razon: eso la vuelve invariante a la submission.",
        "",
        "El codigo dice por que esa es la propiedad que importa, y lo dice dos veces:",
        "",
        "- `metrics/de.py:608-612` (gate de `lfc_nmae`): *\"the gate rules below drop a",
        "  perturbation whose gate is empty, smaller than `min_gate_size`, or whose",
        "  denominator is zero or non-finite -- those score nothing rather than 1.0. **All",
        "  three are real-side decisions, so the surviving set is the same for every",
        "  submission.**\"*",
        "- `metrics/direction.py:296-299` (`_reference_stats`, que produce el `n_conf` de",
        "  `reach`): *\"Never touches the join: `_direction_frame` is an INNER join, so",
        "  anything computed from it is prediction-dependent, and `q`/`N_conf` must not be",
        "  -- **that invariance is what stops a model shrinking its own BUDGET by omitting",
        "  genes**.\"*",
        "",
        "Asi que la hipotesis se prueba en dos mitades, y las dos se miden:",
        "",
        "1. ¿es la misma entre SUBMISSIONS sobre la misma poblacion? (nuestro suelo contra",
        "   el `b` oficial)",
        "2. ¿es la misma entre OBJETOS, incluidas las particiones del anchor, que corren",
        "   sobre media poblacion?",
        "",
        "### 2.2 Mitad 1: el lado real es invariante a la submission. Medido.",
        "",
        "Dos emisiones completamente distintas sobre el MISMO `real.h5ad` de cada linea:",
        "nuestro suelo (control remuestreado, sin escalar) y el `b` oficial (control",
        "escalado al perfil de respuesta promedio). Se comparan los conteos del lado real",
        "perturbacion por perturbacion.",
        "",
    ]
    rows = []
    todo_igual = True
    for ln in LINEAS:
        if ln not in coh:
            continue
        u = coh[ln]["objetos"]["suelo"]["por_pert"]
        b = coh[ln]["objetos"]["b_oficial"]["por_pert"]
        perts = sorted(set(u) | set(b))
        dif_real = [p for p in perts
                    if u.get(p, {}).get("n_DE_real") != b.get(p, {}).get("n_DE_real")]
        dif_conf = [p for p in perts
                    if u.get(p, {}).get("n_conf") != b.get(p, {}).get("n_conf")]
        dif_pred = [p for p in perts
                    if u.get(p, {}).get("n_DE_pred") != b.get(p, {}).get("n_DE_pred")]
        if dif_real or dif_conf:
            todo_igual = False
        rows.append([
            ln, len(perts), len(dif_real), len(dif_conf), len(dif_pred),
            sum(v["n_DE_pred"] for v in u.values()),
            sum(v["n_DE_pred"] for v in b.values()),
            "**si**" if not (dif_real or dif_conf) else "**NO**",
        ])
    tabla(L, ["linea", "perturbaciones", "difieren en n_DE_real",
              "difieren en n_conf", "difieren en n_DE_pred (lado predicho)",
              "suma n_DE_pred suelo", "suma n_DE_pred b",
              "lado real identico"], rows)
    L += [
        ("**El lado real es identico celda por celda en las seis lineas.**"
         if todo_igual else
         "**ATENCION: el lado real NO salio identico. Ver la tabla.**"),
        "",
        "Las tres ultimas columnas estan ahi como control, y hay que leerlas con cuidado:",
        "el lado predicho difiere en 0, 3, 0, 12, 2 y 6 perturbaciones de 33 -- **no en",
        "casi todas**. La razon esta en las dos columnas de suma: las dos emisiones",
        "predicen **casi ningun gen significativo en ningun sitio**. Nuestro suelo suma 0,",
        "1, 0, 3, 0 y 0 genes significativos sobre las 33 perturbaciones de cada linea; el",
        "`b` oficial suma 0, 2, 0, 14, 2 y 7. En A549 y HAP1 las dos suman cero, asi que",
        "ahi el control positivo **no discrimina**: las dos emisiones coinciden en el lado",
        "predicho porque las dos estan vacias, no porque sean la misma. Donde si",
        "discrimina -- HT29, con 12 perturbaciones distintas -- el lado real sigue saliendo",
        "identico celda por celda. Esa es la evidencia que sostiene la invariancia; las",
        "otras cinco lineas la acompanan sin probarla por si solas.",
        "",
        "### 2.3 Mitad 2: pero NO es invariante a la POBLACION",
        "",
        "El `r` del anchor no se calcula sobre la poblacion completa: cada particion puntua",
        "`half_b` contra `half_a` (`anchor.py:120`), media poblacion a cada lado. Y el gate",
        "es del lado real, luego depende de cuantas celulas tiene ese lado.",
        "",
        "`anchor_splits.parquet` estampa `n_perturbations` por metrica y particion. Eso es",
        "el tamano de la cohorte, directamente del artefacto oficial:",
        "",
    ]
    import polars as pl
    rows = []
    for ln in LINEAS:
        d = pl.read_parquet(os.path.join(BUNDLES, f"BUNDLE_{ln}", "anchor_splits.parquet"))
        fila = [ln]
        for m in (PDS, JAC, LFC, FID, REACH):
            v = d.filter(pl.col("metric") == m).sort("split_index")["n_perturbations"].to_list()
            fila.append("/".join(str(x) for x in v))
        # cohorte del run completo, del mismo objeto que puntua: el b oficial
        if ln in coh:
            b = coh[ln]["objetos"]["b_oficial"]["por_pert"]
            n_lfc = sum(1 for v in b.values() if v["lfc_nmae"] is not None)
            n_fid = sum(1 for v in b.values() if v["fid"] is not None)
            n_rch = sum(1 for v in b.values() if v["reach"] is not None)
            fila.append(f"{n_lfc} / {n_fid} / {n_rch}")
        else:
            fila.append("n/d")
        rows.append(fila)
    tabla(L, ["linea", "pds (5 particiones)", "jac", "lfc_nmae", "fid", "reach",
              "run completo: lfc/fid/reach"], rows)
    L += [
        "Dos cosas, las dos medidas:",
        "",
        "- **`lfc_nmae` no se mueve.** Su cohorte es la misma en las 5 particiones Y en el",
        "  run completo, en las seis lineas. Su gate es `>= min_gate_size = 10` genes",
        "  significativos del lado real (`metrics/de.py:601`), y perder la mitad de las",
        "  celulas no cruza ese umbral en ninguna perturbacion de este panel.",
        "- **`fid` y `reach` si se mueven**, y entre particiones de la misma linea. Su",
        "  cohorte es `n_conf > 0`: basta un gen que pierda la significancia al quedarse",
        "  con media poblacion para que la perturbacion salga.",
        "",
        "### 2.4 Las cohortes, por linea",
        "",
        "La cohorte del lado real que la hipotesis propone -- `n_conf > 0`, las",
        "perturbaciones donde la referencia tiene senal DE -- medida sobre la poblacion",
        "completa, con `_reference_stats` del marcador:",
        "",
    ]
    rows = []
    for ln in LINEAS:
        if ln not in coh:
            continue
        u = coh[ln]["objetos"]["suelo"]["por_pert"]
        n_sig = sum(1 for v in u.values() if v["n_DE_real"] > 0)
        n_conf = sum(1 for v in u.values() if (v["n_conf"] or 0) > 0)
        n_gate = sum(1 for v in u.values() if (v["n_conf"] or 0) >= MIN_GATE)
        n_q = sum(1 for v in u.values() if v["q"] is not None)
        nc = [v["n_conf"] or 0 for v in u.values()]
        rows.append([
            ln, N_PANEL, n_sig, n_conf, n_gate, n_q,
            f"{min(nc)} / {int(np.median(nc))} / {max(nc)}",
        ])
    tabla(L, ["linea", "panel", "n_DE_real > 0", "n_conf > 0 (**la cohorte**)",
              "n_conf >= 10 (gate de lfc)", "q definido",
              "n_conf min/mediana/max"], rows)
    L += [
        "`n_DE_real > 0` y `n_conf > 0` dan el mismo numero: los dos cuentan genes",
        "significativos del lado real con el gen propio excluido, y son el mismo conteo",
        "escrito en dos sitios del paquete (`contar_de_por_perturbacion`, que repite las dos",
        "lineas de filtrado de `de_sig_jaccard`, y `_reference_stats`, direction.py:343-350).",
        "Que coincidan no era obligatorio y se comprueba aqui.",
        "",
        "### 2.5 Veredicto de Hueco 2",
        "",
        "**La cohorte del lado real SI es invariante a la submission -- medido, celda por",
        "celda, en las seis lineas -- y NO es invariante a la poblacion.** Las dos mitades",
        "de la hipotesis no caen del mismo lado, y la consecuencia es concreta:",
        "",
        "- comparar `u` contra `b` sobre la cohorte del lado real es legitimo: es la misma",
        "  cohorte, por construccion y por medicion;",
        "- comparar cualquiera de los dos contra `r` sobre \"la\" cohorte no lo es sin mas,",
        "  porque `r` se calcula sobre media poblacion y su cohorte de `fid` y `reach` no es",
        "  la del run completo. En `jac` y `pds` no estorba (33/33 siempre); en `lfc_nmae`",
        "  tampoco (no se movio en ninguna linea ni particion); en `fid` y `reach` si.",
        "",
        "Lo que la cohorte del lado real NO resuelve: la cohorte de `|union| != 0` que la",
        "componente genuina de `jac` necesita. `|union| = |R u P|` incluye el lado predicho,",
        "asi que **esa** cohorte no es invariante a la submission. Medido:",
        "",
    ]
    rows = []
    for ln in LINEAS:
        if ln not in coh:
            continue
        u = coh[ln]["objetos"]["suelo"]["por_pert"]
        b = coh[ln]["objetos"]["b_oficial"]["por_pert"]
        Vu = {p for p, v in u.items() if v["n_DE_union"] == 0}
        Vb = {p for p, v in b.items() if v["n_DE_union"] == 0}
        Er = [s["E_union_vacia"] for s in rgen.get(ln, {}).get("splits", [])]
        cohorte = sum(1 for v in u.values() if (v["n_conf"] or 0) > 0)
        mueven = sorted(p for p, v in b.items()
                        if (v["n_conf"] or 0) == 0 and v["n_DE_union"] > 0)
        rows.append([
            ln, len(Vu), len(Vb), len(Vu & Vb), len(Vu ^ Vb),
            "/".join(str(e) for e in Er) if Er else "n/d",
            "si" if Vu == Vb else "**no**",
            N_PANEL - cohorte,
            ", ".join(mueven) if mueven else "--",
        ])
    tabla(L, ["linea", "E de nuestro suelo", "E del b oficial", "en los dos",
              "en uno solo", "E de las 5 particiones de r", "misma cohorte u/b",
              "33 - cohorte", "celdas que la mueven"], rows)
    L += [
        "**Las dos cohortes estan ligadas por una identidad, y la medicion la fija.** La",
        "columna `E de nuestro suelo` coincide exactamente con `33 - cohorte` en las seis",
        "lineas. No es casualidad: `|union| = |R u P|`, nuestro suelo no predice **ningun**",
        "gen significativo en casi ningun sitio (§2.2: suma 0 en cuatro lineas de seis),",
        "asi que su union es vacia exactamente donde el lado real esta vacio. **La `E` de",
        "nuestro suelo ES el complemento de la cohorte del lado real.**",
        "",
        "Y la columna `celdas que la mueven` nombra lo que rompe esa identidad para el `b`:",
        "las perturbaciones con `n_conf = 0` donde el `b` SI predijo algo significativo.",
        "Son cuatro en total en las seis lineas: **BRD9, IRF7 y RAP1GAP en HT29, y ZBP1 en",
        "K562**. En esas cuatro el `b` pierde el 1.0 gratis de la convencion y se queda con",
        "un jaccard genuino de 0 (interseccion 0 sobre union 2), que es por lo que su `E` es",
        "3 y 7 donde la de nuestro suelo es 6 y 8. Cada una cuesta `1/33 = 0.030303`.",
        "",
        "Eso es la demostracion concreta, con nombre y apellido, de que la cohorte de union",
        "no vacia **la mueve la submission**: cuatro perturbaciones de 198, y mueven el `b`",
        "de dos lineas.",
        "",
        "La cohorte de union no vacia **cambia con la submission** -- es la columna `en uno",
        "solo` -- y cambia tambien de particion a particion del anchor. Con `u` y `b`",
        "genuinos en 0 no estorbaba porque el promedio de un conjunto vacio de jaccards",
        "genuinos es 0 en los dos casos. Con un candidato que SI tenga senal, la media",
        "genuina de `u` y la de `b` se calcularian sobre conjuntos distintos.",
        "",
        "**Eso no lo arregla la cohorte del lado real**: una perturbacion con",
        "`n_conf = 0` puede tener union no vacia (si la submission predice algo",
        "significativo) y una con `n_conf > 0` siempre la tiene. La cohorte del lado real",
        "es un SUBCONJUNTO de la de union no vacia, y la diferencia es exactamente lo que la",
        "submission controla. Lo que la medicion dice es cual de las dos elegir si se quiere",
        "una cohorte no manipulable; no dice que la componente genuina de `jac` quede bien",
        "definida sobre ella. Para eso haria falta decidir que recibe una perturbacion con",
        "`n_conf > 0` cuyo `|union|` sea no vacio solo por el lado predicho, y esa decision",
        "no es de este encargo.",
        "",
    ]


# ------------------------------------------------------- hueco 3
def hueco3(L, coh, nat):
    porb = {r["linea"]: r for r in nat}
    L += [
        "---",
        "",
        "## Hueco 3 — ¿nuestro `b` de `reach` cae sobre NUESTRO `c/N_conf`?",
        "",
        "### 3.1 De donde sale cada pieza",
        "",
        "- **`reach = k* / N_conf`**, `metrics/direction.py:1054`.",
        "- **`N_conf`** lo construye `_reference_stats` desde `prepared.real_df` SOLO",
        "  (`direction.py:293`, `343-350`): genes del lado real con `p_adj < alpha`,",
        "  excluido el gen propio del target. No se persiste en ningun artefacto; aqui se",
        "  obtiene llamando a esa funcion, sin reimplementarla.",
        "- **`c`**: se usa **`c = 0.96`**, de `direction.py:967`: *\"the fitted `c` is ~0.96",
        "  for a coin-flip predictor against ~1.90 for one that ranks a single",
        "  HABITUAL-DIRECTION gene first\"*. Es un ajuste medido por el propio paquete en",
        "  0.13.0 sobre 200 replicas por celda (`direction.py:965-966`), no una constante",
        "  del codigo ejecutable: no hay ningun `c` en el modulo. El punto sin habilidad es",
        f"  `~c/N_conf` (`direction.py:946`), **no** una constante.",
        "",
        "⚠️ El propio codigo acota donde la ley falla, y hay que decirlo antes de usarla:",
        "`E[reach_raw]` medido *\"runs 0.5100 at `N_conf = 1` to 0.0019 at 500\"*",
        "(`direction.py:965-966`). En `N_conf = 1` la ley daria 0.96 y lo medido fue 0.51:",
        "la ley sobreestima donde `N_conf` es diminuto. La tabla trae la distribucion de",
        "`N_conf` para que se vea cuanto peso tiene esa zona en nuestras lineas.",
        "",
        "### 3.2 El azar esperado sobre NUESTRAS perturbaciones",
        "",
        "Para cada linea: la media de `c/N_conf` sobre las perturbaciones que **entran** a",
        "`reach` (las que tienen `reach` no-NaN en el run, que son exactamente las de",
        "`N_conf > 0`, `direction.py:1048-1053`), contra nuestro `b` oficial de `reach`.",
        "",
    ]
    rows = []
    det = []
    c_impls = []
    for ln in LINEAS:
        if ln not in coh:
            continue
        b = coh[ln]["objetos"]["b_oficial"]["por_pert"]
        u = coh[ln]["objetos"]["suelo"]["por_pert"]
        dentro = {p: v for p, v in b.items() if v["reach"] is not None}
        nc = np.array([v["n_conf"] for v in dentro.values()], dtype=float)
        azar = C_FIT / nc
        b_of = porb[ln]["b"][REACH]
        b_run = float(np.mean([v["reach"] for v in dentro.values()]))
        u_run = float(np.mean([u[p]["reach"] for p in dentro
                               if u[p]["reach"] is not None]))
        c_impl = b_of / float(np.mean(1.0 / nc))
        rows.append([
            ln, len(dentro), f"{int(nc.min())}", f"{int(np.median(nc))}",
            f"{int(nc.max())}", int((nc <= 3).sum()),
            f(float(azar.mean())), f(b_of), f(b_run),
            f"{b_of - float(azar.mean()):+.6f}",
            f"{b_of/float(azar.mean()):.3f}x", f"{c_impl:.4f}", f(u_run),
        ])
        det.append((ln, dentro, nc, azar, b_of))
        c_impls.append(c_impl)
    tabla(L, ["linea", "n entran a reach", "N_conf min", "N_conf mediana",
              "N_conf max", "N_conf <= 3", "azar = media(0.96/N_conf)",
              "nuestro b (bundle)", "b recomputado (media por pert)",
              "b - azar", "b / azar", "c implicado por nuestro b",
              "reach de nuestro suelo"], rows)
    L += [
        "La columna `b recomputado` esta como control: es la media de los `reach` por",
        "perturbacion del run del `b`, y debe coincidir con el `b` del bundle, que es",
        "`agg = mean` sobre las que entran. La ultima columna es nuestro suelo sobre la",
        "MISMA cohorte -- §2.2 midio que `n_conf` es identico entre los dos objetos, asi",
        "que las dos medias se calculan sobre las mismas perturbaciones.",
        "",
        "### 3.3 Veredicto de Hueco 3",
        "",
    ]
    razones = [(ln, b_of / float(az.mean())) for ln, _, _, az, b_of in det]
    peor = max(razones, key=lambda t: abs(np.log(t[1])))
    arriba = [ln for ln, r in razones if r > 1.1]
    abajo = [ln for ln, r in razones if r < 0.9]
    cerca = [ln for ln, r in razones if 0.9 <= r <= 1.1]
    azs = np.array([float(az.mean()) for _, _, _, az, _ in det])
    bs = np.array([b_of for _, _, _, _, b_of in det])
    corr_ab = float(np.corrcoef(azs, bs)[0, 1])
    L += [
        "**No coinciden, y la discrepancia no tiene un signo unico.** Linea por linea, el",
        "cociente `b / azar`:",
        "",
    ]
    for ln, r in sorted(razones, key=lambda t: t[1]):
        donde = ("por DEBAJO" if r < 0.95 else
                 ("practicamente ENCIMA" if r <= 1.05 else "por ENCIMA"))
        L.append(f"- {ln}: **{r:.2f}x** -- {donde} del azar")
    L += [
        "",
        f"Asi que **{len([1 for _, r in razones if r > 1.05])} por encima, "
        f"{len([1 for _, r in razones if 0.95 <= r <= 1.05])} encima y "
        f"{len([1 for _, r in razones if r < 0.95])} por debajo.** No hay una direccion",
        "unica que atribuir a nada.",
        "",
        "Pero lo decisivo no es el signo: es la **amplitud**. El azar esperado varia",
        f"**{azs.max()/azs.min():.0f}x** entre lineas ({azs.min():.4f} en MCF7 a "
        f"{azs.max():.4f} en K562), porque la distribucion de `N_conf` lo es. El `b` varia",
        f"solo **{bs.max()/bs.min():.1f}x** ({bs.min():.4f} a {bs.max():.4f}).",
        "",
        f"La correlacion entre los dos es **{corr_ab:+.3f}**, asi que el `b` se mueve en la",
        "MISMA direccion que la ley -- no es independiente de ella. Lo que no hace es",
        "moverse por lo mismo: un objeto que estuviera en el punto sin habilidad tendria que",
        "escalar con `1/N_conf`, y un cambio de 39x en el azar le produce 4x. **La forma de",
        "la dependencia es otra**, y la seccion siguiente mide cual.",
        "",
        "### 3.3bis De que esta hecho nuestro `b` de `reach`",
        "",
        "La tabla por perturbacion dice de donde sale, y no hace falta inventar nada:",
        "",
    ]
    rows = []
    for ln, dentro, nc, az, b_of in det:
        vals = np.array([v["reach"] for v in dentro.values()])
        ncv = np.array([v["n_conf"] for v in dentro.values()])
        unos = vals == 1.0
        resto = float(vals[~unos].mean()) if (~unos).any() else float("nan")
        rows.append([
            ln, len(dentro), int(unos.sum()),
            sorted(int(x) for x in ncv[unos]) if unos.any() else "--",
            f"{vals.sum():.4f}", f"{vals[unos].sum():.4f}",
            f"{vals[unos].sum()/vals.sum():.3f}" if vals.sum() else "n/d",
            f(b_of), f(resto), int((ncv <= 3).sum()), int(ncv.min()),
        ])
    tabla(L, ["linea", "n entran", "perts con reach = 1.0 exacto",
              "sus N_conf", "suma de reach", "suma de los 1.0",
              "fraccion del total", "b", "media sin los 1.0",
              "perts con N_conf <= 3", "N_conf min"], rows)
    L += [
        "**Entre el 53% y el 88% del `b` de `reach` lo ponen entre 4 y 7 perturbaciones que",
        "sacan 1.0 exacto, y todas tienen `N_conf <= 10`** -- la mayoria `<= 3`. Sin ellas,",
        "la media del resto cae a 0.05-0.20.",
        "",
        "El codigo describe exactamente esa banda y dice que en la val oficial A **esta",
        "vacia**: *\"val A's smallest N_conf is 4, so the `N_conf <= 3` band where the head",
        "flip is worth 0.29-0.51 is empty\"* (`direction.py:977-979`). En nuestras lineas no",
        "lo esta: hay 8, 2, 6, 10, 12 y 0 perturbaciones en esa banda. **MCF7 es la unica",
        "con la banda vacia** (`N_conf` minimo 11) y es tambien la unica sin ninguna",
        "perturbacion en 1.0, y la de `b` mas bajo (0.1019).",
        "",
        "Los dos regimenes, dichos aparte porque son distintos:",
        "",
        "- **Lineas con la banda poblada** (las otras cinco): el `b` lo decide el CONTEO de",
        "  perturbaciones que resuelven del todo, no el presupuesto. Por eso es plano.",
        "- **MCF7, banda vacia**: ninguna resuelve del todo, y aun asi el `b` (0.1019) queda",
        "  casi 10x por encima del azar (0.0103). La caida `1/N_conf` de la ley es mucho mas",
        "  empinada que lo que el `b` consigue de verdad sobre perturbaciones con `N_conf`",
        "  de tres digitos.",
        "",
        "En **ninguna** linea el `b` se apoya sobre `c/N_conf`, y en las dos puntas la razon",
        "medida es distinta. No se propone ninguna causa comun: no esta en el codigo.",
        "",
        "### 3.3ter Lo que eso dice",
        "",
        "- El extremo 0 de la escala de `reach` **no** es el piso de azar sobre nuestra",
        f"  poblacion. En {len([1 for _, r in razones if r > 1.05])} lineas el `b` esta por",
        "  encima, asi que un predictor de moneda cobraria `s < 0`; en A549 esta un 6% por",
        "  debajo, asi que ahi cobraria `s > 0`; en K562 cae encima por 1%. Lo que no hay",
        "  es coincidencia sistematica.",
        "- La cifra con que el codigo ilustra la coincidencia (*\"the baseline arm reads",
        "  0.042314 against an ESTIMATED no-skill 0.0378, so the baseline already sits at",
        "  the no-skill point\"*, `direction.py:971-972`) es de la val oficial A, cuya",
        "  distribucion de `N_conf` es otra: *\"val A's smallest N_conf is 4\"*",
        "  (`direction.py:978`). En nuestras lineas la mediana de `N_conf` es de un digito",
        "  en varias y el minimo es 1. La ley `c/N_conf` es justamente la que el codigo",
        "  declara sobreestimada en esa zona, asi que el azar que calculamos con ella es si",
        "  acaso un techo del azar, y nuestro `b` queda aun mas arriba.",
        "- La columna `c implicado` resuelve `media(c/N_conf) = b` por `c`. Es la lectura",
        "  no-inventiva de la discrepancia: nuestro `b` se comporta como un predictor con",
        f"  `c` entre {min(c_impls):.2f} y {max(c_impls):.2f}"
        ", contra el 0.96 de la moneda y el 1.90",
        "  del que *\"ranks a single HABITUAL-DIRECTION gene first\"* (`direction.py:967`).",
        "  El codigo dice de ese segundo predictor que es *\"close to free\"*. **No se afirma",
        "  que el `b` sea ese predictor**: se dice donde cae su `c` implicado y de donde",
        "  salen los dos puntos de referencia.",
        "",
        "### 3.4 Las filas trabajadas de Hueco 3",
        "",
        "Los intermedios crudos: `N_conf` por perturbacion, el azar que implica, y el",
        "`reach` que el `b` saco de verdad.",
        "",
        "**Criterio de seleccion, declarado antes de mirar los valores**: la PRIMERA y la",
        "ULTIMA linea en orden alfabetico -- A549 y MCF7 -- y de cada una **todas** las",
        "perturbaciones que entran a `reach`, sin truncar. Es un criterio posicional y da",
        "los dos conjuntos completos, no una muestra. Las filas van ordenadas por `N_conf`",
        "creciente, que es una funcion del lado real y por tanto invariante a la submission",
        "(§2.2).",
        "",
    ]
    rows = []
    for objetivo in (det[0], det[-1]):
        lnx, dentro, _, _, _ = objetivo
        for p, v in sorted(dentro.items(), key=lambda kv: (kv[1]["n_conf"], kv[0])):
            az = C_FIT / v["n_conf"]
            rows.append([
                lnx, p, v["n_conf"], f"{az:.6f}", f(v["reach"]),
                f"{v['reach'] - az:+.6f}",
                "**si**" if v["reach"] == 1.0 else "",
                f(v["q"], 4) if v["q"] is not None else "n/d",
                v["n_DE_real"], v["n_DE_pred"], v["n_DE_union"],
            ])
    tabla(L, ["linea", "perturbacion", "N_conf", "azar 0.96/N_conf",
              "reach del b", "reach - azar", "reach = 1.0", "q", "n_DE_real",
              "n_DE_pred", "\\|union\\|"], rows)
    n_a = sum(1 for r in rows if r[0] == det[0][0])
    L += [
        f"Son **{len(rows)}** filas: {n_a} de {det[0][0]} y {len(rows)-n_a} de "
        f"{det[-1][0]}, los dos conjuntos enteros.",
        "",
        "Tres cosas se leen directamente en ellas:",
        "",
        "1. La columna `azar` es `0.96/N_conf` celda a celda, y su media sobre cada linea",
        "   es el `azar` de la tabla de §3.2.",
        "2. `n_DE_pred = 0` en todas las filas de A549 -- el `b` no predice ni un gen",
        "   significativo -- y aun asi cinco de ellas sacan `reach = 1.0`. Es la propiedad",
        "   que `descomposicion.md` ya habia explicado: con `universe='adjudicated'` el",
        "   ranking recorre solo los genes significativos de la REFERENCIA",
        "   (`direction.py:980-982`), asi que el conjunto predicho vacio no impide nada.",
        "3. `n_DE_real` y `|union|` coinciden en todas las filas de A549, por lo mismo: la",
        "   union es el lado real solo. Y `n_DE_real` coincide con `N_conf`, que es la",
        "   identidad que §2.4 comprueba.",
        "",
    ]


# ------------------------------------------------------- hueco 4
def hueco4(L, pdsj, coh):
    L += [
        "---",
        "",
        "## Hueco 4 — `pds`: la causa en HAP1 y BXPC3",
        "",
        "### 4.1 La hipotesis, declarada antes de medir",
        "",
        "> El estadistico asigna exactamente 0.5 a una fila empatada",
        "> (`discrimination.py:70`, `tie_policy=\"midrank\"`). Nuestro suelo esta MENOS",
        "> disperso, asi que produciria mas casi-empates y se pegaria mas a 0.5; el `b`, mas",
        "> disperso, se alejaria de 0.5 en una direccion u otra segun la linea.",
        "",
        "### 4.2 La distribucion por perturbacion",
        "",
        "`pds_cosine` emite una fila por perturbacion (`agg=\"mean\"`, `catalog.py:319-324`).",
        "Las 33 de cada linea, para los dos objetos:",
        "",
    ]
    rows = []
    for ln in LINEAS:
        if ln not in pdsj:
            continue
        d = pdsj[ln]
        rows.append([
            ln, d["n"], f(d["media_u"]), f(d["se_u"] * np.sqrt(d["n"])),
            f(d["media_b"]), f(d["se_b"] * np.sqrt(d["n"])),
            f"{d['dif']:+.6f}", f"{d['corr']:+.4f}",
        ])
    tabla(L, ["linea", "n perts", "media suelo", "sd suelo", "media b", "sd b",
              "b - suelo", "corr(suelo, b) por pert"], rows)
    # fracciones exactas en 0.5
    L += [
        "Y la fraccion que cae **exactamente** en 0.5, que es lo que la hipotesis predecia:",
        "",
    ]
    rows = []
    for ln in LINEAS:
        if ln not in coh:
            continue
        u = [v["pds"] for v in coh[ln]["objetos"]["suelo"]["por_pert"].values()
             if v["pds"] is not None]
        b = [v["pds"] for v in coh[ln]["objetos"]["b_oficial"]["por_pert"].values()
             if v["pds"] is not None]
        nu = sum(1 for x in u if x == 0.5)
        nb = sum(1 for x in b if x == 0.5)
        rows.append([ln, len(u), nu, f"{nu/len(u):.4f}", len(b), nb, f"{nb/len(b):.4f}"])
    tabla(L, ["linea", "n suelo", "exactamente 0.5 (suelo)", "fraccion",
              "n b", "exactamente 0.5 (b)", "fraccion"], rows)
    tot_u = sum(r[2] for r in rows)
    tot_b = sum(r[5] for r in rows)
    L += [
        f"**{tot_u} de {sum(r[1] for r in rows)} en nuestro suelo, {tot_b} de "
        f"{sum(r[4] for r in rows)} en el `b`.** El `b` tiene MAS celdas en 0.5 exacto que",
        "nuestro suelo. La hipotesis predecia lo contrario.",
        "",
        "### 4.3 Dos razones por las que el conteo de 0.5 no era la prueba que parecia",
        "",
        "**Primera: 0.5 es un punto de la rejilla.** `vcc2026` usa",
        "`rank_denominator: n-1` (`configs/vcc2026.yaml:33`), asi que con n = 33",
        "perturbaciones el valor es `1 - rango/32` y vive en la rejilla",
        "`{0, 1/32, ..., 32/32}` -- 33 puntos, y `16/32 = 0.5` es uno de ellos. Verificado:",
        "todos los valores de los dos objetos en las seis lineas son multiplos exactos de",
        "1/32. Una fila SIN ningun empate puede caer en 0.5 por rango 16. Contar 0.5",
        "exactos no cuenta empates.",
        "",
        "**Segunda: el unico mecanismo de empate conocido no se activa.** El codigo nombra",
        "el caso (`discrimination.py:63-68`): *\"under `cosine` a zero-norm predicted effect",
        "(a submission pasting the reference control cells for a target) ties the entire row",
        "at distance 1.0\"*, y `midrank` *\"scores every such target exactly 0.5\"*",
        "(`discrimination.py:69`). `descomp_pds.json` ya midio `frac_norma_cero` en los tres",
        "objetos y las seis lineas: **0 en todos**. Ningun delta tiene norma cero, asi que",
        "ninguna fila esta empatada entera. Sobre flotantes continuos, un empate exacto",
        "entre dos distancias cosine distintas tiene probabilidad practicamente nula.",
        "",
        "### 4.4 Entonces de que esta hecha la variacion: el yardstick correcto",
        "",
        "Lo que la tabla de §4.2 muestra es que **ninguno de los dos objetos esta",
        "concentrado en 0.5**. Los dos estan repartidos casi uniformemente sobre la rejilla:",
        "",
        f"- `sd` teorico de la rejilla uniforme de 33 puntos: **{np.std(np.arange(33)/32.0):.6f}**",
        "- `sd` medido: suelo "
        + " / ".join(f"{pdsj[ln]['se_u']*np.sqrt(pdsj[ln]['n']):.4f}" for ln in LINEAS if ln in pdsj)
        + "; `b` "
        + " / ".join(f"{pdsj[ln]['se_b']*np.sqrt(pdsj[ln]['n']):.4f}" for ln in LINEAS if ln in pdsj),
        "",
        "Eso es exactamente lo que el codigo dice que debe pasar con una emision sin",
        "informacion: *\"the ranks are a bijection onto `{0..n-1}`, so the panel mean is",
        "0.5\"* (`discrimination.py:73-74`), y `catalog.py:313-314`: *\"Its comparator sits at",
        "the metric's chance point (a uniform rank is 0.5)\"*. La media de panel sale ~0.5",
        "**no porque las filas se peguen a 0.5, sino porque los rangos se reparten**.",
        "",
        "Y una media de 33 rangos casi uniformes tiene ruido de muestreo propio:",
        "",
        f"```",
        f"SE(media de panel) = sd_rejilla / sqrt(33) = {np.std(np.arange(33)/32.0):.6f} / "
        f"{np.sqrt(33):.6f} = {np.std(np.arange(33)/32.0)/np.sqrt(33):.6f}",
        f"```",
        "",
        "**Ese** es el yardstick de `|b - u|`, no el `sd` entre particiones del anchor. Con",
        "el `SE` de la diferencia calculado sobre las diferencias por perturbacion",
        "observadas (que no asume independencia):",
        "",
    ]
    rows = []
    for ln in LINEAS:
        if ln not in pdsj:
            continue
        d = pdsj[ln]
        rows.append([
            ln, f"{d['dif']:+.6f}", f(d["se_u"], 4), f(d["se_b"], 4),
            f(d["se_dif_indep"], 4), f(d["se_dif_obs"], 4),
            f"**{d['z_obs']:.2f}**", f"{d['anchor_sd']:.4f}",
            f"{abs(d['dif'])/d['anchor_sd']:.2f}",
        ])
    tabla(L, ["linea", "b - u", "SE media u", "SE media b", "SE dif (indep)",
              "SE dif (observado)", "\\|b-u\\| / SE dif",
              "sd entre particiones de r", "\\|b-u\\| / sd del anchor"], rows)
    zs = [pdsj[ln]["z_obs"] for ln in LINEAS if ln in pdsj]
    L += [
        f"**Ninguna linea pasa de {max(zs):.2f} SE.** Las dos que \"excedian\" el `sd` del",
        "anchor caen en "
        + ", ".join(f"{ln} {pdsj[ln]['z_obs']:.2f}" for ln in ("HAP1", "BXPC3") if ln in pdsj)
        + " contra el yardstick correcto. Las dos ultimas columnas muestran la diferencia",
        "entre los dos yardsticks: el del anchor declara "
        + " y ".join(f"{abs(pdsj[ln]['dif'])/pdsj[ln]['anchor_sd']:.1f}"
                     for ln in ("BXPC3", "HAP1") if ln in pdsj)
        + " sigma donde el del estadistico declara "
        + " y ".join(f"{pdsj[ln]['z_obs']:.1f}" for ln in ("BXPC3", "HAP1") if ln in pdsj)
        + ".",
        "",
        "**Por que el `sd` entre particiones del anchor era el yardstick equivocado**, con",
        "los numeros a la vista: las medias de panel de `pds` de las 5 particiones son",
        "",
    ]
    rows = []
    for ln in LINEAS:
        if ln not in pdsj:
            continue
        d = pdsj[ln]
        rows.append([ln, " / ".join(f"{v:.4f}" for v in d["anchor_vals"]),
                     f"{np.mean(d['anchor_vals']):.4f}", f"{d['anchor_sd']:.4f}",
                     f(d["media_u"], 4), f(d["media_b"], 4)])
    tabla(L, ["linea", "pds de las 5 particiones de r", "media", "sd",
              "media u", "media b"], rows)
    L += [
        "Las particiones del anchor son real-contra-real: sus valores por perturbacion se",
        "apilan cerca de 1.0, donde la rejilla deja poco sitio para moverse, y su `sd` de",
        "panel es chico por eso. Es el `sd` de otro regimen. Usarlo para juzgar la",
        "variacion de un objeto sin informacion, cuyos rangos se reparten por toda la",
        "rejilla, subestima el ruido por un factor de 2 a 3.",
        "",
        "### 4.5 Veredicto de Hueco 4",
        "",
        "**La hipotesis no se sostiene, en sus dos mitades.**",
        "",
        "1. El mecanismo de empate que proponia no se activa: ningun delta tiene norma cero",
        "   (`descomp_pds.json`, seis lineas, tres objetos), que es el unico caso de fila",
        "   entera empatada que el codigo nombra.",
        "2. La prediccion direccional sale al reves: el `b` tiene mas celdas en 0.5 exacto",
        f"   ({tot_b}) que nuestro suelo ({tot_u}) -- y de todos modos ese conteo no mide",
        "   empates, porque 0.5 esta en la rejilla.",
        "3. Nuestro suelo no esta \"menos disperso\" en el espacio del estadistico: su `sd`",
        "   por perturbacion es indistinguible del del `b`, y los dos coinciden con el de",
        "   la rejilla uniforme. La menor dispersion que `descomposicion.md` midio estaba en",
        "   el espacio de los DELTAS (distancia cosine entre pares), no en el de los RANGOS,",
        "   y el estadistico puntua rangos.",
        "",
        "**Lo que queda explicado**: el signo cambiante y los dos \"excesos\". Los dos objetos",
        "estan en el punto sin informacion del estadistico; sus medias de panel son medias",
        f"de 33 rangos casi uniformes, con un `SE` propio de ~{np.std(np.arange(33)/32.0)/np.sqrt(33):.3f};",
        "las seis diferencias caben en ese ruido. Nada que atribuir al escalado de",
        "`baseline.py:344`.",
        "",
        "**Lo que queda sin explicar**: nada de esta pregunta. Lo que no se midio aqui es",
        "si un candidato CON habilidad en `pds` tiene ese mismo `SE` -- no lo tendria: el",
        "regimen del anchor muestra que cuando los valores se apilan arriba el `sd` baja.",
        "Eso no era la pregunta y no se estima.",
        "",
        "### 4.6 Las 20 filas trabajadas de Hueco 4",
        "",
        "Criterio declarado antes de mirar los valores: las perturbaciones de **HAP1** -- la",
        "linea del `|b - u|` mas grande, que es la que el encargo nombra -- en orden",
        "alfabetico, las 20 primeras.",
        "",
    ]
    ln = "HAP1"
    rows = []
    if ln in coh:
        u = coh[ln]["objetos"]["suelo"]["por_pert"]
        b = coh[ln]["objetos"]["b_oficial"]["por_pert"]
        for p in sorted(u)[:20]:
            pu, pb = u[p]["pds"], b[p]["pds"]
            rows.append([
                ln, p, f(pu, 6), f"{round(pu*32)}/32" if pu is not None else "n/d",
                f(pb, 6), f"{round(pb*32)}/32" if pb is not None else "n/d",
                f"{pb-pu:+.6f}" if (pu is not None and pb is not None) else "n/d",
                "si" if pu == 0.5 else "", "si" if pb == 0.5 else "",
            ])
    tabla(L, ["linea", "perturbacion", "pds suelo", "rango suelo", "pds b",
              "rango b", "b - suelo", "suelo en 0.5", "b en 0.5"], rows)
    L += [
        f"Son {len(rows)} filas. Las columnas `rango` muestran que los valores son enteros",
        "sobre 32: la rejilla. La media de la columna `pds b` sobre las 33 (no solo estas",
        f"20) es {f(pdsj[ln]['media_b'])}, y la de `pds suelo` {f(pdsj[ln]['media_u'])}; su",
        f"diferencia {pdsj[ln]['dif']:+.6f} es la celda de HAP1 de la tabla de §4.4.",
        "",
    ]


def cierre(L):
    L += [
        "---",
        "",
        "## Lo que no se hizo",
        "",
        "- **No se redacto ningun criterio de adopcion.** No es de este encargo.",
        "- **No se construyo ningun candidato.** Las dos emisiones que se corrieron son",
        "  nuestro suelo (`predict_floor` sin modificar) y la del `b` oficial",
        "  (`cell-eval2 baseline --seed 0`), las dos ya conocidas.",
        "- **No se modifico ninguna metrica ni ningun arnes.** `build_fold_adatas`,",
        "  `predict_floor`, `contar_de_por_perturbacion`, `_disjoint_halves` y",
        "  `_reference_stats` se usaron tal cual. `descomp_rgen.py` -- script de analisis de",
        "  esta sesion, no un arnes -- gano un argumento de linea de comandos y una columna",
        "  de conteos por perturbacion; su metodo no cambio, y §1.1 lo verifica contra los",
        "  valores que ya habia estampado el bundle.",
        "- **No se selecciono ninguna celda por su valor ni por su signo.** Los criterios de",
        "  seleccion que aparecen (las filas trabajadas de §3.4 y §4.6, y las lineas de cada",
        "  tabla) estan declarados antes de los numeros y son posicionales o funcion del",
        "  lado real. §3.4 muestra los dos conjuntos ENTEROS, sin truncar.",
        "- **No se derivo ningun objetivo numerico por linea.**",
        "- No se leyo nada de `~/vcc/benchmark/` ni de",
        "  `/datos/benchmark_VEDADO_techo_metricas/`.",
        "- Sin envios, sin `vcc sample`, sin commits.",
        "",
    ]


def main():
    rgen = json.load(open(RGEN))
    coh = json.load(open(COH))
    pdsj = json.load(open(PDSJ))
    nat = leer_jsonl(NATIVAS)

    faltan = [ln for ln in LINEAS if ln not in rgen or ln not in coh]
    if faltan:
        sys.stderr.write(f"ATENCION: faltan datos de {faltan}\n")

    L = []
    cabecera(L)
    resumen(L)
    verif_b(L, coh)
    hueco1(L, rgen, nat)
    hueco1_span(L, rgen, coh, nat)
    hueco2(L, coh, rgen)
    hueco3(L, coh, nat)
    hueco4(L, pdsj, coh)
    cierre(L)

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    sys.stderr.write(f"escrito {OUT_MD} ({len(L)} lineas)\n")


if __name__ == "__main__":
    main()
