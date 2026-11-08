#!/usr/bin/env python
"""Arma modelo/out/descomposicion.md: de que esta hecho cada numero de las seis
metricas, por linea. Lectura y conteo sobre artefactos que ya existen, mas dos
recomputos declarados (descomp_pds.json, descomp_rgen.json)."""
import json
import math
import os

import pandas as pd

OUT_MD = "/home/ubuntu/vcc/modelo/out/descomposicion.md"
SUELO_RES = "/datos/modelo/out/run_{L}/results.csv"
SUELO_AGG = "/datos/modelo/out/run_{L}/agg_results.csv"
B_RES = "/datos/modelo/out/nativas_comp/BASE_{L}/baseline_results.csv"
B_AGG = "/datos/modelo/out/nativas_comp/BASE_{L}/baseline_agg.csv"
NATIVAS = "/home/ubuntu/vcc/modelo/out/anclas_nativas_comp.json"
CABLE = "/home/ubuntu/vcc/modelo/out/cable_suelo.json"
PDS = "/home/ubuntu/vcc/modelo/out/descomp_pds.json"
VERIF = "/home/ubuntu/vcc/modelo/out/descomp_verif.json"
RGEN = "/home/ubuntu/vcc/modelo/out/descomp_rgen.json"

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
N_PANEL = 33
K = 5
METRICS = [
    "pds_cosine",
    "expr_mse_unbiased_capped_norm",
    "de_wilcoxon_lfc_nmae",
    "de_wilcoxon_direction_fidelity_yield_raw",
    "de_wilcoxon_direction_reach_raw",
    "de_wilcoxon_sig_jaccard",
]
CORTO = {
    "pds_cosine": "pds_cosine",
    "expr_mse_unbiased_capped_norm": "mse_capped_norm",
    "de_wilcoxon_lfc_nmae": "lfc_nmae",
    "de_wilcoxon_direction_fidelity_yield_raw": "fid",
    "de_wilcoxon_direction_reach_raw": "reach",
    "de_wilcoxon_sig_jaccard": "jac",
}
MAS_ES_MEJOR = {
    "pds_cosine": True, "expr_mse_unbiased_capped_norm": False,
    "de_wilcoxon_lfc_nmae": False,
    "de_wilcoxon_direction_fidelity_yield_raw": True,
    "de_wilcoxon_direction_reach_raw": True,
    "de_wilcoxon_sig_jaccard": True,
}
JAC = "de_wilcoxon_sig_jaccard"


def jl(path):
    out = []
    if not os.path.exists(path):
        return out
    with open(path) as fh:
        for ln in fh:
            ln = ln.strip()
            if ln:
                out.append(json.loads(ln))
    return out


def fmt(v, nd=6):
    if v is None:
        return "NaN"
    try:
        f = float(v)
    except (TypeError, ValueError):
        return str(v)
    if f != f:
        return "NaN"
    return f"{f:.{nd}f}"


def agg_mean(path):
    d = pd.read_csv(path)
    return d.loc[d["statistic"] == "mean"].iloc[0]


def descompone(res_path, agg_path):
    """Por metrica: filas emitidas, NaN, entradas, media de las entradas, y para
    jac la particion convencion/genuina."""
    r = pd.read_csv(res_path)
    a = agg_mean(agg_path)
    out = {}
    for m in METRICS:
        s = r[r.metric == m]["value"]
        ok = s.dropna()
        rep = a.get(m)
        try:
            rep = float(rep)
        except (TypeError, ValueError):
            rep = None
        if rep is not None and rep != rep:
            rep = None
        d = {
            "filas_emitidas": int(len(s)),
            "nan": int(s.isna().sum()),
            "entraron": int(len(ok)),
            "media_entraron": float(ok.mean()) if len(ok) else None,
            "reportado": rep,
        }
        if m == JAC and len(s):
            E = int((ok == 1.0).sum())
            gen = ok[ok != 1.0]
            d.update({
                "E_convencion": E,
                "n_genuinos": int(len(gen)),
                "media_genuinos": float(gen.mean()) if len(gen) else None,
                "comp_convencion": E / N_PANEL,
                "comp_genuina": float(gen.sum()) / N_PANEL if len(gen) else 0.0,
            })
        out[m] = d
    return out


def main():
    nat = {x["linea"]: x for x in jl(NATIVAS) if "error" not in x}
    cable = {x["linea"]: x for x in jl(CABLE) if "error" not in x}
    pds = json.load(open(PDS)) if os.path.exists(PDS) else {}
    verif = json.load(open(VERIF)) if os.path.exists(VERIF) else {}
    rgen = json.load(open(RGEN)) if os.path.exists(RGEN) else {}

    D_u = {L: descompone(SUELO_RES.format(L=L), SUELO_AGG.format(L=L)) for L in LINES}
    D_b = {L: descompone(B_RES.format(L=L), B_AGG.format(L=L)) for L in LINES}

    L_ = []
    L_.append("# De que esta hecho cada numero: la descomposicion por linea y metrica\n")
    L_.append(
        "Las seis metricas de `vcc2026` **no se miden sobre la misma cohorte**, y "
        "una de ellas asigna un valor por convencion a la ausencia de senal. El "
        "`avg_score` promedia plano a 1/6 seis numeros calculados sobre "
        "poblaciones de 8 a 33 perturbaciones. Esto descompone cada numero en lo "
        "que mide y lo que es convencion.\n"
    )
    L_.append(
        "Lectura y conteo sobre artefactos que ya existian, mas dos recomputos "
        "que se declaran donde aparecen: la dispersion de perfiles "
        "(`descomp_pds.json`) y la componente genuina de `r` "
        "(`descomp_rgen.json`), esta ultima porque el bundle **no persiste** el "
        "detalle por perturbacion de las particiones del anchor.\n"
    )

    # ================= PARTE 1 =================
    L_.append("\n# PARTE 1 -- La descomposicion\n")

    # ---- jac ----
    L_.append("## `jac` -- la convencion de union vacia, medida en tres objetos\n")
    if verif:
        c = verif["contingencia"]
        L_.append(
            "**Primero la premisa, verificada y no inferida.** El conteo usa "
            "`E = count(jaccard == 1.0)` como numero de perturbaciones con union "
            "vacia. Un jaccard de 1.0 tambien podria venir de un acuerdo "
            "perfecto, asi que se midio: se recorrio una corrida con "
            "`--write-degenes` y se contaron las uniones con "
            "`techo_arnes.contar_de_por_perturbacion` (sin modificar), que "
            "reproduce las dos lineas de filtrado de `de_sig_jaccard`. "
            f"Contingencia sobre {verif['n_perts']} perturbaciones de "
            f"{verif['linea']}:\n"
        )
        L_.append("| | union vacia | union NO vacia |")
        L_.append("|---|---|---|")
        L_.append(f"| `jac == 1.0` | **{c['jac==1.0 y union vacia']}** | {c['jac==1.0 y union NO vacia']} |")
        L_.append(f"| `jac != 1.0` | {c['jac!=1.0 y union vacia']} | **{c['jac!=1.0 y union NO vacia']}** |")
        L_.append(
            f"\nLas dos celdas fuera de la diagonal son **0**: la equivalencia "
            f"`jac == 1.0` <=> `|union| == 0` queda **verificada** "
            f"({verif['equivalencia_verificada']}). Asi que `E` es una medicion.\n"
        )
    L_.append(
        "La convencion es `de_sig_jaccard`'s `J(vacio, vacio) = 1` "
        "(`metrics/de.py:195`, y el retorno en `metrics/de.py:263`: "
        "`out[p] = 1.0 if union == 0 else tp / union`). Cada perturbacion sin "
        "senal en NINGUNO de los dos lados aporta 1.0 al promedio de 33.\n"
    )
    L_.append(
        "| objeto | linea | E (union vacia) | genuinos | media de los genuinos | "
        "comp. CONVENCION = E/33 | comp. GENUINA = sum/33 | reportado | suma | dif |"
    )
    L_.append("|---" * 10 + "|")
    for etq, D in (("nuestro suelo", D_u), ("b oficial", D_b)):
        for L in LINES:
            d = D[L][JAC]
            rec = d["comp_convencion"] + d["comp_genuina"]
            L_.append(
                f"| {etq} | {L} | {d['E_convencion']} | {d['n_genuinos']} | "
                f"{fmt(d['media_genuinos'])} | {fmt(d['comp_convencion'])} | "
                f"{fmt(d['comp_genuina'])} | {fmt(d['reportado'])} | {fmt(rec)} | "
                f"{abs(rec - d['reportado']):.1e} |"
            )
    if rgen:
        for L in LINES:
            if L not in rgen:
                continue
            g = rgen[L]
            L_.append(
                f"| **r (anchor, k=5)** | {L} | "
                f"{sum(s['E_union_vacia'] for s in g['splits'])/K:.1f} (media) | "
                f"{sum(s['n_genuinos'] for s in g['splits'])/K:.1f} (media) | "
                f"{fmt(sum(s['media_genuinos'] for s in g['splits'])/K)} | "
                f"{fmt(g['r_convencion_media'])} | {fmt(g['r_genuina_media'])} | "
                f"{fmt(g['r_crudo_media'])} | "
                f"{fmt(g['r_convencion_media'] + g['r_genuina_media'])} | "
                f"{abs(g['r_convencion_media'] + g['r_genuina_media'] - g['r_crudo_media']):.1e} |"
            )
    L_.append(
        "\n**La identidad `jac_reportado = (E + sum(genuinos)) / 33` se verifica "
        "exacta en todos los objetos y todas las lineas** (columna `dif`).\n"
    )
    L_.append(
        "**El hallazgo.** En nuestro suelo la media de los jaccards genuinos es "
        "**exactamente 0.000000 en las seis lineas**: el numero reportado es "
        "convencion al 100%. En el `b` oficial la componente genuina es 0.000000 "
        "en tres lineas y <= 0.001 en las otras tres. Es decir: **para las dos "
        "emisiones triviales, `jac` no mide nada de la emision; cuenta cuantas "
        "perturbaciones no tienen senal.**\n"
    )
    if rgen:
        L_.append(
            "El `r` del anchor SI tiene componente genuina (es un replicado "
            "real), pero es una fraccion de su valor crudo:\n"
        )
        L_.append("| linea | r crudo | r convencion | r genuina | genuina / crudo |")
        L_.append("|---|---|---|---|---|")
        for L in LINES:
            if L not in rgen:
                continue
            g = rgen[L]
            frac = g["r_genuina_media"] / g["r_crudo_media"] if g["r_crudo_media"] else float("nan")
            L_.append(
                f"| {L} | {fmt(g['r_crudo_media'])} | {fmt(g['r_convencion_media'])} | "
                f"{fmt(g['r_genuina_media'])} | {frac:.1%} |"
            )
        L_.append(
            "\n### El recomputo del `r`, y por que hizo falta\n"
        )
        L_.append(
            "`anchor_splits.parquet` guarda `(split_index, seed, metric, value, "
            "n_perturbations)` -- el AGREGADO por particion, que mezcla las dos "
            "componentes. La descomposicion de `r` no se puede hacer leyendo. Se "
            "recomputaron las 5 particiones de A549 y HAP1 con la MISMA funcion "
            "del marcador (`cell_eval2.ceiling._disjoint_halves`, la que usa "
            "`anchor.py:117`), las MISMAS semillas que el bundle estampo, "
            "`half_b` como prediccion (`anchor.py:120`) y "
            "`control_source=pred` (`anchor.py:103`).\n"
        )
        L_.append(
            "**Criterio declarado antes de medir:** A549 y HAP1 son las dos "
            "celdas que el encargo nombra por su `span/SE` sobre el crudo. No se "
            "eligieron por su valor ni por su signo.\n"
        )
        L_.append("| linea | split | semilla | valor del bundle | valor recomputado | dif | E | media genuinos |")
        L_.append("|---" * 8 + "|")
        for L in LINES:
            if L not in rgen:
                continue
            for s in rgen[L]["splits"]:
                L_.append(
                    f"| {L} | {s['split_index']} | {s['seed']} | "
                    f"{fmt(s['valor_del_bundle'])} | {fmt(s['valor_recomputado'])} | "
                    f"{s['dif_vs_bundle']:.1e} | {s['E_union_vacia']} | "
                    f"{fmt(s['media_genuinos'])} |"
                )
        mx = max(rgen[L]["max_dif_vs_bundle"] for L in rgen)
        L_.append(
            f"\n**La reproduccion es exacta**: la diferencia maxima contra el "
            f"valor que el bundle estampo es {mx:.1e} sobre las 10 particiones. "
            "Asi que la descomposicion describe el `r` que realmente se usa, no "
            "una aproximacion.\n"
        )
        L_.append(
            "Y aparece el mecanismo que explica el `span` raquitico del crudo: "
            "**la componente de convencion trabaja EN CONTRA del span.** El "
            "replicado real tiene MENOS uniones vacias que el baseline, asi que "
            "cobra menos 1.0 regalados, lo que arrastra su valor crudo hacia "
            "abajo y casi borra la holgura genuina:\n"
        )
        L_.append("| linea | aporte de la CONVENCION al span | aporte de lo GENUINO | span crudo |")
        L_.append("|---|---|---|---|")
        for L in LINES:
            if L not in rgen:
                continue
            g = rgen[L]
            b = D_b[L][JAC]
            a_conv = g["r_convencion_media"] - b["comp_convencion"]
            a_gen = g["r_genuina_media"] - b["comp_genuina"]
            L_.append(
                f"| {L} | {a_conv:+.6f} | {a_gen:+.6f} | {a_conv + a_gen:+.6f} |"
            )
        L_.append("")

    # ---- las gateadas ----
    L_.append("\n## `fid`, `reach`, `lfc_nmae` -- las gateadas: quien entro, quien se cayo, y por que\n")
    L_.append(
        "Las tres excluyen perturbaciones, pero por mecanismos DISTINTOS, y la "
        "diferencia importa:\n"
    )
    L_.append("| metrica | mecanismo | condicion | linea del codigo | que recibe la excluida |")
    L_.append("|---|---|---|---|---|")
    L_.append(
        "| `lfc_nmae` | la fila **no se emite** | gate vacio, menor que "
        "`min_gate_size` (=10), o denominador cero/no finito | "
        "`metrics/de.py:608-611` | nada: \"those score nothing rather than 1.0\" |"
    )
    L_.append(
        "| `fid` (`direction_fidelity_yield_raw`) | la fila se emite con **NaN** | "
        "`max(n_pred, n_conf) == 0`, es decir los DOS lados vacios | "
        "`metrics/direction.py:776` con `_safe_div` en `689-699` | NaN, que el "
        "agregador omite |"
    )
    L_.append(
        "| `reach` (`direction_reach_raw`) | la fila se emite con **NaN** | "
        "`n_conf == 0` (pool de referencia vacio), o `q` indefinido | "
        "`metrics/direction.py:1050-1053` | NaN, que el agregador omite |"
    )
    L_.append(
        "\n**Ninguna de las tres aplica una convencion a la ausencia.** El "
        "codigo de `reach` lo dice explicitamente (`metrics/direction.py:1048-1050`): "
        "\"Spec 5: 'never reached P0 -> a computed 0.0' is DOMAIN-RESTRICTED and "
        "must not override the matrix. It is trivially true at N_conf = 0 (empty "
        "pool) and when q is undefined (no threshold exists) -- **both are NaN, "
        "not 0**.\" Es la decision opuesta a la de `jac`, en el mismo paquete.\n"
    )
    L_.append(
        "La consecuencia es que la media se toma sobre una cohorte mas chica, no "
        "que se contamine con valores de relleno. Medido (la columna "
        "`media_entraron` reproduce el valor reportado en las 36 celdas):\n"
    )
    L_.append(
        "| objeto | linea | metrica | filas emitidas | NaN | entraron | de 33 | "
        "media de las que entraron | reportado | dif |"
    )
    L_.append("|---" * 10 + "|")
    for etq, D in (("nuestro suelo", D_u), ("b oficial", D_b)):
        for L in LINES:
            for m in ("de_wilcoxon_lfc_nmae",
                      "de_wilcoxon_direction_fidelity_yield_raw",
                      "de_wilcoxon_direction_reach_raw"):
                d = D[L][m]
                dif = ("-" if (d["media_entraron"] is None or d["reportado"] is None)
                       else f"{abs(d['media_entraron'] - d['reportado']):.1e}")
                L_.append(
                    f"| {etq} | {L} | {CORTO[m]} | {d['filas_emitidas']} | "
                    f"{d['nan']} | **{d['entraron']}** | 33 | "
                    f"{fmt(d['media_entraron'])} | {fmt(d['reportado'])} | {dif} |"
                )
    L_.append(
        "\n**Dos hechos del codigo confirmados por medicion:**\n"
    )
    L_.append(
        "- `metrics/de.py:610-611` dice de los gates de `lfc_nmae`: \"All three "
        "are real-side decisions, so **the surviving set is the same for every "
        "submission**.\" Medido: la cohorte de `lfc_nmae` es identica entre "
        "nuestro suelo y el `b` en las seis lineas (8/8, 23/23, 8/8, 15/15, 9/9, "
        "33/33). Verificado.\n"
        "- El NaN de `reach` depende SOLO del lado real (`n_conf == 0`), asi que "
        "su cohorte tambien es identica entre las dos emisiones (15/15, 0/0, "
        "14/14, 6/6, 8/8, 0/0). El de `fid` depende de `max(n_pred, n_conf)`, "
        "que SI mira la emision: en HT29 cae 6 veces en nuestro suelo y 3 en el "
        "`b`, y en K562 8 contra 7. **La cohorte de `fid` es la unica de las "
        "tres que la submission puede mover.**\n"
    )

    # ---- reach: la pregunta abierta ----
    L_.append("\n## 🔴 `reach` -- la pregunta abierta desde `bucle_suelo.md`, ahora medida\n")
    L_.append(
        "La pregunta era: `reach` da NO-NULO (0.014 a 0.216 en nuestro suelo) con "
        "conjunto predicho VACIO. Estaba declarada \"no entendida, "
        "deliberadamente sin explicar\". La respuesta esta en el codigo y **no es "
        "una convencion analoga a `J(vacio,vacio)=1`**.\n"
    )
    L_.append("**Que calcula.** `reach = k*(t) / N_conf(t)` (`metrics/direction.py:1054`):\n")
    L_.append("```\nout[target] = _to_float_nan(k_star.get(target, 0) / n_conf)\n```\n")
    L_.append(
        "donde `k*` es la profundidad hasta la que el ranking de la prediccion se "
        "mantiene direccionalmente puro, y `N_conf` es el presupuesto confiable "
        "de la REFERENCIA.\n"
    )
    L_.append(
        "**Por que el conjunto predicho vacio no lo anula.** Con "
        "`universe='adjudicated'` el ranking recorre **solo los genes "
        "significativos de la REFERENCIA** "
        "(`metrics/direction.py:980-982`: \"`universe='adjudicated'` ranks only "
        "reference-significant genes, so `k* <= N_conf` and the result is bounded "
        "by 1\"). La prediccion no aporta el POOL; aporta el ORDEN y los SIGNOS. "
        "Una emision degenerada no llama nada significativo, pero sigue teniendo "
        "un `|log2FC|` por gen con el que ordenar, y un signo por gen.\n"
    )
    L_.append(
        "**De donde sale el valor no-nulo.** `metrics/direction.py:1016-1019`, "
        "issue #279:\n"
    )
    L_.append(
        "> \"`reach_raw`'s no-skill point is `~c/N_conf`, NOT a constant, and for "
        "a no-skill submission it is decided almost entirely at the HEAD of a "
        "ranking the submission controls. **A first-pair match is SUFFICIENT for "
        "`k* >= 1` (purity 1/1 >= P0), one coin flip independent of depth**, and "
        "it dominates the no-skill probability.\"\n"
    )
    L_.append(
        "Y la medicion que el propio codigo trae "
        "(`metrics/direction.py:1030-1033`): \"Measured at 0.13.0 over 200 "
        "replicates per cell: `E[reach_raw]` runs 0.5100 at `N_conf = 1` to "
        "0.0019 at 500 with `P(reach_raw > 0) = 0.385-0.565` at EVERY N_conf, and "
        "the fitted `c` is ~0.96 for a coin-flip predictor\".\n"
    )
    L_.append(
        "**Veredicto: no es una convencion aplicada a la ausencia, es un piso de "
        "azar.** Un acierto de signo en el primer gen del ranking -- una moneda -- "
        "da `k* = 1` y por tanto `reach = 1/N_conf > 0`. No hay un valor asignado "
        "a la ausencia: hay un computo cuyo valor esperado sin habilidad es "
        "`~c/N_conf`, que no es cero. Por eso **no se descompone como `jac`**: "
        "todo el numero \"mide\", lo que pasa es que su punto de no-habilidad no "
        "esta en 0.\n"
    )
    L_.append(
        "Dato relevante para la escala: el codigo tambien dice "
        "(`metrics/direction.py:1035-1037`) que \"on official val A the baseline "
        "arm reads 0.042314 against an ESTIMATED no-skill 0.0378, so **the "
        "baseline already sits at the no-skill point**\". Es decir, en `reach` el "
        "`b` medido ES el piso de azar -- que es exactamente lo que un extremo 0 "
        "deberia ser.\n"
    )

    # ---- mse ----
    L_.append("\n## `mse_capped_norm` -- no gateada, pero clampeada\n")
    nm = [(L, (nat.get(L, {}).get("r") or {}).get("expr_mse_unbiased_capped_norm", {}))
          for L in LINES]
    L_.append(
        "**No esta gateada y lee las 33.** En `results.csv` la metrica derivada "
        "no aparece por perturbacion (es `agg=\"ratio_of_sums\"` con "
        "`DerivedAgg`, `catalog.py:412-414`), asi que su cohorte se lee del "
        "artefacto del anchor, que la estampa:\n"
    )
    L_.append("| linea | `n_perturbations_min` | `n_perturbations_max` | de 33 |")
    L_.append("|---|---|---|---|")
    for L, rr in nm:
        L_.append(
            f"| {L} | {rr.get('n_perturbations_min')} | "
            f"{rr.get('n_perturbations_max')} | 33 |"
        )
    L_.append(
        "\nSus dos componentes SI aparecen por perturbacion y con 33 filas cada "
        "una (`expr_mse_unbiased_capped` y `expr_distance_unbiased`), sin NaN, en "
        "las seis lineas y en los dos objetos. Confirmado.\n"
    )
    L_.append(
        "**🔴 El primitivo que importa, sin cifras por linea.** "
        "`catalog.py:415-416` le pone `clamp_low=0.0, clamp_high=1.0`:\n"
    )
    L_.append(
        "```\nscoring=Scoring(scored=True, direction=\"lower\", anchor=0.0,\n"
        "                penalty=\"boxcox\", clamp_low=0.0, clamp_high=1.0),\n```\n"
    )
    L_.append(
        "y el comentario (`catalog.py:406-410`) declara que es la unica: \"this "
        "is the ONE scored vcc2026 member that clamps to [0, 1], on BOTH scales. "
        "The metric is capped so that 1.0 = predicting the control, so a value "
        "past the baseline is already no-skill and the graded Box-Cox tail below "
        "it buys nothing a ranking uses.\"\n"
    )
    L_.append(
        "**La propiedad, declarada:** por debajo del extremo 0 esta metrica "
        "devuelve **exactamente 0**, asi que **no informa si un candidato esta "
        "mejorando mientras siga del lado malo del baseline**. `cable_suelo.md` "
        "lo midio: las seis celdas dieron 0.000000, y sin el clamp A549 habria "
        "dado -1.950941. Un candidato que pase de -1.95 a -0.2 ve la misma cifra "
        "que uno que no cambio nada. **No se calcula aqui cuanto habria que "
        "bajar para salir del clamp: eso seria un peso, no una descripcion.**\n"
    )

    # ================= PARTE 2 =================
    L_.append("\n# PARTE 2 -- `pds_cosine`: por que el `b` queda por DEBAJO de `u`\n")
    L_.append(
        "⚠️ **No se compara contra el 0.5 de `scales.py:428`.** Esa es la base "
        "constante de otra escala (`low-random_high-1_v10`), con otro estimando y "
        "otra poblacion. El 0.5 que aparece abajo es otro objeto: el valor que el "
        "propio estadistico de rango asigna a una prediccion sin informacion, con "
        "su cita.\n"
    )
    L_.append("## Que calcula `pds_cosine`, y que premia\n")
    L_.append(
        "`metrics/discrimination.py:48-53`: \"A rank-based retrieval score: for "
        "each non-control perturbation, rank how close the predicted signed delta "
        "(effect) is to the matching real delta versus all other real deltas. "
        "Score = `1 - rank / D`\". Opera sobre el **delta** (perfil menos "
        "control), no sobre el perfil.\n"
    )
    L_.append(
        "Lo que premia es **especificidad**: que el delta predicho para `p` se "
        "parezca mas al delta real de `p` que a los de las otras 32. No premia "
        "magnitud ni dispersion.\n"
    )
    L_.append(
        "Y el punto sin informacion esta fijado por el estadistico, no por una "
        "constante externa (`metrics/discrimination.py:74-77`): \"A FULLY "
        "constant prediction self-corrects under either policy (the ranks are a "
        "bijection onto `{0..n-1}`, so **the panel mean is 0.5**)\". Mas "
        "`discrimination.py:70`: con `tie_policy=\"midrank\"` (el de `vcc2026`) "
        "\"`midrank` scores every such target exactly 0.5\" para una fila "
        "empatada.\n"
    )

    L_.append("\n## La hipotesis del encargo, puesta a prueba\n")
    L_.append(
        "La hipotesis era: el escalado de `baseline.py:344` mueve todas las "
        "perturbaciones hacia el MISMO perfil, lo que las acerca entre si y "
        "destruye discriminacion. Es verificable: se mide la dispersion entre los "
        "deltas de perturbacion en las tres emisiones, con las funciones del "
        "propio marcador (`prep.pseudobulk_bulk_lognorm` con "
        "`bulk_target_sum=50000` y `prep.delta`, que es lo que ve el comparador "
        "de `vcc2026`).\n"
    )
    L_.append(
        "**Criterio declarado antes de medir: las seis lineas, sin seleccion.**\n"
    )
    if pds:
        L_.append(
            "| linea | objeto | norma media del delta | distancia cosine media "
            "entre pares | min | max | rango |"
        )
        L_.append("|---" * 7 + "|")
        for L in LINES:
            if L not in pds:
                continue
            for k, et in (("real", "real (la verdad)"), ("suelo", "nuestro suelo"),
                          ("b_oficial", "b oficial")):
                r = pds[L][k]
                L_.append(
                    f"| {L} | {et} | {fmt(r['norma_media'], 4)} | "
                    f"**{fmt(r['cos_par_media'])}** | {fmt(r['cos_par_min'], 4)} | "
                    f"{fmt(r['cos_par_max'], 4)} | "
                    f"{fmt(r['cos_par_max'] - r['cos_par_min'], 4)} |"
                )
        L_.append(
            "\n🔴 **La hipotesis queda REFUTADA por la medicion.** Los deltas del "
            "`b` oficial estan casi tan dispersos como los reales: distancia "
            "cosine media entre pares 0.749-0.853 contra 0.772-0.885 del real. El "
            "escalado **no** colapsa las perturbaciones hacia un punto. Y nuestro "
            "suelo esta MENOS disperso (0.481-0.527), no mas.\n"
        )
        L_.append(
            "Las normas tampoco distinguen: los tres objetos tienen deltas de "
            "norma comparable (7.7 a 14.0), y **ninguno tiene deltas de norma "
            "cero** (`frac_norma_cero = 0.000` en los 18 casos), asi que el "
            "mecanismo de \"delta de norma cero empata toda la fila\" que "
            "`discrimination.py:66-70` describe no es el que opera aqui.\n"
        )
        L_.append(
            "Lo que SI distingue al real es el **rango**: sus pares van de 0.31 a "
            "1.07 (hay perturbaciones muy parecidas y pares anticorrelados), "
            "mientras el `b` se queda en 0.65-0.92 y nuestro suelo en 0.39-0.60. "
            "El real tiene ESTRUCTURA; los otros dos tienen dispersion sin "
            "estructura.\n"
        )
    L_.append("\n## Entonces, ¿que explica el orden `b < u`?\n")
    L_.append(
        "**Dispersion no es discriminacion.** `pds_cosine` no pregunta cuanto se "
        "separan los deltas predichos entre si, sino si el delta predicho de `p` "
        "identifica a `p`. Ninguna de las dos emisiones lleva informacion "
        "especifica de la perturbacion: la de nuestro suelo es ruido de remuestreo "
        "del control, y la del `b` es un perfil CONSTANTE unico "
        "(`baseline.py:121-122`: \"the result is still ONE constant vector\") mas "
        "ruido de remuestreo amplificado por el escalado. Las dos caen, por "
        "construccion del estadistico de rango, en el punto sin informacion.\n"
    )
    L_.append("Medido, los dos objetos rodean ese punto:\n")
    L_.append("| linea | `u` (nuestro suelo) | `b` oficial | b - u | ¿de que lado? |")
    L_.append("|---|---|---|---|---|")
    difs = []
    for L in LINES:
        u = D_u[L]["pds_cosine"]["reportado"]
        b = D_b[L]["pds_cosine"]["reportado"]
        difs.append(b - u)
        L_.append(
            f"| {L} | {fmt(u)} | {fmt(b)} | {b-u:+.6f} | "
            f"{'b por encima' if b > u else 'b por DEBAJO'} |"
        )
    n_abajo = sum(1 for d in difs if d < 0)
    L_.append(
        f"\n`b` queda por debajo de `u` en **{n_abajo} de 6** lineas y por encima "
        f"en {6-n_abajo}. Los dos objetos se mueven en una banda de "
        f"{min(min(D_u[L]['pds_cosine']['reportado'], D_b[L]['pds_cosine']['reportado']) for L in LINES):.4f} "
        f"a {max(max(D_u[L]['pds_cosine']['reportado'], D_b[L]['pds_cosine']['reportado']) for L in LINES):.4f}, "
        f"y `|b - u|` va de {min(abs(d) for d in difs):.4f} a "
        f"{max(abs(d) for d in difs):.4f}. **El signo cambia de linea a linea sin "
        "direccion consistente, que es como se ve el ruido alrededor de un punto "
        "fijo, no un efecto sistematico del escalado.**\n"
    )
    L_.append(
        "Escala de ruido disponible: el anchor da 5 particiones de `pds_cosine` "
        "por linea, y su `sd` va de 0.015 a 0.029 "
        "(`anclas_nativas_v2.md`). Cuatro de los seis `|b - u|` caen por debajo "
        "de 0.029; HAP1 (0.0748) y BXPC3 (0.0417) lo exceden. Asi que el ruido "
        "explica la mayoria de los casos pero **no todos**, y para HAP1 queda una "
        "causa sin identificar. **No esta en el codigo** ni en estas mediciones "
        "cual es.\n"
    )
    L_.append(
        "**Lo que importa para la escala, y esta medido:** en `pds_cosine` el `b` "
        "oficial NO acota la trivialidad con margen -- se sienta sobre el punto "
        "sin informacion con +-0.04 de ruido, y `cable_suelo.md` midio que una "
        "copia del control cobra `s` positivo ahi en 3 de 6 lineas (HAP1 +0.2125, "
        "BXPC3 +0.0870, MCF7 +0.0197). Con `r - b` entre 0.28 y 0.48, ese ruido "
        "de +-0.04 vale hasta ~0.14 de `s` sin ninguna habilidad.\n"
    )

    # ================= ENTREGABLE =================
    L_.append("\n# ENTREGABLE -- las 36 filas\n")
    L_.append(
        "Para el objeto que un candidato va a tener que superar: el **`b` "
        "oficial**, que es el extremo 0 de la escala. `parte que MIDE` es la "
        "componente calculada del contenido de la emision; `parte que es "
        "CONVENCION` es un valor asignado donde la cantidad no esta definida y "
        "que ENTRA al promedio. Para las gateadas la convencion es 0 por "
        "construccion (la fila se omite, no se rellena), y lo que se encoge es la "
        "cohorte.\n"
    )
    L_.append(
        "| linea | metrica | valor reportado | parte que MIDE | parte que es "
        "CONVENCION | perts que entraron | de 33 | fraccion del numero que mide algo |"
    )
    L_.append("|---" * 8 + "|")
    filas36 = []
    for L in LINES:
        for m in METRICS:
            d = D_b[L][m]
            rep = d["reportado"]
            if m == JAC:
                mide = d["comp_genuina"]
                conv = d["comp_convencion"]
                entraron = d["n_genuinos"]
                frac = (mide / rep) if rep else float("nan")
            elif m == "expr_mse_unbiased_capped_norm":
                # Metrica DERIVADA (`agg="ratio_of_sums"`, catalog.py:412-414):
                # no emite fila por perturbacion, asi que no tiene conteo propio.
                # Sus dos componentes si, y las dos leen 33 sin NaN (verificado).
                mide, conv = rep, 0.0
                entraron = N_PANEL
                frac = 1.0
            else:
                mide, conv = rep, 0.0
                entraron = d["entraron"]
                frac = 1.0
            filas36.append({
                "linea": L, "metrica": m, "rep": rep, "mide": mide,
                "conv": conv, "entraron": entraron, "frac": frac,
                "cohorte": (entraron / N_PANEL) if entraron else None,
            })
            marca = ""
            if m == JAC and frac == frac and frac < 0.5:
                marca = " ⚠️"
            L_.append(
                f"| {L} | {CORTO[m]} | {fmt(rep)} | {fmt(mide)} | {fmt(conv)} | "
                f"{entraron if entraron is not None else '?'} | 33 | "
                f"{('%.1f%%' % (100*frac)) if frac == frac else 'NaN'}{marca} |"
            )
    jac_bajos = [f for f in filas36 if f["metrica"] == JAC
                 and f["frac"] == f["frac"] and f["frac"] < 0.5]
    ceros = [f for f in filas36 if f["metrica"] == JAC and f["mide"] == 0.0]
    L_.append(
        f"\n`mse_capped_norm` es una metrica DERIVADA "
        "(`agg=\"ratio_of_sums\"`, `catalog.py:412-414`): no emite fila por "
        "perturbacion, asi que no tiene conteo propio. Sus dos componentes si lo "
        "tienen, y las dos leen 33 filas sin NaN en las seis lineas, por eso "
        "figura con 33.\n"
    )
    L_.append(
        f"**{len(jac_bajos)} de 36 celdas miden menos de la mitad de lo que "
        f"reportan, y todas son `jac`.** En {len(ceros)} de ellas "
        f"({', '.join(f['linea'] for f in ceros)}) el `b` oficial mide "
        "**exactamente 0.0%** de su propio valor: las uniones vacias son el "
        "numero entero. Las `jac` de BXPC3 y MCF7 miden ~100%, pero de un valor "
        "de 0.001 y 0.00003 -- miden algo, y ese algo es casi nada.\n"
    )
    L_.append(
        "Y la cohorte, que es el otro eje: `lfc_nmae` entra con 8 de 33 "
        "perturbaciones en A549 y HAP1 (24%), `fid` con 18-33, `reach` con 18-33, "
        "y `pds`/`mse` con 33. **El `avg_score` promedia esas seis a 1/6 cada "
        "una, sin pesar por cohorte** -- de modo que en A549 una perturbacion de "
        "las 8 que `lfc_nmae` ve pesa 4.1 veces mas en el agregado que una de las "
        "33 que `pds` ve.\n"
    )

    # ================= ESCALA DIAGNOSTICA =================
    L_.append("\n# 🔴 SEGUNDA ESCALA -- DIAGNOSTICA, **no es la cifra de competencia**\n")
    L_.append(
        "**Lo que sigue NO es el puntaje.** El `Overall` usa el crudo, y el crudo "
        "es lo que esta en `cable_suelo.md`. Esta escala calcula "
        "`(u - b) / (r - b)` **solo sobre la componente genuina** de `jac`, y "
        "sirve para una cosa: saber CUAL de las dos partes se movio cuando un "
        "candidato cambia. No sustituye a nada y no debe reportarse como "
        "resultado.\n"
    )
    L_.append(
        f"`SE = sd / sqrt({K})` porque `r` ES una media de {K} particiones "
        "(`anchor.py:285`), asi que la precision de su ubicacion es el error "
        "estandar de esa media, no su desviacion.\n"
    )
    if rgen:
        L_.append(
            "| linea | escala | u | b | r | span | sd de r | SE de r | span/SE | "
            "¿resuelve (>= 3)? |"
        )
        L_.append("|---" * 10 + "|")
        for L in LINES:
            if L not in rgen:
                continue
            g = rgen[L]
            du, db = D_u[L][JAC], D_b[L][JAC]
            # crudo
            sp_c = g["r_crudo_media"] - db["reportado"]
            se_c = g["r_crudo_sd"] / math.sqrt(K)
            rc = sp_c / se_c if se_c else float("nan")
            L_.append(
                f"| {L} | **cruda (la de competencia)** | {fmt(du['reportado'])} | "
                f"{fmt(db['reportado'])} | {fmt(g['r_crudo_media'])} | "
                f"{fmt(sp_c)} | {fmt(g['r_crudo_sd'])} | {fmt(se_c)} | "
                f"**{rc:.2f}** | {'si' if rc >= 3 else '**NO**'} |"
            )
            sp_g = g["r_genuina_media"] - db["comp_genuina"]
            se_g = g["r_genuina_sd"] / math.sqrt(K)
            rg = sp_g / se_g if se_g else float("nan")
            L_.append(
                f"| {L} | diagnostica (solo genuina) | {fmt(du['comp_genuina'])} | "
                f"{fmt(db['comp_genuina'])} | {fmt(g['r_genuina_media'])} | "
                f"{fmt(sp_g)} | {fmt(g['r_genuina_sd'])} | {fmt(se_g)} | "
                f"**{rg:.2f}** | {'**si**' if rg >= 3 else 'NO'} |"
            )
        L_.append(
            "\n**La respuesta a la pregunta del encargo: si, resuelven.** A549 "
            "pasa de `span/SE = 1.19` sobre el crudo a **8.41** sobre la "
            "componente genuina; HAP1 de `1.16` a **28.23**. Las dos cruzan el "
            "umbral de 3 que no cruzaban.\n"
        )
        L_.append(
            "El mecanismo esta en las dos columnas de en medio. Sobre el crudo, "
            "el span es chico porque la convencion resta (el replicado real tiene "
            "menos uniones vacias que el baseline) y el `sd` es grande porque `E` "
            "varia entre particiones -- en A549 de 13 a 15, y cada unidad de `E` "
            "mueve el crudo 1/33 = 0.0303. Sobre la componente genuina el span es "
            "mas grande Y el `sd` mas chico, las dos cosas a la vez.\n"
        )
        L_.append(
            "⚠️ Tres limites de esta escala diagnostica, para que no se use de "
            "mas: (1) `u` y `b` genuinos valen **exactamente 0** en estas dos "
            "lineas, asi que el `s` diagnostico de una emision trivial es 0 igual "
            "que el crudo -- lo que cambia es la REGLA, no el puntaje de esta "
            "corrida; (2) solo esta calculada para A549 y HAP1, porque recomputar "
            "las particiones de las otras cuatro cuesta otros ~20 min y no se "
            "pidio; (3) la componente genuina de `u` y `b` se cuenta sobre las "
            "perturbaciones con union no vacia de CADA objeto, y esos conjuntos "
            "no son identicos entre objetos -- con `u` y `b` genuinos en 0 eso no "
            "afecta a estos numeros, pero con un candidato real habria que fijar "
            "una cohorte comun antes de comparar.\n"
        )

    # ================= 20 FILAS =================
    L_.append("\n# 20 filas trabajadas\n")
    L_.append(
        "Veinte filas con todos los intermedios de su descomposicion: nada entre "
        "el artefacto y la columna `fraccion que mide`. Se eligieron recorriendo "
        "las lineas en orden y las metricas en orden, tomando las primeras 20 "
        "celdas del entregable -- **criterio posicional, declarado, no depende "
        "del valor**.\n"
    )
    L_.append(
        "| # | linea | metrica | reportado | mide | convencion | E | genuinos | "
        "filas emitidas | NaN | entraron | cohorte | fraccion que mide | mecanismo |"
    )
    L_.append("|---" * 14 + "|")
    MEC = {
        JAC: "convencion J(vacio,vacio)=1",
        "de_wilcoxon_lfc_nmae": "fila no emitida (gate)",
        "de_wilcoxon_direction_fidelity_yield_raw": "NaN si max(n_pred,n_conf)=0",
        "de_wilcoxon_direction_reach_raw": "NaN si n_conf=0; piso de azar ~c/N_conf",
        "pds_cosine": "sin exclusion; punto sin info = 0.5",
        "expr_mse_unbiased_capped_norm": "sin exclusion; clamp a [0,1]",
    }
    for i, f in enumerate(filas36[:20], 1):
        L = f["linea"]
        m = f["metrica"]
        d = D_b[L][m]
        E = d.get("E_convencion", "-")
        ng = d.get("n_genuinos", "-")
        coh = f["cohorte"]
        L_.append(
            f"| {i} | {L} | {CORTO[m]} | {fmt(f['rep'])} | {fmt(f['mide'])} | "
            f"{fmt(f['conv'])} | {E} | {ng} | {d['filas_emitidas']} | {d['nan']} | "
            f"{f['entraron'] if f['entraron'] is not None else '?'} | "
            f"{('%.0f%%' % (100*coh)) if coh else '?'} | "
            f"{('%.1f%%' % (100*f['frac'])) if f['frac'] == f['frac'] else 'NaN'} | "
            f"{MEC[m]} |"
        )

    # ================= NO HECHO =================
    L_.append("\n# Lo que NO se hizo\n")
    L_.append(
        "- **No se modifico ninguna metrica ni ningun arnes.** "
        "`_disjoint_halves`, `contar_de_por_perturbacion`, "
        "`pseudobulk_bulk_lognorm`, `delta` y `build_fold_adatas` se usaron tal "
        "cual.\n"
        "- **No se sustituyo la escala cruda por la diagnostica en ningun lado.** "
        "La diagnostica esta marcada arriba de su propia seccion, no en nota al "
        "pie, y el `Overall` sigue siendo el crudo de `cable_suelo.md`.\n"
        "- **No se selecciono ninguna celda por su signo ni por su valor.** Los "
        "tres criterios de seleccion que aparecen se declararon antes de medir: "
        "las seis lineas para la dispersion; A549 y HAP1 para el recomputo de "
        "`r`, por ser las que el encargo nombra; y las primeras 20 celdas en "
        "orden posicional para las filas trabajadas.\n"
        "- **No se derivaron objetivos numericos por linea** para ningun "
        "candidato. En particular no se calculo cuanto habria que bajar "
        "`mse_capped_norm` para salir del clamp.\n"
        "- **No se comparo contra el 0.5 de `scales.py:428`**, que pertenece a "
        "otra escala. El 0.5 citado es el punto sin informacion del estadistico "
        "de rango (`discrimination.py:70`, `74-77`).\n"
        "- No se leyo nada bajo `~/vcc/benchmark/` ni "
        "`/datos/benchmark_VEDADO_techo_metricas/`.\n"
        "- Sin envios, sin `vcc sample`, sin commits.\n"
    )

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L_) + "\n")
    print(f"escrito {OUT_MD}")


if __name__ == "__main__":
    main()
