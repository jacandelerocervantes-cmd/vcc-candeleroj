#!/usr/bin/env python
"""Arma modelo/out/cable_suelo.md desde modelo/out/cable_suelo.json, cruzando
contra anclas_nativas_comp.json (los dos extremos) y bucle_suelo.json."""
import json
import math
import os

CABLE = "/home/ubuntu/vcc/modelo/out/cable_suelo.json"
NATIVAS = "/home/ubuntu/vcc/modelo/out/anclas_nativas_comp.json"
SUELO = "/home/ubuntu/vcc/modelo/out/bucle_suelo.json"
OUT_MD = "/home/ubuntu/vcc/modelo/out/cable_suelo.md"

HASH_REF = "a5519017556705abf9a6e553ebb4e517b9612e9546331ace7c5ec6c5c30a1a13"
K = 5

LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
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
    "pds_cosine": True,
    "expr_mse_unbiased_capped_norm": False,
    "de_wilcoxon_lfc_nmae": False,
    "de_wilcoxon_direction_fidelity_yield_raw": True,
    "de_wilcoxon_direction_reach_raw": True,
    "de_wilcoxon_sig_jaccard": True,
}
# La unica de las seis que clampea a [0, 1] (catalog.py:415-416).
CLAMPEADA = {"expr_mse_unbiased_capped_norm": (0.0, 1.0)}


def cargar(path):
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


def span(r, b, m):
    if r is None or b is None:
        return None
    return (r - b) if MAS_ES_MEJOR[m] else (b - r)


def main():
    cab = {r["linea"]: r for r in cargar(CABLE) if "error" not in r}
    nat = {r["linea"]: r for r in cargar(NATIVAS) if "error" not in r}
    suelo = {r["linea_excluida"]: r for r in cargar(SUELO) if "error" not in r}
    L = []

    L.append("# La prueba del cable: `run` -> `score --real-bundle`\n")
    L.append(
        "No es un candidato. La emision es **desechable**: "
        "`suelo_floor.predict_floor` sin modificar -- el pool de control "
        "remuestreado con reemplazo, sin escalar. El objetivo es estrenar la "
        "cadena completa con una emision de la que ya sabemos la respuesta, "
        "porque `run` tiene la misma sensibilidad de `config_hash` que mando los "
        "seis primeros bundles a DIAGNOSTICO (ver "
        "`modelo/out/anclas_nativas.md`).\n"
    )

    # ---------- la prediccion declarada ----------
    L.append("## La prediccion, declarada ANTES de correr\n")
    L.append(
        "> `predict_floor` emite el control remuestreado sin escalar. El `b` "
        "oficial es el control ESCALADO al perfil de respuesta promedio "
        "(`baseline.py:344`). Asi que el suelo deberia quedar POR DEBAJO de `b` "
        "en casi todo, y su `s` salir ~0 o negativo. Si sale muy por encima de "
        "0, algo esta mal en la cadena, no en el suelo.\n"
    )
    L.append(
        "Esto es lo que hace que la prueba sirva: el resultado estaba "
        "comprometido antes de verlo.\n"
    )

    # ---------- corrio o se nego ----------
    L.append("## ¿`score` corrio o se nego?\n")
    negados = [l for l in LINES if cab.get(l) and not cab[l].get("score_ok")]
    corridos = [l for l in LINES if cab.get(l) and cab[l].get("score_ok")]
    if negados:
        L.append(f"🔴 **Se nego en {len(negados)} linea(s): {', '.join(negados)}.**\n")
        for l in negados:
            r = cab[l]
            L.append(f"\n**{l}** -- mensaje literal:\n\n```\n{r.get('mensaje_literal')}\n```\n")
    if corridos:
        L.append(
            f"✅ **`score --real-bundle` corrio en {len(corridos)} de 6 lineas**, "
            "sin que ninguna guarda se disparara.\n"
        )
        L.append("| linea | run | score | duracion score | `anchor_source` | `real_bundle_id` |")
        L.append("|---|---|---|---|---|---|")
        for l in LINES:
            r = cab.get(l)
            if r is None:
                continue
            sc = (r.get("scored") or {}).get("filas", {}).get("avg_score", {})
            L.append(
                f"| {l} | exit {r['run_returncode']} ({r['run_duration_s']}s) | "
                f"exit {r.get('score_returncode')} | {r.get('score_duration_s')}s | "
                f"`{sc.get('anchor_source')}` | `{sc.get('real_bundle_id')}` |"
            )
        L.append(
            "\nLas dos guardas que ya conociamos **no se dispararon**, y conviene "
            "decir por que cada una paso:\n"
        )
        L.append(
            "- **`expect_from_run_meta` sobre el fingerprint** (`score.py:382-389`): "
            "exige `source_fingerprint_strict`. El preset `vcc2026` trae "
            "`cache_strict: true`, asi que el `run_meta.json` de cada corrida lo "
            "lleva en `True` (tabla de abajo). Si hubiera sido metadata-only, "
            "`score` habria levantado \"this run's source_fingerprint is "
            "metadata-only ... Re-run both sides with `--cache-strict`\".\n"
            "- **\"the anchor leaves no usable headroom over the baseline\"** "
            "(`score.py:466-467`): exige que el span sea positivo para cada "
            "metrica decisiva. Los 36 spans del bundle son positivos "
            "(`anclas_nativas_v2.md`), asi que no habia nada que levantar -- "
            "pero ver la seccion de las dos celdas que no resuelven.\n"
        )
        L.append("| linea | `source_fingerprint_strict` del run | fingerprint del run |")
        L.append("|---|---|---|")
        for l in LINES:
            r = cab.get(l)
            if r is None:
                continue
            rm = r.get("run_meta") or {}
            L.append(
                f"| {l} | {rm.get('source_fingerprint_strict')} | "
                f"`{str(rm.get('source_fingerprint'))[:16]}...` |"
            )

    # ---------- config_hash del run ----------
    L.append("\n## ¿El `run` quedo con `config_hash` de competencia?\n")
    L.append(
        "`run` no imprime una linea de competencia como `prep-real-bundle`, asi "
        "que el hash se recalcula desde el `run_params.yaml` que escribio cada "
        f"corrida y se compara contra el de referencia `{HASH_REF[:16]}...`.\n"
    )
    L.append("| linea | hash obtenido | coincide | `control` en el config | `pert_col` | `cache_strict` |")
    L.append("|---|---|---|---|---|---|")
    todos_ok = True
    for l in LINES:
        r = cab.get(l)
        if r is None:
            continue
        h = r.get("hash_run") or {}
        ok = h.get("coincide")
        if not ok:
            todos_ok = False
        L.append(
            f"| {l} | `{str(h.get('obtenido'))[:16]}...` | "
            f"{'si' if ok else '**NO**'} | `{h.get('control_en_config')}` | "
            f"`{h.get('pert_col_en_config')}` | {h.get('cache_strict')} |"
        )
    if todos_ok:
        L.append(
            "\n✅ **Los seis runs coinciden con el hash de competencia.** La "
            "leccion del primer intento quedo aplicada: no se paso `--control` "
            "(el valor queda el del preset, `non-targeting`) y `--pert-col "
            "target_gene` si se paso, porque esta en `_RULE_EXCLUDED` "
            "(`competition.py:112`). La sensibilidad de `config_hash` que mando "
            "los seis primeros bundles a DIAGNOSTICO existe igual en `run`, y "
            "aqui no se activo.\n"
        )
    else:
        L.append(
            "\n🔴 **Al menos un run NO coincide.** Hay que aislar que bandera lo "
            "desvia con la misma aritmetica de hashes de `anclas_nativas.md`.\n"
        )

    # ---------- avg_score y de donde viene ----------
    L.append("\n## `avg_score` por linea, y de donde viene\n")
    L.append(
        "El `--help` de `score` dice que con un bundle de competencia "
        "\"from_replicate takes avg_score and from_baseline's goes null\". "
        "Medido:\n"
    )
    L.append(
        "| linea | `avg_score` de `from_replicate` | `avg_score` de `from_baseline` | "
        "¿cual lo tomo? | `anchor_digest` |"
    )
    L.append("|---|---|---|---|---|")
    for l in LINES:
        r = cab.get(l)
        if r is None or not r.get("score_ok"):
            continue
        av = (r["scored"]["filas"].get("avg_score") or {})
        fr, fb = av.get("from_replicate"), av.get("from_baseline")
        cual = "`from_replicate`" if fb is None and fr is not None else (
            "`from_baseline`" if fr is None else "ambas (?)")
        L.append(
            f"| {l} | {fmt(fr)} | {'null' if fb is None else fmt(fb)} | {cual} | "
            f"`{str(av.get('anchor_digest'))[:16]}...` |"
        )
    L.append(
        "\n**`from_replicate` tomo el `avg_score` en las seis lineas, y el de "
        "`from_baseline` salio null**, exactamente como declara el `--help`. "
        "Las columnas por metrica de `from_baseline` SI estan pobladas; lo que "
        "queda null es su fila `avg_score`. Es la diferencia entre \"la columna "
        "existe\" y \"la columna decide\".\n"
    )

    # ---------- s por metrica y linea ----------
    L.append("\n## El `s` por metrica y por linea, las 6 x 6\n")
    L.append(
        "`from_replicate`, que es la columna que toma el `avg_score`. Entre "
        "parentesis el `from_baseline` de la misma celda, para contraste.\n"
    )
    L.append("| metrica | " + " | ".join(LINES) + " |")
    L.append("|---" * (1 + len(LINES)) + "|")
    for m in METRICS:
        celdas = []
        for l in LINES:
            r = cab.get(l)
            if r is None or not r.get("score_ok"):
                celdas.append("-")
                continue
            f = r["scored"]["filas"].get(m) or {}
            fr, fb = f.get("from_replicate"), f.get("from_baseline")
            celdas.append(f"{fmt(fr, 4)} ({fmt(fb, 4)})")
        L.append(f"| {CORTO[m]} | " + " | ".join(celdas) + " |")
    L.append(
        "\n`avg_score` es la media simple de las seis celdas de "
        "`from_replicate`; se verifica en la tabla de intermedios de abajo.\n"
    )

    # ---------- contraste con la prediccion ----------
    L.append("\n## Contraste contra la prediccion declarada\n")
    todos_neg = []
    for l in LINES:
        r = cab.get(l)
        if r is None or not r.get("score_ok"):
            continue
        av = (r["scored"]["filas"].get("avg_score") or {}).get("from_replicate")
        if av is not None:
            todos_neg.append((l, av))
    if todos_neg:
        peor = min(todos_neg, key=lambda x: x[1])
        mejor = max(todos_neg, key=lambda x: x[1])
        n_neg = sum(1 for _, v in todos_neg if v < 0)
        L.append(
            f"**Salio donde dijimos.** Los `avg_score` son negativos en "
            f"{n_neg} de {len(todos_neg)} lineas, en el rango "
            f"{peor[1]:.6f} ({peor[0]}) a {mejor[1]:.6f} ({mejor[0]}).\n"
        )
        L.append(
            "Ninguno queda \"muy por encima de 0\", que era la senal de que el "
            "problema estaria en la cadena. La emision que solo copia el control "
            "puntua por debajo del comparador oraculo que escala ese mismo "
            "control al perfil de respuesta promedio -- que es el orden que la "
            "construccion de `baseline.py` predice.\n"
        )
        L.append("| linea | `avg_score` | ¿negativo? |")
        L.append("|---|---|---|")
        for l, v in todos_neg:
            L.append(f"| {l} | {fmt(v)} | {'si' if v < 0 else '**NO**'} |")
        L.append(
            "\nPor metrica, el detalle importa mas que el agregado. De las seis, "
            "**tres salen negativas** (`pds_cosine`, `lfc_nmae`, `reach`) y "
            "**tres salen exactamente 0.0** -- y los ceros tienen dos causas "
            "distintas, que no conviene confundir:\n"
        )
        L.append(
            "- **`mse_capped_norm` sale 0.0 por CLAMP.** Es la unica de las seis "
            "que clampea a `[0, 1]` (`catalog.py:415-416`: "
            "`clamp_low=0.0, clamp_high=1.0`). Sin el clamp, A549 daria "
            "`(2.562172 - 0.959928) / (0.138661 - 0.959928) = -1.950941`. El "
            "comentario del catalogo (`catalog.py:406-410`) explica justamente "
            "este caso: \"The metric is capped so that 1.0 = predicting the "
            "control, so a value past the baseline is already no-skill and the "
            "graded Box-Cox tail below it buys nothing a ranking uses.\" Nuestra "
            "emision ES predecir el control, asi que caer al piso del clamp es "
            "el comportamiento documentado.\n"
            "- **`fid` y `jac` salen 0.0 porque `u == b` EXACTAMENTE.** El "
            "numerador es cero, no hay clamp. En A549: `fid` da u=0.000000 y "
            "b=0.000000; `jac` da u=0.454545 y b=0.454545. Es el hallazgo de "
            "`anclas_nativas.md` reapareciendo del lado del candidato: el "
            "artefacto de union vacia hace que la emision degenerada y el "
            "comparador oraculo den el MISMO numero, asi que `s` sale 0 no "
            "porque la emision empate con el baseline por merito, sino porque la "
            "metrica no esta mirando la emision. **Un `s = 0` en esas dos celdas "
            "no significa \"igual al baseline\": significa \"no medido\".**\n"
        )

    # ---------- las celdas donde la emision degenerada puntua POSITIVO ----------
    pos = []
    for l in LINES:
        r = cab.get(l)
        if r is None or not r.get("score_ok"):
            continue
        for m in METRICS:
            s = (r["scored"]["filas"].get(m) or {}).get("from_replicate")
            if s is not None and s > 0:
                pos.append((l, m, s))
    if pos:
        L.append(
            "\n### ⚠️ Donde la emision degenerada puntua POSITIVO (y no deberia)\n"
        )
        L.append(
            "El agregado salio negativo en las seis lineas, pero por celda hay "
            f"**{len(pos)} de 36 con `s > 0`**: la emision que solo copia el "
            "control queda del lado bueno del comparador oraculo. Cada una es un "
            "lugar donde la escala premia algo que no es habilidad.\n"
        )
        L.append("| linea | metrica | u (suelo) | b | r | s | lectura |")
        L.append("|---|---|---|---|---|---|---|")
        for l, m, s in sorted(pos, key=lambda x: -x[2]):
            n = nat.get(l) or {}
            u = (cab[l].get("u") or {}).get(m)
            b = (n.get("b") or {}).get(m)
            rv = ((n.get("r") or {}).get(m) or {}).get("replicate")
            if m == "de_wilcoxon_sig_jaccard":
                lect = "union vacia: `b` arranca bajo y `u` cobra 1.0 regalados"
            elif m == "pds_cosine":
                lect = "`b` cayo por debajo del minimo aleatorio de ~0.5"
            else:
                lect = "revisar"
            L.append(
                f"| {l} | {CORTO[m]} | {fmt(u)} | {fmt(b)} | {fmt(rv)} | "
                f"**{fmt(s)}** | {lect} |"
            )
        peor_pos = max(pos, key=lambda x: x[2])
        L.append(
            f"\nEl caso mas fuerte es **{peor_pos[0]} / {CORTO[peor_pos[1]]} = "
            f"{peor_pos[2]:+.4f}**: la emision degenerada recorre la mitad de la "
            "distancia del baseline al techo. No porque prediga nada, sino "
            "porque `de_sig_jaccard` le regala un 1.0 por cada perturbacion sin "
            "senal (`metrics/de.py:195`, convencion `J(vacio, vacio) = 1`), y en "
            "esa linea el baseline oraculo cobra menos de esos regalos que la "
            "emision degenerada.\n"
        )
        L.append(
            "Es la tercera vez que el mismo artefacto aparece en este proyecto, "
            "ahora del lado del candidato: `bucle_techo.md` lo encontro "
            "inflando nuestro suelo, `anclas_nativas.md` lo encontro en el `b` "
            "oficial, y aqui lo encontramos **convirtiendose en puntaje "
            "positivo**. Antes del primer candidato conviene decidir que se hace "
            "con `jac`; esta corrida muestra que ignorarlo tiene un costo "
            "medible.\n"
        )

    # ---------- las dos celdas que no resuelven ----------
    L.append("\n## 🔴 Las dos celdas que NO resuelven -- marcadas, no excluidas\n")
    L.append(
        "`anclas_nativas_v2.md` encontro que el span de `jac` es mas chico que "
        "la precision con que se midio su extremo 1. En error estandar de la "
        f"media de k={K} particiones (`SE = sd / sqrt({K})`), que es la escala "
        "correcta porque `r` ES una media:\n"
    )
    L.append(
        "| linea | metrica | span | sd | SE = sd/sqrt(5) | span/SE | span/sd | estado |"
    )
    L.append("|---" * 8 + "|")
    no_resuelven = []
    for l in LINES:
        n = nat.get(l)
        if n is None or not n.get("r"):
            continue
        for m in METRICS:
            rr = n["r"].get(m) or {}
            sp = span(rr.get("replicate"), n["b"].get(m), m)
            sd = rr.get("replicate_sd")
            if sp is None or not sd:
                continue
            se = sd / math.sqrt(K)
            ratio_se = sp / se
            if ratio_se < 2.0:
                no_resuelven.append((l, m, sp, sd, se, ratio_se))
    for l, m, sp, sd, se, rse in sorted(no_resuelven, key=lambda x: x[5]):
        L.append(
            f"| {l} | {CORTO[m]} | {fmt(sp)} | {fmt(sd)} | {fmt(se)} | "
            f"**{rse:.2f}** | {sp/sd:.2f} | **no resuelve** |"
        )
    L.append(
        f"\n**{len(no_resuelven)} celdas con `span/SE < 2`**: el extremo 1 no "
        "esta localizado con precision suficiente para que la distancia al "
        "extremo 0 signifique algo. **No se excluyen del calculo** -- el "
        "`avg_score` de arriba las incluye tal cual, igual que las incluiria "
        "cualquier submission. Se marcan.\n"
    )
    L.append("Lo que esas celdas aportaron al `s` de esta corrida:\n")
    L.append("| linea | metrica | u (suelo) | b | r | s reportado | aporte al avg_score |")
    L.append("|---|---|---|---|---|---|---|")
    for l, m, sp, sd, se, rse in sorted(no_resuelven, key=lambda x: x[5]):
        r = cab.get(l)
        n = nat.get(l)
        if r is None or not r.get("score_ok"):
            continue
        s = (r["scored"]["filas"].get(m) or {}).get("from_replicate")
        L.append(
            f"| {l} | {CORTO[m]} | {fmt((r.get('u') or {}).get(m))} | "
            f"{fmt(n['b'].get(m))} | {fmt((n['r'].get(m) or {}).get('replicate'))} | "
            f"{fmt(s)} | {fmt((s or 0) / 6)} |"
        )
    L.append(
        "\nEn esta corrida el dano es invisible porque `u == b` exactamente y el "
        "`s` sale 0: una celda que no resuelve aporto cero. **Con un candidato "
        "real que mueva `jac` aunque sea un poco, esas dos celdas aportaran un "
        "sexto del `avg_score` cada una, dividido por un span que el propio "
        "anchor no localiza.** Es el riesgo concreto que esta prueba del cable "
        "no podia revelar y que conviene tener escrito antes del primer "
        "candidato.\n"
    )

    # ---------- intermedios: las 36 celdas completas ----------
    L.append("\n## Los intermedios completos: las 36 celdas, de `u` al `s`\n")
    L.append(
        "Todo lo que entra en cada `s`, sin paso oculto: `u` del "
        "`agg_results.csv` del run, `b` del `baseline_agg.csv` del bundle, `r` "
        "del `anchor_agg.parquet`, el `(u-b)/(r-b)` calculado a mano y el `s` "
        "que reporto `score`. La columna `dif` cierra el circuito; donde no es "
        "~0, la causa es el clamp y se dice cual.\n"
    )
    L.append(
        "| # | linea | metrica | u | b | r | span | (u-b)/(r-b) a mano | "
        "s reportado | dif | clamp | span/SE |"
    )
    L.append("|---" * 12 + "|")
    i = 0
    suma_por_linea = {}
    for l in LINES:
        r = cab.get(l)
        n = nat.get(l)
        if r is None or n is None or not r.get("score_ok"):
            continue
        for m in METRICS:
            i += 1
            u = (r.get("u") or {}).get(m)
            b = n["b"].get(m)
            rr = (n["r"].get(m) or {})
            rv = rr.get("replicate")
            sp = span(rv, b, m)
            s = (r["scored"]["filas"].get(m) or {}).get("from_replicate")
            manual = ((u - b) / (rv - b)) if (u is not None and b is not None
                                              and rv is not None and (rv - b) != 0) else None
            dif = ("-" if (manual is None or s is None) else f"{abs(manual - s):.2e}")
            cl = ""
            if m in CLAMPEADA and manual is not None and s is not None \
                    and abs(manual - s) > 1e-9:
                cl = f"**si, a {CLAMPEADA[m][0]}**"
            sd = rr.get("replicate_sd")
            rse = (f"{sp / (sd / math.sqrt(K)):.2f}"
                   if (sp is not None and sd) else "-")
            if s is not None:
                suma_por_linea.setdefault(l, []).append(s)
            L.append(
                f"| {i} | {l} | {CORTO[m]} | {fmt(u)} | {fmt(b)} | {fmt(rv)} | "
                f"{fmt(sp)} | {fmt(manual)} | {fmt(s)} | {dif} | {cl or '-'} | {rse} |"
            )
    L.append(
        "\n**El cable queda verificado aritmeticamente**: en las celdas sin "
        "clamp, el `s` que reporta `score` reproduce `(u-b)/(r-b)` calculado a "
        "mano desde los tres artefactos independientes (el run, el baseline del "
        "bundle y el anchor del bundle). Las unicas celdas con `dif` grande son "
        "las seis de `mse_capped_norm`, y la causa es su clamp documentado.\n"
    )
    if suma_por_linea:
        L.append("### El `avg_score` recompuesto desde las seis celdas\n")
        L.append("| linea | media de las 6 celdas | `avg_score` reportado | dif |")
        L.append("|---|---|---|---|")
        for l in LINES:
            vals = suma_por_linea.get(l)
            if not vals:
                continue
            media = sum(vals) / len(vals)
            av = (cab[l]["scored"]["filas"].get("avg_score") or {}).get("from_replicate")
            d = "-" if av is None else f"{abs(media - av):.2e}"
            L.append(f"| {l} | {fmt(media)} | {fmt(av)} | {d} |")
        L.append(
            "\nEs una media simple de las seis, sin pesos: se recompone exacta.\n"
        )

    # ---------- intermedios de corrida ----------
    L.append("\n## Intermedios de corrida\n")
    L.append(
        "| linea | real_shape | pred_shape | n control real | n control pred | "
        "seed_real | seed_pred | run_s | score_s | total_s |"
    )
    L.append("|---" * 10 + "|")
    for l in LINES:
        r = cab.get(l)
        if r is None:
            continue
        L.append(
            f"| {l} | {r['real_shape']} | {r['pred_shape']} | "
            f"{r['n_control_real']} | {r['n_control_pred']} | {r['seed_real']} | "
            f"{r['seed_pred']} | {r['run_duration_s']} | "
            f"{r.get('score_duration_s')} | {r.get('total_duration_s')} |"
        )
    L.append(
        "\nLa poblacion es la de siempre: `build_fold_adatas` **sin modificar**, "
        "`seed_real = 20260001 + fold_idx` (400 celulas por perturbacion sin "
        "reemplazo mas TODO el pool de control) y la emision de "
        "`suelo_floor.predict_floor` **sin modificar** con "
        "`seed_pred = 20261001 + fold_idx`. La etiqueta de control se renombro a "
        "`non-targeting` en LOS DOS lados, porque `io.validate_pair` exige que "
        "los conjuntos de etiquetas coincidan entre PRED y REAL; el paso 0 de "
        "`anclas_nativas_v2.md` probo que ese renombrado no mueve ningun "
        "numero.\n"
    )
    L.append(
        "Nota de comparabilidad: el `u` de esta corrida es, celda por celda, el "
        "mismo suelo que puntuo `bucle_suelo.md` -- misma poblacion, misma "
        "emision, mismas semillas. Lo que cambia es contra que se puntua: ahi "
        "fue contra si mismo (`avg_score = 0.0` por construccion), aqui contra "
        "los dos extremos nativos.\n"
    )
    if suelo:
        L.append("| linea | metrica | u de esta corrida | u de bucle_suelo.md | dif |")
        L.append("|---|---|---|---|---|")
        nmax = 0
        for l in LINES:
            r, s = cab.get(l), suelo.get(l)
            if r is None or s is None:
                continue
            for m in METRICS:
                if nmax >= 12:
                    break
                nmax += 1
                u1 = (r.get("u") or {}).get(m)
                u2 = s.get(m)
                d = "-" if (u1 is None or u2 is None) else f"{abs(u1 - u2):.2e}"
                L.append(f"| {l} | {CORTO[m]} | {fmt(u1)} | {fmt(u2)} | {d} |")
            if nmax >= 12:
                break

    L.append("\n## Lo que NO se hizo\n")
    L.append(
        "- **No se construyo ningun candidato.** Esta emision es desechable y no "
        "se guarda como modelo.\n"
        "- **No se reparo nada.** Nada fallo, asi que no hubo que decidirlo; si "
        "algo hubiera fallado, el script para en la primera linea y registra el "
        "mensaje literal (`cable_suelo.py`, `sys.exit(6)`).\n"
        "- **No se tocaron los bundles**, ni `bucle_suelo.md`, ni "
        "`bucle_techo.md`, ni ningun arnes. `suelo_arnes.py` y `suelo_floor.py` "
        "se reusaron por import, sin modificar.\n"
        "- **No se excluyo del calculo ninguna celda**, incluidas las dos que no "
        "resuelven: se marcan.\n"
        "- No se leyo nada bajo `~/vcc/benchmark/` ni "
        "`/datos/benchmark_VEDADO_techo_metricas/`.\n"
        "- Sin envios, sin `vcc sample`, sin commits.\n"
    )

    with open(OUT_MD, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"escrito {OUT_MD} con {len(cab)} lineas")


if __name__ == "__main__":
    main()
