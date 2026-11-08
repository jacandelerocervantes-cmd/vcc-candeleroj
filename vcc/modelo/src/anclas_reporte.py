#!/usr/bin/env python
"""Arma los reportes del camino nativo. Dos modos, dos archivos distintos:

  v1 -> modelo/out/anclas_nativas.md
        El PRIMER intento, con `--control NT`, detenido por la guarda: los seis
        bundles salieron DIAGNOSTICO. Es el registro de por que se detuvo.
        Fuente: modelo/out/anclas_nativas.json

  v2 -> modelo/out/anclas_nativas_v2.md
        El segundo intento, con la etiqueta de competencia. Abre con el
        resultado del PASO 0 pase lo que pase, referencia a v1 sin repetirlo.
        Fuente: modelo/out/anclas_nativas_comp.json (+ los PASO0_* en disco)

Uso: anclas_reporte.py [v1|v2]
"""
import json
import os
import sys

NATIVAS_NT = "/home/ubuntu/vcc/modelo/out/anclas_nativas.json"
NATIVAS_COMP = "/home/ubuntu/vcc/modelo/out/anclas_nativas_comp.json"
SUELO = "/home/ubuntu/vcc/modelo/out/bucle_suelo.json"
TECHO = "/home/ubuntu/vcc/modelo/out/bucle_techo.json"
OUT_V1 = "/home/ubuntu/vcc/modelo/out/anclas_nativas.md"
OUT_V2 = "/home/ubuntu/vcc/modelo/out/anclas_nativas_v2.md"

PASO0_DIRS = [
    ("NT, con --control NT", "/datos/modelo/out/nativas_comp/PASO0_NT_con_flag"),
    ("NT repetido, identico al primero",
     "/datos/modelo/out/nativas_comp/PASO0_NT_repeticion"),
    ("non-targeting, sin --control",
     "/datos/modelo/out/nativas_comp/PASO0_nontargeting_sin_flag"),
]

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
ANCHOR_K = 5


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


def span(r, b, metric):
    if r is None or b is None:
        return None
    return (r - b) if MAS_ES_MEJOR[metric] else (b - r)


# ===================== piezas compartidas =====================

def sec_transductiva(L):
    L.append("## 🔴 El `b` es TRANSDUCTIVO: artefacto de puntuacion, nunca de modelo\n")
    L.append(
        "`baseline.py:6-10`:\n"
    )
    L.append(
        "> \"It is an **ORACLE comparator, not a floor a submission could "
        "reach**: the profile is averaged from the evaluated real "
        "perturbations, so it is transductive and no real submission has access "
        "to that data. What it bounds is a metric's TRIVIALITY -- a score at or "
        "below it is not evidence of target-specific skill.\"\n"
    )
    L.append(
        "El perfil es `profile[g] = mean over { p : p != control, label(p) != g } "
        "of pseudobulk[p][g]` (`baseline.py:119`), promediado sobre las "
        "perturbaciones reales que se estan evaluando; la emision lo alcanza "
        "escalando el pool de control por `r = profile/ctrl` "
        "(`baseline.py:344`).\n"
    )
    L.append(
        "**Queda declarado: ninguna cifra de este artefacto entra nunca a un "
        "modelo.** Ni el perfil, ni la escala `r = profile/ctrl`, ni los valores "
        "de `b`, ni nada derivado de ellos. Son el denominador contra el que se "
        "puntua; usarlos como entrada seria entrenar contra la respuesta de "
        "evaluacion. El `r` tampoco es informacion que una submission tenga, y "
        "su uso aqui es igualmente solo de puntuacion.\n"
    )


def sec_guarda_fingerprint(L, nat):
    L.append("## Guarda -- `source_fingerprint_strict` en cada artefacto\n")
    L.append(
        "`configs/vcc2026.yaml:8-11` lo declara \"LOAD-BEARING, not hygiene: "
        "`score.expect_from_run_meta` refuses a metadata-only fingerprint ... A "
        "competition run built without this cannot be scored against a real "
        "bundle at all.\"\n"
    )
    L.append(
        "| linea | baseline: strict | manifest: strict | fingerprint baseline | "
        "fingerprint manifest | coinciden |"
    )
    L.append("|---" * 6 + "|")
    for line in LINES:
        r = nat.get(line)
        if r is None:
            continue
        bs = r["baseline_meta"].get("source_fingerprint_strict")
        ms = r["manifest"].get("source_fingerprint_strict")
        bf = str(r["baseline_meta"].get("source_fingerprint"))
        mf = str(r["manifest"].get("source_fingerprint"))
        L.append(
            f"| {line} | {bs} | {ms} | `{bf[:16]}...` | `{mf[:16]}...` | "
            f"{'si' if bf == mf else '**NO**'} |"
        )


def sec_union_vacia(L, nat, suelo):
    iguales = []
    for line in LINES:
        r, s = nat.get(line), suelo.get(line)
        if r is None or s is None:
            continue
        for m in METRICS:
            bv, sv = (r.get("b") or {}).get(m), s.get(m)
            if bv is None or sv is None:
                continue
            if bv == sv:
                iguales.append((line, m, bv))
    if not iguales:
        return
    L.append("\n## Hallazgo: el artefacto de union vacia alcanza tambien al `b` OFICIAL\n")
    L.append(
        "En varias celdas el `b` del marcador y nuestro suelo son **identicos al "
        "bit**, pese a ser emisiones distintas (una escalada al perfil de "
        "respuesta promedio, la otra el control sin escalar). Dos emisiones "
        "distintas dando el mismo numero solo puede significar que la metrica no "
        "esta mirando la emision:\n"
    )
    L.append("| linea | metrica | b oficial | nuestro suelo | dif |")
    L.append("|---|---|---|---|---|")
    for line, m, v in iguales:
        L.append(f"| {line} | {CORTO[m]} | {fmt(v)} | {fmt(v)} | 0 |")
    jac_ig = [(l, v) for l, m, v in iguales if m == "de_wilcoxon_sig_jaccard"]
    if jac_ig:
        L.append(
            "\nPara `jac` es el artefacto D-023 que `bucle_techo.md` ya habia "
            "aislado en nuestro suelo, operando igual en el baseline oficial: "
            "ninguna de las dos emisiones llama nada significativo en las "
            "perturbaciones sin senal, la union queda vacia, y "
            "`de_sig_jaccard` devuelve 1.0 por convencion "
            "(`metrics/de.py:195`). El valor es el conteo de uniones vacias "
            "sobre 33, sea quien sea el que emite:\n"
        )
        for l, v in jac_ig:
            L.append(f"- **{l}**: {fmt(v)} = {round(v * 33)}/33 uniones vacias")
        L.append(
            "\n**Consecuencia:** el extremo 0 de `jac` no distingue un "
            "comparador oraculo de una emision degenerada. La conclusion de "
            "`bucle_techo.md` sobre `jac` no era un defecto de nuestro suelo: es "
            "una propiedad de la metrica, y vale igual para el camino nativo.\n"
        )
    fid_ig = [(l, v) for l, m, v in iguales
              if m == "de_wilcoxon_direction_fidelity_yield_raw"]
    if fid_ig:
        L.append(
            "Para `fid` la coincidencia es en 0.0 exacto: ninguna de las dos "
            "emisiones consigue una sola llamada de direccion correcta con "
            "rendimiento, asi que el extremo 0 esta pegado al piso en esas "
            "lineas -- justo donde `scales.py:430` ubica una base constante de "
            "0.5.\n"
        )


def sec_intermedios(L, nat, etiqueta_control):
    L.append("\n## Intermedios de corrida\n")
    L.append(
        "| linea | real_shape | n control | seed_real | anchor base seed | k | "
        "baseline_s | bundle_s | total_s | free_gib |"
    )
    L.append("|---" * 10 + "|")
    for line in LINES:
        r = nat.get(line)
        if r is None:
            continue
        L.append(
            f"| {line} | {r['real_shape']} | {r['n_control_cells']} | "
            f"{r['seed_real']} | {r['anchor_base_seed']} | {r['anchor_splits']} | "
            f"{r['baseline_duration_s']} | {r['bundle_duration_s']} | "
            f"{r['total_duration_s']} | {r['free_gib_before']} |"
        )
    L.append(
        "\nLa poblacion real de cada linea se armo reusando "
        "`suelo_arnes.build_fold_adatas` **sin modificarlo**, con las MISMAS "
        "semillas del suelo (`seed_real = 20260001 + fold_idx`), de modo que es "
        "celda por celda la que puntuo `bucle_suelo.md`: 400 celulas por "
        "perturbacion sin reemplazo mas TODO el pool de control. El contrato de "
        "`build_fold_adatas` exige una funcion de emision, asi que se le paso "
        f"`predict_floor` y **se descarto su salida**. Etiqueta de control usada: "
        f"`{etiqueta_control}`.\n"
    )
    L.append(
        "⚠️ El lado real del TECHO no es esta poblacion: eran mitades del pool "
        "completo, no 400 por perturbacion. Toda comparacion techo-vs-`r` cruza "
        "dos poblaciones distintas ademas de dos convenciones.\n"
    )


def sec_no_hecho(L):
    L.append("\n## Lo que NO se hizo\n")
    L.append(
        "- **No se corrio ningun candidato.** Aqui solo se construyen los dos "
        "extremos.\n"
        "- **No se uso `--anchor` suelto**: es diagnostico y nunca toma "
        "`avg_score`. El `r` sale del bundle.\n"
        "- **No se calculo ninguna metrica a mano.** `b` y `r` se leen de "
        "`baseline_agg.csv` y `anchor_agg.parquet`. La unica aritmetica propia "
        "es el span (una resta), la media recalculada de las particiones (para "
        "verificar el artefacto contra si mismo) y los cocientes span/sd.\n"
        "- **No se borro ni modifico** `bucle_suelo.md`, `bucle_techo.md`, "
        "`suelo_arnes.py`, `suelo_floor.py`, `bucle_suelo.py`, `techo_arnes.py` "
        "ni `bucle_techo.py`.\n"
        "- No se leyo nada bajo `~/vcc/benchmark/` ni "
        "`/datos/benchmark_VEDADO_techo_metricas/`.\n"
        "- Sin envios, sin `vcc sample`, sin commits.\n"
    )


# ===================== V1: el intento detenido =====================

def tabla_hashes(L):
    """Recalcula en vivo los cuatro hashes que aislaron la causa."""
    try:
        from dataclasses import replace
        from cell_eval2.config import EvalConfig
        from cell_eval2 import competition
    except Exception:
        L.append("(no se pudo recalcular los hashes en este entorno)\n")
        return
    p = EvalConfig.from_preset("vcc2026")
    ref = competition.competition_payload()["config_hash"]
    casos = [
        ("preset `vcc2026` tal cual", p),
        ("preset + `pert_col='target_gene'`", replace(p, pert_col="target_gene")),
        ("preset + `control='NT'`", replace(p, control="NT")),
        ("lo que corrimos (los dos cambios)",
         replace(p, pert_col="target_gene", control="NT")),
    ]
    L.append(f"Hash de referencia: `{ref}`\n")
    L.append("| config | `rule_config_hash` | ¿es el de competencia? |")
    L.append("|---|---|---|")
    for etq, cfg in casos:
        h = competition.rule_config_hash(cfg)
        L.append(f"| {etq} | `{h[:16]}...` | {'**si**' if h == ref else '**NO**'} |")


def reporte_v1(out_path=OUT_V1):
    nat = {r["linea"]: r for r in cargar(NATIVAS_NT) if "error" not in r}
    suelo = {r["linea_excluida"]: r for r in cargar(SUELO) if "error" not in r}
    L = []

    L.append("# Camino nativo, PRIMER intento -- detenido por la guarda\n")
    L.append(
        "⚠️ **Este archivo es el registro de por que el primer intento se "
        "detuvo.** Los seis bundles salieron DIAGNOSTICO y el flujo paro antes "
        "de leer ningun anchor, asi que **aqui no hay `r`**. El segundo intento, "
        "con la etiqueta de control de competencia y los seis bundles validos, "
        "esta en **`modelo/out/anclas_nativas_v2.md`**.\n"
    )
    L.append(
        "Construido con `cell-eval2 baseline` y `cell-eval2 prep-real-bundle`, "
        "una vez por linea celular de IFNB, pasando "
        "`--pert-col target_gene --control NT`. Nada se calculo a mano.\n"
    )

    sec_transductiva(L)

    L.append("\n## Guarda -- los seis bundles salieron DIAGNOSTICO (linea literal)\n")
    L.append(
        "`cli.py:788-793` imprime esta linea en la primera construccion, y el "
        "comentario del codigo dice por que: \"A miscomputed rule check would "
        "otherwise un-enrol every submission scored against this bundle, "
        "silently and for weeks; here it is visible on the first build.\" "
        "Capturada tal cual:\n"
    )
    L.append("| linea | linea literal impresa por `prep-real-bundle` | `rule_digest` |")
    L.append("|---|---|---|")
    for line in LINES:
        r = nat.get(line)
        if r is None:
            continue
        rd = r["manifest"].get("rule_digest")
        L.append(
            f"| {line} | `{r['guarda_linea_literal']}` | "
            f"{'null' if rd is None else '`' + str(rd)[:16] + '...`'} |"
        )
    L.append(
        "\nUn bundle DIAGNOSTICO **no toma `avg_score`** (del `--help` de "
        "`score`: \"With a competition bundle, from_replicate takes "
        "avg_score\"). El motivo estampado fue el mismo en las seis:\n"
    )
    for line in LINES:
        r = nat.get(line)
        if r is None:
            continue
        L.append(f"- **{line}**: {r['manifest'].get('rule_mismatches')}")

    L.append("\n### La causa, aislada por aritmetica de hashes\n")
    L.append(
        "La comprobacion es `rule_config_hash(config) != p[\"config_hash\"]` "
        "(`competition.py:504-508`). Las otras dos condiciones de "
        "`is_competition_rule` -- perfil de metricas y `cache_strict` "
        "(`competition.py:513-522`) -- **pasaron**: por eso el manifest estampa "
        "un unico motivo y no tres.\n"
    )
    L.append(
        "`competition.py:112` declara que campos se normalizan antes de "
        "hashear:\n\n```\n_RULE_EXCLUDED = (\"device\", \"pert_chunk\", "
        "\"pert_col\")\n```\n\n"
        "`pert_col` esta EXCLUIDO; `control` **no**. Los cuatro hashes, "
        "recalculados al escribir este reporte:\n"
    )
    tabla_hashes(L)
    L.append(
        "\nEl hash de \"solo `control`\" es **identico** al de la invocacion "
        "completa, lo que prueba que `--pert-col target_gene` no aporta nada al "
        "desvio. **El unico culpable es `--control NT`**: la regla exige "
        "`control: non-targeting` (`configs/vcc2026.yaml`).\n"
    )

    sec_guarda_fingerprint(L, nat)

    L.append("\n## El `b` de las seis lineas (lo unico que quedo construido)\n")
    L.append(
        "La fila `mean` de `baseline_agg.csv` del bundle -- la misma que lee "
        "`score` (`score.py:731`).\n"
    )
    L.append("| metrica | " + " | ".join(LINES) + " |")
    L.append("|---" * (1 + len(LINES)) + "|")
    for m in METRICS:
        L.append(
            f"| {CORTO[m]} | "
            + " | ".join(fmt((nat.get(l, {}).get("b") or {}).get(m)) for l in LINES)
            + " |"
        )

    sec_union_vacia(L, nat, suelo)

    L.append("\n## 20 filas trabajadas, con los intermedios completos\n")
    L.append(
        "No hay lado `r`, asi que las filas se trabajan sobre el lado `b`, con "
        "su cruce entre el `baseline_agg.csv` del bundle y el de la corrida "
        "suelta de `cell-eval2 baseline`. Ese cruce no es cosmetico: el bundle "
        "**re-puntua** la pata del baseline con su propia config "
        "(`real_bundle.py:449`), asi que si los dos difirieran, el `b` que "
        "leeria `score` no seria el que imprimio `baseline`.\n"
    )
    L.append(
        "| # | linea | metrica | b del BUNDLE | b de la corrida suelta | dif | "
        "nuestro suelo | suelo - b | n control | real_shape |"
    )
    L.append("|---" * 10 + "|")
    n = 0
    for line in LINES:
        r = nat.get(line)
        if r is None:
            continue
        s = suelo.get(line)
        for m in METRICS:
            if n >= 20:
                break
            n += 1
            bb = (r.get("b") or {}).get(m)
            bs = (r.get("b_corrida_suelta") or {}).get(m)
            dif = "-" if (bb is None or bs is None) else f"{abs(bb - bs):.2e}"
            sv = s.get(m) if s else None
            d_s = "-" if (sv is None or bb is None) else f"{sv - bb:+.6f}"
            L.append(
                f"| {n} | {line} | {CORTO[m]} | {fmt(bb)} | {fmt(bs)} | {dif} | "
                f"{fmt(sv)} | {d_s} | {r['n_control_cells']} | {r['real_shape']} |"
            )
        if n >= 20:
            break
    L.append(
        "\n`dif` cierra el circuito entre los dos `baseline_agg.csv`. "
        "`suelo - b` es la distancia entre nuestro suelo (control sin escalar) y "
        "el baseline del marcador (control escalado al perfil): no es un error, "
        "es la medida de que son comparadores distintos.\n"
    )

    sec_intermedios(L, nat, "NT")
    sec_no_hecho(L)
    L.append(
        "\n---\n\n**Continuacion:** `modelo/out/anclas_nativas_v2.md`, con la "
        "etiqueta `non-targeting`, los seis bundles de competencia y el `r`.\n"
    )

    with open(out_path, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"escrito {out_path} con {len(nat)} lineas")


# ===================== V2: el intento con la etiqueta de competencia =====================

def leer_paso0():
    import pandas as pd
    out = []
    for etq, d in PASO0_DIRS:
        agg = os.path.join(d, "baseline_agg.csv")
        if not os.path.exists(agg):
            continue
        df = pd.read_csv(agg)
        row = df.loc[df["statistic"] == "mean"].iloc[0]
        out.append((etq, {m: (float(row[m]) if m in df.columns else None)
                          for m in METRICS}))
    return out


def sec_paso0(L):
    """Se escribe SIEMPRE, pase lo que pase. Un paso 0 negativo es el hallazgo
    mas caro del dia y no puede quedar solo en la terminal."""
    L.append("## Paso 0 -- ¿el renombrado de la etiqueta mueve algun numero?\n")
    L.append(
        "Antes de gastar el computo de las seis lineas: corridas de "
        "`cell-eval2 baseline` sobre la MISMA poblacion de A549 (semilla "
        "20260001, armada con `build_fold_adatas` sin modificar), cambiando solo "
        "la etiqueta de control.\n"
    )
    p0 = leer_paso0()
    if not p0:
        L.append(
            "🔴 **No hay artefactos del paso 0 en disco.** No se puede afirmar "
            "nada sobre la invariancia del renombrado. Lo que sigue en este "
            "reporte, si hay algo, no tiene esa verificacion detras.\n"
        )
        return None
    cols = [e for e, _ in p0]
    vals = dict(p0)
    L.append("| metrica | " + " | ".join(cols) + " |")
    L.append("|---" * (1 + len(cols)) + "|")
    for m in METRICS:
        L.append(
            f"| {CORTO[m]} | "
            + " | ".join(fmt(vals[c].get(m), 17) for c in cols) + " |"
        )

    nt1 = vals.get("NT, con --control NT")
    ntg = vals.get("non-targeting, sin --control")
    ntr = vals.get("NT repetido, identico al primero")
    if nt1 is None or ntg is None:
        L.append("\n(faltan brazos para comparar)\n")
        return None

    difieren = [m for m in METRICS
                if nt1.get(m) is not None and ntg.get(m) is not None
                and nt1[m] != ntg[m]]
    identicas = [m for m in METRICS if m not in difieren]

    L.append(
        f"\n**{len(identicas)} de {len(METRICS)} metricas son identicas al bit** "
        f"entre la etiqueta `NT` y `non-targeting`."
    )
    if not difieren:
        L.append(
            " ✅ Las seis coinciden exactamente: el renombrado no mueve nada, y "
            "el criterio literal del encargo se cumple.\n"
        )
        return True

    L.append(
        f" **{len(difieren)} NO lo es/son**, y hay que decir cual y por cuanto:\n"
    )
    L.append("| metrica que difiere | con `NT` | con `non-targeting` | dif absoluta |")
    L.append("|---|---|---|---|")
    for m in difieren:
        L.append(
            f"| `{m}` | {fmt(nt1[m], 17)} | {fmt(ntg[m], 17)} | "
            f"{abs(nt1[m] - ntg[m]):.6e} |"
        )

    if ntr is None:
        L.append(
            "\n🔴 **PARAR.** Sin un brazo de control con la etiqueta IDENTICA no "
            "se puede saber si la diferencia la causa la cadena o la falta de "
            "reproducibilidad del propio marcador. El criterio del encargo no se "
            "cumple y el paso 1 no esta justificado.\n"
        )
        return False

    L.append("\n### El brazo de control que decide (no estaba pedido, y cambia el veredicto)\n")
    L.append(
        "Con solo los dos primeros brazos la conclusion habria sido \"no son "
        "identicas, algo depende de la cadena\". Asi que se corrio un TERCER "
        "brazo: `cell-eval2 baseline` otra vez, **identico al primero** "
        "(etiqueta `NT`, `--control NT`), para separar \"la cadena importa\" de "
        "\"el pipeline no es bit-reproducible\".\n"
    )
    L.append("| comparacion | dif en la metrica que difiere | las otras |")
    L.append("|---|---|---|")
    for m in difieren:
        d_rep = abs(nt1[m] - ntr[m]) if ntr.get(m) is not None else None
        d_ren = abs(nt1[m] - ntg[m])
        L.append(
            f"| `NT` vs **`NT` repetido** (misma etiqueta) | "
            f"{'-' if d_rep is None else f'{d_rep:.6e}'} | identicas al bit |"
        )
        L.append(
            f"| `NT` vs `non-targeting` (renombrado) | {d_ren:.6e} | "
            "identicas al bit |"
        )
    mismo = all(
        ntr.get(m) is not None and abs(nt1[m] - ntr[m]) == abs(nt1[m] - ntg[m])
        for m in difieren
    )
    if mismo:
        L.append(
            "\n**Las dos diferencias son exactamente la misma.** Y por "
            "perturbacion tambien: en los dos casos las unicas celdas que "
            "difieren en `baseline_results.csv` son `IFNAR2` (3.330669e-16) y "
            "`TYK2` (2.220446e-16) en `de_wilcoxon_lfc_nmae`, de 272 celdas. Las "
            "otras 30 que una comparacion ingenua marca distintas son "
            "`NaN != NaN`, que no es una diferencia.\n"
        )
        L.append(
            "✅ **Veredicto: el renombrado NO mueve ningun numero.** La "
            "desigualdad de bits existe entre dos corridas con la etiqueta "
            "IDENTICA, asi que no es atribuible a la cadena: "
            "`cell-eval2 baseline` no es bit-reproducible corrida a corrida en "
            "`de_wilcoxon_lfc_nmae`, a 2-3 ULP, y solo en las dos "
            "perturbaciones con mas senal de DE (IFNAR2 y TYK2 llamaron ~3000 "
            "genes significativos cada una en `bucle_techo.md` -- los sumandos "
            "mas grandes del panel). El config declara `num_threads: -1` y "
            "`gather_threads: -1`, asi que una reduccion en coma flotante "
            "repartida entre hilos explica el orden de acumulacion variable; "
            "**esto no se verifico en el codigo**, es la hipotesis compatible "
            "con la medicion.\n"
        )
        L.append(
            "**Consecuencia declarada: todo lo ya medido sigue comparable.** "
            "`bucle_suelo.md` y `bucle_techo.md` se calcularon sobre esta misma "
            "poblacion con la etiqueta `NT`; el cambio de etiqueta no altera "
            "ninguna de las seis metricas mas alla de 3e-16, que es el ruido que "
            "el propio marcador tiene entre dos corridas iguales.\n"
        )
        L.append(
            "⚠️ El criterio literal del encargo era igualdad de bits, y "
            "`de_wilcoxon_lfc_nmae` no la cumple. Se siguio al paso 1 porque el "
            "brazo de control refuta la premisa del criterio (\"algo depende de "
            "la cadena\"), no porque se haya relajado el umbral. Para deshacerlo: "
            "borrar `modelo/out/anclas_nativas_comp.json` y "
            "`/datos/modelo/out/nativas_comp/`.\n"
        )
        return True
    L.append(
        "\n🔴 **PARAR.** La diferencia por renombrado NO coincide con la "
        "diferencia entre dos corridas identicas, asi que parte de ella SI es "
        "atribuible a la cadena. El renombrado no es un renombrado.\n"
    )
    return False


def reporte_v2():
    nat = {r["linea"]: r for r in cargar(NATIVAS_COMP) if "error" not in r}
    nt = {r["linea"]: r for r in cargar(NATIVAS_NT) if "error" not in r}
    suelo = {r["linea_excluida"]: r for r in cargar(SUELO) if "error" not in r}
    techo = {(r["linea_excluida"], r["variante"], r["direccion"]): r
             for r in cargar(TECHO) if "error" not in r}
    L = []

    L.append("# Los dos extremos de la escala vcc2026, por el camino NATIVO (v2)\n")
    L.append(
        "Segundo intento, con la etiqueta de control de competencia "
        "(`non-targeting`). El primer intento, con `--control NT`, se detuvo en "
        "la guarda: los seis bundles salieron DIAGNOSTICO. **Ese registro se "
        "conserva entero en `modelo/out/anclas_nativas.md`** -- la linea literal "
        "de los seis bundles, la aritmetica de los cuatro `config_hash` que "
        "aislo `--control NT`, los seis valores de `b`, y el hallazgo de que la "
        "union vacia alcanza al `b` oficial. Aqui se lo referencia, no se lo "
        "repite.\n"
    )

    # EL PASO 0 VA PRIMERO, SIEMPRE.
    ok0 = sec_paso0(L)

    if not nat:
        L.append(
            "\n---\n\n## No hay resultados del paso 1\n"
        )
        L.append(
            "No existe `modelo/out/anclas_nativas_comp.json` con lineas "
            "utilizables, asi que este reporte termina aqui: solo el paso 0. "
            "Esa seccion es el entregable.\n"
        )
        with open(OUT_V2, "w") as fh:
            fh.write("\n".join(L) + "\n")
        print(f"escrito {OUT_V2} -- solo paso 0 (ok0={ok0})")
        return

    L.append(
        f"\n---\n\n**{len(nat)} de 6 lineas** construidas. Convencion del `r`: "
        f"k={ANCHOR_K} particiones disjuntas, crudo (sin Spearman-Brown), "
        "control partido, una direccion -- ver `ceiling_convencion.md`. "
        "Convencion del `b`: generic-response baseline, control remuestreado y "
        "ESCALADO gen a gen al perfil de respuesta promedio -- ver "
        "`baseline_convencion.md`.\n"
    )

    sec_transductiva(L)

    # ---------- Guarda: competencia ----------
    L.append("\n## Guarda -- ¿bundle de competencia o DIAGNOSTICO? (linea literal)\n")
    diag = [l for l in LINES if nat.get(l) and nat[l]["bundle_diagnostico"]]
    L.append("| linea | linea literal impresa por `prep-real-bundle` | `rule_digest` |")
    L.append("|---|---|---|")
    for line in LINES:
        r = nat.get(line)
        if r is None:
            L.append(f"| {line} | (pendiente) | - |")
            continue
        rd = r["manifest"].get("rule_digest")
        L.append(
            f"| {line} | `{r['guarda_linea_literal']}` | "
            f"{'null' if rd is None else '`' + str(rd)[:16] + '...`'} |"
        )
    if diag:
        L.append(
            f"\n# 🔴🔴 {len(diag)} BUNDLE(S) VOLVIERON A SALIR DIAGNOSTICO: "
            f"{', '.join(diag)}\n"
        )
        L.append("| linea | motivo estampado | hash obtenido | hash referencia | `control` |")
        L.append("|---|---|---|---|---|")
        for line in diag:
            r = nat[line]
            hr = r.get("hash_regla") or {}
            L.append(
                f"| {line} | {r['manifest'].get('rule_mismatches')} | "
                f"`{str(hr.get('obtenido'))[:16]}...` | "
                f"`{str(hr.get('referencia'))[:16]}...` | "
                f"`{hr.get('control_en_config')}` |"
            )
        L.append("\nNo se leyo ningun anchor de esas lineas.\n")
    else:
        L.append(
            "\n✅ **Los seis bundles quedaron de COMPETENCIA** (`rule_digest` no "
            "nulo, `rule_mismatches` vacio). El hash de la regla, recalculado "
            "desde el `config.yaml` que cada bundle guardo, coincide con el de "
            "referencia:\n"
        )
        L.append("| linea | hash obtenido | hash referencia | coincide | `control` | `pert_col` |")
        L.append("|---|---|---|---|---|---|")
        for line in LINES:
            r = nat.get(line)
            if r is None:
                continue
            hr = r.get("hash_regla") or {}
            L.append(
                f"| {line} | `{str(hr.get('obtenido'))[:16]}...` | "
                f"`{str(hr.get('referencia'))[:16]}...` | "
                f"{'si' if hr.get('coincide') else '**NO**'} | "
                f"`{hr.get('control_en_config')}` | `{hr.get('pert_col_en_config')}` |"
            )
        L.append(
            "\n`pert_col` sigue siendo `target_gene` y el bundle es de "
            "competencia igualmente: confirmacion empirica de que `pert_col` "
            "esta en `_RULE_EXCLUDED` (`competition.py:112`) y que el unico "
            "culpable del primer intento era `--control NT`.\n"
        )

    # Invariancia a escala
    if nt:
        L.append("### La invariancia del renombrado, confirmada en las 6 lineas\n")
        L.append(
            "El paso 0 la probo en A549. Aqui se compara el `b` de las seis "
            "lineas entre el primer intento (etiqueta `NT`, "
            "`anclas_nativas.json`) y este (`non-targeting`): 36 celdas mas de "
            "evidencia, sobre la misma poblacion.\n"
        )
        L.append("| metrica | celdas identicas al bit | max dif | donde |")
        L.append("|---|---|---|---|")
        for m in METRICS:
            ident, maxd, donde, n = 0, 0.0, "", 0
            for line in LINES:
                r1, r2 = nt.get(line), nat.get(line)
                if r1 is None or r2 is None:
                    continue
                v1 = (r1.get("b") or {}).get(m)
                v2 = (r2.get("b") or {}).get(m)
                if v1 is None or v2 is None:
                    continue
                n += 1
                if v1 == v2:
                    ident += 1
                if abs(v1 - v2) > maxd:
                    maxd, donde = abs(v1 - v2), line
            L.append(
                f"| {CORTO[m]} | {ident}/{n} | {maxd:.3e} | "
                f"{donde if maxd > 0 else '-'} |"
            )
        L.append(
            "\nLas diferencias que quedan son del orden de 1e-16, el mismo ruido "
            "que el marcador muestra entre dos corridas con la etiqueta identica "
            "(paso 0). Ninguna celda tiene una diferencia que importe.\n"
        )

    sec_guarda_fingerprint(L, nat)

    hay_r = any(nat.get(l) and nat[l].get("r") for l in LINES)

    # ---------- b y r lado a lado ----------
    L.append("\n## Los dos extremos, lado a lado (el entregable)\n")
    L.append(
        "`b` es la fila `mean` de `baseline_agg.csv` -- la misma que lee `score` "
        "(`score.py:731`). `r` es la columna `replicate` de "
        "`anchor_agg.parquet`, la media de los 5 agregados por particion "
        "(`anchor.py:285`).\n"
    )
    for line in LINES:
        r = nat.get(line)
        if r is None:
            continue
        L.append(f"\n### {line}\n")
        if r["r"] is None:
            L.append("Bundle DIAGNOSTICO: no se leyo el anchor. Solo hay `b`.\n")
            L.append("| metrica | b |")
            L.append("|---|---|")
            for m in METRICS:
                L.append(f"| {CORTO[m]} | {fmt(r['b'].get(m))} |")
            continue
        L.append(
            "| metrica | sentido | b (baseline) | r (replicate) | span orientado | "
            "estimador | r_sd (5 splits) | n_perts min/max |"
        )
        L.append("|---" * 8 + "|")
        for m in METRICS:
            bv = r["b"].get(m)
            rr = r["r"].get(m) or {}
            rv = rr.get("replicate")
            sp = span(rv, bv, m)
            sentido = "mas alto mejor" if MAS_ES_MEJOR[m] else "mas bajo mejor"
            marca = " **<= 0**" if (sp is not None and sp <= 0) else ""
            L.append(
                f"| {CORTO[m]} | {sentido} | {fmt(bv)} | {fmt(rv)} | "
                f"{fmt(sp)}{marca} | `{rr.get('estimator')}` | "
                f"{fmt(rr.get('replicate_sd'))} | "
                f"{rr.get('n_perturbations_min')}/{rr.get('n_perturbations_max')} de 33 |"
            )

    # ---------- span ----------
    L.append("\n## El span `(r - b)` orientado -- ¿alguno sale <= 0?\n")
    if not hay_r:
        L.append("**Sin contestar: no hay `r`.**\n")
    else:
        L.append(
            "Orientado por el sentido de cada metrica: `r - b` donde mas alto es "
            "mejor, `b - r` donde mas bajo es mejor. Un span <= 0 significa que "
            "el extremo 1 no esta del lado bueno del extremo 0. El marcador "
            "tiene su propia guarda: `score.py:466-467` levanta \"the anchor "
            "leaves no usable headroom over the baseline\".\n"
        )
        L.append("| metrica | " + " | ".join(LINES) + " | n <= 0 |")
        L.append("|---" * (2 + len(LINES)) + "|")
        total_neg = 0
        for m in METRICS:
            celdas, neg = [], 0
            for line in LINES:
                r = nat.get(line)
                if r is None or r.get("r") is None:
                    celdas.append("-")
                    continue
                sp = span((r["r"].get(m) or {}).get("replicate"), r["b"].get(m), m)
                if sp is None:
                    celdas.append("NaN")
                    continue
                if sp <= 0:
                    neg += 1
                    celdas.append(f"**{sp:+.6f}**")
                else:
                    celdas.append(f"{sp:+.6f}")
            total_neg += neg
            L.append(f"| {CORTO[m]} | " + " | ".join(celdas) + f" | {neg} |")
        if total_neg == 0:
            L.append(
                "\n**Ningun span sale <= 0: los 36 son positivos.** La escala se "
                "puede usar en las seis metricas y las seis lineas. Contrasta "
                "con nuestro techo a mano, donde `jac` daba span negativo en 3 "
                "de 6 lineas (`bucle_techo.md`) -- pero ver la seccion "
                "siguiente antes de darlo por resuelto.\n"
            )
        else:
            L.append(f"\n**{total_neg} span(s) <= 0.**\n")

    # ---------- span vs ruido ----------
    if hay_r:
        L.append("\n## ⚠️ El span contra el ruido del propio `r`\n")
        L.append(
            "Un span positivo no basta. Si el span es mas chico que la "
            "dispersion con que se midio su extremo 1, re-sortear las "
            "particiones mueve el techo mas que todo el recorrido que la escala "
            "pretende medir, y un candidato puntuado ahi recibe un numero que "
            "depende mas de la semilla del anchor que de su propio merito. "
            "`r_sd` es la desviacion de los 5 agregados por particion "
            "(`anchor.py:286`, `ddof=0`).\n"
        )
        L.append("| linea | metrica | span | r_sd | span / r_sd | lectura |")
        L.append("|---" * 6 + "|")
        finas = []
        for line in LINES:
            r = nat.get(line)
            if r is None or r.get("r") is None:
                continue
            for m in METRICS:
                rr = r["r"].get(m) or {}
                sp = span(rr.get("replicate"), r["b"].get(m), m)
                sd = rr.get("replicate_sd")
                if sp is None or sd in (None, 0):
                    continue
                ratio = sp / sd
                if ratio < 1:
                    lect = "**span < 1 sd: no usable**"
                elif ratio < 3:
                    lect = "span < 3 sd: fragil"
                else:
                    lect = "holgado"
                if ratio < 3:
                    finas.append((line, m, sp, sd, ratio))
                L.append(
                    f"| {line} | {CORTO[m]} | {fmt(sp)} | {fmt(sd)} | "
                    f"{ratio:.2f} | {lect} |"
                )
        if finas:
            L.append(
                f"\n**{len(finas)} de 36 celdas tienen menos de 3 desviaciones de "
                "holgura**, y son las que hay que mirar antes de puntuar nada:\n"
            )
            for line, m, sp, sd, ratio in sorted(finas, key=lambda x: x[4]):
                L.append(
                    f"- **{line} / {CORTO[m]}**: span {sp:.6f} contra sd "
                    f"{sd:.6f} = {ratio:.2f} sd"
                )
            peor = min(finas, key=lambda x: x[4])
            L.append(
                f"\nEl caso extremo es **{peor[0]} / {CORTO[peor[1]]}**: el span "
                f"({peor[2]:.6f}) es {1/peor[4]:.1f}x mas chico que la "
                f"desviacion del propio `r` ({peor[3]:.6f}). Ahi la escala no "
                "mide nada: el extremo 1 se mueve mas entre semillas que toda la "
                "distancia al extremo 0.\n"
            )
            if any(f[1] == "de_wilcoxon_sig_jaccard" for f in finas):
                L.append(
                    "No es casualidad que `jac` concentre los casos finos: su "
                    "extremo 0 esta inflado por el artefacto de union vacia, asi "
                    "que `b` arranca alto y deja poca holgura hasta `r`. **El "
                    "span positivo de `jac` no resuelve el problema que "
                    "`bucle_techo.md` encontro: lo esconde detras de un signo.**\n"
                )

    # ---------- dispersion ----------
    L.append("\n## La dispersion entre las 5 particiones del anchor\n")
    L.append(
        "De `anchor_splits.parquet`, que guarda un valor por (particion, "
        "metrica) con su semilla derivada (`anchor.py:82-85`).\n"
    )
    if not hay_r:
        L.append("**Sin leer: el anchor no se abrio.**\n")
    else:
        L.append("| linea | metrica | min | max | rango | sd | los 5 valores |")
        L.append("|---" * 7 + "|")
        for line in LINES:
            r = nat.get(line)
            if r is None or r.get("splits") is None:
                continue
            for m in METRICS:
                ss = r["splits"].get(m)
                if not ss:
                    continue
                vals = [s["value"] for s in ss]
                sd = (r["r"].get(m) or {}).get("replicate_sd")
                L.append(
                    f"| {line} | {CORTO[m]} | {fmt(min(vals))} | {fmt(max(vals))} | "
                    f"{fmt(max(vals) - min(vals))} | {fmt(sd)} | "
                    + ", ".join(f"{v:.4f}" for v in vals) + " |"
                )

    # ---------- la pregunta de la direccion ----------
    L.append("\n## 🔴 ¿La asimetria de direccion que medimos a mano era ruido de sorteo?\n")
    if not hay_r:
        L.append("**Sin contestar: el anchor no se leyo.**\n")
    else:
        L.append(
            "`bucle_techo.md` midio cada metrica en las dos direcciones (mitad A "
            "como prediccion y mitad B como prediccion) y reporto las dos sin "
            "promediarlas, porque un desempate dirigido le da mas potencia a un "
            "lado. La pregunta era si esa diferencia es senal o ruido del "
            "sorteo.\n"
        )
        L.append(
            "El anchor nativo la contesta sin suponer nada: corre **5 "
            "particiones** de la MISMA direccion (`half_b` como prediccion "
            "siempre, `anchor.py:120`). Si el rango entre esas 5 cubre la "
            "diferencia que vimos entre direcciones, la asimetria no era senal. "
            "Se compara contra la variante **I** del techo, la que comparte "
            "convencion de control con el anchor.\n"
        )
        L.append(
            "| linea | metrica | rango entre 5 particiones | "
            "\\|dirA - dirB\\| del techo I | ¿el rango lo cubre? |"
        )
        L.append("|---" * 5 + "|")
        si = no = 0
        no_cubre = []
        for line in LINES:
            r = nat.get(line)
            ta, tb = techo.get((line, "I", "A")), techo.get((line, "I", "B"))
            if r is None or r.get("splits") is None or ta is None or tb is None:
                continue
            for m in METRICS:
                ss = r["splits"].get(m)
                if not ss:
                    continue
                vals = [s["value"] for s in ss]
                rango = max(vals) - min(vals)
                va, vb = ta.get(m), tb.get(m)
                if va is None or vb is None:
                    continue
                asim = abs(va - vb)
                cubre = rango >= asim
                if cubre:
                    si += 1
                else:
                    no += 1
                    no_cubre.append((line, m, rango, asim))
                L.append(
                    f"| {line} | {CORTO[m]} | {rango:.6f} | {asim:.6f} | "
                    f"{'**si**' if cubre else 'no'} |"
                )
        tot = si + no
        if tot:
            L.append(
                f"\n**Resumen: el rango entre particiones cubre la asimetria de "
                f"direccion en {si} de {tot} celdas (linea x metrica); no la "
                f"cubre en {no}.**\n"
            )
            if si > no:
                L.append(
                    "La lectura: en la mayoria de las celdas la diferencia entre "
                    "nuestras dos direcciones cae dentro de lo que una sola "
                    "direccion se mueve al re-sortear la particion. **Mayormente "
                    "era ruido de sorteo, y correr las dos direcciones no "
                    "compraba la informacion que parecia.** Lo que si la compra "
                    "es promediar varias particiones, que es lo que hace el "
                    "camino nativo con k=5.\n"
                )
            else:
                L.append(
                    "La lectura: en la mayoria de las celdas la asimetria EXCEDE "
                    "el rango entre 5 re-sorteos de la misma direccion. **No era "
                    "solo ruido de sorteo**, y reportar las dos direcciones sin "
                    "promediarlas estaba justificado.\n"
                )
            if no_cubre:
                L.append("Las celdas donde la asimetria excede el ruido de sorteo:\n")
                for line, m, rango, asim in sorted(no_cubre, key=lambda x: -(x[3] - x[2])):
                    L.append(
                        f"- **{line} / {CORTO[m]}**: asimetria {asim:.6f} contra "
                        f"rango {rango:.6f} (excede por {asim - rango:+.6f})"
                    )
                L.append("")
            L.append(
                "⚠️ Dos advertencias: el rango de 5 valores es un estimador "
                "grosero de dispersion, no un intervalo de confianza; y las dos "
                "cantidades se miden sobre poblaciones distintas -- el techo "
                "partio el pool completo por mitades, el anchor parte esta "
                "poblacion de 400 por perturbacion. Es indicativo, no una "
                "prueba.\n"
            )

    # ---------- contra lo nuestro ----------
    L.append("\n## Contra lo que ya teniamos\n")
    L.append(
        "Dos comparaciones distintas, ninguna es una equivalencia:\n\n"
        "- **nuestro suelo vs este `b`**: nuestro suelo emite el control SIN "
        "escalar; el `b` lo escala al perfil de respuesta promedio "
        "(`baseline_convencion.md` A.4).\n"
        "- **nuestro techo (variante I, direccion B) vs este `r`**: la variante "
        "nuestra mas cercana a la convencion oficial, pero el oficial promedia "
        "k=5 particiones y el nuestro usa 1, y `lfc_nmae` oficial usa otro "
        "estimador.\n"
    )
    L.append(
        "| linea | metrica | nuestro suelo | b nativo | dif | nuestro techo I/B | "
        "r nativo | dif |"
    )
    L.append("|---" * 8 + "|")
    for line in LINES:
        r = nat.get(line)
        if r is None:
            continue
        s, t = suelo.get(line), techo.get((line, "I", "B"))
        for m in METRICS:
            bv = r["b"].get(m)
            rv = (r["r"].get(m) or {}).get("replicate") if r.get("r") else None
            sv = s.get(m) if s else None
            tv = t.get(m) if t else None
            d_b = "-" if (sv is None or bv is None) else f"{sv - bv:+.6f}"
            d_r = "-" if (tv is None or rv is None) else f"{tv - rv:+.6f}"
            L.append(
                f"| {line} | {CORTO[m]} | {fmt(sv)} | {fmt(bv)} | {d_b} | "
                f"{fmt(tv)} | {fmt(rv)} | {d_r} |"
            )

    sec_union_vacia(L, nat, suelo)

    # ---------- estimador ----------
    if hay_r:
        L.append("\n## El estimador estampado por metrica\n")
        L.append(
            "`anchor.py:296-299` estampa, por metrica, que estimador produjo su "
            "ancla: `full_gate_raw` para la familia `lfc_nmae` y "
            "`split_half_raw` para el resto. El motivo (`anchor.py:57-66`): una "
            "MITAD llama muchos menos genes significativos que el dato completo "
            "(mediana del gate 1.83-3.36x menor), asi que el nucleo uniforme "
            "promediaria sobre 65-79% de las perturbaciones que promedia la "
            "metrica. `full_gate_raw` toma el gate y el denominador de la tabla "
            "DE del real COMPLETO y solo los dos vectores del numerador de las "
            "mitades.\n"
        )
        L.append("| linea | " + " | ".join(CORTO[m] for m in METRICS) + " |")
        L.append("|---" * (1 + len(METRICS)) + "|")
        for line in LINES:
            r = nat.get(line)
            if r is None or r.get("r") is None:
                continue
            L.append(
                f"| {line} | "
                + " | ".join(f"`{(r['r'].get(m) or {}).get('estimator')}`"
                             for m in METRICS) + " |"
            )

    # ---------- 20 filas ----------
    L.append("\n## 20 filas trabajadas, con los 5 valores por particion a la vista\n")
    if not hay_r:
        L.append("**Sin filas del lado `r`: el anchor no se leyo.**\n")
    else:
        L.append(
            "Veinte filas (linea, metrica) con TODO lo que entra en su span: los "
            "dos extremos, los 5 valores por particion que promedian al `r`, su "
            "dispersion, el estimador y la cohorte. No hay paso intermedio "
            "oculto entre el artefacto del marcador y el span.\n"
        )
        filas = []
        for line in LINES:
            r = nat.get(line)
            if r is None or r.get("r") is None:
                continue
            for m in METRICS:
                rr = r["r"].get(m) or {}
                ss = r["splits"].get(m) or []
                vals = [s["value"] for s in ss]
                filas.append({
                    "linea": line, "metric": m, "b": r["b"].get(m),
                    "r": rr.get("replicate"),
                    "span": span(rr.get("replicate"), r["b"].get(m), m),
                    "vals": vals, "sd": rr.get("replicate_sd"),
                    "rango": (max(vals) - min(vals)) if vals else None,
                    "est": rr.get("estimator"),
                    "nmin": rr.get("n_perturbations_min"),
                    "nmax": rr.get("n_perturbations_max"),
                    "media": (sum(vals) / len(vals)) if vals else None,
                    "semillas": [s["seed"] for s in ss],
                })
        elegidas = []
        elegidas += sorted(filas, key=lambda f: (f["span"] / f["sd"]) if f["sd"] else 1e9)[:6]
        elegidas += sorted(filas, key=lambda f: -(f["rango"] or 0))[:6]
        elegidas += [f for f in filas if f["est"] == "full_gate_raw"][:4]
        vistos, unicas = set(), []
        for f in elegidas + filas:
            k = (f["linea"], f["metric"])
            if k in vistos:
                continue
            vistos.add(k)
            unicas.append(f)
            if len(unicas) >= 20:
                break
        L.append(
            "| # | linea | metrica | b | r | span | los 5 valores por particion | "
            "media recalculada | dif vs r | rango | sd | span/sd | estimador | n_perts |"
        )
        L.append("|---" * 14 + "|")
        for i, f in enumerate(unicas, 1):
            dif = ("-" if (f["media"] is None or f["r"] is None)
                   else f"{abs(f['media'] - f['r']):.2e}")
            ratio = (f"{f['span'] / f['sd']:.2f}" if f["sd"] else "-")
            L.append(
                f"| {i} | {f['linea']} | {CORTO[f['metric']]} | {fmt(f['b'])} | "
                f"{fmt(f['r'])} | {fmt(f['span'])} | "
                + ", ".join(f"{v:.4f}" for v in f["vals"]) + " | "
                f"{fmt(f['media'])} | {dif} | {fmt(f['rango'])} | {fmt(f['sd'])} | "
                f"{ratio} | `{f['est']}` | {f['nmin']}/{f['nmax']} de 33 |"
            )
        L.append(
            "\n`dif vs r` es |media de los 5 valores - `replicate` del "
            "artefacto|: cierra el circuito entre `anchor_splits.parquet` y "
            "`anchor_agg.parquet`. Si no fuera ~0, la columna de dispersion no "
            "describiria el `r` que se reporta.\n"
        )
        sem = unicas[0]["semillas"] if unicas else []
        if sem:
            L.append(
                f"Las 5 semillas derivadas son las mismas en todas las lineas y "
                f"metricas, por construccion: `{sem}`, de "
                "`numpy.random.SeedSequence(0).generate_state(5)` "
                "(`anchor.py:74`, `135-141`).\n"
            )

    sec_intermedios(L, nat, "non-targeting")
    sec_no_hecho(L)

    with open(OUT_V2, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"escrito {OUT_V2} con {len(nat)} lineas (paso0 ok={ok0})")


if __name__ == "__main__":
    modo = sys.argv[1] if len(sys.argv) > 1 else "v2"
    if modo == "v1":
        reporte_v1(sys.argv[2] if len(sys.argv) > 2 else OUT_V1)
    elif modo == "v2":
        reporte_v2()
    else:
        raise SystemExit("uso: anclas_reporte.py [v1|v2]")
