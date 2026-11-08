#!/usr/bin/env python
"""Escribe modelo/out/dispersion_b.md desde dispersion_b*.json.

Lee y cuenta. No mide nada nuevo, no sustituye el `b` de ningun bundle y no
propone ningun umbral: calcula la resolucion, no decide un corte.
"""
import json
import os
import sys

import numpy as np

OUT_MD = "/home/ubuntu/vcc/modelo/out/dispersion_b.md"
DISP_B = ["/home/ubuntu/vcc/modelo/out/dispersion_b.json",
          "/home/ubuntu/vcc/modelo/out/dispersion_b_p2.json",
          "/home/ubuntu/vcc/modelo/out/dispersion_b_p3.json"]
DISP_U = "/home/ubuntu/vcc/modelo/out/dispersion_u.json"
NATIVAS = "/home/ubuntu/vcc/modelo/out/anclas_nativas_comp.json"

LINEAS = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
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
JAC = "de_wilcoxon_sig_jaccard"
FID = "de_wilcoxon_direction_fidelity_yield_raw"
N_SEM = 8
K_REPLICAS = [1, 2, 4, 8]


def leer_jsonl(path):
    out = []
    if not os.path.exists(path):
        return out
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
    recs = []
    for p in DISP_B:
        recs += leer_jsonl(p)
    idx = {(r["linea"], r["j"]): r for r in recs if not r.get("fallo_en")}
    faltan = [(ln, j) for ln in LINEAS for j in range(N_SEM)
              if (ln, j) not in idx]
    if faltan:
        sys.stderr.write(f"faltan {len(faltan)} corridas: {faltan[:6]}\n")
        sys.exit(1)

    urecs = leer_jsonl(DISP_U)
    uidx = {(r["linea"], r["j"]): r for r in urecs if not r.get("fallo_en")}
    nat = {r["linea"]: r for r in leer_jsonl(NATIVAS)}

    L = [
        "# La dispersion de `b`",
        "",
        "El tercer termino del margen, que hasta hoy faltaba. `dispersion_u.md` §2 dejo",
        "escrito que su margen OMITE el termino de `b` porque `b` era UNA sola",
        "realizacion con `--seed 0`, y que mientras lo fuera todo margen seria cota",
        "inferior. Aqui se estima: `cell-eval2 baseline` con **8 semillas** por linea",
        "sobre la MISMA poblacion real. 48 corridas, ninguna fallida.",
        "",
        "Con este termino, la resolucion de cada celda deja de ser una cota y pasa a ser",
        "un numero.",
        "",
        "## 1. Lo que vario y lo que no, con la prueba",
        "",
        "**Vario solo la semilla del baseline.** Familia nueva `20265000 + 100*fold + j`,",
        "`j = 0..7`. Ninguna colision con `2026000x` (verdad), `2026100x` (suelo),",
        "`2026200x`/`2026300x` (techo), `2026400x` (la dispersion de `u`) ni la `0` del",
        "baseline/anchor.",
        "",
        "**No vario la poblacion real.** `build_fold_adatas` sin modificar con",
        "`seed_real = 20260001 + fold_idx`, etiqueta `non-targeting`. Se comprueba con el",
        "mismo `sha256` de `dispersion_u.md` (matriz `data`/`indices`/`indptr` mas las",
        "etiquetas), y el script aborta si no coincide.",
        "",
        "🔴 **Y la huella es la MISMA que la de `dispersion_u.md`**, no solo consistente",
        "dentro de esta tanda. Eso es lo que permite componer `sd_b` con `sd_s` en un solo",
        "margen: los dos terminos se midieron sobre la misma poblacion, celda por celda.",
        "",
    ]
    rows = []
    for ln in LINEAS:
        hs = {idx[(ln, j)]["huella_real"] for j in range(N_SEM)}
        hp = {idx[(ln, j)]["huella_pred_baseline"] for j in range(N_SEM)}
        hu = {uidx[(ln, j)]["huella_real"] for j in range(N_SEM)
              if (ln, j) in uidx}
        ok = (len(hs) == 1 and len(hu) == 1 and hs == hu)
        hashes = {idx[(ln, j)]["hash_regla"]["obtenido"] for j in range(N_SEM)}
        fps = {idx[(ln, j)]["source_fingerprint_strict"] for j in range(N_SEM)}
        rows.append([
            ln, N_SEM,
            f"{min(idx[(ln, j)]['seed_baseline'] for j in range(N_SEM))}.."
            f"{max(idx[(ln, j)]['seed_baseline'] for j in range(N_SEM))}",
            f"**{len(hs)}**", f"**{len(hp)}**",
            ("**si**" if ok else "**NO**"),
            f"**{len(hashes)}**",
            ("**si**" if fps == {True} else f"{sorted(map(str, fps))}"),
        ])
    tabla(L, ["linea", "corridas", "semillas del baseline",
              "huellas distintas del REAL", "huellas distintas del EMITIDO",
              "huella del real = la de dispersion_u.md",
              "config_hash distintos entre las 8",
              "source_fingerprint_strict"], rows)

    h0 = idx[(LINEAS[0], 0)]["hash_regla"]
    L += [
        "**Una sola huella del real por linea y ocho del h5ad que el baseline emite**",
        "(`--save-pred`): la poblacion esta fija y la emision del baseline cambia en las",
        "ocho. Es la condicion que hace que la dispersion medida tenga una sola fuente.",
        "",
        "**La semilla no entra en el `config_hash`.** `run_params.yaml` del baseline no",
        f"contiene `seed` -- es parametro de la emision, no de `EvalConfig` -- y el hash",
        f"de regla sale identico en las 8 semillas de cada linea: `{h0['obtenido'][:24]}...`.",
        "Ese valor **no** es el de competencia (`a5519017556705ab...`) y no tiene por que",
        "serlo: es el del brazo baseline, y coincide **bit a bit con el del `BASE_<linea>`",
        "oficial** que se uso para construir los bundles. Verificado recalculandolo desde",
        "los dos `run_params.yaml`. Lo que esto descarta es que cambiar la semilla cambie",
        "la configuracion; lo unico que cambia es el sorteo.",
        "",
        "`source_fingerprint_strict` sale `true` en las 48.",
        "",
        "## 2. El margen con sus TRES terminos, derivado",
        "",
        "La escala es `s = (u - b)/(r - b)` (`score.py:395-399`). Las tres cantidades son",
        "mediciones, no constantes. Las derivadas de primer orden, escritas:",
        "",
        "```",
        "s = (u - b) / (r - b)",
        "",
        "ds/du = 1/(r - b)",
        "",
        "ds/dr = -(u - b)/(r - b)^2            = -s/(r - b)",
        "",
        "ds/db = -1/(r - b) + (u - b)/(r - b)^2",
        "      = [-(r - b) + (u - b)] / (r - b)^2",
        "      = (u - r)/(r - b)^2             = (s - 1)/(r - b)",
        "```",
        "",
        "La forma compacta `ds/db = (s - 1)/(r - b)` dice dos cosas que la forma larga",
        "esconde:",
        "",
        "- En `s = 0` (un candidato pegado al baseline) la derivada vale `-1/(r - b)`, el",
        "  mismo tamano que `ds/du` y de signo opuesto. Ahi `b` pesa tanto como la propia",
        "  emision.",
        "- En `s = 1` (un candidato pegado al ancla) la derivada vale **0**. Un candidato",
        "  que alcanza el ancla es insensible a donde cayo el baseline.",
        "",
        "Las 48 corridas de `dispersion_u.md` dan `s` entre **-1.6197 y +0.5012** (el",
        "minimo es el salto de convencion de A549/`jac`, §5.2 de ese reporte); por celda,",
        "las medias corren de -0.878 a +0.480. En los dos casos `|s - 1| >= 0.50`, y llega",
        "a **2.62**. Como `|ds/du| = 1/|r - b|` y `|ds/db| = |s - 1|/|r - b|`, el factor que",
        "los separa es exactamente `|s - 1|`: el termino de `b` pesa **mas** que el de `u`",
        "en toda celda con `s < 0`, que son la mayoria de las medidas.",
        "",
        "Componiendo en cuadratura, con `u`, `r` y `b` tratados como independientes:",
        "",
        "```",
        "margen_celda = sqrt(  sd_s^2",
        "                    + ( |s|     * SE_r / |r - b| )^2",
        "                    + ( |s - 1| * sd_b / |r - b| )^2  )",
        "```",
        "",
        "🔴 **Y la asimetria entre `SE_r` y `sd_b`, que no es un descuido.** Para `r` entra",
        "`SE_r = sd_r/sqrt(5)` porque el artefacto del bundle guarda la MEDIA de 5",
        "particiones: el numero que usa el marcador ya es un promedio. Para `b` entra la",
        "`sd_b` ENTERA, sin dividir por `sqrt(8)`, porque el bundle guarda **un solo**",
        "sorteo con `--seed 0` -- las 8 semillas estiman cuanto se mueve ese sorteo, no lo",
        "promedian. Dividir por `sqrt(8)` seria describir un `b` promediado que no existe",
        "en ningun bundle.",
        "",
        "**Lo que sigue sin estimarse, declarado.** La poblacion real es una sola",
        "realizacion (`seed_real = 20260001 + fold_idx`) y los tres terminos se miden",
        "sobre ella. El margen de arriba es el ruido de los sorteos DENTRO de esa",
        "poblacion; no cubre cambiar la poblacion. Eso sigue siendo cota inferior respecto",
        "al ruido total, y queda dicho.",
        "",
        "## 3. Las 36 celdas: los 8 valores de `b`",
        "",
        "`b` es la columna `mean` de `baseline_agg.csv` de cada corrida. El `b` de semilla",
        "0 -- el que esta DENTRO de los 6 bundles de competencia y que no se sustituye ni",
        "se promedia -- se localiza dentro de la distribucion de las 8, marcado como",
        "`sem0`.",
        "",
    ]

    # ---- recoger celdas
    celdas = []
    for ln in LINEAS:
        for met, corto in METRICAS:
            bs = [idx[(ln, j)]["b"].get(met) for j in range(N_SEM)]
            bs = [x for x in bs if x is not None]
            b0 = nat[ln]["b"].get(met)
            rr = nat[ln]["r"].get(met)
            r = (rr.get("replicate") if isinstance(rr, dict) else rr)
            sp = nat[ln]["splits"].get(met)
            sd_r = float(np.std([s["value"] for s in sp], ddof=0)) if sp else None
            se_r = (sd_r / np.sqrt(5)) if sd_r is not None else None
            # `s` y su sd vienen de dispersion_u.md, misma poblacion
            ss = []
            for j in range(N_SEM):
                rq = uidx.get((ln, j))
                if rq is None:
                    continue
                v = (rq.get("scored") or {}).get("filas", {}).get(met, {})
                if v.get("from_replicate") is not None:
                    ss.append(float(v["from_replicate"]))
            media_b = float(np.mean(bs)) if bs else None
            sd_b1 = float(np.std(bs, ddof=1)) if len(bs) > 1 else None
            sd_b0 = float(np.std(bs, ddof=0)) if bs else None
            # donde cae sem0
            z0 = ((b0 - media_b) / sd_b1
                  if (b0 is not None and media_b is not None and sd_b1) else None)
            pct0 = (100.0 * sum(1 for x in bs if x < b0) / len(bs)
                    if (b0 is not None and bs) else None)
            # el pool de 9, marcado como tal
            pool = bs + ([b0] if b0 is not None else [])
            sd_pool = float(np.std(pool, ddof=1)) if len(pool) > 1 else None
            s_med = float(np.mean(ss)) if ss else None
            sd_s = float(np.std(ss, ddof=1)) if len(ss) > 1 else None
            span = (r - b0) if (r is not None and b0 is not None) else None
            t_r = (abs(s_med) * se_r / abs(span)
                   if (s_med is not None and se_r is not None and span) else None)
            t_b = (abs(s_med - 1.0) * sd_b1 / abs(span)
                   if (s_med is not None and sd_b1 is not None and span) else None)
            m_viejo = (float(np.hypot(sd_s, t_r))
                       if (sd_s is not None and t_r is not None) else None)
            m_nuevo = (float(np.sqrt(sd_s ** 2 + t_r ** 2 + t_b ** 2))
                       if (sd_s is not None and t_r is not None
                           and t_b is not None) else None)
            celdas.append({
                "linea": ln, "metrica": corto, "canon": met,
                "bs": bs, "b0": b0, "media_b": media_b,
                "sd_b1": sd_b1, "sd_b0": sd_b0,
                "rango_b": (max(bs) - min(bs)) if bs else None,
                "min_b": min(bs) if bs else None, "max_b": max(bs) if bs else None,
                "z0": z0, "pct0": pct0, "sd_pool9": sd_pool,
                "r": r, "span": span, "sd_r": sd_r, "se_r": se_r,
                "s_med": s_med, "sd_s": sd_s,
                "t_r": t_r, "t_b": t_b,
                "m_viejo": m_viejo, "m_nuevo": m_nuevo,
                "crece": ((m_nuevo / m_viejo) if (m_viejo and m_nuevo) else None),
            })

    rows = []
    for c in celdas:
        rows.append([
            c["linea"], f"`{c['metrica']}`", SENTIDO[c["canon"]],
            " / ".join(f"{x:.4f}" for x in c["bs"]),
            f(c["media_b"]), f(c["sd_b1"]), f(c["sd_b0"]), f(c["rango_b"]),
            f(c["min_b"], 4), f(c["max_b"], 4),
        ])
    tabla(L, ["linea", "metrica", "sentido", "los 8 valores de b", "media",
              "sd (ddof=1)", "sd (ddof=0)", "rango", "min", "max"], rows)

    L += [
        "### 3.1 Donde cae el `b` de semilla 0 -- el que esta en los bundles",
        "",
        "`z` son desviaciones desde la media de las 8; las columnas de conteo dicen",
        "cuantas de las 8 quedan por debajo, iguales y por encima de `sem0`. Se separan",
        "tres situaciones distintas que el percentil solo confundiria:",
        "",
        "- **fuera**: `sem0` es estrictamente mas extremo que las ocho. El bundle se",
        "  construyo con una realizacion que ningun otro sorteo reprodujo.",
        "- **empate**: `sem0` coincide con el valor de varias de las 8. Suele ser un `b`",
        "  exactamente 0, no una punta.",
        "- **dentro**: hay valores a los dos lados.",
        "",
        "🔴 **Y una columna mas, porque sin ella la `z` enganaria.** Donde las 8 semillas",
        "toman muy pocos valores distintos, la `sd` es casi 0 y la `z` se dispara sin que",
        "eso signifique una cola: significa que `sem0` cayo en **otro valor discreto**. La",
        "columna `valores distintos de las 8` lo marca, y en esas celdas la `z` no se debe",
        "leer como desviaciones de una nube.",
        "",
    ]
    rows = []
    n_punta = 0
    n_disc = 0
    for c in celdas:
        bs = c["bs"]
        b0 = c["b0"]
        n_bajo = sum(1 for x in bs if x < b0) if b0 is not None else None
        n_ig = sum(1 for x in bs if x == b0) if b0 is not None else None
        n_alto = sum(1 for x in bs if x > b0) if b0 is not None else None
        ndist = len(set(bs))
        discreta = ndist <= 2
        if discreta:
            n_disc += 1
        if n_bajo is None:
            sit = "n/d"
        elif n_ig and (n_bajo == 0 or n_alto == 0):
            sit = "empate"
        elif n_bajo == 0 or n_alto == 0:
            sit = "🔴 **fuera**"
            n_punta += 1
        else:
            sit = "dentro"
        c["situacion"] = sit
        c["discreta"] = discreta
        rows.append([
            c["linea"], f"`{c['metrica']}`", f(b0), f(c["media_b"]),
            f(c["sd_b1"]), ndist,
            (f"{c['z0']:+.2f}" if c["z0"] is not None else "n/d"),
            f"{n_bajo}/{n_ig}/{n_alto}",
            f(c["sd_pool9"]), sit,
            ("🔴 **si**" if discreta else ""),
        ])
    tabla(L, ["linea", "metrica", "b de semilla 0", "media de las 8",
              "sd de las 8", "valores distintos de las 8", "z de sem0",
              "debajo/igual/encima", "sd del pool de 9 (marcado)",
              "situacion de sem0", "sd degenerada"], rows)

    zs = [abs(c["z0"]) for c in celdas
          if c["z0"] is not None and not c["discreta"]]
    zdeg = [(c, abs(c["z0"])) for c in celdas
            if c["z0"] is not None and c["discreta"]]
    L += [
        f"**{n_punta} de 36 celdas tienen el `b` del bundle estrictamente fuera** del",
        "rango de las 8 semillas: el bundle se construyo con una realizacion que ningun",
        "otro sorteo reprodujo.",
        "",
        f"En las {len(zs)} celdas con `sd` no degenerada, `|z|` mediano "
        f"{np.median(zs):.2f} y maximo {max(zs):.2f} -- una dispersion normal. Las "
        f"{len(zdeg)} degeneradas son otra cosa y van aparte: "
        + ", ".join(f"{c['linea']}/`{c['metrica']}` (`z` = {z:+.0f}, "
                    f"{len(set(c['bs']))} valor(es) distinto(s) en las 8)"
                    for c, z in sorted(zdeg, key=lambda t: -t[1])[:4])
        + ". Ahi las 8 coinciden en un valor y `sem0` cayo en otro: la `z` es grande",
        "porque el denominador es ~0, no porque `sem0` este lejos en unidades de la",
        "metrica. §6 lo explica: son saltos de la convencion `1/33`.",
        "",
        "La columna `sd del pool de 9` cumple con que el `b` de semilla 0 **entre como una",
        "observacion mas**: es la `sd` de las 8 nuevas junto con esa. Esta al lado, no en",
        "lugar de la de 8, y no se usa para el margen: para el margen entra la `sd` de las",
        "8, porque `sem0` es el valor cuyo error se esta describiendo y no puede ayudar a",
        "estimar su propia desviacion. 🔴 **Ningun bundle se toco y ningun `b` se",
        "sustituyo por una media.**",
        "",
        "### 3.2 Cuanto vale esa punta, en `s` -- DIAGNOSTICO, no cifra de competencia",
        "",
        "🔴 **Marcado arriba y no en una nota: esta tabla NO es la puntuacion.** La `s` de",
        "competencia es la de la columna `s con el b del bundle`, que sale de los 6 bundles",
        "tal como estan y es la que `dispersion_u.md` midio. La segunda columna es un",
        "CONTRAFACTUAL: que `s` habria salido si el baseline hubiera caido en la media de",
        "sus 9 realizaciones en vez de en la que cayo. Sirve para una sola cosa -- poner en",
        "unidades de `s` lo que §3.1 dice en unidades de `b` -- y para nada mas. **Ningun",
        "bundle se modifico y ningun `b` se sustituyo por una media.**",
        "",
        "Se calcula exacto, `s = (u - b)/(r - b)` con las dos `b`, no con la derivada de",
        "§2: donde `|b - b0|` es grande respecto del span la aproximacion de primer orden",
        "se rompe, y precisamente esas son las celdas que interesan. `u` es la media de las",
        "8 emisiones de `dispersion_u.md` y `r` no se toca. El centro es la media de las **8",
        "semillas nuevas**, sin `sem0`, por el mismo motivo que la `sd` de §3.1: `sem0` es",
        "el valor cuyo desplazamiento se esta midiendo y no puede entrar en su propio",
        "centro.",
        "",
        "🔴 **Y el clamp, que si no se aplica da una columna falsa.** `catalog.py:415-416`",
        "acota `expr_mse_unbiased_capped_norm` a `[0, 1]`, y solo esa. Las dos columnas de",
        "`s` llevan el clamp en esas seis celdas, igual que lo lleva el marcador. La columna",
        "`s medida por score` es la que escribio `score --real-bundle` en las 48 corridas de",
        "`dispersion_u.md` (media de las 8), y esta al lado como comprobacion: la `s`",
        "rehecha a mano tiene que reproducirla, y donde no lo hace es porque el clamp",
        "intervino.",
        "",
    ]

    def clampa(canon, s):
        if s is None:
            return None
        if canon == "expr_mse_unbiased_capped_norm":
            return float(min(max(s, 0.0), 1.0))
        return s

    rows = []
    for c in celdas:
        b0, mb, r = c["b0"], c["media_b"], c["r"]
        us, ss = [], []
        for j in range(N_SEM):
            rq = uidx.get((c["linea"], j))
            if rq is None:
                continue
            if rq.get("u", {}).get(c["canon"]) is not None:
                us.append(float(rq["u"][c["canon"]]))
            v = (rq.get("scored") or {}).get("filas", {}).get(c["canon"], {})
            if v.get("from_replicate") is not None:
                ss.append(float(v["from_replicate"]))
        u = float(np.mean(us)) if us else None
        s_med_score = float(np.mean(ss)) if ss else None
        s0 = clampa(c["canon"], (u - b0) / (r - b0)
                    if (u is not None and b0 is not None and r is not None
                        and (r - b0)) else None)
        s1 = clampa(c["canon"], (u - mb) / (r - mb)
                    if (u is not None and mb is not None and r is not None
                        and (r - mb)) else None)
        d = (s1 - s0) if (s0 is not None and s1 is not None) else None
        cambia_signo = (s0 is not None and s1 is not None
                        and s0 * s1 < 0 and abs(s0) > 1e-6 and abs(s1) > 1e-6)
        c["s0_exacto"], c["s_media_b"], c["dif_cf"] = s0, s1, d
        c["s_med_score"] = s_med_score
        c["cambia_signo"] = cambia_signo
        rep = (abs(s0 - s_med_score)
               if (s0 is not None and s_med_score is not None) else None)
        rows.append([
            c["linea"], f"`{c['metrica']}`", f(u, 4), f(b0, 4), f(mb, 4), f(r, 4),
            (f"{s_med_score:+.4f}" if s_med_score is not None else "n/d"),
            (f"{s0:+.4f}" if s0 is not None else "n/d"),
            (f"{rep:.1e}" if rep is not None else "n/d"),
            (f"{s1:+.4f}" if s1 is not None else "n/d"),
            (f"**{d:+.4f}**" if d is not None else "n/d"),
            ("🔴 **cambia de signo**" if cambia_signo else ""),
        ])
    tabla(L, ["linea", "metrica", "u medio", "b del bundle (sem0)",
              "media de las 8", "r", "s medida por score (COMPETENCIA)",
              "s rehecha con el b del bundle", "\\|rehecha - medida\\|",
              "s con la b media (DIAGNOSTICO)", "cuanto se movia",
              "signo"], rows)

    flips = [c for c in celdas if c.get("cambia_signo")]
    grandes = sorted((c for c in celdas
                      if c.get("dif_cf") is not None and c.get("sd_s") is not None
                      and abs(c["dif_cf"]) > c["sd_s"]),
                     key=lambda c: -abs(c["dif_cf"]))
    L += [
        f"**Lo que esto mide.** En **{len(grandes)} de 36 celdas el sorteo del baseline",
        "mueve la `s` mas que toda la dispersion de la emision junta** (`|cuanto se "
        "movia| > sd_s` de §4). Las cinco mayores: "
        + ", ".join(f"{c['linea']}/`{c['metrica']}` ({c['dif_cf']:+.4f} contra "
                    f"sd_s {f(c['sd_s'], 4)})" for c in grandes[:5])
        + ".",
        "",
    ]
    if flips:
        suben = [c for c in flips if c["s0_exacto"] < 0]
        bajan = [c for c in flips if c["s0_exacto"] > 0]
        L += [
            f"🔴 **Y en {len(flips)} celdas la `s` cambia de SIGNO**: "
            + ", ".join(f"{c['linea']}/`{c['metrica']}` "
                        f"({c['s0_exacto']:+.4f} -> {c['s_media_b']:+.4f})"
                        for c in flips)
            + ".",
            "",
            f"El cambio va en las dos direcciones, no en una: en {len(bajan)} celdas "
            f"({', '.join(c['linea'] + '/`' + c['metrica'] + '`' for c in bajan)}) "
            "el `b` del bundle le da a la emision un puntaje positivo que con el `b` medio",
            f"seria negativo, y en {len(suben)} "
            f"({', '.join(c['linea'] + '/`' + c['metrica'] + '`' for c in suben)}) "
            "ocurre lo contrario. No hay sesgo de direccion: hay sensibilidad. **El signo",
            "de esas celdas lo decide el sorteo del baseline, no la emision**, y en las dos",
            "de `pds_cosine` la `s` de competencia esta a menos de 0.04 de cero, que es",
            "menos que lo que el sorteo la mueve.",
            "",
        ]
    ht = [c for c in celdas if c["linea"] == "HT29" and c["metrica"] == "jac"]
    if ht and ht[0]["dif_cf"] is not None:
        c = ht[0]
        L += [
            "**El caso que ya estaba senalado y ahora tiene causa.** `cable_suelo.md`",
            "marco HT29/`jac` como la celda donde una copia del control cobra `s` = +0.5012",
            "por convencion, y `descomposicion.md` y `huecos.md` la persiguieron por el lado",
            f"de la `E`. Aqui se ve de donde sale la mitad: el `b` del bundle "
            f"({f(c['b0'], 4)}) es **el minimo de las nueve realizaciones** "
            f"(media {f(c['media_b'], 4)}, maximo {f(c['max_b'], 4)}), y con la media la "
            f"misma emision cobraria {c['s_media_b']:+.4f} en vez de "
            f"{c['s0_exacto']:+.4f}. **{abs(c['dif_cf']) / abs(c['s0_exacto']) * 100:.0f}%",
            "de ese puntaje es la realizacion concreta del baseline**, no la emision ni la",
            "convencion. Las dos explicaciones no compiten: la convencion pone el `u` arriba",
            "y el sorteo del baseline pone el `b` abajo.",
            "",
        ]

    L += [
        "## 4. El margen completo, y cuanto crecia",
        "",
        "Los tres terminos de §2 celda por celda, contra el margen de `dispersion_u.md`",
        "que omitia el tercero. `sd_s` y `s` vienen de `dispersion_u.json` -- misma",
        "poblacion, verificada en §1.",
        "",
    ]
    rows = []
    for c in celdas:
        rows.append([
            c["linea"], f"`{c['metrica']}`",
            f(c["s_med"]), f(c["span"]),
            (f"{abs((c['s_med'] - 1.0) / c['span']):.4f}"
             if (c["s_med"] is not None and c["span"]) else "n/d"),
            f(c["sd_s"]), f(c["t_r"]), f(c["sd_b1"]), f(c["t_b"]),
            f(c["m_viejo"]), f"**{f(c['m_nuevo'])}**",
            (f"{c['crece']:.2f}x" if c["crece"] else "n/d"),
        ])
    tabla(L, ["linea", "metrica", "s medio", "span r-b", "\\|ds/db\\|",
              "sd_s (termino u)", "termino de r", "sd_b", "termino de b",
              "margen de dispersion_u.md", "margen COMPLETO", "crece"], rows)

    crecen = [c for c in celdas if c["crece"]]
    if crecen:
        pk = max(crecen, key=lambda c: c["crece"])
        dom = [c for c in crecen if c["t_b"] and c["sd_s"] is not None
               and c["t_b"] > c["sd_s"]]
        L += [
            f"**El margen crece en {sum(1 for c in crecen if c['crece'] > 1.001)} de "
            f"{len(crecen)} celdas comparables**, por factores de "
            f"{min(c['crece'] for c in crecen):.2f}x a {pk['crece']:.2f}x. El maximo es "
            f"{pk['linea']}/`{pk['metrica']}`.",
            "",
            f"Las otras {len([c for c in celdas if not c['crece']])} no llevan factor "
            "porque su margen anterior era **exactamente 0** y no se puede dividir por el."
            " Son las `sd_s = 0` de `dispersion_u.md` §5, y ahi el cambio es mayor que en "
            "ninguna otra: pasan de 0 a un margen real de "
            + f"{min(c['m_nuevo'] for c in celdas if not c['crece'] and c['m_nuevo']):.4f}"
            + " a "
            + f"{max(c['m_nuevo'] for c in celdas if not c['crece'] and c['m_nuevo']):.4f}"
            + ". Una celda que parecia sin ruido tiene el ruido entero en `b`.",
            "",
            f"🔴 **Y en {len(dom)} celdas el termino de `b` es el MAYOR de los tres** -- "
            "mayor que la dispersion de la propia emision:",
            "",
        ]
        rows = []
        for c in sorted(dom, key=lambda c: -(c["t_b"] / (c["sd_s"] or 1e-12)))[:12]:
            rows.append([
                c["linea"], f"`{c['metrica']}`", f(c["sd_s"]), f(c["t_r"]),
                f"**{f(c['t_b'])}**",
                (f"{c['t_b'] / c['sd_s']:.1f}x" if c["sd_s"] else "inf"),
            ])
        tabla(L, ["linea", "metrica", "sd_s (u)", "termino de r",
                  "termino de b", "b / u"], rows)
        L += [
            "Es la consecuencia directa de §2: `|ds/db| = |s - 1|/|r - b|` y los `s`",
            "medidos son ~0 o negativos, asi que `|s - 1| >= 1` en casi todas las celdas",
            "mientras `|ds/du| = 1/|r - b|`. Con `sd_b` del orden de `sd_u`, el termino de",
            "`b` sale mayor que el de `u`. **El baseline era la fuente de ruido mas grande",
            "de la escala y estaba sin medir.**",
            "",
        ]

    L += [
        "## 5. La resolucion de cada celda",
        "",
        "El tamano minimo de diferencia entre dos candidatos que la celda distingue. Aqui",
        "`r` y `b` **se cancelan**, y eso cambia el resultado por completo:",
        "",
        "```",
        "s1 - s2 = (u1 - b)/(r - b) - (u2 - b)/(r - b) = (u1 - u2)/(r - b)",
        "```",
        "",
        "Los dos candidatos se puntuan contra el MISMO bundle, asi que `r` y `b` son",
        "comunes: su ruido mueve las dos `s` juntas y no entra en la diferencia. Con `k`",
        "replicas (semillas de emision) por candidato, promediadas:",
        "",
        "```",
        "SE( s1_media - s2_media ) = sd_s * sqrt(2/k)",
        "",
        "resolucion_k = sd_s * sqrt(2/k)",
        "```",
        "",
        "Es el margen de la DIFERENCIA, que es la cantidad que §3.16 compara. No se",
        "multiplica por ningun factor de confianza: el encargo pide la resolucion, no un",
        "corte, y la regla de decision ya esta escrita.",
        "",
        "🔴 **Y la diferencia con §4 que conviene no confundir.** El margen completo de §4",
        "es para situar UNA `s` en la escala absoluta (por ejemplo, decir si un candidato",
        "supera el ancla). La resolucion de aqui es para comparar DOS candidatos. La",
        "segunda es mucho mas pequena que la primera porque no carga ni `r` ni `b`. Las dos",
        "columnas estan para que se use la que corresponda a la pregunta.",
        "",
    ]
    rows = []
    for c in celdas:
        fila = [c["linea"], f"`{c['metrica']}`", f(c["sd_s"])]
        for k in K_REPLICAS:
            fila.append(f(c["sd_s"] * np.sqrt(2.0 / k)) if c["sd_s"] is not None
                        else "n/d")
        fila.append(f"**{f(c['m_nuevo'])}**")
        rows.append(fila)
    tabla(L, ["linea", "metrica", "sd_s"]
          + [f"resolucion k={k}" for k in K_REPLICAS]
          + ["margen absoluto de §4 (otra pregunta)"], rows)

    ceros = [c for c in celdas if c["sd_s"] == 0.0]
    c_mse = [c for c in ceros if c["canon"] == "expr_mse_unbiased_capped_norm"]
    c_otras = [c for c in ceros if c["canon"] != "expr_mse_unbiased_capped_norm"]
    L += [
        f"**{len(ceros)} celdas dan resolucion 0 en todo `k`**, y eso NO quiere decir que",
        "distingan diferencias arbitrariamente pequenas. Son exactamente las de",
        f"`dispersion_u.md` §5: las {len(c_mse)} de `mse_capped_norm`, donde el clamp a",
        "`[0,1]` (`catalog.py:415-416`) devuelve 0 para los dos candidatos si los dos estan",
        "por debajo de `b`, y la celda no distingue **nada**; y "
        + ", ".join(f"{c['linea']}/`{c['metrica']}`" for c in c_otras)
        + ", donde la emision no mueve el `u` crudo (`sd_u = 0`) y el cero si es",
        "estructural. Una resolucion de 0 por ceguera y una por estabilidad se escriben",
        "igual y significan lo contrario; esta columna no las separa y §5 de",
        "`dispersion_u.md` si.",
        "",
        "🔴 **Y lo que `b` le hace a esas mismas celdas.** Su resolucion sigue siendo 0 --",
        "comparar dos candidatos ahi no depende de `b`, que se cancela -- pero su margen",
        "ABSOLUTO deja de ser 0: §4 les da 0.0353 a 0.0891. `dispersion_u.md` escribio 0",
        "en las dos columnas porque solo tenia la primera. Con el tercer termino las dos",
        "dejan de coincidir, y es justo lo que debe pasar: una celda puede ser ciega a la",
        "diferencia entre candidatos y seguir estando mal situada en la escala absoluta.",
        "",
        "## 6. La `E` del baseline, sus cohortes y si alguna `sd` es Bernoulli",
        "",
        "`dispersion_u.md` §5.2 destapo que la `E` del SUELO si se mueve con su semilla, y",
        "que un solo gen predicho sobre una referencia vacia cuesta `1/33` de `jac` crudo,",
        "amplificado por `1/span`. La misma pregunta, ahora para el baseline.",
        "",
        "**Como se mide, y la validacion del proxy.** `baseline_results.csv` da el valor",
        "por perturbacion, asi que las cohortes de `fid`, `reach` y `lfc_nmae` son el",
        "conteo de filas no-NaN, y la `E` se lee como `jac == 1.0`. Que `jac == 1.0`",
        "equivalga a `|union| = 0` es un PROXY, y esta validado antes de usarse: contra los",
        "conteos de union de `huecos_cohorte.json` -- que salieron de `de_real.parquet` y",
        "`de_pred.parquet` con `--write-degenes` -- los dos conjuntos coinciden exactamente",
        "en las seis lineas con semilla 0 (E = 15/0/14/3/7/0). Ademas cada semilla lleva su",
        "guarda: union vacia implica `n_real = n_pred = 0`, luego `max(n_pred, n_conf) = 0`",
        "y `fid` sale NaN (`direction.py:776` con `_safe_div`, 689-699); asi que toda",
        "perturbacion con `jac == 1.0` tiene que tener `fid` NaN.",
        "",
    ]
    rows = []
    viol_tot = 0
    for ln in LINEAS:
        Es = [idx[(ln, j)]["cohortes"]["E_union_vacia"] for j in range(N_SEM)]
        fids = [idx[(ln, j)]["cohortes"]["n_fid"] for j in range(N_SEM)]
        rch = [idx[(ln, j)]["cohortes"]["n_reach"] for j in range(N_SEM)]
        lfc = [idx[(ln, j)]["cohortes"]["n_lfc"] for j in range(N_SEM)]
        conj = [tuple(idx[(ln, j)]["cohortes"]["perts_union_vacia"])
                for j in range(N_SEM)]
        v = sum(len(idx[(ln, j)]["cohortes"]["proxy_violaciones"])
                for j in range(N_SEM))
        viol_tot += v
        rows.append([
            ln, "/".join(map(str, Es)), len(set(Es)), len(set(conj)),
            "/".join(map(str, fids)), "/".join(map(str, rch)),
            "/".join(map(str, lfc)),
            ("**0**" if v == 0 else f"🔴 **{v}**"),
        ])
    tabla(L, ["linea", "E por semilla (j=0..7)", "valores distintos de E",
              "CONJUNTOS distintos de union vacia", "cohorte de fid por semilla",
              "cohorte de reach", "cohorte de lfc_nmae",
              "violaciones de la guarda del proxy"], rows)
    L += [
        f"Violaciones de la guarda en las 48 corridas: **{viol_tot}**. El proxy se sostiene",
        "en todas, no solo en la semilla 0.",
        "",
    ]

    # las perturbaciones que entran y salen de la union vacia del baseline
    L += [
        "### 6.1 Que perturbaciones mueve la semilla del baseline",
        "",
        "`huecos.md` §2.5 nombro BRD9, IRF7 y RAP1GAP (HT29) y ZBP1 (K562) como las que el",
        "`b` de semilla 0 saco de la union vacia. Con 8 semillas se ve si son siempre esas.",
        "El conjunto de referencia es el de `dispersion_u.md`: la union vacia de NUESTRO",
        "suelo, que es `33 - |cohorte del lado real|` y no depende de ninguna emision.",
        "",
    ]
    # el conjunto de referencia: la union vacia del SUELO, que es
    # 33 - |cohorte del lado real| y no depende de ninguna emision
    ref = {}
    for ln in LINEAS:
        conj_suelo = [tuple(uidx[(ln, j)]["cohortes"]["perts_union_vacia"])
                      for j in range(N_SEM) if (ln, j) in uidx]
        ref[ln] = set(max(set(conj_suelo), key=conj_suelo.count)) if conj_suelo \
            else set()
    # lo que saco la semilla 0, de los conteos de union de huecos_cohorte.json
    SEM0 = {"A549": [], "BXPC3": [], "HAP1": [],
            "HT29": ["BRD9", "IRF7", "RAP1GAP"], "K562": ["ZBP1"], "MCF7": []}
    rows = []
    for ln in LINEAS:
        sal = {}
        for j in range(N_SEM):
            s = set(idx[(ln, j)]["cohortes"]["perts_union_vacia"])
            for p in (ref[ln] - s):
                sal.setdefault(p, []).append(j)
        nombres = sorted(set(sal) | set(SEM0[ln]))
        if not nombres:
            rows.append([ln, len(ref[ln]), "-- ninguna --", "0/8", "no", "--"])
            continue
        for p in sorted(nombres, key=lambda q: (-len(sal.get(q, [])), q)):
            js = sal.get(p, [])
            rows.append([
                ln, len(ref[ln]), f"**{p}**", f"{len(js)}/8",
                ("**si**" if p in SEM0[ln] else "no"),
                "/".join(map(str, js)) if js else "--",
            ])
    tabla(L, ["linea", "E del suelo (referencia)",
              "perturbacion que el baseline saca de la union vacia",
              "en cuantas de las 8 semillas nuevas",
              "la saco tambien la semilla 0 del bundle", "en cuales (j)"], rows)

    L += [
        "🔴 **La respuesta: no son siempre esas.** El encargo preguntaba si las cuatro de",
        "`huecos.md` §2.5 son estables o cambian. Cambian, y de dos maneras distintas:",
        "",
        "- **HT29.** La semilla 0 saco **las tres a la vez** (BRD9, IRF7 y RAP1GAP), que es",
        "  `E = 3`, el **minimo** de las nueve realizaciones -- igualado por una sola de las",
        "  8 nuevas, `j=2`, que saca exactamente esas mismas tres. Las otras siete sacan",
        "  menos: `E` = 6/5/6/5/4/4/6. Entre las 8 esos nombres aparecen sueltos --",
        "  RAP1GAP en 3/8, BRD9 en 2/8, IRF7 en 1/8 -- y aparecen ademas **dos que la",
        "  semilla 0 no saco: STAT3 (2/8) y HERC5 (1/8)**. El",
        "  conjunto de candidatas es mayor que lo que revelaba una sola realizacion; lo que",
        "  el sorteo decide es cuales caen.",
        "- **K562.** La semilla 0 saco ZBP1. **Ninguna de las 8 nuevas saca nada**, ZBP1",
        "  incluido: las ocho dan `E = 8`, que es exactamente la `E` del suelo. Ese evento",
        "  es 1 de 9 sorteos.",
        "- **A549 y HAP1.** La semilla 0 no saco ninguna, pero las nuevas si: STAT5A (1/8)",
        "  en A549, BRD9 y MAP3K14 (1/8 cada una) en HAP1.",
        "",
        "**Y una coincidencia que vale anotar sin explicarla.** STAT5A es la unica que el",
        "baseline saca en A549, y es la misma que el SUELO saco en su `j=5`",
        "(`dispersion_u.md` §5.2). Las dos emisiones son distintas y los sorteos",
        "independientes; que la marginal de A549 sea la misma para las dos apunta a una",
        "propiedad de esa perturbacion en esa linea, no de la emision. Queda dicho como",
        "observacion, no como mecanismo: no se midio.",
        "",
        "**Lo que esto le hace a `huecos.md` §2.5.** Esa seccion midio bien lo que midio --",
        "con el `b` del bundle, esas cuatro celdas son las que rompen `E_u = E_b` -- pero",
        "la lectura de que BRD9/IRF7/RAP1GAP/ZBP1 fueran *las* perturbaciones donde el",
        "baseline se despega no se sostiene con nueve realizaciones a la vista. Son las que",
        "se despegaron en una.",
        "",
    ]

    # Bernoulli: sd acumulada de b
    L += [
        "### 6.2 La `sd` acumulada de `b` con k = 2..8: nube o Bernoulli",
        "",
        "`dispersion_u.md` §4 encontro celdas donde la `sd` arranca en 0 y salta: ahi la",
        "dispersion no es una nube que se estreche con mas draws, es un evento discreto. La",
        "misma comprobacion sobre `b`, en el orden `j = 0..7` fijado antes de medir.",
        "",
    ]
    rows = []
    bern = []
    for c in celdas:
        acc = []
        for k in range(2, N_SEM + 1):
            acc.append(float(np.std(c["bs"][:k], ddof=1)))
        arranca_cero = acc[0] == 0.0 and max(acc) > 0
        if arranca_cero:
            bern.append(c)
        rows.append([c["linea"], f"`{c['metrica']}`"]
                    + [("0" if a == 0.0 else f(a)) for a in acc]
                    + ["🔴 **si**" if arranca_cero else "no"])
    tabla(L, ["linea", "metrica"] + [f"k={k}" for k in range(2, N_SEM + 1)]
          + ["arranca en 0 y salta"], rows)
    if bern:
        L += [
            f"**{len(bern)} celdas arrancan en 0 y saltan**: "
            + ", ".join(f"{c['linea']}/`{c['metrica']}`" for c in bern)
            + ". Ahi la `sd_b` de 8 draws estima la probabilidad de un salto, no la",
            "anchura de una nube.",
            "",
            "Donde esto se tiene que leer con cuidado es en el **margen de §4**, que es el",
            "que carga `sd_b` -- no en la resolucion de §5, que solo carga `sd_s` porque `b`",
            "se cancela entre dos candidatos. Y las dos de K562 son el ejemplo mas claro del",
            "problema: su `sd_b` sale ~0.0002 porque las 8 semillas coinciden, asi que su",
            "termino de `b` en §4 es despreciable (0.0006 y 0.0005) -- y sin embargo §3.2",
            "mide que el `b` de semilla 0 mueve su `s` en 0.1751 y 0.1110, con cambio de",
            "signo en las dos. **El margen de primer orden no ve ese salto**: lo estima con",
            "una `sd` que describe ocho draws de acuerdo, mientras el noveno cayo en otro",
            "valor discreto. En esas celdas el contrafactual exacto de §3.2 es la cifra",
            "informativa y el termino de `b` de §4 no lo es.",
            "",
        ]
    else:
        L += [
            "**Ninguna celda arranca en 0 y salta.** A diferencia de `u`, el `b` se mueve de",
            "forma continua en las 36 celdas: la emision dispersa del baseline reescala",
            "todos los genes en cada sorteo, asi que no hay un umbral discreto que cruzar",
            "como el de la significancia DE que producia los saltos de `jac` en",
            "`dispersion_u.md` §5.2.",
            "",
        ]

    # ------------------------------------------------ 7. 20 filas trabajadas
    L += [
        "## 7. Las 20 filas trabajadas",
        "",
        "**Criterio de seleccion, declarado antes de mirar los valores**: las celdas de",
        "`pds_cosine`, `lfc_nmae` y `reach` en las seis lineas = 18 filas, mas las dos de",
        "`jac` de HT29 y K562 -- las dos lineas donde `huecos.md` §2.5 habia identificado,",
        "ANTES de este encargo, que la `E` del `b` no coincide con la del suelo",
        "(BRD9/IRF7/RAP1GAP y ZBP1) = **20**. Es el mismo criterio de `dispersion_u.md`",
        "§7, para que las filas se puedan poner lado a lado.",
        "",
        "Cada fila trae los 8 `b` crudos, el `b` del bundle, y los tres terminos del margen",
        "con su aritmetica a la vista.",
        "",
    ]
    sel = [c for c in celdas if c["metrica"] in ("pds_cosine", "lfc_nmae", "reach")]
    sel += [c for c in celdas if c["metrica"] == "jac"
            and c["linea"] in ("HT29", "K562")]
    rows = []
    for c in sel:
        rows.append([
            c["linea"], f"`{c['metrica']}`",
            " / ".join(f"{x:.4f}" for x in c["bs"]),
            f(c["media_b"]), f(c["sd_b1"]), f(c["b0"]),
            (f"{c['z0']:+.2f}" if c["z0"] is not None else "n/d"),
            f(c["r"]), f(c["span"]), f(c["s_med"]),
            (f"{abs((c['s_med'] - 1.0) / c['span']):.4f}"
             if (c["s_med"] is not None and c["span"]) else "n/d"),
            f(c["sd_s"]), f(c["t_r"]), f(c["t_b"]),
            f"**{f(c['m_nuevo'])}**",
            f(c["sd_s"] * np.sqrt(2.0 / 8) if c["sd_s"] is not None else None),
        ])
    tabla(L, ["linea", "metrica", "los 8 b crudos", "media de b", "sd_b",
              "b del bundle (sem0)", "z de sem0", "r", "span r-b", "s medio",
              "\\|ds/db\\|", "sd_s", "termino de r", "termino de b",
              "margen COMPLETO", "resolucion k=8"], rows)
    L += [
        f"Son **{len(rows)}** filas.",
        "",
        "## Lo que no se hizo",
        "",
        "- **No se reconstruyo ningun bundle.** Los 6 de competencia de",
        "  `/datos/modelo/out/nativas_comp/` estan exactamente como estaban.",
        "- **No se sustituyo el `b` de los bundles por una media.** El `b` de semilla 0",
        "  entra como una observacion mas y marcado (§3.1); el margen usa la `sd` de las 8.",
        "- **No se construyo ningun candidato.** Aqui no hay emision de candidato: lo que",
        "  se emite es el baseline del propio marcador.",
        "- **No se propone ningun umbral.** §4 da el margen y §5 la resolucion con su",
        "  formula; la regla de decision ya esta escrita en §3.16 y no se toca.",
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
