#!/usr/bin/env python
"""Escribe modelo/out/candidato_01.md.

Lee y cuenta. No ajusta el candidato, no propone umbrales, no toca bundles.
"""
import json
import os
import sys

import numpy as np

OUT_MD = "/home/ubuntu/vcc/modelo/out/candidato_01.md"
CAND = ["/home/ubuntu/vcc/modelo/out/candidato_01.json",
        "/home/ubuntu/vcc/modelo/out/candidato_01_p2.json",
        "/home/ubuntu/vcc/modelo/out/candidato_01_p3.json"]
REFUTA = "/home/ubuntu/vcc/modelo/out/cand01_refuta.json"
DISP_U = "/home/ubuntu/vcc/modelo/out/dispersion_u.json"
DISP_B = ["/home/ubuntu/vcc/modelo/out/dispersion_b.json",
          "/home/ubuntu/vcc/modelo/out/dispersion_b_p2.json",
          "/home/ubuntu/vcc/modelo/out/dispersion_b_p3.json"]
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
MSE = "expr_mse_unbiased_capped_norm"
N_SEM = 8
KS = [1, 2, 4, 8]


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
    for p in CAND:
        recs += leer_jsonl(p)
    idx = {(r["linea"], r["j"]): r for r in recs if not r.get("fallo_en")}
    faltan = [(l, j) for l in LINEAS for j in range(N_SEM) if (l, j) not in idx]
    if faltan:
        sys.stderr.write(f"faltan {len(faltan)}: {faltan[:6]}\n")
        sys.exit(1)
    ref = json.load(open(REFUTA))
    nat = {r["linea"]: r for r in leer_jsonl(NATIVAS)}
    uidx = {(r["linea"], r["j"]): r for r in leer_jsonl(DISP_U)
            if not r.get("fallo_en")}
    brecs = []
    for p in DISP_B:
        brecs += leer_jsonl(p)
    bidx = {(r["linea"], r["j"]): r for r in brecs if not r.get("fallo_en")}

    L = [
        "# Candidato 01: respuesta lineal a lo largo de la covarianza del control",
        "",
        "Su trabajo era uno solo: mostrar habilidad **especifica del blanco**, por encima",
        "de un `b` que emite el control escalado a UN perfil de respuesta promedio, el",
        "mismo para las 33 perturbaciones.",
        "",
        "**No la muestra.** La refutacion esta medida en las seis lineas y es el resultado",
        "de este encargo. La medicion completa se corrio igual, por tres razones que se",
        "declaran: el punto 5 (como mueve las cohortes una emision con mucha senal) es",
        "independiente de que el candidato tenga habilidad y cierra el hueco que §3.17",
        "dejo abierto; '¿supera al `b`?' merece quedar registrado con su margen; y la",
        "dispersion de un candidato real no se habia medido nunca. **No se itero ni se",
        "ajusto nada mirando la `s`.**",
        "",
        "## 1. 🔴 LA REFUTACION, primero porque es el resultado",
        "",
        "**La hipotesis.** Para el blanco `g`, los genes que mas covarian con `g` en el",
        "pool de control de ese contexto son los candidatos a moverse cuando `g` se",
        "perturba.",
        "",
        "**La medicion.** Spearman del delta de pseudobulk contra el control, sobre los",
        "genes expresados (>= 1 en CP50k en los DOS controles). Cuatro cantidades, y las",
        "dos ultimas son las que impiden una conclusion falsa:",
        "",
        "| clave | que es | para que |",
        "|---|---|---|",
        "| **A** | mitad-contra-mitad del REAL | el **techo alcanzable** a 400 celulas. Si fuera ~0, el test no tendria potencia y no se podria concluir NADA |",
        "| **B** | real de un blanco vs real de OTRO | la respuesta **generica** compartida. Es la cantidad que el `b` captura |",
        "| **C** | la prediccion vs el real del MISMO blanco | lo que se quiere medir |",
        "| **D** | la prediccion vs el real de OTRO blanco | **control negativo**. `C - D` es lo especifico del blanco que el candidato captura |",
        "",
    ]
    rows = []
    for ln in LINEAS:
        v = ref[ln]
        sg = [x for x in v["concordancia_signo"] if x is not None]
        A = float(np.median(v["A_mitades_real"]))
        B = float(np.median(v["B_cruzado_real"]))
        C = float(np.median(v["C_pred_vs_mismo"]))
        D = float(np.median(v["D_pred_vs_otro"]))
        rows.append([
            ln, v["n_genes_expresados"], v["n_control_cells"],
            f"{A:+.4f}", f"{B:+.4f}", f"**{C:+.4f}**", f"{D:+.4f}",
            f"**{C - D:+.4f}**", f"**{np.mean(sg):.4f}**",
        ])
    tabla(L, ["linea", "genes expresados", "celulas de control",
              "A techo", "B generica", "C candidato", "D control negativo",
              "C - D", "concordancia de signo"], rows)

    As = [float(np.median(ref[l]["A_mitades_real"])) for l in LINEAS]
    Bs = [float(np.median(ref[l]["B_cruzado_real"])) for l in LINEAS]
    Cs = [float(np.median(ref[l]["C_pred_vs_mismo"])) for l in LINEAS]
    Ds = [float(np.median(ref[l]["D_pred_vs_otro"])) for l in LINEAS]
    sgs = [float(np.mean([x for x in ref[l]["concordancia_signo"]
                          if x is not None])) for l in LINEAS]
    L += [
        f"**`A` corre de {min(As):+.4f} a {max(As):+.4f} en las seis lineas: el test tiene",
        "potencia.** A 400 celulas por perturbacion el delta real es reproducible, asi que",
        "un modelo con senal tendria donde mostrarla.",
        "",
        f"**`C` corre de {min(Cs):+.4f} a {max(Cs):+.4f}: cero.** Y la concordancia de signo",
        "-- la fraccion de genes que el lado real mueve de verdad (`|dlog| > 0.25`) en que",
        f"el signo predicho acierta -- va de {min(sgs):.4f} a {max(sgs):.4f}. **0.5 es azar",
        "exacto.** El candidato no acierta el signo de la respuesta mejor que una moneda.",
        "",
        f"**Y `D` descarta la explicacion alternativa.** `C - D` va de {min(c - d for c, d in zip(Cs, Ds)):+.4f} "
        f"a {max(c - d for c, d in zip(Cs, Ds)):+.4f}: el candidato no captura ni lo",
        "especifico del blanco ni la respuesta generica. Si hubiera capturado lo generico,",
        "`D` habria salido alto y `C - D` bajo; salen los dos en cero.",
        "",
        "🔴 **Y el hallazgo que sale de paso y vale mas que el candidato: `B > A` en cinco",
        "de las seis lineas.** La respuesta compartida entre perturbaciones distintas es",
        "mas fuerte que la reproducibilidad de una sola perturbacion consigo misma:",
        "",
    ]
    rows = []
    for ln, a, b in zip(LINEAS, As, Bs):
        rows.append([ln, f"{a:+.4f}", f"{b:+.4f}", f"{b - a:+.4f}",
                     ("**si**" if b > a else "no")])
    tabla(L, ["linea", "A (una perturbacion consigo misma)",
              "B (entre perturbaciones distintas)", "B - A", "B > A"], rows)
    L += [
        "Medido en contexto, sin tocar el `b` ni el anchor. Explica mecanicamente por que",
        "el `b` es un comparador tan duro: **esta capturando la componente dominante de la",
        "varianza**, y la parte especifica del blanco que queda por encima es la minoritaria",
        "incluso en el lado real. Un candidato que quiera superar al `b` tiene que acertar",
        "la componente pequena, no la grande.",
        "",
        "### 1.1 Las dos magnitudes que se midieron y se descartaron ANTES de puntuar",
        "",
        "Ninguna de las dos se puntuo nunca: se descartaron por una propiedad medida de la",
        "EMISION, no por su `s`. Queda escrito para que no se repitan.",
        "",
        "| version | formula | lo que se midio en A549 | por que se descarto |",
        "|---|---|---|---|",
        "| v1 | `beta_h = cov/var(y_g)`, `Delta_y_h = -beta_h * mean(y_g)` (apagar `g` del todo) | `r_min = 0.0000` con `r_max = 1.0000` en **24 de 33** blancos; cientos de genes por blanco a cero exacto (RAP1GAP 381, ZBP1 335, STAT4 332) | `beta_h = rho * sd_h/sd_g` crece sin cota cuando `sd_h >> sd_g`; la extrapolacion sale del soporte |",
        "| v2 | `Delta_y_h = -rho_h * sd(y_h)` (esperanza condicional a una `sd`) | la **mediana** de los genes que responden va a cero; **18 de 33** blancos saturados; `\\|log2FC\\|` de pseudobulk hasta **17.4** | para un gen disperso `mean/sd ~ sqrt(p/(1-p))` = 0.23 a `p = 0.05`, asi que cualquier `\\|rho\\| > 0.23` lo cruza por debajo de cero. log1p de cuentas dispersas esta inflado de ceros, no es una conjunta normal |",
        "| **v3, la que corre** | `r_h = 1 - rho_h`, acotada en `[0, 2]` | **0 de 33** saturados | adimensional y acotada por construccion; la unica cantidad que entra es `rho` en contexto, que es literalmente 'el orden y el signo, nunca el valor' |",
        "",
        "La magnitud queda declarada como **hueco**, no como resultado: que la fraccion de",
        "cambio sea `rho` es una eleccion de escala, no una cifra medida de cuanto baja un",
        "gen cuando se perturba su regulador. Eso exigiria el lado real perturbado, que",
        "esta prohibido.",
        "",
        "## 2. Lo que vario y lo que no, con la prueba",
        "",
        "**Vario solo la semilla de emision.** Familia nueva `20266000 + 100*fold + j`,",
        "`j = 0..7`. Ninguna colision con `2026000x` (verdad), `2026100x` (suelo),",
        "`2026200x`/`2026300x` (techo), `2026400x` (dispersion de `u`), `2026500x`",
        "(dispersion de `b`) ni la `0` del baseline/anchor.",
        "",
        "**No vario la poblacion real**, ni **la estructura del candidato**. Los ratios por",
        "blanco son deterministicos dado el pool de control (`construye_respuesta` no usa",
        "`rng`), asi que `sd_s` mide **solo el muestreo de la emision**, que es lo que se",
        "quiere medir.",
        "",
        "**No transductivo, verificado en el codigo del arnes y no asumido.**",
        "`build_fold_adatas` llama `predict_fn(control_pool, panel + [CONTROL_LABEL],",
        "rng_pred)` (`suelo_arnes.py:176`), y `control_pool` sale de",
        "`gene_local == CONTROL_LABEL` (147-148). La funcion no recibe la identidad de la",
        "linea ni una sola celula perturbada. `var_names` entra por clausura porque la firma",
        "del arnes no lo lleva y sin el no se puede localizar la columna del gen blanco: son",
        "los nombres de las columnas del dataset, no el lado real perturbado. Queda",
        "declarado. **El arnes no se modifico** -- `predict_fn` es un parametro suyo.",
        "",
        "**Ninguna cifra del `b` ni del anchor entra como entrada**: son artefactos de",
        "puntuacion. Nada importado: ni Replogle, ni el benchmark, ni tabla de pesos.",
        "",
    ]
    rows = []
    for ln in LINEAS:
        hr = {idx[(ln, j)]["huella_real"] for j in range(N_SEM)}
        hp = {idx[(ln, j)]["huella_pred"] for j in range(N_SEM)}
        hu = {uidx[(ln, j)]["huella_real"] for j in range(N_SEM)
              if (ln, j) in uidx}
        hh = {idx[(ln, j)]["hash_run"]["coincide"] for j in range(N_SEM)}
        con = [idx[(ln, j)]["contrato_emision"] for j in range(N_SEM)]
        rows.append([
            ln, f"**{len(hr)}**", f"**{len(hp)}**",
            ("**si**" if (len(hu) == 1 and hr == hu) else "**NO**"),
            ("**si**" if hh == {True} else f"**NO** {hh}"),
            ("**si**" if all(c["todo_entero"] for c in con) else "**NO**"),
            ("**no**" if not any(c["hay_negativos"] for c in con) else "**SI**"),
            f"{min(c['nnz'] for c in con):,}-{max(c['nnz'] for c in con):,}",
            f"{max(c['max_valor'] for c in con):.0f}",
        ])
    tabla(L, ["linea", "huellas del REAL", "huellas del PREDICHO",
              "huella del real = dispersion_u.md",
              "config_hash de competencia en las 8", "todo entero",
              "hay negativos", "nnz", "valor max"], rows)
    L += [
        "Una sola huella del real por linea y ocho del predicho en las seis. El",
        "`config_hash` recalculado coincide con `a5519017556705ab...` en las **48**.",
        "",
        "**El contrato de emision se comprueba en cada corrida, no se asume.** Conteos",
        "enteros muestreados: se parte de celulas de control remuestreadas con reemplazo y",
        "se reescala gen a gen con `r_h <= 1: Binomial(x, r)` y `r_h > 1: x + Poisson(x(r-1))`.",
        "Las dos ramas tienen esperanza exactamente `x*r` y las dos devuelven enteros. No",
        "hay ningun perfil esperado suave, que es lo que `paso1_contrato.md` prohibe por el",
        "tope `HARD_CAP_NNZ`.",
        "",
        "Ninguna corrida fallo: 48 de 48 `run` y 48 de 48 `score --real-bundle`.",
        "",
        "## 3. El agregado y las 36 celdas",
        "",
    ]
    rows = []
    for ln in LINEAS:
        a = [idx[(ln, j)]["scored"]["filas"]["avg_score"]["from_replicate"]
             for j in range(N_SEM)]
        au = [uidx[(ln, j)]["scored"]["filas"]["avg_score"]["from_replicate"]
              for j in range(N_SEM) if (ln, j) in uidx]
        rows.append([
            ln, " / ".join(f"{x:+.3f}" for x in a),
            f"**{np.mean(a):+.4f}**", f(float(np.std(a, ddof=1))),
            f(max(a) - min(a)),
            f"{np.mean(au):+.4f}" if au else "n/d",
            f"{np.mean(a) - np.mean(au):+.4f}" if au else "n/d",
        ])
    tabla(L, ["linea", "los 8 avg_score", "media", "sd (ddof=1)", "rango",
              "el suelo (dispersion_u.md)", "candidato - suelo"], rows)

    # ---- celdas
    celdas = []
    for ln in LINEAS:
        for met, corto in METRICAS:
            ss = [idx[(ln, j)]["scored"]["filas"][met]["from_replicate"]
                  for j in range(N_SEM)]
            us = [idx[(ln, j)]["u"][met] for j in range(N_SEM)]
            ssu = [uidx[(ln, j)]["scored"]["filas"][met]["from_replicate"]
                   for j in range(N_SEM) if (ln, j) in uidx]
            b0 = nat[ln]["b"].get(met)
            rr = nat[ln]["r"].get(met)
            r = (rr.get("replicate") if isinstance(rr, dict) else rr)
            sp = nat[ln]["splits"].get(met)
            sd_r = float(np.std([s["value"] for s in sp], ddof=0)) if sp else None
            se_r = (sd_r / np.sqrt(5)) if sd_r is not None else None
            bs = [bidx[(ln, j)]["b"][met] for j in range(N_SEM)
                  if (ln, j) in bidx]
            sd_b = float(np.std(bs, ddof=1)) if len(bs) > 1 else None
            span = (r - b0) if (r is not None and b0 is not None) else None
            s_med = float(np.mean(ss))
            sd_s = float(np.std(ss, ddof=1))
            t_r = (abs(s_med) * se_r / abs(span)
                   if (se_r is not None and span) else None)
            t_b = (abs(s_med - 1.0) * sd_b / abs(span)
                   if (sd_b is not None and span) else None)
            margen = (float(np.sqrt(sd_s ** 2 + t_r ** 2 + t_b ** 2))
                      if (t_r is not None and t_b is not None) else None)
            celdas.append({
                "linea": ln, "metrica": corto, "canon": met,
                "s": ss, "u": us, "s_suelo": ssu,
                "s_med": s_med, "sd_s": sd_s, "sd_s0": float(np.std(ss, ddof=0)),
                "rango": max(ss) - min(ss), "min": min(ss), "max": max(ss),
                "u_med": float(np.mean(us)), "sd_u": float(np.std(us, ddof=1)),
                "b": b0, "r": r, "span": span, "se_r": se_r, "sd_b": sd_b,
                "t_r": t_r, "t_b": t_b, "margen": margen,
                "s_suelo_med": float(np.mean(ssu)) if ssu else None,
            })

    L += [
        "`s` es la columna `from_replicate` de `scored.csv`, que con un bundle de",
        "competencia es la que toma `avg_score`.",
        "",
    ]
    rows = []
    for c in celdas:
        rows.append([
            c["linea"], f"`{c['metrica']}`", SENTIDO[c["canon"]],
            " / ".join(f"{x:+.3f}" for x in c["s"]),
            f"**{c['s_med']:+.4f}**", f(c["sd_s"]), f(c["sd_s0"]),
            f(c["rango"]), f"{c['min']:+.4f}", f"{c['max']:+.4f}",
            (f"{c['s_suelo_med']:+.4f}" if c["s_suelo_med"] is not None else "n/d"),
        ])
    tabla(L, ["linea", "metrica", "sentido", "los 8 valores de s", "media",
              "sd_s (ddof=1)", "sd_s (ddof=0)", "rango", "min", "max",
              "el suelo"], rows)

    # ---- 4. supera al b?
    L += [
        "## 4. 🔴 ¿Supera al `b`? Celda por celda, con el margen ABSOLUTO de los tres terminos",
        "",
        "Superar al `b` es, por definicion de la escala, `s > 0`: `s = (u - b)/(r - b)` vale",
        "exactamente 0 cuando `u = b`. Esta es la pregunta **absoluta**, no una comparacion",
        "entre dos candidatos, asi que lleva el margen completo de §3.19 con sus tres",
        "terminos -- incluido el de `b`, que no se cancela aqui:",
        "",
        "```",
        "margen = sqrt( sd_s^2 + (|s| * SE_r/|r-b|)^2 + (|s-1| * sd_b/|r-b|)^2 )",
        "```",
        "",
        "`sd_s` es la **del candidato**, medida con sus 8 semillas, no la del suelo (punto 1",
        "del encargo). `SE_r = sd_r/sqrt(5)` de `anclas_nativas_v2.md`; `sd_b` de las 8",
        "semillas de `dispersion_b.md`. La columna de veredicto exige `s > margen`.",
        "",
    ]
    rows = []
    gana, pierde, empata = [], [], []
    for c in celdas:
        m = c["margen"]
        s = c["s_med"]
        if m is None:
            ver = "n/d"
        elif s > m:
            ver = "🔴 **SUPERA**"
            gana.append(c)
        elif s < -m:
            ver = "**por debajo**"
            pierde.append(c)
        else:
            ver = "indistinguible de `b`"
            empata.append(c)
        rows.append([
            c["linea"], f"`{c['metrica']}`", f"{s:+.4f}", f(c["sd_s"]),
            f(c["t_r"]), f(c["t_b"]), f"**{f(m)}**",
            (f"{s / m:+.2f}" if m else "n/d"), ver,
        ])
    tabla(L, ["linea", "metrica", "s medio", "termino de u (sd_s)",
              "termino de r", "termino de b", "margen absoluto", "s / margen",
              "veredicto"], rows)
    L += [
        f"**Supera al `b` por mas que su margen en {len(gana)} de 36 celdas**; queda por",
        f"debajo en {len(pierde)}; indistinguible en {len(empata)}.",
        "",
    ]
    if gana:
        L += ["Las que superan:", ""]
        rows = []
        for c in sorted(gana, key=lambda c: -(c["s_med"] / c["margen"])):
            rows.append([
                c["linea"], f"`{c['metrica']}`", f"{c['s_med']:+.4f}",
                f(c["margen"]), f"**{c['s_med'] / c['margen']:+.2f}**",
                f(c["u_med"], 4), f(c["b"], 4), f(c["r"], 4),
                (f"{c['s_suelo_med']:+.4f}" if c["s_suelo_med"] is not None
                 else "n/d"),
            ])
        tabla(L, ["linea", "metrica", "s medio", "margen", "s / margen",
                  "u medio", "b", "r", "el suelo en la misma celda"], rows)
        L += [
            "🔴 **La comprobacion que hace falta para que esa tabla no engañe: ¿superan",
            "tambien al SUELO?** El suelo es una copia remuestreada del control, sin ninguna",
            "estructura. Si el candidato no le gana a eso en una celda, lo que esa celda",
            "premia no es su estructura. Se compara contra el mismo `sd_s` del candidato,",
            "porque las dos emisiones se puntuaron contra el mismo bundle y ahi `r` y `b` se",
            "cancelan (§7):",
            "",
        ]
        rows = []
        for c in sorted(gana, key=lambda c: -(c["s_med"] / c["margen"])):
            d = ((c["s_med"] - c["s_suelo_med"])
                 if c["s_suelo_med"] is not None else None)
            res = c["sd_s"] * np.sqrt(2.0 / N_SEM)
            rows.append([
                c["linea"], f"`{c['metrica']}`", f"{c['s_med']:+.4f}",
                (f"{c['s_suelo_med']:+.4f}" if c["s_suelo_med"] is not None
                 else "n/d"),
                (f"**{d:+.4f}**" if d is not None else "n/d"),
                f(res), (f"{d / res:+.2f}" if (d is not None and res) else "n/d"),
                f(c["b"], 4),
                ("🔴 **b ~ 0**" if (c["b"] is not None and abs(c["b"]) < 0.01)
                 else ""),
                ("gana al suelo" if (d is not None and d > res)
                 else ("**PIERDE con el suelo**" if (d is not None and d < -res)
                       else "indistinguible del suelo")),
            ])
        tabla(L, ["linea", "metrica", "s del candidato", "s del suelo",
                  "diferencia", "resolucion k=8", "dif / resolucion", "b",
                  "el b es ~0", "veredicto contra el suelo"], rows)
        nfid = sum(1 for c in gana if c["metrica"] == "fid")
        nb0 = sum(1 for c in gana if c["b"] is not None and abs(c["b"]) < 0.01)
        L += [
            f"**De las {len(gana)} celdas que superan al `b`, {nfid} son `fid` y {nb0} tienen",
            "`b` practicamente cero.** Y una, HT29/`jac`, **queda por debajo del suelo**: el",
            "suelo cobra +0.4803 ahi y el candidato +0.2988. Esa celda no premia estructura,",
            "premia la convencion `J(vacio, vacio) = 1` -- y `dispersion_b.md` §3.2 midio que",
            "la mitad de ese +0.4803 la pone el `b` del bundle por haber caido en el minimo",
            "de sus nueve realizaciones. El candidato, al predecir genes, PIERDE parte de ese",
            "regalo.",
            "",
            "**El cierre del punto 3, contado.** Celdas que (a) superan al `b` por mas que su",
            "margen, (b) superan tambien al suelo por mas que la resolucion, y (c) lo hacen",
            "contra un `b` que no es ~0: **"
            + str(len([c for c in gana
                       if c["s_suelo_med"] is not None
                       and (c["s_med"] - c["s_suelo_med"])
                       > c["sd_s"] * np.sqrt(2.0 / N_SEM)
                       and c["b"] is not None and abs(c["b"]) >= 0.01]))
            + " de 36**.",
            "",
            "🔴 **Y lo que esas celdas NO demuestran.** §1 midio que la correlacion con el",
            "delta real es cero y la concordancia de signo es 0.5. Una `s` positiva en `fid`",
            "o `pds_cosine` con signo al azar no es habilidad especifica del blanco: es que",
            "esas dos metricas premian propiedades que una emision con dispersion y con",
            "genes significativos satisface sin acertar CUALES. `fid` es",
            "`k/max(n_pred, n_conf)`: el `b` tiene `n_pred = 0` en A549 y HAP1 y cobra 0",
            "exacto, asi que cualquier emision que prediga algo con una fraccion de acierto",
            "no nula lo supera. **Superar a un comparador que vale 0 por no predecir nada no",
            "es evidencia de habilidad**, y es exactamente la lectura que `baseline.py:6-10`",
            "advierte en el otro sentido. La refutacion de §1 manda sobre esta tabla.",
            "",
        ]

    # ---- 5. el clamp de mse
    L += [
        "## 5. ⏸️ ¿Cruzo el clamp de `mse`?",
        "",
        "**No. En ninguna de las seis lineas.** `u` sigue muy por encima de `b` (sentido",
        "bajo = peor), asi que la celda sigue devolviendo exactamente 0 y sigue ciega.",
        "",
    ]
    rows = []
    for ln in LINEAS:
        c = [x for x in celdas if x["linea"] == ln and x["canon"] == MSE][0]
        cu = [x for x in uidx.values() if x["linea"] == ln]
        u_suelo = float(np.mean([x["u"][MSE] for x in cu])) if cu else None
        rows.append([
            ln, " / ".join(f"{x:.3f}" for x in c["u"]), f(c["u_med"]),
            f(c["sd_u"]), f(c["b"]), f(c["r"]),
            ("**no**" if c["u_med"] > c["b"] else "**SI**"),
            f(u_suelo) if u_suelo else "n/d",
            (f"{c['u_med'] - u_suelo:+.4f}" if u_suelo else "n/d"),
            f"{c['s_med']:+.4f}",
        ])
    tabla(L, ["linea", "los 8 u crudos", "u medio", "sd_u", "b", "r",
              "cruzo (u < b)", "u del suelo", "candidato - suelo",
              "s"], rows)
    peor = [x for x in celdas if x["canon"] == MSE]
    L += [
        "La columna `candidato - suelo` dice algo que conviene no pasar por alto: el",
        "candidato esta **peor** que el suelo en `mse` en las seis lineas (`mse` es de",
        "sentido bajo). Mover genes sin acertar cuales empeora el error cuadratico del",
        "pseudobulk respecto a no moverlos. La celda no lo registra -- el clamp devuelve 0",
        "en los dos casos -- pero el `u` crudo si, y es la misma ceguera que",
        "`dispersion_u.md` §5.1 declaro: dos emisiones distintas, las dos por debajo de `b`,",
        "son indistinguibles.",
        "",
        "## 6. 🔴 Las dos cohortes moviles, y la tasa a la que mueve `E`",
        "",
        "Esto cierra el hueco de §3.17. El suelo predecia genes significativos en 0-3",
        "perturbaciones por corrida, asi que su `E` se movia por accidente, una unidad cada",
        "varias semillas. Este candidato predice mucho mas, y aqui se ve a que tasa.",
        "",
    ]
    rows = []
    for ln in LINEAS:
        E = [idx[(ln, j)]["cohortes"]["E_union_vacia"] for j in range(N_SEM)]
        un = [idx[(ln, j)]["cohortes"]["n_union_no_vacia"] for j in range(N_SEM)]
        fid = [idx[(ln, j)]["cohortes"]["n_fid"] for j in range(N_SEM)]
        dep = [idx[(ln, j)]["cohortes"]["suma_n_DE_pred"] for j in range(N_SEM)]
        npr = [idx[(ln, j)]["cohortes"]["n_perts_con_pred_sig"]
               for j in range(N_SEM)]
        inter = [idx[(ln, j)]["cohortes"]["suma_n_DE_interseccion"]
                 for j in range(N_SEM)]
        der = [idx[(ln, j)]["cohortes"]["suma_n_DE_real"] for j in range(N_SEM)]
        Eu = [uidx[(ln, j)]["cohortes"]["E_union_vacia"] for j in range(N_SEM)
              if (ln, j) in uidx]
        dpu = [uidx[(ln, j)]["cohortes"]["suma_n_DE_pred"] for j in range(N_SEM)
               if (ln, j) in uidx]
        rows.append([
            ln, "/".join(map(str, E)), len(set(E)),
            "/".join(map(str, fid)), len(set(fid)),
            f"{min(npr)}-{max(npr)}", f"{min(dep)}-{max(dep)}",
            f"{min(der)}-{max(der)}", f"{min(inter)}-{max(inter)}",
            "/".join(map(str, Eu)) if Eu else "n/d",
            f"{min(dpu)}-{max(dpu)}" if dpu else "n/d",
        ])
    tabla(L, ["linea", "E del candidato por semilla", "valores distintos de E",
              "cohorte de fid por semilla", "valores distintos de fid",
              "perts con algun gen predicho sig", "suma n_DE_pred",
              "suma n_DE_real", "suma n_DE_interseccion",
              "E del SUELO (dispersion_u.md)", "n_DE_pred del suelo"], rows)
    L += [
        "**La tasa, dicha con los numeros delante.** El suelo sumaba 0-3 genes predichos",
        "significativos en todo el panel; este candidato suma decenas. Y el efecto sobre",
        "`E` es inmediato: donde el suelo tenia `E` casi fija, el candidato la baja y la",
        "mueve. La interseccion con el lado real es la columna que dice si eso compra algo.",
        "",
        "🔴 **A549 es el caso que lo paga entero.** Su span en `jac` es 0.0187, menos que una",
        "unidad de convencion (`1/33 = 0.0303`), asi que cada perturbacion que pierde su",
        "`1.0` gratis cuesta 1.62 puntos de `s`. El candidato, prediciendo genes, pierde",
        f"varias: su `jac` sale {[c for c in celdas if c['linea'] == 'A549' and c['canon'] == JAC][0]['s_med']:+.4f} "
        f"con `sd_s` {f([c for c in celdas if c['linea'] == 'A549' and c['canon'] == JAC][0]['sd_s'], 4)}, "
        "y arrastra el `avg_score` de la linea a",
        f"{np.mean([idx[('A549', j)]['scored']['filas']['avg_score']['from_replicate'] for j in range(N_SEM)]):+.4f}. "
        "Es el mecanismo que `dispersion_u.md` §5.2 aislo con un solo gen",
        "(STAT5A, 1 de 8 semillas); aqui se ve a escala. **Predecir señal en una celda con",
        "span pequeno es caro con independencia de si la señal acierta.**",
        "",
        "## 7. La resolucion del candidato",
        "",
        "Para comparar DOS candidatos en la misma celda del mismo bundle, `r` y `b` se",
        "cancelan (`s1 - s2 = (u1 - u2)/(r - b)`) y solo entra `sd_s`. Con `k` replicas:",
        "",
        "```",
        "resolucion_k = sd_s * sqrt(2/k)",
        "```",
        "",
        "Con la `sd_s` **del candidato**, que es lo que el punto 1 pedia y que no se asume",
        "igual a la del suelo -- y no lo es:",
        "",
    ]
    rows = []
    for c in celdas:
        su = [x for x in celdas if False]
        sd_suelo = None
        ssu = c["s_suelo"]
        if len(ssu) > 1:
            sd_suelo = float(np.std(ssu, ddof=1))
        fila = [c["linea"], f"`{c['metrica']}`", f(c["sd_s"]),
                f(sd_suelo) if sd_suelo is not None else "n/d",
                (f"{c['sd_s'] / sd_suelo:.2f}x"
                 if (sd_suelo and sd_suelo > 0) else "n/d")]
        for k in KS:
            fila.append(f(c["sd_s"] * np.sqrt(2.0 / k)))
        fila.append(f"**{f(c['margen'])}**")
        rows.append(fila)
    tabla(L, ["linea", "metrica", "sd_s del CANDIDATO", "sd_s del suelo",
              "razon"] + [f"resolucion k={k}" for k in KS]
          + ["margen absoluto de §4"], rows)
    razones = [(c["sd_s"] / float(np.std(c["s_suelo"], ddof=1)))
               for c in celdas
               if len(c["s_suelo"]) > 1 and np.std(c["s_suelo"], ddof=1) > 0]
    L += [
        f"**La `sd_s` del candidato NO es la del suelo**: la razon va de "
        f"{min(razones):.2f}x a {max(razones):.2f}x en las {len(razones)} celdas",
        "comparables. Asumir la del suelo habria dado un margen equivocado en los dos",
        "sentidos segun la celda. Es la razon por la que el punto 1 pedia medirla.",
        "",
        "## 8. La `sd` acumulada k=2..8, y donde es Bernoulli",
        "",
        "En el orden `j = 0..7` de la familia de semillas, fijado antes de medir.",
        "",
    ]
    rows = []
    bern = []
    for c in celdas:
        acc = [float(np.std(c["s"][:k], ddof=1)) for k in range(2, N_SEM + 1)]
        arr = acc[0] == 0.0 and max(acc) > 0
        salto = (max(acc) / min(x for x in acc if x > 0)
                 if any(x > 0 for x in acc) else None)
        if arr:
            bern.append(c)
        rows.append([c["linea"], f"`{c['metrica']}`"]
                    + [("0" if a == 0 else f(a)) for a in acc]
                    + ["🔴 **si**" if arr else "no"])
    tabla(L, ["linea", "metrica"] + [f"k={k}" for k in range(2, N_SEM + 1)]
          + ["arranca en 0 y salta"], rows)
    if bern:
        L += [
            f"**{len(bern)} celdas arrancan en 0 y saltan**: "
            + ", ".join(f"{c['linea']}/`{c['metrica']}`" for c in bern) + ".",
            "",
            "🔴 **Y en A549/`jac` la `sd` no describe una nube ni una Bernoulli: describe una",
            "RETICULA, y se puede verificar con aritmetica exacta.** Sus 8 valores de `s` son",
            "",
        ]
        cj = [c for c in celdas if c["linea"] == "A549" and c["canon"] == JAC]
        if cj:
            c = cj[0]
            paso = (1.0 / 33.0) / c["span"] if c["span"] else None
            rows = []
            for j in range(N_SEM):
                E = idx[("A549", j)]["cohortes"]["E_union_vacia"]
                rows.append([
                    j, idx[("A549", j)]["seed_emision"], E,
                    f"{15 - E}", f(c['u'][j], 5), f"{c['u'][j] * 33:.0f}/33",
                    f"{c['s'][j]:+.4f}",
                    (f"{-(15 - E) * paso:+.4f}" if paso else "n/d"),
                    (f"{abs(c['s'][j] + (15 - E) * paso):.1e}" if paso else "n/d"),
                ])
            tabla(L, ["j", "semilla", "E del candidato",
                      "unidades perdidas vs el suelo (E=15)", "u de jac",
                      "u en 33avos", "s medido",
                      "-(unidades) x (1/33)/span", "\\|dif\\|"], rows)
            L += [
                f"`(1/33)/span = (1/33)/{f(c['span'], 6)} = {f(paso, 4)}`, y los 8 `s` son "
                "multiplos exactos de esa cantidad con residuo < 1e-15. La celda **solo puede",
                "tomar valores en esa reticula**: no hay nada continuo que promediar. Pedir",
                "mas semillas ahi no estrecha un intervalo, estima la distribucion de cuantas",
                "unidades de convencion pierde la emision. Es la misma aritmetica que",
                "`dispersion_u.md` §5.2 verifico con un solo gen; con un candidato que predice",
                "señal se ve en las 8 semillas a la vez.",
                "",
            ]
    else:
        L += [
            "**Ninguna celda arranca en 0 y salta.** A diferencia del suelo, que tenia 5",
            "(`dispersion_u.md` §4), este candidato mueve genes en las 8 semillas de las seis",
            "lineas, asi que no hay celdas donde el evento ocurra o no ocurra: la dispersion",
            "es continua. Lo que SI queda discreto es el tamano del salto en `jac` de A549,",
            "pero ahi el evento ocurre en las 8, no en una.",
            "",
        ]
    ceros = [c for c in celdas if c["sd_s"] == 0.0]
    L += [
        f"**{len(ceros)} celdas con `sd_s = 0`**, todas de `mse_capped_norm`: el clamp (§5).",
        "No son estabilidad, son ceguera; el `u` crudo se mueve (`sd_u` 0.02-0.07) y el `s`",
        "no.",
        "",
    ]

    # ---- 20 filas trabajadas
    L += [
        "## 9. Las 20 filas trabajadas",
        "",
        "**Criterio de seleccion, declarado antes de mirar los valores**: las celdas de",
        "`pds_cosine`, `fid` y `jac` en las seis lineas = 18 filas, mas `lfc_nmae` en A549",
        "y MCF7, los dos extremos del `span` de esa metrica = **20**. Se eligen `pds`, `fid`",
        "y `jac` porque son las tres donde §4 encuentra celdas que superan al `b`, que es la",
        "pregunta del encargo; no por su valor.",
        "",
        "Cada fila trae los 8 `u` crudos, los 8 `s`, la aritmetica `(u-b)/(r-b)` rehecha a",
        "mano sobre la media, y los tres terminos del margen.",
        "",
    ]
    sel = [c for c in celdas if c["metrica"] in ("pds_cosine", "fid", "jac")]
    sel += [c for c in celdas if c["metrica"] == "lfc_nmae"
            and c["linea"] in ("A549", "MCF7")]
    rows = []
    for c in sel:
        reh = ((c["u_med"] - c["b"]) / c["span"]) if c["span"] else None
        rows.append([
            c["linea"], f"`{c['metrica']}`",
            " / ".join(f"{x:.4f}" for x in c["u"]),
            f(c["u_med"]), f(c["sd_u"]), f(c["b"]), f(c["r"]), f(c["span"]),
            (f"{reh:+.4f}" if reh is not None else "n/d"),
            f"{c['s_med']:+.4f}",
            (f"{abs(reh - c['s_med']):.1e}" if reh is not None else "n/d"),
            f(c["sd_s"]), f(c["t_r"]), f(c["t_b"]), f"**{f(c['margen'])}**",
            f(c["sd_s"] * np.sqrt(2.0 / 8)),
        ])
    tabla(L, ["linea", "metrica", "los 8 u crudos", "u medio", "sd_u", "b", "r",
              "span r-b", "s a mano desde u medio", "s medio medido",
              "\\|dif\\|", "sd_s", "termino de r", "termino de b",
              "margen absoluto", "resolucion k=8"], rows)
    L += [
        f"Son **{len(rows)}** filas.",
        "",
        "## Lo que no se hizo",
        "",
        "- **No se leyo nada del lado real perturbado para construir la prediccion.** El",
        "  arnes lo impide estructuralmente y se verifico en su codigo (§2). No hay fallo",
        "  del arnes que reportar.",
        "- **No se uso ninguna cifra del `b` ni del anchor como entrada.** Entran solo en la",
        "  puntuacion, que es su sitio.",
        "- **No se ajusto el candidato mirando su `s`.** Las dos magnitudes descartadas",
        "  (§1.1) se descartaron por propiedades medidas de la EMISION y nunca se puntuaron.",
        "  La version que corre se midio una vez, con 8 semillas, y no se toco despues.",
        "- **No se itero.** Una iteracion seria un candidato nuevo con su propia medicion.",
        "- **No se proponen umbrales ni se decide adopcion.** §4 aplica la regla ya escrita",
        "  en §3.16 y §3.18; no se añade ninguna.",
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
