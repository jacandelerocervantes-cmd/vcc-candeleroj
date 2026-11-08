#!/usr/bin/env python
"""Escribe modelo/out/sensibilidad_entrega.md de lo que haya en disco.

Degrada: lo que no se corrio sale como "no corrida", no inventa nada.
"""
import glob
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, "/home/ubuntu/vcc/modelo/src/noche")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import comun as C  # noqa: E402
import banco as B  # noqa: E402
import sem as S  # noqa: E402
import sup as U  # noqa: E402
import a_bio as A  # noqa: E402

OUT_MD = "/home/ubuntu/vcc/modelo/out/sensibilidad_entrega.md"
OUT_JSON = "/home/ubuntu/vcc/modelo/out/sensibilidad_entrega.json"
MET = U.MET
CORTO = {"pds_cosine": "pds", "expr_mse_unbiased_capped_norm": "mse",
         "de_wilcoxon_lfc_nmae": "lfc", "de_wilcoxon_direction_fidelity_yield_raw": "fid",
         "de_wilcoxon_direction_reach_raw": "reach", "de_wilcoxon_sig_jaccard": "jac"}


def f4(x):
    return C.f4a(x) if x is not None else "n/d"


def fs(x):
    return C.f4(x) if x is not None else "n/d"


def carga():
    d = {"celdas": {}}
    for f in sorted(glob.glob(os.path.join(B.CELDAS, "*.json"))):
        r = json.load(open(f))
        d["celdas"][os.path.basename(f)[:-5]] = r
    for k, p in (("de", A.OUT_DE), ("rep", A.OUT_REP), ("res", A.OUT_RES),
                 ("sem", os.path.join(U.WORK, "semillas.json"))):
        d[k] = json.load(open(p)) if os.path.exists(p) else None
    return d


def u_de(r, brazo, m):
    b = (r.get("brazos") or {}).get(brazo) or {}
    return (b.get("u") or {}).get(m)


def s_de(r, brazo, m):
    b = (r.get("brazos") or {}).get(brazo) or {}
    s = b.get("s") or {}
    fila = s.get(m) or {}
    return fila.get("from_replicate")


def s_avg(r, brazo):
    b = (r.get("brazos") or {}).get(brazo) or {}
    s = b.get("s") or {}
    return (s.get("avg_score") or {}).get("from_replicate")


# ──────────────────────────────── secciones ────────────────────────────────
def cab(W, d):
    W += ["# El mapa de sensibilidad de la superficie de entrega", "",
          "## 🔴 ENCABEZADO OBLIGATORIO", "",
          "El **brazo B** y el **brazo C** miden **LA SENSIBILIDAD DEL MARCADOR** sobre "
          "**los controles de la entrega**, no la biologia de D/E/F, **que no existe**: "
          "de D/E/F solo hay 18,400 celulas de control por contexto y ninguna "
          "perturbada. La estructura inyectada viene del **brazo A** y es una "
          "**TRANSFERENCIA DE ESTRUCTURA declarada, dada como rango** — tres numeros "
          "de forma (amplitud, fraccion de genes movidos, fraccion comun), nunca "
          "valores de gen ni de perturbacion. **Ninguna cifra de aqui es una "
          "prediccion del score real.**", "",
          "El **brazo A** se mide en Jiang (6 contextos, somera) y VCC 2025 (1 contexto, "
          "profunda), que **no son** la superficie de entrega.", "",
          "Dos referencias externas que fijan la escala, y que no se recalculan aqui:", "",
          "* el **nulo** de `ntc_nulo.md`: a `n = 400` y en estos mismos contextos, "
          "partiendo controles contra controles, el conjunto ED sale **vacio** y las "
          "metricas de DE no registran nada. Cero medido, no supuesto.",
          "* el **azar** de `pds_cosine`: **0.5 exacto**, derivado del codigo en "
          "`pds_objetivo.md` §1.2 (`1 - ((n-1)/2)/(n-1)`), para todo `n`.", ""]
    return W


def sec0(W, d):
    W += ["---", "", "## 0. La construccion, declarada antes de medir", ""]
    W += ["### 0.1 Las tres superficies y que mide cada brazo", ""]
    W += C.tabla(["brazo", "superficie", "n contextos", "que mide", "que NO mide"], [
        ["A", "Jiang IFNB + VCC 2025", "6 + 1", "que estructura tiene el delta REAL",
         "nada de la entrega"],
        ["B", "controles de D/E/F (y A/B/C aparte)", "3 + 3",
         "que registra el marcador con verdad SINTETICA de rejilla declarada",
         "biologia; la rejilla la elegi yo"],
        ["C", "controles de D/E/F", "3",
         "lo mismo con la estructura MEDIDA en A, el techo y el mapa calidad->puntuacion",
         "no transfiere valores, solo forma"]])
    W += [""]
    if d["sem"]:
        sm = d["sem"]
        W += ["### 0.2 Semillas: familia NUEVA, comprobadas EN CODIGO y guardadas enteras", "",
              f"Familia: `{sm['familia']}`. **{sm['n']} semillas**, guardadas enteras en "
              f"`{os.path.join(U.WORK, 'semillas.json')}` y en el json de este reporte "
              f"(`lista_plana`), no solo su conteo. Minimo {min(sm['lista_plana'])}, "
              f"maximo {max(sm['lista_plana'])}.", ""]
        W += C.tabla(["comprobacion", "resultado"], [
            ["choques internos entre familias", len(sm["choques_internos"]) or "**0**"],
            ["dentro de los rangos vetados " + str(sm["rangos_vetados"]),
             len(sm["choques_veto"]) or "**0**"],
            ["contra la enumeracion de `noche/comun.py`", len(sm["choques_comun"]) or "**0**"],
            ["contra PASO 5a (`seleccion/semillas.py`)", len(sm["choques_5a"]) or "**0**"],
            ["contra el suelo de llamada (`sl_semillas.py`)",
             len(sm["choques_suelo_llamada"]) or "**0**"],
            [f"contra el barrido de {sm['n_barrido']} enteros 202[6-8]xxxx en "
             "`modelo/src/**/*.py` y los json de `modelo/out/`",
             len(sm["choques_barrido"]) or "**0**"],
            ["**sin choques**", "**si**" if sm["sin_choques"] else "**NO**"]])
        W += [""]
    cel = d["celdas"]
    una = next(iter(cel.values()), None)
    W += ["### 0.3 La superficie de entrega como banco de sensibilidad", ""]
    if una:
        inv, par = una.get("inventario", {}), una.get("particion", {})
        W += C.tabla(["que", "valor"], [
            ["contextos con verdad en disco", "**ninguno** (D/E/F son solo control)"],
            ["celulas de control por contexto", "18,400 (46 `ntc_id` x 400 exactas)"],
            ["genes", f"{inv.get('shape', ['', ''])[1]:,}" if inv.get("shape") else "n/d"],
            ["nnz del contexto leido", f"{inv.get('nnz', 0):,}"],
            ["X entero", "si" if inv.get("X_entero") else "NO"],
            ["particion POR GUIA (nunca por celula)",
             f"{par.get('n_control', 0):,} celulas de control ({len(par.get('guias_control', []))} guias) "
             f"/ {par.get('n_pp', 0)} pseudo-perturbaciones de 400"],
            ["guia partida", "**no**" if par.get("guia_partida") is False else "n/d"],
            ["regimen de la llamada", "400 contra 9,200 — el brazo `n400_ref_completa` de "
             "`ntc_nulo.md`, el declarado ahi como el que se parece a la entrega"],
            ["etiquetas", "23 nombres de gen del propio reto, sorteados con semilla "
             f"`{par.get('diag_etiquetas', {}).get('semilla_etiquetas')}` entre las "
             f"{par.get('diag_etiquetas', {}).get('n_candidatas_que_son_gen')} del panel "
             "(el marcador exige que resuelvan a gen: con etiquetas inventadas "
             "`baseline` se niega)"],
            ["genes filtrables (CPM control > 5, menos los 23 de etiqueta)",
             f"{inv.get('n_filtrables', 0):,}"],
            ["lado real", f"{una.get('real', {}).get('shape')} filas x genes"],
            ["lado pred", f"{(una.get('brazos', {}).get('suelo', {}).get('pred', {}) or {}).get('shape')}"]])
        W += [""]
    W += ["### 0.4 El dispositivo de inyeccion y la rejilla del brazo B", "",
          "La inyeccion es `cand01_covar._reescala_entero` **sin modificar**: Binomial "
          "si la razon es <= 1, Poisson aditivo si es > 1; la esperanza es `x * r` en "
          "las dos ramas y la salida sigue siendo cuentas enteras. Es el mismo "
          "dispositivo con el que `pds_objetivo.md` §3.1 emitio sus oraculos.", "",
          "🔴 La rejilla del brazo B esta **declarada antes de correr** y es mia — por "
          "eso es un peso inventado, y por eso existe el brazo C:", ""]
    W += C.tabla(["eje", "valores", "que es"], [
        ["amplitud", ", ".join(str(a) for a in B.AMPLITUDES), "`|log2FC|` inyectado"],
        ["fraccion de genes", ", ".join(str(a) for a in B.FRACCIONES),
         "de los filtrables, por pseudo-perturbacion"],
        ["fraccion comun", ", ".join(str(a) for a in B.COMUNES),
         "del movimiento, comun a las 23 (mismo gen, mismo signo)"],
        ["celdas", f"{len(B.REJILLA_B)} (producto)",
         f"centro declarado: celda {B.CENTRO_B} = {B.REJILLA_B[B.CENTRO_B]}"]])
    W += ["", "Y el barrido de calidad del brazo C: `f` en "
          f"{{{', '.join(str(f) for f in B.F_GRID)}}}, en dos variantes declaradas "
          "distintas (**magnitud**: direccion correcta y magnitud encogida; "
          "**identidad**: `floor(f*n)` genes correctos y el resto en genes equivocados "
          "con la misma amplitud). En las dos el eje comun se da por acertado, para que "
          "`f` module solo la componente especifica.", ""]
    return W


def sec_A(W, d):
    W += ["---", "", "## 1. BRAZO A — la biologia: que estructura tiene el delta real", "",
          "🔴 Jiang y VCC 2025. **Ninguna cifra de esta seccion es de la entrega.**", ""]
    de, res = d["de"], d["res"]
    if not de:
        W += ["_(brazo A sin correr: no hay `a_de.json`)_", ""]
        return W
    W += ["### 1.1 Verificacion de las poblaciones y de la puerta", ""]
    fil = []
    for L in C.LINES:
        if L in de.get("jiang", {}):
            r = de["jiang"][L]
            fil.append([f"Jiang {L}", "si" if r["huella_coincide"] else "**NO**",
                        f"{r['perts_con_|ED|_igual_a_p2']}/{r['n_perts']}",
                        C.dur(r["dt_s"])])
    if "v2025" in de:
        r = de["v2025"]
        fil.append(["VCC 2025", "(poblacion de `vcc2025_suelo`, ya verificada por sha256)",
                    f"{r['perts_con_|ED|_igual_a_p2']}/{r['n_perts']}", C.dur(r["dt_s"])])
    W += C.tabla(["superficie", "huella sha256 = la guardada",
                  "perts con |ED| identico al PASO 5a", "duracion"], fil)
    W += [""]
    if res:
        W += ["### 1.2 🔴 A1/A3 como RANGO, y el nucleo-5 aparte", "",
              f"El nucleo son **{', '.join(A.NUCLEO)}**. Criterio: {res['criterio_nucleo']}.", ""]
        enc = ["grupo", "n perts", "|ED| mediana", "|ED| [min, max]", "frac ED mediana",
               "|log2FC| mediana de los sig.", "|log2FC| p90 de los sig."]
        fil = []
        for k, v in res["grupos"].items():
            def g(n, c="mediana"):
                return (v.get(n) or {}).get(c)
            fil.append([f"`{k}`", v["n_perts"],
                        f"{g('n_ed'):.0f}" if g("n_ed") is not None else "n/d",
                        f"[{g('n_ed', 'min'):.0f}, {g('n_ed', 'max'):.0f}]"
                        if g("n_ed", "min") is not None else "n/d",
                        C.fe(g("frac_ed"), 2), f4(g("lfc_mediana_sig")), f4(g("lfc_p90_sig"))])
        W += C.tabla(enc, fil)
        W += ["", "**Fraccion de genes probados que se mueve por encima de cada corte de "
              "`|log2FC|` (mediana del grupo):**", ""]
        enc = ["grupo"] + [f"≥ {c}" for c in A.CORTES]
        fil = []
        for k, v in res["grupos"].items():
            fil.append([f"`{k}`"] + [C.fe((v.get(f"frac_sobre_{c}") or {}).get("mediana"), 2)
                                     for c in A.CORTES])
        W += C.tabla(enc, fil)
        W += [""]
    W += ["### 1.3 🔴 A2 — el reparto COMPARTIDO / ESPECIFICO, ANCLADO en su control", ""]
    rep = d["rep"]
    if not rep:
        W += ["_(sin correr)_", ""]
    else:
        W += ["**(a) ENTRE CONTEXTOS, mismo blanco** — se LEE de `magnitud_centro.md` §4, "
              "no se recalcula. Solo existe en Jiang: 2025 es un contexto y ahi esta "
              "cantidad no se puede medir.", "",
              "| cantidad | agregada (6 lineas) | [min, max] sobre los 15 pares |",
              "|---|---|---|",
              "| delta CRUDO, fraccion ANCLADA | **+0.2596** | [+0.1445, +0.4519] |",
              "| delta ESPECIFICO del blanco, fraccion ANCLADA | **+0.3080** | [+0.1842, +0.5053] |",
              "", "(la observada sin anclar es +0.4212 / +0.4233 y **no es** la fraccion "
              "de centro: lleva dentro el suelo de no-comparticion, que el control "
              "permutado mide en +0.2182 / +0.1666 contra un suelo analitico de 1/6 = "
              "0.1667.)", "",
              "**(b) ENTRE PERTURBACIONES, dentro de un contexto** — la que el brazo C "
              "puede inyectar, medida aqui, con su control de signos ±1 por "
              "perturbacion (el mismo en las dos mitades):", ""]
        fil = []
        for L, v in (rep.get("entre_perts_jiang") or {}).items():
            fil.append([f"Jiang {L}", v["n_perts"], f"{v['n_genes_espacio']:,}",
                        fs(v["frac_obs_mediana"]),
                        fs(np.median([f["control_medio"] for f in v["por_particion"]])),
                        f4(v["por_particion"][0]["suelo_analitico_1_sobre_n"]),
                        f"**{fs(v['frac_anclada_mediana'])}**",
                        f"[{fs(v['frac_anclada_rango'][0])}, {fs(v['frac_anclada_rango'][1])}]"])
        v = rep.get("entre_perts_2025")
        if v:
            fil.append(["VCC 2025", v["n_perts"], f"{v['n_genes_espacio']:,}",
                        fs(np.median([f["frac_obs"] for f in v["por_particion"]])),
                        fs(np.median([f["control_medio"] for f in v["por_particion"]])),
                        f4(v["por_particion"][0]["suelo_analitico_1_sobre_n"]),
                        f"**{fs(v['frac_anclada_mediana'])}**",
                        f"[{fs(v['frac_anclada_rango'][0])}, {fs(v['frac_anclada_rango'][1])}]"])
        W += C.tabla(["superficie", "n perts", "genes del espacio", "frac. observada",
                      "control de signos", "suelo 1/n", "frac. ANCLADA (mediana sobre j)",
                      "[min, max] sobre j"], fil)
        W += [""]
    return W


def sec_A_inter(W, d):
    de = d["de"]
    if not de:
        return W
    W += ["### 1.4 🔴 Los intermedios crudos: perturbacion por perturbacion", "",
          "Jiang, las 33 filas del panel. `|ED|` y, entre parentesis, la mediana de "
          "`|log2FC|` de sus genes significativos. El nucleo-5 en **negrita**.", ""]
    enc = ["#", "pert"] + [L for L in C.LINES if L in de.get("jiang", {})]
    fil = []
    panel = sorted((de["jiang"][next(iter(de["jiang"]))]["por_pert"]).keys()) \
        if de.get("jiang") else []
    for i, g in enumerate(panel, 1):
        row = [i, f"**`{g}`**" if g in A.NUCLEO else f"`{g}`"]
        for L in C.LINES:
            if L not in de.get("jiang", {}):
                continue
            p = de["jiang"][L]["por_pert"].get(g, {})
            med = (p.get("lfc_abs_pct_sig") or {}).get("50")
            row.append(f"{p.get('n_ed', 0)}" + (f" ({med:.2f})" if med else ""))
        fil.append(row)
    W += C.tabla(enc, fil)
    W += [""]
    if "v2025" in de:
        pp = de["v2025"]["por_pert"]
        ks = sorted(pp, key=lambda k: -(pp[k].get("n_ed") or 0))
        sel = [ks[int(round(x))] for x in np.linspace(0, len(ks) - 1, 20)]
        W += ["**VCC 2025**: 20 de sus 100 perturbaciones, equiespaciadas por `|ED|` "
              "(de la mayor a la menor), con sus intermedios enteros.", ""]
        fil = []
        for i, g in enumerate(sel, 1):
            p = pp[g]
            q = p.get("lfc_abs_pct_sig") or {}
            fil.append([i, f"`{g}`", p.get("n_probados"), p.get("n_ed"),
                        C.fe(p.get("frac_ed"), 2), f4(q.get("50")), f4(q.get("90")),
                        f4(p.get("lfc_abs_max_sig")),
                        f4(p.get("signo_neg_sig"))])
        W += C.tabla(["#", "pert", "n probados", "|ED|", "frac ED", "|log2FC| p50 sig",
                      "p90", "max", "frac. de signo negativo"], fil)
        W += [""]
    return W


def filas_celda(r, brazos, met):
    out = []
    for b in brazos:
        if b not in (r.get("brazos") or {}):
            continue
        out.append([b] + [f4(u_de(r, b, m)) for m in met] + [fs(s_avg(r, b))])
    return out


def sec_B(W, d):
    W += ["---", "", "## 2. BRAZO B — que registra el marcador a n=400", "",
          "🔴 Verdad SINTETICA, rejilla mia. Cada celda tiene su propio lado real, su "
          "propio `baseline` y su propio `prep-real-bundle` (el bundle no se puede "
          "compartir entre celdas porque el lado real cambia), y dos brazos de "
          "prediccion contra el mismo bundle: **suelo** (sin estructura) y **oraculo** "
          "(la misma estructura, celulas distintas, ruido nuevo).", ""]
    cel = {k: v for k, v in d["celdas"].items() if "_B" in k}
    if not cel:
        W += ["_(sin celdas del brazo B en disco)_", ""]
        return W
    W += ["### 2.1 Las celdas corridas", ""]
    fil = []
    for k, r in sorted(cel.items()):
        p = r["params"]
        e = r.get("estructura", {})
        fil.append([f"`{k}`", r["ctx"], p["amplitud"], p["frac_genes"], p["frac_comun"],
                    f"{e.get('m_movidos', 0):,}", f"{e.get('n_comun', 0):,}",
                    f"{e.get('n_especifico', 0):,}",
                    "si" if (r.get("hash_bundle") or {}).get("coincide") else
                    ("NO" if r.get("hash_bundle") else "n/d"),
                    C.dur(r.get("dt_celda_s", 0))])
    W += C.tabla(["celda", "ctx", "amplitud", "frac genes", "frac comun", "genes movidos",
                  "comunes", "especificos", "config_hash de competencia", "duracion"], fil)
    W += ["", "### 2.2 🔴 La superficie: las seis metricas, `u` y `s`, celda por celda y "
          "brazo por brazo", "",
          "`u` es la fila `mean` de `agg_results.csv`; `s` es `from_replicate` de "
          "`scored.csv` contra el bundle de esa misma celda. El **recorrido** es "
          "`oraculo − suelo`: lo que el marcador puede registrar para esa estructura.", ""]
    for k, r in sorted(cel.items()):
        p = r["params"]
        W += [f"**`{k}`** — ctx {r['ctx']}, amplitud {p['amplitud']}, "
              f"frac genes {p['frac_genes']}, frac comun {p['frac_comun']} "
              f"({r.get('estructura', {}).get('m_movidos', 0):,} genes movidos)", ""]
        enc = ["cantidad"] + [CORTO[m] for m in MET] + ["avg"]
        fil = []
        for b in ("suelo", "oraculo"):
            if b not in (r.get("brazos") or {}):
                continue
            fil.append([f"`u` {b}"] + [f4(u_de(r, b, m)) for m in MET] + [""])
        for b in ("suelo", "oraculo"):
            if b not in (r.get("brazos") or {}):
                continue
            fil.append([f"`s` {b}"] + [fs(s_de(r, b, m)) for m in MET] + [fs(s_avg(r, b))])
        if "suelo" in (r.get("brazos") or {}) and "oraculo" in (r.get("brazos") or {}):
            fil.append(["**recorrido `u`**"] + [
                fs((u_de(r, "oraculo", m) - u_de(r, "suelo", m))
                   if (u_de(r, "oraculo", m) is not None and u_de(r, "suelo", m) is not None)
                   else None) for m in MET] + [""])
            fil.append(["**recorrido `s`**"] + [
                fs((s_de(r, "oraculo", m) - s_de(r, "suelo", m))
                   if (s_de(r, "oraculo", m) is not None and s_de(r, "suelo", m) is not None)
                   else None) for m in MET] + [
                fs((s_avg(r, "oraculo") - s_avg(r, "suelo"))
                   if (s_avg(r, "oraculo") is not None and s_avg(r, "suelo") is not None)
                   else None)])
        W += C.tabla(enc, fil)
        W += [""]
    W += ["### 2.3 Como se mueve cada metrica al mover UN eje", "",
          "Del subconjunto en cruz: el centro y, desde el, un eje a la vez. Cada fila es "
          "el `u` del ORACULO (lo mejor alcanzable para esa estructura).", ""]
    enc = ["celda", "amplitud", "frac genes", "frac comun"] + [CORTO[m] for m in MET] + ["s avg"]
    fil = []
    for k, r in sorted(cel.items()):
        p = r["params"]
        fil.append([f"`{k}`", p["amplitud"], p["frac_genes"], p["frac_comun"]]
                   + [f4(u_de(r, "oraculo", m)) for m in MET] + [fs(s_avg(r, "oraculo"))])
    W += C.tabla(enc, fil)
    W += ["", "La misma tabla para el SUELO (sin estructura): es el control de que lo que "
          "se mueve arriba es la estructura y no la celda.", ""]
    fil = []
    for k, r in sorted(cel.items()):
        p = r["params"]
        fil.append([f"`{k}`", p["amplitud"], p["frac_genes"], p["frac_comun"]]
                   + [f4(u_de(r, "suelo", m)) for m in MET] + [fs(s_avg(r, "suelo"))])
    W += C.tabla(enc, fil)
    W += [""]
    return W


def sec_C(W, d):
    W += ["---", "", "## 3. BRAZO C — la union, de donde saldrian los pesos", ""]
    res = d["res"]
    cel = {k: v for k, v in d["celdas"].items() if "_C" in k}
    if res and res.get("puntos_C"):
        W += ["### 3.1 Los tres puntos de estructura MEDIDA, y de donde salen", "",
              res["puente_declarado"], ""]
        W += C.tabla(["punto", "amplitud |log2FC|", "frac de genes movidos", "frac comun"],
                     [[p["nombre"], p["amplitud"], p["frac_genes"], p["frac_comun"]]
                      for p in res["puntos_C"]])
        W += [""]
    if not cel:
        W += ["_(sin celdas del brazo C en disco)_", ""]
        return W
    W += ["### 3.2 🔴 C2 — EL TECHO DE LA SUPERFICIE: que puntua un predictor PERFECTO", "",
          "El brazo `mag_f1` es el oraculo exacto: la misma estructura del lado real, "
          "aplicada a celulas distintas del pool de control, con ruido Binomial/Poisson "
          "nuevo. **No tiene por que dar 1.0**: el `n = 400` y la profundidad limitan lo "
          "detectable.", ""]
    enc = ["celda", "ctx", "punto"] + [CORTO[m] for m in MET] + ["avg"]
    fil = []
    for k, r in sorted(cel.items()):
        nm = (r.get("estructura", {}).get("punto_C") or {}).get("nombre", "")
        fil.append([f"`{k}`", r["ctx"], nm] + [f4(u_de(r, "mag_f1", m)) for m in MET] + [""])
        fil.append([f"`{k}` (`s`)", r["ctx"], nm]
                   + [fs(s_de(r, "mag_f1", m)) for m in MET] + [fs(s_avg(r, "mag_f1"))])
    W += C.tabla(enc, fil)
    W += ["", "### 3.3 🔴 C3 — EL MAPA CALIDAD -> PUNTUACION, por metrica", "",
          "`s` (`from_replicate`) de un modelo que acierta una fraccion `f` de la "
          "componente especifica. Dos variantes declaradas distintas.", ""]
    for k, r in sorted(cel.items()):
        W += [f"**`{k}`** — ctx {r['ctx']}, "
              f"estructura {r.get('estructura', {}).get('punto_C')}", ""]
        for modo, pref in (("magnitud", "mag"), ("identidad", "ide")):
            fil = []
            for f in [None] + B.F_GRID:
                b = "suelo" if f is None else f"{pref}_f{f:g}"
                if b not in (r.get("brazos") or {}):
                    continue
                fil.append([("suelo" if f is None else f"f = {f:g}")]
                           + [fs(s_de(r, b, m)) for m in MET] + [fs(s_avg(r, b))]
                           + [f4(u_de(r, b, "pds_cosine"))])
            if fil:
                W += [f"variante **{modo}**", ""]
                W += C.tabla(["calidad"] + [CORTO[m] for m in MET] + ["avg", "`u` pds"], fil)
                W += [""]
    W += ["### 3.4 🔴 C4 — a que parte de la estructura responde cada metrica", "",
          "La tabla, no la eleccion. Por metrica: el recorrido que le abre el brazo B al "
          "mover cada eje, y la pendiente del mapa de calidad del brazo C. **No se eligen "
          "pesos aqui.**", ""]
    celB = {k: v for k, v in d["celdas"].items() if "_B" in k}
    enc = ["metrica", "recorrido `u` (oraculo − suelo) en el centro de la rejilla",
           "pendiente `s` de f=0 a f=1 (magnitud)", "pendiente `s` de f=0 a f=1 (identidad)",
           "`u` del suelo", "techo `u` (f=1)"]
    fil = []
    centro = celB.get(f"D_B{B.CENTRO_B:02d}")
    c0 = next(iter(sorted(cel.values())), None) if cel else None
    for m in MET:
        rec = None
        if centro:
            a, b_ = u_de(centro, "oraculo", m), u_de(centro, "suelo", m)
            rec = None if (a is None or b_ is None) else a - b_
        def pend(pref):
            if not c0:
                return None
            x0, x1 = s_de(c0, f"{pref}_f0", m), s_de(c0, f"{pref}_f1", m)
            return None if (x0 is None or x1 is None) else x1 - x0
        fil.append([f"`{m}`", fs(rec), fs(pend("mag")), fs(pend("ide")),
                    f4(u_de(c0, "suelo", m)) if c0 else "n/d",
                    f4(u_de(c0, "mag_f1", m)) if c0 else "n/d"])
    W += C.tabla(enc, fil)
    W += [""]
    return W


def sec_def(W, d):
    W += ["---", "", "## 4. Defectos de la construccion, como datos sobre la medicion", "",
          "1. **El banco tiene 23 pseudo-perturbaciones, la entrega tiene 300.** Con "
          "18,400 celulas de control por contexto, repartirlas en pool de control y "
          "pseudo-perturbadas es un trueque: 23/23 da el regimen 400 contra 9,200 y 23 "
          "competidores para `pds_cosine` (granularidad 1/22 = 0.0455), no 300 (1/299).",
          "2. **El lado pred se remuestrea del MISMO pool de control que el lado real usa "
          "de referencia.** Es el contrato del arnes (`predict_floor`), pero significa "
          "que las celulas del pred son filas del control del real, con reemplazo.",
          "3. **La rejilla del brazo B la elegi yo.** Es el peso inventado que el encargo "
          "nombra; el brazo C existe para sustituirla por estructura medida, y las dos "
          "cosas no se mezclan en ninguna tabla.",
          "4. **La transferencia del brazo C es de FORMA, y la forma se mide en otras "
          "superficies.** Amplitud, fraccion de genes y fraccion comun vienen de Jiang y "
          "2025; que esos tres numeros sean los de D/E/F no se puede comprobar.",
          "5. **La amplitud inyectada es la MISMA en todos los genes movidos.** El delta "
          "real tiene una distribucion de amplitudes (seccion 1.1 la da): el banco inyecta "
          "su mediana, no su forma.",
          "6. **El eje comun del brazo C se da por acertado** en las dos variantes de `f`. "
          "Si un modelo no acierta el eje comun, su `s` es menor que el de este mapa.",
          "7. **El gen propio de cada etiqueta no se mueve** y los 23 genes-etiqueta salen "
          "del sorteo de genes movidos: la exclusion del panel los quita del espacio. No "
          "se modela eficiencia de knockdown.",
          "8. **Un `s` de aqui no es un `s` de la entrega.** El bundle es el de esta "
          "superficie sintetica, con sus 23 etiquetas; el `b` y el `r` de la entrega real "
          "no existen porque no hay verdad.",
          "9. **El ruido de la inyeccion es Binomial/Poisson sobre una celula de control.** "
          "`sobredispersion.md` ya midio que ese dispositivo no reproduce la "
          "sobredispersion del delta real; lo que el banco mide es sensibilidad a un "
          "desplazamiento de medias, no a un cambio de dispersion.", ""]
    W += ["---", "", "## 5. Las entradas", ""]
    W += C.tabla(["que", "de donde"], [
        ["controles de la entrega", "`/datos/reto/final/context_{D,E,F}.h5ad`, "
         "`/datos/reto/val/context_{A,B,C}.h5ad` (solo lectura)"],
        ["etiquetas", "`pert_counts.csv` de la particion de cada contexto"],
        ["Jiang", "`comun.arma_linea` (= `suelo_arnes.build_fold_adatas`), "
         "`seed_real = 20260001 + fold`"],
        ["VCC 2025", "`v25.selecciona` (20280201 / 20280202)"],
        ["deltas de Jiang para A2b", "`/datos/modelo/out/magnitud_centro/deltas_*.npz`"],
        ["reparto entre contextos", "`modelo/out/magnitud_centro/{calc,ctrl}.json`"],
        ["|ED| de referencia", "`/datos/modelo/out/seleccion/p2_datos.json`"],
        ["codigo", "`modelo/src/sensib/`"],
        ["celdas crudas", f"`{B.CELDAS}/*.json`"]])
    W += ["", f"_Generado {time.strftime('%Y-%m-%d %H:%M:%S')} UTC. "
          f"Celdas en disco: {len(d['celdas'])}._"]
    return W


def main():
    d = carga()
    W = []
    for f in (cab, sec0, sec_A, sec_A_inter, sec_B, sec_C, sec_def):
        try:
            W = f(W, d)
        except Exception as e:
            W += ["", f"_(seccion `{f.__name__}` fallo al escribirse: {e!r})_", ""]
    C.escribe_md(OUT_MD, W)
    resumen = {"celdas": sorted(d["celdas"]), "n_celdas": len(d["celdas"]),
               "semillas": d["sem"], "a_resumen": d["res"],
               "u_s_por_celda": {k: {b: {"u": (v.get("brazos", {}).get(b, {}) or {}).get("u"),
                                         "s": (v.get("brazos", {}).get(b, {}) or {}).get("s")}
                                     for b in (v.get("brazos") or {})}
                                 for k, v in d["celdas"].items()}}
    C.escribe_json(OUT_JSON, resumen)
    print(f"escrito {OUT_MD}: {len(W)} lineas; {len(d['celdas'])} celdas")


if __name__ == "__main__":
    main()
