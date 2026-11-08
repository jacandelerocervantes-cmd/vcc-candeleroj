#!/usr/bin/env python
"""Escribe modelo/out/suelo_llamada.md a partir de suelo_llamada.json."""
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sl_analisis as A  # noqa: E402

OUT_MD = "/home/ubuntu/vcc/modelo/out/suelo_llamada.md"
LINES = A.LINES


def tabla(enc, filas):
    o = ["| " + " | ".join(map(str, enc)) + " |", "|" + "|".join(["---"] * len(enc)) + "|"]
    for f in filas:
        o.append("| " + " | ".join("" if c is None else str(c) for c in f) + " |")
    return o


def f4(x, s=False):
    if x is None:
        return "n/d"
    return f"{x:+.4f}" if s else f"{x:.4f}"


def rs(r, s=False):
    """mediana [p10, p90] [min, max] (n)"""
    if not r or r.get("n", 0) == 0:
        return ["n/d"] * 4
    return [f"**{f4(r['mediana'], s)}**", f"[{f4(r['p10'], s)}, {f4(r['p90'], s)}]",
            f"[{f4(r['min'], s)}, {f4(r['max'], s)}]", r["n"]]


def main():
    r = A.main()
    W = ["# El suelo de reproducibilidad de la llamada, dentro de linea", "",
         "## 🔴 CÓMO SE LEE", "",
         "* Esto mide **Jiang** (6 lineas, IFNB, 33 perturbaciones) y **VCC 2025** (`test`, un contexto). De aqui "
         "solo viaja el **MECANISMO**; ninguna cifra es de D/E/F.",
         "* **MCF7 aporta 12,867 de los 24,410 genes de R(g,L) en Jiang** y domina cualquier agregado: todo se da "
         "**por linea, por par y por perturbacion**, con rango. Ninguna media sola.",
         "* Los dos brazos llaman con **n/2 celulas perturbadas y n/2 de control, disjuntas**. La unica diferencia "
         "entre DENTRO y ENTRE es si la segunda mitad es de la misma linea o de otra.",
         "* Formato: datos. Sin veredicto, sin conclusiones.", "",
         "## 0. Construccion, declarada antes de medir", ""]
    W += tabla(["# ", "cantidad", "valor"], [
        [1, "poblacion Jiang", "`comun.arma_linea`, `seed_real = 20260001 + indice`, la de `venn_ed`; 400 celulas por perturbacion y TODO el pool de control"],
        [2, "poblacion 2025", "`v25.selecciona` (20280201 / 20280202), la de `vcc2025_suelo`: min(400, disponibles) por pert, 18,400 controles"],
        [3, "particion", "POR CELULA; una rng por (superficie, j): primero el control, luego cada perturbacion en orden de panel; cada grupo se permuta y se parte en dos"],
        [4, "🔴 control", "PARTIDO en dos mitades disjuntas: la mitad A de perturbadas se llama contra cA, la B contra cB"],
        [5, "desempate de impares", "ALTERNADO: la celula sobrante va a A si j es par y a B si j es impar (control y perturbaciones truncadas de 2025)"],
        [6, "particiones", "8 por superficie, j = 0..7"],
        [7, "la llamada", "`cell-eval2 run -ar A -ap B --preset vcc2026 --profile vcc2026 --set control_source=pred --write-degenes`, luego `prepare_de` + `p_adj < 0.05` + `_exclude_own_gene` (lo de venn_ed)"],
        [8, "🔴 por que `control_source=pred`", "el preset puro tiene `control_source = real`: el lado pred se llamaria contra el control del REAL (`run.py`, el combinado de referencia), compartiendo su ruido. Con `pred`, cada mitad usa su medio control. Se asserta en `run_params.yaml` de cada corrida"],
        [9, "verificacion", "en j = 0 de cada superficie, corrida con los papeles cambiados: cada mitad por el lado REAL (preset intacto en ese lado) tiene que dar el mismo conjunto, gen a gen y con el mismo p_adj, que por el lado pred (§0.2)"],
        [10, "brazo DENTRO", "J(R(L,A_j), R(L,B_j)); f = |X \\ Y|/|X| en los dos sentidos, media de los definidos"],
        [11, "brazo ENTRE", "J(R(L1,A_j), R(L2,B_j)) y J(R(L1,B_j), R(L2,A_j)), media de los definidos; f(L1→L2) igual"],
        [12, "cabecera", "dJ = ½(J_DENTRO(L1) + J_DENTRO(L2)) − J_ENTRE(L1,L2), por (par, pert, j); por (par, pert): media sobre j y rango sobre j; df(L1→L2) = f_ENTRE(L1→L2) − f_DENTRO(L1)"],
        [13, "union vacia", "J indefinida; convencion (excl) la excluye — la de cabecera —, (marc) la pone en 1.0 como el marcador"],
        [14, "igualada", "recorte de venn_ed §5: por (pert, j), m = menor tamano NO VACIO de los 12 conjuntos-mitad de Jiang; cada uno a sus m de p_adj mas baja; empates por clave aleatoria; vacios siguen vacios. 2025: m = menor no vacio de sus 2 mitades"],
        [15, "tamaño a n entero", "Jiang: `venn_ed.json` (CLI); 2025: `p2_datos.json` del paso 5a (`compute_de`, el camino que reprodujo el CLI 198/198 en Jiang)"],
    ])
    sem = r["semillas"]
    W += ["", "### 0.1 Semillas: familia NUEVA, comprobadas EN CODIGO y guardadas enteras", "",
          "🔴 **fold se numera desde 1** (A549 = 1 … MCF7 = 6, VCC 2025 = 7). El rango vetado por el encargo, "
          "«2027000x–2028200x», llega hasta **20282009**; con fold = 0 la formula `20282000 + 100*fold + j` daria "
          "20282000–20282007, **dentro** del veto. Con fold desde 1 se cumplen la formula y el veto a la vez.", ""]
    W += tabla(["uso", "semillas enteras"], [[k, ", ".join(map(str, v))] for k, v in sem["semillas_enteras"].items()])
    W += ["", f"Rangos vetados comprobados: {sem['rangos_vetados']}. Choques: internos {len(sem['choques_internos'])}, "
          f"con el veto {len(sem['choques_veto'])}, con `noche/comun.py` {len(sem['choques_comun'])}, con PASO 5a "
          f"{len(sem['choques_5a'])}, con el barrido ({sem['n_barrido']} enteros del codigo y de los JSON parseados; maximo "
          f"previo {sem['max_barrido']}) {len(sem['choques_barrido'])}. **Sin choques: {'si' if sem['sin_choques'] else 'NO'}.** "
          "Guardadas en `/datos/modelo/out/suelo_llamada/semillas.json` y en `suelo_llamada.json`.", ""]
    W += ["### 0.2 Verificacion de poblaciones y del camino de la llamada", ""]
    filas = []
    for t in LINES + ["v2025"]:
        v = r["verif"][t]; m = r["meta"][t]
        hu = m.get("coincide_venn", m.get("coincide"))
        filas.append([t, "si" if hu else "**NO**", m["n_control"], "si" if v["A_real_igual_A_pred"] else "**NO**",
                      "si" if v["B_real_igual_B_pred"] else "**NO**", v["genes_distintos_A"], v["genes_distintos_B"]])
    W += tabla(["superficie", "huella sha256 = la guardada", "n control", "mitad A: lado real = lado pred (j=0)",
                "mitad B: lado real = lado pred (j=0)", "genes distintos A", "genes distintos B"], filas)
    W += ["", "**Cortes del control (cA / cB) por j**", ""]
    W += tabla(["superficie"] + [f"j={j}" for j in range(8)],
               [[t] + [f"{a}/{b}" for a, b in (r["cortes"][t][str(j)] if str(j) in r["cortes"][t] else r["cortes"][t][j] for j in range(8))]
                for t in LINES + ["v2025"]])
    imp = r["cortes_pert_impar_2025"]
    W += ["", f"Perturbaciones de 2025 con n impar ({len(imp)}), cortes A/B por j: " +
          "; ".join(f"`{g}` " + " ".join(f"{a}/{b}" for a, b in v) for g, v in list(imp.items())[:15]) +
          (" …" if len(imp) > 15 else "") + ". Las de Jiang son todas 200/200.", ""]

    for var, etq in (("jiang_crudo", "CRUDA"), ("jiang_igualado", "IGUALADA EN POTENCIA")):
        J = r[var]
        W += ["---", "", f"## 1. Solape de los conjuntos ED, los dos brazos — Jiang, {etq}", "",
              "### 1.1 Brazo DENTRO, por linea (33 perts × 8 particiones = 264 celdas)", ""]
        W += tabla(["linea", "J DENTRO mediana", "[p10, p90]", "[min, max]", "n definidas", "# union vacia",
                    "J (marc) mediana", "f DENTRO mediana", "f [min, max]"],
                   [[L] + rs(J["dentro_linea"][L]["Jd"]) + [J["dentro_linea"][L]["n_celdas_union_vacia"],
                     f4(J["dentro_linea"][L]["Jd_marc"].get("mediana")), f4(J["dentro_linea"][L]["fd"].get("mediana")),
                     f"[{f4(J['dentro_linea'][L]['fd'].get('min'))}, {f4(J['dentro_linea'][L]['fd'].get('max'))}]"]
                    for L in LINES])
        W += ["", "### 1.2 Brazo ENTRE, por par (264 celdas por par)", ""]
        filas = []
        for k, v in J["entre_par"].items():
            a, b = k.split("–")
            filas.append([k] + rs(v["Je"]) + [f4(J["fe_opar"][f"{a}→{b}"].get("mediana")),
                                              f4(J["fe_opar"][f"{b}→{a}"].get("mediana"))])
        W += tabla(["par", "J ENTRE mediana", "[p10, p90]", "[min, max]", "n definidas", "f ENTRE (1º→2º) mediana",
                    "f ENTRE (2º→1º) mediana"], filas)
        W += ["", f"## 2. 🔴 LA CIFRA DE CABECERA: J DENTRO − J ENTRE — Jiang, {etq}", "",
              "Cada valor de la distribucion es una celda (par, pert) con dJ = media sobre las 8 particiones. "
              "≈ 0: lo especifico de contexto es ruido de llamada. Grande: es contexto. (El texto lo dice el encargo; aqui solo los numeros.)", "",
              f"Celdas (par, pert) definidas: {J['n_celdas_par_g_definidas']} de {J['n_celdas_par_g']}.", "",
              "### 2.1 Por par", ""]
        W += tabla(["par", "dJ mediana", "[p10, p90]", "[min, max]", "# perts", "mediana del ancho del rango sobre j",
                    "dJ (marc) mediana", "dJ (marc) [min, max]"],
                   [[k] + rs(v, True) + [f4(v.get("ancho_rango_j_mediana")), f4(J["dJm_par"][k].get("mediana"), True),
                                         f"[{f4(J['dJm_par'][k].get('min'), True)}, {f4(J['dJm_par'][k].get('max'), True)}]"]
                    for k, v in J["dJ_par"].items()])
        W += ["", "### 2.2 Por linea (sus 5 pares × 33 perts)", ""]
        W += tabla(["linea", "dJ mediana", "[p10, p90]", "[min, max]", "# celdas", "df (L→otra) mediana", "df [min, max]"],
                   [[L] + rs(J["dJ_linea"][L], True) + [f4(J["df_linea"][L].get("mediana"), True),
                     f"[{f4(J['df_linea'][L].get('min'), True)}, {f4(J['df_linea'][L].get('max'), True)}]"] for L in LINES])
        W += ["", "### 2.3 Por perturbacion (sobre los 15 pares) y J DENTRO de cada linea (media sobre j)", ""]
        W += tabla(["# ", "pert", "dJ mediana", "dJ [min, max]", "# pares"] + [f"Jd {L}" for L in LINES],
                   [[i + 1, f"`{g}`", f4(v.get("mediana"), True), f"[{f4(v.get('min'), True)}, {f4(v.get('max'), True)}]",
                     v.get("n", 0)] + [f4(v["Jd_por_linea"][L]) for L in LINES]
                    for i, (g, v) in enumerate(J["dJ_pert"].items())])
        W += ["", "### 2.4 df = f ENTRE(L1→L2) − f DENTRO(L1), por par ordenado", ""]
        W += tabla(["par", "df mediana", "[p10, p90]", "[min, max]", "# perts"],
                   [[k] + rs(v, True) for k, v in J["df_opar"].items()])
        W += ["", f"🔴 **Intermedios ({etq}): 20 celdas (par, pert)** equiespaciadas por dJ", ""]
        W += tabla(["# ", "par", "pert", "Jd L1", "Jd L2", "Je", "dJ media j", "dJ min j", "dJ max j", "# j", "|R(L1,A)| j=0", "|R(L2,B)| j=0"],
                   [[i + 1, x["par"], f"`{x['pert']}`", f4(x["Jd_L1"]), f4(x["Jd_L2"]), f4(x["Je"]), f4(x["dJ_media_j"], True),
                     f4(x["dJ_min_j"], True), f4(x["dJ_max_j"], True), x["n_j"], x["tam_L1_A_j0"], x["tam_L2_B_j0"]]
                    for i, x in enumerate(J["intermedios"])])
        if J.get("intermedios_j"):
            ij = J["intermedios_j"]
            W += ["", f"**Intermedios por particion de la celda mediana** ({ij['par']}, `{ij['pert']}`)", ""]
            W += tabla(["j", "|L1A|", "|L1B|", "|L2A|", "|L2B|", "Jd L1", "Jd L2", "Je", "dJ"],
                       [[x["j"], x["|L1A|"], x["|L1B|"], x["|L2A|"], x["|L2B|"], f4(x["Jd_L1"]), f4(x["Jd_L2"]), f4(x["Je"]),
                         f4(x["dJ"], True)] for x in ij["filas"]])
        W += [""]
    W += ["---", "", "## 3. El recorte igualado: m por perturbacion (las 8 particiones)", ""]
    W += tabla(["# ", "pert", "m por j", "m min", "m max"],
               [[i + 1, f"`{g}`", " ".join(map(str, v)), min(v), max(v)] for i, (g, v) in enumerate(r["m_recorte_jiang"].items())])
    W += ["", "---", "", "## 4. Tamaño de los conjuntos ED a n/2 contra n entero", ""]
    T = r["tamanos"]
    W += tabla(["superficie", "Σ|R| a n entero", "Σ|R| a n/2: media de 16 mitades", "[min, max]", "n/2 ÷ n",
                "control n / mitad (j=0)", "UMI mediana/celula de control"],
               [[t, f"{T[t]['full']:,}", f"{T[t]['half_media']:,.1f}", f"[{T[t]['half_min']:,}, {T[t]['half_max']:,}]",
                 f4(T[t]["ratio_media"]), f"{sum(T[t]['cortes_j0'])} / {T[t]['cortes_j0'][0]}", f"{r['profundidad'][t]:,.0f}"]
                for t in LINES + ["v2025"]])
    W += ["", "🔴 **Intermedios Jiang: |R| a n entero / mediana de las 16 mitades, por perturbacion y linea**", ""]
    W += tabla(["# ", "pert"] + LINES,
               [[i + 1, f"`{g}`"] + [f"{T[L]['por_pert'][g]['full']} / {T[L]['por_pert'][g]['half_mediana']:.0f}" for L in LINES]
                for i, g in enumerate(r["panel"])])
    pp = sorted(r["perts25"], key=lambda g: T["v2025"]["por_pert"][g]["full"])
    sel = [pp[i] for i in np.linspace(0, len(pp) - 1, 20).round().astype(int)]
    W += ["", "🔴 **Intermedios 2025: 20 perturbaciones equiespaciadas por |R| a n entero**", ""]
    W += tabla(["# ", "pert", "|R| n entero", "|R| n/2 mediana", "[min, max] de las 16 mitades", "J DENTRO media j", "J [min, max] sobre j"],
               [[i + 1, f"`{g}`", T["v2025"]["por_pert"][g]["full"], f"{T['v2025']['por_pert'][g]['half_mediana']:.0f}",
                 f"[{T['v2025']['por_pert'][g]['half_min']}, {T['v2025']['por_pert'][g]['half_max']}]",
                 f4(r["v2025_crudo"]["por_pert"][g]["Jd"][0] if r["v2025_crudo"]["por_pert"][g]["Jd"] else None),
                 (f"[{f4(r['v2025_crudo']['por_pert'][g]['Jd'][1])}, {f4(r['v2025_crudo']['por_pert'][g]['Jd'][2])}]"
                  if r["v2025_crudo"]["por_pert"][g]["Jd"] else "n/d")] for i, g in enumerate(sel)])
    W += ["", "---", "", "## 5. 🔴 El brazo DENTRO en VCC 2025, al lado del de Jiang", "",
          "Un solo contexto: no hay brazo ENTRE. Las dos cifras, lado a lado, con la profundidad y el n de cada una.", ""]
    filas = []
    for t in LINES:
        dc, di = r["jiang_crudo"]["dentro_linea"][t], r["jiang_igualado"]["dentro_linea"][t]
        filas.append([f"Jiang {t}", f"{r['profundidad'][t]:,.0f}", T[t]["cortes_j0"][0], 200] + rs(dc["Jd"])[:3] +
                     [f4(dc["fd"].get("mediana")), f4(di["Jd"].get("mediana")), f"{T[t]['half_media']:,.0f}"])
    for var in ("v2025",):
        c, i_ = r["v2025_crudo"], r["v2025_igualado"]
        filas.append(["**VCC 2025**", f"{r['profundidad']['v2025']:,.0f}", T["v2025"]["cortes_j0"][0], "≤ 200"] + rs(c["Jd"])[:3] +
                      [f4(c["fd"].get("mediana")), f4(i_["Jd"].get("mediana")), f"{T['v2025']['half_media']:,.0f}"])
    W += tabla(["superficie", "UMI mediana control", "control / mitad", "perturbadas / mitad", "J DENTRO mediana (cruda)",
                "[p10, p90]", "[min, max]", "f DENTRO mediana", "J DENTRO mediana (igualada)", "Σ|R| por mitad"], filas)
    c5 = r["v2025_crudo"]
    W += ["", f"2025, sobre las {len(r['perts25'])} perturbaciones × 8 particiones: celdas con union vacia "
          f"{c5['n_celdas_union_vacia']}; J (marc) mediana {f4(c5['Jd_marc'].get('mediana'))}; J DENTRO por perturbacion "
          f"(media sobre j): mediana {f4(c5['Jd_pert'].get('mediana'))}, rango [{f4(c5['Jd_pert'].get('min'))}, {f4(c5['Jd_pert'].get('max'))}].", ""]
    W += ["---", "", "## Defectos propios", "",
          "1. **n/2 no es n.** Los dos brazos llaman a n/2 y son comparables entre si, pero ninguno es la llamada a n entero "
          "del paso 5a: el residuo 0.9484 de 5a se midio a n entero y esta cifra no lo corrige directamente (§4 da cuanto se encoge).",
          "2. **El brazo ENTRE en Jiang mezcla dos cosas de la otra linea**: contexto biologico Y otra potencia (otro n de control, "
          "otra profundidad). La version igualada solo recorta el TAMAÑO del conjunto, no iguala n ni profundidad de las celdas.",
          "3. **La convencion de la union vacia mueve J mucho** donde los conjuntos son chicos (A549, HAP1, K562 con perts "
          "de 0–3 genes): por eso van las dos.",
          "4. **dJ promedia dos J DENTRO** (de L1 y de L2). Si una linea reproduce mucho mejor que la otra, el promedio lo esconde; "
          "las tablas 1.1 y 2.3 dan cada una.",
          "5. **El recorte igualado de Jiang usa el menor de 12 conjuntos**: en las perts donde una mitad de una linea "
          "llama 1 gen, todo se recorta a 1 y J solo toma 0 o 1.",
          "6. **2025 a n entero salio de `compute_de`**, no del CLI; el camino es el que reprodujo el CLI 198/198 en Jiang, "
          "pero en 2025 no se cruzo.",
          "7. **Las 8 particiones no son independientes** (las mismas celdas se reparten 8 veces): el rango sobre j "
          "mide la variabilidad de la particion, no un intervalo de confianza.",
          "8. **`control_source=pred` es una modificacion del preset**, aplicada solo al lado pred; la verificacion de §0.2 "
          "la cruza contra el lado real (preset intacto) en j = 0 de cada superficie, no en las 8 particiones.",
          "9. **2025 tiene 15 perturbaciones truncadas** (< 400 celulas), asi que sus mitades tienen < 200.", "",
          f"_Generado {time.strftime('%F %T')} UTC. Tiempo de sondas acumulado: {r['dt_cli_total_s']/60:.1f} min. "
          "Codigo: `modelo/src/suelo_llamada/`; conjuntos: `/datos/modelo/out/suelo_llamada/sets/`._"]
    open(OUT_MD + ".tmp", "w").write("\n".join(W) + "\n")
    os.replace(OUT_MD + ".tmp", OUT_MD)


if __name__ == "__main__":
    main()
