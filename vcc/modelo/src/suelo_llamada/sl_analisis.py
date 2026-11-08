#!/usr/bin/env python
"""Analisis del SUELO DE LA LLAMADA. Lee los conjuntos de `sl_sondas.py`.

MEDIDAS, DECLARADAS ANTES DE MEDIR
  J(X,Y) = |X∩Y|/|X∪Y|; union vacia -> indefinida. Dos convenciones: (excl) se
    excluye la celda; (marc) vale 1.0, la del marcador. La de cabecera es (excl).
  f(X->Y) = |X \\ Y| / |X|; X vacio -> indefinida.
  DENTRO(L,g,j): Jd = J(R(L,A_j), R(L,B_j)); fd = media de f(A->B) y f(B->A) definidas.
  ENTRE(L1,L2,g,j), par no ordenado: Je = media de J(R(L1,A_j),R(L2,B_j)) y
    J(R(L1,B_j),R(L2,A_j)) definidas. Ordenado: fe(L1->L2) = media de
    f(R(L1,A_j)->R(L2,B_j)) y f(R(L1,B_j)->R(L2,A_j)).
  CABECERA: dJ(L1,L2,g,j) = ½(Jd(L1)+Jd(L2)) − Je(L1,L2) (todas definidas);
    df(L1->L2,g,j) = fe(L1->L2) − fd(L1). Por (par, g): media sobre j y rango
    [min, max] sobre j. Luego distribuciones por par, por linea y por perturbacion.
  IGUALADA (recorte de venn_ed §5): por (g, j), m = el menor tamano NO VACIO de
    los 12 conjuntos-mitad (6 lineas x 2 mitades); cada uno se recorta a sus m
    de p_adj mas baja, empates por clave aleatoria (semilla 20282800+j); los
    vacios siguen vacios. En 2025 (un contexto): m = menor no vacio de las 2
    mitades, semilla 20282810+j.
"""
import json
import os
import sys
from itertools import combinations, permutations

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sl_semillas as S  # noqa: E402

SETS = "/datos/modelo/out/suelo_llamada/sets"
OUT_JSON = "/home/ubuntu/vcc/modelo/out/suelo_llamada.json"
LINES = S.LINES
JJ = range(S.J)


def jac(a, b, conv="excl"):
    if not a and not b:
        return None if conv == "excl" else 1.0
    return len(a & b) / len(a | b)


def fno(a, b):
    return None if not a else len(a - b) / len(a)


def mean_def(xs):
    v = [x for x in xs if x is not None]
    return float(np.mean(v)) if v else None


def carga():
    D = {}
    for t in LINES + ["v2025"]:
        D[t] = {}
        for j in JJ:
            p = os.path.join(SETS, f"{t}_j{j}.json")
            if os.path.exists(p):
                D[t][j] = json.load(open(p))
    return D


def recorte(D, tags, perts, seed_base):
    """Devuelve {tag: {j: {"A": {g: set}, "B": {...}}}} recortado, y m por (g,j)."""
    out = {t: {j: {"A": {}, "B": {}} for j in JJ} for t in tags}
    ms = {}
    for j in JJ:
        rng = np.random.default_rng(seed_base + j)
        for g in perts:
            tam = [len(D[t][j][h].get(g, {})) for t in tags for h in "AB"]
            nv = [x for x in tam if x > 0]
            m = min(nv) if nv else 0
            ms[(g, j)] = m
            for t in tags:
                for h in "AB":
                    pv = D[t][j][h].get(g, {})
                    gs = sorted(pv)
                    if not gs or m == 0:
                        out[t][j][h][g] = set(); continue
                    k = rng.random(len(gs))
                    o = sorted(range(len(gs)), key=lambda i: (pv[gs[i]], k[i]))
                    out[t][j][h][g] = {gs[i] for i in o[:m]}
    return out, ms


def crudo(D, tags, perts):
    return {t: {j: {h: {g: set(D[t][j][h].get(g, {})) for g in perts} for h in "AB"} for j in JJ} for t in tags}


def medidas_jiang(R, perts):
    """Devuelve dicts por (L,g,j) y (L1,L2,g,j)."""
    Jd, fd, Jdm = {}, {}, {}
    for L in LINES:
        for g in perts:
            for j in JJ:
                A, B = R[L][j]["A"][g], R[L][j]["B"][g]
                Jd[(L, g, j)] = jac(A, B); Jdm[(L, g, j)] = jac(A, B, "marc")
                fd[(L, g, j)] = mean_def([fno(A, B), fno(B, A)])
    Je, Jem, fe = {}, {}, {}
    for L1, L2 in combinations(LINES, 2):
        for g in perts:
            for j in JJ:
                a1, b1 = R[L1][j]["A"][g], R[L1][j]["B"][g]
                a2, b2 = R[L2][j]["A"][g], R[L2][j]["B"][g]
                Je[(L1, L2, g, j)] = mean_def([jac(a1, b2), jac(b1, a2)])
                Jem[(L1, L2, g, j)] = mean_def([jac(a1, b2, "marc"), jac(b1, a2, "marc")])
                fe[(L1, L2, g, j)] = mean_def([fno(a1, b2), fno(b1, a2)])
                fe[(L2, L1, g, j)] = mean_def([fno(a2, b1), fno(b2, a1)])
    return Jd, Jdm, fd, Je, Jem, fe


def cabecera(Jd, Jdm, fd, Je, Jem, fe, perts):
    dJ, dJm, df = {}, {}, {}
    for L1, L2 in combinations(LINES, 2):
        for g in perts:
            for j in JJ:
                a, b, e = Jd[(L1, g, j)], Jd[(L2, g, j)], Je[(L1, L2, g, j)]
                dJ[(L1, L2, g, j)] = None if None in (a, b, e) else 0.5 * (a + b) - e
                dJm[(L1, L2, g, j)] = 0.5 * (Jdm[(L1, g, j)] + Jdm[(L2, g, j)]) - Jem[(L1, L2, g, j)]
    for L1, L2 in permutations(LINES, 2):
        for g in perts:
            for j in JJ:
                x, y = fe[(L1, L2, g, j)], fd[(L1, g, j)]
                df[(L1, L2, g, j)] = None if None in (x, y) else x - y
    return dJ, dJm, df


def por_celda(d, keys, perts):
    """{(key..., g): (media_j, min_j, max_j, n_j)} sobre j definidos."""
    out = {}
    for k in keys:
        for g in perts:
            v = [d[k + (g, j)] for j in JJ if d.get(k + (g, j)) is not None]
            out[k + (g,)] = (float(np.mean(v)), float(min(v)), float(max(v)), len(v)) if v else None
    return out


def resumen(vals):
    v = [x for x in vals if x is not None]
    if not v:
        return {"n": 0}
    v = np.array(v)
    return {"n": int(v.size), "mediana": float(np.median(v)), "p10": float(np.percentile(v, 10)),
            "p90": float(np.percentile(v, 90)), "min": float(v.min()), "max": float(v.max()),
            "media": float(v.mean())}


def analiza_jiang(R, perts):
    Jd, Jdm, fd, Je, Jem, fe = medidas_jiang(R, perts)
    dJ, dJm, df = cabecera(Jd, Jdm, fd, Je, Jem, fe, perts)
    pares = list(combinations(LINES, 2)); opares = list(permutations(LINES, 2))
    cJd = por_celda(Jd, [(L,) for L in LINES], perts)
    cfd = por_celda(fd, [(L,) for L in LINES], perts)
    cJe = por_celda(Je, pares, perts)
    cfe = por_celda(fe, opares, perts)
    cdJ = por_celda(dJ, pares, perts)
    cdJm = por_celda(dJm, pares, perts)
    cdf = por_celda(df, opares, perts)
    m = lambda c: None if c is None else c[0]
    res = {"dentro_linea": {}, "entre_par": {}, "dJ_par": {}, "dJm_par": {}, "dJ_linea": {}, "dJ_pert": {},
           "df_opar": {}, "df_linea": {}, "fe_opar": {}}
    for L in LINES:
        res["dentro_linea"][L] = {"Jd": resumen([Jd[(L, g, j)] for g in perts for j in JJ]),
                                  "Jd_marc": resumen([Jdm[(L, g, j)] for g in perts for j in JJ]),
                                  "fd": resumen([fd[(L, g, j)] for g in perts for j in JJ]),
                                  "n_celdas_union_vacia": sum(1 for g in perts for j in JJ if Jd[(L, g, j)] is None)}
        res["dJ_linea"][L] = resumen([m(cdJ[p + (g,)]) for p in pares if L in p for g in perts])
        res["df_linea"][L] = resumen([m(cdf[(L, L2, g)]) for L2 in LINES if L2 != L for g in perts])
    for p in pares:
        k = f"{p[0]}–{p[1]}"
        res["entre_par"][k] = {"Je": resumen([Je[p + (g, j)] for g in perts for j in JJ])}
        res["dJ_par"][k] = resumen([m(cdJ[p + (g,)]) for g in perts])
        res["dJ_par"][k]["ancho_rango_j_mediana"] = float(np.median([cdJ[p + (g,)][2] - cdJ[p + (g,)][1]
                                                                      for g in perts if cdJ[p + (g,)]])) if any(cdJ[p + (g,)] for g in perts) else None
        res["dJm_par"][k] = resumen([m(cdJm[p + (g,)]) for g in perts])
    for p in opares:
        k = f"{p[0]}→{p[1]}"
        res["df_opar"][k] = resumen([m(cdf[p + (g,)]) for g in perts])
        res["fe_opar"][k] = resumen([fe[p + (g, j)] for g in perts for j in JJ])
    for g in perts:
        res["dJ_pert"][g] = resumen([m(cdJ[p + (g,)]) for p in pares])
        res["dJ_pert"][g]["Jd_por_linea"] = {L: m(cJd[(L, g)]) for L in LINES}
        res["dJ_pert"][g]["Jd_rango_j_por_linea"] = {L: None if cJd[(L, g)] is None else [cJd[(L, g)][1], cJd[(L, g)][2]] for L in LINES}
    # intermedios: 20 celdas (par, g) equiespaciadas por dJ
    cel = [(p, g, cdJ[p + (g,)]) for p in pares for g in perts if cdJ[p + (g,)] is not None]
    cel.sort(key=lambda x: x[2][0])
    sel = [cel[i] for i in np.linspace(0, len(cel) - 1, 20).round().astype(int)] if cel else []
    res["intermedios"] = [{"par": f"{p[0]}–{p[1]}", "pert": g, "Jd_L1": m(cJd[(p[0], g)]), "Jd_L2": m(cJd[(p[1], g)]),
                           "Je": m(cJe[p + (g,)]), "dJ_media_j": c[0], "dJ_min_j": c[1], "dJ_max_j": c[2], "n_j": c[3],
                           "tam_L1_A_j0": len(R[p[0]][0]["A"][g]), "tam_L2_B_j0": len(R[p[1]][0]["B"][g])}
                          for p, g, c in sel]
    # intermedios de una celda por particion: la de mediana dJ
    if cel:
        p, g, _ = cel[len(cel) // 2]
        res["intermedios_j"] = {"par": f"{p[0]}–{p[1]}", "pert": g, "filas": [
            {"j": j, "|L1A|": len(R[p[0]][j]["A"][g]), "|L1B|": len(R[p[0]][j]["B"][g]),
             "|L2A|": len(R[p[1]][j]["A"][g]), "|L2B|": len(R[p[1]][j]["B"][g]),
             "Jd_L1": Jd[(p[0], g, j)], "Jd_L2": Jd[(p[1], g, j)], "Je": Je[p + (g, j)], "dJ": dJ[p + (g, j)]}
            for j in JJ]}
    res["n_celdas_par_g"] = len(pares) * len(perts)
    res["n_celdas_par_g_definidas"] = len(cel)
    return res


def analiza_2025(R, perts):
    Jd, Jdm, fd = {}, {}, {}
    for g in perts:
        for j in JJ:
            A, B = R["v2025"][j]["A"][g], R["v2025"][j]["B"][g]
            Jd[(g, j)] = jac(A, B); Jdm[(g, j)] = jac(A, B, "marc"); fd[(g, j)] = mean_def([fno(A, B), fno(B, A)])
    porg = {}
    for g in perts:
        v = [Jd[(g, j)] for j in JJ if Jd[(g, j)] is not None]
        w = [fd[(g, j)] for j in JJ if fd[(g, j)] is not None]
        porg[g] = {"Jd": (float(np.mean(v)), min(v), max(v), len(v)) if v else None,
                   "fd": (float(np.mean(w)), min(w), max(w)) if w else None,
                   "tam_medio": float(np.mean([len(R["v2025"][j][h][g]) for j in JJ for h in "AB"]))}
    return {"Jd": resumen([Jd[k] for k in Jd]), "Jd_marc": resumen([Jdm[k] for k in Jdm]),
            "fd": resumen([fd[k] for k in fd]), "n_celdas_union_vacia": sum(1 for k in Jd if Jd[k] is None),
            "Jd_pert": resumen([x["Jd"][0] if x["Jd"] else None for x in porg.values()]), "por_pert": porg}


def tamanos(D, full, tags, perts_of):
    out = {}
    for t in tags:
        perts = perts_of[t]
        sums = [sum(len(D[t][j][h].get(g, {})) for g in perts) for j in JJ for h in "AB"]
        F = sum(len(full[t].get(g, [])) for g in perts)
        porp = {}
        for g in perts:
            hs = [len(D[t][j][h].get(g, {})) for j in JJ for h in "AB"]
            porp[g] = {"full": len(full[t].get(g, [])), "half_mediana": float(np.median(hs)),
                       "half_min": int(min(hs)), "half_max": int(max(hs))}
        out[t] = {"full": F, "half_media": float(np.mean(sums)), "half_min": int(min(sums)), "half_max": int(max(sums)),
                  "ratio_media": float(np.mean(sums)) / F if F else None,
                  "cortes_j0": D[t][0]["cortes"].get("non-targeting"), "por_pert": porp}
    return out


def main():
    D = carga()
    faltan = [(t, j) for t in D for j in JJ if j not in D[t]]
    assert not faltan, faltan
    venn = json.load(open("/home/ubuntu/vcc/modelo/out/venn_ed.json"))
    p2 = json.load(open("/datos/modelo/out/seleccion/p2_datos.json"))
    panel = venn["panel"]
    perts25 = p2["v2025"]["perts"]
    res = {"semillas": S.comprueba(), "verif": {t: D[t][0].get("verif") for t in D},
           "meta": {t: json.load(open(os.path.join(SETS, f"{t}_meta.json"))) for t in D}}
    res["cortes"] = {t: {j: D[t][j]["cortes"]["non-targeting"] for j in JJ} for t in D}
    res["cortes_pert_impar_2025"] = {g: [D["v2025"][j]["cortes"][g] for j in JJ]
                                     for g in perts25 if D["v2025"][0]["cortes"][g][0] != D["v2025"][0]["cortes"][g][1]
                                     or D["v2025"][1]["cortes"][g][0] != D["v2025"][1]["cortes"][g][1]}
    Rc = crudo(D, LINES, panel)
    Ri, ms = recorte(D, LINES, panel, S.SEMILLAS["desempate_recorte_jiang"][0])
    res["jiang_crudo"] = analiza_jiang(Rc, panel)
    res["jiang_igualado"] = analiza_jiang(Ri, panel)
    res["m_recorte_jiang"] = {g: [ms[(g, j)] for j in JJ] for g in panel}
    R5c = crudo(D, ["v2025"], perts25)
    R5i, ms5 = recorte(D, ["v2025"], perts25, S.SEMILLAS["desempate_recorte_2025"][0])
    res["v2025_crudo"] = analiza_2025(R5c, perts25)
    res["v2025_igualado"] = analiza_2025(R5i, perts25)
    full = {L: venn["ed"][L] for L in LINES}
    full["v2025"] = p2["v2025"]["ed"]
    po = {L: panel for L in LINES}; po["v2025"] = perts25
    res["tamanos"] = tamanos(D, full, LINES + ["v2025"], po)
    res["perts25"] = perts25; res["panel"] = panel
    res["profundidad"] = {L: p2["jiang"]["crudo"][L]["L_mediana"] for L in LINES}
    res["profundidad"]["v2025"] = p2["v2025"]["L_mediana"]
    res["dt_cli_total_s"] = sum(D[t][j]["dt_total_s"] for t in D for j in JJ)
    json.dump(res, open(OUT_JSON, "w"), default=lambda o: list(o) if isinstance(o, (set, tuple)) else str(o))
    return res


if __name__ == "__main__":
    main()
