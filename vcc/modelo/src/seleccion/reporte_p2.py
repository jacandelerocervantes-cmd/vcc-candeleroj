#!/usr/bin/env python
"""Tablas de la PARTE 2 (Jiang y VCC 2025, superficies con verdad). Lee
`p2_datos.json` y `p2_mecanismo.json`. Ninguna tabla de aqui lleva cifras de
la Parte 1."""
import json

import numpy as np

from reporte_p1 import tabla

DAT = "/datos/modelo/out/seleccion/p2_datos.json"
MEC = "/datos/modelo/out/seleccion/p2_mecanismo.json"
LINES = ["A549", "BXPC3", "HAP1", "HT29", "K562", "MCF7"]
CLAVES = [("p1", "solo_p1", "acum_p1"), ("p2", "solo_p2", "acum_p1+p2"),
          ("p3", "solo_p3", "acum_p1+p2+p3"), ("p4", "solo_p4", "acum_p1+p2+p3+p4")]
NOMBRE = {"p1": "p1 pasa el filtro de L2", "p2": "p2 CPM medio (deciles)",
          "p3": "p3 tasa de deteccion (deciles)", "p4": "p4 potencia (log n, log UMI mediana)"}


def anc(o, nul):
    m = float(np.mean(nul))
    return (o - m) / (1 - m) if m < 1 else float("nan")


def nstat(nul):
    return f"{np.mean(nul):.4f} [{min(nul):.4f}, {max(nul):.4f}]"


def main():
    d = json.load(open(DAT)); m = json.load(open(MEC))
    J = d["jiang"]; U = d["umbral"]
    W = ["## PARTE 2 — JIANG Y VCC 2025 (CON VERDAD): EL MECANISMO", "",
         "🔴 De esta parte viaja solo el MECANISMO. Ninguna fraccion de aqui es una cifra de D/E/F.", ""]
    # ------------- construccion verificada
    W += ["### 2.0 Las poblaciones, verificadas", ""]
    filas = []
    for i, L in enumerate(LINES):
        c, g = J["crudo"][L], J["igualado"][L]
        filas.append([i + 1, f"**{L}**", "si" if c["huella_coincide"] else "**NO**", f"{c['perts_identicas_venn']}/33",
                      f"{c['n_control']:,}", f"{c['L_mediana']:,.0f}", f"{c['n_pasa']:,}",
                      sum(len(v) for v in c["ed"].values()),
                      f"{g['p_adelgaza']:.4f}", f"{g['n_control']:,}", f"{g['L_mediana']:,.0f}", f"{g['n_pasa']:,}",
                      sum(len(v) for v in g["ed"].values()), f"{g['semilla_control']} / {g['semilla_adelgaza']}"])
    W += tabla(["# ", "linea", "huella real = venn_ed", "R(g,L) = venn_ed gen a gen", "n control", "UMI mediana",
                "# pasan filtro", "Σ|R| crudo", "p adelgaza", "n control ig.", "UMI mediana ig.", "# pasan ig.",
                "Σ|R| igualado", "semillas ig. (control / adelgaza)"], filas)
    v = d["v2025"]
    W += ["", f"Igualacion: n_min = {J['igualacion']['n_min']:,} controles, D_min = {J['igualacion']['D_min']:,.0f} UMI "
          f"mediana/celula. Genes del espacio Jiang: {len(J['genes']):,}.", "",
          f"VCC 2025 (`test`): huella del real = la de `real_2025_test_meta.json`: "
          f"**{'si' if v['huella_coincide'] else 'NO'}**; control {v['n_control']:,} celulas, UMI mediana "
          f"{v['L_mediana']:,.0f}; filtrables {v['n_pasa']:,} de {len(v['genes']):,}; Σ|R(g)| = "
          f"{sum(len(x) for x in v['ed'].values()):,} sobre {len(v['perts'])} perturbaciones.", "",
          f"Chequeo de construccion: genes de R(g,L1) que NO pasan el filtro de L1 — crudo "
          f"**{m['crudo']['viol_filtro_L1']}**, igualado **{m['igualado']['viol_filtro_L1']}** (tiene que ser 0).", ""]
    for var, etq in (("crudo", "CRUDO"), ("igualado", "IGUALADO EN POTENCIA")):
        r = m[var]
        W += [f"Unidades (g, gen, L1, L2) — {etq}: **{r['n_filas']:,}** filas, de ellas **{r['n_y1']:,}** "
              f"especificas (gen en R(g,L1), no en R(g,L2)) = {100*r['frac_y1']:.2f} %.", ""]
    # ------------- 2.1
    W += ["### 2.1 De los genes en R(g,L1) y no en R(g,L2): que fraccion NO PASA el filtro de L2", "",
          "Sin biologia: un gen que no pasa el filtro de L2 no se PRUEBA en L2, asi que no puede estar en R(g,L2).", ""]
    for var, etq in (("crudo", "CRUDO"), ("igualado", "IGUALADO")):
        r = m[var]
        nul = [x["A_det_p1"] for x in r["nulo"]]
        o = r["obs"]["A_det_p1"]
        W += [f"**{etq}** — agregado: **{100*o:.2f} %** de {r['n_y1']:,}; nulo (20 permutaciones del contexto de "
              f"las estadisticas de control) {100*np.mean(nul):.2f} % [{100*min(nul):.2f}, {100*max(nul):.2f}]; "
              f"anclado **{anc(o, nul):.4f}**.", ""]
        filas = []
        for i, p in enumerate(r["por_par"]):
            filas.append([i + 1, p["L1"], p["L2"], f"{p['n_filas']:,}", f"{p['n_especificos']:,}", f"{p['n_no_pasan_L2']:,}",
                          f"{100*p['n_no_pasan_L2']/p['n_especificos']:.1f} %" if p["n_especificos"] else "—"])
        W += tabla(["# ", "L1", "L2", "# en R(g,L1)", "# no en R(g,L2)", "# de esos que no pasan filtro L2", "fraccion"], filas)
        W += [""]
    W += ["🔴 **Intermedios (CRUDO): 20 filas especificas**, equiespaciadas por CPM del control de L2", ""]
    W += tabla(["# ", "pert", "gen", "L1", "L2", "CPM ctrl L1", "CPM ctrl L2", "pasa L2", "det ctrl L2", "n ctrl L2", "UMI med L2"],
               [[i + 1, f"`{x['pert']}`", f"`{x['gen']}`", x["L1"], x["L2"], f"{x['cpm_L1']:.3f}", f"{x['cpm_L2']:.3f}",
                 "si" if x["cpm_L2"] > U else "no", f"{x['det_L2']:.4f}", x["n_L2"], f"{x['Lmed_L2']:.0f}"]
                for i, x in enumerate(m["crudo"]["intermedios"])])
    W += [""]
    # ------------- 2.2
    W += ["### 2.2 La atribucion: predictores del control de L2, orden DECLARADO p1 → p2 → p3 → p4", "",
          "D² = 1 − devianza(modelo)/devianza(solo intercepto), logistica. Las dos columnas: cada predictor SOLO, y "
          "el ACUMULADO en el orden declarado. Nulo: media [min, max] de 20 sorteos. Anclado = (obs − media nulo)/(1 − media nulo).", ""]
    for var, etq in (("crudo", "CRUDO"), ("igualado", "IGUALADO EN POTENCIA")):
        r = m[var]
        filas = []
        prev = 0.0
        for p, ks, ka in CLAVES:
            os_, oa = r["obs"][ks], r["obs"][ka]
            ns, na = [x[ks] for x in r["nulo"]], [x[ka] for x in r["nulo"]]
            nota = "" if (p != "p4" or r["usar_p4"]) else " (constante tras igualar: omitido)"
            filas.append([NOMBRE[p] + nota, f"**{os_:.4f}**", nstat(ns), f"{anc(os_, ns):.4f}",
                          f"**{oa:.4f}**", nstat(na), f"{anc(oa, na):.4f}", f"{oa - prev:+.4f}"])
            prev = oa
        W += [f"**{etq}**", ""]
        W += tabla(["predictor", "SOLO D²", "SOLO nulo", "SOLO anclado", "ACUM D²", "ACUM nulo", "ACUM anclado",
                    "incremento ACUM"], filas)
        W += [""]
    W += ["🔴 **Intermedios del nulo (CRUDO): los 20 sorteos**", ""]
    W += tabla(["k", "semilla", "permutacion (L2 → π(L2), indices en " + "/".join(LINES) + ")", "# puntos fijos",
                "A det p1", "SOLO p1", "SOLO p2", "SOLO p3", "SOLO p4", "ACUM p1..p4"],
               [[i, x["semilla"], str(x["perm"]), x["fijos"], f"{x['A_det_p1']:.4f}", f"{x['solo_p1']:.4f}",
                 f"{x['solo_p2']:.4f}", f"{x['solo_p3']:.4f}", f"{x['solo_p4']:.4f}", f"{x['acum_p1+p2+p3+p4']:.4f}"]
                for i, x in enumerate(m["crudo"]["nulo"])])
    W += ["", f"Bordes de los deciles (CRUDO): p2 log10(CPM+1e-3) = {[round(x, 3) for x in m['crudo']['bordes']['p2']]}; "
          f"p3 deteccion = {[round(x, 4) for x in m['crudo']['bordes']['p3']]}.", ""]
    # ------------- 2.3
    W += ["### 2.3 🔴 EL RESIDUO: la especificidad de contexto que NO explica ninguna cantidad del control", ""]
    filas = []
    for var, etq in (("crudo", "CRUDO"), ("igualado", "IGUALADO EN POTENCIA")):
        r = m[var]
        k = "acum_p1+p2+p3+p4"
        o = r["obs"][k]; nul = [x[k] for x in r["nulo"]]
        oA = r["obs"]["A_det_p1"]; nA = [x["A_det_p1"] for x in r["nulo"]]
        filas.append([etq, f"{r['n_y1']:,}", f"**{1-o:.4f}**", f"**{1-anc(o, nul):.4f}**",
                      f"{1-oA:.4f}", f"{1-anc(oA, nA):.4f}"])
    W += tabla(["R(g,L)", "# especificas", "**residuo B = 1 − D² acum (crudo)**", "**residuo B anclado = 1 − D² acum anclado**",
                "residuo A = 1 − frac. que no pasa filtro L2", "residuo A anclado"], filas)
    W += [""]
    # ------------- 2.4
    W += ["### 2.4 El universo U(g)", "",
          "Jiang: U(g) = ∪_L R(g,L) sobre las 6 lineas. VCC 2025: un solo contexto, asi que se da |R(g)| / |filtrables|.", ""]
    uj = m["universo"]
    for var in ("crudo", "igualado"):
        pp = uj[var]["por_pert"]
        us = [x["U"] for x in pp]
        W += [f"**Jiang {var.upper()}**: |U(g)| rango **[{min(us)}, {max(us)}]**, mediana {int(np.median(us))}; "
              f"filtrables en al menos una linea: {uj[var]['filtrables_alguna_linea']:,} (|U(g)|/eso: "
              f"[{min(us)/uj[var]['filtrables_alguna_linea']:.4f}, {max(us)/uj[var]['filtrables_alguna_linea']:.4f}]).", ""]
    pc, pi_ = uj["crudo"]["por_pert"], uj["igualado"]["por_pert"]
    W += ["🔴 **Intermedios: las 33 perturbaciones**", ""]
    W += tabla(["# ", "pert", "|U| crudo"] + [f"|R| {L}" for L in LINES] + ["en las 6", "# lineas no vacias", "|U| igualado"],
               [[i + 1, f"`{a['pert']}`", f"**{a['U']}**"] + [a["por_linea"][L] for L in LINES] +
                [a["en_6"], a["n_lineas_que_aportan"], b["U"]] for i, (a, b) in enumerate(zip(pc, pi_))])
    v5 = uj["v2025"]
    fr = [x["frac_de_filtrables"] for x in v5["por_pert"]]
    Rs = [x["R"] for x in v5["por_pert"]]
    W += ["", f"**VCC 2025**: filtrables {v5['n_filtrables']:,}. |R(g)| rango [{min(Rs)}, {max(Rs)}], mediana {int(np.median(Rs))}; "
          f"fraccion de los filtrables en R(g): rango **[{min(fr):.4f}, {max(fr):.4f}]**, mediana {np.median(fr):.4f}. "
          f"Perts con R vacio: {sum(1 for x in Rs if x == 0)}. Genes de R(g) fuera de los filtrables (tiene que ser 0): "
          f"{sum(x['R_fuera_filtrables'] for x in v5['por_pert'])}. ∪_g R(g) sobre las {len(Rs)} perts: {v5['union_R']:,} = "
          f"{v5['union_R_frac']:.4f} de los filtrables.", ""]
    o = sorted(v5["por_pert"], key=lambda x: x["R"])
    sel = [o[j] for j in np.linspace(0, len(o) - 1, 20).round().astype(int)]
    W += ["🔴 **Intermedios 2025: 20 perturbaciones equiespaciadas por |R(g)|**", ""]
    W += tabla(["# ", "pert", "celulas", "|R(g)|", "fraccion de filtrables"],
               [[i + 1, f"`{x['pert']}`", x["n_celulas"], x["R"], f"{x['frac_de_filtrables']:.4f}"] for i, x in enumerate(sel)])
    W += [""]
    return W, d, m


if __name__ == "__main__":
    W, _, _ = main()
    print("\n".join(W))
