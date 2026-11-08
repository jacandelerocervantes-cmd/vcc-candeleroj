#!/usr/bin/env python
"""Reporte `magnitud_centro.md` -- la fraccion CENTRO."""
import json, os
import numpy as np
OUTDIR = "/home/ubuntu/vcc/modelo/out/magnitud_centro"
VAR = [("delta_crudo", "el delta CRUDO"),
       ("delta_especifico_del_blanco", "el delta ESPECIFICO del blanco")]


def main():
    d = json.load(open(os.path.join(OUTDIR, "calc.json")))
    c = json.load(open(os.path.join(OUTDIR, "ctrl.json")))
    fc = d["fraccion_centro"]
    L = []; W = L.append

    W("## 4. 🔴 La fraccion CENTRO: `V_sh / V_tot` entre contextos")
    W("")
    W("### 4.1 La identidad, y por que el producto cruzado de mitades")
    W("")
    W("La identidad de `b_contra_a.md` §5, con el indice corrido sobre **LINEAS** en")
    W("vez de sobre blancos. Para el blanco `g` y un conjunto de lineas `S`, con")
    W("`h1_L = delta_g,L,A` y `h2_L = delta_g,L,B`:")
    W("")
    W("```")
    W("V_tot(g) = (1/|S|) Σ_L <h1_L, h2_L> / G      senal reproducible total")
    W("V_sh(g)  = <m1, m2> / G                      el CENTRO (m = media sobre lineas)")
    W("V_sp(g)  = (1/|S|) Σ_L <h1_L-m1, h2_L-m2> / G   ESPECIFICA de la linea")
    W("```")
    W("")
    W("`V_tot = V_sh + V_sp` es una **identidad algebraica exacta**, no un ajuste: los")
    W("terminos cruzados se anulan porque `Σ_L (h1_L - m1) = 0`. Y el producto")
    W("**cruzado de mitades** es lo que quita el ruido **sin modelarlo**: su esperanza")
    W("es el producto de las senales verdaderas, porque el ruido de las dos mitades es")
    W("independiente (celulas disjuntas y medios controles disjuntos). La varianza de")
    W("**una** mitad llevaria el ruido dentro.")
    W("")
    W("**El residuo de la identidad, medido** (si no fuera ~0 la descomposicion estaria mal):")
    W("")
    W("| # | variante | conjunto | residuo abs. max | residuo relativo max |")
    W("|---|---|---|---|---|")
    i = 0
    for k, nom in VAR:
        for etiq in ("las_6_lineas",):
            i += 1
            v = fc[k][etiq]
            W(f"| {i} | {nom} | las 6 lineas | {v['residuo_identidad_max']:.2e} | "
              f"{v['residuo_identidad_rel_max']:.2e} |")
        pares = [e for e in fc[k] if e != "las_6_lineas"]
        rm = max(fc[k][e]["residuo_identidad_max"] for e in pares)
        rr = max(fc[k][e]["residuo_identidad_rel_max"] for e in pares)
        i += 1
        W(f"| {i} | {nom} | peor de los 15 pares | {rm:.2e} | {rr:.2e} |")
    W("")
    W("### 4.2 🔴 El control negativo que la fraccion necesita, y por que")
    W("")
    W("`V_sh` mide la norma reproducible de la **media empirica** de los deltas sobre")
    W("las lineas. Si las lineas no compartieran NADA del delta especifico del blanco,")
    W("esa media **no seria cero**: la media de `n` vectores mutuamente incorrelados")
    W("conserva ~`1/n` de la norma cuadrada tipica, y **el producto cruzado de mitades**")
    W("**no lo quita**, porque el efecto especifico de linea es la misma cosa real en")
    W("las dos mitades (celulas disjuntas de la misma linea, no ruido independiente).")
    W("")
    W(f"Asi que `V_sh/V_tot` tiene un **suelo de no-comparticion**, y se **MIDE** en vez")
    W("de suponerse: para cada blanco `g` se toma, en cada linea, el delta de un blanco")
    W("DISTINTO por una permutacion **independiente por linea y sin punto fijo**. Eso")
    W("destruye la comparticion atribuible al blanco y **deja intacto el eje comun al**")
    W(f"**IFNB**. {c['n_permutaciones']} sorteos, semillas")
    W(f"`{c['semilla_perm_base']} + k`.")
    W("")
    for k, nom in VAR:
        W(f"**TABLA 4.{1 if k == 'delta_crudo' else 2} — la fraccion CENTRO sobre {nom}.**")
        W("")
        W("| # | conjunto | lineas | `V_tot` | `V_sh` | `V_sp` | "
          "**frac. observada** | **IC 95 % (bootstrap sobre perts)** | "
          "**control permutado** (sd) | suelo analitico `1/n` | "
          "**frac. ANCLADA en el control** | frac. anclada en `1/n` |")
        W("|---|---|---|---|---|---|---|---|---|---|---|---|")
        orden = ["las_6_lineas"] + sorted(e for e in fc[k] if e != "las_6_lineas")
        for i, etiq in enumerate(orden, 1):
            v = fc[k][etiq]; cv = c["variantes"][k][etiq]
            b = v["bootstrap_sobre_perturbaciones"]
            ic = (f"[{b['ic95'][0]:+.4f}, {b['ic95'][1]:+.4f}]" if b["ic95"] else "n/d")
            W(f"| {i} | **{etiq.replace('_', ' ')}** | {v['n_lineas']} | "
              f"{v['V_tot_media']:.6f} | {v['V_sh_media']:.6f} | {v['V_sp_media']:.6f} | "
              f"**{v['frac_centro_agregada']:+.4f}** | {ic} | "
              f"{cv['frac_control_permutado_media']:+.4f} "
              f"({cv['frac_control_permutado_sd']:.4f}) | "
              f"{cv['suelo_analitico_1_sobre_n']:.4f} | "
              f"**{cv['frac_anclada_en_el_control']:+.4f}** | "
              f"{cv['frac_anclada_en_1_sobre_n']:+.4f} |")
        W("")
    W("### 4.3 Los rangos sobre los 15 pares, y la cifra agregada")
    W("")
    W("| # | variante | cantidad | agregada (las 6 lineas) | **[min, max] sobre los 15 pares** |")
    W("|---|---|---|---|---|")
    i = 0
    for k, nom in VAR:
        pares = sorted(e for e in fc[k] if e != "las_6_lineas")
        for lab, get in (("frac. observada",
                          lambda e: fc[k][e]["frac_centro_agregada"]),
                         ("control permutado",
                          lambda e: c["variantes"][k][e]["frac_control_permutado_media"]),
                         ("**frac. ANCLADA en el control**",
                          lambda e: c["variantes"][k][e]["frac_anclada_en_el_control"])):
            i += 1
            vs = [get(e) for e in pares]
            ag = get("las_6_lineas")
            W(f"| {i} | {nom} | {lab} | **{ag:+.4f}** | "
              f"**[{min(vs):+.4f}, {max(vs):+.4f}]** |")
    W("")
    k = "delta_especifico_del_blanco"
    pares = sorted(e for e in fc[k] if e != "las_6_lineas")
    anc = [c["variantes"][k][e]["frac_anclada_en_el_control"] for e in pares]
    ancA = c["variantes"][k]["las_6_lineas"]["frac_anclada_en_el_control"]
    obsA = fc[k]["las_6_lineas"]["frac_centro_agregada"]
    ctlA = c["variantes"][k]["las_6_lineas"]["frac_control_permutado_media"]
    kc = "delta_crudo"
    obsC = fc[kc]["las_6_lineas"]["frac_centro_agregada"]
    ctlC = c["variantes"][kc]["las_6_lineas"]["frac_control_permutado_media"]
    ancC = c["variantes"][kc]["las_6_lineas"]["frac_anclada_en_el_control"]
    W("🔴 **Lo que estas tablas dicen, como aritmetica:**")
    W("")
    W("* 🔴 **El control permutado separa las dos variantes, y eso es la TRAMPA 2**")
    W("  **medida.** Sobre el delta **CRUDO** y las seis lineas el control queda en")
    W(f"  **{ctlC:+.4f}**, **por encima** del suelo analitico `1/6 = 0.1667`: esa")
    W(f"  diferencia de **{ctlC - 1 / 6:+.4f}** es el eje comun al IFNB, que sobrevive a")
    W("  la permutacion de blancos porque es compartido entre blancos Y entre lineas.")
    W("  Sobre el delta **ESPECIFICO del blanco** el control cae a")
    W(f"  **{ctlA:+.4f}**, que es `1/6` **a tres decimales** -- restar la media sobre")
    W("  blancos dentro de cada linea le quito al vector exactamente lo que el control")
    W("  permutado detectaba. Las dos mediciones se cruzan y dan lo mismo por dos")
    W("  caminos.")
    W(f"* La fraccion OBSERVADA es **{obsC:+.4f}** (delta crudo) y **{obsA:+.4f}**")
    W("  (delta especifico), pero **ninguna de las dos es la fraccion de centro**:")
    W("  llevan dentro el suelo de no-comparticion. Las comparables son las ANCLADAS:")
    W(f"  **{ancC:+.4f}** y **{ancA:+.4f}** sobre las seis lineas.")
    W(f"* Sobre los 15 pares, la fraccion anclada va de **{min(anc):+.4f} a")
    W(f"  {max(anc):+.4f}** (delta especifico del blanco). El par mas compartido es")
    W(f"  `{[e for e in pares if c['variantes'][k][e]['frac_anclada_en_el_control'] == max(anc)][0]}`")
    W(f"  y el menos `{[e for e in pares if c['variantes'][k][e]['frac_anclada_en_el_control'] == min(anc)][0]}`.")
    W("* 🔴 **La fraccion OBSERVADA por PAR es sistematicamente mayor que la de las**")
    W("  **seis lineas, y es aritmetica del estimador, no biologia**: con `n = 2` la")
    W("  media empirica conserva la mitad de la norma por si sola (suelo 0.50) y con")
    W("  `n = 6` solo un sexto (0.1667). Por eso las observadas de la TABLA 4.2 **no se**")
    W("  **pueden comparar entre conjuntos de distinto `n`** y las ancladas si.")
    conj = ["las_6_lineas"] + pares
    marg = []
    for e in conj:
        lo = fc[k][e]["bootstrap_sobre_perturbaciones"]["ic95"][0]
        ct = c["variantes"][k][e]["frac_control_permutado_media"]
        marg.append((lo - ct, e, lo, ct))
    ok = sum(1 for m, *_ in marg if m > 0)
    peor = min(marg)
    W("* **El IC del bootstrap se compara con el control DE SU PROPIO conjunto**, no")
    W("  con el de otro (los suelos difieren con `n`). El extremo inferior del IC supera")
    W(f"  al control permutado de su conjunto en **{ok} de {len(conj)}**. El margen mas")
    W(f"  estrecho es `{peor[1]}`: IC inferior **{peor[2]:+.4f}** contra control")
    W(f"  **{peor[3]:+.4f}**, margen **{peor[0]:+.4f}**.")
    W("")
    W("### 4.4 Por perturbacion")
    W("")
    W("**TABLA 4.3 — `V_tot`, `V_sh`, `V_sp` y la fraccion por perturbacion, sobre las")
    W("seis lineas, media de las 8 particiones.** Variante del delta especifico del")
    W("blanco. Ordenada por `V_tot`.")
    W("")
    W("| # | pert | `V_tot` | `V_sh` | `V_sp` | **frac. centro** | sd sobre particiones |")
    W("|---|---|---|---|---|---|---|")
    rows = sorted(fc[k]["las_6_lineas"]["por_pert"], key=lambda r: -r["V_tot"])
    for i, r in enumerate(rows, 1):
        fcv = (f"**{r['frac_centro']:+.4f}**" if r["frac_centro"] is not None else "n/d")
        W(f"| {i} | `{r['pert']}` | {r['V_tot']:.6f} | {r['V_sh']:.6f} | "
          f"{r['V_sp']:.6f} | {fcv} | "
          f"{r['frac_centro_sd_sobre_particiones']:.4f} |")
    W("")
    fr = [r["frac_centro"] for r in rows if r["frac_centro"] is not None]
    W(f"Rango de la fraccion por perturbacion: **[{min(fr):+.4f}, {max(fr):+.4f}]**,")
    W(f"con **{sum(1 for x in fr if x < 0)} de {len(fr)}** negativas (una fraccion")
    W("negativa sale cuando `V_sh` y `V_tot` tienen signos distintos, que a")
    W("`V_tot` pequeno es ruido de estimacion, no una cantidad interpretable).")
    W("")
    W("---")
    W("")
    return L


if __name__ == "__main__":
    open("/tmp/claude-1000/-home-ubuntu-vcc/52816e9e-e42f-4af3-b6ec-69a5b07f6ce9/scratchpad/mcrep_d.txt",
         "w").write("\n".join(main()))
    print("ok")
