#!/usr/bin/env python
"""Reporte `magnitud_centro.md` -- signo y magnitud, POR SEPARADO."""
import json, os
import numpy as np
OUTDIR = "/home/ubuntu/vcc/modelo/out/magnitud_centro"
NUCLEO = ["IFNAR1", "IFNAR2", "JAK1", "STAT1", "TYK2"]


def main():
    d = json.load(open(os.path.join(OUTDIR, "calc.json")))
    pares = d["por_par"]; pp = d["por_par_y_pert"]
    L = []; W = L.append

    W("## 5. Signo y magnitud, POR SEPARADO")
    W("")
    W("🔴 **Nunca una resumiendo a la otra.** La regla del artefacto las distingue, asi")
    W("que van en columnas distintas y con rangos distintos:")
    W("")
    W("| | pregunta | con que se mide |")
    W("|---|---|---|")
    W("| **ORDEN / SIGNO** | ¿transfiere el orden? | **concordancia de signo** (azar = "
      "0.5 exactamente) y **Spearman** entre contextos, los dos con mitades cruzadas |")
    W("| **MAGNITUD** | ¿transfiere la magnitud? | **pendiente** de la regresion entre "
      "contextos (`<dA,dB>/<dA,dA>`) y su **R²** |")
    W("")
    W("La pendiente se atenua por la fiabilidad del **PREDICTOR solo**, no por la raiz")
    W("del producto, asi que su desatenuada es `pendiente / fiab_L1`. Se da cruda y")
    W("desatenuada.")
    W("")
    W("**TABLA 5.1 — signo y magnitud por par de lineas, media sobre 8 x 33 celdas.**")
    W("")
    W("| # | par | **concordancia de signo** | exceso sobre el azar | **Spearman** | "
      "**pendiente cruda** | **pendiente desat.** | **R²** |")
    W("|---|---|---|---|---|---|---|---|")
    for i, p in enumerate(pares, 1):
        sg = p["concordancia_signo"]["media"]
        pd_ = p["pendiente_desat"]["media"]
        W(f"| {i} | **{p['par']}** | **{sg:.4f}** | **{sg - 0.5:+.4f}** | "
          f"**{p['spearman']['media']:+.4f}** | "
          f"**{p['pendiente']['media']:+.4f}** | "
          f"{'**%+.4f**' % pd_ if pd_ is not None else 'n/d'} | "
          f"**{p['R2']['media']:.5f}** |")
    W("")
    W("**Los rangos sobre los 15 pares:**")
    W("")
    W("| # | cantidad | **[min, max]** | media | referencia |")
    W("|---|---|---|---|---|")
    campos = [("concordancia de signo", "concordancia_signo", "azar = **0.5000**"),
              ("Spearman", "spearman", "ausencia = **0.0000**"),
              ("pendiente cruda", "pendiente", "transferencia perfecta = **1.0000**"),
              ("pendiente desatenuada", "pendiente_desat",
               "transferencia perfecta = **1.0000**"),
              ("R²", "R2", "transferencia perfecta = **1.0000**")]
    for i, (lab, k, ref) in enumerate(campos, 1):
        vs = [p[k]["media"] for p in pares if p[k]["media"] is not None]
        W(f"| {i} | {lab} | **[{min(vs):+.4f}, {max(vs):+.4f}]** | "
          f"{float(np.mean(vs)):+.4f} | {ref} |")
    W("")
    sg = [p["concordancia_signo"]["media"] for p in pares]
    sp = [p["spearman"]["media"] for p in pares]
    pe = [p["pendiente"]["media"] for p in pares]
    pdv = [p["pendiente_desat"]["media"] for p in pares
           if p["pendiente_desat"]["media"] is not None]
    r2 = [p["R2"]["media"] for p in pares]
    W("🔴 **Lo que estas filas dicen, como aritmetica:**")
    W("")
    W(f"* **Signo**: la concordancia va de **{min(sg):.4f} a {max(sg):.4f}**, es decir")
    W(f"  **{min(sg) - 0.5:+.4f} a {max(sg) - 0.5:+.4f}** sobre el azar, y esta por")
    W(f"  encima de 0.5 en **los 15 pares de 15**. El Spearman va de")
    W(f"  **{min(sp):+.4f} a {max(sp):+.4f}**, positivo en los 15.")
    W(f"* **Magnitud**: la pendiente cruda va de **{min(pe):+.4f} a {max(pe):+.4f}** --")
    W("  lejos de 1 en los 15 pares -- y el **R² va de")
    W(f"  {min(r2):.5f} a {max(r2):.5f}**, es decir que entre el **0.15 % y el")
    W(f"  {max(r2) * 100:.2f} %** de la varianza del delta de un contexto se explica")
    W("  desde el otro.")
    W(f"* 🔴 **La pendiente desatenuada es INESTABLE**: va de **{min(pdv):+.4f} a")
    W(f"  {max(pdv):+.4f}**, y en el par `"
      f"{[p['par'] for p in pares if p['pendiente_desat']['media'] == min(pdv)][0]}` sale")
    W("  **negativa**. Dividir por una fiabilidad de centesimas es lo que produce eso.")
    W("  No se recorta ni se omite.")
    W("")
    W("**TABLA 5.2 — signo y magnitud partidos por FIABILIDAD**, con el mismo corte de")
    W("la seccion 3.3.")
    W("")
    W("| # | subconjunto | n perts | n filas | `fiab` media | **concordancia de signo** | "
      "**Spearman** | **pendiente cruda** | **pendiente desat.** | **R²** |")
    W("|---|---|---|---|---|---|---|---|---|---|")
    for i, (nom, sel) in enumerate((("fiable (`>= 0.15`)", lambda r: r["pert"] in NUCLEO),
                                    ("el resto", lambda r: r["pert"] not in NUCLEO),
                                    ("todas", lambda r: True)), 1):
        rs = [r for r in pp if sel(r)]
        fi = float(np.mean([0.5 * (r["fiab_L1"] + r["fiab_L2"]) for r in rs]))
        pdd = [r["pendiente_desat"] for r in rs if r["pendiente_desat"] is not None
               and np.isfinite(r["pendiente_desat"])]
        W(f"| {i} | {nom} | {len({r['pert'] for r in rs})} | {len(rs)} | {fi:.4f} | "
          f"**{np.mean([r['concordancia_signo'] for r in rs]):.4f}** | "
          f"**{np.mean([r['spearman'] for r in rs]):+.4f}** | "
          f"**{np.mean([r['pendiente'] for r in rs]):+.4f}** | "
          f"**{np.mean(pdd):+.4f}** | "
          f"**{np.mean([r['R2'] for r in rs]):.5f}** |")
    W("")
    rsN = [r for r in pp if r["pert"] in NUCLEO]
    rsR = [r for r in pp if r["pert"] not in NUCLEO]
    sgN = np.mean([r["concordancia_signo"] for r in rsN])
    sgR = np.mean([r["concordancia_signo"] for r in rsR])
    W("🔴 **Y aqui aparece algo que hay que mirar antes de leer el signo: la**")
    W("**concordancia es MAS ALTA en el subconjunto POCO fiable.**")
    W(f"{sgR:.4f} en las 28 contra **{sgN:.4f}** en las cinco fiables. Si la")
    W("concordancia de signo midiera transferencia de la perturbacion, tendria que ir")
    W("al contrario. Asi que se midio de donde sale, en vez de suponerlo.")
    W("")
    sg = json.load(open(os.path.join(OUTDIR, "signo.json")))
    W("### 5.1 🔴 De donde sale la concordancia de signo: el bloque de ceros")
    W("")
    W("**La sospecha, declarada antes de medir.** Un gen sin expresion en la mitad")
    W("perturbada da `delta = log1p(0) - log1p(ctrl) = -log1p(ctrl)`: negativo, y casi")
    W("el mismo en las dos mitades y en las dos lineas. Esos genes concuerdan en signo")
    W("**por la estructura del espacio**, no por la perturbacion.")
    W("")
    W("**TABLA 5.3 — concordancia de signo por DECIL de `|delta|` del lado predictor**,")
    W("agregada sobre los 15 pares x 33 perts x 8 particiones.")
    W("")
    W("| decil de `|delta|` | n (pares gen-celda) | **concordancia** | nucleo (5 perts) | resto (28 perts) |")
    W("|---|---|---|---|---|")
    for r in sg["por_decil_de_abs_delta"]:
        W(f"| {r['decil']} | {r['n']:,} | **{r['concordancia']:.4f}** | "
          f"{r['concordancia_nucleo']:.4f} | {r['concordancia_resto']:.4f} |")
    W("")
    d3 = sg["por_decil_de_abs_delta"][2]; d4 = sg["por_decil_de_abs_delta"][3]
    d10 = sg["por_decil_de_abs_delta"][9]; d1 = sg["por_decil_de_abs_delta"][0]
    f = sg["frac_de_la_concordancia_que_es_ambos_negativos"]
    W("🔴 **La concordancia NO es monotona en `|delta|`, y eso la delata:**")
    W("")
    W(f"* Los deciles **3 y 4** concentran **{d3['n'] + d4['n']:,}** pares gen-celda")
    W(f"  -- los `|delta|` estan EMPATADOS en bloques enormes, asi que los cortes por")
    W(f"  cuantil colapsan -- y ahi la concordancia sube a **{d3['concordancia']:.4f}**")
    W(f"  y **{d4['concordancia']:.4f}**, con el **resto** de perturbaciones todavia mas")
    W(f"  alto ({d3['concordancia_resto']:.4f} y {d4['concordancia_resto']:.4f}) que el")
    W("  nucleo. Ese es el bloque de ceros estructurales.")
    W(f"* En el decil **10** (los `|delta|` mas grandes, donde esta la senal) la")
    W(f"  concordancia cae a **{d10['concordancia']:.4f}** en conjunto, y solo en el")
    W(f"  **nucleo** llega a **{d10['concordancia_nucleo']:.4f}** contra")
    W(f"  **{d10['concordancia_resto']:.4f}** en el resto.")
    W(f"* En el decil **1** la concordancia esta **por debajo** del azar")
    W(f"  (**{d1['concordancia']:.4f}**).")
    W(f"* De los genes que concuerdan en signo, una media de **{f['media'] * 100:.1f} %**")
    W(f"  (rango {f['min'] * 100:.1f}-{f['max'] * 100:.1f} %) lo hacen siendo **los dos**")
    W("  **negativos**, que es la direccion que impone el cero.")
    W("")
    W("**Asi que el 0.5501 de la TABLA 5.1 no es transferencia de signo**: es en su")
    W("mayor parte el bloque empatado de genes sin expresion. La cifra de signo que")
    W("queda en pie, medida donde hay senal, es la del decil 10:")
    W(f"**{d10['concordancia']:.4f}** en conjunto y **{d10['concordancia_nucleo']:.4f}**")
    W("en las cinco perturbaciones fiables, contra un azar de 0.5000.")
    W("")
    z = json.load(open(os.path.join(OUTDIR, "zeros.json")))
    W("### 5.2 🔴 El mismo control aplicado al PEARSON: ahi NO contamina")
    W("")
    W("La misma pregunta sobre `G_mismo`: ¿su `+0.02` sale del bloque de ceros? Se")
    W("recalculo quitando los genes cuyo `|delta|` cae en el **bloque modal** (el valor")
    W("de `|delta|` mas repetido, que es el del cero) y restringiendo al **decil")
    W("superior** de `|delta|`.")
    W("")
    W("**TABLA 5.4 — `G_mismo` y `G_otro` en los tres espacios de genes.**")
    W(f"`G_otro` aqui usa **un** blanco fijo y declarado (`q = (g+7) mod 33`), no la")
    W("media sobre los 32, para que las tres columnas sean comparables entre si.")
    W("")
    W("| # | subconjunto | filas | frac. de genes en el bloque modal | "
      "**`G_mismo` todos** | **`G_mismo` sin el bloque** | **`G_mismo` decil 10** | "
      "`G_otro` todos | `G_otro` sin el bloque | `G_otro` decil 10 | "
      "**dif. todos** | **dif. decil 10** |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, (nom, d) in enumerate(z["subconjuntos"].items(), 1):
        W(f"| {i} | {nom} | {d['n_filas']:,} | "
          f"{d['frac_genes_en_el_bloque_modal']:.4f} | "
          f"**{d['mismo_todos']['media']:+.4f}** | "
          f"**{d['mismo_sin_bloque']['media']:+.4f}** | "
          f"**{d['mismo_decil_sup']['media']:+.4f}** | "
          f"{d['otro_todos']['media']:+.4f} | {d['otro_sin_bloque']['media']:+.4f} | "
          f"{d['otro_decil_sup']['media']:+.4f} | "
          f"**{d['dif_todos']:+.4f}** | **{d['dif_decil_sup']:+.4f}** |")
    W("")
    t = z["subconjuntos"]["todas"]; n5 = z["subconjuntos"]["nucleo (5)"]
    W("🔴 **El Pearson se comporta al CONTRARIO que el signo, y eso es el control que**")
    W("**le faltaba a `G_mismo`:**")
    W("")
    W(f"* El **{t['frac_genes_en_el_bloque_modal'] * 100:.2f} %** de los genes esta en el")
    W("  bloque modal del cero.")
    W(f"* Quitarlo **no mueve `G_mismo`**: {t['mismo_todos']['media']:+.4f} ->")
    W(f"  {t['mismo_sin_bloque']['media']:+.4f}. En el nucleo,")
    W(f"  {n5['mismo_todos']['media']:+.4f} -> {n5['mismo_sin_bloque']['media']:+.4f}.")
    W(f"* Restringir al **decil superior de `|delta|`** lo SUBE:")
    W(f"  {t['mismo_todos']['media']:+.4f} -> **{t['mismo_decil_sup']['media']:+.4f}** en")
    W(f"  conjunto y {n5['mismo_todos']['media']:+.4f} -> "
      f"**{n5['mismo_decil_sup']['media']:+.4f}** en el nucleo.")
    W(f"* `G_otro` se queda en el ruido en los tres espacios "
      f"({t['otro_todos']['media']:+.4f}, {t['otro_sin_bloque']['media']:+.4f}, "
      f"{t['otro_decil_sup']['media']:+.4f}).")
    W("")
    W("**Asi que las dos medidas de \"orden\" no se comportan igual y no se pueden**")
    W("**resumir juntas**: la concordancia de signo esta dominada por el bloque de")
    W("ceros y baja a casi el azar cuando se la saca de ahi; el Pearson entre contextos")
    W("no depende del bloque y se concentra donde el delta es grande.")
    W("")
    W("🔴 **El Spearman NO comparte ese defecto en la misma medida** y va al revés del")
    W("signo: **+%.4f** en el nucleo contra **+%.4f** en el resto, porque los empates se"
      % (np.mean([r["spearman"] for r in rsN]), np.mean([r["spearman"] for r in rsR])))
    W("promedian y un bloque empatado recibe todos el mismo rango, que no aporta orden.")
    W("")
    W(f"Y el `R²` pasa de **{np.mean([r['R2'] for r in rsN]):.5f}** en el nucleo a")
    W(f"**{np.mean([r['R2'] for r in rsR]):.5f}** en el resto: la magnitud explica")
    W(f"**menos del {max(np.mean([r['R2'] for r in rsN]), np.mean([r['R2'] for r in rsR])) * 100:.1f} %**")
    W("de la varianza en los dos subconjuntos.")
    W("")
    W("---")
    W("")
    return L


if __name__ == "__main__":
    open("/tmp/claude-1000/-home-ubuntu-vcc/52816e9e-e42f-4af3-b6ec-69a5b07f6ce9/scratchpad/mcrep_e.txt",
         "w").write("\n".join(main()))
    print("ok")
