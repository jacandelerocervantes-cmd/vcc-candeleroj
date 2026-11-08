#!/usr/bin/env python
"""Reporte `magnitud_centro.md` -- TRAMPA 1: fiabilidad y desatenuacion."""
import json, os
import numpy as np
OUTDIR = "/home/ubuntu/vcc/modelo/out/magnitud_centro"
ALTA = 0.15          # umbral de fiabilidad para el subconjunto de sensibilidad


def main():
    d = json.load(open(os.path.join(OUTDIR, "calc.json")))
    panel = d["panel"]
    L = []; W = L.append
    fl = d["fiabilidad"]["por_linea"]
    fp = d["fiabilidad"]["por_linea_y_pert"]
    lineas = sorted(fl)

    W("## 2. TRAMPA 1 — la atenuacion, y las tres columnas")
    W("")
    W("### 2.1 🔴 La fiabilidad DENTRO de linea: la cota de todo lo demas")
    W("")
    W("`fiab[L,g] = corr(delta_g,L,A , delta_g,L,B)` sobre los 33,992 genes, con las")
    W("dos mitades de celulas disjuntas y los dos medios controles disjuntos. Es el")
    W("techo de cualquier correlacion entre contextos: la desatenuada divide por")
    W("`sqrt(fiab_L1 * fiab_L2)`, asi que **si la fiabilidad es baja, la desatenuada**")
    W("**es inestable**, y si es <= 0 la raiz no existe.")
    W("")
    W("**TABLA 2.1 — la fiabilidad por linea, sobre las 8 x 33 = 264 celdas de cada una.**")
    W("")
    W("| # | linea | media | mediana | sd | min | max | **celdas con `fiab <= 0`** |")
    W("|---|---|---|---|---|---|---|---|")
    for i, l in enumerate(lineas, 1):
        v = fl[l]
        W(f"| {i} | **{l}** | {v['media']:.4f} | {v['mediana']:.4f} | {v['sd']:.4f} | "
          f"{v['min']:+.4f} | {v['max']:.4f} | "
          f"**{v['n_no_positivas']}/{v['n_celdas']}** |")
    g = d["fiabilidad"]["global"]
    W(f"| — | **las 6** | **{g['media']:.4f}** | — | — | {g['min']:+.4f} | "
      f"{g['max']:.4f} | **{g['n_no_positivas']}/{g['n_celdas']}** |")
    W("")
    W(f"🔴 **La fiabilidad media es {g['media']:.4f}**: el delta de una mitad de 200")
    W("celulas apenas se reproduce contra la otra mitad. Y")
    W(f"**{g['n_no_positivas']} de {g['n_celdas']} celdas (linea, particion, blanco)")
    W("tienen `fiab <= 0`**, donde la desatenuacion no esta definida.")
    W("")
    W("**TABLA 2.2 — la fiabilidad por PERTURBACION, media sobre las 8 particiones.**")
    W("Ordenada por la media sobre las seis lineas. Esta es la tabla que acota todas")
    W("las demas cifras del reporte.")
    W("")
    W("| # | pert | " + " | ".join(f"**{l}**" for l in lineas) +
      " | **media 6 lineas** | min | max | n lineas con `<= 0` |")
    W("|---|---|" + "---|" * (len(lineas) + 4))
    filas = []
    for p in panel:
        v = [fp[l][p]["media"] for l in lineas]
        filas.append((float(np.mean(v)), p, v))
    for i, (m, p, v) in enumerate(sorted(filas, reverse=True), 1):
        W(f"| {i} | `{p}` | " + " | ".join(f"{x:+.4f}" for x in v) +
          f" | **{m:+.4f}** | {min(v):+.4f} | {max(v):+.4f} | "
          f"{sum(1 for x in v if x <= 0)} |")
    W("")
    altas = [p for m, p, v in filas if m >= ALTA]
    bajas = [p for m, p, v in filas if m < ALTA]
    W(f"🔴 **La fiabilidad esta CONCENTRADA en {len(altas)} de las 33 perturbaciones.**")
    W(f"Con media sobre las seis lineas `>= {ALTA}`:")
    W(f"{', '.join('`%s`' % x for x in altas)} -- que son el nucleo del receptor de")
    W("interferon y la via JAK-STAT. Las otras")
    W(f"**{len(bajas)}** quedan por debajo, y")
    W(f"**{sum(1 for m, p, v in filas if min(v) < 0.05)} de 33** caen por debajo de")
    W("0.05 en **alguna** linea. El rango de la media por perturbacion es")
    W(f"**[{min(m for m, _, _ in filas):+.4f}, {max(m for m, _, _ in filas):+.4f}]**.")
    W("")
    W("Eso significa que casi todo lo que sigue se mide sobre deltas que **apenas se**")
    W("**reproducen consigo mismos**. No se corrige: se reporta, y §3.3 da la misma")
    W("medicion restringida al subconjunto fiable para que se vea cuanto cambia.")
    W("")
    W("---")
    W("")
    return L, altas, bajas


if __name__ == "__main__":
    t = main()
    open("/tmp/claude-1000/-home-ubuntu-vcc/52816e9e-e42f-4af3-b6ec-69a5b07f6ce9/scratchpad/mcrep_b.txt",
         "w").write("\n".join(t[0]))
    print("ok, altas:", t[1])
