#!/usr/bin/env python
"""Reporte `magnitud_centro.md` -- TRAMPA 2: MISMO / OTRO blanco, y por par."""
import json, os
import numpy as np
OUTDIR = "/home/ubuntu/vcc/modelo/out/magnitud_centro"
ALTA = 0.15
NUCLEO = ["IFNAR1", "IFNAR2", "JAK1", "STAT1", "TYK2"]


def ra(v):
    v = [x for x in v if x is not None and np.isfinite(x)]
    return (min(v), max(v), float(np.mean(v))) if v else (None, None, None)


def main():
    d = json.load(open(os.path.join(OUTDIR, "calc.json")))
    pares = d["por_par"]; pp = d["por_par_y_pert"]
    L = []; W = L.append

    W("## 3. TRAMPA 2 — el eje compartido: MISMO blanco contra OTRO blanco")
    W("")
    W("**Por que sin este control la cifra no significa nada.** Las 33 perturbaciones")
    W("son todas de la via del interferon, asi que \"compartido entre contextos\" se")
    W("confunde con \"compartido entre perturbaciones\". Es la forma `G_mismo` /")
    W("`G_otro` de `gwps_jiang_datos.md` §1.5, aplicada entre lineas:")
    W("")
    W("| | que es |")
    W("|---|---|")
    W("| **`G_mismo`** | MISMO blanco, linea distinta -> **el centro que se busca** |")
    W("| **`G_otro`** | OTRO blanco, linea distinta -> **la respuesta comun al IFNB** "
      "(media sobre los 32 `q != g`) |")
    W("| **diferencia** | 🔴 la parte **especifica del blanco** que cruza contextos |")
    W("")
    W("Las dos se calculan con **mitades cruzadas** y simetrizadas: media de")
    W("`corr(A de L1, B de L2)` y `corr(B de L1, A de L2)`. Las dos usan celulas")
    W("disjuntas y medios controles disjuntos, asi que el ruido no se comparte en")
    W("ninguna.")
    W("")
    W("### 3.1 Por par de lineas: las tres columnas, cruda y desatenuada")
    W("")
    W("**TABLA 3.1 — los 15 pares. `n/d` = celdas donde `fiab <= 0` en alguna de las")
    W("dos lineas y la raiz no existe.**")
    W("")
    W("| # | par | `fiab` L1 | `fiab` L2 | **`G_mismo` cruda** | **`G_otro` cruda** | "
      "**dif. cruda** | **`G_mismo` desat.** | **`G_otro` desat.** | **dif. desat.** | "
      "`G_otro` sd sobre `q` | `G_otro` [min, max] sobre `q` | celdas `n/d` | "
      "**celdas desat. > 1** |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, p in enumerate(pares, 1):
        W(f"| {i} | **{p['par']}** | {p['fiab_L1']['media']:.4f} | "
          f"{p['fiab_L2']['media']:.4f} | "
          f"**{p['mismo_cruda']['media']:+.4f}** | {p['otro_cruda']['media']:+.4f} | "
          f"**{p['dif_cruda']['media']:+.4f}** | "
          f"**{p['mismo_desat']['media']:+.4f}** | {p['otro_desat']['media']:+.4f} | "
          f"**{p['dif_desat']['media']:+.4f}** | {p['otro_sd_media']:.4f} | "
          f"[{p['otro_min_global']:+.4f}, {p['otro_max_global']:+.4f}] | "
          f"{p['n_desat_nd']}/264 | **{p['n_desat_mayor_1']}** |")
    W("")
    W("**Los rangos sobre los 15 pares:**")
    W("")
    W("| # | cantidad | **[min, max]** | media |")
    W("|---|---|---|---|")
    campos = [("`G_mismo` cruda", "mismo_cruda"), ("`G_otro` cruda", "otro_cruda"),
              ("**diferencia cruda**", "dif_cruda"),
              ("`G_mismo` desatenuada", "mismo_desat"),
              ("`G_otro` desatenuada", "otro_desat"),
              ("**diferencia desatenuada**", "dif_desat")]
    for i, (lab, k) in enumerate(campos, 1):
        mn, mx, me = ra([p[k]["media"] for p in pares])
        W(f"| {i} | {lab} | **[{mn:+.4f}, {mx:+.4f}]** | {me:+.4f} |")
    W("")
    nd = sum(p["n_desat_nd"] for p in pares)
    g1 = sum(p["n_desat_mayor_1"] for p in pares)
    W("🔴 **Lo que estas filas dicen, como aritmetica:**")
    W("")
    mn, mx, _ = ra([p["mismo_cruda"]["media"] for p in pares])
    W(f"* La correlacion CRUDA entre contextos con el MISMO blanco va de")
    W(f"  **{mn:+.4f} a {mx:+.4f}**. Es casi nada en valor absoluto.")
    mno, mxo, _ = ra([p["otro_cruda"]["media"] for p in pares])
    W(f"* La de OTRO blanco va de **{mno:+.4f} a {mxo:+.4f}**, asi que el eje comun al")
    W("  IFNB **no es cero** -- la TRAMPA 2 es real -- pero es **menor** que la del")
    W("  mismo blanco en los 15 pares de 15.")
    mnd, mxd, _ = ra([p["dif_cruda"]["media"] for p in pares])
    W(f"* La DIFERENCIA, que es la parte especifica del blanco que cruza contextos, va")
    W(f"  de **{mnd:+.4f} a {mxd:+.4f}** y es **positiva en los 15 pares**.")
    mnx, mxx, _ = ra([p["mismo_desat"]["media"] for p in pares])
    W(f"* Desatenuada, `G_mismo` sube a **[{mnx:+.4f}, {mxx:+.4f}]**: la atenuacion era")
    W("  de un orden de magnitud, que es lo que cabia esperar de una fiabilidad de")
    W("  ~0.09.")
    W(f"* 🔴 **{nd} de {264 * len(pares):,} celdas salen `n/d`** porque `fiab <= 0` en")
    W(f"  una de las dos lineas, y **{g1} celdas dan una desatenuada > 1**. "
      f"**No se recortan** (lo prohibe el encargo y seria esconder la inestabilidad):")
    W("  una correlacion desatenuada por encima de 1 es la senal de que el denominador")
    W("  esta estimado con mas ruido que el numerador.")
    W("")
    W("### 3.2 Por par y por perturbacion")
    W("")
    W(f"**TABLA 3.2 — las {len(pp)} filas (15 pares x 33 perturbaciones), media sobre")
    W("las 8 particiones.** Ordenada por `G_mismo` desatenuada descendente.")
    W("")
    W("| # | par | pert | `fiab` L1 | `fiab` L2 | `G_mismo` cruda (sd) | "
      "`G_mismo` desat. | `G_otro` cruda | `G_otro` desat. | dif. cruda | dif. desat. |")
    W("|---|---|---|---|---|---|---|---|---|---|---|")
    ordenado = sorted(pp, key=lambda r: -(r["mismo_desat"] if r["mismo_desat"]
                                          is not None and np.isfinite(r["mismo_desat"])
                                          else -9e9))
    for i, r in enumerate(ordenado, 1):
        md = (f"{r['mismo_desat']:+.4f}" if r["mismo_desat"] is not None
              and np.isfinite(r["mismo_desat"]) else "**n/d**")
        od = (f"{r['otro_desat']:+.4f}" if r["otro_desat"] is not None
              and np.isfinite(r["otro_desat"]) else "n/d")
        dd = (f"{r['dif_desat']:+.4f}" if r["dif_desat"] is not None
              and np.isfinite(r["dif_desat"]) else "n/d")
        W(f"| {i} | {r['par']} | `{r['pert']}` | {r['fiab_L1']:+.4f} | "
          f"{r['fiab_L2']:+.4f} | {r['mismo_cruda']:+.4f} ({r['mismo_cruda_sd']:.4f}) | "
          f"{md} | {r['otro_cruda']:+.4f} | {od} | {r['dif_cruda']:+.4f} | {dd} |")
    W("")
    W("### 3.3 🔴 La misma medicion, partida por FIABILIDAD")
    W("")
    W(f"El subconjunto fiable son las **{len(NUCLEO)}** perturbaciones con fiabilidad")
    W(f"media sobre las seis lineas `>= {ALTA}` (TABLA 2.2):")
    W(f"{', '.join('`%s`' % x for x in NUCLEO)}. Las otras 28 van aparte. **No se**")
    W("**elige entre las dos: las dos estan.**")
    W("")
    W("**TABLA 3.3 — las tres columnas, por subconjunto de fiabilidad.**")
    W("")
    W("| # | subconjunto | n perts | n filas | `fiab` media | **`G_mismo` cruda** | "
      "**`G_otro` cruda** | **dif. cruda** | **`G_mismo` desat.** | **dif. desat.** | "
      "filas `n/d` |")
    W("|---|---|---|---|---|---|---|---|---|---|---|")
    for i, (nom, sel) in enumerate((("fiable (`>= %.2f`)" % ALTA,
                                     lambda r: r["pert"] in NUCLEO),
                                    ("el resto", lambda r: r["pert"] not in NUCLEO),
                                    ("todas", lambda r: True)), 1):
        rs = [r for r in pp if sel(r)]
        fi = float(np.mean([0.5 * (r["fiab_L1"] + r["fiab_L2"]) for r in rs]))
        mc = float(np.mean([r["mismo_cruda"] for r in rs]))
        oc = float(np.mean([r["otro_cruda"] for r in rs]))
        dc = float(np.mean([r["dif_cruda"] for r in rs]))
        mdv = [r["mismo_desat"] for r in rs if r["mismo_desat"] is not None
               and np.isfinite(r["mismo_desat"])]
        ddv = [r["dif_desat"] for r in rs if r["dif_desat"] is not None
               and np.isfinite(r["dif_desat"])]
        W(f"| {i} | {nom} | {len({r['pert'] for r in rs})} | {len(rs)} | {fi:.4f} | "
          f"**{mc:+.4f}** | {oc:+.4f} | **{dc:+.4f}** | "
          f"**{np.mean(mdv):+.4f}** | **{np.mean(ddv):+.4f}** | "
          f"{len(rs) - len(mdv)} |")
    W("")
    rsN = [r for r in pp if r["pert"] in NUCLEO]
    rsR = [r for r in pp if r["pert"] not in NUCLEO]
    W("🔴 **La lectura es aritmetica, no interpretacion:** en el subconjunto fiable la")
    W(f"`G_mismo` cruda es **{np.mean([r['mismo_cruda'] for r in rsN]):+.4f}** y en el")
    W(f"resto **{np.mean([r['mismo_cruda'] for r in rsR]):+.4f}**, una razon de")
    W(f"**x{np.mean([r['mismo_cruda'] for r in rsN]) / np.mean([r['mismo_cruda'] for r in rsR]):.2f}**.")
    W("Casi toda la correlacion entre contextos del reporte sale de las cinco")
    W("perturbaciones del nucleo del receptor de interferon; las otras 28 aportan una")
    W("correlacion cruda de centesimas sobre una fiabilidad de centesimas.")
    W("")
    W("---")
    W("")
    return L


if __name__ == "__main__":
    open("/tmp/claude-1000/-home-ubuntu-vcc/52816e9e-e42f-4af3-b6ec-69a5b07f6ce9/scratchpad/mcrep_c.txt",
         "w").write("\n".join(main()))
    print("ok")
