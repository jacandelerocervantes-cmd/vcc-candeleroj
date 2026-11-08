#!/usr/bin/env python
"""Reporte `ntc_nulo.md` -- PARTE 1: el nulo, por contexto."""
import json, os
import numpy as np
OUTDIR = "/home/ubuntu/vcc/modelo/out/ntc_nulo"
BRAZOS = ["completo", "n400_400", "n400_ref_completa"]
ETQ = {"completo": "`completo` 9,200 v 9,200",
       "n400_400": "`n400_400` 400 v 400",
       "n400_ref_completa": "`n400_ref_completa` 400 v 9,200"}


def main():
    p1 = json.load(open(os.path.join(OUTDIR, "p1.json")))
    L = []; W = L.append
    cors = p1["corridas"]
    parts = p1["particiones"]
    ctxs = sorted({r["contexto"] for r in cors})

    W("## PARTE 1 — El nulo, por contexto")
    W("")
    W("### 1.1 La construccion, declarada antes de medir")
    W("")
    W("* **La unidad es la GUIA, no la celula.** Las 400 celulas de una guia van")
    W("  enteras a un lado. Partir celulas DENTRO de una guia subestimaria el nulo")
    W("  porque las dos mitades compartirian el efecto de esa guia; aqui no se hace")
    W("  nunca, y se verifica particion por particion.")
    W("* **Las 46 guias se parten en 23 contra 23.** `rng.permutation(46)` sobre las")
    W("  46 categorias de `ntc_id` ORDENADAS alfabeticamente; los primeros 23 son la")
    W("  mitad de **REFERENCIA**, los ultimos 23 la mitad **BLANCO**.")
    W("* **Semilla**: `20274000 + 100*idx_contexto + j`, con")
    W("  `idx` A=0, B=1, C=2, D=3, E=4, F=5, y `j = 0..7` -- **ocho particiones por**")
    W("  **contexto**, 48 en total.")
    W("* **La puerta es la del preset `vcc2026`**, la misma que usa el marcador, con")
    W("  cada parametro leido del preset (tabla de la seccion 0).")
    W("")
    W("**Los tres brazos de `n`, declarados aparte y no mezclados:**")
    W("")
    W("| brazo | `n` referencia | `n` blanco | que es |")
    W("|---|---|---|---|")
    W("| `completo` | 9,200 | 9,200 | las 23 guias enteras de cada lado |")
    W("| `n400_400` | 400 | 400 | el `n` POR PERTURBACION del reto, a los dos lados. "
      "**Es el brazo que pide el encargo.** |")
    W("| `n400_ref_completa` | 9,200 | 400 | 🔴 **AÑADIDO y declarado.** El reto no "
      "compara 400 contra 400: compara las 400 celulas de una perturbacion contra el "
      "pool de control del contexto (`manifest.json`: `138,400 = 300x400 + 18,400`). "
      "Este brazo es el que se parece al regimen de la entrega. |")
    W("")
    W("El submuestreo a 400 es **sin reemplazo** sobre las 9,200 celulas de la mitad,")
    W("con el mismo `rng` de la particion y en orden declarado (referencia primero,")
    W("blanco despues). **No parte guias**: las guias siguen enteras de un lado, y el")
    W("submuestreo solo adelgaza celulas dentro de la mitad que ya les toco.")
    W("")
    W("🔴 **La asimetria de la puerta, declarada como propiedad, no como ajuste.** El")
    W("filtro de genes del preset (`filter_gene_min_cpm_cell = 5.0`) lo decide **el")
    W("grupo de REFERENCIA solo** (`_apply_cpm_filter`, #351), y BH se recalcula sobre")
    W("los genes que quedan. Asi que el conjunto probado `m` depende de que mitad hace")
    W("de referencia. Por eso `m` varia entre particiones y se reporta en cada corrida.")
    W("")

    # -------------------------------------------------- 1.2 particiones
    W("### 1.2 Las 48 particiones: ninguna parte una guia")
    W("")
    W("**TABLA 1.1 — las 48 particiones, verificadas una por una.** `guia partida` es")
    W("falso en las 48; los conjuntos de guias son disjuntos y su union son las 46; y")
    W("las 46 guias tienen 400 celulas cada una a los dos lados.")
    W("")
    W("| # | ctx | j | semilla | guias ref | guias blanco | disjuntas | union = 46 | "
      "celulas ref | celulas blanco | 400/guia ref | 400/guia blanco | **guia partida** |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, r in enumerate(parts, 1):
        W(f"| {i} | **{r['contexto']}** | {r['j']} | `{r['semilla']}` | "
          f"{r['n_guias_ref']} | {r['n_guias_tgt']} | "
          f"{'si' if r['disjuntas'] else '**no**'} | "
          f"{'si' if r['union_46'] else '**no**'} | {r['celulas_ref']:,} | "
          f"{r['celulas_tgt']:,} | {'si' if r['celulas_por_guia_ref_ok'] else '**no**'} | "
          f"{'si' if r['celulas_por_guia_tgt_ok'] else '**no**'} | "
          f"**{'SI' if r['guia_partida'] else 'no'}** |")
    W("")
    W("Ejemplo de una particion escrita entera, para que se pueda reproducir")
    W(f"(contexto {parts[0]['contexto']}, j={parts[0]['j']}, semilla "
      f"`{parts[0]['semilla']}`):")
    W("")
    W(f"* mitad REFERENCIA (23 guias): {', '.join('`%s`' % g for g in parts[0]['guias_ref'])}")
    W(f"* mitad BLANCO (23 guias): {', '.join('`%s`' % g for g in parts[0]['guias_tgt'])}")
    W("")

    # -------------------------------------------------- 1.3 resumen
    W("### 1.3 🔴 Cuantos genes llama significativos, y su distribucion sobre particiones")
    W("")
    W("**TABLA 1.2 — resumen por contexto y brazo, sobre las 8 particiones.**")
    W("`rango` es `[min, max]` de `n_sig` sobre las 8 particiones de ese contexto.")
    W("")
    W("| # | ctx | brazo | corridas | **`n_sig` media** | **`n_sig` sd** | mediana | "
      "**rango [min, max]** | corridas con `n_sig>0` | `m` probados [min, max] | "
      "`p_crit = 0.05/m` | `p_min` mediana | genes llamados (union) |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    i = 0
    resumen = {}
    for c in ctxs:
        for b in BRAZOS:
            rs = [r for r in cors if r["contexto"] == c and r["brazo"] == b]
            if not rs:
                continue
            i += 1
            ns = np.array([r["n_sig"] for r in rs])
            m = np.array([r["n_probados"] for r in rs])
            pm = np.array([r["p_min"] for r in rs], dtype=float)
            un = sorted({g for r in rs for g in (r["genes_sig"] or [])})
            resumen[(c, b)] = (ns, m, un)
            W(f"| {i} | **{c}** | {ETQ[b]} | {len(rs)} | **{ns.mean():.3f}** | "
              f"**{ns.std(ddof=1):.3f}** | {np.median(ns):.1f} | "
              f"**[{ns.min()}, {ns.max()}]** | {int((ns > 0).sum())}/{len(rs)} | "
              f"[{m.min():,}, {m.max():,}] | "
              f"[{0.05 / m.max():.2e}, {0.05 / m.min():.2e}] | "
              f"{np.median(pm):.3e} | {len(un)} |")
    W("")
    W("**El rango por contexto, en una sola linea por brazo** (sobre los 6 contextos")
    W("x 8 particiones de cada brazo):")
    W("")
    W("| brazo | corridas | `n_sig` media | `n_sig` sd | rango global [min, max] | "
      "corridas con `n_sig>0` | suma de `n_sig` | genes distintos llamados |")
    W("|---|---|---|---|---|---|---|---|")
    for b in BRAZOS:
        rs = [r for r in cors if r["brazo"] == b]
        if not rs:
            continue
        ns = np.array([r["n_sig"] for r in rs])
        un = sorted({g for r in rs for g in (r["genes_sig"] or [])})
        W(f"| {ETQ[b]} | {len(rs)} | {ns.mean():.4f} | {ns.std(ddof=1):.4f} | "
          f"**[{ns.min()}, {ns.max()}]** | {int((ns > 0).sum())}/{len(rs)} | "
          f"{int(ns.sum())} | {len(un)} |")
    W("")
    W("**Los genes que el nulo llego a llamar** (con en cuantas corridas, de las de")
    W("ese brazo):")
    W("")
    W("| brazo | gen | corridas que lo llaman | ctx y j |")
    W("|---|---|---|---|")
    hay = False
    for b in BRAZOS:
        rs = [r for r in cors if r["brazo"] == b]
        cnt = {}
        for r in rs:
            for g in (r["genes_sig"] or []):
                cnt.setdefault(g, []).append(f"{r['contexto']}{r['j']}")
        for g, w in sorted(cnt.items(), key=lambda kv: -len(kv[1])):
            hay = True
            W(f"| {ETQ[b]} | `{g}` | {len(w)} | {', '.join(w)} |")
    if not hay:
        W("| — | — | 0 | ninguna corrida de ningun brazo llamo ningun gen |")
    W("")
    # recurrencia: el nulo NO es intercambiable entre genes
    rec = {}
    for b in BRAZOS:
        cnt = {}
        for r in [x for x in cors if x["brazo"] == b]:
            for g in (r["genes_sig"] or []):
                cnt.setdefault(g, []).append(f"{r['contexto']}{r['j']}")
        rec[b] = {g: w for g, w in cnt.items() if len(w) > 1}
    if any(rec.values()):
        W("🔴 **El nulo NO es intercambiable entre genes: hay recurrencia.** Si los")
        W("genes llamados fueran falsos positivos independientes, con tan pocas")
        W("llamadas en total la probabilidad de que un mismo gen salga varias veces")
        W("seria minuscula. Los que repiten:")
        W("")
        W("| brazo | gen | corridas que lo llaman | de | ctx y j |")
        W("|---|---|---|---|---|")
        for b in BRAZOS:
            n_b = len([x for x in cors if x["brazo"] == b])
            for g, w in sorted(rec[b].items(), key=lambda kv: -len(kv[1])):
                W(f"| {ETQ[b]} | **`{g}`** | **{len(w)}** | {n_b} | {', '.join(w)} |")
        W("")
        W("Esa recurrencia es la huella de la estructura de guia que la seccion 2.3")
        W("mide por separado: un gen cuyo nivel depende de QUE guia lleva la celula")
        W("sale llamado cada vez que la particion separa las guias que lo mueven.")
        W("Cuantifica tambien el limite de la propia construccion: las 8 particiones")
        W("de un contexto NO son 8 muestras independientes del nulo, porque comparten")
        W("las mismas 46 guias.")
        W("")
    return L, cors, parts, ctxs


if __name__ == "__main__":
    L, *_ = main()
    open("/tmp/claude-1000/-home-ubuntu-vcc/52816e9e-e42f-4af3-b6ec-69a5b07f6ce9/scratchpad/rep_c.txt",
         "w").write("\n".join(L))
    print("ok, lineas:", len(L))
