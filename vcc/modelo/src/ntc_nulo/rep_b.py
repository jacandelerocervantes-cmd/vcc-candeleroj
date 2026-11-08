#!/usr/bin/env python
"""Reporte `ntc_nulo.md` -- seccion 0.c (nnz de control y el Poisson aditivo)."""
import json, os
OUTDIR = "/home/ubuntu/vcc/modelo/out/ntc_nulo"


def main():
    p0 = json.load(open(os.path.join(OUTDIR, "p0.json")))
    K = p0["constantes"]
    nn = p0["c_nnz_por_contexto"]
    pc = p0["c_poisson"]
    L = []; W = L.append

    W("### 0.c — El `nnz` de las celulas de control, y el Poisson aditivo")
    W("")
    W("**El tope y de donde sale.** `HARD_CAP_NNZ = 4,750,000,000` entradas")
    W("almacenadas por entrega (`vcc/sizing.py` del CLI instalado, ya confirmado")
    W("leyendo el paquete en `paso1_contrato.md` punto 42). Una entrega lleva")
    W("`300 x 400 x 3 = 360,000` filas, asi que el presupuesto por celula es")
    W(f"`4.75e9 / 360,000 = {K['tope_por_celula']:,.4f}`.")
    W("")
    W("🔴 **Cual de las dos cifras ATA.** El tope es sobre el TOTAL, no por")
    W("celula. `13,194.44` es el presupuesto por celula **si se gasta parejo**, asi")
    W("que la cantidad que ata es la **MEDIA** (total = media x 360,000); el maximo")
    W("de una celula sola puede pasarse de 13,194 sin violar nada mientras el total")
    W("aguante. Las dos columnas de margen van separadas por eso.")
    W("")
    W("**TABLA 0.4 — `nnz` por celula de control, por contexto. 18,400 celulas cada uno.**")
    W("")
    W("| # | ctx | celulas | nnz media | nnz mediana | nnz p99 | nnz max | nnz min | nnz sd | "
      "frac tope (media) | margen (media) | frac tope (max) | margen (max) |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, r in enumerate(nn, 1):
        W(f"| {i} | **{r['contexto']}** | {r['n_celulas']:,} | {r['nnz_media']:,.1f} | "
          f"{r['nnz_mediana']:,.1f} | {r['nnz_p99']:,.1f} | {r['nnz_max']:,} | "
          f"{r['nnz_min']:,} | {r['nnz_sd']:,.1f} | {r['frac_tope_media']:.4f} | "
          f"x{r['margen_sobre_media']:.3f} | {r['frac_tope_max']:.4f} | "
          f"x{r['margen_sobre_max']:.3f} |")
    W("")
    W("**El margen, dicho como lo gasta una entrega.** Una entrega que reusa filas")
    W("de control tal cual (el contrato ya probado) gasta, en esperanza,")
    W("`360,000 x media_nnz` de los tres contextos que use:")
    W("")
    W("| particion | contextos | media nnz/celula de los 3 | nnz total previsto | frac del tope | margen |")
    W("|---|---|---|---|---|---|")
    ref = {}
    for part, letras in (("val", "ABC"), ("final", "DEF")):
        ms = [r["nnz_media"] for r in nn if r["contexto"] in letras]
        mm = sum(ms) / 3.0
        tot = mm * K["filas_entrega"]
        ref[part] = tot
        W(f"| `{part}` | {', '.join(letras)} | {mm:,.2f} | {tot:,.0f} | "
          f"{tot / K['HARD_CAP_NNZ']:.4f} | x{K['HARD_CAP_NNZ'] / tot:.3f} |")
    W("")
    W("🔴 **Cruce independiente.** `paso1_contrato.md` punto 42 MIDIO el nnz real de")
    W("una carga trivial sobre `val` (A/B/C) construida de verdad:")
    W("**2,084,329,279**. La proyeccion de arriba desde la media de `nnz` por celula")
    W(f"da **{ref['val']:,.0f}** -- difieren en")
    W(f"{abs(ref['val'] - 2084329279) / 2084329279 * 100:.3f} %. Son dos caminos")
    W("distintos a la misma cifra (uno construyo la matriz, el otro suma `indptr`).")
    W("")
    W("#### El Poisson aditivo de `cand01_covar._reescala_entero`")
    W("")
    W("**La construccion, declarada.** `_reescala_entero(bloque, rat, rng)` reescala")
    W("un bloque CSR de cuentas enteras por el factor por gen `rat`:")
    W("")
    W("```")
    W("r_h <= 1:  nuevo_h ~ Binomial(x_h, r_h)            (adelgazamiento)")
    W("r_h  > 1:  nuevo_h = x_h + Poisson(x_h * (r_h-1))  (engorde)")
    W("```")
    W("")
    W(f"Se aplico a una muestra de **400 celulas de control** del contexto")
    W(f"**{pc['contexto']}**, sorteadas sin reemplazo con semilla fija")
    W(f"`{pc['semilla_muestra']}`, con las **dos** emisiones de densidad que se")
    W("venian usando, y tres blancos cada una:")
    W("")
    W("* `cand01_covar.construye_respuesta` -- ratios estimados EN CONTEXTO sobre el")
    W("  pool de 18,400 celulas del propio contexto (`r = 1 - rho`, acotado en [0,2]).")
    W("* `gwps_emision.construye_ratios` -- ratios EXTERNOS del perfil de k562_gwps")
    W("  (`r = 1 + s`, `s` del rango normalizado, acotado en [0,2]).")
    W("")
    W(f"El bloque de partida tiene **{pc['nnz_antes']:,}** entradas almacenadas,")
    W(f"**{pc['nnz_antes_por_celula_media']:,.1f}** por celula.")
    W("")
    W("**TABLA 0.5 — los seis brazos, `nnz` antes y despues. ES VALIDEZ, NO ESTILO.**")
    W("")
    W("| # | emision | blanco | genes `r>1` | genes `r<1` | genes `r==1` | `r` min | `r` max | "
      "nnz despues | **delta nnz** | nnz/celula media | nnz max de una celula | frac tope (max) | "
      "**revienta el tope** |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, b in enumerate(pc["brazos"], 1):
        em = ("cand01 (en contexto)" if "cand01" in b["emision"]
              else "gwps (externa)")
        W(f"| {i} | {em} | `{b['blanco']}` | {b['n_genes_r_mayor_1']:,} | "
          f"{b['n_genes_r_menor_1']:,} | {b['n_genes_r_igual_1']:,} | "
          f"{b['r_min']:.4f} | {b['r_max']:.4f} | {b['nnz_despues']:,} | "
          f"**{b['delta_nnz']:+,}** | {b['nnz_despues_por_celula_media']:,.1f} | "
          f"{b['nnz_max_por_celula']:,} | {b['frac_tope_max']:.4f} | "
          f"{'**SI**' if b['revienta_tope'] else 'no'} |")
    W("")
    W("**Y la suma de cuentas, que es otra cantidad** (el tope de `nnz` no la mira,")
    W("pero la clausula de 1,000,000 por celula si):")
    W("")
    W("| # | emision | blanco | suma antes | suma despues | razon |")
    W("|---|---|---|---|---|---|")
    for i, b in enumerate(pc["brazos"], 1):
        em = "cand01" if "cand01" in b["emision"] else "gwps"
        W(f"| {i} | {em} | `{b['blanco']}` | {b['suma_antes']:,.0f} | "
          f"{b['suma_despues']:,.0f} | {b['suma_despues'] / b['suma_antes']:.4f} |")
    W("")
    W("**La cota, y por que es estructural y no un ajuste.** `_reescala_entero`")
    W("opera sobre `b.data` -- las entradas YA GUARDADAS del CSR -- y luego llama")
    W("`eliminate_zeros()`. El brazo Poisson escribe `x + Poisson(...)` **solo donde")
    W("`x` ya era distinto de cero**, asi que no puede crear un no-cero nuevo; el")
    W("brazo binomial si puede destruirlos. De ahi `nnz_despues <= nnz_antes`,")
    W("siempre, sea cual sea `rat`.")
    W("")
    W(f"Medido: el `delta nnz` es **negativo en los 6 de 6 brazos** "
      f"(`{p0['c_poisson']['cota_estructural']['verificado_en_todos_los_brazos']}`), de")
    W(f"**{min(b['delta_nnz'] for b in pc['brazos']):+,}** a")
    W(f"**{max(b['delta_nnz'] for b in pc['brazos']):+,}**. El `nnz` maximo de una")
    W(f"celula despues del reescalado llega a")
    W(f"**{max(b['nnz_max_por_celula'] for b in pc['brazos']):,}**, que es")
    W(f"**{max(b['frac_tope_max'] for b in pc['brazos']):.4f}** del presupuesto por")
    W("celula. **El Poisson aditivo no lo revienta**, y no por poco: no puede")
    W("subirlo en absoluto.")
    W("")
    W("**Resultado de la Parte 0: las tres pasaron. No se disparo ninguna parada.**")
    W("")
    W("---")
    W("")
    return L


if __name__ == "__main__":
    open("/tmp/claude-1000/-home-ubuntu-vcc/52816e9e-e42f-4af3-b6ec-69a5b07f6ce9/scratchpad/rep_b.txt",
         "w").write("\n".join(main()))
    print("ok")
