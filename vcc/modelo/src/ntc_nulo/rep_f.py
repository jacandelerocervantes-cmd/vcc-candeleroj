#!/usr/bin/env python
"""Reporte `ntc_nulo.md` -- PARTE 3: contra el regimen de Jiang."""
import json, os
import numpy as np
OUTDIR = "/home/ubuntu/vcc/modelo/out/ntc_nulo"


def main():
    p3 = json.load(open(os.path.join(OUTDIR, "p3.json")))
    L = []; W = L.append

    W("## PARTE 3 — Contra el regimen de Jiang")
    W("")
    W("**La referencia que se compara.** `control_permutado.md` §7 MIDIO, en las seis")
    W("lineas de Jiang, el numero de perturbaciones con DE REAL vacio: **A549 15,")
    W("BXPC3 0, HAP1 14, HT29 6, K562 8, MCF7 0**, de 33 -- el rango 0 a 15 de 33 que")
    W("cita el encargo. En el reto, `n = 400` es fijo por `manifest.json`")
    W("(`cells_per_pert: 400`) y el knockdown es >=80 % por criterio de seleccion.")
    W("")
    W("**La construccion, declarada.** Las dos cantidades se miden con la MISMA")
    W("definicion en los dos lados:")
    W("")
    W("* **profundidad por celula** = suma de cuentas de la fila (UMIs). En Jiang se")
    W("  cruzo contra `obs/nCount_RNA` (coincide fila por fila en las seis lineas);")
    W("  en el reto se suma del propio CSR.")
    W("* **`nnz` por celula** = `diff(indptr)`. En Jiang se cruzo contra")
    W(f"  `obs/nFeature_RNA` (identico: `{p3['jiang_nfeature_igual_indptr']}`).")
    W("* **panel probable** = nº de genes cuyo **CPM medio por celula en el pool de")
    W("  control** supera `5.0` -- la MISMA cantidad con la que `_apply_cpm_filter`")
    W("  del preset decide que genes entran en BH. Es el `m` de BH.")
    W("")
    W("De Jiang se leyeron **solo las celulas `NT`** de cada linea. El lado")
    W("perturbado no se toco.")
    W("")
    W("**TABLA 3.1 — profundidad y `nnz` de los controles: reto contra Jiang, por linea.**")
    W("")
    W("| # | fuente | unidad | celulas de control | genes del panel | prof. mediana | "
      "prof. media | prof. p1 | prof. p99 | prof. sd | `nnz` mediana | `nnz` media | "
      "`nnz` p99 | `nnz` max | `nnz` sd | **panel probable (CPM>5)** | frac del panel | "
      "`nnz`/panel |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    i = 0
    for k, fuente in (("reto", "reto 2026"), ("jiang", "Jiang IFNB")):
        for r in p3[k]:
            i += 1
            n = r.get("n_celulas_control") or r.get("n_celulas_control_NT")
            W(f"| {i} | {fuente} | **{r['unidad']}** | {n:,} | {r['n_genes_panel']:,} | "
              f"{r['profundidad_mediana']:,.0f} | {r['profundidad_media']:,.0f} | "
              f"{r['profundidad_p1']:,.0f} | {r['profundidad_p99']:,.0f} | "
              f"{r['profundidad_sd']:,.0f} | {r['nnz_mediana']:,.0f} | "
              f"{r['nnz_media']:,.0f} | {r['nnz_p99']:,.0f} | {r['nnz_max']:,.0f} | "
              f"{r['nnz_sd']:,.0f} | **{r['n_genes_cpm_mayor_5']:,}** | "
              f"{r['frac_panel_probable']:.4f} | {r['nnz_sobre_panel']:.4f} |")
    W("")
    rt, ji = p3["reto"], p3["jiang"]
    pr_r = np.array([r["profundidad_mediana"] for r in rt])
    pr_j = np.array([r["profundidad_mediana"] for r in ji])
    nz_r = np.array([r["nnz_mediana"] for r in rt])
    nz_j = np.array([r["nnz_mediana"] for r in ji])
    pa_r = np.array([r["n_genes_cpm_mayor_5"] for r in rt])
    pa_j = np.array([r["n_genes_cpm_mayor_5"] for r in ji])
    nc_r = np.array([r["n_celulas_control"] for r in rt])
    nc_j = np.array([r["n_celulas_control_NT"] for r in ji])

    W("### 3.1 🔴 Las cantidades, lado a lado. NO se concluye si la verdad vacia transfiere.")
    W("")
    W("**TABLA 3.2 — las razones entre los dos regimenes, calculadas de la TABLA 3.1.**")
    W("")
    W("| # | cantidad | reto (6 contextos) [min, max] | Jiang (6 lineas) [min, max] | "
      "razon reto/Jiang [min, max] |")
    W("|---|---|---|---|---|")
    filas = [
        ("celulas de control por unidad", nc_r, nc_j, "{:,.0f}"),
        ("profundidad mediana (UMIs/celula)", pr_r, pr_j, "{:,.0f}"),
        ("`nnz` mediana por celula", nz_r, nz_j, "{:,.0f}"),
        ("genes del panel entregado", np.array([r["n_genes_panel"] for r in rt]),
         np.array([r["n_genes_panel"] for r in ji]), "{:,.0f}"),
        ("**panel probable (el `m` de BH)**", pa_r, pa_j, "{:,.0f}"),
    ]
    for i, (lab, a, b, fmt) in enumerate(filas, 1):
        W(f"| {i} | {lab} | [{fmt.format(a.min())}, {fmt.format(a.max())}] | "
          f"[{fmt.format(b.min())}, {fmt.format(b.max())}] | "
          f"[x{a.min() / b.max():.2f}, x{a.max() / b.min():.2f}] |")
    W("")
    W("**Y la heterogeneidad DENTRO de cada lado**, que es la cantidad que el encargo")
    W("nombra como causa del rango 0-15 en Jiang:")
    W("")
    W("| # | cantidad | reto: max/min entre contextos | Jiang: max/min entre lineas |")
    W("|---|---|---|---|")
    het = [("celulas de control", nc_r, nc_j),
           ("profundidad mediana", pr_r, pr_j),
           ("`nnz` mediana", nz_r, nz_j),
           ("panel probable", pa_r, pa_j)]
    for i, (lab, a, b) in enumerate(het, 1):
        W(f"| {i} | {lab} | **x{a.max() / a.min():.3f}** | **x{b.max() / b.min():.3f}** |")
    W("")
    W("**Las cifras que quedan en la mesa, sin lectura:**")
    W("")
    W(f"1. El `n` de control por unidad: el reto lleva **18,400 en los seis**")
    W(f"   contextos, exactamente iguales. Jiang lleva de **{nc_j.min():,}** (A549) a")
    W(f"   **{nc_j.max():,}** (BXPC3) celulas `NT`, una dispersion de")
    W(f"   **x{nc_j.max() / nc_j.min():.2f}** entre lineas. El pool de control del reto")
    W(f"   es de **x{nc_r.min() / nc_j.max():.1f} a x{nc_r.max() / nc_j.min():.1f}** el")
    W("   de una linea de Jiang.")
    W(f"2. La profundidad por celula: el reto esta en")
    W(f"   **{pr_r.min():,.0f}-{pr_r.max():,.0f}** UMIs de mediana, Jiang en")
    W(f"   **{pr_j.min():,.0f}-{pr_j.max():,.0f}**. El reto es")
    W(f"   **x{pr_r.min() / pr_j.max():.2f} a x{pr_r.max() / pr_j.min():.2f}** mas profundo.")
    W(f"3. El `nnz` por celula: **{nz_r.min():,.0f}-{nz_r.max():,.0f}** en el reto")
    W(f"   contra **{nz_j.min():,.0f}-{nz_j.max():,.0f}** en Jiang, razon")
    W(f"   **x{nz_r.min() / nz_j.max():.2f} a x{nz_r.max() / nz_j.min():.2f}**.")
    W(f"4. 🔴 **El panel probable va al reves.** Aunque el reto es mas profundo, su")
    W(f"   panel entregado tiene 18,533 genes y Jiang 34,025, asi que los genes que")
    W(f"   pasan el filtro CPM>5 -- el `m` sobre el que corre BH -- son")
    W(f"   **{pa_r.min():,}-{pa_r.max():,}** en el reto y **{pa_j.min():,}-{pa_j.max():,}**")
    W(f"   en Jiang: **x{pa_r.min() / pa_j.max():.3f} a x{pa_r.max() / pa_j.min():.3f}**.")
    W(f"   Jiang prueba MAS genes, asi que su escalon de BH de rango 1 es MAS exigente")
    W(f"   (`0.05/{pa_j.max():,} = {0.05 / pa_j.max():.3e}` contra")
    W(f"   `0.05/{pa_r.min():,} = {0.05 / pa_r.min():.3e}`).")
    W(f"5. La heterogeneidad dentro del lado: en el reto la profundidad mediana varia")
    W(f"   **x{pr_r.max() / pr_r.min():.3f}** entre contextos y el `n` de control")
    W(f"   **x1.000**; en Jiang la profundidad varia **x{pr_j.max() / pr_j.min():.3f}** y")
    W(f"   el `n` **x{nc_j.max() / nc_j.min():.3f}**.")
    W("")
    W("🔴 **No se declara si la verdad vacia de Jiang transfiere ni si no**. Las")
    W("cantidades de arriba son lo que se midio. Lo que falta para cerrar esa")
    W("pregunta -- los tamanos de efecto reales de las 300 perturbaciones del reto --")
    W("esta en el lado perturbado, que no existe en 2026; la seccion 2.4 lo detalla.")
    W("")
    W("---")
    W("")
    return L


if __name__ == "__main__":
    L = main()
    open("/tmp/claude-1000/-home-ubuntu-vcc/52816e9e-e42f-4af3-b6ec-69a5b07f6ce9/scratchpad/rep_f.txt",
         "w").write("\n".join(L))
    print("ok, lineas:", len(L))
