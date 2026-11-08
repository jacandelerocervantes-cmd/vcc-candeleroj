#!/usr/bin/env python
"""Reporte `ntc_nulo.md` -- PARTE 2: lo que el nulo contesta."""
import json, os
import numpy as np
OUTDIR = "/home/ubuntu/vcc/modelo/out/ntc_nulo"
ETQ = {"completo": "`completo` 9,200 v 9,200",
       "n400_400": "`n400_400` 400 v 400",
       "n400_ref_completa": "`n400_ref_completa` 400 v 9,200"}
BRAZOS = ["completo", "n400_400", "n400_ref_completa"]


def main():
    p2 = json.load(open(os.path.join(OUTDIR, "p2ab.json")))
    p2c = json.load(open(os.path.join(OUTDIR, "p2c.json")))
    p1 = json.load(open(os.path.join(OUTDIR, "p1.json")))
    L = []; W = L.append
    A = p2["alpha"]

    W("## PARTE 2 — Lo que el nulo contesta")
    W("")
    W("### 2.1 🔴 La tasa de falsos positivos del conjunto ED en contexto, con `n = 400`")
    W("")
    W("**Dos cantidades distintas que no hay que confundir.**")
    W("")
    W("* **`n_sig` por corrida**: cuantos genes pasa la puerta en UNA comparacion")
    W("  control-contra-control. Es la cantidad comparable con el `n_DE` de una")
    W("  perturbacion, y es la que contesta \"cuantos genes pasa la puerta sobre puro")
    W("  ruido\".")
    W("* **FPR POR GEN**: `suma(n_sig) / suma(m)`. Es la tasa por prueba. Por")
    W(f"  construccion de BH se espera MUY por debajo de {A}: {A} es el nivel de FDR")
    W("  del conjunto llamado, no la tasa por prueba.")
    W("")
    W("**TABLA 2.1 — la tasa de falsos positivos, por brazo.**")
    W("")
    W("| # | brazo | `n` ref | `n` blanco | corridas | **`n_sig` media** | **`n_sig` sd** | "
      "mediana | **[min, max]** | corridas con `n_sig>0` | suma `n_sig` | `m` medio | "
      "**FPR por gen** | `p_crit` [min, max] | `p` min sobre corridas |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, b in enumerate(BRAZOS, 1):
        d = p2["brazos"].get(b)
        if not d:
            continue
        W(f"| {i} | {ETQ[b]} | {d['n_ref']:,} | {d['n_tgt']:,} | {d['n_corridas']} | "
          f"**{d['n_sig_media']:.4f}** | **{d['n_sig_sd']:.4f}** | "
          f"{d['n_sig_mediana']:.1f} | **[{d['n_sig_min']}, {d['n_sig_max']}]** | "
          f"{d['n_corridas_con_alguno']}/{d['n_corridas']} | {d['n_sig_total']} | "
          f"{d['m_media']:,.1f} | **{d['fpr_por_gen']:.3e}** | "
          f"[{d['p_crit_rango'][0]:.3e}, {d['p_crit_rango'][1]:.3e}] | "
          f"{d['p_min_sobre_corridas_min']:.3e} |")
    W("")
    W("**TABLA 2.2 — `n_sig` por contexto y brazo, con sus 8 valores crudos.**")
    W("")
    W("| # | brazo | ctx | los 8 valores de `n_sig` (j=0..7) | media | sd | **[min, max]** |")
    W("|---|---|---|---|---|---|---|")
    i = 0
    for b in BRAZOS:
        d = p2["brazos"].get(b)
        if not d:
            continue
        for c, v in sorted(d["n_sig_por_contexto"].items()):
            i += 1
            W(f"| {i} | {ETQ[b]} | **{c}** | `{v['valores']}` | {v['media']:.3f} | "
              f"{v['sd']:.3f} | **[{v['min']}, {v['max']}]** |")
    W("")
    n4 = p2["brazos"]["n400_400"]
    nr = p2["brazos"]["n400_ref_completa"]
    W("**Lo que dicen esas filas, como aritmetica:**")
    W("")
    W(f"* Con **`n = 400` a los dos lados**: `n_sig` vale **{n4['n_sig_min']} en las")
    W(f"  {n4['n_corridas']} corridas** ({n4['n_corridas_con_alguno']} de")
    W(f"  {n4['n_corridas']} con alguno), media **{n4['n_sig_media']:.4f}**, sd")
    W(f"  **{n4['n_sig_sd']:.4f}**, FPR por gen **{n4['fpr_por_gen']:.3e}** sobre")
    W(f"  {n4['n_sig_total'] if n4['n_sig_total'] else 0} genes llamados de")
    W(f"  {int(n4['m_media'] * n4['n_corridas']):,} pruebas.")
    W(f"* Con **400 contra el pool de 9,200** (el brazo que se parece a la entrega):")
    W(f"  media **{nr['n_sig_media']:.4f}**, sd **{nr['n_sig_sd']:.4f}**, rango")
    W(f"  **[{nr['n_sig_min']}, {nr['n_sig_max']}]**, "
      f"{nr['n_corridas_con_alguno']}/{nr['n_corridas']} corridas con alguno.")
    co = p2["brazos"]["completo"]
    W(f"* Con **9,200 contra 9,200**: media **{co['n_sig_media']:.4f}**, sd")
    W(f"  **{co['n_sig_sd']:.4f}**, rango **[{co['n_sig_min']}, {co['n_sig_max']}]**,")
    W(f"  {co['n_corridas_con_alguno']}/{co['n_corridas']} corridas con alguno.")
    W("")
    W("**El diagnostico de uniformidad de la seccion 1.4** cruza con esto: la fraccion")
    W(f"de `p` crudas bajo 0.05 promedia 0.0501-0.0539 segun el brazo (esperado 0.05")
    W("bajo un nulo uniforme) y el KS contra la uniforme promedia 0.013-0.014.")
    W("")

    # ------------------------------------------------ 2.2 piso
    W("### 2.2 🔴 El piso de deteccion, calculado del propio nulo")
    W("")
    W("**La construccion, declarada.** La puerta llama un gen cuando `p_adj < 0.05`.")
    W("Para una senal **AISLADA** (el resto del panel nulo), el escalon de BH de rango")
    W("1 exige `p <= alpha/m`. Ese `p_crit` se calcula de cada corrida (`m` = genes")
    W("probados tras el filtro CPM), no se importa de ninguna formula.")
    W("")
    W("El piso depende FUERTE del nivel de expresion del gen, asi que se mide por")
    W("banda de CPM medio de la referencia. Para cada banda se dan:")
    W("")
    W("* el `|log2FC|` MINIMO que el nulo alcanzo en cada umbral de `p` -- lectura")
    W("  directa, sin modelo;")
    W("* el **techo de ruido**: el `|log2FC|` maximo que el nulo PRODUCE en esa banda;")
    W("* el **alcance**: `n_filas x p_crit` dice cuantas filas del nulo se esperan por")
    W("  debajo de `p_crit`. Si es << 1, el piso en `p_crit` **no se puede leer")
    W("  directo**.")
    W("")
    for b in ("n400_400", "n400_ref_completa"):
        d = p2["brazos"].get(b, {}).get("piso_deteccion")
        if not d:
            continue
        W(f"#### Brazo {ETQ[b]}")
        W("")
        W(f"Pool de **{d['n_filas_pool']:,} filas** "
          f"({p2['brazos'][b]['n_archivos_pool']} corridas x ~{d['m_media']:,.0f} genes).")
        W(f"`p_crit = 0.05/{d['m_media']:,.1f} = {d['p_crit']:.4e}`. "
          f"`p` minimo del pool = **{d['p_min_pool']:.4e}**.")
        W("")
        W(f"🔴 **El alcance.** Filas esperadas bajo `p_crit`: "
          f"**{d['filas_esperadas_bajo_p_crit']:.2f}**. Observadas: "
          f"**{d['filas_observadas_bajo_p_crit']}**. `p_crit` "
          f"{'**SI** cae' if d['p_crit_dentro_del_alcance'] else '**NO** cae'} dentro "
          f"del alcance directo de este nulo.")
        if d.get("lfc_de_las_filas_bajo_p_crit"):
            W("")
            W(f"Las filas del nulo que SI bajaron de `p_crit` llevan "
              f"`|log2FC|` = {d['lfc_de_las_filas_bajo_p_crit']} con CPM de referencia "
              f"{d['cpm_de_las_filas_bajo_p_crit']}.")
        W("")
        W("**TABLA 2.3%s — el piso de deteccion por banda de expresion.**" %
          ("a" if b == "n400_400" else "b"))
        W("")
        W("`F(t) = min{|log2FC| : p <= t}` sobre el pool, que es **no creciente en")
        W("`t`**. Como `p_crit` es menor que cada ancla alcanzada, **cada ancla da una")
        W("cota INFERIOR medida** del piso. Una fila del nulo que SI cruzo `p_crit` da")
        W("una cota SUPERIOR medida.")
        W("")
        W("| # | banda CPM ref | filas | `|lfc|` mediana | **techo de ruido (`|lfc|` max)** | "
          "`p` min | `F(1e-2)` (n) | `F(1e-3)` (n) | `F(1e-4)` (n) | `F(1e-5)` (n) | "
          "**cota INFERIOR medida** | **cota SUPERIOR medida** | filas esp. bajo `p_crit` | "
          "filas obs. | **extrapolacion a `p_crit`** | pend./decada | decadas extrap. |")
        W("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for i, bb in enumerate(d["bandas"], 1):
            lo, hi = bb["banda"]
            nom = f"({lo}, {hi}]" if hi else f"> {lo}"
            if bb.get("insuficiente"):
                W(f"| {i} | {nom} | {bb['n_filas']:,} | — | — | — | — | — | — | — | — "
                  f"| — | — | — | *menos de 50 filas: no se mide* | — | — |")
                continue
            anc = bb["anclas"]
            def c(k):
                a = anc.get(k, {})
                return (f"{a['lfc_min']:.4f} ({a['n']})" if a.get("lfc_min") is not None
                        else "— (0)")
            ci, cs = bb.get("cota_inferior_medida"), bb.get("cota_superior_medida")
            cis = (f"**>= {ci['lfc_min']:.4f}** (de `p<={ci['umbral_p_alcanzado']:.0e}`, "
                   f"n={ci['n_filas_en_el_umbral']})" if ci else "**n/d**")
            css = (f"**<= {cs['lfc_min_de_las_filas_que_cruzaron']:.4f}** "
                   f"(n={cs['n_filas']})" if cs else "**n/d** -- ninguna fila cruzo")
            ex = bb.get("extrapolacion_declarada")
            W(f"| {i} | {nom} | {bb['n_filas']:,} | {bb['lfc_mediana_del_nulo']:.4f} | "
              f"**{bb['lfc_max_del_nulo']:.4f}** | {bb['p_min_de_la_banda']:.2e} | "
              f"{c('1e-02')} | {c('1e-03')} | {c('1e-04')} | {c('1e-05')} | {cis} | "
              f"{css} | {bb['filas_esperadas_bajo_p_crit']:.3f} | "
              f"{bb['filas_observadas_bajo_p_crit']} | "
              + (f"`|lfc|` = **{ex['lfc_en_p_crit']:.4f}** (EXTRAPOLADO) | "
                 f"{ex['pendiente_por_decada']:.4f} | {ex['decadas_de_extrapolacion']:.2f} |"
                 if ex else "n/d (menos de 3 anclas) | — | — |"))
        W("")
    W("🔴 **Lo que se puede y lo que no se puede calcular de este nulo.**")
    W("")
    d4 = p2["brazos"]["n400_400"]["piso_deteccion"]
    dr = p2["brazos"]["n400_ref_completa"]["piso_deteccion"]
    bb4 = [x for x in d4["bandas"] if not x.get("insuficiente")]
    bbr = [x for x in dr["bandas"] if not x.get("insuficiente")]
    W("1. **SE PUEDE, y esta medido: una cota INFERIOR del piso en cada banda.**")
    W("   `F` es no creciente y `p_crit` esta por debajo de cada ancla alcanzada, asi")
    W("   que `F(ancla)` acota el piso por abajo sin ningun modelo. Con `n = 400` a")
    W("   los dos lados va de")
    W(f"   **{min(x['cota_inferior_medida']['lfc_min'] for x in bb4 if x.get('cota_inferior_medida')):.3f}**")
    W(f"   (banda de CPM mas alto) a")
    W(f"   **{max(x['cota_inferior_medida']['lfc_min'] for x in bb4 if x.get('cota_inferior_medida')):.3f}**")
    W("   (banda de CPM mas bajo) en `|log2FC|`.")
    cs_ok = [x for x in bbr if x.get("cota_superior_medida")]
    if cs_ok:
        x = cs_ok[0]
        W("2. **SE PUEDE, en UNA banda de UN brazo: un intervalo con las dos cotas**")
        W("   **medidas.** En el brazo 400 v 9,200, banda")
        W(f"   CPM ({x['banda'][0]}, {x['banda'][1]}], una fila del nulo cruzo `p_crit`")
        W(f"   (`p = {dr['p_min_pool']:.3e}`), asi que el piso queda **entre")
        W(f"   {x['cota_inferior_medida']['lfc_min']:.3f} y "
          f"{x['cota_superior_medida']['lfc_min_de_las_filas_que_cruzaron']:.3f}** en")
        W(f"   `|log2FC|`, las dos cotas medidas. La extrapolacion de esa misma banda da")
        W(f"   {x['extrapolacion_declarada']['lfc_en_p_crit']:.3f}, dentro del intervalo.")
    else:
        W("2. En ninguna banda de ningun brazo se midieron las dos cotas a la vez.")
    W("3. 🔴 **NO SE PUEDE leer directo el piso EXACTO en `p_crit`.** Las filas")
    W("   esperadas por debajo de `p_crit` son ~0.01-0.64 por banda (columna de la")
    W("   tabla): el nulo, como se construyo, no llega ahi. Donde se da una cifra, es")
    W("   una **extrapolacion DECLARADA** -- `F(t)` lineal en `log10(t)` ajustada a las")
    W("   anclas MEDIDAS de esa banda -- y va marcada `(EXTRAPOLADO)` con su pendiente")
    W("   por decada y cuantas decadas se extrapola (0.29 donde se alcanzo `p<=1e-5`,")
    W("   1.29 donde solo se alcanzo `p<=1e-4`). **Es un rango, no un numero medido.**")
    W("   La forma log-lineal **es** una suposicion importada; el nulo solo da los puntos.")
    W("4. **Medido directo: el techo de ruido por banda** -- el `|log2FC|` maximo que")
    W("   el nulo PRODUCE. 🔴 En la banda de CPM mas bajo llega a")
    W(f"   **{max(x['lfc_max_del_nulo'] for x in bb4):.2f}** con `n = 400`: un")
    W("   `|log2FC|` enorme es compatible con puro ruido cuando el gen tiene ceros en")
    W("   un lado (`epsilon = 1e-9` contra una media de CPM positiva). **Ningun umbral")
    W("   de `|log2FC|` solo separa senal de ruido a `n = 400`**; la potencia de la")
    W("   puerta viene del test de rangos, no del tamano de efecto. Esa cola gruesa es")
    W("   tambien la razon de que la extrapolacion de las dos bandas mas bajas tenga la")
    W("   pendiente por decada mas alta (~0.30 contra ~0.05 en la banda mas alta).")
    W("")
    # ------------------------------------------------ 2.3 guia a guia
    W("### 2.3 La variacion guia a guia: cuanto del nulo lo explica la identidad de la guia")
    W("")
    W("**La construccion, declarada.** Dos mediciones, ninguna parte una guia:")
    W("")
    W("1. **Particion que AGRUPA guias contra particiones ALEATORIAS de guias.** La")
    W("   particion aleatoria de la Parte 1 reparte las 46 guias al azar. La que")
    W("   AGRUPA las ordena por su propio perfil y corta por la mitad: las 23 de un")
    W("   extremo contra las 23 del otro. Es el peor caso que la estructura de guia")
    W("   permite **sin partir ninguna guia**.")
    W("   * perfil de guia = media del CPM por celula sobre sus 400 celulas,")
    W("     restringido al panel probable, `log1p`, centrado por gen;")
    W("   * se ordena por la **primera componente principal** de esa matriz 46 x G;")
    W("   * mitad REFERENCIA = las 23 de PC1 mas baja, BLANCO = las 23 mas alta;")
    W("   * **determinista: no lleva semilla**. Solo el submuestreo a 400 la lleva")
    W("     (`20274500 + idx_contexto`).")
    W("2. **Descomposicion de varianza por gen** (ANOVA de un factor, 46 guias x 400")
    W("   celulas), en `log1p(CPM)` -- el espacio que ve el Wilcoxon antes de pasar a")
    W("   rangos:")
    W("   `MSB = SSB/45`, `MSW = SSW/(18400-46)`, `F = MSB/MSW`,")
    W("   `ICC = (MSB - MSW)/(MSB + 399*MSW)`.")
    W("   `F` dice si la identidad de la guia mueve la media del gen mas de lo que la")
    W("   mueve el muestreo de 400 celulas; `ICC` dice QUE FRACCION de la varianza")
    W("   por celula se lleva la guia.")
    W("")
    W("**TABLA 2.4 — la particion que agrupa guias contra las 8 aleatorias, por contexto.**")
    W("")
    W("| # | ctx | brazo | `n_sig` AGRUPA guias | `n_sig` aleatorias: media (sd) | "
      "aleatorias [min, max] | `m` AGRUPA | `p` min AGRUPA | `p` min aleatorias (mediana) | "
      "`|lfc|` max AGRUPA | genes llamados por AGRUPA |")
    W("|---|---|---|---|---|---|---|---|---|---|---|")
    i = 0
    cors = p1["corridas"]
    for cc in p2c["contextos"]:
        c = cc["contexto"]
        for cr in cc["corridas"]:
            i += 1
            b = cr["brazo"]
            rs = [r for r in cors if r["contexto"] == c and r["brazo"] == b]
            ns = np.array([r["n_sig"] for r in rs]) if rs else np.array([np.nan])
            pm = np.array([r["p_min"] for r in rs], dtype=float) if rs else np.array([np.nan])
            if cr.get("genes_sig"):
                gl = ", ".join("`%s`" % g for g in cr["genes_sig"])
            elif cr["n_sig"] > 0:
                gl = f"*{cr['n_sig']:,} genes: lista omitida (>200)*"
            else:
                gl = "—"
            W(f"| {i} | **{c}** | {ETQ[b]} | **{cr['n_sig']}** | "
              f"{np.mean(ns):.3f} ({np.std(ns, ddof=1):.3f}) | "
              f"[{int(np.min(ns))}, {int(np.max(ns))}] | {cr['n_probados']:,} | "
              f"{cr['p_min']:.3e} | {np.median(pm):.3e} | {cr['abs_lfc_max']:.4f} | {gl} |")
    W("")
    W("**TABLA 2.5 — la descomposicion de varianza por gen, por contexto.**")
    W("`n` genes evaluados = los del panel probable con `MSW > 0`.")
    W("")
    W("| # | ctx | panel probable | genes evaluados | **genes con guia significativa (BH 0.05)** | "
      "**frac** | `F` mediana | `F` p99 | `F` max | **`ICC` mediana** | `ICC` p75 | "
      "`ICC` p99 | `ICC` max | frac `ICC<0` | PC1 % var | gl entre | gl dentro |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, cc in enumerate(p2c["contextos"], 1):
        a = cc["anova"]
        W(f"| {i} | **{cc['contexto']}** | {cc['n_panel_probable']:,} | "
          f"{a['n_genes_evaluados']:,} | **{a['n_genes_guia_significativa_BH']:,}** | "
          f"**{a['frac_genes_guia_significativa']:.5f}** | {a['F_mediana']:.4f} | "
          f"{a['F_p99']:.4f} | {a['F_max']:.4f} | **{a['ICC_mediana']:.6f}** | "
          f"{a['ICC_p75']:.6f} | {a['ICC_p99']:.6f} | {a['ICC_max']:.6f} | "
          f"{a['ICC_frac_negativo']:.4f} | "
          f"{cc['pc1_varianza_explicada'] * 100:.2f} | {a['df_entre']} | "
          f"{a['df_dentro']:,} |")
    W("")
    def agr(brazo):
        return {cc["contexto"]: [cr["n_sig"] for cr in cc["corridas"]
                                 if cr["brazo"] == brazo][0]
                for cc in p2c["contextos"]}
    def ale(brazo):
        o = {}
        for cc in p2c["contextos"]:
            v = [r["n_sig"] for r in cors
                 if r["contexto"] == cc["contexto"] and r["brazo"] == brazo]
            o[cc["contexto"]] = v
        return o
    W("**TABLA 2.6 — el contraste directo: AGRUPA contra aleatorias, por brazo y contexto.**")
    W("")
    W("| # | brazo | ctx | **`n_sig` AGRUPA** | `n_sig` aleatorias (los 8) | "
      "aleatorias media | **AGRUPA - media aleatoria** |")
    W("|---|---|---|---|---|---|---|")
    i = 0
    for b in BRAZOS:
        A_, R_ = agr(b), ale(b)
        for c in sorted(A_):
            i += 1
            mu = float(np.mean(R_[c]))
            W(f"| {i} | {ETQ[b]} | **{c}** | **{A_[c]:,}** | `{R_[c]}` | {mu:.3f} | "
              f"**{A_[c] - mu:+,.3f}** |")
    W("")
    W("**Lo que esas tres tablas dicen, como aritmetica:**")
    W("")
    Ac, Rc = agr("completo"), ale("completo")
    A4, R4 = agr("n400_400"), ale("n400_400")
    W(f"* En el brazo `completo`, la particion que AGRUPA guias llama")
    W(f"  **{', '.join('%s=%s' % (c, f'{Ac[c]:,}') for c in sorted(Ac))}** genes,")
    W(f"  contra una media de **{np.mean([v for c in Rc for v in Rc[c]]):.4f}** (rango")
    W(f"  [{min(v for c in Rc for v in Rc[c])}, {max(v for c in Rc for v in Rc[c])}]) de")
    W("  las 48 particiones aleatorias. 🔴 **El contexto D llega a")
    W(f"  {Ac['D']:,} genes de {[cr['n_probados'] for cc in p2c['contextos'] if cc['contexto'] == 'D' for cr in cc['corridas'] if cr['brazo'] == 'completo'][0]:,}**")
    W("  probados, y es tambien el contexto con mas varianza entre guias capturada por")
    W(f"  PC1 ({[c['pc1_varianza_explicada'] for c in p2c['contextos'] if c['contexto'] == 'D'][0] * 100:.2f} %).")
    W(f"* En el brazo `n400_400`, donde las 48 aleatorias dan **0 en 48**, la particion")
    W(f"  que AGRUPA da **{', '.join('%s=%d' % (c, A4[c]) for c in sorted(A4))}**.")
    W(f"* El `F` mediano por gen esta en")
    W(f"  [{min(c['anova']['F_mediana'] for c in p2c['contextos']):.4f}, "
      f"{max(c['anova']['F_mediana'] for c in p2c['contextos']):.4f}] sobre los seis")
    W("  contextos. Un `F` mediano de 1 es exactamente \"la guia no explica nada mas")
    W("  que el muestreo\".")
    W(f"* Los genes con efecto de guia significativo a BH 0.05 van de")
    W(f"  **{min(c['anova']['n_genes_guia_significativa_BH'] for c in p2c['contextos']):,}**")
    W(f"  a **{max(c['anova']['n_genes_guia_significativa_BH'] for c in p2c['contextos']):,}**")
    W(f"  por contexto, es decir")
    W(f"  **{min(c['anova']['frac_genes_guia_significativa'] for c in p2c['contextos']) * 100:.3f}-"
      f"{max(c['anova']['frac_genes_guia_significativa'] for c in p2c['contextos']) * 100:.3f} %**")
    W("  del panel probable.")
    W(f"* El `ICC` mediano va de")
    W(f"  **{min(c['anova']['ICC_mediana'] for c in p2c['contextos']):.6f}** a")
    W(f"  **{max(c['anova']['ICC_mediana'] for c in p2c['contextos']):.6f}** -- la guia se")
    W("  lleva una fraccion diminuta de la varianza POR CELULA de un gen tipico -- y su")
    W(f"  maximo por contexto llega a")
    W(f"  **{max(c['anova']['ICC_max'] for c in p2c['contextos']):.4f}**.")
    dD = [c["anova"] for c in p2c["contextos"] if c["contexto"] == "D"][0]
    W("* 🔴 **El contexto D tiene el `F` MAXIMO mas bajo de los seis y a la vez las**")
    W("  **1,937 llamadas.** `F_max` en D vale")
    W(f"  **{dD['F_max']:.4f}** y su `ICC_max` **{dD['ICC_max']:.6f}**, los mas bajos de")
    W(f"  la TABLA 2.5 (los otros contextos llegan a `F_max` de "
      f"{min(c['anova']['F_max'] for c in p2c['contextos'] if c['contexto'] != 'D'):.1f}-"
      f"{max(c['anova']['F_max'] for c in p2c['contextos'] if c['contexto'] != 'D'):.1f}).")
    W("  Es decir: en D **ningun gen suelto** tiene un efecto de guia grande, pero PC1")
    W(f"  captura **{[c['pc1_varianza_explicada'] for c in p2c['contextos'] if c['contexto'] == 'D'][0] * 100:.2f} %**")
    W("  de la varianza entre guias, el doble que en cualquier otro contexto. La")
    W("  estructura de guia de D es **ancha y de poca amplitud** -- un desplazamiento")
    W("  pequeno repartido sobre muchos genes -- y eso es exactamente lo que un corte a")
    W("  lo largo de PC1 convierte en llamadas masivas, mientras que un ANOVA gen a gen")
    W("  no lo ve en ningun gen en particular. Las dos tablas miran la misma estructura")
    W("  por ejes distintos.")
    W("* 🔴 **Las dos mediciones no se contradicen, y la aritmetica de por que es**")
    W("  **esta:** el `ICC` por gen es diminuto, pero la puerta no compara celulas")
    W("  sueltas: compara las MEDIAS de dos grupos de 9,200. El error de una media de")
    W("  9,200 celulas es ~96 veces menor que la `sd` por celula, asi que un desplazami"
      "ento")
    W("  de medias que vale una fraccion minuscula de la varianza por celula sigue")
    W("  siendo muchas veces su propio error estandar. Esa es la razon de que `ICC`")
    W("  mediano ~1e-4 convivan con 1,937 genes llamados cuando las guias se agrupan.")
    W("")
    # ------------------------------------------------ 2.4
    W("### 2.4 🔴 Cuantas de las 300 perturbaciones esperarian conjunto ED vacio por potencia")
    W("")
    W("**El nulo NO permite contestarlo, y no se estima.**")
    W("")
    W("La razon es estructural, no de precision. El nulo mide una sola de las dos")
    W("mitades de la pregunta:")
    W("")
    W("| | lo que hace falta | lo que el nulo da |")
    W("|---|---|---|")
    W("| **el umbral** | el tamano de efecto que la puerta exige a `n = 400` con esta "
      "profundidad | **SI** -- las secciones 2.1 y 2.2, con su rango y su "
      "extrapolacion declarada |")
    W("| **la distribucion** | el tamano de efecto REAL de cada una de las 300 "
      "perturbaciones en cada contexto | **NO** -- esta en el lado perturbado |")
    W("")
    W("Contar cuantas perturbaciones caen por debajo del umbral exige las dos. El")
    W("lado perturbado de 2026 **no existe**: `manifest.json` declara")
    W("`ground_truth_cells: 138,400` por contexto, pero lo que se entrega son las")
    W("18,400 celulas de control (TABLA 0.2); las 120,000 celulas perturbadas no")
    W("estan en `/datos/reto`. Y el de 2025/Jiang esta fuera de este encargo.")
    W("")
    W("Sustituir la distribucion que falta por la de Jiang seria **suponer que el")
    W("regimen transfiere**, que es exactamente lo que la Parte 3 se niega a declarar")
    W("en cualquiera de los dos sentidos. Asi que aqui no hay cifra, y no se da una.")
    W("")
    W("Lo unico que el nulo acota sin verdad es la direccion contraria: **cuantas")
    W("perturbaciones esperarian ED vacio POR RUIDO DEL CONTROL**, que es")
    W(f"`n_sig` del brazo correspondiente -- **{n4['n_sig_media']:.4f}** genes de media")
    W("a `n = 400` contra `n = 400`. Eso es el piso de falsos positivos, no el conteo")
    W("de verdad vacia.")
    W("")
    W("---")
    W("")
    return L


if __name__ == "__main__":
    L = main()
    open("/tmp/claude-1000/-home-ubuntu-vcc/52816e9e-e42f-4af3-b6ec-69a5b07f6ce9/scratchpad/rep_e.txt",
         "w").write("\n".join(L))
    print("ok, lineas:", len(L))
