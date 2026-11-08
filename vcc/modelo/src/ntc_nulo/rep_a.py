#!/usr/bin/env python
"""Generador del reporte `ntc_nulo.md` -- PARTE 0 y encabezado."""
import json, os

OUTDIR = "/home/ubuntu/vcc/modelo/out/ntc_nulo"
J = lambda n: json.load(open(os.path.join(OUTDIR, n)))


def f(x, n=4):
    if x is None:
        return "n/d"
    if isinstance(x, bool):
        return "si" if x else "**no**"
    if isinstance(x, int):
        return f"{x:,}"
    return f"{x:,.{n}f}"


def e(x, n=3):
    return "n/d" if x is None else f"{x:.{n}e}"


def main():
    p0, p0b, p1 = J("p0.json"), J("p0b.json"), J("p1.json")
    L = []
    W = L.append

    W("# El nulo de `ntc_id` en contexto, y las tres verificaciones de validez")
    W("")
    W("## 0. Que es este reporte y que no")
    W("")
    W("Esto PRESENTA DATOS. No lleva veredicto, no concluye, y no va dirigido a")
    W("ningun candidato siguiente. Cada construccion se declara **antes** de la")
    W("medicion que la usa.")
    W("")
    W("**La cantidad que se mide.** Cada `context_{A..F}.h5ad` trae 18,400 celulas")
    W("NO DIRIGIDAS y la columna `ntc_id` nombra cual de 46 guias lleva cada una,")
    W("400 por guia. Esa es estructura de replica DENTRO del control, y existe en")
    W("A/B/C y en D/E/F. Partir las 46 guias en dos mitades y correr la puerta de")
    W("significancia del preset `vcc2026` mitad contra mitad da la distribucion")
    W("NULA de expresion diferencial aparente entre dos grupos de celulas no")
    W("dirigidas, en el contexto que se entrega, con su profundidad y su `n`. Es")
    W("independiente de la convencion del marcador **por construccion**: son")
    W("controles contra controles, no hay perturbacion que acertar.")
    W("")
    W("**Lo que se leyo**: el documento de formatos del reto")
    W("(`modelo/out/paso1_contrato.md`), `modelo/out/control_permutado.md` §6-§7,")
    W("y el codigo de la puerta (`cell_eval2/configs/vcc2026.yaml`,")
    W("`cell_eval2/de_compute.py`, `cell_eval2/de.py`).")
    W("")
    W("**Lo que NO se toco**: el lado perturbado -- no existe en 2026, y el de")
    W("2025/Jiang no entra aqui; de Jiang solo se leyeron las celulas `NT`. No se")
    W("leyo nada de `~/vcc/benchmark/` ni de `/datos/benchmark_VEDADO_*`. Sin")
    W("envios, sin `vcc sample`, sin commits, y `vcc prep` no se corrio en absoluto.")
    W("")
    W("**La puerta, leida del propio preset** (ningun numero se escribio a mano:")
    W("todos salen de `EvalConfig.from_preset(\"vcc2026\")`):")
    W("")
    W("| parametro | valor |")
    W("|---|---|")
    for k, v in p1["puerta"].items():
        W(f"| `{k}` | `{v}` |")
    W("")
    W("El camino, el mismo que recorre el marcador en su lado de CPU:")
    W("`compute_de` (normalize_total a `target_sum` -> medias aritmeticas -> log2FC")
    W("con `epsilon` -> `log1p` -> MWU del backend -> `_apply_cpm_filter`, que")
    W("conserva los genes cuyo **CPM medio en el grupo de REFERENCIA** supera")
    W("`filter_gene_min_cpm_cell` y **recalcula BH por blanco** sobre las filas que")
    W("quedan) -> `prep_de_side` (`nan_lfc_policy='mask'`, `min_abs_log2fc=0.0`,")
    W("que es no-op) -> conjunto DE = filas con `p_adj < 0.05`.")
    W("")
    W("**Verificacion de que el arnes CORRE la puerta y no la reimplementa.** Se")
    W("comprobo, sobre una corrida (contexto D, semilla `20274300`, 400 v 400), que")
    W("las cuatro piezas que definen el conjunto DE salen exactas")
    W("(`modelo/out/ntc_nulo/verif_puerta.log`):")
    W("")
    W("| # | lo que define la puerta | esperado | obtenido | igual |")
    W("|---|---|---|---|---|")
    W("| 1 | conjunto probado = genes con CPM medio de la REFERENCIA > 5.0 | 9,619 | "
      "9,619 | si |")
    W("| 2 | `p_adj` = Benjamini-Hochberg de `p_value` DENTRO de ese conjunto | — | "
      "`max|dif| = 0.000e+00` | si |")
    W("| 3 | conjunto DE = filas con `p_adj < 0.05` | 0 | 0 | si |")
    W("| 4 | `log2FC = log2((media_CPM_blanco + eps)/(media_CPM_ref + eps))` | — | "
      "`max|dif| = 1.289e-15` | si |")
    W("")
    W("🔴 **El backend resuelto aqui es `pdex`, no `gpudge`**: `de.backend='auto'`")
    W("y la maquina no tiene CUDA. El propio `cell_eval2` lo avisa (\"DE numbers")
    W("differ between engines\"). Todas las cifras de este reporte son de `pdex`.")
    W("")
    W("---")
    W("")

    # ---------------------------------------------------------------- 0a
    W("## PARTE 0 — Las tres verificaciones de validez")
    W("")
    W("Las tres se hicieron **primero** y son bloqueantes. Ninguna se disparo.")
    W("")
    W("### 0.a — `controls-final` esta bajado y completo")
    W("")
    a = p0["a_validez_datos"]
    W("**TABLA 0.1 — los 14 archivos de `/datos/reto`, con `sha256` de cada uno.**")
    W("")
    W("| # | particion | archivo | bytes | sha256 |")
    W("|---|---|---|---|---|")
    for i, r in enumerate(a["archivos"], 1):
        W(f"| {i} | `{r['particion']}` | `{r['archivo']}` | {r['bytes']:,} | "
          f"`{r['sha256']}` |")
    W("")
    W("Los tres `final/context_{D,E,F}.h5ad`, `final/gene_names.csv` y")
    W("`final/pert_counts.csv` estan los cinco, mas el `vcc_2026_controls_final.zip`")
    W("de origen. `val/` esta completo igual.")
    W("")
    W("**Y el contenido de los dos `.zip` de origen, listado**: cada uno trae")
    W("**exactamente 6 archivos** y los tamanos coinciden byte a byte con los ya")
    W("extraidos (p. ej. `context_D.h5ad` = 217,126,334 en el zip y en disco):")
    W("")
    W("| zip | fecha | archivos | contenido |")
    W("|---|---|---|---|")
    W("| `val/vcc_2026_controls.zip` | 2026-08-19 | 6 | `context_{A,B,C}.h5ad`, "
      "`gene_names.csv`, `pert_counts.csv`, `manifest.json` |")
    W("| `final/vcc_2026_controls_final.zip` | 2026-10-09 | 6 | `context_{D,E,F}.h5ad`, "
      "`gene_names.csv`, `pert_counts.csv`, `manifest.json` |")
    W("")
    W("🔴 **Ni el zip ni el directorio traen una sola celula perturbada.** Los 6")
    W("archivos de cada lado son los tres contextos de CONTROL mas los tres de panel.")
    W("Esa es la verificacion literal de la afirmacion que la seccion 2.4 necesita.")
    W("")
    W("**TABLA 0.2 — forma y conteos de los seis contextos, leidos del h5ad.**")
    W("")
    W("| # | ctx | particion | celulas | genes | nnz almacenado | `nnz==indptr[-1]` | "
      "`target_gene` | `context` | `ntc_id` cats | celulas/guia | suma |")
    W("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for i, r in enumerate(a["contextos"], 1):
        pg = ("400 en las 46" if r["ntc_id_todas_400"]
              else f"**min {r['ntc_id_celulas_min']} max {r['ntc_id_celulas_max']}**")
        W(f"| {i} | **{r['contexto']}** | `{r['particion']}` | {r['n_celulas']:,} | "
          f"{r['n_genes']:,} | {r['nnz_total']:,} | "
          f"{f(r['nnz_consistente_con_indptr'])} | "
          f"`{r['target_gene_categorias'][0]}` | `{r['context_categorias'][0]}` | "
          f"{r['ntc_id_n_categorias']} | {pg} | {r['ntc_id_suma']:,} |")
    W("")
    W(f"Los 46 `ntc_id` son **el mismo conjunto en los seis** archivos: "
      f"`{f(a['ntc_id_identicos_los_6'])}`. 46 x 400 = 18,400 en cada contexto, sin")
    W("excepcion, y `18400 = ntc_id_suma` cierra el conteo.")
    W("")
    W("**`manifest.json` esta JUNTO a `pert_counts.csv`** en las dos particiones --")
    W("que es de donde `vcc prep` lee los contextos:")
    W("")
    W("| particion | `manifest.json` existe | `pert_counts.csv` existe | mismo directorio |")
    W("|---|---|---|---|")
    for part, v in sorted(a["manifest"].items()):
        W(f"| `{part}` | {f(v['manifest_existe'])} | {f(v['pert_counts_existe'])} | "
          f"{f(v['mismo_directorio'])} |")
    W("")
    W("Y lo que dice el de `final/`, verbatim en sus campos:")
    W("")
    mf = a["manifest"]["final"]["contenido"]
    W("| campo | valor |")
    W("|---|---|")
    for k in sorted(mf):
        if k == "per_context":
            for c in sorted(mf[k]):
                W(f"| `per_context.{c}` | `{json.dumps(mf[k][c], sort_keys=True)}` |")
        else:
            W(f"| `{k}` | `{json.dumps(mf[k])}` |")
    W("")
    W("🔴 `ground_truth_cells = 138,400` por contexto y `control_cells = 18,400`:")
    W("**138,400 = 300 x 400 + 18,400**. El lado real contra el que se mide una")
    W("perturbacion es el pool de control de 18,400 celulas, no otras 400. Esa")
    W("aritmetica es la razon del tercer brazo de `n` de la Parte 1.")
    W("")

    # ---------------------------------------------------------------- 0b
    W("### 0.b — Cobertura real POR SIMBOLO contra gwps")
    W("")
    W("🔴 **Correccion de lectura, declarada.** `gene_names.csv` trae **cabecera**")
    W("(`gene_name`), asi que contarlo por lineas da 18,534. Los simbolos son")
    W("**18,533**, y la cuenta de abajo los lee con `csv.reader` saltando la")
    W("cabecera. Se cruzo ademas contra `var/_index` del h5ad, elemento por")
    W("elemento.")
    W("")
    W("gwps = `/datos/ref/replogle/k562_gwps_raw_bulk.h5ad`. Dos conjuntos")
    W("distintos que no hay que confundir: los genes **MEDIDOS** (las 8,248")
    W("columnas de `var`) y los blancos con **PERFIL** (las filas de `obs`,")
    W("nombradas `*_GEN_CONSTRUCCION`).")
    W("")
    g = p0b["gwps"]
    W("| cantidad de gwps | valor |")
    W("|---|---|")
    W(f"| filas totales | {g['filas']:,} |")
    W(f"| filas `core_control` | {g['filas_core_control']:,} |")
    W(f"| columnas de `var` | {g['var_columnas']:,} |")
    W(f"| nombres de gen unicos | {g['var_nombres_unicos']:,} |")
    W(f"| nombres duplicados (se excluyen enteros) | "
      f"{', '.join('`%s`' % x for x in g['var_nombres_duplicados'])} |")
    W(f"| blancos con perfil, unicos | {g['blancos_con_perfil_unicos']:,} |")
    W("")
    W("**TABLA 0.3 — cobertura por simbolo, las dos particiones.**")
    W("")
    W("| # | cantidad | `val` | `final` |")
    W("|---|---|---|---|")
    filas = [
        ("simbolos en `gene_names.csv`", "gene_names_n_simbolos"),
        ("simbolos unicos", "gene_names_unicos"),
        ("`n_vars` del h5ad", "n_genes_h5ad"),
        ("orden identico a `var/_index`, elemento por elemento", "orden_identico_a_var_index"),
        ("**de los 18,533: cuantos estan MEDIDOS en gwps**", "en_gwps_medidos"),
        ("idem, excluyendo los 2 nombres duplicados", "en_gwps_medidos_sin_duplicados"),
        ("de los 18,533: cuantos NO estan en gwps", "fuera_de_gwps"),
        ("cobertura %", "cobertura_pct"),
        ("perturbaciones en `pert_counts.csv`", "pert_n"),
        ("**de las 300: cuantas tienen PERFIL en gwps**", "pert_con_perfil_en_gwps"),
        ("de las 300: cuantas NO tienen perfil", "pert_sin_perfil_en_gwps"),
        ("las 300 estan todas en `gene_names.csv`", "pert_todas_en_gene_names"),
    ]
    for i, (lab, k) in enumerate(filas, 1):
        W(f"| {i} | {lab} | {f(p0b['val'][k])} | {f(p0b['final'][k])} |")
    W("")
    for part in ("val", "final"):
        lst = p0b[part]["pert_lista_sin_perfil"]
        W(f"**Las {len(lst)} perturbaciones de `{part}/pert_counts.csv` SIN perfil en gwps**")
        W(f"(lista completa): {', '.join('`%s`' % x for x in lst)}.")
        W("")
    W("---")
    W("")
    return L


if __name__ == "__main__":
    open("/tmp/claude-1000/-home-ubuntu-vcc/52816e9e-e42f-4af3-b6ec-69a5b07f6ce9/scratchpad/rep_a.txt",
         "w").write("\n".join(main()))
    print("ok")
