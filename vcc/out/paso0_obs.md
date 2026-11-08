# Paso 0 de T13 — `.obs` fuente por fuente (modo backed, sin tocar `X`)

**Método.** `h5py.File(path, "r")`, leyendo solo el grupo `/obs` (y el `attrs`/`shape` de
`/X` para contexto, nunca sus datos). Ningún archivo se cargó completo. Entorno aislado
(`uv venv` en `~/obs_venv`, interpretado con `python -I`), fuera de `/datos`, script en
`~/obs_scripts/inspect_h5ad_obs.py`, paths pasados como argumento — nunca se ejecutó nada
desde dentro de `/datos`. No se asumió ningún nombre de columna: se enumeraron los
presentes y luego se buscaron candidatas por substring (`well, lane, batch, lib, pool,
replicate, guide, sgrna, grna, construct, plate, run, flowcell, hash, sample, gem,
transduc, moi, …`).

---

## 0 · Hallazgo que cambia el plan de T13

Los seis `.h5ad` de `/datos/reto/{val,final}` — los que el inventario de archivos del
paso 0 anterior listó sin abrir — **no contienen perturbaciones**. Cada uno tiene
`n_obs = 18400`, y las tres columnas de `obs` (`context`, `ntc_id`, `target_gene`) son
categóricas con `target_gene` colapsado a una sola categoría: `"non-targeting"`. Los
`cell_barcode` (vía `_index`, codificado `nullable-string-array`) siguen el patrón
`"{contexto}_ctrl_000000"` — son los **18,400 controles** que el `manifest.json` ya
declaraba (`control_cells: 18400`), no los `ground_truth_cells: 138400` por contexto.

**Consecuencia:** T13 (réplicas por línea×perturbación) no se puede medir sobre los
archivos del reto — ahí solo hay controles. Tiene que correr sobre las fuentes de
referencia (`nadig`, `replogle`), que sí traen perturbaciones a nivel de célula.

---

## 1 · `/datos/reto/{val,final}/context_{A..F}.h5ad` — solo controles

Idéntico esquema en los seis (difieren solo en qué letra de contexto llevan). `shape`
reportado es `[n_obs, n_vars]` leído del atributo de `X`, no de sus datos.

| archivo | n_obs | shape |
|---|---|---|
| `val/context_A.h5ad` | 18400 | [18400, 18533] |
| `val/context_B.h5ad` | 18400 | [18400, 18533] |
| `val/context_C.h5ad` | 18400 | [18400, 18533] |
| `final/context_D.h5ad` | 18400 | [18400, 18533] |
| `final/context_E.h5ad` | 18400 | [18400, 18533] |
| `final/context_F.h5ad` | 18400 | [18400, 18533] |

Columnas de `obs` (iguales en los seis, solo cambia el valor de `context`):

| columna | dtype | nº únicos presentes | ejemplos |
|---|---|---|---|
| `_index` | `nullable-string-array<object>` | 18400 (no recalculado, es el barcode) | `A_ctrl_000000`, `A_ctrl_000001`, `A_ctrl_000002`, `A_ctrl_000003`, `A_ctrl_000004` |
| `context` | `categorical<object>` | 1 | `"A"` (o la letra correspondiente) |
| `ntc_id` | `categorical<object>` | 46 | `non-targeting-11`, `-12`, `-13`, `-15`, `-17` |
| `target_gene` | `categorical<object>` | 1 | `"non-targeting"` |

Ninguna candidata a (a)/(b)/(c): no hay perturbaciones que comparar. `ntc_id` con 46
únicos **confirma el lado "documento" (46) del cotejo pendiente de `PLAN.md` §8.1** —
sigue sin resolver el lado "100 (el reto)".

---

## 2 · `/datos/ref/nadig/{jurkat,hepg2}_raw_sc.h5ad`

| archivo | n_obs | shape |
|---|---|---|
| `jurkat_raw_sc.h5ad` | 262,956 | [262956, 8882] |
| `hepg2_raw_sc.h5ad` | 145,473 | [145473, 9624] |

Columnas completas de `obs` (idénticas en ambos archivos; cifras de Jurkat / HepG2):

| columna | dtype | únicos (Jurkat / HepG2) | candidata | ejemplos |
|---|---|---|---|---|
| `cell_barcode` | `object` (string) | 262956 / 145473 | no | `AAACCCAAGAAACTGT-27`, `AAACCCAAGAAATCCA-12`, … |
| `UMI_count` | `float32` | 32846 / 38700 | no | 1751.0, 1752.0, 1753.0, 1754.0, 1755.0 |
| `mitopercent` | `float32` | 245638 / 142677 | no | 0.0, 7.50e-05, 8.93e-05, … |
| **`gem_group`** | `int64` | **55 / 56** | **sí** | 1, 2, 3, 4, 5 |
| **`z_gemgroup_UMI`** | `float32` | 222301 / 136061 | sí (score derivado de `gem_group`, no un eje nuevo) | -2.591, -2.586, -2.585, … |
| `gene` | legacy-categorical (`__categories/gene`), `int16` | 2394 / 2394 | no (es el gen blanco, no la guía) | `AAAS`, `AAMP`, `AAR2`, `AARS`, `AARS2` |
| `gene_id` | legacy-categorical, `int16` | 2392 / 2392 | no | `ENSG00000001084`, `ENSG00000001497`, … |
| `gene_transcript` | legacy-categorical, `int16` | 2679 / 2662 | no | `2_AAAS_P1P2_ENSG00000094914`, … |
| **`sgID_AB`** | legacy-categorical, `int16` | **2679 / 2662** | **sí** | `AAAS_-_53715438.23-P1P2\|AAAS_+_53715355.23-P1P2`, … |
| `transcript` | legacy-categorical, `int8` | 62 / 62 | no | `ENST00000248058.1`, … |

(`legacy-categorical` = codificación anndata <0.7: enteros en `/obs/<col>` + etiquetas en
`/obs/__categories/<col>`, no el formato `categories`/`codes` moderno — se decodificó
igual, pero vale la nota porque un lector que asuma el formato nuevo lo leería mal.)

**Respuesta a lo que T13 pide, para esta fuente:**

- **(a) pozo/transducción aparte:** candidata = `gem_group` (55–56 valores, enteros
  secuenciales 1..N). Es el identificador de "GEM well" de 10x — cada valor es una
  carga/cápsula de Chromium independiente. **No puedo confirmar desde la columna sola
  si cada `gem_group` es también una transducción independiente**, o si algunas
  comparten pool de transducción y solo difieren en la carga del chip — esa semántica
  vive en el método del paper, no en el dato. El sufijo `-N` del `cell_barcode` coincide
  en rango con `gem_group` (ambos van hasta 55/56), así que **probablemente** sea el
  mismo eje expuesto dos veces (convención de `cellranger aggr`), no dos ejes
  independientes — queda sin confirmar contra la fuente primaria del paper.
- **(b) mismo pozo resecuenciado:** **no encontrada.** No hay columna de lane/flowcell/
  run/secuenciación separada de `gem_group`. Si existe, no está en `obs` de este
  archivo.
- **(c) guía distinta al mismo gen:** **sí, resuelta.** `sgID_AB` (2679/2662 valores)
  es mayor que `gene` (2394/2394) en los dos archivos — hay genes con más de un par de
  guías. La comparación correcta es `(gene, sgID_AB)`: mismo `gene` + distinto `sgID_AB`
  = guías distintas al mismo blanco → **no es réplica, es benchmark** (regla (c) de
  `PLAN.md`).

---

## 3 · `/datos/ref/replogle/{k562_essential,k562_gwps,rpe1}_raw_sc.h5ad`

Mismo esquema que Nadig (misma convención de pipeline/lab):

| archivo | n_obs | `gem_group` únicos | `gene` únicos | `sgID_AB` únicos |
|---|---|---|---|---|
| `k562_essential_raw_sc.h5ad` | 310,385 | 48 | 2058 | 2273 |
| `k562_gwps_raw_sc.h5ad` | 1,989,578 | 267 | 9867 | 11187 |
| `rpe1_raw_sc.h5ad` | 247,914 | 56 | 2394 | 2662 |

Columnas (idénticas a las de Nadig más una): `cell_barcode`, `UMI_count`, `mitopercent`,
**`gem_group`**, `z_gemgroup_UMI`, `gene`, `gene_id`, `gene_transcript`, **`sgID_AB`**,
`transcript`, y además:

| columna | dtype | candidata | notas |
|---|---|---|---|
| **`core_adjusted_UMI_count`** | `float32` | sí | UMI ajustado — probablemente por `gem_group`/`core_scale_factor`, no un eje nuevo |
| **`core_scale_factor`** | `float32` | sí | nº de únicos = nº de `gem_group` en cada archivo (48/267/56) — **es un factor por `gem_group`**, no columna independiente |

5 ejemplos de `sgID_AB` (k562_gwps, el más grande): `A1BG_+_58858964.23-P1\|A1BG_-_58858788.23-P1`,
`A1BG_-_58864840.23-P2\|A1BG_-_58864822.23-P2`, `AAAS_-_53715438.23-P1P2\|AAAS_+_53715355.23-P1P2`,
`AACS_+_125549983.23-P1P2\|AACS_-_125550169.23-P1P2`, `AAED1_-_99417574.23-P1P2\|AAED1_+_99417525.23-P1P2`.
Nótese `A1BG` con dos entradas (`-P1`, `-P2`) — **dos guías distintas al mismo gen**,
exactamente el caso (c).

**(a)/(b)/(c) para Replogle: mismas conclusiones que Nadig** — `gem_group` es la única
candidata a (a), (b) no tiene columna propia, (c) se resuelve con `(gene, sgID_AB)`.
`core_scale_factor` confirma que `core_scale_factor` y (probablemente)
`core_adjusted_UMI_count` son derivados de `gem_group`, no ejes nuevos — su nº de
únicos coincide exactamente con el nº de `gem_group`.

---

## 4 · `/datos/ref/replogle/{k562_essential,k562_gwps,rpe1}_raw_bulk.h5ad`

Ya son pseudobulk — **una fila por perturbación**, no por célula. `gene_transcript` es
único en las tres (`nunique == n_obs`: 2285/11258/2679). **No hay `gem_group` ni
`sgID_AB`** — la agregación a bulk ya colapsó ese eje. Columnas completas:
`TE_ratio`, `UMI_count_unfiltered`, `anderson_darling_counts`, `cnv_score_z`,
`control_expr`, `core_control` (bool), `energy_test_p_value`, `fold_expr`,
`gene_transcript`, `mann_whitney_counts`, `mean_leverage_score`, `mitopercent`,
`num_cells_filtered`, `num_cells_unfiltered`, `pct_expr`, `std_leverage_score`,
`z_gemgroup_UMI` (residual del nombre, pero aquí es un valor por perturbación, no por
célula/gem_group).

**Para T13: esta fuente no sirve** — no tiene la granularidad celular que (a)/(b)/(c)
necesitan. Solo el `_sc` correspondiente (§3) la tiene.

---

## 5 · `/datos/ref/orion/{hct116,hek293t}.parquet` — grano real

Esquema leído de los metadatos del parquet (`pyarrow.parquet.ParquetFile.schema_arrow`
+ `metadata.num_rows`, **sin leer ninguna fila**):

| archivo | filas | columnas |
|---|---|---|
| `hct116.parquet` | 102,228,672 | `target, feature, target_mean, ref_mean, target_ncells, ref_ncells, log2_fold_change, p_value, Ueffect, p_adj` |
| `hek293t.parquet` | 53,646,719 | (idénticas) |

**Grano real: `(target, feature)`** — una fila por par (gen perturbado, gen medido),
ya como tabla de DE agregada (incluye `target_ncells`/`ref_ncells`, es decir, ya pasó
por pseudobulk). **No hay columna de pozo, transducción, guía ni lane** — ese nivel se
perdió en la agregación. Orion **no aporta a T13**; es insumo para las métricas de
fidelidad/dirección (familia `fid`/`reach`), no para contar réplicas.

---

## 6 · `/datos/ref/feng/*.tsv.gz` — grano real

| archivo | columnas | grano |
|---|---|---|
| `GenomeWideScreen_LFC_byGene.tsv.gz` | `Target, Expressed_Gene_Symbol, Expressed_Gene_Ens_ID, lfc, pval_adj` | `(Target, Expressed_Gene)` |
| `TargetedScreen_LFC_byGene.tsv.gz` | igual | `(Target, Expressed_Gene)` |
| `TargetedScreen_LFC_byGene-perLine.tsv.gz` | `Target, Expressed_Gene_Symbol, Expressed_Gene_Ens_ID, Cell_Line, wt_expr, lfc, pval_adj` | `(Target, Expressed_Gene, Cell_Line)` |

Igual que Orion: ya es una tabla de LFC agregada por gen perturbado × gen medido (y por
línea en la versión `-perLine`). **Sin columna de pozo/transducción/guía** — no aporta
a T13.

---

## 7 · `/datos/ref/otros/{GSE225775,GSE139944}_RAW.tar` — grano real, sin extraer

Se listó el contenido (`tar -tvf`, sin extraer) y se leyeron 2–3 archivos de anotación
por `tar -xO -f <tar> <miembro> | zcat | head`, sin desempaquetar el `.tar` completo.

**`GSE225775_RAW.tar`** (31 miembros) — **grano real: un `GSM` (muestra GEO) = una
corrida/sublibrería sci-Plex-GxE**. Dentro de cada `GSM`:
- `*_cell.annotations.txt.gz`: sin cabecera, 2 columnas — `(barcode_RT+Lig, etiqueta_de_muestra)`,
  ej. `01B_A01_RT_BC_10_Lig_BC_12 → sci3_A172_MMR_HPRT1_CROPseq`. Es una etiqueta de
  **pool/sublibrería**, no una columna de pozo per se.
- `*_hashTable.out.txt.gz`: sin cabecera, 5 columnas — `(etiqueta_muestra, barcode,
  etiqueta_hash_oligo, count, rank)`. La tercera columna **sí codifica el pozo** como
  substring, ej. `plate16_H2_A172_CRISPRi_HPRT1_0_dmso_96` → placa 16, pozo H2, línea
  A172, condición. **No es una columna estructurada de "well"; es un string a parsear.**
- `*_sgRNATable_out.rds.gz` / `*_gRNATable_reads.out.txt.gz`: existen (nivel de guía),
  **no se abrieron** — quedan fuera de este paso.

**`GSE139944_RAW.tar`** (32 miembros) — **grano real: un `GSM` = una corrida
sci-Plex**. `*_cell.annotations.txt.gz` igual que arriba (barcode → etiqueta de
muestra, ej. `sciChem19`). `*_pData.txt.gz` (volcado de `data.frame` de R, separado por
espacio y comillas) trae `Cell, sample, Size_Factor, n.umi, hash_umis, pval, qval,
top_to_second_best_ratio, top_oligo` — **`top_oligo` codifica el pozo+línea** tras
deconvolución de hash, ej. `"NIH3T3_E07"`.

**Para (a)/(b)/(c) en estos dos GSE: no resuelto.** El pozo vive embebido en un string
de hash (`hashTable`/`top_oligo`), no en una columna categórica limpia, y las tablas de
guía no se abrieron — no alcanza para distinguir pozo-aparte de resecuenciado, ni para
confirmar guía-distinta-al-mismo-gen en estos dos GSE. Queda abierto.

---

## 8 · `/datos/ref/otros/*.rds` (9 archivos) — grano real

**No se pudo inspeccionar: no hay R instalado en esta máquina** (`which R`/`Rscript`
no devolvieron nada). Son objetos Seurat, formato binario de R (serialización `RDS`),
no legible con herramientas de Python sin un puente (`rds2py`, `pyreadr` con
limitaciones, o un `Rscript` real). Por el `A_readme.txt` del propio depósito (Zenodo
14518762), el grano **esperado** es una célula por fila en el assay de cada objeto
(`Seurat_object_{IFNB,IFNG,TNFA,TGFB1,INS}_Perturb_seq.rds`), más los genelists y el
bulk RNA-seq — pero esto es lo que **dice el README**, no lo que confirmé abriendo el
archivo. Queda abierto; cerrarlo cuesta instalar R (o un puente RDS→Python) y sí entra
en el alcance de "costaría cerrarlo" de `registro/limitaciones.md` §4.

---

## 9 · Resumen para T13

| fuente | sirve para T13 (nivel célula) | (a) pozo/transducción | (b) resecuenciado | (c) guía≠gen |
|---|---|---|---|---|
| `reto/*.h5ad` | ❌ solo controles | — | — | — |
| `nadig/*_sc.h5ad` | ✅ | `gem_group` (sin confirmar transducción) | no encontrada | ✅ `(gene, sgID_AB)` |
| `replogle/*_sc.h5ad` | ✅ | `gem_group` (ídem) | no encontrada | ✅ `(gene, sgID_AB)` |
| `replogle/*_bulk.h5ad` | ❌ ya agregado | — | — | — |
| `orion/*.parquet` | ❌ ya agregado | — | — | — |
| `feng/*.tsv.gz` | ❌ ya agregado | — | — | — |
| `otros/GSE*.tar` | parcial, sin abrir tablas de guía | pozo embebido en string de hash, no columna | no resuelto | no resuelto |
| `otros/*.rds` | sin inspeccionar (sin R) | — | — | — |
