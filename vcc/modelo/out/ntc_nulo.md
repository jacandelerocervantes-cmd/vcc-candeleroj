# El nulo de `ntc_id` en contexto, y las tres verificaciones de validez

## 0. Que es este reporte y que no

Esto PRESENTA DATOS. No lleva veredicto, no concluye, y no va dirigido a
ningun candidato siguiente. Cada construccion se declara **antes** de la
medicion que la usa.

**La cantidad que se mide.** Cada `context_{A..F}.h5ad` trae 18,400 celulas
NO DIRIGIDAS y la columna `ntc_id` nombra cual de 46 guias lleva cada una,
400 por guia. Esa es estructura de replica DENTRO del control, y existe en
A/B/C y en D/E/F. Partir las 46 guias en dos mitades y correr la puerta de
significancia del preset `vcc2026` mitad contra mitad da la distribucion
NULA de expresion diferencial aparente entre dos grupos de celulas no
dirigidas, en el contexto que se entrega, con su profundidad y su `n`. Es
independiente de la convencion del marcador **por construccion**: son
controles contra controles, no hay perturbacion que acertar.

**Lo que se leyo**: el documento de formatos del reto
(`modelo/out/paso1_contrato.md`), `modelo/out/control_permutado.md` §6-§7,
y el codigo de la puerta (`cell_eval2/configs/vcc2026.yaml`,
`cell_eval2/de_compute.py`, `cell_eval2/de.py`).

**Lo que NO se toco**: el lado perturbado -- no existe en 2026, y el de
2025/Jiang no entra aqui; de Jiang solo se leyeron las celulas `NT`. No se
leyo nada de `~/vcc/benchmark/` ni de `/datos/benchmark_VEDADO_*`. Sin
envios, sin `vcc sample`, sin commits, y `vcc prep` no se corrio en absoluto.

**La puerta, leida del propio preset** (ningun numero se escribio a mano:
todos salen de `EvalConfig.from_preset("vcc2026")`):

| parametro | valor |
|---|---|
| `control` | `non-targeting` |
| `de.backend` | `auto` |
| `de.clip_value` | `None` |
| `de.epsilon` | `1e-09` |
| `de.fdr_scope` | `per_pert` |
| `de.mean_calc` | `arithmetic` |
| `de.method` | `wilcoxon` |
| `de.min_abs_log2fc` | `0.0` |
| `de.nan_lfc_policy` | `mask` |
| `de.p_adj_threshold` | `0.05` |
| `de.sort_by` | `abs_log2_fold_change` |
| `filter_gene_min_cpm_cell` | `5.0` |
| `input_type` | `counts` |
| `num_threads` | `-1` |
| `pert_col` | `target` |
| `preset` | `vcc2026` |
| `target_sum` | `1000000.0` |
| `version` | `v2` |

El camino, el mismo que recorre el marcador en su lado de CPU:
`compute_de` (normalize_total a `target_sum` -> medias aritmeticas -> log2FC
con `epsilon` -> `log1p` -> MWU del backend -> `_apply_cpm_filter`, que
conserva los genes cuyo **CPM medio en el grupo de REFERENCIA** supera
`filter_gene_min_cpm_cell` y **recalcula BH por blanco** sobre las filas que
quedan) -> `prep_de_side` (`nan_lfc_policy='mask'`, `min_abs_log2fc=0.0`,
que es no-op) -> conjunto DE = filas con `p_adj < 0.05`.

**Verificacion de que el arnes CORRE la puerta y no la reimplementa.** Se
comprobo, sobre una corrida (contexto D, semilla `20274300`, 400 v 400), que
las cuatro piezas que definen el conjunto DE salen exactas
(`modelo/out/ntc_nulo/verif_puerta.log`):

| # | lo que define la puerta | esperado | obtenido | igual |
|---|---|---|---|---|
| 1 | conjunto probado = genes con CPM medio de la REFERENCIA > 5.0 | 9,619 | 9,619 | si |
| 2 | `p_adj` = Benjamini-Hochberg de `p_value` DENTRO de ese conjunto | — | `max|dif| = 0.000e+00` | si |
| 3 | conjunto DE = filas con `p_adj < 0.05` | 0 | 0 | si |
| 4 | `log2FC = log2((media_CPM_blanco + eps)/(media_CPM_ref + eps))` | — | `max|dif| = 1.289e-15` | si |

🔴 **El backend resuelto aqui es `pdex`, no `gpudge`**: `de.backend='auto'`
y la maquina no tiene CUDA. El propio `cell_eval2` lo avisa ("DE numbers
differ between engines"). Todas las cifras de este reporte son de `pdex`.

---

## PARTE 0 — Las tres verificaciones de validez

Las tres se hicieron **primero** y son bloqueantes. Ninguna se disparo.

### 0.a — `controls-final` esta bajado y completo

**TABLA 0.1 — los 14 archivos de `/datos/reto`, con `sha256` de cada uno.**

| # | particion | archivo | bytes | sha256 |
|---|---|---|---|---|
| 1 | `val` | `context_A.h5ad` | 224,973,316 | `f22e71968487d3da9402769dcd0c73a47851fad4cd7dd25510283c923031fe3e` |
| 2 | `val` | `context_B.h5ad` | 210,966,111 | `41557555b3febf6db19fd930b4dd7db88941840077bf985d14dfff959c782380` |
| 3 | `val` | `context_C.h5ad` | 226,056,691 | `e090e7c7cfdd0006f689214ba38ba0cd2669578e89c9bc325b65f0bd940c02e0` |
| 4 | `val` | `gene_names.csv` | 119,295 | `25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201` |
| 5 | `val` | `manifest.json` | 715 | `2c2b5d425df817ca35e9edf08282d12494ce42913ba824375bf0f98a1cc676f5` |
| 6 | `val` | `pert_counts.csv` | 1,906 | `f57edd7b912ebd718efc7ee9d0f334772513e7cc418d133ce525470e373b3276` |
| 7 | `val` | `vcc_2026_controls.zip` | 662,118,680 | `329a22bc29cac4f43ad0b846b9bb1a383f7a945c83a8cb0450c465f6dc5cc4b5` |
| 8 | `final` | `context_D.h5ad` | 217,126,334 | `97dfff32570748677fef6ee732955f71af5e7596a3747a11aac36903206a0d13` |
| 9 | `final` | `context_E.h5ad` | 213,808,470 | `1942a0b6e00331da4503a9f67101a002564675bc775c2b8534a9628941436e93` |
| 10 | `final` | `context_F.h5ad` | 212,299,193 | `5d835a036045fdf3849b5c3f8633282498aba75038ef77c999e596815a2ee38c` |
| 11 | `final` | `gene_names.csv` | 119,295 | `25bfa66715e186bebabce7ac788bbcea47e2bf59ca70be1f8f3a06f2f0e47201` |
| 12 | `final` | `manifest.json` | 717 | `7fd5a2ea490c40b6d156dd189db2c8cbddedabab13dec1c27b69089a04b34e28` |
| 13 | `final` | `pert_counts.csv` | 1,913 | `6a472047d58ddb2a59b2345f962dabaf4120276c381a5e97ec19ea03f640b118` |
| 14 | `final` | `vcc_2026_controls_final.zip` | 643,356,568 | `d625f47f56a6b582d009f88ed5ca9473de0a19260c36cdf0a71c6bb67a59b2e9` |

Los tres `final/context_{D,E,F}.h5ad`, `final/gene_names.csv` y
`final/pert_counts.csv` estan los cinco, mas el `vcc_2026_controls_final.zip`
de origen. `val/` esta completo igual.

**Y el contenido de los dos `.zip` de origen, listado**: cada uno trae
**exactamente 6 archivos** y los tamanos coinciden byte a byte con los ya
extraidos (p. ej. `context_D.h5ad` = 217,126,334 en el zip y en disco):

| zip | fecha | archivos | contenido |
|---|---|---|---|
| `val/vcc_2026_controls.zip` | 2026-08-19 | 6 | `context_{A,B,C}.h5ad`, `gene_names.csv`, `pert_counts.csv`, `manifest.json` |
| `final/vcc_2026_controls_final.zip` | 2026-10-09 | 6 | `context_{D,E,F}.h5ad`, `gene_names.csv`, `pert_counts.csv`, `manifest.json` |

🔴 **Ni el zip ni el directorio traen una sola celula perturbada.** Los 6
archivos de cada lado son los tres contextos de CONTROL mas los tres de panel.
Esa es la verificacion literal de la afirmacion que la seccion 2.4 necesita.

**TABLA 0.2 — forma y conteos de los seis contextos, leidos del h5ad.**

| # | ctx | particion | celulas | genes | nnz almacenado | `nnz==indptr[-1]` | `target_gene` | `context` | `ntc_id` cats | celulas/guia | suma |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **A** | `val` | 18,400 | 18,533 | 109,910,478 | si | `non-targeting` | `A` | 46 | 400 en las 46 | 18,400 |
| 2 | **B** | `val` | 18,400 | 18,533 | 101,546,777 | si | `non-targeting` | `B` | 46 | 400 en las 46 | 18,400 |
| 3 | **C** | `val` | 18,400 | 18,533 | 108,213,983 | si | `non-targeting` | `C` | 46 | 400 en las 46 | 18,400 |
| 4 | **D** | `final` | 18,400 | 18,533 | 104,487,981 | si | `non-targeting` | `D` | 46 | 400 en las 46 | 18,400 |
| 5 | **E** | `final` | 18,400 | 18,533 | 101,608,841 | si | `non-targeting` | `E` | 46 | 400 en las 46 | 18,400 |
| 6 | **F** | `final` | 18,400 | 18,533 | 100,855,632 | si | `non-targeting` | `F` | 46 | 400 en las 46 | 18,400 |

Los 46 `ntc_id` son **el mismo conjunto en los seis** archivos: `si`. 46 x 400 = 18,400 en cada contexto, sin
excepcion, y `18400 = ntc_id_suma` cierra el conteo.

**`manifest.json` esta JUNTO a `pert_counts.csv`** en las dos particiones --
que es de donde `vcc prep` lee los contextos:

| particion | `manifest.json` existe | `pert_counts.csv` existe | mismo directorio |
|---|---|---|---|
| `final` | si | si | si |
| `val` | si | si | si |

Y lo que dice el de `final/`, verbatim en sus campos:

| campo | valor |
|---|---|
| `cells_per_pert` | `400` |
| `context_col` | `"context"` |
| `contexts` | `["D", "E", "F"]` |
| `control_label` | `"non-targeting"` |
| `n_constructs` | `300` |
| `n_genes` | `18533` |
| `panel_id` | `"vcc2026-test-1"` |
| `partition` | `"test"` |
| `per_context.D` | `{"control_cells": 18400, "ground_truth_cells": 138400, "n_ntc_ids": 46, "n_perturbations": 300}` |
| `per_context.E` | `{"control_cells": 18400, "ground_truth_cells": 138400, "n_ntc_ids": 46, "n_perturbations": 300}` |
| `per_context.F` | `{"control_cells": 18400, "ground_truth_cells": 138400, "n_ntc_ids": 46, "n_perturbations": 300}` |
| `pert_col` | `"target_gene"` |
| `season` | `"2026"` |

🔴 `ground_truth_cells = 138,400` por contexto y `control_cells = 18,400`:
**138,400 = 300 x 400 + 18,400**. El lado real contra el que se mide una
perturbacion es el pool de control de 18,400 celulas, no otras 400. Esa
aritmetica es la razon del tercer brazo de `n` de la Parte 1.

### 0.b — Cobertura real POR SIMBOLO contra gwps

🔴 **Correccion de lectura, declarada.** `gene_names.csv` trae **cabecera**
(`gene_name`), asi que contarlo por lineas da 18,534. Los simbolos son
**18,533**, y la cuenta de abajo los lee con `csv.reader` saltando la
cabecera. Se cruzo ademas contra `var/_index` del h5ad, elemento por
elemento.

gwps = `/datos/ref/replogle/k562_gwps_raw_bulk.h5ad`. Dos conjuntos
distintos que no hay que confundir: los genes **MEDIDOS** (las 8,248
columnas de `var`) y los blancos con **PERFIL** (las filas de `obs`,
nombradas `*_GEN_CONSTRUCCION`).

| cantidad de gwps | valor |
|---|---|
| filas totales | 11,258 |
| filas `core_control` | 514 |
| columnas de `var` | 8,248 |
| nombres de gen unicos | 8,246 |
| nombres duplicados (se excluyen enteros) | `HSPA14`, `TBCE` |
| blancos con perfil, unicos | 9,867 |

**TABLA 0.3 — cobertura por simbolo, las dos particiones.**

| # | cantidad | `val` | `final` |
|---|---|---|---|
| 1 | simbolos en `gene_names.csv` | 18,533 | 18,533 |
| 2 | simbolos unicos | 18,533 | 18,533 |
| 3 | `n_vars` del h5ad | 18,533 | 18,533 |
| 4 | orden identico a `var/_index`, elemento por elemento | si | si |
| 5 | **de los 18,533: cuantos estan MEDIDOS en gwps** | 7,681 | 7,681 |
| 6 | idem, excluyendo los 2 nombres duplicados | 7,679 | 7,679 |
| 7 | de los 18,533: cuantos NO estan en gwps | 10,852 | 10,852 |
| 8 | cobertura % | 41.4450 | 41.4450 |
| 9 | perturbaciones en `pert_counts.csv` | 300 | 300 |
| 10 | **de las 300: cuantas tienen PERFIL en gwps** | 272 | 277 |
| 11 | de las 300: cuantas NO tienen perfil | 28 | 23 |
| 12 | las 300 estan todas en `gene_names.csv` | si | si |

**Las 28 perturbaciones de `val/pert_counts.csv` SIN perfil en gwps**
(lista completa): `ABCD1`, `CAPRIN2`, `EPHB2`, `HEATR6`, `KIF21B`, `MAPK7`, `MLKL`, `NICN1`, `PARP3`, `PBLD`, `PDK3`, `PHF11`, `PSMB9`, `RAB11FIP5`, `SEMA4F`, `SLC44A1`, `SNN`, `SOCS5`, `SPATA20`, `TAF4`, `TBC1D19`, `TDRD7`, `TESK2`, `TMEM104`, `TMEM79`, `TTBK2`, `TXNDC16`, `ZC2HC1A`.

**Las 23 perturbaciones de `final/pert_counts.csv` SIN perfil en gwps**
(lista completa): `ANG`, `ARSA`, `ATP10D`, `CAPRIN2`, `HCFC2`, `HEATR6`, `ITPRIPL1`, `LYSMD4`, `MAPK7`, `MLKL`, `NEK3`, `PARP3`, `PDK3`, `PSMB9`, `RAB26`, `SEMA4F`, `SIK2`, `STK10`, `TBC1D19`, `TESK2`, `TMEM79`, `UBTD2`, `ZDHHC13`.

---

### 0.c — El `nnz` de las celulas de control, y el Poisson aditivo

**El tope y de donde sale.** `HARD_CAP_NNZ = 4,750,000,000` entradas
almacenadas por entrega (`vcc/sizing.py` del CLI instalado, ya confirmado
leyendo el paquete en `paso1_contrato.md` punto 42). Una entrega lleva
`300 x 400 x 3 = 360,000` filas, asi que el presupuesto por celula es
`4.75e9 / 360,000 = 13,194.4444`.

🔴 **Cual de las dos cifras ATA.** El tope es sobre el TOTAL, no por
celula. `13,194.44` es el presupuesto por celula **si se gasta parejo**, asi
que la cantidad que ata es la **MEDIA** (total = media x 360,000); el maximo
de una celula sola puede pasarse de 13,194 sin violar nada mientras el total
aguante. Las dos columnas de margen van separadas por eso.

**TABLA 0.4 — `nnz` por celula de control, por contexto. 18,400 celulas cada uno.**

| # | ctx | celulas | nnz media | nnz mediana | nnz p99 | nnz max | nnz min | nnz sd | frac tope (media) | margen (media) | frac tope (max) | margen (max) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **A** | 18,400 | 5,973.4 | 6,147.0 | 8,003.0 | 8,453 | 2,205 | 1,185.5 | 0.4527 | x2.209 | 0.6406 | x1.561 |
| 2 | **B** | 18,400 | 5,518.8 | 5,756.0 | 7,178.0 | 7,649 | 580 | 1,189.4 | 0.4183 | x2.391 | 0.5797 | x1.725 |
| 3 | **C** | 18,400 | 5,881.2 | 6,006.0 | 7,944.0 | 8,482 | 1,937 | 1,144.4 | 0.4457 | x2.243 | 0.6428 | x1.556 |
| 4 | **D** | 18,400 | 5,678.7 | 5,891.5 | 7,673.0 | 8,146 | 787 | 1,255.4 | 0.4304 | x2.323 | 0.6174 | x1.620 |
| 5 | **E** | 18,400 | 5,522.2 | 5,634.0 | 7,606.0 | 8,107 | 1,607 | 1,140.5 | 0.4185 | x2.389 | 0.6144 | x1.628 |
| 6 | **F** | 18,400 | 5,481.3 | 5,630.0 | 7,423.0 | 7,914 | 1,571 | 1,107.7 | 0.4154 | x2.407 | 0.5998 | x1.667 |

**El margen, dicho como lo gasta una entrega.** Una entrega que reusa filas
de control tal cual (el contrato ya probado) gasta, en esperanza,
`360,000 x media_nnz` de los tres contextos que use:

| particion | contextos | media nnz/celula de los 3 | nnz total previsto | frac del tope | margen |
|---|---|---|---|---|---|
| `val` | A, B, C | 5,791.15 | 2,084,812,422 | 0.4389 | x2.278 |
| `final` | D, E, F | 5,560.73 | 2,001,863,830 | 0.4214 | x2.373 |

🔴 **Cruce independiente.** `paso1_contrato.md` punto 42 MIDIO el nnz real de
una carga trivial sobre `val` (A/B/C) construida de verdad:
**2,084,329,279**. La proyeccion de arriba desde la media de `nnz` por celula
da **2,084,812,422** -- difieren en
0.023 %. Son dos caminos
distintos a la misma cifra (uno construyo la matriz, el otro suma `indptr`).

#### El Poisson aditivo de `cand01_covar._reescala_entero`

**La construccion, declarada.** `_reescala_entero(bloque, rat, rng)` reescala
un bloque CSR de cuentas enteras por el factor por gen `rat`:

```
r_h <= 1:  nuevo_h ~ Binomial(x_h, r_h)            (adelgazamiento)
r_h  > 1:  nuevo_h = x_h + Poisson(x_h * (r_h-1))  (engorde)
```

Se aplico a una muestra de **400 celulas de control** del contexto
**D**, sorteadas sin reemplazo con semilla fija
`20274999`, con las **dos** emisiones de densidad que se
venian usando, y tres blancos cada una:

* `cand01_covar.construye_respuesta` -- ratios estimados EN CONTEXTO sobre el
  pool de 18,400 celulas del propio contexto (`r = 1 - rho`, acotado en [0,2]).
* `gwps_emision.construye_ratios` -- ratios EXTERNOS del perfil de k562_gwps
  (`r = 1 + s`, `s` del rango normalizado, acotado en [0,2]).

El bloque de partida tiene **2,232,478** entradas almacenadas,
**5,581.2** por celula.

**TABLA 0.5 — los seis brazos, `nnz` antes y despues. ES VALIDEZ, NO ESTILO.**

| # | emision | blanco | genes `r>1` | genes `r<1` | genes `r==1` | `r` min | `r` max | nnz despues | **delta nnz** | nnz/celula media | nnz max de una celula | frac tope (max) | **revienta el tope** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | cand01 (en contexto) | `ABCC1` | 3,335 | 4,610 | 10,588 | 0.0000 | 1.1994 | 2,215,403 | **-17,075** | 5,538.5 | 8,115 | 0.6150 | no |
| 2 | cand01 (en contexto) | `ACLY` | 3,080 | 6,464 | 8,989 | 0.0000 | 1.2929 | 2,181,770 | **-50,708** | 5,454.4 | 8,062 | 0.6110 | no |
| 3 | cand01 (en contexto) | `ADNP` | 1,230 | 6,256 | 11,047 | 0.0000 | 1.1297 | 2,202,676 | **-29,802** | 5,506.7 | 8,116 | 0.6151 | no |
| 4 | gwps (externa) | `ABCC1` | 2,003 | 5,676 | 10,854 | 0.0000 | 2.0000 | 1,813,867 | **-418,611** | 4,534.7 | 7,147 | 0.5417 | no |
| 5 | gwps (externa) | `ACLY` | 3,101 | 4,578 | 10,854 | 0.0000 | 2.0000 | 1,892,506 | **-339,972** | 4,731.3 | 7,279 | 0.5517 | no |
| 6 | gwps (externa) | `ADNP` | 2,348 | 5,331 | 10,854 | 0.0000 | 1.9996 | 1,831,501 | **-400,977** | 4,578.8 | 7,144 | 0.5414 | no |

**Y la suma de cuentas, que es otra cantidad** (el tope de `nnz` no la mira,
pero la clausula de 1,000,000 por celula si):

| # | emision | blanco | suma antes | suma despues | razon |
|---|---|---|---|---|---|
| 1 | cand01 | `ABCC1` | 8,156,733 | 8,252,150 | 1.0117 |
| 2 | cand01 | `ACLY` | 8,156,733 | 7,593,674 | 0.9310 |
| 3 | cand01 | `ADNP` | 8,156,733 | 7,833,194 | 0.9603 |
| 4 | gwps | `ABCC1` | 8,156,733 | 6,446,926 | 0.7904 |
| 5 | gwps | `ACLY` | 8,156,733 | 7,454,481 | 0.9139 |
| 6 | gwps | `ADNP` | 8,156,733 | 6,642,554 | 0.8144 |

**La cota, y por que es estructural y no un ajuste.** `_reescala_entero`
opera sobre `b.data` -- las entradas YA GUARDADAS del CSR -- y luego llama
`eliminate_zeros()`. El brazo Poisson escribe `x + Poisson(...)` **solo donde
`x` ya era distinto de cero**, asi que no puede crear un no-cero nuevo; el
brazo binomial si puede destruirlos. De ahi `nnz_despues <= nnz_antes`,
siempre, sea cual sea `rat`.

Medido: el `delta nnz` es **negativo en los 6 de 6 brazos** (`True`), de
**-418,611** a
**-17,075**. El `nnz` maximo de una
celula despues del reescalado llega a
**8,116**, que es
**0.6151** del presupuesto por
celula. **El Poisson aditivo no lo revienta**, y no por poco: no puede
subirlo en absoluto.

**Resultado de la Parte 0: las tres pasaron. No se disparo ninguna parada.**

---

## PARTE 1 — El nulo, por contexto

### 1.1 La construccion, declarada antes de medir

* **La unidad es la GUIA, no la celula.** Las 400 celulas de una guia van
  enteras a un lado. Partir celulas DENTRO de una guia subestimaria el nulo
  porque las dos mitades compartirian el efecto de esa guia; aqui no se hace
  nunca, y se verifica particion por particion.
* **Las 46 guias se parten en 23 contra 23.** `rng.permutation(46)` sobre las
  46 categorias de `ntc_id` ORDENADAS alfabeticamente; los primeros 23 son la
  mitad de **REFERENCIA**, los ultimos 23 la mitad **BLANCO**.
* **Semilla**: `20274000 + 100*idx_contexto + j`, con
  `idx` A=0, B=1, C=2, D=3, E=4, F=5, y `j = 0..7` -- **ocho particiones por**
  **contexto**, 48 en total.
* **La puerta es la del preset `vcc2026`**, la misma que usa el marcador, con
  cada parametro leido del preset (tabla de la seccion 0).

**Los tres brazos de `n`, declarados aparte y no mezclados:**

| brazo | `n` referencia | `n` blanco | que es |
|---|---|---|---|
| `completo` | 9,200 | 9,200 | las 23 guias enteras de cada lado |
| `n400_400` | 400 | 400 | el `n` POR PERTURBACION del reto, a los dos lados. **Es el brazo que pide el encargo.** |
| `n400_ref_completa` | 9,200 | 400 | 🔴 **AÑADIDO y declarado.** El reto no compara 400 contra 400: compara las 400 celulas de una perturbacion contra el pool de control del contexto (`manifest.json`: `138,400 = 300x400 + 18,400`). Este brazo es el que se parece al regimen de la entrega. |

El submuestreo a 400 es **sin reemplazo** sobre las 9,200 celulas de la mitad,
con el mismo `rng` de la particion y en orden declarado (referencia primero,
blanco despues). **No parte guias**: las guias siguen enteras de un lado, y el
submuestreo solo adelgaza celulas dentro de la mitad que ya les toco.

🔴 **La asimetria de la puerta, declarada como propiedad, no como ajuste.** El
filtro de genes del preset (`filter_gene_min_cpm_cell = 5.0`) lo decide **el
grupo de REFERENCIA solo** (`_apply_cpm_filter`, #351), y BH se recalcula sobre
los genes que quedan. Asi que el conjunto probado `m` depende de que mitad hace
de referencia. Por eso `m` varia entre particiones y se reporta en cada corrida.

### 1.2 Las 48 particiones: ninguna parte una guia

**TABLA 1.1 — las 48 particiones, verificadas una por una.** `guia partida` es
falso en las 48; los conjuntos de guias son disjuntos y su union son las 46; y
las 46 guias tienen 400 celulas cada una a los dos lados.

| # | ctx | j | semilla | guias ref | guias blanco | disjuntas | union = 46 | celulas ref | celulas blanco | 400/guia ref | 400/guia blanco | **guia partida** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **A** | 0 | `20274000` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 2 | **A** | 1 | `20274001` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 3 | **A** | 2 | `20274002` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 4 | **A** | 3 | `20274003` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 5 | **A** | 4 | `20274004` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 6 | **A** | 5 | `20274005` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 7 | **A** | 6 | `20274006` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 8 | **A** | 7 | `20274007` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 9 | **B** | 0 | `20274100` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 10 | **B** | 1 | `20274101` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 11 | **B** | 2 | `20274102` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 12 | **B** | 3 | `20274103` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 13 | **B** | 4 | `20274104` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 14 | **B** | 5 | `20274105` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 15 | **B** | 6 | `20274106` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 16 | **B** | 7 | `20274107` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 17 | **C** | 0 | `20274200` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 18 | **C** | 1 | `20274201` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 19 | **C** | 2 | `20274202` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 20 | **C** | 3 | `20274203` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 21 | **C** | 4 | `20274204` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 22 | **C** | 5 | `20274205` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 23 | **C** | 6 | `20274206` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 24 | **C** | 7 | `20274207` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 25 | **D** | 0 | `20274300` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 26 | **D** | 1 | `20274301` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 27 | **D** | 2 | `20274302` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 28 | **D** | 3 | `20274303` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 29 | **D** | 4 | `20274304` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 30 | **D** | 5 | `20274305` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 31 | **D** | 6 | `20274306` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 32 | **D** | 7 | `20274307` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 33 | **E** | 0 | `20274400` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 34 | **E** | 1 | `20274401` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 35 | **E** | 2 | `20274402` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 36 | **E** | 3 | `20274403` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 37 | **E** | 4 | `20274404` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 38 | **E** | 5 | `20274405` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 39 | **E** | 6 | `20274406` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 40 | **E** | 7 | `20274407` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 41 | **F** | 0 | `20274500` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 42 | **F** | 1 | `20274501` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 43 | **F** | 2 | `20274502` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 44 | **F** | 3 | `20274503` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 45 | **F** | 4 | `20274504` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 46 | **F** | 5 | `20274505` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 47 | **F** | 6 | `20274506` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |
| 48 | **F** | 7 | `20274507` | 23 | 23 | si | si | 9,200 | 9,200 | si | si | **no** |

Ejemplo de una particion escrita entera, para que se pueda reproducir
(contexto A, j=0, semilla `20274000`):

* mitad REFERENCIA (23 guias): `non-targeting-12`, `non-targeting-21`, `non-targeting-22`, `non-targeting-24`, `non-targeting-3`, `non-targeting-36`, `non-targeting-37`, `non-targeting-39`, `non-targeting-40`, `non-targeting-44`, `non-targeting-52`, `non-targeting-59`, `non-targeting-65`, `non-targeting-68`, `non-targeting-70`, `non-targeting-73`, `non-targeting-77`, `non-targeting-8`, `non-targeting-80`, `non-targeting-82`, `non-targeting-85`, `non-targeting-89`, `non-targeting-98`
* mitad BLANCO (23 guias): `non-targeting-11`, `non-targeting-13`, `non-targeting-15`, `non-targeting-17`, `non-targeting-20`, `non-targeting-25`, `non-targeting-28`, `non-targeting-29`, `non-targeting-30`, `non-targeting-45`, `non-targeting-53`, `non-targeting-55`, `non-targeting-58`, `non-targeting-61`, `non-targeting-7`, `non-targeting-71`, `non-targeting-72`, `non-targeting-74`, `non-targeting-83`, `non-targeting-88`, `non-targeting-91`, `non-targeting-93`, `non-targeting-99`

### 1.3 🔴 Cuantos genes llama significativos, y su distribucion sobre particiones

**TABLA 1.2 — resumen por contexto y brazo, sobre las 8 particiones.**
`rango` es `[min, max]` de `n_sig` sobre las 8 particiones de ese contexto.

| # | ctx | brazo | corridas | **`n_sig` media** | **`n_sig` sd** | mediana | **rango [min, max]** | corridas con `n_sig>0` | `m` probados [min, max] | `p_crit = 0.05/m` | `p_min` mediana | genes llamados (union) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **A** | `completo` 9,200 v 9,200 | 8 | **0.500** | **0.756** | 0.0 | **[0, 2]** | 3/8 | [9,924, 9,934] | [5.03e-06, 5.04e-06] | 2.728e-05 | 4 |
| 2 | **A** | `n400_400` 400 v 400 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [9,909, 9,936] | [5.03e-06, 5.05e-06] | 1.245e-04 | 0 |
| 3 | **A** | `n400_ref_completa` 400 v 9,200 | 8 | **0.125** | **0.354** | 0.0 | **[0, 1]** | 1/8 | [9,924, 9,934] | [5.03e-06, 5.04e-06] | 7.133e-05 | 1 |
| 4 | **B** | `completo` 9,200 v 9,200 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [9,627, 9,636] | [5.19e-06, 5.19e-06] | 3.997e-05 | 0 |
| 5 | **B** | `n400_400` 400 v 400 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [9,590, 9,642] | [5.19e-06, 5.21e-06] | 8.673e-05 | 0 |
| 6 | **B** | `n400_ref_completa` 400 v 9,200 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [9,627, 9,636] | [5.19e-06, 5.19e-06] | 6.418e-05 | 0 |
| 7 | **C** | `completo` 9,200 v 9,200 | 8 | **0.500** | **0.756** | 0.0 | **[0, 2]** | 3/8 | [10,117, 10,131] | [4.94e-06, 4.94e-06] | 3.208e-05 | 2 |
| 8 | **C** | `n400_400` 400 v 400 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [10,066, 10,127] | [4.94e-06, 4.97e-06] | 5.909e-05 | 0 |
| 9 | **C** | `n400_ref_completa` 400 v 9,200 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [10,117, 10,131] | [4.94e-06, 4.94e-06] | 1.145e-04 | 0 |
| 10 | **D** | `completo` 9,200 v 9,200 | 8 | **1.125** | **2.475** | 0.0 | **[0, 7]** | 2/8 | [9,625, 9,637] | [5.19e-06, 5.19e-06] | 4.740e-05 | 9 |
| 11 | **D** | `n400_400` 400 v 400 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [9,599, 9,657] | [5.18e-06, 5.21e-06] | 8.534e-05 | 0 |
| 12 | **D** | `n400_ref_completa` 400 v 9,200 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [9,625, 9,637] | [5.19e-06, 5.19e-06] | 7.662e-05 | 0 |
| 13 | **E** | `completo` 9,200 v 9,200 | 8 | **0.875** | **1.126** | 0.5 | **[0, 3]** | 4/8 | [9,737, 9,753] | [5.13e-06, 5.14e-06] | 7.250e-06 | 6 |
| 14 | **E** | `n400_400` 400 v 400 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [9,731, 9,764] | [5.12e-06, 5.14e-06] | 1.177e-04 | 0 |
| 15 | **E** | `n400_ref_completa` 400 v 9,200 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [9,737, 9,753] | [5.13e-06, 5.14e-06] | 4.721e-05 | 0 |
| 16 | **F** | `completo` 9,200 v 9,200 | 8 | **1.250** | **1.035** | 1.0 | **[0, 3]** | 6/8 | [9,530, 9,542] | [5.24e-06, 5.25e-06] | 3.214e-06 | 7 |
| 17 | **F** | `n400_400` 400 v 400 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [9,504, 9,547] | [5.24e-06, 5.26e-06] | 4.050e-05 | 0 |
| 18 | **F** | `n400_ref_completa` 400 v 9,200 | 8 | **0.000** | **0.000** | 0.0 | **[0, 0]** | 0/8 | [9,530, 9,542] | [5.24e-06, 5.25e-06] | 9.524e-05 | 0 |

**El rango por contexto, en una sola linea por brazo** (sobre los 6 contextos
x 8 particiones de cada brazo):

| brazo | corridas | `n_sig` media | `n_sig` sd | rango global [min, max] | corridas con `n_sig>0` | suma de `n_sig` | genes distintos llamados |
|---|---|---|---|---|---|---|---|
| `completo` 9,200 v 9,200 | 48 | 0.7083 | 1.2709 | **[0, 7]** | 18/48 | 34 | 22 |
| `n400_400` 400 v 400 | 48 | 0.0000 | 0.0000 | **[0, 0]** | 0/48 | 0 | 0 |
| `n400_ref_completa` 400 v 9,200 | 48 | 0.0208 | 0.1443 | **[0, 1]** | 1/48 | 1 | 1 |

**Los genes que el nulo llego a llamar** (con en cuantas corridas, de las de
ese brazo):

| brazo | gen | corridas que lo llaman | ctx y j |
|---|---|---|---|
| `completo` 9,200 v 9,200 | `CST3` | 9 | A7, C4, C6, C7, E3, E6, F2, F3, F4 |
| `completo` 9,200 v 9,200 | `NUP37` | 4 | A0, E1, F0, F5 |
| `completo` 9,200 v 9,200 | `TPBG` | 2 | E4, F0 |
| `completo` 9,200 v 9,200 | `IRX5` | 1 | A0 |
| `completo` 9,200 v 9,200 | `SLC12A7` | 1 | A3 |
| `completo` 9,200 v 9,200 | `MYO18A` | 1 | C6 |
| `completo` 9,200 v 9,200 | `ITPR1` | 1 | D2 |
| `completo` 9,200 v 9,200 | `RET` | 1 | D2 |
| `completo` 9,200 v 9,200 | `C14orf93` | 1 | D6 |
| `completo` 9,200 v 9,200 | `CEP72` | 1 | D6 |
| `completo` 9,200 v 9,200 | `MCCC2` | 1 | D6 |
| `completo` 9,200 v 9,200 | `GRHPR` | 1 | D6 |
| `completo` 9,200 v 9,200 | `RASA1` | 1 | D6 |
| `completo` 9,200 v 9,200 | `TAF7` | 1 | D6 |
| `completo` 9,200 v 9,200 | `PSMB9` | 1 | D6 |
| `completo` 9,200 v 9,200 | `RRM2` | 1 | E1 |
| `completo` 9,200 v 9,200 | `MAPKAPK3` | 1 | E3 |
| `completo` 9,200 v 9,200 | `IRS1` | 1 | E3 |
| `completo` 9,200 v 9,200 | `BOD1L1` | 1 | F1 |
| `completo` 9,200 v 9,200 | `ABI3BP` | 1 | F1 |
| `completo` 9,200 v 9,200 | `SELENOW` | 1 | F3 |
| `completo` 9,200 v 9,200 | `PBX2` | 1 | F3 |
| `n400_ref_completa` 400 v 9,200 | `TCF19` | 1 | A3 |

🔴 **El nulo NO es intercambiable entre genes: hay recurrencia.** Si los
genes llamados fueran falsos positivos independientes, con tan pocas
llamadas en total la probabilidad de que un mismo gen salga varias veces
seria minuscula. Los que repiten:

| brazo | gen | corridas que lo llaman | de | ctx y j |
|---|---|---|---|---|
| `completo` 9,200 v 9,200 | **`CST3`** | **9** | 48 | A7, C4, C6, C7, E3, E6, F2, F3, F4 |
| `completo` 9,200 v 9,200 | **`NUP37`** | **4** | 48 | A0, E1, F0, F5 |
| `completo` 9,200 v 9,200 | **`TPBG`** | **2** | 48 | E4, F0 |

Esa recurrencia es la huella de la estructura de guia que la seccion 2.3
mide por separado: un gen cuyo nivel depende de QUE guia lleva la celula
sale llamado cada vez que la particion separa las guias que lo mueven.
Cuantifica tambien el limite de la propia construccion: las 8 particiones
de un contexto NO son 8 muestras independientes del nulo, porque comparten
las mismas 46 guias.

### 1.4 🔴 Los intermedios crudos: las corridas una por una

Ninguna cifra de las tablas de resumen es un agregado sin sus filas. Aqui
estan las **144 corridas** completas.

Columnas: `m` = genes probados tras el filtro CPM del preset; `p_crit` =
`0.05/m`, el escalon de BH de rango 1, que es lo que una senal AISLADA tiene
que batir; `frac p<0.05` = fraccion de los `m` genes con `p` cruda bajo 0.05,
que bajo un nulo uniforme vale ~0.05; `KS` = estadistico de
Kolmogorov-Smirnov de las `p` crudas contra la uniforme (diagnostico del nulo,
**no** parte de la puerta); `|lfc|` en el espacio log2 del propio marcador.

**TABLA 1.3 — las 144 corridas crudas.**

| # | ctx | j | semilla | brazo | n ref | n blanco | `m` | **`n_sig`** | `p_adj` min | `p` min | `p_crit` | `p` esp. min unif. | frac `p<0.05` | frac `p<0.01` | `p` mediana | KS unif. | `|lfc|` mediana | `|lfc|` p99 | `|lfc|` max | `lfc` NaN | seg |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **A** | 0 | `20274000` | `completo` | 9,200 | 9,200 | 9,926 | **2** | 0.0327 | 3.296e-06 | 5.037e-06 | 1.007e-04 | 0.0666 | 0.0160 | 0.4824 | 0.0286 | 0.0165 | 0.1192 | 0.3497 | 0 | 19.1 |
| 2 | **A** | 0 | `20274000` | `n400_400` | 400 | 400 | 9,930 | **0** | 0.6867 | 3.620e-04 | 5.035e-06 | 1.007e-04 | 0.0488 | 0.0111 | 0.4968 | 0.0079 | 0.0771 | 0.5528 | 1.4433 | 0 | 0.5 |
| 3 | **A** | 0 | `20274000` | `n400_ref` | 9,200 | 400 | 9,926 | **0** | 0.8303 | 1.327e-04 | 5.037e-06 | 1.007e-04 | 0.0526 | 0.0110 | 0.5061 | 0.0084 | 0.0546 | 0.3882 | 1.0029 | 0 | 6.1 |
| 4 | **A** | 1 | `20274001` | `completo` | 9,200 | 9,200 | 9,928 | **0** | 0.1798 | 2.668e-05 | 5.036e-06 | 1.007e-04 | 0.0518 | 0.0113 | 0.5011 | 0.0063 | 0.0160 | 0.1143 | 0.2437 | 0 | 12.0 |
| 5 | **A** | 1 | `20274001` | `n400_400` | 400 | 400 | 9,919 | **0** | 0.7208 | 1.304e-04 | 5.041e-06 | 1.008e-04 | 0.0442 | 0.0084 | 0.5145 | 0.0171 | 0.0754 | 0.5569 | 1.8148 | 0 | 0.5 |
| 6 | **A** | 1 | `20274001` | `n400_ref` | 9,200 | 400 | 9,928 | **0** | 0.9983 | 2.004e-04 | 5.036e-06 | 1.007e-04 | 0.0433 | 0.0075 | 0.5095 | 0.0141 | 0.0543 | 0.3915 | 0.9705 | 0 | 6.1 |
| 7 | **A** | 2 | `20274002` | `completo` | 9,200 | 9,200 | 9,927 | **0** | 0.5013 | 5.050e-05 | 5.037e-06 | 1.007e-04 | 0.0571 | 0.0143 | 0.4847 | 0.0174 | 0.0164 | 0.1202 | 0.2780 | 0 | 12.0 |
| 8 | **A** | 2 | `20274002` | `n400_400` | 400 | 400 | 9,909 | **0** | 0.6217 | 1.186e-04 | 5.046e-06 | 1.009e-04 | 0.0489 | 0.0105 | 0.5092 | 0.0102 | 0.0756 | 0.5450 | 1.2061 | 0 | 0.5 |
| 9 | **A** | 2 | `20274002` | `n400_ref` | 9,200 | 400 | 9,927 | **0** | 0.4827 | 4.862e-05 | 5.037e-06 | 1.007e-04 | 0.0531 | 0.0109 | 0.4945 | 0.0080 | 0.0565 | 0.3970 | 0.9061 | 0 | 6.1 |
| 10 | **A** | 3 | `20274003` | `completo` | 9,200 | 9,200 | 9,934 | **1** | 0.0157 | 1.580e-06 | 5.033e-06 | 1.007e-04 | 0.0626 | 0.0133 | 0.4800 | 0.0254 | 0.0172 | 0.1190 | 0.2658 | 0 | 12.6 |
| 11 | **A** | 3 | `20274003` | `n400_400` | 400 | 400 | 9,923 | **0** | 0.9224 | 1.371e-04 | 5.039e-06 | 1.008e-04 | 0.0529 | 0.0094 | 0.4905 | 0.0117 | 0.0807 | 0.5921 | 1.4561 | 0 | 0.5 |
| 12 | **A** | 3 | `20274003` | `n400_ref` | 9,200 | 400 | 9,934 | **1** | 0.0060 | 6.087e-07 | 5.033e-06 | 1.007e-04 | 0.0676 | 0.0165 | 0.4716 | 0.0383 | 0.0594 | 0.4198 | 1.0515 | 0 | 6.4 |
| 13 | **A** | 4 | `20274004` | `completo` | 9,200 | 9,200 | 9,928 | **0** | 0.4645 | 1.110e-04 | 5.036e-06 | 1.007e-04 | 0.0550 | 0.0127 | 0.4958 | 0.0154 | 0.0163 | 0.1133 | 0.2644 | 0 | 12.6 |
| 14 | **A** | 4 | `20274004` | `n400_400` | 400 | 400 | 9,926 | **0** | 0.4782 | 5.418e-05 | 5.037e-06 | 1.007e-04 | 0.0455 | 0.0083 | 0.5084 | 0.0127 | 0.0754 | 0.5746 | 1.5170 | 0 | 0.5 |
| 15 | **A** | 4 | `20274004` | `n400_ref` | 9,200 | 400 | 9,928 | **0** | 0.8820 | 2.243e-04 | 5.036e-06 | 1.007e-04 | 0.0424 | 0.0087 | 0.5076 | 0.0110 | 0.0559 | 0.4031 | 1.0584 | 0 | 6.4 |
| 16 | **A** | 5 | `20274005` | `completo` | 9,200 | 9,200 | 9,924 | **0** | 0.1876 | 3.490e-05 | 5.038e-06 | 1.008e-04 | 0.0665 | 0.0167 | 0.4716 | 0.0332 | 0.0166 | 0.1190 | 0.2485 | 0 | 12.5 |
| 17 | **A** | 5 | `20274005` | `n400_400` | 400 | 400 | 9,917 | **0** | 0.3165 | 3.622e-05 | 5.042e-06 | 1.008e-04 | 0.0536 | 0.0121 | 0.4888 | 0.0160 | 0.0806 | 0.5394 | 1.3100 | 0 | 0.5 |
| 18 | **A** | 5 | `20274005` | `n400_ref` | 9,200 | 400 | 9,924 | **0** | 0.5125 | 6.443e-05 | 5.038e-06 | 1.008e-04 | 0.0538 | 0.0101 | 0.4929 | 0.0085 | 0.0566 | 0.3828 | 0.7422 | 0 | 6.2 |
| 19 | **A** | 6 | `20274006` | `completo` | 9,200 | 9,200 | 9,924 | **0** | 0.2173 | 2.788e-05 | 5.038e-06 | 1.008e-04 | 0.0557 | 0.0121 | 0.4882 | 0.0139 | 0.0166 | 0.1142 | 0.2632 | 0 | 12.1 |
| 20 | **A** | 6 | `20274006` | `n400_400` | 400 | 400 | 9,922 | **0** | 0.1739 | 1.753e-05 | 5.039e-06 | 1.008e-04 | 0.0509 | 0.0099 | 0.5017 | 0.0076 | 0.0759 | 0.5648 | 2.8370 | 0 | 0.5 |
| 21 | **A** | 6 | `20274006` | `n400_ref` | 9,200 | 400 | 9,924 | **0** | 0.4182 | 7.824e-05 | 5.038e-06 | 1.008e-04 | 0.0527 | 0.0114 | 0.4901 | 0.0143 | 0.0559 | 0.4004 | 0.7516 | 0 | 6.2 |
| 22 | **A** | 7 | `20274007` | `completo` | 9,200 | 9,200 | 9,929 | **1** | 0.0044 | 4.462e-07 | 5.036e-06 | 1.007e-04 | 0.0538 | 0.0110 | 0.4969 | 0.0124 | 0.0161 | 0.1117 | 0.2332 | 0 | 12.3 |
| 23 | **A** | 7 | `20274007` | `n400_400` | 400 | 400 | 9,936 | **0** | 0.6622 | 1.547e-04 | 5.032e-06 | 1.006e-04 | 0.0500 | 0.0114 | 0.4894 | 0.0113 | 0.0771 | 0.5640 | 1.8571 | 0 | 0.5 |
| 24 | **A** | 7 | `20274007` | `n400_ref` | 9,200 | 400 | 9,929 | **0** | 0.1797 | 1.810e-05 | 5.036e-06 | 1.007e-04 | 0.0559 | 0.0122 | 0.4885 | 0.0160 | 0.0564 | 0.4081 | 0.7824 | 0 | 6.2 |
| 25 | **B** | 0 | `20274100` | `completo` | 9,200 | 9,200 | 9,632 | **0** | 0.3141 | 3.261e-05 | 5.191e-06 | 1.038e-04 | 0.0517 | 0.0123 | 0.5053 | 0.0093 | 0.0187 | 0.1536 | 0.6960 | 0 | 11.1 |
| 26 | **B** | 0 | `20274100` | `n400_400` | 400 | 400 | 9,593 | **0** | 0.5961 | 1.049e-04 | 5.212e-06 | 1.042e-04 | 0.0472 | 0.0107 | 0.5065 | 0.0069 | 0.0915 | 0.7665 | 5.9884 | 0 | 0.4 |
| 27 | **B** | 0 | `20274100` | `n400_ref` | 9,200 | 400 | 9,632 | **0** | 0.0920 | 9.552e-06 | 5.191e-06 | 1.038e-04 | 0.0580 | 0.0124 | 0.4821 | 0.0239 | 0.0693 | 0.5372 | 2.0098 | 0 | 5.7 |
| 28 | **B** | 1 | `20274101` | `completo` | 9,200 | 9,200 | 9,636 | **0** | 0.4443 | 9.189e-05 | 5.189e-06 | 1.038e-04 | 0.0581 | 0.0115 | 0.4858 | 0.0171 | 0.0188 | 0.1493 | 0.4648 | 0 | 11.1 |
| 29 | **B** | 1 | `20274101` | `n400_400` | 400 | 400 | 9,642 | **0** | 0.0962 | 1.435e-05 | 5.186e-06 | 1.037e-04 | 0.0617 | 0.0155 | 0.4902 | 0.0201 | 0.0877 | 0.8095 | 32.8837 | 0 | 0.5 |
| 30 | **B** | 1 | `20274101` | `n400_ref` | 9,200 | 400 | 9,636 | **0** | 0.1207 | 1.253e-05 | 5.189e-06 | 1.038e-04 | 0.0569 | 0.0133 | 0.5019 | 0.0093 | 0.0597 | 0.4830 | 4.1428 | 0 | 5.6 |
| 31 | **B** | 2 | `20274102` | `completo` | 9,200 | 9,200 | 9,636 | **0** | 0.0794 | 8.237e-06 | 5.189e-06 | 1.038e-04 | 0.0439 | 0.0074 | 0.5063 | 0.0091 | 0.0186 | 0.1433 | 0.7470 | 0 | 11.2 |
| 32 | **B** | 2 | `20274102` | `n400_400` | 400 | 400 | 9,619 | **0** | 0.4369 | 4.542e-05 | 5.198e-06 | 1.040e-04 | 0.0495 | 0.0100 | 0.4922 | 0.0124 | 0.0933 | 0.7818 | 32.4881 | 0 | 0.4 |
| 33 | **B** | 2 | `20274102` | `n400_ref` | 9,200 | 400 | 9,636 | **0** | 0.9969 | 2.868e-04 | 5.189e-06 | 1.038e-04 | 0.0465 | 0.0089 | 0.5061 | 0.0119 | 0.0636 | 0.4961 | 1.4790 | 0 | 5.6 |
| 34 | **B** | 3 | `20274103` | `completo` | 9,200 | 9,200 | 9,632 | **0** | 0.0703 | 7.302e-06 | 5.191e-06 | 1.038e-04 | 0.0545 | 0.0115 | 0.4852 | 0.0182 | 0.0192 | 0.1442 | 0.4423 | 0 | 11.2 |
| 35 | **B** | 3 | `20274103` | `n400_400` | 400 | 400 | 9,590 | **0** | 0.9406 | 1.535e-04 | 5.214e-06 | 1.043e-04 | 0.0429 | 0.0089 | 0.5171 | 0.0192 | 0.0844 | 0.6469 | 3.8063 | 0 | 0.4 |
| 36 | **B** | 3 | `20274103` | `n400_ref` | 9,200 | 400 | 9,632 | **0** | 0.9674 | 1.246e-04 | 5.191e-06 | 1.038e-04 | 0.0422 | 0.0097 | 0.5116 | 0.0193 | 0.0594 | 0.4951 | 2.5156 | 0 | 5.6 |
| 37 | **B** | 4 | `20274104` | `completo` | 9,200 | 9,200 | 9,627 | **0** | 0.7662 | 2.677e-04 | 5.194e-06 | 1.039e-04 | 0.0521 | 0.0098 | 0.4943 | 0.0078 | 0.0189 | 0.1469 | 1.0831 | 0 | 11.2 |
| 38 | **B** | 4 | `20274104` | `n400_400` | 400 | 400 | 9,630 | **0** | 0.6603 | 6.857e-05 | 5.192e-06 | 1.038e-04 | 0.0505 | 0.0088 | 0.4990 | 0.0058 | 0.0847 | 0.6999 | 32.4846 | 0 | 0.5 |
| 39 | **B** | 4 | `20274104` | `n400_ref` | 9,200 | 400 | 9,627 | **0** | 0.4825 | 5.012e-05 | 5.194e-06 | 1.039e-04 | 0.0495 | 0.0105 | 0.5003 | 0.0064 | 0.0624 | 0.4923 | 4.0310 | 0 | 5.7 |
| 40 | **B** | 5 | `20274105` | `completo` | 9,200 | 9,200 | 9,633 | **0** | 0.1641 | 1.704e-05 | 5.190e-06 | 1.038e-04 | 0.0457 | 0.0093 | 0.5054 | 0.0151 | 0.0185 | 0.1473 | 0.4947 | 0 | 11.1 |
| 41 | **B** | 5 | `20274105` | `n400_400` | 400 | 400 | 9,623 | **0** | 0.9942 | 2.149e-04 | 5.196e-06 | 1.039e-04 | 0.0440 | 0.0080 | 0.5082 | 0.0145 | 0.0868 | 0.7500 | 33.1910 | 0 | 0.4 |
| 42 | **B** | 5 | `20274105` | `n400_ref` | 9,200 | 400 | 9,633 | **0** | 0.7638 | 2.385e-04 | 5.190e-06 | 1.038e-04 | 0.0435 | 0.0094 | 0.5126 | 0.0183 | 0.0631 | 0.4868 | 32.4531 | 0 | 5.6 |
| 43 | **B** | 6 | `20274106` | `completo` | 9,200 | 9,200 | 9,631 | **0** | 0.4149 | 4.733e-05 | 5.192e-06 | 1.038e-04 | 0.0515 | 0.0107 | 0.4985 | 0.0040 | 0.0185 | 0.1443 | 0.6715 | 0 | 11.2 |
| 44 | **B** | 6 | `20274106` | `n400_400` | 400 | 400 | 9,637 | **0** | 0.7098 | 4.001e-04 | 5.188e-06 | 1.038e-04 | 0.0563 | 0.0125 | 0.4923 | 0.0135 | 0.0872 | 0.7063 | 32.4402 | 0 | 0.5 |
| 45 | **B** | 6 | `20274106` | `n400_ref` | 9,200 | 400 | 9,631 | **0** | 0.2725 | 4.260e-05 | 5.192e-06 | 1.038e-04 | 0.0536 | 0.0110 | 0.4974 | 0.0072 | 0.0612 | 0.4581 | 1.8749 | 0 | 5.7 |
| 46 | **B** | 7 | `20274107` | `completo` | 9,200 | 9,200 | 9,635 | **0** | 0.5498 | 5.706e-05 | 5.189e-06 | 1.038e-04 | 0.0541 | 0.0097 | 0.4938 | 0.0157 | 0.0194 | 0.1496 | 0.5676 | 0 | 11.2 |
| 47 | **B** | 7 | `20274107` | `n400_400` | 400 | 400 | 9,624 | **0** | 0.2523 | 2.621e-05 | 5.195e-06 | 1.039e-04 | 0.0715 | 0.0177 | 0.4516 | 0.0509 | 0.0892 | 0.7387 | 3.9362 | 0 | 0.5 |
| 48 | **B** | 7 | `20274107` | `n400_ref` | 9,200 | 400 | 9,635 | **0** | 0.4608 | 7.823e-05 | 5.189e-06 | 1.038e-04 | 0.0646 | 0.0156 | 0.4734 | 0.0301 | 0.0680 | 0.5541 | 1.7714 | 0 | 5.6 |
| 49 | **C** | 0 | `20274200` | `completo` | 9,200 | 9,200 | 10,123 | **0** | 0.6831 | 1.800e-04 | 4.939e-06 | 9.878e-05 | 0.0563 | 0.0122 | 0.4806 | 0.0230 | 0.0169 | 0.1176 | 0.5880 | 0 | 11.8 |
| 50 | **C** | 0 | `20274200` | `n400_400` | 400 | 400 | 10,120 | **0** | 0.3843 | 3.798e-05 | 4.941e-06 | 9.880e-05 | 0.0510 | 0.0121 | 0.4929 | 0.0119 | 0.0829 | 0.5823 | 3.2592 | 0 | 0.5 |
| 51 | **C** | 0 | `20274200` | `n400_ref` | 9,200 | 400 | 10,123 | **0** | 0.8447 | 8.344e-05 | 4.939e-06 | 9.878e-05 | 0.0460 | 0.0086 | 0.4997 | 0.0074 | 0.0574 | 0.4095 | 1.3169 | 0 | 6.2 |
| 52 | **C** | 1 | `20274201` | `completo` | 9,200 | 9,200 | 10,124 | **0** | 0.2701 | 6.336e-05 | 4.939e-06 | 9.877e-05 | 0.0498 | 0.0117 | 0.5026 | 0.0090 | 0.0166 | 0.1153 | 0.7871 | 0 | 12.0 |
| 53 | **C** | 1 | `20274201` | `n400_400` | 400 | 400 | 10,096 | **0** | 0.7955 | 7.880e-05 | 4.952e-06 | 9.904e-05 | 0.0467 | 0.0085 | 0.5085 | 0.0139 | 0.0761 | 0.5681 | 1.7734 | 0 | 0.5 |
| 54 | **C** | 1 | `20274201` | `n400_ref` | 9,200 | 400 | 10,124 | **0** | 0.8100 | 1.266e-04 | 4.939e-06 | 9.877e-05 | 0.0478 | 0.0103 | 0.5035 | 0.0104 | 0.0570 | 0.3907 | 2.6352 | 0 | 6.0 |
| 55 | **C** | 2 | `20274202` | `completo` | 9,200 | 9,200 | 10,117 | **0** | 0.4996 | 1.320e-04 | 4.942e-06 | 9.883e-05 | 0.0484 | 0.0097 | 0.5064 | 0.0110 | 0.0166 | 0.1163 | 0.2405 | 0 | 11.8 |
| 56 | **C** | 2 | `20274202` | `n400_400` | 400 | 400 | 10,066 | **0** | 0.1730 | 1.718e-05 | 4.967e-06 | 9.933e-05 | 0.0472 | 0.0099 | 0.5011 | 0.0088 | 0.0806 | 0.5529 | 2.7004 | 0 | 0.5 |
| 57 | **C** | 2 | `20274202` | `n400_ref` | 9,200 | 400 | 10,117 | **0** | 0.5207 | 2.263e-04 | 4.942e-06 | 9.883e-05 | 0.0497 | 0.0105 | 0.5151 | 0.0166 | 0.0563 | 0.4000 | 2.0024 | 0 | 6.0 |
| 58 | **C** | 3 | `20274203` | `completo` | 9,200 | 9,200 | 10,124 | **0** | 0.4479 | 5.504e-05 | 4.939e-06 | 9.877e-05 | 0.0520 | 0.0094 | 0.4967 | 0.0062 | 0.0169 | 0.1175 | 0.8278 | 0 | 11.8 |
| 59 | **C** | 3 | `20274203` | `n400_400` | 400 | 400 | 10,102 | **0** | 0.3633 | 4.603e-05 | 4.950e-06 | 9.898e-05 | 0.0475 | 0.0101 | 0.5009 | 0.0062 | 0.0774 | 0.5558 | 2.3272 | 0 | 0.5 |
| 60 | **C** | 3 | `20274203` | `n400_ref` | 9,200 | 400 | 10,124 | **0** | 0.9999 | 1.902e-04 | 4.939e-06 | 9.877e-05 | 0.0491 | 0.0091 | 0.5057 | 0.0117 | 0.0556 | 0.4083 | 1.3058 | 0 | 5.9 |
| 61 | **C** | 4 | `20274204` | `completo` | 9,200 | 9,200 | 10,121 | **1** | 0.0004 | 4.363e-08 | 4.940e-06 | 9.879e-05 | 0.0497 | 0.0107 | 0.5035 | 0.0088 | 0.0165 | 0.1179 | 0.6944 | 0 | 11.8 |
| 62 | **C** | 4 | `20274204` | `n400_400` | 400 | 400 | 10,111 | **0** | 0.9821 | 9.713e-05 | 4.945e-06 | 9.889e-05 | 0.0415 | 0.0074 | 0.5068 | 0.0134 | 0.0783 | 0.5552 | 1.1992 | 0 | 0.5 |
| 63 | **C** | 4 | `20274204` | `n400_ref` | 9,200 | 400 | 10,121 | **0** | 0.9974 | 4.472e-04 | 4.940e-06 | 9.879e-05 | 0.0443 | 0.0075 | 0.5012 | 0.0088 | 0.0566 | 0.3950 | 0.9610 | 0 | 6.0 |
| 64 | **C** | 5 | `20274205` | `completo` | 9,200 | 9,200 | 10,130 | **0** | 0.0923 | 9.115e-06 | 4.936e-06 | 9.871e-05 | 0.0513 | 0.0114 | 0.4867 | 0.0156 | 0.0169 | 0.1176 | 0.4405 | 0 | 11.9 |
| 65 | **C** | 5 | `20274205` | `n400_400` | 400 | 400 | 10,127 | **0** | 0.6358 | 7.216e-05 | 4.937e-06 | 9.874e-05 | 0.0536 | 0.0119 | 0.4887 | 0.0148 | 0.0788 | 0.5435 | 2.6201 | 0 | 0.5 |
| 66 | **C** | 5 | `20274205` | `n400_ref` | 9,200 | 400 | 10,130 | **0** | 0.1601 | 2.741e-05 | 4.936e-06 | 9.871e-05 | 0.0492 | 0.0102 | 0.5052 | 0.0065 | 0.0554 | 0.4068 | 1.1291 | 0 | 6.0 |
| 67 | **C** | 6 | `20274206` | `completo` | 9,200 | 9,200 | 10,130 | **2** | 0.0008 | 7.907e-08 | 4.936e-06 | 9.871e-05 | 0.0530 | 0.0113 | 0.4986 | 0.0086 | 0.0168 | 0.1221 | 0.8575 | 0 | 11.8 |
| 68 | **C** | 6 | `20274206` | `n400_400` | 400 | 400 | 10,120 | **0** | 0.7130 | 3.584e-04 | 4.941e-06 | 9.880e-05 | 0.0490 | 0.0098 | 0.5020 | 0.0049 | 0.0802 | 0.5742 | 2.4903 | 0 | 0.5 |
| 69 | **C** | 6 | `20274206` | `n400_ref` | 9,200 | 400 | 10,130 | **0** | 0.1859 | 1.836e-05 | 4.936e-06 | 9.871e-05 | 0.0433 | 0.0085 | 0.5111 | 0.0172 | 0.0568 | 0.4001 | 1.9113 | 0 | 6.0 |
| 70 | **C** | 7 | `20274207` | `completo` | 9,200 | 9,200 | 10,131 | **1** | 0.0046 | 4.537e-07 | 4.935e-06 | 9.870e-05 | 0.0530 | 0.0113 | 0.4936 | 0.0120 | 0.0167 | 0.1207 | 0.9799 | 0 | 11.9 |
| 71 | **C** | 7 | `20274207` | `n400_400` | 400 | 400 | 10,108 | **0** | 0.0570 | 5.641e-06 | 4.947e-06 | 9.892e-05 | 0.0466 | 0.0074 | 0.5058 | 0.0088 | 0.0787 | 0.5683 | 2.6656 | 0 | 0.5 |
| 72 | **C** | 7 | `20274207` | `n400_ref` | 9,200 | 400 | 10,131 | **0** | 0.8714 | 1.025e-04 | 4.935e-06 | 9.870e-05 | 0.0494 | 0.0101 | 0.5064 | 0.0115 | 0.0558 | 0.4202 | 1.9883 | 0 | 6.1 |
| 73 | **D** | 0 | `20274300` | `completo` | 9,200 | 9,200 | 9,636 | **0** | 0.7381 | 1.672e-04 | 5.189e-06 | 1.038e-04 | 0.0568 | 0.0123 | 0.4730 | 0.0284 | 0.0178 | 0.1336 | 0.4984 | 0 | 11.5 |
| 74 | **D** | 0 | `20274300` | `n400_400` | 400 | 400 | 9,619 | **0** | 0.3758 | 1.058e-04 | 5.198e-06 | 1.040e-04 | 0.0609 | 0.0125 | 0.4749 | 0.0289 | 0.0840 | 0.5927 | 4.8852 | 0 | 0.5 |
| 75 | **D** | 0 | `20274300` | `n400_ref` | 9,200 | 400 | 9,636 | **0** | 0.9984 | 2.145e-04 | 5.189e-06 | 1.038e-04 | 0.0440 | 0.0093 | 0.5043 | 0.0116 | 0.0578 | 0.4261 | 2.0273 | 0 | 5.8 |
| 76 | **D** | 1 | `20274301` | `completo` | 9,200 | 9,200 | 9,630 | **0** | 0.2675 | 6.352e-05 | 5.192e-06 | 1.038e-04 | 0.0526 | 0.0105 | 0.4945 | 0.0084 | 0.0175 | 0.1273 | 0.9641 | 0 | 11.4 |
| 77 | **D** | 1 | `20274301` | `n400_400` | 400 | 400 | 9,618 | **0** | 0.0548 | 5.700e-06 | 5.199e-06 | 1.040e-04 | 0.0459 | 0.0073 | 0.4939 | 0.0119 | 0.0812 | 0.6168 | 6.7964 | 0 | 0.5 |
| 78 | **D** | 1 | `20274301` | `n400_ref` | 9,200 | 400 | 9,630 | **0** | 0.8316 | 8.636e-05 | 5.192e-06 | 1.038e-04 | 0.0442 | 0.0084 | 0.5134 | 0.0176 | 0.0579 | 0.4271 | 4.3087 | 0 | 5.8 |
| 79 | **D** | 2 | `20274302` | `completo` | 9,200 | 9,200 | 9,633 | **2** | 0.0414 | 7.012e-06 | 5.190e-06 | 1.038e-04 | 0.0618 | 0.0122 | 0.4753 | 0.0283 | 0.0184 | 0.1253 | 0.3731 | 0 | 11.4 |
| 80 | **D** | 2 | `20274302` | `n400_400` | 400 | 400 | 9,642 | **0** | 0.8055 | 8.354e-05 | 5.186e-06 | 1.037e-04 | 0.0419 | 0.0079 | 0.5170 | 0.0231 | 0.0806 | 0.6196 | 32.6137 | 0 | 0.5 |
| 81 | **D** | 2 | `20274302` | `n400_ref` | 9,200 | 400 | 9,633 | **0** | 0.9990 | 1.700e-04 | 5.190e-06 | 1.038e-04 | 0.0407 | 0.0075 | 0.5240 | 0.0275 | 0.0570 | 0.4146 | 0.9879 | 0 | 5.8 |
| 82 | **D** | 3 | `20274303` | `completo` | 9,200 | 9,200 | 9,637 | **0** | 0.3190 | 4.792e-05 | 5.188e-06 | 1.038e-04 | 0.0566 | 0.0113 | 0.4863 | 0.0174 | 0.0176 | 0.1342 | 0.6823 | 0 | 11.3 |
| 83 | **D** | 3 | `20274303` | `n400_400` | 400 | 400 | 9,615 | **0** | 0.9995 | 3.601e-04 | 5.200e-06 | 1.040e-04 | 0.0433 | 0.0079 | 0.5211 | 0.0287 | 0.0822 | 0.6089 | 2.2319 | 0 | 0.5 |
| 84 | **D** | 3 | `20274303` | `n400_ref` | 9,200 | 400 | 9,637 | **0** | 0.5944 | 6.168e-05 | 5.188e-06 | 1.038e-04 | 0.0612 | 0.0119 | 0.4772 | 0.0293 | 0.0627 | 0.4311 | 0.9042 | 0 | 5.9 |
| 85 | **D** | 4 | `20274304` | `completo` | 9,200 | 9,200 | 9,635 | **0** | 0.2920 | 4.688e-05 | 5.189e-06 | 1.038e-04 | 0.0622 | 0.0147 | 0.4734 | 0.0295 | 0.0187 | 0.1390 | 0.4307 | 0 | 11.5 |
| 86 | **D** | 4 | `20274304` | `n400_400` | 400 | 400 | 9,657 | **0** | 0.9949 | 5.590e-04 | 5.178e-06 | 1.035e-04 | 0.0447 | 0.0075 | 0.5093 | 0.0186 | 0.0861 | 0.6685 | 3.3767 | 0 | 0.4 |
| 87 | **D** | 4 | `20274304` | `n400_ref` | 9,200 | 400 | 9,635 | **0** | 0.2486 | 2.580e-05 | 5.189e-06 | 1.038e-04 | 0.0449 | 0.0080 | 0.5103 | 0.0125 | 0.0619 | 0.4347 | 1.3805 | 0 | 5.8 |
| 88 | **D** | 5 | `20274305` | `completo` | 9,200 | 9,200 | 9,631 | **0** | 0.0696 | 7.231e-06 | 5.192e-06 | 1.038e-04 | 0.0609 | 0.0113 | 0.4792 | 0.0250 | 0.0182 | 0.1254 | 0.5021 | 0 | 11.4 |
| 89 | **D** | 5 | `20274305` | `n400_400` | 400 | 400 | 9,599 | **0** | 0.8364 | 8.713e-05 | 5.209e-06 | 1.042e-04 | 0.0507 | 0.0091 | 0.5100 | 0.0138 | 0.0797 | 0.5800 | 3.8508 | 0 | 0.5 |
| 90 | **D** | 5 | `20274305` | `n400_ref` | 9,200 | 400 | 9,631 | **0** | 0.9259 | 1.773e-04 | 5.192e-06 | 1.038e-04 | 0.0507 | 0.0087 | 0.5013 | 0.0089 | 0.0575 | 0.4252 | 1.0168 | 0 | 5.9 |
| 91 | **D** | 6 | `20274306` | `completo` | 9,200 | 9,200 | 9,625 | **7** | 0.0120 | 1.248e-06 | 5.195e-06 | 1.039e-04 | 0.0845 | 0.0224 | 0.4523 | 0.0583 | 0.0178 | 0.1365 | 0.4259 | 0 | 11.4 |
| 92 | **D** | 6 | `20274306` | `n400_400` | 400 | 400 | 9,609 | **0** | 0.2087 | 2.172e-05 | 5.203e-06 | 1.041e-04 | 0.0567 | 0.0129 | 0.4823 | 0.0220 | 0.0833 | 0.5879 | 2.2923 | 0 | 0.5 |
| 93 | **D** | 6 | `20274306` | `n400_ref` | 9,200 | 400 | 9,625 | **0** | 0.0764 | 7.935e-06 | 5.195e-06 | 1.039e-04 | 0.0464 | 0.0098 | 0.5011 | 0.0083 | 0.0591 | 0.4281 | 0.8945 | 0 | 5.8 |
| 94 | **D** | 7 | `20274307` | `completo` | 9,200 | 9,200 | 9,632 | **0** | 0.7706 | 8.001e-05 | 5.191e-06 | 1.038e-04 | 0.0422 | 0.0071 | 0.5049 | 0.0184 | 0.0172 | 0.1192 | 0.4142 | 0 | 11.4 |
| 95 | **D** | 7 | `20274307` | `n400_400` | 400 | 400 | 9,640 | **0** | 0.5748 | 5.962e-05 | 5.187e-06 | 1.037e-04 | 0.0457 | 0.0093 | 0.5052 | 0.0124 | 0.0810 | 0.6623 | 5.1090 | 0 | 0.5 |
| 96 | **D** | 7 | `20274307` | `n400_ref` | 9,200 | 400 | 9,632 | **0** | 0.4862 | 6.689e-05 | 5.191e-06 | 1.038e-04 | 0.0487 | 0.0096 | 0.4988 | 0.0049 | 0.0591 | 0.4126 | 0.9788 | 0 | 5.9 |
| 97 | **E** | 0 | `20274400` | `completo` | 9,200 | 9,200 | 9,752 | **0** | 0.1012 | 1.038e-05 | 5.127e-06 | 1.025e-04 | 0.0517 | 0.0109 | 0.4964 | 0.0073 | 0.0171 | 0.1253 | 1.1609 | 0 | 11.1 |
| 98 | **E** | 0 | `20274400` | `n400_400` | 400 | 400 | 9,755 | **0** | 0.2777 | 6.215e-05 | 5.126e-06 | 1.025e-04 | 0.0513 | 0.0103 | 0.4855 | 0.0156 | 0.0828 | 0.5830 | 3.8046 | 0 | 0.4 |
| 99 | **E** | 0 | `20274400` | `n400_ref` | 9,200 | 400 | 9,752 | **0** | 0.2757 | 4.150e-05 | 5.127e-06 | 1.025e-04 | 0.0489 | 0.0107 | 0.4963 | 0.0091 | 0.0596 | 0.4316 | 1.8369 | 0 | 5.6 |
| 100 | **E** | 1 | `20274401` | `completo` | 9,200 | 9,200 | 9,745 | **2** | 0.0240 | 4.119e-06 | 5.131e-06 | 1.026e-04 | 0.0758 | 0.0167 | 0.4595 | 0.0459 | 0.0182 | 0.1233 | 1.3157 | 0 | 11.1 |
| 101 | **E** | 1 | `20274401` | `n400_400` | 400 | 400 | 9,736 | **0** | 0.4142 | 1.481e-04 | 5.136e-06 | 1.027e-04 | 0.0490 | 0.0108 | 0.5074 | 0.0125 | 0.0817 | 0.5930 | 2.3648 | 0 | 0.4 |
| 102 | **E** | 1 | `20274401` | `n400_ref` | 9,200 | 400 | 9,745 | **0** | 0.5158 | 5.293e-05 | 5.131e-06 | 1.026e-04 | 0.0546 | 0.0109 | 0.4894 | 0.0150 | 0.0595 | 0.4340 | 1.9250 | 0 | 5.6 |
| 103 | **E** | 2 | `20274402` | `completo` | 9,200 | 9,200 | 9,749 | **0** | 0.3534 | 4.636e-05 | 5.129e-06 | 1.026e-04 | 0.0592 | 0.0140 | 0.4784 | 0.0272 | 0.0180 | 0.1221 | 0.7531 | 0 | 11.1 |
| 104 | **E** | 2 | `20274402` | `n400_400` | 400 | 400 | 9,745 | **0** | 0.9985 | 3.600e-04 | 5.131e-06 | 1.026e-04 | 0.0437 | 0.0087 | 0.5213 | 0.0240 | 0.0811 | 0.5744 | 2.3338 | 0 | 0.4 |
| 105 | **E** | 2 | `20274402` | `n400_ref` | 9,200 | 400 | 9,749 | **0** | 0.9994 | 1.374e-04 | 5.129e-06 | 1.026e-04 | 0.0432 | 0.0070 | 0.5220 | 0.0237 | 0.0580 | 0.4058 | 1.0828 | 0 | 5.7 |
| 106 | **E** | 3 | `20274403` | `completo` | 9,200 | 9,200 | 9,753 | **3** | 0.0080 | 8.182e-07 | 5.127e-06 | 1.025e-04 | 0.0659 | 0.0143 | 0.4867 | 0.0207 | 0.0179 | 0.1243 | 0.8744 | 0 | 11.1 |
| 107 | **E** | 3 | `20274403` | `n400_400` | 400 | 400 | 9,731 | **0** | 0.4456 | 8.738e-05 | 5.138e-06 | 1.028e-04 | 0.0552 | 0.0122 | 0.4844 | 0.0197 | 0.0834 | 0.5862 | 1.8732 | 0 | 0.4 |
| 108 | **E** | 3 | `20274403` | `n400_ref` | 9,200 | 400 | 9,753 | **0** | 0.0615 | 6.308e-06 | 5.127e-06 | 1.025e-04 | 0.0533 | 0.0117 | 0.4949 | 0.0092 | 0.0585 | 0.4253 | 1.7774 | 0 | 5.6 |
| 109 | **E** | 4 | `20274404` | `completo` | 9,200 | 9,200 | 9,743 | **1** | 0.0109 | 1.115e-06 | 5.132e-06 | 1.026e-04 | 0.0503 | 0.0089 | 0.5023 | 0.0066 | 0.0170 | 0.1253 | 1.1248 | 0 | 11.1 |
| 110 | **E** | 4 | `20274404` | `n400_400` | 400 | 400 | 9,764 | **0** | 0.0928 | 9.504e-06 | 5.121e-06 | 1.024e-04 | 0.0519 | 0.0118 | 0.4924 | 0.0116 | 0.0837 | 0.6101 | 4.2209 | 0 | 0.5 |
| 111 | **E** | 4 | `20274404` | `n400_ref` | 9,200 | 400 | 9,743 | **0** | 0.3571 | 3.665e-05 | 5.132e-06 | 1.026e-04 | 0.0631 | 0.0129 | 0.4785 | 0.0297 | 0.0626 | 0.4485 | 2.4526 | 0 | 5.6 |
| 112 | **E** | 5 | `20274405` | `completo` | 9,200 | 9,200 | 9,737 | **0** | 0.1535 | 1.577e-05 | 5.135e-06 | 1.027e-04 | 0.0536 | 0.0122 | 0.4968 | 0.0094 | 0.0171 | 0.1187 | 0.7342 | 0 | 11.1 |
| 113 | **E** | 5 | `20274405` | `n400_400` | 400 | 400 | 9,732 | **0** | 0.9990 | 2.654e-04 | 5.138e-06 | 1.027e-04 | 0.0465 | 0.0084 | 0.5062 | 0.0095 | 0.0795 | 0.5828 | 5.2059 | 0 | 0.4 |
| 114 | **E** | 5 | `20274405` | `n400_ref` | 9,200 | 400 | 9,737 | **0** | 0.0914 | 2.580e-05 | 5.135e-06 | 1.027e-04 | 0.0481 | 0.0096 | 0.5039 | 0.0087 | 0.0564 | 0.3972 | 2.8055 | 0 | 5.7 |
| 115 | **E** | 6 | `20274406` | `completo` | 9,200 | 9,200 | 9,747 | **1** | 0.0000 | 4.946e-09 | 5.130e-06 | 1.026e-04 | 0.0544 | 0.0137 | 0.4906 | 0.0162 | 0.0176 | 0.1259 | 1.0868 | 0 | 11.1 |
| 116 | **E** | 6 | `20274406` | `n400_400` | 400 | 400 | 9,760 | **0** | 0.6955 | 1.779e-04 | 5.123e-06 | 1.024e-04 | 0.0516 | 0.0106 | 0.5067 | 0.0092 | 0.0821 | 0.6340 | 32.5648 | 0 | 0.5 |
| 117 | **E** | 6 | `20274406` | `n400_ref` | 9,200 | 400 | 9,747 | **0** | 0.6246 | 2.324e-04 | 5.130e-06 | 1.026e-04 | 0.0489 | 0.0112 | 0.5027 | 0.0129 | 0.0574 | 0.4119 | 2.1209 | 0 | 5.7 |
| 118 | **E** | 7 | `20274407` | `completo` | 9,200 | 9,200 | 9,746 | **0** | 0.0870 | 1.540e-05 | 5.130e-06 | 1.026e-04 | 0.0559 | 0.0115 | 0.4953 | 0.0082 | 0.0176 | 0.1219 | 1.2827 | 0 | 11.1 |
| 119 | **E** | 7 | `20274407` | `n400_400` | 400 | 400 | 9,732 | **0** | 0.6512 | 6.691e-05 | 5.138e-06 | 1.027e-04 | 0.0467 | 0.0103 | 0.5076 | 0.0155 | 0.0789 | 0.5859 | 1.4517 | 0 | 0.4 |
| 120 | **E** | 7 | `20274407` | `n400_ref` | 9,200 | 400 | 9,746 | **0** | 0.7992 | 8.200e-05 | 5.130e-06 | 1.026e-04 | 0.0469 | 0.0078 | 0.5162 | 0.0189 | 0.0576 | 0.4324 | 1.5991 | 0 | 5.6 |
| 121 | **F** | 0 | `20274500` | `completo` | 9,200 | 9,200 | 9,542 | **2** | 0.0277 | 3.973e-06 | 5.240e-06 | 1.048e-04 | 0.0582 | 0.0127 | 0.4861 | 0.0152 | 0.0177 | 0.1317 | 1.0473 | 0 | 11.2 |
| 122 | **F** | 0 | `20274500` | `n400_400` | 400 | 400 | 9,541 | **0** | 0.2813 | 2.948e-05 | 5.241e-06 | 1.048e-04 | 0.0488 | 0.0086 | 0.5096 | 0.0144 | 0.0780 | 0.6392 | 3.6293 | 0 | 0.4 |
| 123 | **F** | 0 | `20274500` | `n400_ref` | 9,200 | 400 | 9,542 | **0** | 0.4923 | 1.687e-04 | 5.240e-06 | 1.048e-04 | 0.0506 | 0.0115 | 0.5103 | 0.0164 | 0.0579 | 0.4547 | 1.8246 | 0 | 5.6 |
| 124 | **F** | 1 | `20274501` | `completo` | 9,200 | 9,200 | 9,538 | **2** | 0.0454 | 6.172e-06 | 5.242e-06 | 1.048e-04 | 0.0606 | 0.0140 | 0.4869 | 0.0167 | 0.0176 | 0.1274 | 0.7177 | 0 | 11.1 |
| 125 | **F** | 1 | `20274501` | `n400_400` | 400 | 400 | 9,504 | **0** | 0.8942 | 3.949e-04 | 5.261e-06 | 1.052e-04 | 0.0477 | 0.0092 | 0.5045 | 0.0106 | 0.0823 | 0.6103 | 2.2451 | 0 | 0.4 |
| 126 | **F** | 1 | `20274501` | `n400_ref` | 9,200 | 400 | 9,538 | **0** | 0.7799 | 1.503e-04 | 5.242e-06 | 1.048e-04 | 0.0478 | 0.0094 | 0.5034 | 0.0100 | 0.0601 | 0.4485 | 1.4290 | 0 | 5.6 |
| 127 | **F** | 2 | `20274502` | `completo` | 9,200 | 9,200 | 9,538 | **1** | 0.0067 | 7.042e-07 | 5.242e-06 | 1.048e-04 | 0.0528 | 0.0116 | 0.4899 | 0.0131 | 0.0178 | 0.1231 | 0.6703 | 0 | 11.1 |
| 128 | **F** | 2 | `20274502` | `n400_400` | 400 | 400 | 9,529 | **0** | 0.2556 | 2.682e-05 | 5.247e-06 | 1.049e-04 | 0.0471 | 0.0070 | 0.5102 | 0.0145 | 0.0817 | 0.6183 | 3.2270 | 0 | 0.5 |
| 129 | **F** | 2 | `20274502` | `n400_ref` | 9,200 | 400 | 9,538 | **0** | 0.2242 | 2.351e-05 | 5.242e-06 | 1.048e-04 | 0.0410 | 0.0085 | 0.5068 | 0.0161 | 0.0582 | 0.4437 | 1.4874 | 0 | 5.7 |
| 130 | **F** | 3 | `20274503` | `completo` | 9,200 | 9,200 | 9,539 | **3** | 0.0076 | 7.933e-07 | 5.242e-06 | 1.048e-04 | 0.0569 | 0.0142 | 0.4864 | 0.0187 | 0.0176 | 0.1279 | 0.4486 | 0 | 11.6 |
| 131 | **F** | 3 | `20274503` | `n400_400` | 400 | 400 | 9,509 | **0** | 0.1185 | 1.246e-05 | 5.258e-06 | 1.052e-04 | 0.0490 | 0.0085 | 0.5096 | 0.0118 | 0.0825 | 0.6073 | 1.9750 | 0 | 0.5 |
| 132 | **F** | 3 | `20274503` | `n400_ref` | 9,200 | 400 | 9,539 | **0** | 0.5134 | 1.890e-04 | 5.242e-06 | 1.048e-04 | 0.0449 | 0.0097 | 0.5144 | 0.0210 | 0.0603 | 0.4534 | 2.6332 | 0 | 5.9 |
| 133 | **F** | 4 | `20274504` | `completo` | 9,200 | 9,200 | 9,532 | **1** | 0.0000 | 2.339e-12 | 5.245e-06 | 1.049e-04 | 0.0538 | 0.0128 | 0.4931 | 0.0111 | 0.0177 | 0.1237 | 0.4131 | 0 | 11.2 |
| 134 | **F** | 4 | `20274504` | `n400_400` | 400 | 400 | 9,547 | **0** | 0.4048 | 5.153e-05 | 5.237e-06 | 1.047e-04 | 0.0555 | 0.0119 | 0.4979 | 0.0084 | 0.0839 | 0.6409 | 33.3059 | 0 | 0.5 |
| 135 | **F** | 4 | `20274504` | `n400_ref` | 9,200 | 400 | 9,532 | **0** | 0.6125 | 6.426e-05 | 5.245e-06 | 1.049e-04 | 0.0522 | 0.0111 | 0.5053 | 0.0095 | 0.0583 | 0.4479 | 1.6427 | 0 | 5.6 |
| 136 | **F** | 5 | `20274505` | `completo` | 9,200 | 9,200 | 9,530 | **1** | 0.0234 | 2.455e-06 | 5.247e-06 | 1.049e-04 | 0.0654 | 0.0174 | 0.4825 | 0.0252 | 0.0182 | 0.1237 | 0.4501 | 0 | 11.6 |
| 137 | **F** | 5 | `20274505` | `n400_400` | 400 | 400 | 9,515 | **0** | 0.3748 | 5.700e-05 | 5.255e-06 | 1.051e-04 | 0.0521 | 0.0102 | 0.4918 | 0.0115 | 0.0827 | 0.5877 | 2.5421 | 0 | 0.5 |
| 138 | **F** | 5 | `20274505` | `n400_ref` | 9,200 | 400 | 9,530 | **0** | 0.6032 | 8.623e-05 | 5.247e-06 | 1.049e-04 | 0.0558 | 0.0102 | 0.4856 | 0.0187 | 0.0592 | 0.4183 | 1.2237 | 0 | 5.7 |
| 139 | **F** | 6 | `20274506` | `completo` | 9,200 | 9,200 | 9,541 | **0** | 0.0866 | 1.513e-05 | 5.241e-06 | 1.048e-04 | 0.0562 | 0.0139 | 0.4920 | 0.0179 | 0.0181 | 0.1241 | 0.4226 | 0 | 11.0 |
| 140 | **F** | 6 | `20274506` | `n400_400` | 400 | 400 | 9,546 | **0** | 0.9973 | 5.069e-04 | 5.238e-06 | 1.047e-04 | 0.0441 | 0.0083 | 0.5094 | 0.0129 | 0.0795 | 0.6299 | 5.1221 | 0 | 0.4 |
| 141 | **F** | 6 | `20274506` | `n400_ref` | 9,200 | 400 | 9,541 | **0** | 0.9509 | 1.042e-04 | 5.241e-06 | 1.048e-04 | 0.0496 | 0.0085 | 0.4866 | 0.0163 | 0.0608 | 0.4246 | 0.8024 | 0 | 5.7 |
| 142 | **F** | 7 | `20274507` | `completo` | 9,200 | 9,200 | 9,539 | **0** | 0.1896 | 1.987e-05 | 5.242e-06 | 1.048e-04 | 0.0491 | 0.0125 | 0.5005 | 0.0062 | 0.0169 | 0.1288 | 0.8689 | 0 | 11.0 |
| 143 | **F** | 7 | `20274507` | `n400_400` | 400 | 400 | 9,528 | **0** | 0.2027 | 2.128e-05 | 5.248e-06 | 1.049e-04 | 0.0609 | 0.0133 | 0.4905 | 0.0164 | 0.0840 | 0.6330 | 33.3823 | 0 | 0.5 |
| 144 | **F** | 7 | `20274507` | `n400_ref` | 9,200 | 400 | 9,539 | **0** | 0.3801 | 4.847e-05 | 5.242e-06 | 1.048e-04 | 0.0463 | 0.0084 | 0.5042 | 0.0113 | 0.0590 | 0.4473 | 2.6169 | 0 | 5.6 |

**Los diagnosticos del nulo, agregados por brazo** (de las mismas filas de
arriba, para poder cruzarlos):

| brazo | corridas | frac `p<0.05` media (sd) | frac `p<0.01` media (sd) | `p` mediana media | KS media (sd) | `p` min minimo | `|lfc|` mediana media | `|lfc|` max maximo | `lfc` NaN total |
|---|---|---|---|---|---|---|---|---|---|
| `completo` | 48 | 0.0561 (0.0075) | 0.0123 (0.0027) | 0.4897 | 0.0171 (0.0105) | 2.339e-12 | 0.0175 | 1.3157 | 0 |
| `n400_400` | 48 | 0.0498 (0.0058) | 0.0101 (0.0022) | 0.5002 | 0.0145 (0.0076) | 5.641e-06 | 0.0817 | 33.3823 | 0 |
| `n400_ref` | 48 | 0.0498 (0.0062) | 0.0101 (0.0019) | 0.5009 | 0.0146 (0.0073) | 6.087e-07 | 0.0588 | 32.4531 | 0 |

---

## PARTE 2 — Lo que el nulo contesta

### 2.1 🔴 La tasa de falsos positivos del conjunto ED en contexto, con `n = 400`

**Dos cantidades distintas que no hay que confundir.**

* **`n_sig` por corrida**: cuantos genes pasa la puerta en UNA comparacion
  control-contra-control. Es la cantidad comparable con el `n_DE` de una
  perturbacion, y es la que contesta "cuantos genes pasa la puerta sobre puro
  ruido".
* **FPR POR GEN**: `suma(n_sig) / suma(m)`. Es la tasa por prueba. Por
  construccion de BH se espera MUY por debajo de 0.05: 0.05 es el nivel de FDR
  del conjunto llamado, no la tasa por prueba.

**TABLA 2.1 — la tasa de falsos positivos, por brazo.**

| # | brazo | `n` ref | `n` blanco | corridas | **`n_sig` media** | **`n_sig` sd** | mediana | **[min, max]** | corridas con `n_sig>0` | suma `n_sig` | `m` medio | **FPR por gen** | `p_crit` [min, max] | `p` min sobre corridas |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `completo` 9,200 v 9,200 | 9,200 | 9,200 | 48 | **0.7083** | **1.2709** | 0.0 | **[0, 7]** | 18/48 | 34 | 9,766.9 | **7.252e-05** | [4.935e-06, 5.247e-06] | 2.339e-12 |
| 2 | `n400_400` 400 v 400 | 400 | 400 | 48 | **0.0000** | **0.0000** | 0.0 | **[0, 0]** | 0/48 | 0 | 9,757.6 | **0.000e+00** | [4.937e-06, 5.261e-06] | 5.641e-06 |
| 3 | `n400_ref_completa` 400 v 9,200 | 9,200 | 400 | 48 | **0.0208** | **0.1443** | 0.0 | **[0, 1]** | 1/48 | 1 | 9,766.9 | **2.133e-06** | [4.935e-06, 5.247e-06] | 6.087e-07 |

**TABLA 2.2 — `n_sig` por contexto y brazo, con sus 8 valores crudos.**

| # | brazo | ctx | los 8 valores de `n_sig` (j=0..7) | media | sd | **[min, max]** |
|---|---|---|---|---|---|---|
| 1 | `completo` 9,200 v 9,200 | **A** | `[2, 0, 0, 1, 0, 0, 0, 1]` | 0.500 | 0.756 | **[0, 2]** |
| 2 | `completo` 9,200 v 9,200 | **B** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |
| 3 | `completo` 9,200 v 9,200 | **C** | `[0, 0, 0, 0, 1, 0, 2, 1]` | 0.500 | 0.756 | **[0, 2]** |
| 4 | `completo` 9,200 v 9,200 | **D** | `[0, 0, 2, 0, 0, 0, 7, 0]` | 1.125 | 2.475 | **[0, 7]** |
| 5 | `completo` 9,200 v 9,200 | **E** | `[0, 2, 0, 3, 1, 0, 1, 0]` | 0.875 | 1.126 | **[0, 3]** |
| 6 | `completo` 9,200 v 9,200 | **F** | `[2, 2, 1, 3, 1, 1, 0, 0]` | 1.250 | 1.035 | **[0, 3]** |
| 7 | `n400_400` 400 v 400 | **A** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |
| 8 | `n400_400` 400 v 400 | **B** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |
| 9 | `n400_400` 400 v 400 | **C** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |
| 10 | `n400_400` 400 v 400 | **D** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |
| 11 | `n400_400` 400 v 400 | **E** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |
| 12 | `n400_400` 400 v 400 | **F** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |
| 13 | `n400_ref_completa` 400 v 9,200 | **A** | `[0, 0, 0, 1, 0, 0, 0, 0]` | 0.125 | 0.354 | **[0, 1]** |
| 14 | `n400_ref_completa` 400 v 9,200 | **B** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |
| 15 | `n400_ref_completa` 400 v 9,200 | **C** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |
| 16 | `n400_ref_completa` 400 v 9,200 | **D** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |
| 17 | `n400_ref_completa` 400 v 9,200 | **E** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |
| 18 | `n400_ref_completa` 400 v 9,200 | **F** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | 0.000 | **[0, 0]** |

**Lo que dicen esas filas, como aritmetica:**

* Con **`n = 400` a los dos lados**: `n_sig` vale **0 en las
  48 corridas** (0 de
  48 con alguno), media **0.0000**, sd
  **0.0000**, FPR por gen **0.000e+00** sobre
  0 genes llamados de
  468,363 pruebas.
* Con **400 contra el pool de 9,200** (el brazo que se parece a la entrega):
  media **0.0208**, sd **0.1443**, rango
  **[0, 1]**, 1/48 corridas con alguno.
* Con **9,200 contra 9,200**: media **0.7083**, sd
  **1.2709**, rango **[0, 7]**,
  18/48 corridas con alguno.

**El diagnostico de uniformidad de la seccion 1.4** cruza con esto: la fraccion
de `p` crudas bajo 0.05 promedia 0.0501-0.0539 segun el brazo (esperado 0.05
bajo un nulo uniforme) y el KS contra la uniforme promedia 0.013-0.014.

### 2.2 🔴 El piso de deteccion, calculado del propio nulo

**La construccion, declarada.** La puerta llama un gen cuando `p_adj < 0.05`.
Para una senal **AISLADA** (el resto del panel nulo), el escalon de BH de rango
1 exige `p <= alpha/m`. Ese `p_crit` se calcula de cada corrida (`m` = genes
probados tras el filtro CPM), no se importa de ninguna formula.

El piso depende FUERTE del nivel de expresion del gen, asi que se mide por
banda de CPM medio de la referencia. Para cada banda se dan:

* el `|log2FC|` MINIMO que el nulo alcanzo en cada umbral de `p` -- lectura
  directa, sin modelo;
* el **techo de ruido**: el `|log2FC|` maximo que el nulo PRODUCE en esa banda;
* el **alcance**: `n_filas x p_crit` dice cuantas filas del nulo se esperan por
  debajo de `p_crit`. Si es << 1, el piso en `p_crit` **no se puede leer
  directo**.

#### Brazo `n400_400` 400 v 400

Pool de **468,363 filas** (48 corridas x ~9,758 genes).
`p_crit = 0.05/9,757.6 = 5.1242e-06`. `p` minimo del pool = **5.6412e-06**.

🔴 **El alcance.** Filas esperadas bajo `p_crit`: **2.40**. Observadas: **0**. `p_crit` **NO** cae dentro del alcance directo de este nulo.

**TABLA 2.3a — el piso de deteccion por banda de expresion.**

`F(t) = min{|log2FC| : p <= t}` sobre el pool, que es **no creciente en
`t`**. Como `p_crit` es menor que cada ancla alcanzada, **cada ancla da una
cota INFERIOR medida** del piso. Una fila del nulo que SI cruzo `p_crit` da
una cota SUPERIOR medida.

| # | banda CPM ref | filas | `|lfc|` mediana | **techo de ruido (`|lfc|` max)** | `p` min | `F(1e-2)` (n) | `F(1e-3)` (n) | `F(1e-4)` (n) | `F(1e-5)` (n) | **cota INFERIOR medida** | **cota SUPERIOR medida** | filas esp. bajo `p_crit` | filas obs. | **extrapolacion a `p_crit`** | pend./decada | decadas extrap. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | (5, 10] | 43,151 | 0.2204 | **33.1910** | 1.32e-05 | 0.1269 (448) | 0.1878 (60) | 0.7266 (7) | — (0) | **>= 0.7266** (de `p<=1e-04`, n=7) | **n/d** -- ninguna fila cruzo | 0.221 | 0 | `|lfc|` = **1.0339** (EXTRAPOLADO) | 0.2998 | 1.29 |
| 2 | (10, 20] | 61,444 | 0.1554 | **33.3823** | 2.17e-05 | 0.0980 (674) | 0.4021 (46) | 0.6786 (3) | — (0) | **>= 0.6786** (de `p<=1e-04`, n=3) | **n/d** -- ninguna fila cruzo | 0.315 | 0 | `|lfc|` = **1.0578** (EXTRAPOLADO) | 0.2903 | 1.29 |
| 3 | (20, 50] | 123,921 | 0.1038 | **5.0902** | 5.64e-06 | 0.0705 (1224) | 0.1998 (134) | 0.4294 (19) | 0.5352 (2) | **>= 0.5352** (de `p<=1e-05`, n=2) | **n/d** -- ninguna fila cruzo | 0.635 | 0 | `|lfc|` = **0.5994** (EXTRAPOLADO) | 0.1624 | 0.29 |
| 4 | (50, 100] | 105,965 | 0.0723 | **3.9241** | 1.43e-05 | 0.0161 (1057) | 0.1412 (107) | 0.3216 (11) | — (0) | **>= 0.3216** (de `p<=1e-04`, n=11) | **n/d** -- ninguna fila cruzo | 0.543 | 0 | `|lfc|` = **0.5094** (EXTRAPOLADO) | 0.1527 | 1.29 |
| 5 | (100, 300] | 107,406 | 0.0505 | **1.8864** | 9.50e-06 | 0.0723 (1087) | 0.1430 (101) | 0.2107 (6) | 0.2627 (1) | **>= 0.2627** (de `p<=1e-05`, n=1) | **n/d** -- ninguna fila cruzo | 0.550 | 0 | `|lfc|` = **0.2866** (EXTRAPOLADO) | 0.0639 | 0.29 |
| 6 | (300, 1000] | 24,023 | 0.0355 | **1.1821** | 2.13e-05 | 0.0724 (214) | 0.1231 (19) | 0.1674 (3) | — (0) | **>= 0.1674** (de `p<=1e-04`, n=3) | **n/d** -- ninguna fila cruzo | 0.123 | 0 | `|lfc|` = **0.2298** (EXTRAPOLADO) | 0.0475 | 1.29 |
| 7 | > 1000 | 2,453 | 0.0288 | **0.3250** | 8.29e-04 | 0.0586 (22) | 0.0808 (1) | — (0) | — (0) | **>= 0.0808** (de `p<=1e-03`, n=1) | **n/d** -- ninguna fila cruzo | 0.013 | 0 | n/d (menos de 3 anclas) | — | — |

#### Brazo `n400_ref_completa` 400 v 9,200

Pool de **468,812 filas** (48 corridas x ~9,767 genes).
`p_crit = 0.05/9,766.9 = 5.1193e-06`. `p` minimo del pool = **6.0868e-07**.

🔴 **El alcance.** Filas esperadas bajo `p_crit`: **2.40**. Observadas: **1**. `p_crit` **SI** cae dentro del alcance directo de este nulo.

Las filas del nulo que SI bajaron de `p_crit` llevan `|log2FC|` = [0.3573] con CPM de referencia [131.64].

**TABLA 2.3b — el piso de deteccion por banda de expresion.**

`F(t) = min{|log2FC| : p <= t}` sobre el pool, que es **no creciente en
`t`**. Como `p_crit` es menor que cada ancla alcanzada, **cada ancla da una
cota INFERIOR medida** del piso. Una fila del nulo que SI cruzo `p_crit` da
una cota SUPERIOR medida.

| # | banda CPM ref | filas | `|lfc|` mediana | **techo de ruido (`|lfc|` max)** | `p` min | `F(1e-2)` (n) | `F(1e-3)` (n) | `F(1e-4)` (n) | `F(1e-5)` (n) | **cota INFERIOR medida** | **cota SUPERIOR medida** | filas esp. bajo `p_crit` | filas obs. | **extrapolacion a `p_crit`** | pend./decada | decadas extrap. |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | (5, 10] | 43,340 | 0.1576 | **32.4531** | 2.74e-05 | 0.0105 (437) | 0.3515 (41) | 0.6143 (6) | — (0) | **>= 0.6143** (de `p<=1e-04`, n=6) | **n/d** -- ninguna fila cruzo | 0.222 | 0 | `|lfc|` = **1.0171** (EXTRAPOLADO) | 0.3019 | 1.29 |
| 2 | (10, 20] | 61,364 | 0.1117 | **4.1428** | 7.93e-06 | 0.1250 (642) | 0.2580 (65) | 0.4007 (4) | 0.7721 (1) | **>= 0.7721** (de `p<=1e-05`, n=1) | **n/d** -- ninguna fila cruzo | 0.314 | 0 | `|lfc|` = **0.7621** (EXTRAPOLADO) | 0.2084 | 0.29 |
| 3 | (20, 50] | 124,135 | 0.0747 | **2.2221** | 1.25e-05 | 0.0136 (1239) | 0.1854 (139) | 0.2463 (11) | — (0) | **>= 0.2463** (de `p<=1e-04`, n=11) | **n/d** -- ninguna fila cruzo | 0.635 | 0 | `|lfc|` = **0.4149** (EXTRAPOLADO) | 0.1163 | 1.29 |
| 4 | (50, 100] | 105,876 | 0.0520 | **0.7761** | 1.81e-05 | 0.0056 (1086) | 0.0930 (110) | 0.1564 (10) | — (0) | **>= 0.1564** (de `p<=1e-04`, n=10) | **n/d** -- ninguna fila cruzo | 0.542 | 0 | `|lfc|` = **0.2578** (EXTRAPOLADO) | 0.0754 | 1.29 |
| 5 | (100, 300] | 107,661 | 0.0365 | **0.9790** | 6.09e-07 | 0.0269 (1090) | 0.0269 (106) | 0.1579 (13) | 0.2565 (3) | **>= 0.2565** (de `p<=1e-05`, n=3) | **<= 0.3573** (n=1) | 0.551 | 1 | `|lfc|` = **0.2639** (EXTRAPOLADO) | 0.0820 | 0.29 |
| 6 | (300, 1000] | 23,989 | 0.0256 | **0.8337** | 2.58e-05 | 0.0617 (234) | 0.0877 (23) | 0.0938 (2) | — (0) | **>= 0.0938** (de `p<=1e-04`, n=2) | **n/d** -- ninguna fila cruzo | 0.123 | 0 | `|lfc|` = **0.1178** (EXTRAPOLADO) | 0.0160 | 1.29 |
| 7 | > 1000 | 2,447 | 0.0207 | **0.2568** | 1.39e-03 | 0.0481 (16) | — (0) | — (0) | — (0) | **>= 0.0481** (de `p<=1e-02`, n=16) | **n/d** -- ninguna fila cruzo | 0.013 | 0 | n/d (menos de 3 anclas) | — | — |

🔴 **Lo que se puede y lo que no se puede calcular de este nulo.**

1. **SE PUEDE, y esta medido: una cota INFERIOR del piso en cada banda.**
   `F` es no creciente y `p_crit` esta por debajo de cada ancla alcanzada, asi
   que `F(ancla)` acota el piso por abajo sin ningun modelo. Con `n = 400` a
   los dos lados va de
   **0.081**
   (banda de CPM mas alto) a
   **0.727**
   (banda de CPM mas bajo) en `|log2FC|`.
2. **SE PUEDE, en UNA banda de UN brazo: un intervalo con las dos cotas**
   **medidas.** En el brazo 400 v 9,200, banda
   CPM (100, 300], una fila del nulo cruzo `p_crit`
   (`p = 6.087e-07`), asi que el piso queda **entre
   0.257 y 0.357** en
   `|log2FC|`, las dos cotas medidas. La extrapolacion de esa misma banda da
   0.264, dentro del intervalo.
3. 🔴 **NO SE PUEDE leer directo el piso EXACTO en `p_crit`.** Las filas
   esperadas por debajo de `p_crit` son ~0.01-0.64 por banda (columna de la
   tabla): el nulo, como se construyo, no llega ahi. Donde se da una cifra, es
   una **extrapolacion DECLARADA** -- `F(t)` lineal en `log10(t)` ajustada a las
   anclas MEDIDAS de esa banda -- y va marcada `(EXTRAPOLADO)` con su pendiente
   por decada y cuantas decadas se extrapola (0.29 donde se alcanzo `p<=1e-5`,
   1.29 donde solo se alcanzo `p<=1e-4`). **Es un rango, no un numero medido.**
   La forma log-lineal **es** una suposicion importada; el nulo solo da los puntos.
4. **Medido directo: el techo de ruido por banda** -- el `|log2FC|` maximo que
   el nulo PRODUCE. 🔴 En la banda de CPM mas bajo llega a
   **33.38** con `n = 400`: un
   `|log2FC|` enorme es compatible con puro ruido cuando el gen tiene ceros en
   un lado (`epsilon = 1e-9` contra una media de CPM positiva). **Ningun umbral
   de `|log2FC|` solo separa senal de ruido a `n = 400`**; la potencia de la
   puerta viene del test de rangos, no del tamano de efecto. Esa cola gruesa es
   tambien la razon de que la extrapolacion de las dos bandas mas bajas tenga la
   pendiente por decada mas alta (~0.30 contra ~0.05 en la banda mas alta).

### 2.3 La variacion guia a guia: cuanto del nulo lo explica la identidad de la guia

**La construccion, declarada.** Dos mediciones, ninguna parte una guia:

1. **Particion que AGRUPA guias contra particiones ALEATORIAS de guias.** La
   particion aleatoria de la Parte 1 reparte las 46 guias al azar. La que
   AGRUPA las ordena por su propio perfil y corta por la mitad: las 23 de un
   extremo contra las 23 del otro. Es el peor caso que la estructura de guia
   permite **sin partir ninguna guia**.
   * perfil de guia = media del CPM por celula sobre sus 400 celulas,
     restringido al panel probable, `log1p`, centrado por gen;
   * se ordena por la **primera componente principal** de esa matriz 46 x G;
   * mitad REFERENCIA = las 23 de PC1 mas baja, BLANCO = las 23 mas alta;
   * **determinista: no lleva semilla**. Solo el submuestreo a 400 la lleva
     (`20274500 + idx_contexto`).
2. **Descomposicion de varianza por gen** (ANOVA de un factor, 46 guias x 400
   celulas), en `log1p(CPM)` -- el espacio que ve el Wilcoxon antes de pasar a
   rangos:
   `MSB = SSB/45`, `MSW = SSW/(18400-46)`, `F = MSB/MSW`,
   `ICC = (MSB - MSW)/(MSB + 399*MSW)`.
   `F` dice si la identidad de la guia mueve la media del gen mas de lo que la
   mueve el muestreo de 400 celulas; `ICC` dice QUE FRACCION de la varianza
   por celula se lleva la guia.

**TABLA 2.4 — la particion que agrupa guias contra las 8 aleatorias, por contexto.**

| # | ctx | brazo | `n_sig` AGRUPA guias | `n_sig` aleatorias: media (sd) | aleatorias [min, max] | `m` AGRUPA | `p` min AGRUPA | `p` min aleatorias (mediana) | `|lfc|` max AGRUPA | genes llamados por AGRUPA |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **A** | `completo` 9,200 v 9,200 | **4** | 0.500 (0.756) | [0, 2] | 9,924 | 1.890e-07 | 2.728e-05 | 0.2447 | `SMARCD1`, `METTL8`, `ATP6V1E1`, `CARF` |
| 2 | **A** | `n400_400` 400 v 400 | **0** | 0.000 (0.000) | [0, 0] | 9,936 | 1.449e-04 | 1.245e-04 | 32.4436 | — |
| 3 | **A** | `n400_ref_completa` 400 v 9,200 | **0** | 0.125 (0.354) | [0, 1] | 9,924 | 1.803e-04 | 7.133e-05 | 0.7263 | — |
| 4 | **B** | `completo` 9,200 v 9,200 | **8** | 0.000 (0.000) | [0, 0] | 9,639 | 5.469e-08 | 3.997e-05 | 1.0599 | `ZC3H3`, `PLCB4`, `CST3`, `CEP76`, `FST`, `WTAP`, `CITED2`, `TMEM182` |
| 5 | **B** | `n400_400` 400 v 400 | **0** | 0.000 (0.000) | [0, 0] | 9,621 | 1.955e-04 | 8.673e-05 | 32.7386 | — |
| 6 | **B** | `n400_ref_completa` 400 v 9,200 | **1** | 0.000 (0.000) | [0, 0] | 9,639 | 3.272e-06 | 6.418e-05 | 5.9564 | `ATXN1` |
| 7 | **C** | `completo` 9,200 v 9,200 | **0** | 0.500 (0.756) | [0, 2] | 10,134 | 1.253e-05 | 3.208e-05 | 1.0993 | — |
| 8 | **C** | `n400_400` 400 v 400 | **0** | 0.000 (0.000) | [0, 0] | 10,132 | 8.394e-05 | 5.909e-05 | 2.9742 | — |
| 9 | **C** | `n400_ref_completa` 400 v 9,200 | **0** | 0.000 (0.000) | [0, 0] | 10,134 | 9.385e-05 | 1.145e-04 | 2.5490 | — |
| 10 | **D** | `completo` 9,200 v 9,200 | **1937** | 1.125 (2.475) | [0, 7] | 9,622 | 1.102e-12 | 4.740e-05 | 0.4853 | *1,937 genes: lista omitida (>200)* |
| 11 | **D** | `n400_400` 400 v 400 | **1** | 0.000 (0.000) | [0, 0] | 9,587 | 6.106e-07 | 8.534e-05 | 2.6891 | `LRRC58` |
| 12 | **D** | `n400_ref_completa` 400 v 9,200 | **0** | 0.000 (0.000) | [0, 0] | 9,622 | 6.939e-05 | 7.662e-05 | 0.9016 | — |
| 13 | **E** | `completo` 9,200 v 9,200 | **15** | 0.875 (1.126) | [0, 3] | 9,742 | 1.924e-06 | 7.250e-06 | 0.4459 | `STPG1`, `DEPDC1B`, `OXCT1`, `NCOA1`, `ZW10`, `CST3`, `TBP`, `LAMA4`, `LDLR`, `PRMT7`, `CRTC2`, `TMEM154`, `SMAGP`, `GSPT2`, `MAP3K5` |
| 14 | **E** | `n400_400` 400 v 400 | **1** | 0.000 (0.000) | [0, 0] | 9,736 | 1.053e-06 | 1.177e-04 | 2.9059 | `SKP2` |
| 15 | **E** | `n400_ref_completa` 400 v 9,200 | **0** | 0.000 (0.000) | [0, 0] | 9,742 | 2.525e-05 | 4.721e-05 | 1.6745 | — |
| 16 | **F** | `completo` 9,200 v 9,200 | **12** | 1.250 (1.035) | [0, 3] | 9,544 | 2.639e-07 | 3.214e-06 | 0.5469 | `TFB1M`, `ITGA6`, `SERPINE1`, `LYST`, `COMMD10`, `TNIK`, `ADAMTS5`, `SLC27A4`, `SCN9A`, `MARCH3`, `ECSCR`, `RTEL1` |
| 17 | **F** | `n400_400` 400 v 400 | **0** | 0.000 (0.000) | [0, 0] | 9,508 | 5.218e-05 | 4.050e-05 | 1.7619 | — |
| 18 | **F** | `n400_ref_completa` 400 v 9,200 | **0** | 0.000 (0.000) | [0, 0] | 9,544 | 4.230e-05 | 9.524e-05 | 1.8065 | — |

**TABLA 2.5 — la descomposicion de varianza por gen, por contexto.**
`n` genes evaluados = los del panel probable con `MSW > 0`.

| # | ctx | panel probable | genes evaluados | **genes con guia significativa (BH 0.05)** | **frac** | `F` mediana | `F` p99 | `F` max | **`ICC` mediana** | `ICC` p75 | `ICC` p99 | `ICC` max | frac `ICC<0` | PC1 % var | gl entre | gl dentro |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **A** | 9,929 | 9,929 | **162** | **0.01632** | 1.0327 | 1.9688 | 18.6673 | **0.000082** | 0.000485 | 0.002416 | 0.042300 | 0.4393 | 7.80 | 45 | 18,354 |
| 2 | **B** | 9,626 | 9,626 | **11** | **0.00114** | 0.9909 | 1.6036 | 18.1370 | **-0.000023** | 0.000347 | 0.001507 | 0.041082 | 0.5192 | 9.08 | 45 | 18,354 |
| 3 | **C** | 10,124 | 10,124 | **6** | **0.00059** | 0.9882 | 1.5598 | 25.2563 | **-0.000030** | 0.000341 | 0.001398 | 0.057174 | 0.5218 | 6.64 | 45 | 18,354 |
| 4 | **D** | 9,633 | 9,633 | **59** | **0.00612** | 1.0826 | 1.8225 | 3.0679 | **0.000207** | 0.000656 | 0.002052 | 0.005143 | 0.3676 | 15.96 | 45 | 18,354 |
| 5 | **E** | 9,743 | 9,743 | **41** | **0.00421** | 1.0253 | 1.6890 | 36.5471 | **0.000063** | 0.000453 | 0.001720 | 0.081615 | 0.4560 | 7.95 | 45 | 18,354 |
| 6 | **F** | 9,537 | 9,537 | **38** | **0.00398** | 1.0309 | 1.7025 | 38.3222 | **0.000077** | 0.000476 | 0.001753 | 0.085343 | 0.4477 | 9.20 | 45 | 18,354 |

**TABLA 2.6 — el contraste directo: AGRUPA contra aleatorias, por brazo y contexto.**

| # | brazo | ctx | **`n_sig` AGRUPA** | `n_sig` aleatorias (los 8) | aleatorias media | **AGRUPA - media aleatoria** |
|---|---|---|---|---|---|---|
| 1 | `completo` 9,200 v 9,200 | **A** | **4** | `[2, 0, 0, 1, 0, 0, 0, 1]` | 0.500 | **+3.500** |
| 2 | `completo` 9,200 v 9,200 | **B** | **8** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+8.000** |
| 3 | `completo` 9,200 v 9,200 | **C** | **0** | `[0, 0, 0, 0, 1, 0, 2, 1]` | 0.500 | **-0.500** |
| 4 | `completo` 9,200 v 9,200 | **D** | **1,937** | `[0, 0, 2, 0, 0, 0, 7, 0]` | 1.125 | **+1,935.875** |
| 5 | `completo` 9,200 v 9,200 | **E** | **15** | `[0, 2, 0, 3, 1, 0, 1, 0]` | 0.875 | **+14.125** |
| 6 | `completo` 9,200 v 9,200 | **F** | **12** | `[2, 2, 1, 3, 1, 1, 0, 0]` | 1.250 | **+10.750** |
| 7 | `n400_400` 400 v 400 | **A** | **0** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+0.000** |
| 8 | `n400_400` 400 v 400 | **B** | **0** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+0.000** |
| 9 | `n400_400` 400 v 400 | **C** | **0** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+0.000** |
| 10 | `n400_400` 400 v 400 | **D** | **1** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+1.000** |
| 11 | `n400_400` 400 v 400 | **E** | **1** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+1.000** |
| 12 | `n400_400` 400 v 400 | **F** | **0** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+0.000** |
| 13 | `n400_ref_completa` 400 v 9,200 | **A** | **0** | `[0, 0, 0, 1, 0, 0, 0, 0]` | 0.125 | **-0.125** |
| 14 | `n400_ref_completa` 400 v 9,200 | **B** | **1** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+1.000** |
| 15 | `n400_ref_completa` 400 v 9,200 | **C** | **0** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+0.000** |
| 16 | `n400_ref_completa` 400 v 9,200 | **D** | **0** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+0.000** |
| 17 | `n400_ref_completa` 400 v 9,200 | **E** | **0** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+0.000** |
| 18 | `n400_ref_completa` 400 v 9,200 | **F** | **0** | `[0, 0, 0, 0, 0, 0, 0, 0]` | 0.000 | **+0.000** |

**Lo que esas tres tablas dicen, como aritmetica:**

* En el brazo `completo`, la particion que AGRUPA guias llama
  **A=4, B=8, C=0, D=1,937, E=15, F=12** genes,
  contra una media de **0.7083** (rango
  [0, 7]) de
  las 48 particiones aleatorias. 🔴 **El contexto D llega a
  1,937 genes de 9,622**
  probados, y es tambien el contexto con mas varianza entre guias capturada por
  PC1 (15.96 %).
* En el brazo `n400_400`, donde las 48 aleatorias dan **0 en 48**, la particion
  que AGRUPA da **A=0, B=0, C=0, D=1, E=1, F=0**.
* El `F` mediano por gen esta en
  [0.9882, 1.0826] sobre los seis
  contextos. Un `F` mediano de 1 es exactamente "la guia no explica nada mas
  que el muestreo".
* Los genes con efecto de guia significativo a BH 0.05 van de
  **6**
  a **162**
  por contexto, es decir
  **0.059-1.632 %**
  del panel probable.
* El `ICC` mediano va de
  **-0.000030** a
  **0.000207** -- la guia se
  lleva una fraccion diminuta de la varianza POR CELULA de un gen tipico -- y su
  maximo por contexto llega a
  **0.0853**.
* 🔴 **El contexto D tiene el `F` MAXIMO mas bajo de los seis y a la vez las**
  **1,937 llamadas.** `F_max` en D vale
  **3.0679** y su `ICC_max` **0.005143**, los mas bajos de
  la TABLA 2.5 (los otros contextos llegan a `F_max` de 18.1-38.3).
  Es decir: en D **ningun gen suelto** tiene un efecto de guia grande, pero PC1
  captura **15.96 %**
  de la varianza entre guias, el doble que en cualquier otro contexto. La
  estructura de guia de D es **ancha y de poca amplitud** -- un desplazamiento
  pequeno repartido sobre muchos genes -- y eso es exactamente lo que un corte a
  lo largo de PC1 convierte en llamadas masivas, mientras que un ANOVA gen a gen
  no lo ve en ningun gen en particular. Las dos tablas miran la misma estructura
  por ejes distintos.
* 🔴 **Las dos mediciones no se contradicen, y la aritmetica de por que es**
  **esta:** el `ICC` por gen es diminuto, pero la puerta no compara celulas
  sueltas: compara las MEDIAS de dos grupos de 9,200. El error de una media de
  9,200 celulas es ~96 veces menor que la `sd` por celula, asi que un desplazamiento
  de medias que vale una fraccion minuscula de la varianza por celula sigue
  siendo muchas veces su propio error estandar. Esa es la razon de que `ICC`
  mediano ~1e-4 convivan con 1,937 genes llamados cuando las guias se agrupan.

### 2.4 🔴 Cuantas de las 300 perturbaciones esperarian conjunto ED vacio por potencia

**El nulo NO permite contestarlo, y no se estima.**

La razon es estructural, no de precision. El nulo mide una sola de las dos
mitades de la pregunta:

| | lo que hace falta | lo que el nulo da |
|---|---|---|
| **el umbral** | el tamano de efecto que la puerta exige a `n = 400` con esta profundidad | **SI** -- las secciones 2.1 y 2.2, con su rango y su extrapolacion declarada |
| **la distribucion** | el tamano de efecto REAL de cada una de las 300 perturbaciones en cada contexto | **NO** -- esta en el lado perturbado |

Contar cuantas perturbaciones caen por debajo del umbral exige las dos. El
lado perturbado de 2026 **no existe**: `manifest.json` declara
`ground_truth_cells: 138,400` por contexto, pero lo que se entrega son las
18,400 celulas de control (TABLA 0.2); las 120,000 celulas perturbadas no
estan en `/datos/reto`. Y el de 2025/Jiang esta fuera de este encargo.

Sustituir la distribucion que falta por la de Jiang seria **suponer que el
regimen transfiere**, que es exactamente lo que la Parte 3 se niega a declarar
en cualquiera de los dos sentidos. Asi que aqui no hay cifra, y no se da una.

Lo unico que el nulo acota sin verdad es la direccion contraria: **cuantas
perturbaciones esperarian ED vacio POR RUIDO DEL CONTROL**, que es
`n_sig` del brazo correspondiente -- **0.0000** genes de media
a `n = 400` contra `n = 400`. Eso es el piso de falsos positivos, no el conteo
de verdad vacia.

---

## PARTE 3 — Contra el regimen de Jiang

**La referencia que se compara.** `control_permutado.md` §7 MIDIO, en las seis
lineas de Jiang, el numero de perturbaciones con DE REAL vacio: **A549 15,
BXPC3 0, HAP1 14, HT29 6, K562 8, MCF7 0**, de 33 -- el rango 0 a 15 de 33 que
cita el encargo. En el reto, `n = 400` es fijo por `manifest.json`
(`cells_per_pert: 400`) y el knockdown es >=80 % por criterio de seleccion.

**La construccion, declarada.** Las dos cantidades se miden con la MISMA
definicion en los dos lados:

* **profundidad por celula** = suma de cuentas de la fila (UMIs). En Jiang se
  cruzo contra `obs/nCount_RNA` (coincide fila por fila en las seis lineas);
  en el reto se suma del propio CSR.
* **`nnz` por celula** = `diff(indptr)`. En Jiang se cruzo contra
  `obs/nFeature_RNA` (identico: `True`).
* **panel probable** = nº de genes cuyo **CPM medio por celula en el pool de
  control** supera `5.0` -- la MISMA cantidad con la que `_apply_cpm_filter`
  del preset decide que genes entran en BH. Es el `m` de BH.

De Jiang se leyeron **solo las celulas `NT`** de cada linea. El lado
perturbado no se toco.

**TABLA 3.1 — profundidad y `nnz` de los controles: reto contra Jiang, por linea.**

| # | fuente | unidad | celulas de control | genes del panel | prof. mediana | prof. media | prof. p1 | prof. p99 | prof. sd | `nnz` mediana | `nnz` media | `nnz` p99 | `nnz` max | `nnz` sd | **panel probable (CPM>5)** | frac del panel | `nnz`/panel |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | reto 2026 | **A** | 18,400 | 18,533 | 20,109 | 21,134 | 4,470 | 47,490 | 9,705 | 6,147 | 5,973 | 8,003 | 8,453 | 1,186 | **9,929** | 0.5357 | 0.3223 |
| 2 | reto 2026 | **B** | 18,400 | 18,533 | 19,946 | 19,990 | 1,007 | 38,812 | 7,701 | 5,756 | 5,519 | 7,178 | 7,649 | 1,189 | **9,626** | 0.5194 | 0.2978 |
| 3 | reto 2026 | **C** | 18,400 | 18,533 | 20,034 | 21,157 | 4,685 | 46,623 | 9,256 | 6,006 | 5,881 | 7,944 | 8,482 | 1,144 | **10,124** | 0.5463 | 0.3173 |
| 4 | reto 2026 | **D** | 18,400 | 18,533 | 20,236 | 21,250 | 2,023 | 48,413 | 10,059 | 5,892 | 5,679 | 7,673 | 8,146 | 1,255 | **9,633** | 0.5198 | 0.3064 |
| 5 | reto 2026 | **E** | 18,400 | 18,533 | 19,681 | 20,921 | 4,117 | 46,983 | 9,529 | 5,634 | 5,522 | 7,606 | 8,107 | 1,140 | **9,743** | 0.5257 | 0.2980 |
| 6 | reto 2026 | **F** | 18,400 | 18,533 | 19,845 | 21,004 | 4,266 | 47,787 | 9,762 | 5,630 | 5,481 | 7,423 | 7,914 | 1,108 | **9,537** | 0.5146 | 0.2958 |
| 7 | Jiang IFNB | **A549** | 1,379 | 34,025 | 3,319 | 4,056 | 784 | 14,747 | 2,999 | 1,952 | 2,088 | 4,972 | 6,480 | 958 | **11,733** | 0.3448 | 0.0614 |
| 8 | Jiang IFNB | **BXPC3** | 3,228 | 34,025 | 4,754 | 5,798 | 919 | 21,997 | 4,343 | 2,710 | 2,839 | 6,283 | 9,272 | 1,203 | **12,005** | 0.3528 | 0.0834 |
| 9 | Jiang IFNB | **HAP1** | 3,172 | 34,025 | 5,231 | 6,248 | 1,063 | 22,340 | 4,338 | 3,062 | 3,181 | 6,688 | 9,482 | 1,256 | **12,041** | 0.3539 | 0.0935 |
| 10 | Jiang IFNB | **HT29** | 2,736 | 34,025 | 4,396 | 5,285 | 894 | 19,104 | 3,762 | 2,493 | 2,597 | 5,739 | 8,048 | 1,106 | **11,181** | 0.3286 | 0.0763 |
| 11 | Jiang IFNB | **K562** | 1,398 | 34,025 | 4,302 | 5,201 | 801 | 16,212 | 3,456 | 2,482 | 2,582 | 5,395 | 6,743 | 1,101 | **11,750** | 0.3453 | 0.0759 |
| 12 | Jiang IFNB | **MCF7** | 2,669 | 34,025 | 6,203 | 8,745 | 1,141 | 43,669 | 8,636 | 3,128 | 3,388 | 8,239 | 9,959 | 1,593 | **11,534** | 0.3390 | 0.0996 |

### 3.1 🔴 Las cantidades, lado a lado. NO se concluye si la verdad vacia transfiere.

**TABLA 3.2 — las razones entre los dos regimenes, calculadas de la TABLA 3.1.**

| # | cantidad | reto (6 contextos) [min, max] | Jiang (6 lineas) [min, max] | razon reto/Jiang [min, max] |
|---|---|---|---|---|
| 1 | celulas de control por unidad | [18,400, 18,400] | [1,379, 3,228] | [x5.70, x13.34] |
| 2 | profundidad mediana (UMIs/celula) | [19,681, 20,236] | [3,319, 6,203] | [x3.17, x6.10] |
| 3 | `nnz` mediana por celula | [5,630, 6,147] | [1,952, 3,128] | [x1.80, x3.15] |
| 4 | genes del panel entregado | [18,533, 18,533] | [34,025, 34,025] | [x0.54, x0.54] |
| 5 | **panel probable (el `m` de BH)** | [9,537, 10,124] | [11,181, 12,041] | [x0.79, x0.91] |

**Y la heterogeneidad DENTRO de cada lado**, que es la cantidad que el encargo
nombra como causa del rango 0-15 en Jiang:

| # | cantidad | reto: max/min entre contextos | Jiang: max/min entre lineas |
|---|---|---|---|
| 1 | celulas de control | **x1.000** | **x2.341** |
| 2 | profundidad mediana | **x1.028** | **x1.869** |
| 3 | `nnz` mediana | **x1.092** | **x1.602** |
| 4 | panel probable | **x1.062** | **x1.077** |

**Las cifras que quedan en la mesa, sin lectura:**

1. El `n` de control por unidad: el reto lleva **18,400 en los seis**
   contextos, exactamente iguales. Jiang lleva de **1,379** (A549) a
   **3,228** (BXPC3) celulas `NT`, una dispersion de
   **x2.34** entre lineas. El pool de control del reto
   es de **x5.7 a x13.3** el
   de una linea de Jiang.
2. La profundidad por celula: el reto esta en
   **19,681-20,236** UMIs de mediana, Jiang en
   **3,319-6,203**. El reto es
   **x3.17 a x6.10** mas profundo.
3. El `nnz` por celula: **5,630-6,147** en el reto
   contra **1,952-3,128** en Jiang, razon
   **x1.80 a x3.15**.
4. 🔴 **El panel probable va al reves.** Aunque el reto es mas profundo, su
   panel entregado tiene 18,533 genes y Jiang 34,025, asi que los genes que
   pasan el filtro CPM>5 -- el `m` sobre el que corre BH -- son
   **9,537-10,124** en el reto y **11,181-12,041**
   en Jiang: **x0.792 a x0.905**.
   Jiang prueba MAS genes, asi que su escalon de BH de rango 1 es MAS exigente
   (`0.05/12,041 = 4.152e-06` contra
   `0.05/9,537 = 5.243e-06`).
5. La heterogeneidad dentro del lado: en el reto la profundidad mediana varia
   **x1.028** entre contextos y el `n` de control
   **x1.000**; en Jiang la profundidad varia **x1.869** y
   el `n` **x2.341**.

🔴 **No se declara si la verdad vacia de Jiang transfiere ni si no**. Las
cantidades de arriba son lo que se midio. Lo que falta para cerrar esa
pregunta -- los tamanos de efecto reales de las 300 perturbaciones del reto --
esta en el lado perturbado, que no existe en 2026; la seccion 2.4 lo detalla.

---

## Defectos de la construccion, como datos sobre la medicion

1. 🔴 **El `n` del nulo no es el `n` de la entrega, en ninguno de los tres**
   **brazos.** El reto compara las 400 celulas de una perturbacion contra el
   pool de control de **18,400** celulas del contexto. Al partir el control en
   dos mitades, la referencia maxima disponible es **9,200**, la mitad. El brazo
   `n400_ref_completa` es el mas cercano (400 contra 9,200) y aun asi su
   referencia es la mitad de la real. El brazo `n400_400` que pide el encargo
   esta mas lejos todavia. Un nulo con referencia de 18,400 **exigiria no**
   **partir el control**, y entonces no habria dos grupos independientes de
   celulas no dirigidas que comparar: es una limitacion del diseno, no un
   parametro que se pueda subir.
2. 🔴 **La puerta es ASIMETRICA y eso no se cuantifico.** El filtro de genes
   del preset lo decide **el grupo de referencia solo** (`_apply_cpm_filter`,
   #351) y BH se recalcula sobre los genes que quedan. Asi que intercambiar las
   dos mitades da un `m` distinto y puede dar un `n_sig` distinto. **Aqui
   siempre hizo de referencia la primera mitad de la permutacion**, y la
   permutacion inversa NO se corrio. El efecto se ve en el rango de `m` de la
   TABLA 1.2 pero no esta aislado.
   Para dimensionarlo: `m` varia de **9,504** a **10,131** sobre las
   144 corridas, un rango de 627 genes (**6.42 %** de la media).
3. 🔴 **El backend es `pdex`, no `gpudge`.** `de.backend='auto'` y la maquina no
   tiene CUDA, asi que la puerta corrio en el camino de CPU. `cell_eval2` avisa
   explicitamente que "DE numbers differ between engines". Si el marcador
   oficial puntua con `gpudge`, las cifras de `n_sig` de este reporte son del
   motor equivocado, y **cuanto difieren no se midio**.
4. **Ocho particiones por contexto es poco para la cola.** Con 8 particiones el
   `sd` de `n_sig` por contexto se estima con 7 grados de libertad, y en los
   brazos donde `n_sig` es 0 en las 8 el `sd` sale exactamente 0 -- que es el
   valor de la muestra, no una afirmacion de que la varianza sea nula. El
   numero 8 y las semillas se fijaron **antes** de la primera corrida (las da
   el encargo) y no se cambiaron al ver ningun valor.
5. 🔴 **Ni las 8 particiones ni los 6 contextos son muestras independientes**
   **del nulo.** Las 8 particiones de un contexto reparten LAS MISMAS 46 guias,
   asi que comparten toda la estructura de guia; y la TABLA 0.2 verifico que los
   46 `ntc_id` son **el mismo conjunto en los seis contextos**. La recurrencia de
   la seccion 1.3 lo hace visible: `CST3` sale llamado en 9 de las 48 corridas
   del brazo `completo`, repartidas por A, C, E y F. Asi que el `sd` de `n_sig`
   sobre particiones **no es** el error de un muestreo independiente, y los 48
   valores de un brazo **no son** 48 replicas. Es la limitacion mas fuerte de
   todo el reporte y no tiene arreglo con estos datos: el reto entrega esas 46
   guias y no mas.
6. **23 contra 23 es una sola forma de partir.** No se probaron mitades
   desiguales (por ejemplo 5 guias contra 41), que es la forma que mas se
   parece a comparar una perturbacion contra el pool. Queda sin medir.
7. 🔴 **El submuestreo a 400 reparte las 400 celulas entre las 23 guias de la**
   **mitad, no las concentra.** 400 celulas de 23 guias son ~17 por guia. Una
   perturbacion del reto lleva sus 400 celulas de **un solo** constructo. Asi
   que el brazo de `n = 400` promedia el efecto de guia sobre 23 guias donde la
   entrega lo tendria sobre uno. Eso **subestima** la componente de guia del
   nulo a `n = 400`, y la direccion del sesgo esta declarada pero no cuantificada.
8. **El piso de deteccion de la seccion 2.2 es, en las bandas estrechas, una**
   **extrapolacion.** Donde `filas_esperadas_bajo_p_crit << 1` la cifra sale de
   un ajuste log-lineal sobre los tramos del propio nulo, no de una lectura
   directa; va marcada `(EXTRAPOLADO)` con su residuo en dex. La forma
   log-lineal **es** una suposicion importada; el nulo solo da los puntos.
9. **El `|log2FC|` del nulo tiene cola gruesa por ceros, no por efecto.** Con
   `epsilon = 1e-9` y 400 celulas, un gen con CPM medio de referencia por
   encima de 5 pero cero cuentas en el lado blanco da
   `|log2FC|` de hasta **33.38** sin que su `p`
   se mueva. Las tablas lo reportan como techo de ruido, pero significa que
   cualquier resumen del nulo en unidades de `|log2FC|` esta dominado por esa
   cola y no por la escala tipica.
10. **La particion que agrupa guias usa PC1, que es UNA eleccion.** PC1 maximiza
   la varianza entre guias en el espacio `log1p(CPM)` del panel probable, pero
   no es el corte que maximiza `n_sig`. Un corte elegido mirando `n_sig`
   seria circular y no se hizo. Asi que la TABLA 2.4 es **una** particion
   adversaria, no el peor caso.
   PC1 explica 6.64-15.96 % de la
   varianza entre guias, asi que el resto de la estructura de guia queda fuera
   del corte.
11. **El ANOVA de la seccion 2.3 mide la MEDIA, la puerta mide RANGOS.** `F` e
    `ICC` se calculan sobre `log1p(CPM)` por celula; el Wilcoxon del preset
    convierte eso a rangos antes de decidir. Las dos cantidades miran la misma
    estructura de guia, pero no son la misma estadistica, y un `ICC` pequeno no
    se traduce en un `n_sig` por una formula que este reporte mida.
12. **La Parte 0c probo 3 blancos por emision, no 300.** La cota de `nnz` es
    estructural (el brazo Poisson no puede crear no-ceros) y por eso 6 brazos
    alcanzan para verificarla; pero la MAGNITUD del `delta nnz` y de la suma de
    cuentas se midio en 6 pares (emision, blanco) y una sola muestra de 400
    celulas de un solo contexto (`D`, semilla `20274999`).
13. **La cobertura de la Parte 0b es de SIMBOLOS, no de identidad de gen.** Se
    cruzaron nombres de gen con nombres de gen. Un simbolo que cambio de
    nomenclatura entre gwps y el panel del reto cuenta como ausente aunque sea
    el mismo gen. No se uso ningun mapa de alias, y cuantos casos hay no se midio.

## Lo que no se midio

* La permutacion inversa de cada particion (intercambiar referencia y blanco).
* Mitades desiguales de guias, ni particiones de mas de dos grupos.
* El mismo nulo con `gpudge` en vez de `pdex`, ni cuanto difieren los dos motores.
* El nulo con el submuestreo a 400 **concentrado en una sola guia**, que es la
  forma de la entrega (una perturbacion = un constructo).
* Las otras cinco metricas del preset: este reporte solo toca la puerta de
  significancia (el conjunto DE). `pds_cosine`, `expr_mse_unbiased_capped_norm`,
  `de_wilcoxon_lfc_nmae`, `de_wilcoxon_direction_fidelity_yield_raw` y
  `de_wilcoxon_direction_reach_raw` no se corrieron aqui.
* El lado perturbado, en 2026 (no existe) y en 2025/Jiang (fuera del encargo).
* Ningun `vcc prep`, ni con `--dry-run`: no hizo falta para nada de esto.

## Las entradas

| que | ruta |
|---|---|
| controles del reto, `val` | `/datos/reto/val/context_{A,B,C}.h5ad` |
| controles del reto, `final` | `/datos/reto/final/context_{D,E,F}.h5ad` |
| paneles y manifiestos | `/datos/reto/{val,final}/{gene_names,pert_counts}.csv`, `manifest.json` |
| gwps | `/datos/ref/replogle/k562_gwps_raw_bulk.h5ad` |
| Jiang IFNB (solo celulas `NT`) | `/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad` |
| la puerta | `cell_eval2/configs/vcc2026.yaml`, `de_compute.py`, `de.py` |
| interprete | `/home/ubuntu/mf/envs/vcc/bin/python -I` (el del marcador) |

| salida | ruta |
|---|---|
| Parte 0a y 0c | `modelo/out/ntc_nulo/p0.json` |
| Parte 0b (corregida) | `modelo/out/ntc_nulo/p0b.json` |
| Parte 1, las corridas crudas | `modelo/out/ntc_nulo/p1.json` |
| Parte 2a y 2b | `modelo/out/ntc_nulo/p2ab.json` |
| Parte 2c | `modelo/out/ntc_nulo/p2c.json` |
| Parte 3 | `modelo/out/ntc_nulo/p3.json` |
| verificacion de la puerta | `modelo/out/ntc_nulo/verif_puerta.log` |
| vectores crudos por corrida (`p`, `log2FC`, CPM ref) | `modelo/out/ntc_nulo/crudo_{ctx}_{j}_{brazo}.npz` |
