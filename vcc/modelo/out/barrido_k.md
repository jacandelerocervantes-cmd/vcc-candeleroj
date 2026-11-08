# El barrido de k sobre VCC 2025, blancos retenidos
## 🔴 CÓMO SE LEE
Este documento mide **la sensibilidad del marcador al factor de des-encogimiento k** sobre una sub-superficie de VCC 2025 con verdad, en blancos que el modelo NO vio. No es una prediccion del score real de la entrega: la linea es otra, el panel es otro (18,080 genes de 2025 contra 18,533 del reto) y el control es otro. Lo unico que podria cruzar de aqui a una reemision es **k, un numero**.
🔴 **No se ha enviado nada** y no se ha tocado D/E/F: alli no hay verdad y no se puede puntuar nada. El `.vcc` de respaldo esta intacto (seccion 0.3).
🔴 **Sin veredicto.** No se elige k aqui. Las tablas se presentan con sus controles y sus defectos.
_Brazos cerrados: 6 de 6 — k=1.0, k=3.4, 8 rasgos, k=1.0, k=2.0, k=1.5, k=2.5._ 
---

## 0. Verificaciones previas
### 0.1 La CLI al dia
|  | version |
|---|---|
| antes de `uv tool upgrade vcc-cli` | `vcc 0.3.0` |
| despues | `vcc 0.3.0` — la misma |

La version no cambio; el arbol de dependencias si: pandas 3.1.0→3.1.1, pydantic 2.14.1→2.14.3, zarr 3.4.1→3.4.2, scverse-misc 0.1.7→0.1.8, tzdata 2026.5→2026.6, python-dateutil 3.0.0→4.0.0 (mayor), six retirado. Por el salto mayor de `python-dateutil` se revalido el `.h5ad` ya empaquetado con la CLI nueva: ✓ dry run; 360,000 celulas; 18,533 genes; D/E/F 120,000 cada uno; targets: verified against the official list
### 0.2 El paquete es el FINAL, no el de validacion
| campo | valor |
|---|---|
| `partition` | `test` |
| `panel_id` | `vcc2026-test-1` |
| `contextos` | `D, E, F` |
| `pert_counts.csv` | `301 lineas (cabecera + 300)` |
| `gene_names.csv` | `18,534 lineas (cabecera + 18,533)` |
| `val_hermano` | `/datos/reto/val/ con "partition": "val", "panel_id": "vcc2026-val-1", contextos A/B/C` |

### 0.3 El respaldo, intocable
|  |  |
|---|---|
| ruta | `/datos/modelo/out/entrega/prediction.vcc` |
| bytes | **3,186,462,720** |
| sha256 | `66149bf6794c9f5d8cc062599eff9856451a62456e5a3f5ad710749656c5bff1` |
| mtime | 2026-11-05 19:12:05 UTC |

`vcc prep` corrio **en firme** (`-f`, sin `--dry-run`): el log del paso 5 dice `✓ prep complete` con `output: …/prediction.vcc`, frente al `✓ dry run — no .vcc written` del paso 4. El sha256 se verifico DESPUES de la revalidacion de 0.1: identico. No esta bajo `~/vcc/out` (vive en `/datos`, otro sistema de archivos, y no es enlace). `VCC_TOKEN` no esta definido en el entorno ni su nombre aparece en ningun archivo de `~/vcc` ni en los dotfiles. Sin commits: `~/vcc` no es repo git.
---

## 1. El banco, declarado antes de correr
barrido de k sobre VCC 2025, blancos retenidos del pliegue 0 de la MISMA particion por blanco de la CV; el modelo se reajusta sin ellos. No se toca D/E/F. Lo unico que viaja es k.
| que | valor |
|---|---|
| donde | VCC 2025, pliegue 0 de la MISMA particion por blanco de la CV |
| particion | leida de `entrega_ajuste.json` (semilla 20284950), no reconstruida |
| blancos retenidos | 60 |
| reajuste | sobre los otros 4 pliegues: ningun blanco puntuado estuvo en su ajuste |
| alpha | 1000.0 |
| celulas por etiqueta | 400 |
| control del banco | 9201 de 114528 disponibles en 2025 |
| regimen de DE | 400 contra 9,200 — el que `ntc_nulo` declaro como el de entrega |
| forma del lado real | [32363, 18080] |
| nnz del lado real | 282,452,518 |
| recorte en log2 | ±2.0 aplicado DESPUES de multiplicar por k |
| sin modificar | `suelo_arnes`, `cand01_covar._reescala_entero` (via `sup.aplica`), marcador oficial preset/profile `vcc2026` |
| orden declarado | k=1.0, k=3.4, 8 rasgos, k=1.0, k=2.0, k=1.5, k=2.5 |
| por que ese orden | extremos y ablacion primero, para que una tabla parcial abarque el rango |

**El bundle reproduce el `config_hash` de competicion**: `a5519017556705abf9a6e553ebb4e517b9612e9546331ace7c5ec6c5c30a1a13` = referencia. Las seis estan medidas con las reglas oficiales.
### 1.1 Semillas: familia 20284000+, comprobadas EN CODIGO y guardadas enteras
| uso | semilla |
|---|---|
| `muestra_real` | 20284600 |
| `control` | 20284601 |
| `pred_k1.0` | 20284700 |
| `pred_k3.4` | 20284701 |
| `pred_r8_k1.0` | 20284702 |
| `pred_k2.0` | 20284703 |
| `pred_k1.5` | 20284704 |
| `pred_k2.5` | 20284705 |

8 semillas nuevas contra **19,028** enumeradas (rangos vetados mas las de la entrega): **choques ninguno**.
### 1.2 Lo que el marcador avisa de ESTE lado real
- `de_lfc_nmae: omitted 1 perturbation(s) for an empty gate (e.g. PCCB). The gate is real-side only, so this set is the same for every submission.`
- `de_lfc_nmae: omitted 11 perturbation(s) for fewer than 10 gated genes (e.g. CD276). The gate is real-side only, so this set is the same for every submission.`
- `pred DE table omits 130 of the 76015 reference-significant (target, feature) pairs (0.171%), across 30 of 60 target(s); worst is 'DNMT1' with 12.` El marcador anade que el efecto NO es de un solo signo (su issue #291).

Escala del lado real: **76,015 pares (blanco, gen) significativos sobre 60 blancos**, unos **1,267 por blanco** segun el aviso; la media medida en la seccion 6 es 1,950 porque alli la puerta es `p_adj` sin el recorte de `N_conf`. Son dos conteos distintos del mismo lado real, no una contradiccion.
---

## 2. 1-bis — los rasgos de contexto y la ablacion de 8 rasgos
🔴 **Correccion de premisa, verificada en el codigo.** `profundidad_mediana` y `n_control` **no estan entre los 10 rasgos**. La lista es:

`gwps_perfil`, `gwps_rango`, `gwps_signo`, `tiene_gwps_par`, `log10_cpm_h`, `det_h`, `log_fano_h`, `pasa_filtro_h`, `g_es_TF`, `log10_cpm_g`

profundidad_mediana y n_control NO entran en matriz_rasgos; el modelo empaquetado ya es el que no los usa. Dentro de UNA linea ambos son escalares unicos: su sd sobre las filas del ajuste es 0 exactamente, asi que no se pueden estandarizar ni identificar. Es decir: **el modelo empaquetado ya es el que no los usa**, y quitarlos de 10 daria los mismos 10.
| cantidad | valor |
|---|---|
| sd de `profundidad_mediana` sobre las filas del ajuste | 0.0 (escalar unico) |
| sd de `n_control` sobre las filas del ajuste | 0.0 (escalar unico) |
| n de control de 2025 | 114,528 |

La ablacion que si se puede hacer, **criterio declarado antes de ajustar**: los 10 menos los dos rasgos del BLANCO g.

|  | rasgos |
|---|---|
| fuera (2) | `g_es_TF`, `log10_cpm_g` |
| dentro (8) | `gwps_perfil`, `gwps_rango`, `gwps_signo`, `tiene_gwps_par`, `log10_cpm_h`, `det_h`, `log_fano_h`, `pasa_filtro_h` |

### 2.1 La CV partida por blanco, los dos modelos
**Pearson por blanco retenido** (los 300 blancos, cada uno predicho por un modelo que no lo vio):

| cantidad | 10 rasgos (empaquetado) | 8 rasgos |
|---|---|---|
| mediana | +0.0507 | +0.0505 |
| media | +0.0628 | +0.0627 |
| p10_p90 | [+0.0012, +0.1455] | [+0.0016, +0.1453] |
| min_max | [-0.0972, +0.3103] | [-0.0972, +0.3087] |
| frac_mayor_0 | +0.9000 | +0.9000 |
| frac_mayor_0.1 | +0.2200 | +0.2200 |
| n | 300 | 300 |

**R² por blanco retenido** (los 300 blancos, cada uno predicho por un modelo que no lo vio):

| cantidad | 10 rasgos (empaquetado) | 8 rasgos |
|---|---|---|
| mediana | -0.0092 | -0.0091 |
| media | -0.0152 | -0.0151 |
| p10_p90 | [-0.0478, +0.0099] | [-0.0463, +0.0097] |
| min_max | [-0.1473, +0.0224] | [-0.1496, +0.0238] |
| frac_mayor_0 | +0.3067 | +0.2933 |
| frac_mayor_0.1 | +0.0000 | +0.0000 |
| n | 300 | 300 |

🔴 **El R² por blanco es negativo en los dos modelos**, con el Pearson positivo: la funcion acierta algo de la FORMA del perfil y falla la ESCALA. Es la cifra que motiva este barrido. En el informe de la entrega el R² publicado era el global (+0.0052), no el de cada blanco.
### 2.2 Intermedio crudo: 20 blancos retenidos, los DOS modelos lado a lado
| # | blanco | Pearson 10r | Pearson 8r | R² 10r | R² 8r |
|---|---|---|---|---|---|
| 1 | `SH3BP4` | -0.0972 | -0.0972 | -0.0708 | -0.0726 |
| 2 | `MAT2A` | -0.0154 | -0.0153 | -0.0726 | -0.0757 |
| 3 | `PBX3` | +0.0016 | +0.0019 | -0.0405 | -0.0367 |
| 4 | `MKI67` | +0.0086 | +0.0085 | -0.0372 | -0.0362 |
| 5 | `MPC2` | +0.0178 | +0.0175 | -0.0388 | -0.0419 |
| 6 | `RRM2` | +0.0241 | +0.0237 | -0.0258 | -0.0247 |
| 7 | `HAT1` | +0.0284 | +0.0281 | -0.0191 | -0.0190 |
| 8 | `NDUFB6` | +0.0343 | +0.0342 | -0.0167 | -0.0195 |
| 9 | `LAP3` | +0.0405 | +0.0407 | -0.0199 | -0.0208 |
| 10 | `DHCR24` | +0.0472 | +0.0467 | -0.0128 | -0.0157 |
| 11 | `CENPO` | +0.0534 | +0.0537 | -0.0093 | -0.0085 |
| 12 | `KRT18` | +0.0604 | +0.0602 | -0.0173 | -0.0156 |
| 13 | `PRR12` | +0.0671 | +0.0669 | -0.0039 | -0.0027 |
| 14 | `PHF14` | +0.0787 | +0.0775 | -0.0021 | -0.0003 |
| 15 | `CTSV` | +0.0896 | +0.0890 | +0.0058 | +0.0054 |
| 16 | `SMARCA5` | +0.1022 | +0.1026 | +0.0073 | +0.0073 |
| 17 | `MED12` | +0.1161 | +0.1164 | +0.0070 | +0.0041 |
| 18 | `SUPV3L1` | +0.1402 | +0.1399 | +0.0135 | +0.0135 |
| 19 | `DNAJA3` | +0.1853 | +0.1846 | +0.0149 | +0.0148 |
| 20 | `KDM1A` | +0.3103 | +0.3087 | +0.0082 | +0.0077 |

_De los 300 blancos retenidos, equiespaciados por su Pearson del modelo de 10 rasgos. Cada blanco predicho por un modelo que no lo vio._

Diferencia |Pearson 10r − Pearson 8r| por blanco: mediana 0.00025, maxima 0.00213, sobre 300 blancos. Es el tamano del efecto de quitar los dos rasgos del blanco g, blanco a blanco.
---

## 3. Las seis metricas por brazo — `u` (fila `mean` de `agg_results.csv`)
| metrica | k=1.0 | k=3.4 | 8 rasgos, k=1.0 | k=2.0 | k=1.5 | k=2.5 |
|---|---|---|---|---|---|---|
| `pds_cosine` | 0.6305 | 0.6542 | 0.6455 | 0.6280 | 0.6271 | 0.6229 |
| `expr_mse_unbiased_capped_norm` | 0.9922 | 1.0509 | 0.9959 | 1.0076 | 0.9964 | 1.0222 |
| `de_wilcoxon_lfc_nmae` | 0.9909 | 0.9761 | 0.9915 | 0.9860 | 0.9895 | 0.9842 |
| `de_wilcoxon_direction_fidelity_yield_raw` | 0.0023 | 0.2461 | 0.0002 | 0.0631 | 0.0218 | 0.1216 |
| `de_wilcoxon_direction_reach_raw` | 0.0398 | 0.0509 | 0.0448 | 0.0541 | 0.0658 | 0.0819 |
| `de_wilcoxon_sig_jaccard` | 0.0333 | 0.0272 | 0.0334 | 0.0183 | 0.0181 | 0.0073 |

## 4. `s` — normalizado contra el bundle real
**`from_replicate`:**

| metrica | k=1.0 | k=3.4 | 8 rasgos, k=1.0 | k=2.0 | k=1.5 | k=2.5 |
|---|---|---|---|---|---|---|
| `pds_cosine` | +0.2793 | +0.3294 | +0.3109 | +0.2739 | +0.2721 | +0.2632 |
| `expr_mse_unbiased_capped_norm` | +0.0000 | +0.0000 | +0.0000 | +0.0000 | +0.0000 | +0.0000 |
| `de_wilcoxon_lfc_nmae` | -0.1025 | -0.0770 | -0.1035 | -0.0940 | -0.1001 | -0.0909 |
| `de_wilcoxon_direction_fidelity_yield_raw` | -0.2035 | +0.1327 | -0.2064 | -0.1197 | -0.1765 | -0.0389 |
| `de_wilcoxon_direction_reach_raw` | -0.0802 | -0.0677 | -0.0746 | -0.0641 | -0.0508 | -0.0326 |
| `de_wilcoxon_sig_jaccard` | +0.0598 | +0.0454 | +0.0599 | +0.0243 | +0.0239 | -0.0015 |
| `avg_score` | -0.0079 | +0.0605 | -0.0023 | +0.0034 | -0.0052 | +0.0165 |

**`from_baseline`:**

| metrica | k=1.0 | k=3.4 | 8 rasgos, k=1.0 | k=2.0 | k=1.5 | k=2.5 |
|---|---|---|---|---|---|---|
| `pds_cosine` | +0.2635 | +0.3108 | +0.2934 | +0.2584 | +0.2568 | +0.2483 |
| `expr_mse_unbiased_capped_norm` | +0.0000 | +0.0000 | +0.0000 | +0.0000 | +0.0000 | +0.0000 |
| `de_wilcoxon_lfc_nmae` | -0.0640 | -0.0481 | -0.0646 | -0.0587 | -0.0625 | -0.0568 |
| `de_wilcoxon_direction_fidelity_yield_raw` | -0.1735 | +0.1132 | -0.1761 | -0.1021 | -0.1506 | -0.0332 |
| `de_wilcoxon_direction_reach_raw` | -0.0796 | -0.0672 | -0.0740 | -0.0636 | -0.0504 | -0.0323 |
| `de_wilcoxon_sig_jaccard` | +0.0256 | +0.0194 | +0.0256 | +0.0104 | +0.0102 | -0.0006 |
| `avg_score` | n/d | n/d | n/d | n/d | n/d | n/d |

## 5. La emision de cada brazo
|  | k=1.0 | k=3.4 | 8 rasgos, k=1.0 | k=2.0 | k=1.5 | k=2.5 |
|---|---|---|---|---|---|---|
| modelo | 10_rasgos | 10_rasgos | 8_rasgos | 10_rasgos | 10_rasgos | 10_rasgos |
| k | 1.0 | 3.4 | 1.0 | 2.0 | 1.5 | 2.5 |
| sd de la razon en log2 | 0.0538 | 0.1828 | 0.0539 | 0.1075 | 0.0806 | 0.1344 |
| sd de la razon | 0.0391 | 0.1503 | 0.0392 | 0.0822 | 0.0602 | 0.1055 |
| razon [min, max] | [0.923, 1.166] | [0.763, 1.686] | [0.925, 1.163] | [0.853, 1.360] | [0.887, 1.259] | [0.819, 1.468] |
| recortados a ±2 | 0 | 0 | 0 | 0 | 0 | 0 |
| semilla | 20284700 | 20284701 | 20284702 | 20284703 | 20284704 | 20284705 |
| segundos | 371.7 | 372.5 | 373.7 | 376.3 | 429.5 | 453.4 |

_Referencia: la sd de la `y` observada en 2025 es 0.1772; la de la entrega empaquetada, 0.0138 en muestra y 0.052 por contexto._
---

## 6. 🔴 La tabla de cuatro celdas, nunca en una sola direccion
Puerta: `p_adj` de los `de_real.parquet` / `de_pred.parquet` que escribe `run --write-degenes`, o sea la del marcador, con el gen propio fuera.

| media por blanco | k=1.0 | k=3.4 | 8 rasgos, k=1.0 | k=2.0 | k=1.5 | k=2.5 |
|---|---|---|---|---|---|---|
| real ∩ pred | 0.033 | 78.133 | 0.117 | 4.150 | 0.633 | 17.117 |
| solo real | 1,950.133 | 1,872.033 | 1,950.050 | 1,946.017 | 1,949.533 | 1,933.050 |
| solo pred (falsos positivos) | 0.100 | 188.067 | 0.083 | 7.650 | 1.267 | 33.617 |
| ambos vacios (de 60 blancos) | 2.000 | 0.000 | 2.000 | 1.000 | 1.000 | 0.000 |
| genes llamados, pred | 0.133 | 266.200 | 0.200 | 11.800 | 1.900 | 50.733 |
| genes llamados, real | 1,950.167 | 1,950.167 | 1,950.167 | 1,950.167 | 1,950.167 | 1,950.167 |

|  | k=1.0 | k=3.4 | 8 rasgos, k=1.0 | k=2.0 | k=1.5 | k=2.5 |
|---|---|---|---|---|---|---|
| precision = ∩ / llamados por pred | 25.0% | 29.4% | 58.3% | 35.2% | 33.3% | 33.7% |
| cobertura = ∩ / llamados por real | 0.0% | 4.0% | 0.0% | 0.2% | 0.0% | 0.9% |

**Y la descomposicion del Jaccard**, para no confundir solape con vacio:

|  | k=1.0 | k=3.4 | 8 rasgos, k=1.0 | k=2.0 | k=1.5 | k=2.5 |
|---|---|---|---|---|---|---|
| `de_wilcoxon_sig_jaccard` medido | 0.0333 | 0.0272 | 0.0334 | 0.0183 | 0.0181 | 0.0073 |
| blancos con ambos conjuntos vacios | 2/60 | 0/60 | 2/60 | 1/60 | 1/60 | 0/60 |
| lo que aportan esos vacios (1.0 cada uno) | 0.0333 | 0.0000 | 0.0333 | 0.0167 | 0.0167 | 0.0000 |

### 6.1 k=1.0 — 20 blancos crudos, de mayor a menor interseccion
| # | blanco | n real | n pred | ∩ | solo real | solo pred |
|---|---|---|---|---|---|---|
| 1 | `TFAM` | 2,811 | 1 | 1 | 2,810 | 0 |
| 2 | `ADNP` | 4,740 | 0 | 0 | 4,740 | 0 |
| 3 | `ATF4` | 1,584 | 0 | 0 | 1,584 | 0 |
| 4 | `CLDN6` | 24 | 1 | 0 | 24 | 1 |
| 5 | `CTNNB1` | 3,718 | 0 | 0 | 3,718 | 0 |
| 6 | `DOT1L` | 5,534 | 0 | 0 | 5,534 | 0 |
| 7 | `EID2` | 13 | 0 | 0 | 13 | 0 |
| 8 | `HDAC8` | 273 | 0 | 0 | 273 | 0 |
| 9 | `HSP90B1` | 3,683 | 0 | 0 | 3,683 | 0 |
| 10 | `LZTR1` | 3,141 | 0 | 0 | 3,141 | 0 |
| 11 | `METTL3` | 8,449 | 0 | 0 | 8,449 | 0 |
| 12 | `MTFR1` | 7 | 0 | 0 | 7 | 0 |
| 13 | `OXCT1` | 2,499 | 0 | 0 | 2,499 | 0 |
| 14 | `SDC1` | 698 | 0 | 0 | 698 | 0 |
| 15 | `SLC25A3` | 957 | 0 | 0 | 957 | 0 |
| 16 | `STAT6` | 1,062 | 1 | 0 | 1,062 | 1 |
| 17 | `TARBP2` | 173 | 0 | 0 | 173 | 0 |
| 18 | `TMSB4X` | 3 | 0 | 0 | 3 | 0 |
| 19 | `WFS1` | 6 | 0 | 0 | 6 | 0 |
| 20 | `ZNF92` | 15 | 0 | 0 | 15 | 0 |

### 6.2 k=3.4 — 20 blancos crudos, de mayor a menor interseccion
| # | blanco | n real | n pred | ∩ | solo real | solo pred |
|---|---|---|---|---|---|---|
| 1 | `METTL3` | 8,449 | 469 | 415 | 8,034 | 54 |
| 2 | `ADNP` | 4,740 | 414 | 216 | 4,524 | 198 |
| 3 | `TAF13` | 4,474 | 355 | 199 | 4,275 | 156 |
| 4 | `NGRN` | 1,766 | 462 | 149 | 1,617 | 313 |
| 5 | `MED15` | 2,038 | 498 | 145 | 1,893 | 353 |
| 6 | `SOX4` | 2,981 | 309 | 125 | 2,856 | 184 |
| 7 | `TFAM` | 2,811 | 273 | 122 | 2,689 | 151 |
| 8 | `ATF4` | 1,584 | 372 | 108 | 1,476 | 264 |
| 9 | `HSBP1` | 1,646 | 313 | 73 | 1,573 | 240 |
| 10 | `SDC1` | 698 | 284 | 56 | 642 | 228 |
| 11 | `COX4I1` | 610 | 287 | 42 | 568 | 245 |
| 12 | `TCF7L2` | 790 | 215 | 23 | 767 | 192 |
| 13 | `TARBP2` | 173 | 243 | 13 | 160 | 230 |
| 14 | `SIX4` | 268 | 312 | 8 | 260 | 304 |
| 15 | `NREP` | 2,707 | 15 | 6 | 2,701 | 9 |
| 16 | `CD276` | 26 | 237 | 1 | 25 | 236 |
| 17 | `TMSB4X` | 3 | 320 | 1 | 2 | 319 |
| 18 | `DZIP3` | 3 | 7 | 0 | 3 | 7 |
| 19 | `MTFR1` | 7 | 267 | 0 | 7 | 267 |
| 20 | `ZNF286A` | 3 | 219 | 0 | 3 | 219 |

### 6.3 8 rasgos, k=1.0 — 20 blancos crudos, de mayor a menor interseccion
| # | blanco | n real | n pred | ∩ | solo real | solo pred |
|---|---|---|---|---|---|---|
| 1 | `DNMT1` | 6,665 | 6 | 5 | 6,660 | 1 |
| 2 | `ACAT2` | 3,139 | 0 | 0 | 3,139 | 0 |
| 3 | `ARID1A` | 2,636 | 0 | 0 | 2,636 | 0 |
| 4 | `CLDN6` | 24 | 0 | 0 | 24 | 0 |
| 5 | `CTNNB1` | 3,718 | 0 | 0 | 3,718 | 0 |
| 6 | `DPH2` | 5,853 | 0 | 0 | 5,853 | 0 |
| 7 | `EPHX1` | 40 | 0 | 0 | 40 | 0 |
| 8 | `HIRA` | 200 | 0 | 0 | 200 | 0 |
| 9 | `KDR` | 11 | 0 | 0 | 11 | 0 |
| 10 | `MBTPS1` | 4,454 | 0 | 0 | 4,454 | 0 |
| 11 | `MTA1` | 683 | 0 | 0 | 683 | 0 |
| 12 | `NREP` | 2,707 | 0 | 0 | 2,707 | 0 |
| 13 | `PLCB3` | 4 | 0 | 0 | 4 | 0 |
| 14 | `SIX4` | 268 | 1 | 0 | 268 | 1 |
| 15 | `SOX4` | 2,981 | 0 | 0 | 2,981 | 0 |
| 16 | `TAF13` | 4,474 | 0 | 0 | 4,474 | 0 |
| 17 | `TET1` | 440 | 0 | 0 | 440 | 0 |
| 18 | `TSFM` | 2,521 | 0 | 0 | 2,521 | 0 |
| 19 | `WFS1` | 6 | 0 | 0 | 6 | 0 |
| 20 | `ZNF92` | 15 | 0 | 0 | 15 | 0 |

### 6.4 k=2.0 — 20 blancos crudos, de mayor a menor interseccion
| # | blanco | n real | n pred | ∩ | solo real | solo pred |
|---|---|---|---|---|---|---|
| 1 | `STRAP` | 6,562 | 57 | 44 | 6,518 | 13 |
| 2 | `DPH2` | 5,853 | 26 | 20 | 5,833 | 6 |
| 3 | `DOT1L` | 5,534 | 12 | 10 | 5,524 | 2 |
| 4 | `SIN3B` | 7,059 | 9 | 8 | 7,051 | 1 |
| 5 | `MLLT1` | 1,064 | 35 | 7 | 1,057 | 28 |
| 6 | `STAT6` | 1,062 | 42 | 5 | 1,057 | 37 |
| 7 | `COX4I1` | 610 | 8 | 4 | 606 | 4 |
| 8 | `MBTPS1` | 4,454 | 3 | 3 | 4,451 | 0 |
| 9 | `C1QBP` | 2,246 | 10 | 2 | 2,244 | 8 |
| 10 | `USF2` | 6,412 | 3 | 2 | 6,410 | 1 |
| 11 | `NREP` | 2,707 | 4 | 1 | 2,706 | 3 |
| 12 | `ZNF593` | 583 | 7 | 1 | 582 | 6 |
| 13 | `CLDN6` | 24 | 3 | 0 | 24 | 3 |
| 14 | `DZIP3` | 3 | 0 | 0 | 3 | 0 |
| 15 | `GSTZ1` | 1 | 1 | 0 | 1 | 1 |
| 16 | `MED25` | 2,907 | 0 | 0 | 2,907 | 0 |
| 17 | `PCCB` | 0 | 3 | 0 | 0 | 3 |
| 18 | `SOX4` | 2,981 | 0 | 0 | 2,981 | 0 |
| 19 | `VCL` | 194 | 3 | 0 | 194 | 3 |
| 20 | `ZNF92` | 15 | 3 | 0 | 15 | 3 |

### 6.5 k=1.5 — 20 blancos crudos, de mayor a menor interseccion
| # | blanco | n real | n pred | ∩ | solo real | solo pred |
|---|---|---|---|---|---|---|
| 1 | `ACAT2` | 3,139 | 22 | 10 | 3,129 | 12 |
| 2 | `HSP90B1` | 3,683 | 3 | 3 | 3,680 | 0 |
| 3 | `STAT6` | 1,062 | 21 | 2 | 1,060 | 19 |
| 4 | `DOT1L` | 5,534 | 1 | 1 | 5,533 | 0 |
| 5 | `MBTPS1` | 4,454 | 4 | 1 | 4,453 | 3 |
| 6 | `TSFM` | 2,521 | 5 | 1 | 2,520 | 4 |
| 7 | `ATF4` | 1,584 | 0 | 0 | 1,584 | 0 |
| 8 | `CLDN6` | 24 | 0 | 0 | 24 | 0 |
| 9 | `DCK` | 1,127 | 2 | 0 | 1,127 | 2 |
| 10 | `EPHX1` | 40 | 0 | 0 | 40 | 0 |
| 11 | `KDR` | 11 | 0 | 0 | 11 | 0 |
| 12 | `MED15` | 2,038 | 3 | 0 | 2,038 | 3 |
| 13 | `MTA1` | 683 | 1 | 0 | 683 | 1 |
| 14 | `OXCT1` | 2,499 | 0 | 0 | 2,499 | 0 |
| 15 | `SIN3B` | 7,059 | 0 | 0 | 7,059 | 0 |
| 16 | `SMARCB1` | 4,627 | 0 | 0 | 4,627 | 0 |
| 17 | `TCF7L2` | 790 | 0 | 0 | 790 | 0 |
| 18 | `TMSB4X` | 3 | 0 | 0 | 3 | 0 |
| 19 | `WFS1` | 6 | 0 | 0 | 6 | 0 |
| 20 | `ZNF92` | 15 | 0 | 0 | 15 | 0 |

### 6.6 k=2.5 — 20 blancos crudos, de mayor a menor interseccion
| # | blanco | n real | n pred | ∩ | solo real | solo pred |
|---|---|---|---|---|---|---|
| 1 | `TAF13` | 4,474 | 168 | 97 | 4,377 | 71 |
| 2 | `DPH2` | 5,853 | 73 | 58 | 5,795 | 15 |
| 3 | `MBTPS1` | 4,454 | 67 | 41 | 4,413 | 26 |
| 4 | `NGRN` | 1,766 | 102 | 38 | 1,728 | 64 |
| 5 | `DOT1L` | 5,534 | 45 | 33 | 5,501 | 12 |
| 6 | `LZTR1` | 3,141 | 46 | 26 | 3,115 | 20 |
| 7 | `ACAT2` | 3,139 | 40 | 20 | 3,119 | 20 |
| 8 | `MED25` | 2,907 | 57 | 19 | 2,888 | 38 |
| 9 | `MLLT1` | 1,064 | 69 | 16 | 1,048 | 53 |
| 10 | `HSP90B1` | 3,683 | 31 | 14 | 3,669 | 17 |
| 11 | `COX4I1` | 610 | 54 | 10 | 600 | 44 |
| 12 | `KLF10` | 319 | 32 | 4 | 315 | 28 |
| 13 | `HDAC8` | 273 | 31 | 2 | 271 | 29 |
| 14 | `TARBP2` | 173 | 10 | 2 | 171 | 8 |
| 15 | `SIX4` | 268 | 29 | 1 | 267 | 28 |
| 16 | `CLDN6` | 24 | 10 | 0 | 24 | 10 |
| 17 | `EID2` | 13 | 34 | 0 | 13 | 34 |
| 18 | `MTFR1` | 7 | 48 | 0 | 7 | 48 |
| 19 | `TMSB4X` | 3 | 73 | 0 | 3 | 73 |
| 20 | `ZNF92` | 15 | 37 | 0 | 15 | 37 |

---

## 7. Defectos propios
1. **Es una sub-superficie de 2025, no la de entrega.** 60 blancos de 300, control submuestreado a 9,200 de 114,528, panel de 18,080 genes. La eleccion del submuestreo fue por coste del bundle (39.6 min medidos), no por biologia.
2. **Un solo pliegue.** Los seis brazos se puntuan sobre los mismos 60 blancos retenidos; no hay repeticion por pliegue, asi que no hay barra de error entre pliegues para las seis metricas, solo la de la CV de la seccion 2.
3. **Una sola semilla de emision por brazo.** La variabilidad de `_reescala_entero` no esta medida: brazos distintos usan semillas distintas, asi que parte de la diferencia entre columnas es ruido de emision y no efecto de k.
4. **12 de los 60 blancos no entran en `lfc_nmae`** (1 por puerta vacia, 11 por menos de 10 genes con puerta), y el conjunto es el mismo para todos los brazos porque la puerta es del lado real.
5. **El aviso del issue #291** afecta a las metricas de direccion de todos los brazos por igual y su efecto no es de un solo signo.
6. **k es un parametro post-hoc.** Multiplicar el delta no cambia la funcion ni su R²; mueve la escala de la emision. Nada de lo medido aqui dice que el k que mejor puntue en 2025 sea el que mejor puntue en D/E/F — son lineas distintas.

_Generado 2026-11-05 22:47:49 UTC desde `barrido_k.json`. Reejecutable: `reporte_bk.py`._
