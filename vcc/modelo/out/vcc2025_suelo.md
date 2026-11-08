# VCC 2025 como superficie: el SUELO con 8 semillas y las seis metricas

## 🔴 0. Que es esto y que no

El **suelo** de la superficie 2025: `suelo_floor.predict_floor` sin modificar -- 400 celulas remuestreadas **con reemplazo** del pool de control, por etiqueta, sin tocar ningun conteo. Es la **misma emision** que fijo el suelo de Jiang, sobre la otra superficie.

Ocho semillas, `20280300 + j`, para que el suelo tenga dispersion medida y no un solo punto -- igual que `dispersion_u`.

## 1. La construccion, declarada antes de medir

| #  | cantidad | valor |
|---|---|---|
| 1 | split | `test` (CON verdad) |
| 2 | perturbaciones | 100 |
| 3 | celulas por perturbacion | `min(400, disponibles)` |
| 4 | 🔴 perturbaciones truncadas (< 400 celulas) | **15** de 100 |
| 5 | celulas perturbadas | 37,157 |
| 6 | pool de control `non-targeting` | 38,176 |
| 7 | control usado | 18,400 (el mismo numero que el reto 2026 lleva por contexto) |
| 8 | forma del h5ad real | 55,557 x 18,080 |
| 9 | `nnz` | 485,028,540 (485.0 M) |
| 10 | semillas | `20280201` (celulas/pert), `20280202` (control) -- familia NUEVA 20280000+ |
| 11 | huella sha256 del real | 816d6c0bec7ccdca79beaef6024dc24a... |
| 12 | dtype | el X de 2025 es `float32` con valores enteros; se convierte a `int64`, como el de Jiang |

| #  | cantidad | valor |
|---|---|---|
| 1 | emision | `suelo_floor.predict_floor`, SIN MODIFICAR |
| 2 | celulas emitidas por etiqueta | 400 (con reemplazo del pool de control) |
| 3 | etiquetas emitidas | 101 (100 perturbaciones + el control) |
| 4 | semillas | `20280300 + j`, `j = 0..7` -- familia NUEVA, 8 semillas, 0 choques |
| 5 | puntuacion `u` | `cell-eval2 run --preset vcc2026 --profile vcc2026 --write-degenes`, fila `mean` de `agg_results.csv` |
| 6 | puntuacion `s` | `cell-eval2 score --real-bundle` contra el bundle de 4b |
| 7 | bundle | `/datos/modelo/out/noche/b4/BUNDLE_2025_test` |

## 2. Las seis metricas sobre el suelo 2025

| corta | metrica | n semillas | **`u` media** | `u` sd | **rango de `u`** | **`s` media** | rango de `s` |
|---|---|---|---|---|---|---|---|
| `pds` | pds_cosine | 8 | +0.4913 | +0.0329 | [+0.4407, +0.5414] | -0.0141 | [-0.1200, +0.0906] |
| `mse` | expr_mse_unbiased_capped_norm | 8 | +1.0048 | +0.0022 | [+1.0017, +1.0079] | +0.0000 | [+0.0000, +0.0000] |
| `lfc` | de_wilcoxon_lfc_nmae | 8 | +1.0015 | +0.0026 | [+0.9979, +1.0059] | -0.0938 | [-0.1011, -0.0877] |
| `fid` | de_wilcoxon_direction_fidelity_yield_raw | 8 | +0.0012 | +0.0023 | [+0.0000, +0.0068] | -0.3158 | [-0.3177, -0.3074] |
| `reach` | de_wilcoxon_direction_reach_raw | 8 | +0.0510 | +0.0125 | [+0.0313, +0.0700] | -0.0339 | [-0.0557, -0.0130] |
| `jac` | de_wilcoxon_sig_jaccard | 8 | +0.0288 | +0.0035 | [+0.0200, +0.0300] | +0.0481 | [+0.0287, +0.0509] |

🔴 **Rango, no solo media.** Con 8 semillas el rango es la cantidad que dice cuanta resolucion tiene la superficie.

## 3. 🔴 Los intermedios crudos: las 8 corridas una por una

| j | semilla | huella pred (12) | forma pred | nnz pred | run ok | duracion del run | score ok | `config_hash` coincide | duracion total |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 20280300 | a21ff3bf1899... | 40,400x18,080 | 353.5 M | si | 16.2 min (973 s) | si | True | 16.5 min (991 s) |
| 1 | 20280301 | 21a216a71151... | 40,400x18,080 | 353.3 M | si | 16.2 min (973 s) | si | True | 16.5 min (987 s) |
| 2 | 20280302 | 74c260201286... | 40,400x18,080 | 353.9 M | si | 16.2 min (972 s) | si | True | 16.4 min (986 s) |
| 3 | 20280303 | eb96dadd53d2... | 40,400x18,080 | 353.3 M | si | 16.2 min (971 s) | si | True | 16.4 min (986 s) |
| 4 | 20280304 | 22eac20edff2... | 40,400x18,080 | 353.6 M | si | 16.2 min (971 s) | si | True | 16.4 min (986 s) |
| 5 | 20280305 | 0ec3e8dff36c... | 40,400x18,080 | 353.3 M | si | 16.2 min (972 s) | si | True | 16.4 min (986 s) |
| 6 | 20280306 | fc7970aba0c2... | 40,400x18,080 | 353.5 M | si | 16.2 min (972 s) | si | True | 16.4 min (987 s) |
| 7 | 20280307 | 7c684e4edaba... | 40,400x18,080 | 353.7 M | si | 16.2 min (972 s) | si | True | 16.5 min (987 s) |

### 3.1 El `u` crudo, corrida por corrida

| j | semilla | `pds` | `mse` | `lfc` | `fid` | `reach` | `jac` |
|---|---|---|---|---|---|---|---|
| 0 | 20280300 | +0.4407 | +1.0061 | +1.0015 | +0.0013 | +0.0413 | +0.0300 |
| 1 | 20280301 | +0.4754 | +1.0079 | +1.0001 | +0.0000 | +0.0632 | +0.0300 |
| 2 | 20280302 | +0.4929 | +1.0056 | +0.9999 | +0.0000 | +0.0450 | +0.0300 |
| 3 | 20280303 | +0.5414 | +1.0026 | +0.9979 | +0.0000 | +0.0568 | +0.0300 |
| 4 | 20280304 | +0.5247 | +1.0017 | +1.0045 | +0.0005 | +0.0313 | +0.0300 |
| 5 | 20280305 | +0.4844 | +1.0038 | +1.0013 | +0.0011 | +0.0459 | +0.0200 |
| 6 | 20280306 | +0.4628 | +1.0070 | +1.0010 | +0.0001 | +0.0700 | +0.0300 |
| 7 | 20280307 | +0.5082 | +1.0040 | +1.0059 | +0.0068 | +0.0550 | +0.0300 |

### 3.2 El `s` crudo, corrida por corrida

| j | semilla | `pds` | `mse` | `lfc` | `fid` | `reach` | `jac` | **avg_score** |
|---|---|---|---|---|---|---|---|---|
| 0 | 20280300 | -0.1200 | +0.0000 | -0.0938 | -0.3157 | -0.0447 | +0.0509 | -0.0872 |
| 1 | 20280301 | -0.0475 | +0.0000 | -0.0915 | -0.3177 | -0.0205 | +0.0509 | -0.0710 |
| 2 | 20280302 | -0.0108 | +0.0000 | -0.0911 | -0.3177 | -0.0406 | +0.0509 | -0.0682 |
| 3 | 20280303 | +0.0906 | +0.0000 | -0.0877 | -0.3177 | -0.0276 | +0.0509 | -0.0486 |
| 4 | 20280304 | +0.0558 | +0.0000 | -0.0988 | -0.3168 | -0.0557 | +0.0509 | -0.0608 |
| 5 | 20280305 | -0.0285 | +0.0000 | -0.0935 | -0.3160 | -0.0396 | +0.0287 | -0.0748 |
| 6 | 20280306 | -0.0737 | +0.0000 | -0.0928 | -0.3176 | -0.0130 | +0.0509 | -0.0744 |
| 7 | 20280307 | +0.0211 | +0.0000 | -0.1011 | -0.3074 | -0.0296 | +0.0509 | -0.0610 |

### 3.3 Las cohortes moviles, corrida por corrida

| j | perts puntuadas | `E_union_vacia` | perts con pred significativo | suma `n_DE_pred` | suma `n_DE_real` | **perts con DE REAL vacio** |
|---|---|---|---|---|---|---|
| 0 | 100 | 3 | 6 | 7 | 173473 | 3 |
| 1 | 100 | 3 | 3 | 4 | 173473 | 3 |
| 2 | 100 | 3 | 4 | 4 | 173473 | 3 |
| 3 | 100 | 3 | 1 | 1 | 173473 | 3 |
| 4 | 100 | 3 | 2 | 2 | 173473 | 3 |
| 5 | 100 | 2 | 7 | 8 | 173473 | 3 |
| 6 | 100 | 3 | 6 | 6 | 173473 | 3 |
| 7 | 100 | 3 | 10 | 10 | 173473 | 3 |

🔴 La ultima columna es del lado real, asi que **no varia con la semilla del suelo**: es una propiedad de la superficie 2025, igual que el `control_permutado.md` §7 lo es de Jiang.

## 4. 🔴 Lo que esta superficie NO puede dar

| #  | lo que no da | por que |
|---|---|---|
| 1 | **Dominancia entre contextos** | Es UN solo contexto. Las 8 semillas dan dispersion de la EMISION, no de contexto. |
| 2 | **Rango entre contextos** | El rango de esta tabla es sobre semillas del suelo, no sobre contextos. En Jiang hay 6 lineas x 8 semillas y los dos rangos se distinguen; aqui solo existe el segundo. |
| 3 | Transferencia del suelo a 2026 | Panel y profundidad difieren; el valor del suelo es de esta superficie. |

## 5. Defectos de la construccion

| #  | defecto | que implica |
|---|---|---|
| 1 | El lado real es UNA realizacion (una semilla de celulas por perturbacion). | Todas las 8 corridas comparten el mismo lado real, a proposito: la dispersion medida es solo de la emision. |
| 2 | 15 de 100 perturbaciones llevan < 400 celulas reales. | Su potencia es menor y arrastra las metricas de DE. La tabla por perturbacion esta en el json. |
| 3 | El suelo emite 400 celulas por etiqueta aunque el real tenga menos. | Los dos lados no tienen el mismo `n` en esas 15 perturbaciones. |
| 4 | El `s` depende del bundle de 4b. | Si 4b fallo, el `s` no existe y queda declarado; el `u` no depende del bundle. |
| 5 | El control emitido tambien es control remuestreado. | Igual que en Jiang: el suelo predice el control con el propio control. |

*Corrida: 2026-11-05 07:10:21 -> 2026-11-05 09:22:02 UTC, duracion 131.7 min (7901 s).*
