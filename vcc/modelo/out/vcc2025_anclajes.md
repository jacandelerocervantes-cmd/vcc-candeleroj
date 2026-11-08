# VCC 2025 como superficie: los anclajes nativos (`b`, `r`, `config_hash`)

## 🔴 0. Que es esto y que no

Los dos extremos de la escala `vcc2026` construidos por el camino **NATIVO del marcador** sobre la superficie 2025, con los **mismos armadores de comando** que se usaron en Jiang (`anclas_nativas.cmd_baseline` y `anclas_nativas.cmd_bundle`, sin modificar).

🔴 **`b` y `r` son informacion de ORACULO** (`baseline.py:6-10`): el perfil del baseline se promedia de las perturbaciones reales evaluadas. Son artefactos de PUNTUACION. **Ninguna de sus cifras puede entrar en un modelo.**

## 1. La superficie, declarada antes de medir

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

## 2. Las dos invocaciones y sus guardas

| comando | ok | exitcode | duracion | error (200 car.) |
|---|---|---|---|---|
| `cell-eval2 baseline` | si | 0 | 17.0 min (1020 s) |  |
| `cell-eval2 prep-real-bundle` | si | 0 | 93.2 min (5595 s) |  |

**Comandos literales:**

```
/home/ubuntu/celleval2_venv/bin/cell-eval2 baseline -ar /datos/modelo/out/noche/b4/real_2025_test.h5ad --preset vcc2026 --profile vcc2026 --pert-col target_gene --emit dispersed --seed 0 --save-pred /datos/modelo/out/noche/b4/baseline_pred_2025_test.h5ad -o /datos/modelo/out/noche/b4/BASELINE_2025_test
```

```
/home/ubuntu/celleval2_venv/bin/cell-eval2 prep-real-bundle --real /datos/modelo/out/noche/b4/real_2025_test.h5ad --baseline /datos/modelo/out/noche/b4/baseline_pred_2025_test.h5ad --preset vcc2026 --profile vcc2026 --pert-col target_gene --anchor-base-seed 0 --anchor-splits 5 -o /datos/modelo/out/noche/b4/BUNDLE_2025_test
```

## 3. 🔴 El `config_hash`, verificado contra `competition_payload()`

| #  | cantidad | valor |
|---|---|---|
| 1 | `rule_config_hash` obtenido del `config.yaml` del bundle | `a5519017556705abf9a6e553ebb4e517b9612e9546331ace7c5ec6c5c30a1a13` |
| 2 | `config_hash` de referencia (`competition_payload()`) | `a5519017556705abf9a6e553ebb4e517b9612e9546331ace7c5ec6c5c30a1a13` |
| 3 | **coincide** | **si** |
| 4 | `control` en el config | `non-targeting` |
| 5 | `pert_col` en el config | `target_gene` |
| 6 | bundle diagnostico |  |
| 7 | `anchor_base_seed` / `anchor_splits` | 0 / 5 (convencion del marcador, `competition.py:38,41`) |

## 4. El `b` por metrica: los dos lados

| corta | metrica | **`b` en VCC 2025 (test)** | rango del `b` en las 6 lineas de Jiang | media en Jiang |
|---|---|---|---|---|
| `pds` | pds_cosine | +0.4981 | [+0.4375, +0.5208] | +0.4878 |
| `mse` | expr_mse_unbiased_capped_norm | +0.9623 | [+0.8244, +1.1733] | +1.0035 |
| `lfc` | de_wilcoxon_lfc_nmae | +0.9450 | [+0.8148, +0.8973] | +0.8587 |
| `fid` | de_wilcoxon_direction_fidelity_yield_raw | +0.2106 | [+0.0000, +0.1142] | +0.0321 |
| `reach` | de_wilcoxon_direction_reach_raw | +0.0818 | [+0.1019, +0.4068] | +0.3083 |
| `jac` | de_wilcoxon_sig_jaccard | +0.0071 | [+0.0000, +0.4545] | +0.1972 |

🔴 Los dos lados **no son comparables como numeros**: distinta superficie, distinto panel (18,080 contra 34,025 genes), distinta profundidad (~54,000 contra 3,319-6,203 UMI/celula). La columna de Jiang esta como **control de que la maquinaria corrio igual**, no como comparacion.

## 5. El `r` por metrica (`replicate`, el anclaje de arriba)

| metrica | **`replicate`** | sd | min | max | n perts min | n perts max | estimador |
|---|---|---|---|---|---|---|---|
| `pds` | +0.9762 | +0.0031 | +0.9719 | +0.9814 | 100 | 100 | split_half_raw |
| `mse` | +0.0017 | +0.0037 | -0.0027 | +0.0067 |  |  | split_half_raw |
| `lfc` | +0.3424 | +0.0081 | +0.3357 | +0.3579 | 85 | 85 | full_gate_raw |
| `fid` | +0.8737 | +0.0158 | +0.8532 | +0.8983 | 97 | 98 | split_half_raw |
| `reach` | +0.9875 | +0.0048 | +0.9819 | +0.9937 | 96 | 97 | split_half_raw |
| `jac` | +0.4567 | +0.0092 | +0.4421 | +0.4706 | 100 | 100 | split_half_raw |

### 5.1 🔴 Los intermedios crudos: los 5 splits de anclaje, uno por uno

| #  | metrica | split | semilla | valor | n perturbaciones |
|---|---|---|---|---|---|
| 1 | `pds` | 0 | 2968811710 | +0.9719 | 100 |
| 2 | `pds` | 1 | 3677149159 | +0.9768 | 100 |
| 3 | `pds` | 2 | 745650761 | +0.9754 | 100 |
| 4 | `pds` | 3 | 2884920346 | +0.9814 | 100 |
| 5 | `pds` | 4 | 2642120001 | +0.9757 | 100 |
| 6 | `mse` | 0 | 2968811710 | +0.0024 |  |
| 7 | `mse` | 1 | 3677149159 | -0.0027 |  |
| 8 | `mse` | 2 | 745650761 | +0.0067 |  |
| 9 | `mse` | 3 | 2884920346 | -0.0025 |  |
| 10 | `mse` | 4 | 2642120001 | +0.0043 |  |
| 11 | `lfc` | 0 | 2968811710 | +0.3394 | 85 |
| 12 | `lfc` | 1 | 3677149159 | +0.3363 | 85 |
| 13 | `lfc` | 2 | 745650761 | +0.3357 | 85 |
| 14 | `lfc` | 3 | 2884920346 | +0.3428 | 85 |
| 15 | `lfc` | 4 | 2642120001 | +0.3579 | 85 |
| 16 | `fid` | 0 | 2968811710 | +0.8532 | 97 |
| 17 | `fid` | 1 | 3677149159 | +0.8616 | 97 |
| 18 | `fid` | 2 | 745650761 | +0.8819 | 98 |
| 19 | `fid` | 3 | 2884920346 | +0.8983 | 98 |
| 20 | `fid` | 4 | 2642120001 | +0.8735 | 98 |
| 21 | `reach` | 0 | 2968811710 | +0.9819 | 97 |
| 22 | `reach` | 1 | 3677149159 | +0.9924 | 97 |
| 23 | `reach` | 2 | 745650761 | +0.9867 | 96 |
| 24 | `reach` | 3 | 2884920346 | +0.9937 | 96 |
| 25 | `reach` | 4 | 2642120001 | +0.9829 | 97 |
| 26 | `jac` | 0 | 2968811710 | +0.4589 | 100 |
| 27 | `jac` | 1 | 3677149159 | +0.4706 | 100 |
| 28 | `jac` | 2 | 745650761 | +0.4583 | 100 |
| 29 | `jac` | 3 | 2884920346 | +0.4534 | 100 |
| 30 | `jac` | 4 | 2642120001 | +0.4421 | 100 |

## 6. 🔴 Lo que esta superficie NO puede dar

| #  | lo que no da | por que |
|---|---|---|
| 1 | **Dominancia entre contextos** | Es UN solo contexto. `b` y `r` son de un punto, no de un rango. |
| 2 | **Rango entre contextos** | En Jiang hay 6 `b` y 6 `r` y por eso hay rango. Aqui hay uno. |
| 3 | Un `b` comparable con el de Jiang | Distinto panel y distinta profundidad; la comparacion es de maquinaria, no de valor. |

## 7. Defectos de la construccion

| #  | defecto | que implica |
|---|---|---|
| 1 | `b` y `r` son transductivos. | Artefactos de puntuacion. No entran en ningun modelo. |
| 2 | 15 de 100 perturbaciones llevan < 400 celulas. | La potencia por perturbacion no es uniforme; el `n_perturbations` de cada split de anclaje lo refleja. |
| 3 | El control se submuestreo a 18,400 de 38,176. | Semilla declarada. El `b` depende del pool de control, asi que depende de ese submuestreo. |
| 4 | Una sola realizacion de la poblacion real. | No hay barra de error sobre la eleccion de las celulas por perturbacion. |
| 5 | `--control` se omite a proposito. | Para que el valor quede el del preset y el `config_hash` no se desvie -- el defecto que `anclas_nativas.py` documenta de Jiang. |

*Corrida: 2026-11-05 05:19:38 -> 2026-11-05 07:10:20 UTC, duracion 110.7 min (6642 s).*
