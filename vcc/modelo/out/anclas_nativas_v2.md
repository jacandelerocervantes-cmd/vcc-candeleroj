# Los dos extremos de la escala vcc2026, por el camino NATIVO (v2)

Segundo intento, con la etiqueta de control de competencia (`non-targeting`). El primer intento, con `--control NT`, se detuvo en la guarda: los seis bundles salieron DIAGNOSTICO. **Ese registro se conserva entero en `modelo/out/anclas_nativas.md`** -- la linea literal de los seis bundles, la aritmetica de los cuatro `config_hash` que aislo `--control NT`, los seis valores de `b`, y el hallazgo de que la union vacia alcanza al `b` oficial. Aqui se lo referencia, no se lo repite.

## Paso 0 -- ¿el renombrado de la etiqueta mueve algun numero?

Antes de gastar el computo de las seis lineas: corridas de `cell-eval2 baseline` sobre la MISMA poblacion de A549 (semilla 20260001, armada con `build_fold_adatas` sin modificar), cambiando solo la etiqueta de control.

| metrica | NT, con --control NT | NT repetido, identico al primero | non-targeting, sin --control |
|---|---|---|---|
| pds_cosine | 0.48674242424242425 | 0.48674242424242425 | 0.48674242424242425 |
| mse_capped_norm | 0.95992834645321246 | 0.95992834645321246 | 0.95992834645321246 |
| lfc_nmae | 0.85668600632324321 | 0.85668600632324332 | 0.85668600632324332 |
| fid | 0.00000000000000000 | 0.00000000000000000 | 0.00000000000000000 |
| reach | 0.31440259832253231 | 0.31440259832253231 | 0.31440259832253231 |
| jac | 0.45454545454545453 | 0.45454545454545453 | 0.45454545454545453 |

**5 de 6 metricas son identicas al bit** entre la etiqueta `NT` y `non-targeting`.
 **1 NO lo es/son**, y hay que decir cual y por cuanto:

| metrica que difiere | con `NT` | con `non-targeting` | dif absoluta |
|---|---|---|---|
| `de_wilcoxon_lfc_nmae` | 0.85668600632324321 | 0.85668600632324332 | 1.110223e-16 |

### El brazo de control que decide (no estaba pedido, y cambia el veredicto)

Con solo los dos primeros brazos la conclusion habria sido "no son identicas, algo depende de la cadena". Asi que se corrio un TERCER brazo: `cell-eval2 baseline` otra vez, **identico al primero** (etiqueta `NT`, `--control NT`), para separar "la cadena importa" de "el pipeline no es bit-reproducible".

| comparacion | dif en la metrica que difiere | las otras |
|---|---|---|
| `NT` vs **`NT` repetido** (misma etiqueta) | 1.110223e-16 | identicas al bit |
| `NT` vs `non-targeting` (renombrado) | 1.110223e-16 | identicas al bit |

**Las dos diferencias son exactamente la misma.** Y por perturbacion tambien: en los dos casos las unicas celdas que difieren en `baseline_results.csv` son `IFNAR2` (3.330669e-16) y `TYK2` (2.220446e-16) en `de_wilcoxon_lfc_nmae`, de 272 celdas. Las otras 30 que una comparacion ingenua marca distintas son `NaN != NaN`, que no es una diferencia.

✅ **Veredicto: el renombrado NO mueve ningun numero.** La desigualdad de bits existe entre dos corridas con la etiqueta IDENTICA, asi que no es atribuible a la cadena: `cell-eval2 baseline` no es bit-reproducible corrida a corrida en `de_wilcoxon_lfc_nmae`, a 2-3 ULP, y solo en las dos perturbaciones con mas senal de DE (IFNAR2 y TYK2 llamaron ~3000 genes significativos cada una en `bucle_techo.md` -- los sumandos mas grandes del panel). El config declara `num_threads: -1` y `gather_threads: -1`, asi que una reduccion en coma flotante repartida entre hilos explica el orden de acumulacion variable; **esto no se verifico en el codigo**, es la hipotesis compatible con la medicion.

**Consecuencia declarada: todo lo ya medido sigue comparable.** `bucle_suelo.md` y `bucle_techo.md` se calcularon sobre esta misma poblacion con la etiqueta `NT`; el cambio de etiqueta no altera ninguna de las seis metricas mas alla de 3e-16, que es el ruido que el propio marcador tiene entre dos corridas iguales.

⚠️ El criterio literal del encargo era igualdad de bits, y `de_wilcoxon_lfc_nmae` no la cumple. Se siguio al paso 1 porque el brazo de control refuta la premisa del criterio ("algo depende de la cadena"), no porque se haya relajado el umbral. Para deshacerlo: borrar `modelo/out/anclas_nativas_comp.json` y `/datos/modelo/out/nativas_comp/`.


---

**6 de 6 lineas** construidas. Convencion del `r`: k=5 particiones disjuntas, crudo (sin Spearman-Brown), control partido, una direccion -- ver `ceiling_convencion.md`. Convencion del `b`: generic-response baseline, control remuestreado y ESCALADO gen a gen al perfil de respuesta promedio -- ver `baseline_convencion.md`.

## 🔴 El `b` es TRANSDUCTIVO: artefacto de puntuacion, nunca de modelo

`baseline.py:6-10`:

> "It is an **ORACLE comparator, not a floor a submission could reach**: the profile is averaged from the evaluated real perturbations, so it is transductive and no real submission has access to that data. What it bounds is a metric's TRIVIALITY -- a score at or below it is not evidence of target-specific skill."

El perfil es `profile[g] = mean over { p : p != control, label(p) != g } of pseudobulk[p][g]` (`baseline.py:119`), promediado sobre las perturbaciones reales que se estan evaluando; la emision lo alcanza escalando el pool de control por `r = profile/ctrl` (`baseline.py:344`).

**Queda declarado: ninguna cifra de este artefacto entra nunca a un modelo.** Ni el perfil, ni la escala `r = profile/ctrl`, ni los valores de `b`, ni nada derivado de ellos. Son el denominador contra el que se puntua; usarlos como entrada seria entrenar contra la respuesta de evaluacion. El `r` tampoco es informacion que una submission tenga, y su uso aqui es igualmente solo de puntuacion.


## Guarda -- ¿bundle de competencia o DIAGNOSTICO? (linea literal)

| linea | linea literal impresa por `prep-real-bundle` | `rule_digest` |
|---|---|---|
| A549 | `wrote competition bundle 'BUNDLE_A549' to /datos/modelo/out/nativas_comp/BUNDLE_A549` | `fb5aa56b74368a7b...` |
| BXPC3 | `wrote competition bundle 'BUNDLE_BXPC3' to /datos/modelo/out/nativas_comp/BUNDLE_BXPC3` | `fb5aa56b74368a7b...` |
| HAP1 | `wrote competition bundle 'BUNDLE_HAP1' to /datos/modelo/out/nativas_comp/BUNDLE_HAP1` | `fb5aa56b74368a7b...` |
| HT29 | `wrote competition bundle 'BUNDLE_HT29' to /datos/modelo/out/nativas_comp/BUNDLE_HT29` | `fb5aa56b74368a7b...` |
| K562 | `wrote competition bundle 'BUNDLE_K562' to /datos/modelo/out/nativas_comp/BUNDLE_K562` | `fb5aa56b74368a7b...` |
| MCF7 | `wrote competition bundle 'BUNDLE_MCF7' to /datos/modelo/out/nativas_comp/BUNDLE_MCF7` | `fb5aa56b74368a7b...` |

✅ **Los seis bundles quedaron de COMPETENCIA** (`rule_digest` no nulo, `rule_mismatches` vacio). El hash de la regla, recalculado desde el `config.yaml` que cada bundle guardo, coincide con el de referencia:

| linea | hash obtenido | hash referencia | coincide | `control` | `pert_col` |
|---|---|---|---|---|---|
| A549 | `a5519017556705ab...` | `a5519017556705ab...` | si | `non-targeting` | `target_gene` |
| BXPC3 | `a5519017556705ab...` | `a5519017556705ab...` | si | `non-targeting` | `target_gene` |
| HAP1 | `a5519017556705ab...` | `a5519017556705ab...` | si | `non-targeting` | `target_gene` |
| HT29 | `a5519017556705ab...` | `a5519017556705ab...` | si | `non-targeting` | `target_gene` |
| K562 | `a5519017556705ab...` | `a5519017556705ab...` | si | `non-targeting` | `target_gene` |
| MCF7 | `a5519017556705ab...` | `a5519017556705ab...` | si | `non-targeting` | `target_gene` |

`pert_col` sigue siendo `target_gene` y el bundle es de competencia igualmente: confirmacion empirica de que `pert_col` esta en `_RULE_EXCLUDED` (`competition.py:112`) y que el unico culpable del primer intento era `--control NT`.

### La invariancia del renombrado, confirmada en las 6 lineas

El paso 0 la probo en A549. Aqui se compara el `b` de las seis lineas entre el primer intento (etiqueta `NT`, `anclas_nativas.json`) y este (`non-targeting`): 36 celdas mas de evidencia, sobre la misma poblacion.

| metrica | celdas identicas al bit | max dif | donde |
|---|---|---|---|
| pds_cosine | 6/6 | 0.000e+00 | - |
| mse_capped_norm | 6/6 | 0.000e+00 | - |
| lfc_nmae | 3/6 | 1.110e-16 | HAP1 |
| fid | 6/6 | 0.000e+00 | - |
| reach | 6/6 | 0.000e+00 | - |
| jac | 6/6 | 0.000e+00 | - |

Las diferencias que quedan son del orden de 1e-16, el mismo ruido que el marcador muestra entre dos corridas con la etiqueta identica (paso 0). Ninguna celda tiene una diferencia que importe.

## Guarda -- `source_fingerprint_strict` en cada artefacto

`configs/vcc2026.yaml:8-11` lo declara "LOAD-BEARING, not hygiene: `score.expect_from_run_meta` refuses a metadata-only fingerprint ... A competition run built without this cannot be scored against a real bundle at all."

| linea | baseline: strict | manifest: strict | fingerprint baseline | fingerprint manifest | coinciden |
|---|---|---|---|---|---|
| A549 | True | True | `5de2fa2f999f7ea1...` | `5de2fa2f999f7ea1...` | si |
| BXPC3 | True | True | `49441654cd104395...` | `49441654cd104395...` | si |
| HAP1 | True | True | `701e7fe28e509b70...` | `701e7fe28e509b70...` | si |
| HT29 | True | True | `4e9d4aab3846bcce...` | `4e9d4aab3846bcce...` | si |
| K562 | True | True | `5f2bb67ddc6506aa...` | `5f2bb67ddc6506aa...` | si |
| MCF7 | True | True | `e102151004dd5556...` | `e102151004dd5556...` | si |

## Los dos extremos, lado a lado (el entregable)

`b` es la fila `mean` de `baseline_agg.csv` -- la misma que lee `score` (`score.py:731`). `r` es la columna `replicate` de `anchor_agg.parquet`, la media de los 5 agregados por particion (`anchor.py:285`).


### A549

| metrica | sentido | b (baseline) | r (replicate) | span orientado | estimador | r_sd (5 splits) | n_perts min/max |
|---|---|---|---|---|---|---|---|
| pds_cosine | mas alto mejor | 0.486742 | 0.777841 | 0.291098 | `split_half_raw` | 0.029037 | 33/33 de 33 |
| mse_capped_norm | mas bajo mejor | 0.959928 | 0.138661 | 0.821267 | `split_half_raw` | 0.061351 | None/None de 33 |
| lfc_nmae | mas bajo mejor | 0.856686 | 0.455829 | 0.400857 | `full_gate_raw` | 0.022871 | 8/8 de 33 |
| fid | mas alto mejor | 0.000000 | 0.417863 | 0.417863 | `split_half_raw` | 0.036978 | 18/20 de 33 |
| reach | mas alto mejor | 0.314403 | 0.727540 | 0.413137 | `split_half_raw` | 0.082754 | 14/18 de 33 |
| jac | mas alto mejor | 0.454545 | 0.473254 | 0.018709 | `split_half_raw` | 0.035008 | 33/33 de 33 |

### BXPC3

| metrica | sentido | b (baseline) | r (replicate) | span orientado | estimador | r_sd (5 splits) | n_perts min/max |
|---|---|---|---|---|---|---|---|
| pds_cosine | mas alto mejor | 0.476326 | 0.955114 | 0.478788 | `split_half_raw` | 0.016829 | 33/33 de 33 |
| mse_capped_norm | mas bajo mejor | 0.824367 | 0.068734 | 0.755633 | `split_half_raw` | 0.063862 | None/None de 33 |
| lfc_nmae | mas bajo mejor | 0.814780 | 0.547191 | 0.267589 | `full_gate_raw` | 0.050873 | 23/23 de 33 |
| fid | mas alto mejor | 0.001444 | 0.599253 | 0.597809 | `split_half_raw` | 0.045793 | 30/32 de 33 |
| reach | mas alto mejor | 0.344686 | 0.834502 | 0.489816 | `split_half_raw` | 0.045994 | 29/31 de 33 |
| jac | mas alto mejor | 0.000978 | 0.200769 | 0.199792 | `split_half_raw` | 0.022772 | 33/33 de 33 |

### HAP1

| metrica | sentido | b (baseline) | r (replicate) | span orientado | estimador | r_sd (5 splits) | n_perts min/max |
|---|---|---|---|---|---|---|---|
| pds_cosine | mas alto mejor | 0.437500 | 0.789583 | 0.352083 | `split_half_raw` | 0.015251 | 33/33 de 33 |
| mse_capped_norm | mas bajo mejor | 1.003465 | 0.185635 | 0.817830 | `split_half_raw` | 0.091747 | None/None de 33 |
| lfc_nmae | mas bajo mejor | 0.897264 | 0.391590 | 0.505674 | `full_gate_raw` | 0.040663 | 8/8 de 33 |
| fid | mas alto mejor | 0.000000 | 0.546846 | 0.546846 | `split_half_raw` | 0.094299 | 19/26 de 33 |
| reach | mas alto mejor | 0.354552 | 0.725538 | 0.370986 | `split_half_raw` | 0.076025 | 16/20 de 33 |
| jac | mas alto mejor | 0.424242 | 0.468156 | 0.043914 | `split_half_raw` | 0.084462 | 33/33 de 33 |

### HT29

| metrica | sentido | b (baseline) | r (replicate) | span orientado | estimador | r_sd (5 splits) | n_perts min/max |
|---|---|---|---|---|---|---|---|
| pds_cosine | mas alto mejor | 0.501894 | 0.783523 | 0.281629 | `split_half_raw` | 0.017083 | 33/33 de 33 |
| mse_capped_norm | mas bajo mejor | 1.137395 | 0.132822 | 1.004573 | `split_half_raw` | 0.095598 | None/None de 33 |
| lfc_nmae | mas bajo mejor | 0.864894 | 0.483357 | 0.381536 | `full_gate_raw` | 0.013128 | 15/15 de 33 |
| fid | mas alto mejor | 0.114189 | 0.431793 | 0.317604 | `split_half_raw` | 0.065468 | 25/30 de 33 |
| reach | mas alto mejor | 0.327624 | 0.618619 | 0.290995 | `split_half_raw` | 0.040970 | 23/25 de 33 |
| jac | mas alto mejor | 0.091090 | 0.272116 | 0.181026 | `split_half_raw` | 0.058478 | 33/33 de 33 |

### K562

| metrica | sentido | b (baseline) | r (replicate) | span orientado | estimador | r_sd (5 splits) | n_perts min/max |
|---|---|---|---|---|---|---|---|
| pds_cosine | mas alto mejor | 0.520833 | 0.822727 | 0.301894 | `split_half_raw` | 0.017558 | 33/33 de 33 |
| mse_capped_norm | mas bajo mejor | 1.173286 | 0.126447 | 1.046839 | `split_half_raw` | 0.105315 | None/None de 33 |
| lfc_nmae | mas bajo mejor | 0.846558 | 0.497343 | 0.349215 | `full_gate_raw` | 0.008552 | 9/9 de 33 |
| fid | mas alto mejor | 0.076923 | 0.511705 | 0.434782 | `split_half_raw` | 0.044777 | 18/23 de 33 |
| reach | mas alto mejor | 0.406814 | 0.809303 | 0.402489 | `split_half_raw` | 0.054991 | 13/23 de 33 |
| jac | mas alto mejor | 0.212121 | 0.489996 | 0.277874 | `split_half_raw` | 0.067887 | 33/33 de 33 |

### MCF7

| metrica | sentido | b (baseline) | r (replicate) | span orientado | estimador | r_sd (5 splits) | n_perts min/max |
|---|---|---|---|---|---|---|---|
| pds_cosine | mas alto mejor | 0.503788 | 0.936174 | 0.432386 | `split_half_raw` | 0.016343 | 33/33 de 33 |
| mse_capped_norm | mas bajo mejor | 0.922312 | 0.104958 | 0.817354 | `split_half_raw` | 0.079845 | None/None de 33 |
| lfc_nmae | mas bajo mejor | 0.872196 | 0.663824 | 0.208372 | `full_gate_raw` | 0.011592 | 33/33 de 33 |
| fid | mas alto mejor | 0.000173 | 0.462605 | 0.462432 | `split_half_raw` | 0.090324 | 31/33 de 33 |
| reach | mas alto mejor | 0.101892 | 0.676395 | 0.574503 | `split_half_raw` | 0.063574 | 30/32 de 33 |
| jac | mas alto mejor | 0.000029 | 0.072243 | 0.072214 | `split_half_raw` | 0.024237 | 33/33 de 33 |

## El span `(r - b)` orientado -- ¿alguno sale <= 0?

Orientado por el sentido de cada metrica: `r - b` donde mas alto es mejor, `b - r` donde mas bajo es mejor. Un span <= 0 significa que el extremo 1 no esta del lado bueno del extremo 0. El marcador tiene su propia guarda: `score.py:466-467` levanta "the anchor leaves no usable headroom over the baseline".

| metrica | A549 | BXPC3 | HAP1 | HT29 | K562 | MCF7 | n <= 0 |
|---|---|---|---|---|---|---|---|
| pds_cosine | +0.291098 | +0.478788 | +0.352083 | +0.281629 | +0.301894 | +0.432386 | 0 |
| mse_capped_norm | +0.821267 | +0.755633 | +0.817830 | +1.004573 | +1.046839 | +0.817354 | 0 |
| lfc_nmae | +0.400857 | +0.267589 | +0.505674 | +0.381536 | +0.349215 | +0.208372 | 0 |
| fid | +0.417863 | +0.597809 | +0.546846 | +0.317604 | +0.434782 | +0.462432 | 0 |
| reach | +0.413137 | +0.489816 | +0.370986 | +0.290995 | +0.402489 | +0.574503 | 0 |
| jac | +0.018709 | +0.199792 | +0.043914 | +0.181026 | +0.277874 | +0.072214 | 0 |

**Ningun span sale <= 0: los 36 son positivos.** La escala se puede usar en las seis metricas y las seis lineas. Contrasta con nuestro techo a mano, donde `jac` daba span negativo en 3 de 6 lineas (`bucle_techo.md`) -- pero ver la seccion siguiente antes de darlo por resuelto.


## ⚠️ El span contra el ruido del propio `r`

Un span positivo no basta. Si el span es mas chico que la dispersion con que se midio su extremo 1, re-sortear las particiones mueve el techo mas que todo el recorrido que la escala pretende medir, y un candidato puntuado ahi recibe un numero que depende mas de la semilla del anchor que de su propio merito. `r_sd` es la desviacion de los 5 agregados por particion (`anchor.py:286`, `ddof=0`).

| linea | metrica | span | r_sd | span / r_sd | lectura |
|---|---|---|---|---|---|
| A549 | pds_cosine | 0.291098 | 0.029037 | 10.03 | holgado |
| A549 | mse_capped_norm | 0.821267 | 0.061351 | 13.39 | holgado |
| A549 | lfc_nmae | 0.400857 | 0.022871 | 17.53 | holgado |
| A549 | fid | 0.417863 | 0.036978 | 11.30 | holgado |
| A549 | reach | 0.413137 | 0.082754 | 4.99 | holgado |
| A549 | jac | 0.018709 | 0.035008 | 0.53 | **span < 1 sd: no usable** |
| BXPC3 | pds_cosine | 0.478788 | 0.016829 | 28.45 | holgado |
| BXPC3 | mse_capped_norm | 0.755633 | 0.063862 | 11.83 | holgado |
| BXPC3 | lfc_nmae | 0.267589 | 0.050873 | 5.26 | holgado |
| BXPC3 | fid | 0.597809 | 0.045793 | 13.05 | holgado |
| BXPC3 | reach | 0.489816 | 0.045994 | 10.65 | holgado |
| BXPC3 | jac | 0.199792 | 0.022772 | 8.77 | holgado |
| HAP1 | pds_cosine | 0.352083 | 0.015251 | 23.09 | holgado |
| HAP1 | mse_capped_norm | 0.817830 | 0.091747 | 8.91 | holgado |
| HAP1 | lfc_nmae | 0.505674 | 0.040663 | 12.44 | holgado |
| HAP1 | fid | 0.546846 | 0.094299 | 5.80 | holgado |
| HAP1 | reach | 0.370986 | 0.076025 | 4.88 | holgado |
| HAP1 | jac | 0.043914 | 0.084462 | 0.52 | **span < 1 sd: no usable** |
| HT29 | pds_cosine | 0.281629 | 0.017083 | 16.49 | holgado |
| HT29 | mse_capped_norm | 1.004573 | 0.095598 | 10.51 | holgado |
| HT29 | lfc_nmae | 0.381536 | 0.013128 | 29.06 | holgado |
| HT29 | fid | 0.317604 | 0.065468 | 4.85 | holgado |
| HT29 | reach | 0.290995 | 0.040970 | 7.10 | holgado |
| HT29 | jac | 0.181026 | 0.058478 | 3.10 | holgado |
| K562 | pds_cosine | 0.301894 | 0.017558 | 17.19 | holgado |
| K562 | mse_capped_norm | 1.046839 | 0.105315 | 9.94 | holgado |
| K562 | lfc_nmae | 0.349215 | 0.008552 | 40.83 | holgado |
| K562 | fid | 0.434782 | 0.044777 | 9.71 | holgado |
| K562 | reach | 0.402489 | 0.054991 | 7.32 | holgado |
| K562 | jac | 0.277874 | 0.067887 | 4.09 | holgado |
| MCF7 | pds_cosine | 0.432386 | 0.016343 | 26.46 | holgado |
| MCF7 | mse_capped_norm | 0.817354 | 0.079845 | 10.24 | holgado |
| MCF7 | lfc_nmae | 0.208372 | 0.011592 | 17.98 | holgado |
| MCF7 | fid | 0.462432 | 0.090324 | 5.12 | holgado |
| MCF7 | reach | 0.574503 | 0.063574 | 9.04 | holgado |
| MCF7 | jac | 0.072214 | 0.024237 | 2.98 | span < 3 sd: fragil |

**3 de 36 celdas tienen menos de 3 desviaciones de holgura**, y son las que hay que mirar antes de puntuar nada:

- **HAP1 / jac**: span 0.043914 contra sd 0.084462 = 0.52 sd
- **A549 / jac**: span 0.018709 contra sd 0.035008 = 0.53 sd
- **MCF7 / jac**: span 0.072214 contra sd 0.024237 = 2.98 sd

El caso extremo es **HAP1 / jac**: el span (0.043914) es 1.9x mas chico que la desviacion del propio `r` (0.084462). Ahi la escala no mide nada: el extremo 1 se mueve mas entre semillas que toda la distancia al extremo 0.

No es casualidad que `jac` concentre los casos finos: su extremo 0 esta inflado por el artefacto de union vacia, asi que `b` arranca alto y deja poca holgura hasta `r`. **El span positivo de `jac` no resuelve el problema que `bucle_techo.md` encontro: lo esconde detras de un signo.**


## La dispersion entre las 5 particiones del anchor

De `anchor_splits.parquet`, que guarda un valor por (particion, metrica) con su semilla derivada (`anchor.py:82-85`).

| linea | metrica | min | max | rango | sd | los 5 valores |
|---|---|---|---|---|---|---|
| A549 | pds_cosine | 0.724432 | 0.809659 | 0.085227 | 0.029037 | 0.8097, 0.7803, 0.7964, 0.7784, 0.7244 |
| A549 | mse_capped_norm | 0.033782 | 0.220074 | 0.186292 | 0.061351 | 0.1222, 0.0338, 0.1666, 0.1507, 0.2201 |
| A549 | lfc_nmae | 0.420630 | 0.486412 | 0.065782 | 0.022871 | 0.4206, 0.4404, 0.4671, 0.4646, 0.4864 |
| A549 | fid | 0.380629 | 0.466681 | 0.086052 | 0.036978 | 0.3812, 0.4038, 0.3806, 0.4667, 0.4570 |
| A549 | reach | 0.595830 | 0.820338 | 0.224508 | 0.082754 | 0.8056, 0.7347, 0.6813, 0.8203, 0.5958 |
| A549 | jac | 0.428102 | 0.516056 | 0.087954 | 0.035008 | 0.4909, 0.4354, 0.5161, 0.4281, 0.4959 |
| BXPC3 | pds_cosine | 0.929924 | 0.981061 | 0.051136 | 0.016829 | 0.9527, 0.9489, 0.9631, 0.9811, 0.9299 |
| BXPC3 | mse_capped_norm | -0.042190 | 0.152455 | 0.194646 | 0.063862 | 0.0964, 0.0811, 0.0560, -0.0422, 0.1525 |
| BXPC3 | lfc_nmae | 0.487332 | 0.617535 | 0.130203 | 0.050873 | 0.5844, 0.5541, 0.6175, 0.4926, 0.4873 |
| BXPC3 | fid | 0.521234 | 0.638259 | 0.117025 | 0.045793 | 0.6383, 0.6349, 0.5726, 0.6292, 0.5212 |
| BXPC3 | reach | 0.756822 | 0.897207 | 0.140385 | 0.045994 | 0.8972, 0.8207, 0.8475, 0.7568, 0.8503 |
| BXPC3 | jac | 0.163549 | 0.233679 | 0.070131 | 0.022772 | 0.1635, 0.2073, 0.1934, 0.2337, 0.2059 |
| HAP1 | pds_cosine | 0.765152 | 0.804924 | 0.039773 | 0.015251 | 0.7652, 0.7926, 0.7803, 0.8049, 0.8049 |
| HAP1 | mse_capped_norm | 0.077379 | 0.329965 | 0.252586 | 0.091747 | 0.3300, 0.1773, 0.2391, 0.1044, 0.0774 |
| HAP1 | lfc_nmae | 0.335637 | 0.430683 | 0.095046 | 0.040663 | 0.3356, 0.3492, 0.4267, 0.4156, 0.4307 |
| HAP1 | fid | 0.389352 | 0.684577 | 0.295225 | 0.094299 | 0.5641, 0.5620, 0.5342, 0.6846, 0.3894 |
| HAP1 | reach | 0.583333 | 0.808824 | 0.225490 | 0.076025 | 0.7458, 0.8088, 0.7278, 0.7619, 0.5833 |
| HAP1 | jac | 0.335480 | 0.567021 | 0.231541 | 0.084462 | 0.5241, 0.5670, 0.4075, 0.5067, 0.3355 |
| HT29 | pds_cosine | 0.767045 | 0.816288 | 0.049242 | 0.017083 | 0.8163, 0.7822, 0.7756, 0.7670, 0.7765 |
| HT29 | mse_capped_norm | -0.018699 | 0.249714 | 0.268413 | 0.095598 | -0.0187, 0.2115, 0.1417, 0.2497, 0.0799 |
| HT29 | lfc_nmae | 0.465065 | 0.504335 | 0.039270 | 0.013128 | 0.5043, 0.4815, 0.4765, 0.4651, 0.4894 |
| HT29 | fid | 0.367506 | 0.552873 | 0.185367 | 0.065468 | 0.5529, 0.3836, 0.4393, 0.4157, 0.3675 |
| HT29 | reach | 0.572941 | 0.675687 | 0.102746 | 0.040970 | 0.6121, 0.5729, 0.6549, 0.5775, 0.6757 |
| HT29 | jac | 0.167906 | 0.325295 | 0.157389 | 0.058478 | 0.3229, 0.2927, 0.2518, 0.3253, 0.1679 |
| K562 | pds_cosine | 0.803977 | 0.854167 | 0.050189 | 0.017558 | 0.8248, 0.8542, 0.8040, 0.8220, 0.8087 |
| K562 | mse_capped_norm | -0.030981 | 0.282206 | 0.313187 | 0.105315 | 0.2822, 0.0739, 0.1205, 0.1866, -0.0310 |
| K562 | lfc_nmae | 0.486471 | 0.508574 | 0.022102 | 0.008552 | 0.5086, 0.4908, 0.4950, 0.4865, 0.5059 |
| K562 | fid | 0.454968 | 0.571023 | 0.116055 | 0.044777 | 0.5418, 0.4550, 0.4650, 0.5710, 0.5258 |
| K562 | reach | 0.731624 | 0.892148 | 0.160524 | 0.054991 | 0.7917, 0.8921, 0.7316, 0.8458, 0.7853 |
| K562 | jac | 0.427257 | 0.573647 | 0.146390 | 0.067887 | 0.4334, 0.5736, 0.4435, 0.4273, 0.5721 |
| MCF7 | pds_cosine | 0.906250 | 0.952652 | 0.046402 | 0.016343 | 0.9527, 0.9489, 0.9384, 0.9347, 0.9062 |
| MCF7 | mse_capped_norm | -0.043988 | 0.192601 | 0.236588 | 0.079845 | 0.1262, -0.0440, 0.1055, 0.1926, 0.1444 |
| MCF7 | lfc_nmae | 0.642556 | 0.677673 | 0.035117 | 0.011592 | 0.6777, 0.6687, 0.6658, 0.6426, 0.6644 |
| MCF7 | fid | 0.354472 | 0.601739 | 0.247267 | 0.090324 | 0.4785, 0.6017, 0.5039, 0.3744, 0.3545 |
| MCF7 | reach | 0.574213 | 0.756330 | 0.182116 | 0.063574 | 0.7277, 0.6747, 0.7563, 0.6491, 0.5742 |
| MCF7 | jac | 0.052904 | 0.114673 | 0.061769 | 0.024237 | 0.0529, 0.1147, 0.0554, 0.0843, 0.0540 |

## 🔴 ¿La asimetria de direccion que medimos a mano era ruido de sorteo?

`bucle_techo.md` midio cada metrica en las dos direcciones (mitad A como prediccion y mitad B como prediccion) y reporto las dos sin promediarlas, porque un desempate dirigido le da mas potencia a un lado. La pregunta era si esa diferencia es senal o ruido del sorteo.

El anchor nativo la contesta sin suponer nada: corre **5 particiones** de la MISMA direccion (`half_b` como prediccion siempre, `anchor.py:120`). Si el rango entre esas 5 cubre la diferencia que vimos entre direcciones, la asimetria no era senal. Se compara contra la variante **I** del techo, la que comparte convencion de control con el anchor.

| linea | metrica | rango entre 5 particiones | \|dirA - dirB\| del techo I | ¿el rango lo cubre? |
|---|---|---|---|---|
| A549 | pds_cosine | 0.085227 | 0.059659 | **si** |
| A549 | mse_capped_norm | 0.186292 | 0.072305 | **si** |
| A549 | lfc_nmae | 0.065782 | 0.019438 | **si** |
| A549 | fid | 0.086052 | 0.067379 | **si** |
| A549 | reach | 0.224508 | 0.026710 | **si** |
| A549 | jac | 0.087954 | 0.000000 | **si** |
| BXPC3 | pds_cosine | 0.051136 | 0.005682 | **si** |
| BXPC3 | mse_capped_norm | 0.194646 | 0.007328 | **si** |
| BXPC3 | lfc_nmae | 0.130203 | 0.078577 | **si** |
| BXPC3 | fid | 0.117025 | 0.012418 | **si** |
| BXPC3 | reach | 0.140385 | 0.024119 | **si** |
| BXPC3 | jac | 0.070131 | 0.000000 | **si** |
| HAP1 | pds_cosine | 0.039773 | 0.017992 | **si** |
| HAP1 | mse_capped_norm | 0.252586 | 0.025853 | **si** |
| HAP1 | lfc_nmae | 0.095046 | 0.052027 | **si** |
| HAP1 | fid | 0.295225 | 0.024769 | **si** |
| HAP1 | reach | 0.225490 | 0.045222 | **si** |
| HAP1 | jac | 0.231541 | 0.000000 | **si** |
| HT29 | pds_cosine | 0.049242 | 0.067235 | no |
| HT29 | mse_capped_norm | 0.268413 | 0.031551 | **si** |
| HT29 | lfc_nmae | 0.039270 | 0.088958 | no |
| HT29 | fid | 0.185367 | 0.123307 | **si** |
| HT29 | reach | 0.102746 | 0.163658 | no |
| HT29 | jac | 0.157389 | 0.000000 | **si** |
| K562 | pds_cosine | 0.050189 | 0.046402 | **si** |
| K562 | mse_capped_norm | 0.313187 | 0.024249 | **si** |
| K562 | lfc_nmae | 0.022102 | 0.173434 | no |
| K562 | fid | 0.116055 | 0.116286 | no |
| K562 | reach | 0.160524 | 0.080307 | **si** |
| K562 | jac | 0.146390 | 0.000000 | **si** |
| MCF7 | pds_cosine | 0.046402 | 0.002841 | **si** |
| MCF7 | mse_capped_norm | 0.236588 | 0.012288 | **si** |
| MCF7 | lfc_nmae | 0.035117 | 0.162889 | no |
| MCF7 | fid | 0.247267 | 0.333888 | no |
| MCF7 | reach | 0.182116 | 0.187406 | no |
| MCF7 | jac | 0.061769 | 0.000000 | **si** |

**Resumen: el rango entre particiones cubre la asimetria de direccion en 28 de 36 celdas (linea x metrica); no la cubre en 8.**

La lectura: en la mayoria de las celdas la diferencia entre nuestras dos direcciones cae dentro de lo que una sola direccion se mueve al re-sortear la particion. **Mayormente era ruido de sorteo, y correr las dos direcciones no compraba la informacion que parecia.** Lo que si la compra es promediar varias particiones, que es lo que hace el camino nativo con k=5.

Las celdas donde la asimetria excede el ruido de sorteo:

- **K562 / lfc_nmae**: asimetria 0.173434 contra rango 0.022102 (excede por +0.151331)
- **MCF7 / lfc_nmae**: asimetria 0.162889 contra rango 0.035117 (excede por +0.127772)
- **MCF7 / fid**: asimetria 0.333888 contra rango 0.247267 (excede por +0.086622)
- **HT29 / reach**: asimetria 0.163658 contra rango 0.102746 (excede por +0.060912)
- **HT29 / lfc_nmae**: asimetria 0.088958 contra rango 0.039270 (excede por +0.049688)
- **HT29 / pds_cosine**: asimetria 0.067235 contra rango 0.049242 (excede por +0.017992)
- **MCF7 / reach**: asimetria 0.187406 contra rango 0.182116 (excede por +0.005290)
- **K562 / fid**: asimetria 0.116286 contra rango 0.116055 (excede por +0.000232)

⚠️ Dos advertencias: el rango de 5 valores es un estimador grosero de dispersion, no un intervalo de confianza; y las dos cantidades se miden sobre poblaciones distintas -- el techo partio el pool completo por mitades, el anchor parte esta poblacion de 400 por perturbacion. Es indicativo, no una prueba.


## Contra lo que ya teniamos

Dos comparaciones distintas, ninguna es una equivalencia:

- **nuestro suelo vs este `b`**: nuestro suelo emite el control SIN escalar; el `b` lo escala al perfil de respuesta promedio (`baseline_convencion.md` A.4).
- **nuestro techo (variante I, direccion B) vs este `r`**: la variante nuestra mas cercana a la convencion oficial, pero el oficial promedia k=5 particiones y el nuestro usa 1, y `lfc_nmae` oficial usa otro estimador.

| linea | metrica | nuestro suelo | b nativo | dif | nuestro techo I/B | r nativo | dif |
|---|---|---|---|---|---|---|---|
| A549 | pds_cosine | 0.474432 | 0.486742 | -0.012311 | 0.782197 | 0.777841 | +0.004356 |
| A549 | mse_capped_norm | 2.562172 | 0.959928 | +1.602244 | 0.094039 | 0.138661 | -0.044621 |
| A549 | lfc_nmae | 1.036018 | 0.856686 | +0.179332 | 0.473392 | 0.455829 | +0.017563 |
| A549 | fid | 0.000000 | 0.000000 | +0.000000 | 0.748485 | 0.417863 | +0.330621 |
| A549 | reach | 0.144540 | 0.314403 | -0.169863 | 0.922787 | 0.727540 | +0.195247 |
| A549 | jac | 0.454545 | 0.454545 | +0.000000 | 0.562245 | 0.473254 | +0.088992 |
| BXPC3 | pds_cosine | 0.517992 | 0.476326 | +0.041667 | 0.991477 | 0.955114 | +0.036364 |
| BXPC3 | mse_capped_norm | 1.332759 | 0.824367 | +0.508392 | 0.021559 | 0.068734 | -0.047176 |
| BXPC3 | lfc_nmae | 0.996965 | 0.814780 | +0.182185 | 0.628513 | 0.547191 | +0.081322 |
| BXPC3 | fid | 0.000026 | 0.001444 | -0.001418 | 0.778034 | 0.599253 | +0.178781 |
| BXPC3 | reach | 0.142889 | 0.344686 | -0.201797 | 0.943964 | 0.834502 | +0.109462 |
| BXPC3 | jac | 0.000000 | 0.000978 | -0.000978 | 0.211960 | 0.200769 | +0.011191 |
| HAP1 | pds_cosine | 0.512311 | 0.437500 | +0.074811 | 0.910985 | 0.789583 | +0.121402 |
| HAP1 | mse_capped_norm | 1.868728 | 1.003465 | +0.865262 | 0.025049 | 0.185635 | -0.160586 |
| HAP1 | lfc_nmae | 1.000821 | 0.897264 | +0.103558 | 0.362917 | 0.391590 | -0.028672 |
| HAP1 | fid | 0.000000 | 0.000000 | +0.000000 | 0.646667 | 0.546846 | +0.099821 |
| HAP1 | reach | 0.192408 | 0.354552 | -0.162144 | 0.845436 | 0.725538 | +0.119898 |
| HAP1 | jac | 0.424242 | 0.424242 | +0.000000 | 0.318423 | 0.468156 | -0.149733 |
| HT29 | pds_cosine | 0.490530 | 0.501894 | -0.011364 | 0.840909 | 0.783523 | +0.057386 |
| HT29 | mse_capped_norm | 1.931565 | 1.137395 | +0.794170 | 0.015348 | 0.132822 | -0.117474 |
| HT29 | lfc_nmae | 1.044314 | 0.864894 | +0.179420 | 0.642498 | 0.483357 | +0.159141 |
| HT29 | fid | 0.009259 | 0.114189 | -0.104930 | 0.407195 | 0.431793 | -0.024598 |
| HT29 | reach | 0.216033 | 0.327624 | -0.111592 | 0.704372 | 0.618619 | +0.085754 |
| HT29 | jac | 0.181818 | 0.091090 | +0.090728 | 0.133187 | 0.272116 | -0.138929 |
| K562 | pds_cosine | 0.495265 | 0.520833 | -0.025568 | 0.809659 | 0.822727 | -0.013068 |
| K562 | mse_capped_norm | 2.913227 | 1.173286 | +1.739941 | 0.153676 | 0.126447 | +0.027229 |
| K562 | lfc_nmae | 1.028846 | 0.846558 | +0.182289 | 0.518946 | 0.497343 | +0.021603 |
| K562 | fid | 0.000000 | 0.076923 | -0.076923 | 0.519768 | 0.511705 | +0.008063 |
| K562 | reach | 0.200448 | 0.406814 | -0.206366 | 0.743649 | 0.809303 | -0.065654 |
| K562 | jac | 0.242424 | 0.212121 | +0.030303 | 0.403280 | 0.489996 | -0.086716 |
| MCF7 | pds_cosine | 0.512311 | 0.503788 | +0.008523 | 0.985795 | 0.936174 | +0.049621 |
| MCF7 | mse_capped_norm | 1.774874 | 0.922312 | +0.852562 | 0.053492 | 0.104958 | -0.051466 |
| MCF7 | lfc_nmae | 1.056808 | 0.872196 | +0.184612 | 0.776169 | 0.663824 | +0.112345 |
| MCF7 | fid | 0.000000 | 0.000173 | -0.000173 | 0.350310 | 0.462605 | -0.112296 |
| MCF7 | reach | 0.013563 | 0.101892 | -0.088330 | 0.673735 | 0.676395 | -0.002660 |
| MCF7 | jac | 0.000000 | 0.000029 | -0.000029 | 0.095610 | 0.072243 | +0.023367 |

## Hallazgo: el artefacto de union vacia alcanza tambien al `b` OFICIAL

En varias celdas el `b` del marcador y nuestro suelo son **identicos al bit**, pese a ser emisiones distintas (una escalada al perfil de respuesta promedio, la otra el control sin escalar). Dos emisiones distintas dando el mismo numero solo puede significar que la metrica no esta mirando la emision:

| linea | metrica | b oficial | nuestro suelo | dif |
|---|---|---|---|---|
| A549 | fid | 0.000000 | 0.000000 | 0 |
| A549 | jac | 0.454545 | 0.454545 | 0 |
| HAP1 | fid | 0.000000 | 0.000000 | 0 |
| HAP1 | jac | 0.424242 | 0.424242 | 0 |

Para `jac` es el artefacto D-023 que `bucle_techo.md` ya habia aislado en nuestro suelo, operando igual en el baseline oficial: ninguna de las dos emisiones llama nada significativo en las perturbaciones sin senal, la union queda vacia, y `de_sig_jaccard` devuelve 1.0 por convencion (`metrics/de.py:195`). El valor es el conteo de uniones vacias sobre 33, sea quien sea el que emite:

- **A549**: 0.454545 = 15/33 uniones vacias
- **HAP1**: 0.424242 = 14/33 uniones vacias

**Consecuencia:** el extremo 0 de `jac` no distingue un comparador oraculo de una emision degenerada. La conclusion de `bucle_techo.md` sobre `jac` no era un defecto de nuestro suelo: es una propiedad de la metrica, y vale igual para el camino nativo.

Para `fid` la coincidencia es en 0.0 exacto: ninguna de las dos emisiones consigue una sola llamada de direccion correcta con rendimiento, asi que el extremo 0 esta pegado al piso en esas lineas -- justo donde `scales.py:430` ubica una base constante de 0.5.


## El estimador estampado por metrica

`anchor.py:296-299` estampa, por metrica, que estimador produjo su ancla: `full_gate_raw` para la familia `lfc_nmae` y `split_half_raw` para el resto. El motivo (`anchor.py:57-66`): una MITAD llama muchos menos genes significativos que el dato completo (mediana del gate 1.83-3.36x menor), asi que el nucleo uniforme promediaria sobre 65-79% de las perturbaciones que promedia la metrica. `full_gate_raw` toma el gate y el denominador de la tabla DE del real COMPLETO y solo los dos vectores del numerador de las mitades.

| linea | pds_cosine | mse_capped_norm | lfc_nmae | fid | reach | jac |
|---|---|---|---|---|---|---|
| A549 | `split_half_raw` | `split_half_raw` | `full_gate_raw` | `split_half_raw` | `split_half_raw` | `split_half_raw` |
| BXPC3 | `split_half_raw` | `split_half_raw` | `full_gate_raw` | `split_half_raw` | `split_half_raw` | `split_half_raw` |
| HAP1 | `split_half_raw` | `split_half_raw` | `full_gate_raw` | `split_half_raw` | `split_half_raw` | `split_half_raw` |
| HT29 | `split_half_raw` | `split_half_raw` | `full_gate_raw` | `split_half_raw` | `split_half_raw` | `split_half_raw` |
| K562 | `split_half_raw` | `split_half_raw` | `full_gate_raw` | `split_half_raw` | `split_half_raw` | `split_half_raw` |
| MCF7 | `split_half_raw` | `split_half_raw` | `full_gate_raw` | `split_half_raw` | `split_half_raw` | `split_half_raw` |

## 20 filas trabajadas, con los 5 valores por particion a la vista

Veinte filas (linea, metrica) con TODO lo que entra en su span: los dos extremos, los 5 valores por particion que promedian al `r`, su dispersion, el estimador y la cohorte. No hay paso intermedio oculto entre el artefacto del marcador y el span.

| # | linea | metrica | b | r | span | los 5 valores por particion | media recalculada | dif vs r | rango | sd | span/sd | estimador | n_perts |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | HAP1 | jac | 0.424242 | 0.468156 | 0.043914 | 0.5241, 0.5670, 0.4075, 0.5067, 0.3355 | 0.468156 | 0.00e+00 | 0.231541 | 0.084462 | 0.52 | `split_half_raw` | 33/33 de 33 |
| 2 | A549 | jac | 0.454545 | 0.473254 | 0.018709 | 0.4909, 0.4354, 0.5161, 0.4281, 0.4959 | 0.473254 | 0.00e+00 | 0.087954 | 0.035008 | 0.53 | `split_half_raw` | 33/33 de 33 |
| 3 | MCF7 | jac | 0.000029 | 0.072243 | 0.072214 | 0.0529, 0.1147, 0.0554, 0.0843, 0.0540 | 0.072243 | 0.00e+00 | 0.061769 | 0.024237 | 2.98 | `split_half_raw` | 33/33 de 33 |
| 4 | HT29 | jac | 0.091090 | 0.272116 | 0.181026 | 0.3229, 0.2927, 0.2518, 0.3253, 0.1679 | 0.272116 | 0.00e+00 | 0.157389 | 0.058478 | 3.10 | `split_half_raw` | 33/33 de 33 |
| 5 | K562 | jac | 0.212121 | 0.489996 | 0.277874 | 0.4334, 0.5736, 0.4435, 0.4273, 0.5721 | 0.489996 | 1.11e-16 | 0.146390 | 0.067887 | 4.09 | `split_half_raw` | 33/33 de 33 |
| 6 | HT29 | fid | 0.114189 | 0.431793 | 0.317604 | 0.5529, 0.3836, 0.4393, 0.4157, 0.3675 | 0.431793 | 0.00e+00 | 0.185367 | 0.065468 | 4.85 | `split_half_raw` | 25/30 de 33 |
| 7 | K562 | mse_capped_norm | 1.173286 | 0.126447 | 1.046839 | 0.2822, 0.0739, 0.1205, 0.1866, -0.0310 | 0.126447 | 0.00e+00 | 0.313187 | 0.105315 | 9.94 | `split_half_raw` | None/None de 33 |
| 8 | HAP1 | fid | 0.000000 | 0.546846 | 0.546846 | 0.5641, 0.5620, 0.5342, 0.6846, 0.3894 | 0.546846 | 0.00e+00 | 0.295225 | 0.094299 | 5.80 | `split_half_raw` | 19/26 de 33 |
| 9 | HT29 | mse_capped_norm | 1.137395 | 0.132822 | 1.004573 | -0.0187, 0.2115, 0.1417, 0.2497, 0.0799 | 0.132822 | 0.00e+00 | 0.268413 | 0.095598 | 10.51 | `split_half_raw` | None/None de 33 |
| 10 | HAP1 | mse_capped_norm | 1.003465 | 0.185635 | 0.817830 | 0.3300, 0.1773, 0.2391, 0.1044, 0.0774 | 0.185635 | 0.00e+00 | 0.252586 | 0.091747 | 8.91 | `split_half_raw` | None/None de 33 |
| 11 | MCF7 | fid | 0.000173 | 0.462605 | 0.462432 | 0.4785, 0.6017, 0.5039, 0.3744, 0.3545 | 0.462605 | 0.00e+00 | 0.247267 | 0.090324 | 5.12 | `split_half_raw` | 31/33 de 33 |
| 12 | MCF7 | mse_capped_norm | 0.922312 | 0.104958 | 0.817354 | 0.1262, -0.0440, 0.1055, 0.1926, 0.1444 | 0.104958 | 0.00e+00 | 0.236588 | 0.079845 | 10.24 | `split_half_raw` | None/None de 33 |
| 13 | A549 | lfc_nmae | 0.856686 | 0.455829 | 0.400857 | 0.4206, 0.4404, 0.4671, 0.4646, 0.4864 | 0.455829 | 0.00e+00 | 0.065782 | 0.022871 | 17.53 | `full_gate_raw` | 8/8 de 33 |
| 14 | BXPC3 | lfc_nmae | 0.814780 | 0.547191 | 0.267589 | 0.5844, 0.5541, 0.6175, 0.4926, 0.4873 | 0.547191 | 0.00e+00 | 0.130203 | 0.050873 | 5.26 | `full_gate_raw` | 23/23 de 33 |
| 15 | HAP1 | lfc_nmae | 0.897264 | 0.391590 | 0.505674 | 0.3356, 0.3492, 0.4267, 0.4156, 0.4307 | 0.391590 | 0.00e+00 | 0.095046 | 0.040663 | 12.44 | `full_gate_raw` | 8/8 de 33 |
| 16 | HT29 | lfc_nmae | 0.864894 | 0.483357 | 0.381536 | 0.5043, 0.4815, 0.4765, 0.4651, 0.4894 | 0.483357 | 0.00e+00 | 0.039270 | 0.013128 | 29.06 | `full_gate_raw` | 15/15 de 33 |
| 17 | A549 | pds_cosine | 0.486742 | 0.777841 | 0.291098 | 0.8097, 0.7803, 0.7964, 0.7784, 0.7244 | 0.777841 | 0.00e+00 | 0.085227 | 0.029037 | 10.03 | `split_half_raw` | 33/33 de 33 |
| 18 | A549 | mse_capped_norm | 0.959928 | 0.138661 | 0.821267 | 0.1222, 0.0338, 0.1666, 0.1507, 0.2201 | 0.138661 | 2.78e-17 | 0.186292 | 0.061351 | 13.39 | `split_half_raw` | None/None de 33 |
| 19 | A549 | fid | 0.000000 | 0.417863 | 0.417863 | 0.3812, 0.4038, 0.3806, 0.4667, 0.4570 | 0.417863 | 0.00e+00 | 0.086052 | 0.036978 | 11.30 | `split_half_raw` | 18/20 de 33 |
| 20 | A549 | reach | 0.314403 | 0.727540 | 0.413137 | 0.8056, 0.7347, 0.6813, 0.8203, 0.5958 | 0.727540 | 0.00e+00 | 0.224508 | 0.082754 | 4.99 | `split_half_raw` | 14/18 de 33 |

`dif vs r` es |media de los 5 valores - `replicate` del artefacto|: cierra el circuito entre `anchor_splits.parquet` y `anchor_agg.parquet`. Si no fuera ~0, la columna de dispersion no describiria el `r` que se reporta.

Las 5 semillas derivadas son las mismas en todas las lineas y metricas, por construccion: `[2968811710, 3677149159, 745650761, 2884920346, 2642120001]`, de `numpy.random.SeedSequence(0).generate_state(5)` (`anchor.py:74`, `135-141`).


## Intermedios de corrida

| linea | real_shape | n control | seed_real | anchor base seed | k | baseline_s | bundle_s | total_s | free_gib |
|---|---|---|---|---|---|---|---|---|---|
| A549 | [14579, 34025] | 1379 | 20260001 | 0 | 5 | 62.6 | 236.52 | 300.71 | 36.12 |
| BXPC3 | [16428, 34025] | 3228 | 20260002 | 0 | 5 | 84.68 | 361.95 | 449.45 | 36.12 |
| HAP1 | [16372, 34025] | 3172 | 20260003 | 0 | 5 | 87.73 | 381.49 | 471.95 | 36.12 |
| HT29 | [15936, 34025] | 2736 | 20260004 | 0 | 5 | 77.4 | 324.0 | 403.65 | 36.12 |
| K562 | [14598, 34025] | 1398 | 20260005 | 0 | 5 | 65.9 | 257.82 | 325.35 | 36.12 |
| MCF7 | [15869, 34025] | 2669 | 20260006 | 0 | 5 | 85.43 | 366.59 | 454.81 | 36.12 |

La poblacion real de cada linea se armo reusando `suelo_arnes.build_fold_adatas` **sin modificarlo**, con las MISMAS semillas del suelo (`seed_real = 20260001 + fold_idx`), de modo que es celda por celda la que puntuo `bucle_suelo.md`: 400 celulas por perturbacion sin reemplazo mas TODO el pool de control. El contrato de `build_fold_adatas` exige una funcion de emision, asi que se le paso `predict_floor` y **se descarto su salida**. Etiqueta de control usada: `non-targeting`.

⚠️ El lado real del TECHO no es esta poblacion: eran mitades del pool completo, no 400 por perturbacion. Toda comparacion techo-vs-`r` cruza dos poblaciones distintas ademas de dos convenciones.


## Lo que NO se hizo

- **No se corrio ningun candidato.** Aqui solo se construyen los dos extremos.
- **No se uso `--anchor` suelto**: es diagnostico y nunca toma `avg_score`. El `r` sale del bundle.
- **No se calculo ninguna metrica a mano.** `b` y `r` se leen de `baseline_agg.csv` y `anchor_agg.parquet`. La unica aritmetica propia es el span (una resta), la media recalculada de las particiones (para verificar el artefacto contra si mismo) y los cocientes span/sd.
- **No se borro ni modifico** `bucle_suelo.md`, `bucle_techo.md`, `suelo_arnes.py`, `suelo_floor.py`, `bucle_suelo.py`, `techo_arnes.py` ni `bucle_techo.py`.
- No se leyo nada bajo `~/vcc/benchmark/` ni `/datos/benchmark_VEDADO_techo_metricas/`.
- Sin envios, sin `vcc sample`, sin commits.

