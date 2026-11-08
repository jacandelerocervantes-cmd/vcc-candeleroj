# Camino nativo, PRIMER intento -- detenido por la guarda

⚠️ **Este archivo es el registro de por que el primer intento se detuvo.** Los seis bundles salieron DIAGNOSTICO y el flujo paro antes de leer ningun anchor, asi que **aqui no hay `r`**. El segundo intento, con la etiqueta de control de competencia y los seis bundles validos, esta en **`modelo/out/anclas_nativas_v2.md`**.

Construido con `cell-eval2 baseline` y `cell-eval2 prep-real-bundle`, una vez por linea celular de IFNB, pasando `--pert-col target_gene --control NT`. Nada se calculo a mano.

## 🔴 El `b` es TRANSDUCTIVO: artefacto de puntuacion, nunca de modelo

`baseline.py:6-10`:

> "It is an **ORACLE comparator, not a floor a submission could reach**: the profile is averaged from the evaluated real perturbations, so it is transductive and no real submission has access to that data. What it bounds is a metric's TRIVIALITY -- a score at or below it is not evidence of target-specific skill."

El perfil es `profile[g] = mean over { p : p != control, label(p) != g } of pseudobulk[p][g]` (`baseline.py:119`), promediado sobre las perturbaciones reales que se estan evaluando; la emision lo alcanza escalando el pool de control por `r = profile/ctrl` (`baseline.py:344`).

**Queda declarado: ninguna cifra de este artefacto entra nunca a un modelo.** Ni el perfil, ni la escala `r = profile/ctrl`, ni los valores de `b`, ni nada derivado de ellos. Son el denominador contra el que se puntua; usarlos como entrada seria entrenar contra la respuesta de evaluacion. El `r` tampoco es informacion que una submission tenga, y su uso aqui es igualmente solo de puntuacion.


## Guarda -- los seis bundles salieron DIAGNOSTICO (linea literal)

`cli.py:788-793` imprime esta linea en la primera construccion, y el comentario del codigo dice por que: "A miscomputed rule check would otherwise un-enrol every submission scored against this bundle, silently and for weeks; here it is visible on the first build." Capturada tal cual:

| linea | linea literal impresa por `prep-real-bundle` | `rule_digest` |
|---|---|---|
| A549 | `wrote DIAGNOSTIC bundle 'BUNDLE_A549' to /datos/modelo/out/nativas/BUNDLE_A549 -- it will NOT take avg_score: the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs` | null |
| BXPC3 | `wrote DIAGNOSTIC bundle 'BUNDLE_BXPC3' to /datos/modelo/out/nativas/BUNDLE_BXPC3 -- it will NOT take avg_score: the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs` | null |
| HAP1 | `wrote DIAGNOSTIC bundle 'BUNDLE_HAP1' to /datos/modelo/out/nativas/BUNDLE_HAP1 -- it will NOT take avg_score: the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs` | null |
| HT29 | `wrote DIAGNOSTIC bundle 'BUNDLE_HT29' to /datos/modelo/out/nativas/BUNDLE_HT29 -- it will NOT take avg_score: the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs` | null |
| K562 | `wrote DIAGNOSTIC bundle 'BUNDLE_K562' to /datos/modelo/out/nativas/BUNDLE_K562 -- it will NOT take avg_score: the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs` | null |
| MCF7 | `wrote DIAGNOSTIC bundle 'BUNDLE_MCF7' to /datos/modelo/out/nativas/BUNDLE_MCF7 -- it will NOT take avg_score: the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs` | null |

Un bundle DIAGNOSTICO **no toma `avg_score`** (del `--help` de `score`: "With a competition bundle, from_replicate takes avg_score"). El motivo estampado fue el mismo en las seis:

- **A549**: ["the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs"]
- **BXPC3**: ["the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs"]
- **HAP1**: ["the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs"]
- **HT29**: ["the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs"]
- **K562**: ["the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs"]
- **MCF7**: ["the run config is not the competition preset's (diff it against `--preset vcc2026`); its config_hash differs"]

### La causa, aislada por aritmetica de hashes

La comprobacion es `rule_config_hash(config) != p["config_hash"]` (`competition.py:504-508`). Las otras dos condiciones de `is_competition_rule` -- perfil de metricas y `cache_strict` (`competition.py:513-522`) -- **pasaron**: por eso el manifest estampa un unico motivo y no tres.

`competition.py:112` declara que campos se normalizan antes de hashear:

```
_RULE_EXCLUDED = ("device", "pert_chunk", "pert_col")
```

`pert_col` esta EXCLUIDO; `control` **no**. Los cuatro hashes, recalculados al escribir este reporte:

Hash de referencia: `a5519017556705abf9a6e553ebb4e517b9612e9546331ace7c5ec6c5c30a1a13`

| config | `rule_config_hash` | ¿es el de competencia? |
|---|---|---|
| preset `vcc2026` tal cual | `a5519017556705ab...` | **si** |
| preset + `pert_col='target_gene'` | `a5519017556705ab...` | **si** |
| preset + `control='NT'` | `d4b0b3f893ed007f...` | **NO** |
| lo que corrimos (los dos cambios) | `d4b0b3f893ed007f...` | **NO** |

El hash de "solo `control`" es **identico** al de la invocacion completa, lo que prueba que `--pert-col target_gene` no aporta nada al desvio. **El unico culpable es `--control NT`**: la regla exige `control: non-targeting` (`configs/vcc2026.yaml`).

## Guarda -- `source_fingerprint_strict` en cada artefacto

`configs/vcc2026.yaml:8-11` lo declara "LOAD-BEARING, not hygiene: `score.expect_from_run_meta` refuses a metadata-only fingerprint ... A competition run built without this cannot be scored against a real bundle at all."

| linea | baseline: strict | manifest: strict | fingerprint baseline | fingerprint manifest | coinciden |
|---|---|---|---|---|---|
| A549 | True | True | `1739af26959640c8...` | `1739af26959640c8...` | si |
| BXPC3 | True | True | `66fa0f672abcbdfa...` | `66fa0f672abcbdfa...` | si |
| HAP1 | True | True | `7994bd8f30ad225d...` | `7994bd8f30ad225d...` | si |
| HT29 | True | True | `a8527a6a46389086...` | `a8527a6a46389086...` | si |
| K562 | True | True | `4998eb921788c280...` | `4998eb921788c280...` | si |
| MCF7 | True | True | `c2d5b2eede847fee...` | `c2d5b2eede847fee...` | si |

## El `b` de las seis lineas (lo unico que quedo construido)

La fila `mean` de `baseline_agg.csv` del bundle -- la misma que lee `score` (`score.py:731`).

| metrica | A549 | BXPC3 | HAP1 | HT29 | K562 | MCF7 |
|---|---|---|---|---|---|---|
| pds_cosine | 0.486742 | 0.476326 | 0.437500 | 0.501894 | 0.520833 | 0.503788 |
| mse_capped_norm | 0.959928 | 0.824367 | 1.003465 | 1.137395 | 1.173286 | 0.922312 |
| lfc_nmae | 0.856686 | 0.814780 | 0.897264 | 0.864894 | 0.846558 | 0.872196 |
| fid | 0.000000 | 0.001444 | 0.000000 | 0.114189 | 0.076923 | 0.000173 |
| reach | 0.314403 | 0.344686 | 0.354552 | 0.327624 | 0.406814 | 0.101892 |
| jac | 0.454545 | 0.000978 | 0.424242 | 0.091090 | 0.212121 | 0.000029 |

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


## 20 filas trabajadas, con los intermedios completos

No hay lado `r`, asi que las filas se trabajan sobre el lado `b`, con su cruce entre el `baseline_agg.csv` del bundle y el de la corrida suelta de `cell-eval2 baseline`. Ese cruce no es cosmetico: el bundle **re-puntua** la pata del baseline con su propia config (`real_bundle.py:449`), asi que si los dos difirieran, el `b` que leeria `score` no seria el que imprimio `baseline`.

| # | linea | metrica | b del BUNDLE | b de la corrida suelta | dif | nuestro suelo | suelo - b | n control | real_shape |
|---|---|---|---|---|---|---|---|---|---|
| 1 | A549 | pds_cosine | 0.486742 | 0.486742 | 0.00e+00 | 0.474432 | -0.012311 | 1379 | [14579, 34025] |
| 2 | A549 | mse_capped_norm | 0.959928 | 0.959928 | 0.00e+00 | 2.562172 | +1.602244 | 1379 | [14579, 34025] |
| 3 | A549 | lfc_nmae | 0.856686 | 0.856686 | 0.00e+00 | 1.036018 | +0.179332 | 1379 | [14579, 34025] |
| 4 | A549 | fid | 0.000000 | 0.000000 | 0.00e+00 | 0.000000 | +0.000000 | 1379 | [14579, 34025] |
| 5 | A549 | reach | 0.314403 | 0.314403 | 0.00e+00 | 0.144540 | -0.169863 | 1379 | [14579, 34025] |
| 6 | A549 | jac | 0.454545 | 0.454545 | 0.00e+00 | 0.454545 | +0.000000 | 1379 | [14579, 34025] |
| 7 | BXPC3 | pds_cosine | 0.476326 | 0.476326 | 0.00e+00 | 0.517992 | +0.041667 | 3228 | [16428, 34025] |
| 8 | BXPC3 | mse_capped_norm | 0.824367 | 0.824367 | 0.00e+00 | 1.332759 | +0.508392 | 3228 | [16428, 34025] |
| 9 | BXPC3 | lfc_nmae | 0.814780 | 0.814780 | 2.22e-16 | 0.996965 | +0.182185 | 3228 | [16428, 34025] |
| 10 | BXPC3 | fid | 0.001444 | 0.001444 | 0.00e+00 | 0.000026 | -0.001418 | 3228 | [16428, 34025] |
| 11 | BXPC3 | reach | 0.344686 | 0.344686 | 0.00e+00 | 0.142889 | -0.201797 | 3228 | [16428, 34025] |
| 12 | BXPC3 | jac | 0.000978 | 0.000978 | 0.00e+00 | 0.000000 | -0.000978 | 3228 | [16428, 34025] |
| 13 | HAP1 | pds_cosine | 0.437500 | 0.437500 | 0.00e+00 | 0.512311 | +0.074811 | 3172 | [16372, 34025] |
| 14 | HAP1 | mse_capped_norm | 1.003465 | 1.003465 | 0.00e+00 | 1.868728 | +0.865262 | 3172 | [16372, 34025] |
| 15 | HAP1 | lfc_nmae | 0.897264 | 0.897264 | 1.11e-16 | 1.000821 | +0.103558 | 3172 | [16372, 34025] |
| 16 | HAP1 | fid | 0.000000 | 0.000000 | 0.00e+00 | 0.000000 | +0.000000 | 3172 | [16372, 34025] |
| 17 | HAP1 | reach | 0.354552 | 0.354552 | 0.00e+00 | 0.192408 | -0.162144 | 3172 | [16372, 34025] |
| 18 | HAP1 | jac | 0.424242 | 0.424242 | 0.00e+00 | 0.424242 | +0.000000 | 3172 | [16372, 34025] |
| 19 | HT29 | pds_cosine | 0.501894 | 0.501894 | 0.00e+00 | 0.490530 | -0.011364 | 2736 | [15936, 34025] |
| 20 | HT29 | mse_capped_norm | 1.137395 | 1.137395 | 0.00e+00 | 1.931565 | +0.794170 | 2736 | [15936, 34025] |

`dif` cierra el circuito entre los dos `baseline_agg.csv`. `suelo - b` es la distancia entre nuestro suelo (control sin escalar) y el baseline del marcador (control escalado al perfil): no es un error, es la medida de que son comparadores distintos.


## Intermedios de corrida

| linea | real_shape | n control | seed_real | anchor base seed | k | baseline_s | bundle_s | total_s | free_gib |
|---|---|---|---|---|---|---|---|---|---|
| A549 | [14579, 34025] | 1379 | 20260001 | 0 | 5 | 62.33 | 236.7 | 300.6 | 36.13 |
| BXPC3 | [16428, 34025] | 3228 | 20260002 | 0 | 5 | 84.49 | 361.65 | 448.95 | 36.13 |
| HAP1 | [16372, 34025] | 3172 | 20260003 | 0 | 5 | 87.8 | 380.56 | 471.1 | 36.13 |
| HT29 | [15936, 34025] | 2736 | 20260004 | 0 | 5 | 77.27 | 320.32 | 399.81 | 36.13 |
| K562 | [14598, 34025] | 1398 | 20260005 | 0 | 5 | 65.97 | 255.27 | 322.88 | 36.13 |
| MCF7 | [15869, 34025] | 2669 | 20260006 | 0 | 5 | 85.36 | 369.08 | 457.22 | 36.13 |

La poblacion real de cada linea se armo reusando `suelo_arnes.build_fold_adatas` **sin modificarlo**, con las MISMAS semillas del suelo (`seed_real = 20260001 + fold_idx`), de modo que es celda por celda la que puntuo `bucle_suelo.md`: 400 celulas por perturbacion sin reemplazo mas TODO el pool de control. El contrato de `build_fold_adatas` exige una funcion de emision, asi que se le paso `predict_floor` y **se descarto su salida**. Etiqueta de control usada: `NT`.

⚠️ El lado real del TECHO no es esta poblacion: eran mitades del pool completo, no 400 por perturbacion. Toda comparacion techo-vs-`r` cruza dos poblaciones distintas ademas de dos convenciones.


## Lo que NO se hizo

- **No se corrio ningun candidato.** Aqui solo se construyen los dos extremos.
- **No se uso `--anchor` suelto**: es diagnostico y nunca toma `avg_score`. El `r` sale del bundle.
- **No se calculo ninguna metrica a mano.** `b` y `r` se leen de `baseline_agg.csv` y `anchor_agg.parquet`. La unica aritmetica propia es el span (una resta), la media recalculada de las particiones (para verificar el artefacto contra si mismo) y los cocientes span/sd.
- **No se borro ni modifico** `bucle_suelo.md`, `bucle_techo.md`, `suelo_arnes.py`, `suelo_floor.py`, `bucle_suelo.py`, `techo_arnes.py` ni `bucle_techo.py`.
- No se leyo nada bajo `~/vcc/benchmark/` ni `/datos/benchmark_VEDADO_techo_metricas/`.
- Sin envios, sin `vcc sample`, sin commits.


---

**Continuacion:** `modelo/out/anclas_nativas_v2.md`, con la etiqueta `non-targeting`, los seis bundles de competencia y el `r`.

