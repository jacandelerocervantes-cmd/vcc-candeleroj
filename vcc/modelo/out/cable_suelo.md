# La prueba del cable: `run` -> `score --real-bundle`

No es un candidato. La emision es **desechable**: `suelo_floor.predict_floor` sin modificar -- el pool de control remuestreado con reemplazo, sin escalar. El objetivo es estrenar la cadena completa con una emision de la que ya sabemos la respuesta, porque `run` tiene la misma sensibilidad de `config_hash` que mando los seis primeros bundles a DIAGNOSTICO (ver `modelo/out/anclas_nativas.md`).

## La prediccion, declarada ANTES de correr

> `predict_floor` emite el control remuestreado sin escalar. El `b` oficial es el control ESCALADO al perfil de respuesta promedio (`baseline.py:344`). Asi que el suelo deberia quedar POR DEBAJO de `b` en casi todo, y su `s` salir ~0 o negativo. Si sale muy por encima de 0, algo esta mal en la cadena, no en el suelo.

Esto es lo que hace que la prueba sirva: el resultado estaba comprometido antes de verlo.

## ¿`score` corrio o se nego?

✅ **`score --real-bundle` corrio en 6 de 6 lineas**, sin que ninguna guarda se disparara.

| linea | run | score | duracion score | `anchor_source` | `real_bundle_id` |
|---|---|---|---|---|---|
| A549 | exit 0 (60.88s) | exit 0 | 3.05s | `real_bundle` | `BUNDLE_A549` |
| BXPC3 | exit 0 (81.05s) | exit 0 | 3.04s | `real_bundle` | `BUNDLE_BXPC3` |
| HAP1 | exit 0 (84.4s) | exit 0 | 3.04s | `real_bundle` | `BUNDLE_HAP1` |
| HT29 | exit 0 (74.37s) | exit 0 | 3.06s | `real_bundle` | `BUNDLE_HT29` |
| K562 | exit 0 (64.38s) | exit 0 | 3.04s | `real_bundle` | `BUNDLE_K562` |
| MCF7 | exit 0 (82.42s) | exit 0 | 3.05s | `real_bundle` | `BUNDLE_MCF7` |

Las dos guardas que ya conociamos **no se dispararon**, y conviene decir por que cada una paso:

- **`expect_from_run_meta` sobre el fingerprint** (`score.py:382-389`): exige `source_fingerprint_strict`. El preset `vcc2026` trae `cache_strict: true`, asi que el `run_meta.json` de cada corrida lo lleva en `True` (tabla de abajo). Si hubiera sido metadata-only, `score` habria levantado "this run's source_fingerprint is metadata-only ... Re-run both sides with `--cache-strict`".
- **"the anchor leaves no usable headroom over the baseline"** (`score.py:466-467`): exige que el span sea positivo para cada metrica decisiva. Los 36 spans del bundle son positivos (`anclas_nativas_v2.md`), asi que no habia nada que levantar -- pero ver la seccion de las dos celdas que no resuelven.

| linea | `source_fingerprint_strict` del run | fingerprint del run |
|---|---|---|
| A549 | True | `5de2fa2f999f7ea1...` |
| BXPC3 | True | `49441654cd104395...` |
| HAP1 | True | `701e7fe28e509b70...` |
| HT29 | True | `4e9d4aab3846bcce...` |
| K562 | True | `5f2bb67ddc6506aa...` |
| MCF7 | True | `e102151004dd5556...` |

## ¿El `run` quedo con `config_hash` de competencia?

`run` no imprime una linea de competencia como `prep-real-bundle`, asi que el hash se recalcula desde el `run_params.yaml` que escribio cada corrida y se compara contra el de referencia `a5519017556705ab...`.

| linea | hash obtenido | coincide | `control` en el config | `pert_col` | `cache_strict` |
|---|---|---|---|---|---|
| A549 | `a5519017556705ab...` | si | `non-targeting` | `target_gene` | True |
| BXPC3 | `a5519017556705ab...` | si | `non-targeting` | `target_gene` | True |
| HAP1 | `a5519017556705ab...` | si | `non-targeting` | `target_gene` | True |
| HT29 | `a5519017556705ab...` | si | `non-targeting` | `target_gene` | True |
| K562 | `a5519017556705ab...` | si | `non-targeting` | `target_gene` | True |
| MCF7 | `a5519017556705ab...` | si | `non-targeting` | `target_gene` | True |

✅ **Los seis runs coinciden con el hash de competencia.** La leccion del primer intento quedo aplicada: no se paso `--control` (el valor queda el del preset, `non-targeting`) y `--pert-col target_gene` si se paso, porque esta en `_RULE_EXCLUDED` (`competition.py:112`). La sensibilidad de `config_hash` que mando los seis primeros bundles a DIAGNOSTICO existe igual en `run`, y aqui no se activo.


## `avg_score` por linea, y de donde viene

El `--help` de `score` dice que con un bundle de competencia "from_replicate takes avg_score and from_baseline's goes null". Medido:

| linea | `avg_score` de `from_replicate` | `avg_score` de `from_baseline` | ¿cual lo tomo? | `anchor_digest` |
|---|---|---|---|---|
| A549 | -0.150136 | null | `from_replicate` | `b3ba695211d5a87d...` |
| BXPC3 | -0.168844 | null | `from_replicate` | `95391b425f312b9b...` |
| HAP1 | -0.071562 | null | `from_replicate` | `836b6f2fb6cb9b4b...` |
| HT29 | -0.120547 | null | `from_replicate` | `5f02c9373faf9953...` |
| K562 | -0.197881 | null | `from_replicate` | `bdb5448803c232ef...` |
| MCF7 | -0.170131 | null | `from_replicate` | `ec71abf164e77242...` |

**`from_replicate` tomo el `avg_score` en las seis lineas, y el de `from_baseline` salio null**, exactamente como declara el `--help`. Las columnas por metrica de `from_baseline` SI estan pobladas; lo que queda null es su fila `avg_score`. Es la diferencia entre "la columna existe" y "la columna decide".


## El `s` por metrica y por linea, las 6 x 6

`from_replicate`, que es la columna que toma el `avg_score`. Entre parentesis el `from_baseline` de la misma celda, para contraste.

| metrica | A549 | BXPC3 | HAP1 | HT29 | K562 | MCF7 |
|---|---|---|---|---|---|---|
| pds_cosine | -0.0423 (-0.0240) | 0.0870 (0.0796) | 0.2125 (0.1330) | -0.0403 (-0.0228) | -0.0847 (-0.0534) | 0.0197 (0.0172) |
| mse_capped_norm | 0.0000 (0.0000) | 0.0000 (0.0000) | 0.0000 (0.0000) | 0.0000 (0.0000) | 0.0000 (0.0000) | 0.0000 (0.0000) |
| lfc_nmae | -0.4474 (-0.2093) | -0.6808 (-0.2236) | -0.2048 (-0.1154) | -0.4703 (-0.2074) | -0.5220 (-0.2153) | -0.8860 (-0.2117) |
| fid | 0.0000 (0.0000) | -0.0024 (-0.0014) | 0.0000 (0.0000) | -0.3304 (-0.1185) | -0.1769 (-0.0833) | -0.0004 (-0.0002) |
| reach | -0.4112 (-0.2478) | -0.4120 (-0.3079) | -0.4371 (-0.2512) | -0.3835 (-0.1660) | -0.5127 (-0.3479) | -0.1538 (-0.0984) |
| jac | 0.0000 (0.0000) | -0.0049 (-0.0010) | 0.0000 (0.0000) | 0.5012 (0.0998) | 0.1091 (0.0385) | -0.0004 (-0.0000) |

`avg_score` es la media simple de las seis celdas de `from_replicate`; se verifica en la tabla de intermedios de abajo.


## Contraste contra la prediccion declarada

**Salio donde dijimos.** Los `avg_score` son negativos en 6 de 6 lineas, en el rango -0.197881 (K562) a -0.071562 (HAP1).

Ninguno queda "muy por encima de 0", que era la senal de que el problema estaria en la cadena. La emision que solo copia el control puntua por debajo del comparador oraculo que escala ese mismo control al perfil de respuesta promedio -- que es el orden que la construccion de `baseline.py` predice.

| linea | `avg_score` | ¿negativo? |
|---|---|---|
| A549 | -0.150136 | si |
| BXPC3 | -0.168844 | si |
| HAP1 | -0.071562 | si |
| HT29 | -0.120547 | si |
| K562 | -0.197881 | si |
| MCF7 | -0.170131 | si |

Por metrica, el detalle importa mas que el agregado. De las seis, **tres salen negativas** (`pds_cosine`, `lfc_nmae`, `reach`) y **tres salen exactamente 0.0** -- y los ceros tienen dos causas distintas, que no conviene confundir:

- **`mse_capped_norm` sale 0.0 por CLAMP.** Es la unica de las seis que clampea a `[0, 1]` (`catalog.py:415-416`: `clamp_low=0.0, clamp_high=1.0`). Sin el clamp, A549 daria `(2.562172 - 0.959928) / (0.138661 - 0.959928) = -1.950941`. El comentario del catalogo (`catalog.py:406-410`) explica justamente este caso: "The metric is capped so that 1.0 = predicting the control, so a value past the baseline is already no-skill and the graded Box-Cox tail below it buys nothing a ranking uses." Nuestra emision ES predecir el control, asi que caer al piso del clamp es el comportamiento documentado.
- **`fid` y `jac` salen 0.0 porque `u == b` EXACTAMENTE.** El numerador es cero, no hay clamp. En A549: `fid` da u=0.000000 y b=0.000000; `jac` da u=0.454545 y b=0.454545. Es el hallazgo de `anclas_nativas.md` reapareciendo del lado del candidato: el artefacto de union vacia hace que la emision degenerada y el comparador oraculo den el MISMO numero, asi que `s` sale 0 no porque la emision empate con el baseline por merito, sino porque la metrica no esta mirando la emision. **Un `s = 0` en esas dos celdas no significa "igual al baseline": significa "no medido".**


### ⚠️ Donde la emision degenerada puntua POSITIVO (y no deberia)

El agregado salio negativo en las seis lineas, pero por celda hay **5 de 36 con `s > 0`**: la emision que solo copia el control queda del lado bueno del comparador oraculo. Cada una es un lugar donde la escala premia algo que no es habilidad.

| linea | metrica | u (suelo) | b | r | s | lectura |
|---|---|---|---|---|---|---|
| HT29 | jac | 0.181818 | 0.091090 | 0.272116 | **0.501189** | union vacia: `b` arranca bajo y `u` cobra 1.0 regalados |
| HAP1 | pds_cosine | 0.512311 | 0.437500 | 0.789583 | **0.212480** | `b` cayo por debajo del minimo aleatorio de ~0.5 |
| K562 | jac | 0.242424 | 0.212121 | 0.489996 | **0.109053** | union vacia: `b` arranca bajo y `u` cobra 1.0 regalados |
| BXPC3 | pds_cosine | 0.517992 | 0.476326 | 0.955114 | **0.087025** | `b` cayo por debajo del minimo aleatorio de ~0.5 |
| MCF7 | pds_cosine | 0.512311 | 0.503788 | 0.936174 | **0.019711** | `b` cayo por debajo del minimo aleatorio de ~0.5 |

El caso mas fuerte es **HT29 / jac = +0.5012**: la emision degenerada recorre la mitad de la distancia del baseline al techo. No porque prediga nada, sino porque `de_sig_jaccard` le regala un 1.0 por cada perturbacion sin senal (`metrics/de.py:195`, convencion `J(vacio, vacio) = 1`), y en esa linea el baseline oraculo cobra menos de esos regalos que la emision degenerada.

Es la tercera vez que el mismo artefacto aparece en este proyecto, ahora del lado del candidato: `bucle_techo.md` lo encontro inflando nuestro suelo, `anclas_nativas.md` lo encontro en el `b` oficial, y aqui lo encontramos **convirtiendose en puntaje positivo**. Antes del primer candidato conviene decidir que se hace con `jac`; esta corrida muestra que ignorarlo tiene un costo medible.


## 🔴 Las dos celdas que NO resuelven -- marcadas, no excluidas

`anclas_nativas_v2.md` encontro que el span de `jac` es mas chico que la precision con que se midio su extremo 1. En error estandar de la media de k=5 particiones (`SE = sd / sqrt(5)`), que es la escala correcta porque `r` ES una media:

| linea | metrica | span | sd | SE = sd/sqrt(5) | span/SE | span/sd | estado |
|---|---|---|---|---|---|---|---|
| HAP1 | jac | 0.043914 | 0.084462 | 0.037772 | **1.16** | 0.52 | **no resuelve** |
| A549 | jac | 0.018709 | 0.035008 | 0.015656 | **1.19** | 0.53 | **no resuelve** |

**2 celdas con `span/SE < 2`**: el extremo 1 no esta localizado con precision suficiente para que la distancia al extremo 0 signifique algo. **No se excluyen del calculo** -- el `avg_score` de arriba las incluye tal cual, igual que las incluiria cualquier submission. Se marcan.

Lo que esas celdas aportaron al `s` de esta corrida:

| linea | metrica | u (suelo) | b | r | s reportado | aporte al avg_score |
|---|---|---|---|---|---|---|
| HAP1 | jac | 0.424242 | 0.424242 | 0.468156 | 0.000000 | 0.000000 |
| A549 | jac | 0.454545 | 0.454545 | 0.473254 | 0.000000 | 0.000000 |

En esta corrida el dano es invisible porque `u == b` exactamente y el `s` sale 0: una celda que no resuelve aporto cero. **Con un candidato real que mueva `jac` aunque sea un poco, esas dos celdas aportaran un sexto del `avg_score` cada una, dividido por un span que el propio anchor no localiza.** Es el riesgo concreto que esta prueba del cable no podia revelar y que conviene tener escrito antes del primer candidato.


## Los intermedios completos: las 36 celdas, de `u` al `s`

Todo lo que entra en cada `s`, sin paso oculto: `u` del `agg_results.csv` del run, `b` del `baseline_agg.csv` del bundle, `r` del `anchor_agg.parquet`, el `(u-b)/(r-b)` calculado a mano y el `s` que reporto `score`. La columna `dif` cierra el circuito; donde no es ~0, la causa es el clamp y se dice cual.

| # | linea | metrica | u | b | r | span | (u-b)/(r-b) a mano | s reportado | dif | clamp | span/SE |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | A549 | pds_cosine | 0.474432 | 0.486742 | 0.777841 | 0.291098 | -0.042290 | -0.042290 | 0.00e+00 | - | 22.42 |
| 2 | A549 | mse_capped_norm | 2.562172 | 0.959928 | 0.138661 | 0.821267 | -1.950941 | 0.000000 | 1.95e+00 | **si, a 0.0** | 29.93 |
| 3 | A549 | lfc_nmae | 1.036018 | 0.856686 | 0.455829 | 0.400857 | -0.447372 | -0.447372 | 1.11e-16 | - | 39.19 |
| 4 | A549 | fid | 0.000000 | 0.000000 | 0.417863 | 0.417863 | 0.000000 | 0.000000 | 0.00e+00 | - | 25.27 |
| 5 | A549 | reach | 0.144540 | 0.314403 | 0.727540 | 0.413137 | -0.411154 | -0.411154 | 0.00e+00 | - | 11.16 |
| 6 | A549 | jac | 0.454545 | 0.454545 | 0.473254 | 0.018709 | 0.000000 | 0.000000 | 0.00e+00 | - | 1.19 |
| 7 | BXPC3 | pds_cosine | 0.517992 | 0.476326 | 0.955114 | 0.478788 | 0.087025 | 0.087025 | 0.00e+00 | - | 63.61 |
| 8 | BXPC3 | mse_capped_norm | 1.332759 | 0.824367 | 0.068734 | 0.755633 | -0.672803 | 0.000000 | 6.73e-01 | **si, a 0.0** | 26.46 |
| 9 | BXPC3 | lfc_nmae | 0.996965 | 0.814780 | 0.547191 | 0.267589 | -0.680838 | -0.680838 | 0.00e+00 | - | 11.76 |
| 10 | BXPC3 | fid | 0.000026 | 0.001444 | 0.599253 | 0.597809 | -0.002371 | -0.002371 | 0.00e+00 | - | 29.19 |
| 11 | BXPC3 | reach | 0.142889 | 0.344686 | 0.834502 | 0.489816 | -0.411985 | -0.411985 | 0.00e+00 | - | 23.81 |
| 12 | BXPC3 | jac | 0.000000 | 0.000978 | 0.200769 | 0.199792 | -0.004893 | -0.004893 | 0.00e+00 | - | 19.62 |
| 13 | HAP1 | pds_cosine | 0.512311 | 0.437500 | 0.789583 | 0.352083 | 0.212480 | 0.212480 | 0.00e+00 | - | 51.62 |
| 14 | HAP1 | mse_capped_norm | 1.868728 | 1.003465 | 0.185635 | 0.817830 | -1.057997 | 0.000000 | 1.06e+00 | **si, a 0.0** | 19.93 |
| 15 | HAP1 | lfc_nmae | 1.000821 | 0.897264 | 0.391590 | 0.505674 | -0.204791 | -0.204791 | 2.78e-17 | - | 27.81 |
| 16 | HAP1 | fid | 0.000000 | 0.000000 | 0.546846 | 0.546846 | 0.000000 | 0.000000 | 0.00e+00 | - | 12.97 |
| 17 | HAP1 | reach | 0.192408 | 0.354552 | 0.725538 | 0.370986 | -0.437061 | -0.437061 | 0.00e+00 | - | 10.91 |
| 18 | HAP1 | jac | 0.424242 | 0.424242 | 0.468156 | 0.043914 | 0.000000 | 0.000000 | 0.00e+00 | - | 1.16 |
| 19 | HT29 | pds_cosine | 0.490530 | 0.501894 | 0.783523 | 0.281629 | -0.040350 | -0.040350 | 0.00e+00 | - | 36.86 |
| 20 | HT29 | mse_capped_norm | 1.931565 | 1.137395 | 0.132822 | 1.004573 | -0.790555 | 0.000000 | 7.91e-01 | **si, a 0.0** | 23.50 |
| 21 | HT29 | lfc_nmae | 1.044314 | 0.864894 | 0.483357 | 0.381536 | -0.470257 | -0.470257 | 1.11e-16 | - | 64.98 |
| 22 | HT29 | fid | 0.009259 | 0.114189 | 0.431793 | 0.317604 | -0.330379 | -0.330379 | 0.00e+00 | - | 10.85 |
| 23 | HT29 | reach | 0.216033 | 0.327624 | 0.618619 | 0.290995 | -0.383483 | -0.383483 | 0.00e+00 | - | 15.88 |
| 24 | HT29 | jac | 0.181818 | 0.091090 | 0.272116 | 0.181026 | 0.501189 | 0.501189 | 0.00e+00 | - | 6.92 |
| 25 | K562 | pds_cosine | 0.495265 | 0.520833 | 0.822727 | 0.301894 | -0.084693 | -0.084693 | 0.00e+00 | - | 38.45 |
| 26 | K562 | mse_capped_norm | 2.913227 | 1.173286 | 0.126447 | 1.046839 | -1.662090 | 0.000000 | 1.66e+00 | **si, a 0.0** | 22.23 |
| 27 | K562 | lfc_nmae | 1.028846 | 0.846558 | 0.497343 | 0.349215 | -0.521996 | -0.521996 | 2.22e-16 | - | 91.31 |
| 28 | K562 | fid | 0.000000 | 0.076923 | 0.511705 | 0.434782 | -0.176923 | -0.176923 | 0.00e+00 | - | 21.71 |
| 29 | K562 | reach | 0.200448 | 0.406814 | 0.809303 | 0.402489 | -0.512725 | -0.512725 | 0.00e+00 | - | 16.37 |
| 30 | K562 | jac | 0.242424 | 0.212121 | 0.489996 | 0.277874 | 0.109053 | 0.109053 | 0.00e+00 | - | 9.15 |
| 31 | MCF7 | pds_cosine | 0.512311 | 0.503788 | 0.936174 | 0.432386 | 0.019711 | 0.019711 | 0.00e+00 | - | 59.16 |
| 32 | MCF7 | mse_capped_norm | 1.774874 | 0.922312 | 0.104958 | 0.817354 | -1.043076 | 0.000000 | 1.04e+00 | **si, a 0.0** | 22.89 |
| 33 | MCF7 | lfc_nmae | 1.056808 | 0.872196 | 0.663824 | 0.208372 | -0.885975 | -0.885975 | 0.00e+00 | - | 40.19 |
| 34 | MCF7 | fid | 0.000000 | 0.000173 | 0.462605 | 0.462432 | -0.000374 | -0.000374 | 0.00e+00 | - | 11.45 |
| 35 | MCF7 | reach | 0.013563 | 0.101892 | 0.676395 | 0.574503 | -0.153750 | -0.153750 | 0.00e+00 | - | 20.21 |
| 36 | MCF7 | jac | 0.000000 | 0.000029 | 0.072243 | 0.072214 | -0.000399 | -0.000399 | 0.00e+00 | - | 6.66 |

**El cable queda verificado aritmeticamente**: en las celdas sin clamp, el `s` que reporta `score` reproduce `(u-b)/(r-b)` calculado a mano desde los tres artefactos independientes (el run, el baseline del bundle y el anchor del bundle). Las unicas celdas con `dif` grande son las seis de `mse_capped_norm`, y la causa es su clamp documentado.

### El `avg_score` recompuesto desde las seis celdas

| linea | media de las 6 celdas | `avg_score` reportado | dif |
|---|---|---|---|
| A549 | -0.150136 | -0.150136 | 0.00e+00 |
| BXPC3 | -0.168844 | -0.168844 | 0.00e+00 |
| HAP1 | -0.071562 | -0.071562 | 0.00e+00 |
| HT29 | -0.120547 | -0.120547 | 0.00e+00 |
| K562 | -0.197881 | -0.197881 | 2.78e-17 |
| MCF7 | -0.170131 | -0.170131 | 0.00e+00 |

Es una media simple de las seis, sin pesos: se recompone exacta.


## Intermedios de corrida

| linea | real_shape | pred_shape | n control real | n control pred | seed_real | seed_pred | run_s | score_s | total_s |
|---|---|---|---|---|---|---|---|---|---|
| A549 | [14579, 34025] | [13600, 34025] | 1379 | 400 | 20260001 | 20261001 | 60.88 | 3.05 | 67.45 |
| BXPC3 | [16428, 34025] | [13600, 34025] | 3228 | 400 | 20260002 | 20261002 | 81.05 | 3.04 | 87.59 |
| HAP1 | [16372, 34025] | [13600, 34025] | 3172 | 400 | 20260003 | 20261003 | 84.4 | 3.04 | 90.88 |
| HT29 | [15936, 34025] | [13600, 34025] | 2736 | 400 | 20260004 | 20261004 | 74.37 | 3.06 | 80.29 |
| K562 | [14598, 34025] | [13600, 34025] | 1398 | 400 | 20260005 | 20261005 | 64.38 | 3.04 | 69.68 |
| MCF7 | [15869, 34025] | [13600, 34025] | 2669 | 400 | 20260006 | 20261006 | 82.42 | 3.05 | 89.01 |

La poblacion es la de siempre: `build_fold_adatas` **sin modificar**, `seed_real = 20260001 + fold_idx` (400 celulas por perturbacion sin reemplazo mas TODO el pool de control) y la emision de `suelo_floor.predict_floor` **sin modificar** con `seed_pred = 20261001 + fold_idx`. La etiqueta de control se renombro a `non-targeting` en LOS DOS lados, porque `io.validate_pair` exige que los conjuntos de etiquetas coincidan entre PRED y REAL; el paso 0 de `anclas_nativas_v2.md` probo que ese renombrado no mueve ningun numero.

Nota de comparabilidad: el `u` de esta corrida es, celda por celda, el mismo suelo que puntuo `bucle_suelo.md` -- misma poblacion, misma emision, mismas semillas. Lo que cambia es contra que se puntua: ahi fue contra si mismo (`avg_score = 0.0` por construccion), aqui contra los dos extremos nativos.

| linea | metrica | u de esta corrida | u de bucle_suelo.md | dif |
|---|---|---|---|---|
| A549 | pds_cosine | 0.474432 | 0.474432 | 0.00e+00 |
| A549 | mse_capped_norm | 2.562172 | 2.562172 | 0.00e+00 |
| A549 | lfc_nmae | 1.036018 | 1.036018 | 0.00e+00 |
| A549 | fid | 0.000000 | 0.000000 | 0.00e+00 |
| A549 | reach | 0.144540 | 0.144540 | 0.00e+00 |
| A549 | jac | 0.454545 | 0.454545 | 0.00e+00 |
| BXPC3 | pds_cosine | 0.517992 | 0.517992 | 0.00e+00 |
| BXPC3 | mse_capped_norm | 1.332759 | 1.332759 | 0.00e+00 |
| BXPC3 | lfc_nmae | 0.996965 | 0.996965 | 1.11e-16 |
| BXPC3 | fid | 0.000026 | 0.000026 | 0.00e+00 |
| BXPC3 | reach | 0.142889 | 0.142889 | 0.00e+00 |
| BXPC3 | jac | 0.000000 | 0.000000 | 0.00e+00 |

## Lo que NO se hizo

- **No se construyo ningun candidato.** Esta emision es desechable y no se guarda como modelo.
- **No se reparo nada.** Nada fallo, asi que no hubo que decidirlo; si algo hubiera fallado, el script para en la primera linea y registra el mensaje literal (`cable_suelo.py`, `sys.exit(6)`).
- **No se tocaron los bundles**, ni `bucle_suelo.md`, ni `bucle_techo.md`, ni ningun arnes. `suelo_arnes.py` y `suelo_floor.py` se reusaron por import, sin modificar.
- **No se excluyo del calculo ninguna celda**, incluidas las dos que no resuelven: se marcan.
- No se leyo nada bajo `~/vcc/benchmark/` ni `/datos/benchmark_VEDADO_techo_metricas/`.
- Sin envios, sin `vcc sample`, sin commits.

