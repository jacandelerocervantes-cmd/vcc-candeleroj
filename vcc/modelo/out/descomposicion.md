# De que esta hecho cada numero: la descomposicion por linea y metrica

Las seis metricas de `vcc2026` **no se miden sobre la misma cohorte**, y una de ellas asigna un valor por convencion a la ausencia de senal. El `avg_score` promedia plano a 1/6 seis numeros calculados sobre poblaciones de 8 a 33 perturbaciones. Esto descompone cada numero en lo que mide y lo que es convencion.

Lectura y conteo sobre artefactos que ya existian, mas dos recomputos que se declaran donde aparecen: la dispersion de perfiles (`descomp_pds.json`) y la componente genuina de `r` (`descomp_rgen.json`), esta ultima porque el bundle **no persiste** el detalle por perturbacion de las particiones del anchor.


# PARTE 1 -- La descomposicion

## `jac` -- la convencion de union vacia, medida en tres objetos

**Primero la premisa, verificada y no inferida.** El conteo usa `E = count(jaccard == 1.0)` como numero de perturbaciones con union vacia. Un jaccard de 1.0 tambien podria venir de un acuerdo perfecto, asi que se midio: se recorrio una corrida con `--write-degenes` y se contaron las uniones con `techo_arnes.contar_de_por_perturbacion` (sin modificar), que reproduce las dos lineas de filtrado de `de_sig_jaccard`. Contingencia sobre 33 perturbaciones de A549:

| | union vacia | union NO vacia |
|---|---|---|
| `jac == 1.0` | **15** | 0 |
| `jac != 1.0` | 0 | **18** |

Las dos celdas fuera de la diagonal son **0**: la equivalencia `jac == 1.0` <=> `|union| == 0` queda **verificada** (True). Asi que `E` es una medicion.

La convencion es `de_sig_jaccard`'s `J(vacio, vacio) = 1` (`metrics/de.py:195`, y el retorno en `metrics/de.py:263`: `out[p] = 1.0 if union == 0 else tp / union`). Cada perturbacion sin senal en NINGUNO de los dos lados aporta 1.0 al promedio de 33.

| objeto | linea | E (union vacia) | genuinos | media de los genuinos | comp. CONVENCION = E/33 | comp. GENUINA = sum/33 | reportado | suma | dif |
|---|---|---|---|---|---|---|---|---|---|
| nuestro suelo | A549 | 15 | 18 | 0.000000 | 0.454545 | 0.000000 | 0.454545 | 0.454545 | 0.0e+00 |
| nuestro suelo | BXPC3 | 0 | 33 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.0e+00 |
| nuestro suelo | HAP1 | 14 | 19 | 0.000000 | 0.424242 | 0.000000 | 0.424242 | 0.424242 | 0.0e+00 |
| nuestro suelo | HT29 | 6 | 27 | 0.000000 | 0.181818 | 0.000000 | 0.181818 | 0.181818 | 0.0e+00 |
| nuestro suelo | K562 | 8 | 25 | 0.000000 | 0.242424 | 0.000000 | 0.242424 | 0.242424 | 0.0e+00 |
| nuestro suelo | MCF7 | 0 | 33 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.0e+00 |
| b oficial | A549 | 15 | 18 | 0.000000 | 0.454545 | 0.000000 | 0.454545 | 0.454545 | 0.0e+00 |
| b oficial | BXPC3 | 0 | 33 | 0.000978 | 0.000000 | 0.000978 | 0.000978 | 0.000978 | 0.0e+00 |
| b oficial | HAP1 | 14 | 19 | 0.000000 | 0.424242 | 0.000000 | 0.424242 | 0.424242 | 0.0e+00 |
| b oficial | HT29 | 3 | 30 | 0.000199 | 0.090909 | 0.000181 | 0.091090 | 0.091090 | 1.4e-17 |
| b oficial | K562 | 7 | 26 | 0.000000 | 0.212121 | 0.000000 | 0.212121 | 0.212121 | 0.0e+00 |
| b oficial | MCF7 | 0 | 33 | 0.000029 | 0.000000 | 0.000029 | 0.000029 | 0.000029 | 0.0e+00 |
| **r (anchor, k=5)** | A549 | 13.8 (media) | 19.2 (media) | 0.095526 | 0.418182 | 0.055072 | 0.473254 | 0.473254 | 0.0e+00 |
| **r (anchor, k=5)** | HAP1 | 11.2 (media) | 21.8 (media) | 0.198131 | 0.339394 | 0.128762 | 0.468156 | 0.468156 | 0.0e+00 |

**La identidad `jac_reportado = (E + sum(genuinos)) / 33` se verifica exacta en todos los objetos y todas las lineas** (columna `dif`).

**El hallazgo.** En nuestro suelo la media de los jaccards genuinos es **exactamente 0.000000 en las seis lineas**: el numero reportado es convencion al 100%. En el `b` oficial la componente genuina es 0.000000 en tres lineas y <= 0.001 en las otras tres. Es decir: **para las dos emisiones triviales, `jac` no mide nada de la emision; cuenta cuantas perturbaciones no tienen senal.**

El `r` del anchor SI tiene componente genuina (es un replicado real), pero es una fraccion de su valor crudo:

| linea | r crudo | r convencion | r genuina | genuina / crudo |
|---|---|---|---|---|
| A549 | 0.473254 | 0.418182 | 0.055072 | 11.6% |
| HAP1 | 0.468156 | 0.339394 | 0.128762 | 27.5% |

### El recomputo del `r`, y por que hizo falta

`anchor_splits.parquet` guarda `(split_index, seed, metric, value, n_perturbations)` -- el AGREGADO por particion, que mezcla las dos componentes. La descomposicion de `r` no se puede hacer leyendo. Se recomputaron las 5 particiones de A549 y HAP1 con la MISMA funcion del marcador (`cell_eval2.ceiling._disjoint_halves`, la que usa `anchor.py:117`), las MISMAS semillas que el bundle estampo, `half_b` como prediccion (`anchor.py:120`) y `control_source=pred` (`anchor.py:103`).

**Criterio declarado antes de medir:** A549 y HAP1 son las dos celdas que el encargo nombra por su `span/SE` sobre el crudo. No se eligieron por su valor ni por su signo.

| linea | split | semilla | valor del bundle | valor recomputado | dif | E | media genuinos |
|---|---|---|---|---|---|---|---|
| A549 | 0 | 2968811710 | 0.490868 | 0.490868 | 0.0e+00 | 14 | 0.115717 |
| A549 | 1 | 3677149159 | 0.435372 | 0.435372 | 0.0e+00 | 13 | 0.068364 |
| A549 | 2 | 745650761 | 0.516056 | 0.516056 | 1.1e-16 | 15 | 0.112769 |
| A549 | 3 | 2884920346 | 0.428102 | 0.428102 | 0.0e+00 | 13 | 0.056368 |
| A549 | 4 | 2642120001 | 0.495873 | 0.495873 | 1.1e-16 | 14 | 0.124410 |
| HAP1 | 0 | 2968811710 | 0.524090 | 0.524090 | 0.0e+00 | 13 | 0.214749 |
| HAP1 | 1 | 3677149159 | 0.567021 | 0.567021 | 0.0e+00 | 14 | 0.247984 |
| HAP1 | 2 | 745650761 | 0.407504 | 0.407504 | 5.6e-17 | 9 | 0.185318 |
| HAP1 | 3 | 2884920346 | 0.506686 | 0.506686 | 0.0e+00 | 13 | 0.186031 |
| HAP1 | 4 | 2642120001 | 0.335480 | 0.335480 | 0.0e+00 | 7 | 0.156571 |

**La reproduccion es exacta**: la diferencia maxima contra el valor que el bundle estampo es 1.1e-16 sobre las 10 particiones. Asi que la descomposicion describe el `r` que realmente se usa, no una aproximacion.

Y aparece el mecanismo que explica el `span` raquitico del crudo: **la componente de convencion trabaja EN CONTRA del span.** El replicado real tiene MENOS uniones vacias que el baseline, asi que cobra menos 1.0 regalados, lo que arrastra su valor crudo hacia abajo y casi borra la holgura genuina:

| linea | aporte de la CONVENCION al span | aporte de lo GENUINO | span crudo |
|---|---|---|---|
| A549 | -0.036364 | +0.055072 | +0.018709 |
| HAP1 | -0.084848 | +0.128762 | +0.043914 |


## `fid`, `reach`, `lfc_nmae` -- las gateadas: quien entro, quien se cayo, y por que

Las tres excluyen perturbaciones, pero por mecanismos DISTINTOS, y la diferencia importa:

| metrica | mecanismo | condicion | linea del codigo | que recibe la excluida |
|---|---|---|---|---|
| `lfc_nmae` | la fila **no se emite** | gate vacio, menor que `min_gate_size` (=10), o denominador cero/no finito | `metrics/de.py:608-611` | nada: "those score nothing rather than 1.0" |
| `fid` (`direction_fidelity_yield_raw`) | la fila se emite con **NaN** | `max(n_pred, n_conf) == 0`, es decir los DOS lados vacios | `metrics/direction.py:776` con `_safe_div` en `689-699` | NaN, que el agregador omite |
| `reach` (`direction_reach_raw`) | la fila se emite con **NaN** | `n_conf == 0` (pool de referencia vacio), o `q` indefinido | `metrics/direction.py:1050-1053` | NaN, que el agregador omite |

**Ninguna de las tres aplica una convencion a la ausencia.** El codigo de `reach` lo dice explicitamente (`metrics/direction.py:1048-1050`): "Spec 5: 'never reached P0 -> a computed 0.0' is DOMAIN-RESTRICTED and must not override the matrix. It is trivially true at N_conf = 0 (empty pool) and when q is undefined (no threshold exists) -- **both are NaN, not 0**." Es la decision opuesta a la de `jac`, en el mismo paquete.

La consecuencia es que la media se toma sobre una cohorte mas chica, no que se contamine con valores de relleno. Medido (la columna `media_entraron` reproduce el valor reportado en las 36 celdas):

| objeto | linea | metrica | filas emitidas | NaN | entraron | de 33 | media de las que entraron | reportado | dif |
|---|---|---|---|---|---|---|---|---|---|
| nuestro suelo | A549 | lfc_nmae | 8 | 0 | **8** | 33 | 1.036018 | 1.036018 | 0.0e+00 |
| nuestro suelo | A549 | fid | 33 | 15 | **18** | 33 | 0.000000 | 0.000000 | 0.0e+00 |
| nuestro suelo | A549 | reach | 33 | 15 | **18** | 33 | 0.144540 | 0.144540 | 2.8e-17 |
| nuestro suelo | BXPC3 | lfc_nmae | 23 | 0 | **23** | 33 | 0.996965 | 0.996965 | 1.1e-16 |
| nuestro suelo | BXPC3 | fid | 33 | 0 | **33** | 33 | 0.000026 | 0.000026 | 0.0e+00 |
| nuestro suelo | BXPC3 | reach | 33 | 0 | **33** | 33 | 0.142889 | 0.142889 | 0.0e+00 |
| nuestro suelo | HAP1 | lfc_nmae | 8 | 0 | **8** | 33 | 1.000821 | 1.000821 | 0.0e+00 |
| nuestro suelo | HAP1 | fid | 33 | 14 | **19** | 33 | 0.000000 | 0.000000 | 0.0e+00 |
| nuestro suelo | HAP1 | reach | 33 | 14 | **19** | 33 | 0.192408 | 0.192408 | 2.8e-17 |
| nuestro suelo | HT29 | lfc_nmae | 15 | 0 | **15** | 33 | 1.044314 | 1.044314 | 2.2e-16 |
| nuestro suelo | HT29 | fid | 33 | 6 | **27** | 33 | 0.009259 | 0.009259 | 0.0e+00 |
| nuestro suelo | HT29 | reach | 33 | 6 | **27** | 33 | 0.216033 | 0.216033 | 5.6e-17 |
| nuestro suelo | K562 | lfc_nmae | 9 | 0 | **9** | 33 | 1.028846 | 1.028846 | 0.0e+00 |
| nuestro suelo | K562 | fid | 33 | 8 | **25** | 33 | 0.000000 | 0.000000 | 0.0e+00 |
| nuestro suelo | K562 | reach | 33 | 8 | **25** | 33 | 0.200448 | 0.200448 | 0.0e+00 |
| nuestro suelo | MCF7 | lfc_nmae | 33 | 0 | **33** | 33 | 1.056808 | 1.056808 | 2.2e-16 |
| nuestro suelo | MCF7 | fid | 33 | 0 | **33** | 33 | 0.000000 | 0.000000 | 0.0e+00 |
| nuestro suelo | MCF7 | reach | 33 | 0 | **33** | 33 | 0.013563 | 0.013563 | 0.0e+00 |
| b oficial | A549 | lfc_nmae | 8 | 0 | **8** | 33 | 0.856686 | 0.856686 | 1.1e-16 |
| b oficial | A549 | fid | 33 | 15 | **18** | 33 | 0.000000 | 0.000000 | 0.0e+00 |
| b oficial | A549 | reach | 33 | 15 | **18** | 33 | 0.314403 | 0.314403 | 0.0e+00 |
| b oficial | BXPC3 | lfc_nmae | 23 | 0 | **23** | 33 | 0.814780 | 0.814780 | 0.0e+00 |
| b oficial | BXPC3 | fid | 33 | 0 | **33** | 33 | 0.001444 | 0.001444 | 0.0e+00 |
| b oficial | BXPC3 | reach | 33 | 0 | **33** | 33 | 0.344686 | 0.344686 | 0.0e+00 |
| b oficial | HAP1 | lfc_nmae | 8 | 0 | **8** | 33 | 0.897264 | 0.897264 | 0.0e+00 |
| b oficial | HAP1 | fid | 33 | 14 | **19** | 33 | 0.000000 | 0.000000 | 0.0e+00 |
| b oficial | HAP1 | reach | 33 | 14 | **19** | 33 | 0.354552 | 0.354552 | 0.0e+00 |
| b oficial | HT29 | lfc_nmae | 15 | 0 | **15** | 33 | 0.864894 | 0.864894 | 0.0e+00 |
| b oficial | HT29 | fid | 33 | 3 | **30** | 33 | 0.114189 | 0.114189 | 0.0e+00 |
| b oficial | HT29 | reach | 33 | 6 | **27** | 33 | 0.327624 | 0.327624 | 0.0e+00 |
| b oficial | K562 | lfc_nmae | 9 | 0 | **9** | 33 | 0.846558 | 0.846558 | 1.1e-16 |
| b oficial | K562 | fid | 33 | 7 | **26** | 33 | 0.076923 | 0.076923 | 0.0e+00 |
| b oficial | K562 | reach | 33 | 8 | **25** | 33 | 0.406814 | 0.406814 | 0.0e+00 |
| b oficial | MCF7 | lfc_nmae | 33 | 0 | **33** | 33 | 0.872196 | 0.872196 | 1.1e-16 |
| b oficial | MCF7 | fid | 33 | 0 | **33** | 33 | 0.000173 | 0.000173 | 0.0e+00 |
| b oficial | MCF7 | reach | 33 | 0 | **33** | 33 | 0.101892 | 0.101892 | 2.8e-17 |

**Dos hechos del codigo confirmados por medicion:**

- `metrics/de.py:610-611` dice de los gates de `lfc_nmae`: "All three are real-side decisions, so **the surviving set is the same for every submission**." Medido: la cohorte de `lfc_nmae` es identica entre nuestro suelo y el `b` en las seis lineas (8/8, 23/23, 8/8, 15/15, 9/9, 33/33). Verificado.
- El NaN de `reach` depende SOLO del lado real (`n_conf == 0`), asi que su cohorte tambien es identica entre las dos emisiones (15/15, 0/0, 14/14, 6/6, 8/8, 0/0). El de `fid` depende de `max(n_pred, n_conf)`, que SI mira la emision: en HT29 cae 6 veces en nuestro suelo y 3 en el `b`, y en K562 8 contra 7. **La cohorte de `fid` es la unica de las tres que la submission puede mover.**


## 🔴 `reach` -- la pregunta abierta desde `bucle_suelo.md`, ahora medida

La pregunta era: `reach` da NO-NULO (0.014 a 0.216 en nuestro suelo) con conjunto predicho VACIO. Estaba declarada "no entendida, deliberadamente sin explicar". La respuesta esta en el codigo y **no es una convencion analoga a `J(vacio,vacio)=1`**.

**Que calcula.** `reach = k*(t) / N_conf(t)` (`metrics/direction.py:1054`):

```
out[target] = _to_float_nan(k_star.get(target, 0) / n_conf)
```

donde `k*` es la profundidad hasta la que el ranking de la prediccion se mantiene direccionalmente puro, y `N_conf` es el presupuesto confiable de la REFERENCIA.

**Por que el conjunto predicho vacio no lo anula.** Con `universe='adjudicated'` el ranking recorre **solo los genes significativos de la REFERENCIA** (`metrics/direction.py:980-982`: "`universe='adjudicated'` ranks only reference-significant genes, so `k* <= N_conf` and the result is bounded by 1"). La prediccion no aporta el POOL; aporta el ORDEN y los SIGNOS. Una emision degenerada no llama nada significativo, pero sigue teniendo un `|log2FC|` por gen con el que ordenar, y un signo por gen.

**De donde sale el valor no-nulo.** `metrics/direction.py:1016-1019`, issue #279:

> "`reach_raw`'s no-skill point is `~c/N_conf`, NOT a constant, and for a no-skill submission it is decided almost entirely at the HEAD of a ranking the submission controls. **A first-pair match is SUFFICIENT for `k* >= 1` (purity 1/1 >= P0), one coin flip independent of depth**, and it dominates the no-skill probability."

Y la medicion que el propio codigo trae (`metrics/direction.py:1030-1033`): "Measured at 0.13.0 over 200 replicates per cell: `E[reach_raw]` runs 0.5100 at `N_conf = 1` to 0.0019 at 500 with `P(reach_raw > 0) = 0.385-0.565` at EVERY N_conf, and the fitted `c` is ~0.96 for a coin-flip predictor".

**Veredicto: no es una convencion aplicada a la ausencia, es un piso de azar.** Un acierto de signo en el primer gen del ranking -- una moneda -- da `k* = 1` y por tanto `reach = 1/N_conf > 0`. No hay un valor asignado a la ausencia: hay un computo cuyo valor esperado sin habilidad es `~c/N_conf`, que no es cero. Por eso **no se descompone como `jac`**: todo el numero "mide", lo que pasa es que su punto de no-habilidad no esta en 0.

Dato relevante para la escala: el codigo tambien dice (`metrics/direction.py:1035-1037`) que "on official val A the baseline arm reads 0.042314 against an ESTIMATED no-skill 0.0378, so **the baseline already sits at the no-skill point**". Es decir, en `reach` el `b` medido ES el piso de azar -- que es exactamente lo que un extremo 0 deberia ser.


## `mse_capped_norm` -- no gateada, pero clampeada

**No esta gateada y lee las 33.** En `results.csv` la metrica derivada no aparece por perturbacion (es `agg="ratio_of_sums"` con `DerivedAgg`, `catalog.py:412-414`), asi que su cohorte se lee del artefacto del anchor, que la estampa:

| linea | `n_perturbations_min` | `n_perturbations_max` | de 33 |
|---|---|---|---|
| A549 | None | None | 33 |
| BXPC3 | None | None | 33 |
| HAP1 | None | None | 33 |
| HT29 | None | None | 33 |
| K562 | None | None | 33 |
| MCF7 | None | None | 33 |

Sus dos componentes SI aparecen por perturbacion y con 33 filas cada una (`expr_mse_unbiased_capped` y `expr_distance_unbiased`), sin NaN, en las seis lineas y en los dos objetos. Confirmado.

**🔴 El primitivo que importa, sin cifras por linea.** `catalog.py:415-416` le pone `clamp_low=0.0, clamp_high=1.0`:

```
scoring=Scoring(scored=True, direction="lower", anchor=0.0,
                penalty="boxcox", clamp_low=0.0, clamp_high=1.0),
```

y el comentario (`catalog.py:406-410`) declara que es la unica: "this is the ONE scored vcc2026 member that clamps to [0, 1], on BOTH scales. The metric is capped so that 1.0 = predicting the control, so a value past the baseline is already no-skill and the graded Box-Cox tail below it buys nothing a ranking uses."

**La propiedad, declarada:** por debajo del extremo 0 esta metrica devuelve **exactamente 0**, asi que **no informa si un candidato esta mejorando mientras siga del lado malo del baseline**. `cable_suelo.md` lo midio: las seis celdas dieron 0.000000, y sin el clamp A549 habria dado -1.950941. Un candidato que pase de -1.95 a -0.2 ve la misma cifra que uno que no cambio nada. **No se calcula aqui cuanto habria que bajar para salir del clamp: eso seria un peso, no una descripcion.**


# PARTE 2 -- `pds_cosine`: por que el `b` queda por DEBAJO de `u`

⚠️ **No se compara contra el 0.5 de `scales.py:428`.** Esa es la base constante de otra escala (`low-random_high-1_v10`), con otro estimando y otra poblacion. El 0.5 que aparece abajo es otro objeto: el valor que el propio estadistico de rango asigna a una prediccion sin informacion, con su cita.

## Que calcula `pds_cosine`, y que premia

`metrics/discrimination.py:48-53`: "A rank-based retrieval score: for each non-control perturbation, rank how close the predicted signed delta (effect) is to the matching real delta versus all other real deltas. Score = `1 - rank / D`". Opera sobre el **delta** (perfil menos control), no sobre el perfil.

Lo que premia es **especificidad**: que el delta predicho para `p` se parezca mas al delta real de `p` que a los de las otras 32. No premia magnitud ni dispersion.

Y el punto sin informacion esta fijado por el estadistico, no por una constante externa (`metrics/discrimination.py:74-77`): "A FULLY constant prediction self-corrects under either policy (the ranks are a bijection onto `{0..n-1}`, so **the panel mean is 0.5**)". Mas `discrimination.py:70`: con `tie_policy="midrank"` (el de `vcc2026`) "`midrank` scores every such target exactly 0.5" para una fila empatada.


## La hipotesis del encargo, puesta a prueba

La hipotesis era: el escalado de `baseline.py:344` mueve todas las perturbaciones hacia el MISMO perfil, lo que las acerca entre si y destruye discriminacion. Es verificable: se mide la dispersion entre los deltas de perturbacion en las tres emisiones, con las funciones del propio marcador (`prep.pseudobulk_bulk_lognorm` con `bulk_target_sum=50000` y `prep.delta`, que es lo que ve el comparador de `vcc2026`).

**Criterio declarado antes de medir: las seis lineas, sin seleccion.**

| linea | objeto | norma media del delta | distancia cosine media entre pares | min | max | rango |
|---|---|---|---|---|---|---|
| A549 | real (la verdad) | 11.8535 | **0.784571** | 0.3254 | 1.0012 | 0.6758 |
| A549 | nuestro suelo | 14.0450 | **0.480861** | 0.4026 | 0.5641 | 0.1616 |
| A549 | b oficial | 11.4206 | **0.748589** | 0.6460 | 0.8232 | 0.1772 |
| BXPC3 | real (la verdad) | 10.3084 | **0.865822** | 0.3120 | 1.0676 | 0.7556 |
| BXPC3 | nuestro suelo | 11.8508 | **0.498443** | 0.4550 | 0.5452 | 0.0902 |
| BXPC3 | b oficial | 9.3891 | **0.804155** | 0.7154 | 0.8777 | 0.1623 |
| HAP1 | real (la verdad) | 8.8310 | **0.885348** | 0.3876 | 1.0319 | 0.6443 |
| HAP1 | nuestro suelo | 10.8103 | **0.526843** | 0.4670 | 0.5947 | 0.1277 |
| HAP1 | b oficial | 8.5245 | **0.852660** | 0.7856 | 0.9175 | 0.1319 |
| HT29 | real (la verdad) | 9.3300 | **0.856785** | 0.3563 | 0.9584 | 0.6021 |
| HT29 | nuestro suelo | 11.6900 | **0.492605** | 0.4446 | 0.5483 | 0.1037 |
| HT29 | b oficial | 9.1193 | **0.830077** | 0.7619 | 0.8981 | 0.1362 |
| K562 | real (la verdad) | 10.2526 | **0.772358** | 0.4561 | 0.9110 | 0.4549 |
| K562 | nuestro suelo | 12.3826 | **0.481366** | 0.3913 | 0.5367 | 0.1453 |
| K562 | b oficial | 10.0912 | **0.752123** | 0.6571 | 0.8286 | 0.1715 |
| MCF7 | real (la verdad) | 7.9252 | **0.788039** | 0.5023 | 0.9793 | 0.4770 |
| MCF7 | nuestro suelo | 9.5723 | **0.502960** | 0.4325 | 0.5614 | 0.1288 |
| MCF7 | b oficial | 7.7583 | **0.754239** | 0.6470 | 0.8367 | 0.1897 |

🔴 **La hipotesis queda REFUTADA por la medicion.** Los deltas del `b` oficial estan casi tan dispersos como los reales: distancia cosine media entre pares 0.749-0.853 contra 0.772-0.885 del real. El escalado **no** colapsa las perturbaciones hacia un punto. Y nuestro suelo esta MENOS disperso (0.481-0.527), no mas.

Las normas tampoco distinguen: los tres objetos tienen deltas de norma comparable (7.7 a 14.0), y **ninguno tiene deltas de norma cero** (`frac_norma_cero = 0.000` en los 18 casos), asi que el mecanismo de "delta de norma cero empata toda la fila" que `discrimination.py:66-70` describe no es el que opera aqui.

Lo que SI distingue al real es el **rango**: sus pares van de 0.31 a 1.07 (hay perturbaciones muy parecidas y pares anticorrelados), mientras el `b` se queda en 0.65-0.92 y nuestro suelo en 0.39-0.60. El real tiene ESTRUCTURA; los otros dos tienen dispersion sin estructura.


## Entonces, ¿que explica el orden `b < u`?

**Dispersion no es discriminacion.** `pds_cosine` no pregunta cuanto se separan los deltas predichos entre si, sino si el delta predicho de `p` identifica a `p`. Ninguna de las dos emisiones lleva informacion especifica de la perturbacion: la de nuestro suelo es ruido de remuestreo del control, y la del `b` es un perfil CONSTANTE unico (`baseline.py:121-122`: "the result is still ONE constant vector") mas ruido de remuestreo amplificado por el escalado. Las dos caen, por construccion del estadistico de rango, en el punto sin informacion.

Medido, los dos objetos rodean ese punto:

| linea | `u` (nuestro suelo) | `b` oficial | b - u | ¿de que lado? |
|---|---|---|---|---|
| A549 | 0.474432 | 0.486742 | +0.012311 | b por encima |
| BXPC3 | 0.517992 | 0.476326 | -0.041667 | b por DEBAJO |
| HAP1 | 0.512311 | 0.437500 | -0.074811 | b por DEBAJO |
| HT29 | 0.490530 | 0.501894 | +0.011364 | b por encima |
| K562 | 0.495265 | 0.520833 | +0.025568 | b por encima |
| MCF7 | 0.512311 | 0.503788 | -0.008523 | b por DEBAJO |

`b` queda por debajo de `u` en **3 de 6** lineas y por encima en 3. Los dos objetos se mueven en una banda de 0.4375 a 0.5208, y `|b - u|` va de 0.0085 a 0.0748. **El signo cambia de linea a linea sin direccion consistente, que es como se ve el ruido alrededor de un punto fijo, no un efecto sistematico del escalado.**

Escala de ruido disponible: el anchor da 5 particiones de `pds_cosine` por linea, y su `sd` va de 0.015 a 0.029 (`anclas_nativas_v2.md`). Cuatro de los seis `|b - u|` caen por debajo de 0.029; HAP1 (0.0748) y BXPC3 (0.0417) lo exceden. Asi que el ruido explica la mayoria de los casos pero **no todos**, y para HAP1 queda una causa sin identificar. **No esta en el codigo** ni en estas mediciones cual es.

**Lo que importa para la escala, y esta medido:** en `pds_cosine` el `b` oficial NO acota la trivialidad con margen -- se sienta sobre el punto sin informacion con +-0.04 de ruido, y `cable_suelo.md` midio que una copia del control cobra `s` positivo ahi en 3 de 6 lineas (HAP1 +0.2125, BXPC3 +0.0870, MCF7 +0.0197). Con `r - b` entre 0.28 y 0.48, ese ruido de +-0.04 vale hasta ~0.14 de `s` sin ninguna habilidad.


# ENTREGABLE -- las 36 filas

Para el objeto que un candidato va a tener que superar: el **`b` oficial**, que es el extremo 0 de la escala. `parte que MIDE` es la componente calculada del contenido de la emision; `parte que es CONVENCION` es un valor asignado donde la cantidad no esta definida y que ENTRA al promedio. Para las gateadas la convencion es 0 por construccion (la fila se omite, no se rellena), y lo que se encoge es la cohorte.

| linea | metrica | valor reportado | parte que MIDE | parte que es CONVENCION | perts que entraron | de 33 | fraccion del numero que mide algo |
|---|---|---|---|---|---|---|---|
| A549 | pds_cosine | 0.486742 | 0.486742 | 0.000000 | 33 | 33 | 100.0% |
| A549 | mse_capped_norm | 0.959928 | 0.959928 | 0.000000 | 33 | 33 | 100.0% |
| A549 | lfc_nmae | 0.856686 | 0.856686 | 0.000000 | 8 | 33 | 100.0% |
| A549 | fid | 0.000000 | 0.000000 | 0.000000 | 18 | 33 | 100.0% |
| A549 | reach | 0.314403 | 0.314403 | 0.000000 | 18 | 33 | 100.0% |
| A549 | jac | 0.454545 | 0.000000 | 0.454545 | 18 | 33 | 0.0% ⚠️ |
| BXPC3 | pds_cosine | 0.476326 | 0.476326 | 0.000000 | 33 | 33 | 100.0% |
| BXPC3 | mse_capped_norm | 0.824367 | 0.824367 | 0.000000 | 33 | 33 | 100.0% |
| BXPC3 | lfc_nmae | 0.814780 | 0.814780 | 0.000000 | 23 | 33 | 100.0% |
| BXPC3 | fid | 0.001444 | 0.001444 | 0.000000 | 33 | 33 | 100.0% |
| BXPC3 | reach | 0.344686 | 0.344686 | 0.000000 | 33 | 33 | 100.0% |
| BXPC3 | jac | 0.000978 | 0.000978 | 0.000000 | 33 | 33 | 100.0% |
| HAP1 | pds_cosine | 0.437500 | 0.437500 | 0.000000 | 33 | 33 | 100.0% |
| HAP1 | mse_capped_norm | 1.003465 | 1.003465 | 0.000000 | 33 | 33 | 100.0% |
| HAP1 | lfc_nmae | 0.897264 | 0.897264 | 0.000000 | 8 | 33 | 100.0% |
| HAP1 | fid | 0.000000 | 0.000000 | 0.000000 | 19 | 33 | 100.0% |
| HAP1 | reach | 0.354552 | 0.354552 | 0.000000 | 19 | 33 | 100.0% |
| HAP1 | jac | 0.424242 | 0.000000 | 0.424242 | 19 | 33 | 0.0% ⚠️ |
| HT29 | pds_cosine | 0.501894 | 0.501894 | 0.000000 | 33 | 33 | 100.0% |
| HT29 | mse_capped_norm | 1.137395 | 1.137395 | 0.000000 | 33 | 33 | 100.0% |
| HT29 | lfc_nmae | 0.864894 | 0.864894 | 0.000000 | 15 | 33 | 100.0% |
| HT29 | fid | 0.114189 | 0.114189 | 0.000000 | 30 | 33 | 100.0% |
| HT29 | reach | 0.327624 | 0.327624 | 0.000000 | 27 | 33 | 100.0% |
| HT29 | jac | 0.091090 | 0.000181 | 0.090909 | 30 | 33 | 0.2% ⚠️ |
| K562 | pds_cosine | 0.520833 | 0.520833 | 0.000000 | 33 | 33 | 100.0% |
| K562 | mse_capped_norm | 1.173286 | 1.173286 | 0.000000 | 33 | 33 | 100.0% |
| K562 | lfc_nmae | 0.846558 | 0.846558 | 0.000000 | 9 | 33 | 100.0% |
| K562 | fid | 0.076923 | 0.076923 | 0.000000 | 26 | 33 | 100.0% |
| K562 | reach | 0.406814 | 0.406814 | 0.000000 | 25 | 33 | 100.0% |
| K562 | jac | 0.212121 | 0.000000 | 0.212121 | 26 | 33 | 0.0% ⚠️ |
| MCF7 | pds_cosine | 0.503788 | 0.503788 | 0.000000 | 33 | 33 | 100.0% |
| MCF7 | mse_capped_norm | 0.922312 | 0.922312 | 0.000000 | 33 | 33 | 100.0% |
| MCF7 | lfc_nmae | 0.872196 | 0.872196 | 0.000000 | 33 | 33 | 100.0% |
| MCF7 | fid | 0.000173 | 0.000173 | 0.000000 | 33 | 33 | 100.0% |
| MCF7 | reach | 0.101892 | 0.101892 | 0.000000 | 33 | 33 | 100.0% |
| MCF7 | jac | 0.000029 | 0.000029 | 0.000000 | 33 | 33 | 100.0% |

`mse_capped_norm` es una metrica DERIVADA (`agg="ratio_of_sums"`, `catalog.py:412-414`): no emite fila por perturbacion, asi que no tiene conteo propio. Sus dos componentes si lo tienen, y las dos leen 33 filas sin NaN en las seis lineas, por eso figura con 33.

**4 de 36 celdas miden menos de la mitad de lo que reportan, y todas son `jac`.** En 3 de ellas (A549, HAP1, K562) el `b` oficial mide **exactamente 0.0%** de su propio valor: las uniones vacias son el numero entero. Las `jac` de BXPC3 y MCF7 miden ~100%, pero de un valor de 0.001 y 0.00003 -- miden algo, y ese algo es casi nada.

Y la cohorte, que es el otro eje: `lfc_nmae` entra con 8 de 33 perturbaciones en A549 y HAP1 (24%), `fid` con 18-33, `reach` con 18-33, y `pds`/`mse` con 33. **El `avg_score` promedia esas seis a 1/6 cada una, sin pesar por cohorte** -- de modo que en A549 una perturbacion de las 8 que `lfc_nmae` ve pesa 4.1 veces mas en el agregado que una de las 33 que `pds` ve.


# 🔴 SEGUNDA ESCALA -- DIAGNOSTICA, **no es la cifra de competencia**

**Lo que sigue NO es el puntaje.** El `Overall` usa el crudo, y el crudo es lo que esta en `cable_suelo.md`. Esta escala calcula `(u - b) / (r - b)` **solo sobre la componente genuina** de `jac`, y sirve para una cosa: saber CUAL de las dos partes se movio cuando un candidato cambia. No sustituye a nada y no debe reportarse como resultado.

`SE = sd / sqrt(5)` porque `r` ES una media de 5 particiones (`anchor.py:285`), asi que la precision de su ubicacion es el error estandar de esa media, no su desviacion.

| linea | escala | u | b | r | span | sd de r | SE de r | span/SE | ¿resuelve (>= 3)? |
|---|---|---|---|---|---|---|---|---|---|
| A549 | **cruda (la de competencia)** | 0.454545 | 0.454545 | 0.473254 | 0.018709 | 0.035008 | 0.015656 | **1.19** | **NO** |
| A549 | diagnostica (solo genuina) | 0.000000 | 0.000000 | 0.055072 | 0.055072 | 0.014645 | 0.006549 | **8.41** | **si** |
| HAP1 | **cruda (la de competencia)** | 0.424242 | 0.424242 | 0.468156 | 0.043914 | 0.084462 | 0.037772 | **1.16** | **NO** |
| HAP1 | diagnostica (solo genuina) | 0.000000 | 0.000000 | 0.128762 | 0.128762 | 0.010201 | 0.004562 | **28.23** | **si** |

**La respuesta a la pregunta del encargo: si, resuelven.** A549 pasa de `span/SE = 1.19` sobre el crudo a **8.41** sobre la componente genuina; HAP1 de `1.16` a **28.23**. Las dos cruzan el umbral de 3 que no cruzaban.

El mecanismo esta en las dos columnas de en medio. Sobre el crudo, el span es chico porque la convencion resta (el replicado real tiene menos uniones vacias que el baseline) y el `sd` es grande porque `E` varia entre particiones -- en A549 de 13 a 15, y cada unidad de `E` mueve el crudo 1/33 = 0.0303. Sobre la componente genuina el span es mas grande Y el `sd` mas chico, las dos cosas a la vez.

⚠️ Tres limites de esta escala diagnostica, para que no se use de mas: (1) `u` y `b` genuinos valen **exactamente 0** en estas dos lineas, asi que el `s` diagnostico de una emision trivial es 0 igual que el crudo -- lo que cambia es la REGLA, no el puntaje de esta corrida; (2) solo esta calculada para A549 y HAP1, porque recomputar las particiones de las otras cuatro cuesta otros ~20 min y no se pidio; (3) la componente genuina de `u` y `b` se cuenta sobre las perturbaciones con union no vacia de CADA objeto, y esos conjuntos no son identicos entre objetos -- con `u` y `b` genuinos en 0 eso no afecta a estos numeros, pero con un candidato real habria que fijar una cohorte comun antes de comparar.


# 20 filas trabajadas

Veinte filas con todos los intermedios de su descomposicion: nada entre el artefacto y la columna `fraccion que mide`. Se eligieron recorriendo las lineas en orden y las metricas en orden, tomando las primeras 20 celdas del entregable -- **criterio posicional, declarado, no depende del valor**.

| # | linea | metrica | reportado | mide | convencion | E | genuinos | filas emitidas | NaN | entraron | cohorte | fraccion que mide | mecanismo |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | A549 | pds_cosine | 0.486742 | 0.486742 | 0.000000 | - | - | 33 | 0 | 33 | 100% | 100.0% | sin exclusion; punto sin info = 0.5 |
| 2 | A549 | mse_capped_norm | 0.959928 | 0.959928 | 0.000000 | - | - | 0 | 0 | 33 | 100% | 100.0% | sin exclusion; clamp a [0,1] |
| 3 | A549 | lfc_nmae | 0.856686 | 0.856686 | 0.000000 | - | - | 8 | 0 | 8 | 24% | 100.0% | fila no emitida (gate) |
| 4 | A549 | fid | 0.000000 | 0.000000 | 0.000000 | - | - | 33 | 15 | 18 | 55% | 100.0% | NaN si max(n_pred,n_conf)=0 |
| 5 | A549 | reach | 0.314403 | 0.314403 | 0.000000 | - | - | 33 | 15 | 18 | 55% | 100.0% | NaN si n_conf=0; piso de azar ~c/N_conf |
| 6 | A549 | jac | 0.454545 | 0.000000 | 0.454545 | 15 | 18 | 33 | 0 | 18 | 55% | 0.0% | convencion J(vacio,vacio)=1 |
| 7 | BXPC3 | pds_cosine | 0.476326 | 0.476326 | 0.000000 | - | - | 33 | 0 | 33 | 100% | 100.0% | sin exclusion; punto sin info = 0.5 |
| 8 | BXPC3 | mse_capped_norm | 0.824367 | 0.824367 | 0.000000 | - | - | 0 | 0 | 33 | 100% | 100.0% | sin exclusion; clamp a [0,1] |
| 9 | BXPC3 | lfc_nmae | 0.814780 | 0.814780 | 0.000000 | - | - | 23 | 0 | 23 | 70% | 100.0% | fila no emitida (gate) |
| 10 | BXPC3 | fid | 0.001444 | 0.001444 | 0.000000 | - | - | 33 | 0 | 33 | 100% | 100.0% | NaN si max(n_pred,n_conf)=0 |
| 11 | BXPC3 | reach | 0.344686 | 0.344686 | 0.000000 | - | - | 33 | 0 | 33 | 100% | 100.0% | NaN si n_conf=0; piso de azar ~c/N_conf |
| 12 | BXPC3 | jac | 0.000978 | 0.000978 | 0.000000 | 0 | 33 | 33 | 0 | 33 | 100% | 100.0% | convencion J(vacio,vacio)=1 |
| 13 | HAP1 | pds_cosine | 0.437500 | 0.437500 | 0.000000 | - | - | 33 | 0 | 33 | 100% | 100.0% | sin exclusion; punto sin info = 0.5 |
| 14 | HAP1 | mse_capped_norm | 1.003465 | 1.003465 | 0.000000 | - | - | 0 | 0 | 33 | 100% | 100.0% | sin exclusion; clamp a [0,1] |
| 15 | HAP1 | lfc_nmae | 0.897264 | 0.897264 | 0.000000 | - | - | 8 | 0 | 8 | 24% | 100.0% | fila no emitida (gate) |
| 16 | HAP1 | fid | 0.000000 | 0.000000 | 0.000000 | - | - | 33 | 14 | 19 | 58% | 100.0% | NaN si max(n_pred,n_conf)=0 |
| 17 | HAP1 | reach | 0.354552 | 0.354552 | 0.000000 | - | - | 33 | 14 | 19 | 58% | 100.0% | NaN si n_conf=0; piso de azar ~c/N_conf |
| 18 | HAP1 | jac | 0.424242 | 0.000000 | 0.424242 | 14 | 19 | 33 | 0 | 19 | 58% | 0.0% | convencion J(vacio,vacio)=1 |
| 19 | HT29 | pds_cosine | 0.501894 | 0.501894 | 0.000000 | - | - | 33 | 0 | 33 | 100% | 100.0% | sin exclusion; punto sin info = 0.5 |
| 20 | HT29 | mse_capped_norm | 1.137395 | 1.137395 | 0.000000 | - | - | 0 | 0 | 33 | 100% | 100.0% | sin exclusion; clamp a [0,1] |

# Lo que NO se hizo

- **No se modifico ninguna metrica ni ningun arnes.** `_disjoint_halves`, `contar_de_por_perturbacion`, `pseudobulk_bulk_lognorm`, `delta` y `build_fold_adatas` se usaron tal cual.
- **No se sustituyo la escala cruda por la diagnostica en ningun lado.** La diagnostica esta marcada arriba de su propia seccion, no en nota al pie, y el `Overall` sigue siendo el crudo de `cable_suelo.md`.
- **No se selecciono ninguna celda por su signo ni por su valor.** Los tres criterios de seleccion que aparecen se declararon antes de medir: las seis lineas para la dispersion; A549 y HAP1 para el recomputo de `r`, por ser las que el encargo nombra; y las primeras 20 celdas en orden posicional para las filas trabajadas.
- **No se derivaron objetivos numericos por linea** para ningun candidato. En particular no se calculo cuanto habria que bajar `mse_capped_norm` para salir del clamp.
- **No se comparo contra el 0.5 de `scales.py:428`**, que pertenece a otra escala. El 0.5 citado es el punto sin informacion del estadistico de rango (`discrimination.py:70`, `74-77`).
- No se leyo nada bajo `~/vcc/benchmark/` ni `/datos/benchmark_VEDADO_techo_metricas/`.
- Sin envios, sin `vcc sample`, sin commits.

