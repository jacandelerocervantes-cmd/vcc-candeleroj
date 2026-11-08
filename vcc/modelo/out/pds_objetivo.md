# `pds_cosine` como objetivo: la escala contra el azar y de que esta hecha

## 🔴 0. Que es este reporte, que no es, y la advertencia que no se puede omitir

**Presenta datos.** Declara la construccion antes de medir, da las tablas con
sus controles y reporta los defectos de la propia construccion. **No lleva
veredicto, ni conclusiones, ni nada dirigido al siguiente candidato.**

🔴 **LOS CUATRO BRAZOS DIAGNOSTICOS NO SON CANDIDATOS.** Los tres oraculos
(§3.1, §3.2, §3.3) leen el **LADO REAL PERTURBADO de la MISMA poblacion que se
puntua**: son **TRANSDUCTIVOS** y van al mismo estante que el `b` oficial, que
tambien lo es. El cuarto (§3.4) no es transductivo pero es igual de diagnostico.
**SON UNA MEDICION DE LA METRICA, NUNCA UNA EMISION QUE PUEDA ENVIARSE.**
Ninguna de sus cifras puede entrar jamas en una submission ni en un rasgo de un
candidato, ni como valor, ni como umbral, ni como punto de partida. Lo que miden
es **la metrica**, no un modelo.

🔴 **LA PREGUNTA.** `pds_cosine` es la unica de las seis que pregunta *"¿es el
blanco correcto entre 33?"*. `control_permutado.md` §4.1 midio que en `fid` el
brazo permutado cobra el **99.19 %-100.06 %** del brazo real, con `u` entre
0.4303 y 0.4459 sobre un azar de 0.50: ahi la metrica paga **emitir**, no
acertar. De `pds_cosine` ya estaban medidos su punto sin informacion y su `b`.
**Lo que faltaba es donde cae `r`**, y de que esta hecho el numero.

---

## 1. El punto sin informacion, DERIVADO DEL CODIGO

### 1.1 Que estadistico es, sobre cuantos competidores

| que | valor | de donde (archivo:linea) |
|---|---|---|
| estadistico | **rango 0-based del propio delta entre los `n` deltas reales**, convertido con `score = 1 - rango/D` | `metrics/discrimination.py:259` y `:304` (`out[str(p)] = 1.0 if D <= 0 else 1.0 - rank / D`) |
| competidores | `n` = **numero total de perturbaciones no control del lado REAL** = **33** | `metrics/discrimination.py:216` (`n = real_keys.size`, con el comentario "denominator uses the FULL real count") |
| denominador | `D = n - 1` = **32** en el preset de competencia | `metrics/discrimination.py:217` (`D = n if rank_denominator == "n" else n - 1`) + `config.py:248-249` (preset `v2`: `"distance": "cosine", "rank_denominator": "n-1", "tie_policy": "midrank"`) |
| empates | **midrank**: `rango = n_menores + (n_iguales - 1)/2` | `metrics/discrimination.py:332` |
| distancia | **cosine**, invariante de escala | `config.py:248`, `catalog.py:319-320` |
| exclusion | `exclude_target_gene=True`, `exclusion_scope="panel"`: **los 33 genes blanco del panel salen del espacio de rasgos una sola vez**, antes de calcular cualquier distancia | `config.py:100`, `config.py:118`, `metrics/discrimination.py:250-259` |
| agregado del panel | **media** sobre las perturbaciones | `catalog.py:321` (`agg="mean"`) |
| politica de puntuacion | `BOUNDED_UNFLOORED` (sin recorte en 0), `anchor = 1.0` | `catalog.py:320`, `scoring.py:203` |

### 1.2 La esperanza sin habilidad: `0.5` EXACTO, y por que no depende de `n`

Sin habilidad el rango del propio delta es **uniforme sobre `{0, 1, ..., n-1}`**,
asi que `E[rango] = (n-1)/2 = 16`, y

```
E[score] = 1 - E[rango]/D = 1 - ((n-1)/2)/(n-1) = 1/2   EXACTO, para todo n
```

El `0.5` no es una aproximacion ni depende de que `n` sea 33: sale de que
`D = n-1` es exactamente el rango maximo posible. El propio catalogo lo dice en
`catalog.py:313-314`: *"Its comparator sits at the metric's chance point (a
uniform rank is 0.5)"*.

**Comprobado numericamente llamando a la funcion del marcador** (no
reimplementandola): `cell_eval2.metrics.discrimination.discrimination_score`
con `n = 33`, `distance="cosine"`, `rank_denominator="n-1"`,
`tie_policy="midrank"`, `control_source="real"`:

| comprobacion | resultado |
|---|---|
| aritmetica: `1 - 16/32` | **0.50000000000000000** |
| prediccion **constante** (todas las filas iguales) | **0.50000000000000000** exacto -- los rangos son una biyeccion sobre `{0..n-1}`, documentado en `metrics/discrimination.py:72-74` |
| prediccion **independiente** del lado real, 200 paneles de 33 | **0.499527**, `sd` entre paneles 0.049967, `SE` **0.003533** (0.13 `SE` del 0.5) |
| prediccion = lado real (oraculo) | **1.00000000000000000** exacto |
| granularidad de un valor por perturbacion | `1/D` = **0.031250** |

🔴 **El punto sin informacion es `0.5`, derivado y comprobado, no supuesto.**

---

## 2. De donde sale cada cifra ya medida

| cantidad | de donde | n |
|---|---|---|
| `b` oficial | `cell-eval2 baseline`, fila `mean` de `baseline_agg.csv` (`score.py:731`), via `anclas_nativas_v2.md` | 1 por linea |
| `r` oficial | columna `replicate` de `anchor_agg.parquet`, media de 5 particiones (`anchor.py:285`), via `anclas_nativas_v2.md` | 5 particiones |
| `sd_b` | `dispersion_b{,_p2,_p3}.json`, 8 realizaciones de `b` sobre la misma poblacion real | 8 por linea |
| `sd_u` del SUELO | `dispersion_u.json`, 8 semillas de emision del control remuestreado | 8 por linea |
| `r_sd` | `anclas_nativas_v2.md`, `sd` entre las 5 particiones de `split_half_raw` | 5 por linea |
| SUELO / gwps REAL / gwps PERMUTADO | `dispersion_u.json`, `gwps_via.json`, `control_permutado.json`, **ya medidos** | 8 por linea cada uno |
| `pds` por perturbacion de esos tres brazos | `results.csv` de las 144 corridas conservadas en `/datos/modelo/out/{dispu,gwps_via,gwps_perm}` | 33 filas por corrida |

🔴 **No se cruzan las tres escalas del paquete.** Todo `s` de este reporte es
`from_replicate`, que es la columna que toma el `avg_score`. No se mezcla con
`from_baseline` ni con la replica #276 ni con las constantes del paquete.

---

## 3. Los cuatro brazos diagnosticos: la construccion, declarada antes de medir

### 3.0 Lo que NO cambia en ninguno de los cuatro

| | |
|---|---|
| arnes | `suelo_arnes.build_fold_adatas` **SIN MODIFICAR**, contrato `predict_fn(control_pool_X, labels, rng)` |
| muestreo entero | `cand01_covar._reescala_entero` **SIN MODIFICAR** (Binomial para `r<=1`, Poisson aditivo para `r>1`; esperanza `x*r` y enteros en las dos ramas) |
| poblacion real | la **MISMA** de siempre, `seed_real = 20260001 + fold`, comprobada por sha256 en cada corrida |
| bundle | el **MISMO** de competencia, `/datos/modelo/out/nativas_comp/BUNDLE_<linea>` |
| preset / profile | `vcc2026` / `vcc2026`, `config_hash` verificado contra `competition_payload()` en cada corrida |
| escala de `s` | `from_replicate`, una sola |

### 3.1 ORACULO COMPLETO, dos realizaciones

🔴 **Transductivo**: lee las medias del lado real perturbado.

| | |
|---|---|
| **1a, por COPIA** | el bloque predicho de `h` **ES** el bloque real perturbado de `h` (las mismas 400 celulas). Reproduce el delta real **exactamente**, asi que es la **cota superior** y la **comprobacion de cordura**: `pds_cosine` tiene que dar **1.0**. Semilla `20273900 + fold` (la emision no consume el rng para las 33; el bloque de control predicho si se remuestrea) |
| **1b, por RAZON** | `r_g = media_real_pert_g / media_control_g` en **todos** los genes, realizado con el mismo `_reescala_entero` de los brazos 2 y 3. Es el **denominador comparable** de la fraccion de cola: 1b y los brazos 2/3 se diferencian **solo en que genes** llevan razon, no en como se realiza. Semilla `20275000 + 100*fold + j`, `j = 0..3` |

La distancia entre 1a y 1b mide lo que cuesta **la realizacion por conteos
enteros**, separado de lo que cuesta **recortar genes**.

### 3.2 ORACULO DE COLA

La **misma razon** del brazo 1b, pero **solo en los genes del conjunto ED real**
de cada perturbacion; `r = 1` en todo el resto. Semilla
`20276000 + 100*fold + j`, `j = 0..3`.

🔴 **El conjunto ED real se obtiene con la maquinaria del propio marcador**, no
con una reimplementacion: una **sonda** por linea corre `cell-eval2 run` con
`--write-degenes` sobre la poblacion real (que depende solo de `seed_real`, asi
que una sonda por linea alcanza) y de ahi se saca `R(h)` con
`cell_eval2.de.prepare_de` + el filtro `p_adj < 0.05` + `_exclude_own_gene`,
exactamente las mismas tres lineas que usa `techo_arnes.contar_de_por_perturbacion`.

🔴 **En las perturbaciones de ED real VACIO este brazo no emite NADA** (`r = 1`
en todos los genes): ahi es el SUELO por construccion. Son 15, 0, 14, 6, 8 y 0
perturbaciones de 33 (`control_permutado.md` §7), asi que la media del panel de
este brazo **esta diluida hacia el azar por construccion**. Por eso la fraccion
se da **dos veces**: sobre las 33 (la TABLA B3) y **solo donde la cola existe** (la TABLA B3b).

### 3.3 ORACULO DE BULTO

El **complemento exacto** del 3.2: la misma razon con los genes del ED real
puestos a `r = 1`. Semilla `20277000 + 100*fold + j`, `j = 0..3`.

🔴 Los dos brazos **particionan los genes**, no el valor de la metrica:
`pds_cosine` no es aditivo sobre genes, asi que **no hay ninguna razon para que
cola + bulto de la metrica sumen el completo**, y la TABLA B3 reporta la suma medida
precisamente para mostrar cuanto se aparta de 1.

### 3.4 SIGNOS ALEATORIOS DENSOS

La **misma densidad** y el **mismo multiconjunto de `|s|`** que la via gwps, con
el **signo de cada gen sorteado**:

```
r_gwps = 1 + s        (la regla de escala de tres_vias_rango.md §2.2, sin tocar)
r_signos = 1 + |s| * sigma,   sigma = ±1 iid con p = 1/2,   recortado a [0, 2]
```

| | |
|---|---|
| semilla de emision | `20273000 + 100*fold + j`, `j = 0..7` (**la familia que el encargo fija**) |
| semilla del sorteo de signos | `20274000 + j`, **la misma en las seis lineas** para que las lineas sean comparables |
| los 6 replegados | `DTX3L`, `H1F0`, `HERC5`, `HERC6`, `MAP3K14`, `ZBP1` siguen en `r = 1`, igual que en gwps REAL y PERMUTADO |

Es el brazo que en `fid` mostro que la metrica paga **emitir**; aqui dice si
`pds_cosine` hace lo mismo.

### 3.5 Las semillas, y que ninguna colisiona

🔴 Despues del defecto §10.8 de `control_permutado.md` (donde la familia de
emision choco con la del desarreglo en MCF7), aqui la ausencia de colision **se
comprueba en codigo antes de la primera corrida** y se guarda en
`/datos/modelo/out/pds_obj/semillas.json`:

| familia | valores | choques |
|---|---|---|
| copia | `20273900 + fold` | - |
| completo | `20275000 + 100*fold + j`, `j=0..3` | - |
| cola | `20276000 + 100*fold + j`, `j=0..3` | - |
| bulto | `20277000 + 100*fold + j`, `j=0..3` | - |
| signos (emision) | `20273000 + 100*fold + j`, `j=0..7` | - |
| signos (sorteo) | `20274000 + j`, `j=0..7` | - |

**134 semillas distintas, 0 choques internos, 0 choques con las familias
2026000x-2026600x, 2027100x, 2027200x y 2027250x ya gastadas.** La semilla de la
**poblacion real** (`20260001 + fold`) **si** se reusa a proposito -- tiene que
ser la misma que en el suelo, gwps REAL y gwps PERMUTADO para que las tablas
sean comparables -- y eso se verifica por sha256 en cada corrida, no por semilla.

---

## 4. PARTE A: la escala de `pds_cosine` contra el azar

### TABLA A1: la escala de `pds_cosine` contra el azar

| linea | azar (derivado) | `b` oficial | `b` - azar | `r` oficial | `r` - azar | span `r-b` | `sd_b` (8 real.) | `sd_u` suelo (8 sem.) | `r_sd` (5 part.) | `resolucion_8` b vs u |
|---|---|---|---|---|---|---|---|---|---|---|
| A549 | 0.5000 | 0.486742 | **-0.013258** | 0.777841 | **+0.277841** | 0.291099 | 0.031003 | 0.082684 | 0.029037 | 0.031220 |
| BXPC3 | 0.5000 | 0.476326 | **-0.023674** | 0.955114 | **+0.455114** | 0.478788 | 0.022594 | 0.042707 | 0.016829 | 0.017082 |
| HAP1 | 0.5000 | 0.437500 | **-0.062500** | 0.789583 | **+0.289583** | 0.352083 | 0.034716 | 0.053724 | 0.015251 | 0.022615 |
| HT29 | 0.5000 | 0.501894 | **+0.001894** | 0.783523 | **+0.283523** | 0.281629 | 0.034730 | 0.050099 | 0.017083 | 0.021552 |
| K562 | 0.5000 | 0.520833 | **+0.020833** | 0.822727 | **+0.322727** | 0.301894 | 0.033133 | 0.044209 | 0.017558 | 0.019533 |
| MCF7 | 0.5000 | 0.503788 | **+0.003788** | 0.936174 | **+0.436174** | 0.432386 | 0.032763 | 0.054783 | 0.016343 | 0.022568 |

**Conteos.** `r` < azar en **0/6**. `b` < azar en **3/6** (A549, BXPC3, HAP1).
`r` - azar va de **+0.277841** a **+0.455114**; `b` - azar de **-0.062500** a **+0.020833**.

**El `b` contra el azar con su propia resolucion** (`SE_b = sd_b/sqrt(8)`, 8 realizaciones de `b`):

| linea | `b` de 8 real. | `sd_b` | `SE_b` | `b` - azar | \|`b`-azar\| > `SE_b` | `u` suelo de 8 sem. | `u` - azar | `SE_u` | \|`u`-azar\| > `SE_u` |
|---|---|---|---|---|---|---|---|---|---|
| A549 | 0.481179 | 0.031003 | 0.010961 | -0.018821 | **si** | 0.503433 | +0.003433 | 0.029233 | no |
| BXPC3 | 0.522609 | 0.022594 | 0.007988 | +0.022609 | **si** | 0.527225 | +0.027225 | 0.015099 | **si** |
| HAP1 | 0.501657 | 0.034716 | 0.012274 | +0.001657 | no | 0.537405 | +0.037405 | 0.018994 | **si** |
| HT29 | 0.510653 | 0.034730 | 0.012279 | +0.010653 | no | 0.504380 | +0.004380 | 0.017713 | no |
| K562 | 0.497869 | 0.033133 | 0.011714 | -0.002131 | no | 0.510417 | +0.010417 | 0.015630 | no |
| MCF7 | 0.491596 | 0.032763 | 0.011583 | -0.008404 | no | 0.526397 | +0.026397 | 0.019369 | **si** |

### 🔴 TABLA A2: donde cae el `b` DEL BUNDLE dentro de sus 8 realizaciones

El `b` que define la escala `s` es **una sola realizacion**, la del bundle. `dispersion_b` midio 8 realizaciones de `b` sobre la misma poblacion real. Esta tabla sitúa la del bundle dentro de esa distribucion. 🔴 Es la cantidad que explica §6.1: donde el `b` del bundle es un valor bajo de su propia distribucion, el suelo cobra `s > 0` sin ninguna habilidad.

| linea | `b` del bundle | media de los 8 `b` | `sd_b` | **`z` del `b` del bundle** | `b` bundle - azar | `s` del SUELO (8 sem.) |
|---|---|---|---|---|---|---|
| A549 | 0.486742 | 0.481179 | 0.031003 | **+0.18** | -0.013258 | +0.0573 |
| BXPC3 | 0.476326 | 0.522609 | 0.022594 | **-2.05** | -0.023674 | +0.1063 |
| HAP1 | 0.437500 | 0.501657 | 0.034716 | **-1.85** | -0.062500 | +0.2838 |
| HT29 | 0.501894 | 0.510653 | 0.034730 | **-0.25** | +0.001894 | +0.0088 |
| K562 | 0.520833 | 0.497869 | 0.033133 | **+0.69** | +0.020833 | -0.0345 |
| MCF7 | 0.503788 | 0.491596 | 0.032763 | **+0.37** | +0.003788 | +0.0523 |

---

## 5. PARTE B: de que esta hecho `pds_cosine`, y si se alcanza desde la cola

### TABLA B0: de que esta hecho el delta real que lee `pds_cosine`

Medido en el MISMO espacio de la metrica (`bulk_lognorm`, `bulk_target_sum = 50000.0` de `run_params.yaml`) y con los 33 genes blanco del panel FUERA del espacio de rasgos, que es lo que hace `exclusion_scope="panel"`. `pds_cosine` usa distancia **cosine**, asi que la cantidad que reparte es la **norma al cuadrado** del delta.

| linea | genes en el espacio | del panel fuera | perts con ED vacio | **frac. de \|\|delta\|\|² en la COLA (media)** | [min, max] | frac. en el BULTO | `cos`(completo, cola) medio | `cos`(completo, bulto) medio |
|---|---|---|---|---|---|---|---|---|
| A549 | 33992 | 33 | 15 | **0.0348** | [0.0000, 0.4069] | 0.9652 | 0.1686 | 0.9810 |
| BXPC3 | 33992 | 33 | 0 | **0.0848** | [0.0062, 0.5592] | 0.9152 | 0.2402 | 0.9537 |
| HAP1 | 33992 | 33 | 14 | **0.0530** | [0.0000, 0.4275] | 0.9470 | 0.2088 | 0.9706 |
| HT29 | 33992 | 33 | 6 | **0.0676** | [0.0000, 0.5052] | 0.9324 | 0.1904 | 0.9623 |
| K562 | 33992 | 33 | 8 | **0.0381** | [0.0000, 0.3481] | 0.9619 | 0.1533 | 0.9798 |
| MCF7 | 33992 | 33 | 0 | **0.1279** | [0.0067, 0.4359] | 0.8721 | 0.3299 | 0.9319 |

**La cola lleva entre 3.48% y 12.79% de la norma al cuadrado del delta real; el bulto el resto.**

### TABLA B1: `pds_cosine` crudo (`u`) por brazo y por linea

| linea | **copia** | **completo** | **cola** | **bulto** | **signos** | SUELO | gwps REAL | gwps PERM |
|---|---|---|---|---|---|---|---|---|
| A549 | **1.0000** | **1.0000** ± 0.0000 | **0.6399** ± 0.0308 | **1.0000** ± 0.0000 | **0.4838** ± 0.0286 | 0.5034 ± 0.0827 | 0.5240 ± 0.0281 | 0.4722 ± 0.0266 |
| BXPC3 | **1.0000** | **1.0000** ± 0.0000 | **0.8198** ± 0.0360 | **1.0000** ± 0.0000 | **0.5123** ± 0.0312 | 0.5272 ± 0.0427 | 0.4906 ± 0.0157 | 0.4912 ± 0.0404 |
| HAP1 | **1.0000** | **1.0000** ± 0.0000 | **0.6520** ± 0.0486 | **1.0000** ± 0.0000 | **0.5099** ± 0.0185 | 0.5374 ± 0.0537 | 0.4734 ± 0.0289 | 0.5026 ± 0.0514 |
| HT29 | **1.0000** | **1.0000** ± 0.0000 | **0.7038** ± 0.0234 | **1.0000** ± 0.0000 | **0.4738** ± 0.0302 | 0.5044 ± 0.0501 | 0.4903 ± 0.0193 | 0.4957 ± 0.0415 |
| K562 | **1.0000** | **1.0000** ± 0.0000 | **0.6113** ± 0.0643 | **1.0000** ± 0.0000 | **0.5374** ± 0.0587 | 0.5104 ± 0.0442 | 0.5200 ± 0.0222 | 0.4995 ± 0.0481 |
| MCF7 | **1.0000** | **1.0000** ± 0.0000 | **0.9605** ± 0.0095 | **1.0000** ± 0.0000 | **0.4967** ± 0.0413 | 0.5264 ± 0.0548 | 0.5107 ± 0.0274 | 0.5098 ± 0.0472 |

### TABLA B2: `s` de `pds_cosine` (`from_replicate`) por brazo y por linea

| linea | **copia** | **completo** | **cola** | **bulto** | **signos** | SUELO | gwps REAL | gwps PERM |
|---|---|---|---|---|---|---|---|---|
| A549 | **+1.7632** | **+1.7632** ± 0.0000 | **+0.5262** ± 0.1060 | **+1.7632** ± 0.0000 | **-0.0102** ± 0.0983 | +0.0573 ± 0.2840 | +0.1281 ± 0.0965 | -0.0500 ± 0.0915 |
| BXPC3 | **+1.0938** | **+1.0938** ± 0.0000 | **+0.7175** ± 0.0753 | **+1.0938** ± 0.0000 | **+0.0752** ± 0.0651 | +0.1063 ± 0.0892 | +0.0299 ± 0.0327 | +0.0312 ± 0.0843 |
| HAP1 | **+1.5976** | **+1.5976** ± 0.0000 | **+0.6092** ± 0.1382 | **+1.5976** ± 0.0000 | **+0.2058** ± 0.0526 | +0.2838 ± 0.1526 | +0.1019 ± 0.0821 | +0.1849 ± 0.1460 |
| HT29 | **+1.7687** | **+1.7687** ± 0.0000 | **+0.7170** ± 0.0832 | **+1.7687** ± 0.0000 | **-0.0996** ± 0.1074 | +0.0088 ± 0.1779 | -0.0412 ± 0.0685 | -0.0219 ± 0.1475 |
| K562 | **+1.5872** | **+1.5872** ± 0.0000 | **+0.2996** ± 0.2131 | **+1.5872** ± 0.0000 | **+0.0549** ± 0.1946 | -0.0345 ± 0.1464 | -0.0027 ± 0.0734 | -0.0706 ± 0.1592 |
| MCF7 | **+1.1476** | **+1.1476** ± 0.0000 | **+1.0562** ± 0.0220 | **+1.1476** ± 0.0000 | **-0.0164** ± 0.0955 | +0.0523 ± 0.1267 | +0.0159 ± 0.0634 | +0.0140 ± 0.1092 |

### 🔴 TABLA B3: la fraccion que alcanza el ORACULO DE COLA

Tres escalas, declaradas. La **anclada en el azar** es la que usa el punto sin informacion derivado en §1 (`0.5`), y es la que responde "¿se puede mover `pds_cosine` tocando solo la cola?".

| linea | `u` completo (razon) | `u` cola | `u` bulto | `u` azar | **(cola-azar)/(completo-azar)** | **(bulto-azar)/(completo-azar)** | suma cola+bulto | `s` cola / `s` completo |
|---|---|---|---|---|---|---|---|---|
| A549 | 1.0000 | 0.6399 | 1.0000 | 0.5000 | **+0.2798** | **+1.0000** | +1.2798 | +0.2984 |
| BXPC3 | 1.0000 | 0.8198 | 1.0000 | 0.5000 | **+0.6397** | **+1.0000** | +1.6397 | +0.6560 |
| HAP1 | 1.0000 | 0.6520 | 1.0000 | 0.5000 | **+0.3040** | **+1.0000** | +1.3040 | +0.3813 |
| HT29 | 1.0000 | 0.7038 | 1.0000 | 0.5000 | **+0.4077** | **+1.0000** | +1.4077 | +0.4054 |
| K562 | 1.0000 | 0.6113 | 1.0000 | 0.5000 | **+0.2225** | **+1.0000** | +1.2225 | +0.1887 |
| MCF7 | 1.0000 | 0.9605 | 1.0000 | 0.5000 | **+0.9209** | **+1.0000** | +1.9209 | +0.9203 |

**Rango de la fraccion de cola anclada en el azar sobre las seis lineas: [+0.2225, +0.9209].**

### 🔴 TABLA B3b: la misma fraccion SOLO sobre las perturbaciones con ED real NO VACIO

En las perturbaciones de ED real vacio el ORACULO DE COLA no emite NADA (`r = 1` en todos los genes), asi que ahi es el suelo por construccion y arrastra la media de la TABLA B3 hacia el azar. Aqui la misma fraccion contada solo donde la cola existe.

| linea | n con ED lleno | `pds` completo (ED lleno) | `pds` cola (ED lleno) | `pds` bulto (ED lleno) | **(cola-azar)/(completo-azar)** | `pds` cola (ED vacio) | `pds` suelo (ED vacio) |
|---|---|---|---|---|---|---|---|
| A549 | 18 | 1.0000 | 0.7613 | 1.0000 | **+0.5226** | 0.4943 | 0.5156 |
| BXPC3 | 33 | 1.0000 | 0.8198 | 1.0000 | **+0.6397** | n/d | n/d |
| HAP1 | 19 | 1.0000 | 0.7627 | 1.0000 | **+0.5255** | 0.5017 | 0.5689 |
| HT29 | 27 | 1.0000 | 0.7610 | 1.0000 | **+0.5220** | 0.4466 | 0.4525 |
| K562 | 25 | 1.0000 | 0.6541 | 1.0000 | **+0.3081** | 0.4775 | 0.4302 |
| MCF7 | 33 | 1.0000 | 0.9605 | 1.0000 | **+0.9209** | n/d | n/d |

**Rango de la fraccion de cola restringida sobre las seis lineas: [+0.3081, +0.9209].**

### TABLA B4: el rango `[min, max]` de `pds_cosine` por brazo

| brazo | `u` [min, max] | `s` [min, max] | n lineas | n corridas |
|---|---|---|---|---|
| ORACULO COMPLETO (copia) | [1.0000, 1.0000] | [+1.0938, +1.7687] | 6 | 6 |
| ORACULO COMPLETO (razon) | [1.0000, 1.0000] | [+1.0938, +1.7687] | 6 | 24 |
| ORACULO DE COLA | [0.6113, 0.9605] | [+0.2996, +1.0562] | 6 | 24 |
| ORACULO DE BULTO | [1.0000, 1.0000] | [+1.0938, +1.7687] | 6 | 24 |
| SIGNOS ALEATORIOS DENSOS | [0.4738, 0.5374] | [-0.0996, +0.2058] | 6 | 48 |
| SUELO | [0.5034, 0.5374] | [-0.0345, +0.2838] | 6 | 48 |
| REAL | [0.4734, 0.5240] | [-0.0412, +0.1281] | 6 | 48 |
| PERMUTADO | [0.4722, 0.5098] | [-0.0706, +0.1849] | 6 | 48 |

### TABLA B5: cada brazo contra el azar, con su resolucion

| linea | brazo | `u` | `sd` | n | `u` - azar | `SE = sd/sqrt(n)` | \|`u`-azar\| > `SE` |
|---|---|---|---|---|---|---|---|
| A549 | ORACULO COMPLETO (copia) | 1.0000 | 0.0000 | 1 | +0.5000 | n/d | n/d (n=1) |
| A549 | ORACULO COMPLETO (razon) | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| A549 | ORACULO DE COLA | 0.6399 | 0.0308 | 4 | +0.1399 | 0.0154 | **si** |
| A549 | ORACULO DE BULTO | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| A549 | SIGNOS ALEATORIOS DENSOS | 0.4838 | 0.0286 | 8 | -0.0162 | 0.0101 | **si** |
| A549 | SUELO | 0.5034 | 0.0827 | 8 | +0.0034 | 0.0292 | no |
| A549 | REAL | 0.5240 | 0.0281 | 8 | +0.0240 | 0.0099 | **si** |
| A549 | PERMUTADO | 0.4722 | 0.0266 | 8 | -0.0278 | 0.0094 | **si** |
| BXPC3 | ORACULO COMPLETO (copia) | 1.0000 | 0.0000 | 1 | +0.5000 | n/d | n/d (n=1) |
| BXPC3 | ORACULO COMPLETO (razon) | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| BXPC3 | ORACULO DE COLA | 0.8198 | 0.0360 | 4 | +0.3198 | 0.0180 | **si** |
| BXPC3 | ORACULO DE BULTO | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| BXPC3 | SIGNOS ALEATORIOS DENSOS | 0.5123 | 0.0312 | 8 | +0.0123 | 0.0110 | **si** |
| BXPC3 | SUELO | 0.5272 | 0.0427 | 8 | +0.0272 | 0.0151 | **si** |
| BXPC3 | REAL | 0.4906 | 0.0157 | 8 | -0.0094 | 0.0055 | **si** |
| BXPC3 | PERMUTADO | 0.4912 | 0.0404 | 8 | -0.0088 | 0.0143 | no |
| HAP1 | ORACULO COMPLETO (copia) | 1.0000 | 0.0000 | 1 | +0.5000 | n/d | n/d (n=1) |
| HAP1 | ORACULO COMPLETO (razon) | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| HAP1 | ORACULO DE COLA | 0.6520 | 0.0486 | 4 | +0.1520 | 0.0243 | **si** |
| HAP1 | ORACULO DE BULTO | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| HAP1 | SIGNOS ALEATORIOS DENSOS | 0.5099 | 0.0185 | 8 | +0.0099 | 0.0065 | **si** |
| HAP1 | SUELO | 0.5374 | 0.0537 | 8 | +0.0374 | 0.0190 | **si** |
| HAP1 | REAL | 0.4734 | 0.0289 | 8 | -0.0266 | 0.0102 | **si** |
| HAP1 | PERMUTADO | 0.5026 | 0.0514 | 8 | +0.0026 | 0.0182 | no |
| HT29 | ORACULO COMPLETO (copia) | 1.0000 | 0.0000 | 1 | +0.5000 | n/d | n/d (n=1) |
| HT29 | ORACULO COMPLETO (razon) | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| HT29 | ORACULO DE COLA | 0.7038 | 0.0234 | 4 | +0.2038 | 0.0117 | **si** |
| HT29 | ORACULO DE BULTO | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| HT29 | SIGNOS ALEATORIOS DENSOS | 0.4738 | 0.0302 | 8 | -0.0262 | 0.0107 | **si** |
| HT29 | SUELO | 0.5044 | 0.0501 | 8 | +0.0044 | 0.0177 | no |
| HT29 | REAL | 0.4903 | 0.0193 | 8 | -0.0097 | 0.0068 | **si** |
| HT29 | PERMUTADO | 0.4957 | 0.0415 | 8 | -0.0043 | 0.0147 | no |
| K562 | ORACULO COMPLETO (copia) | 1.0000 | 0.0000 | 1 | +0.5000 | n/d | n/d (n=1) |
| K562 | ORACULO COMPLETO (razon) | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| K562 | ORACULO DE COLA | 0.6113 | 0.0643 | 4 | +0.1113 | 0.0322 | **si** |
| K562 | ORACULO DE BULTO | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| K562 | SIGNOS ALEATORIOS DENSOS | 0.5374 | 0.0587 | 8 | +0.0374 | 0.0208 | **si** |
| K562 | SUELO | 0.5104 | 0.0442 | 8 | +0.0104 | 0.0156 | no |
| K562 | REAL | 0.5200 | 0.0222 | 8 | +0.0200 | 0.0078 | **si** |
| K562 | PERMUTADO | 0.4995 | 0.0481 | 8 | -0.0005 | 0.0170 | no |
| MCF7 | ORACULO COMPLETO (copia) | 1.0000 | 0.0000 | 1 | +0.5000 | n/d | n/d (n=1) |
| MCF7 | ORACULO COMPLETO (razon) | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| MCF7 | ORACULO DE COLA | 0.9605 | 0.0095 | 4 | +0.4605 | 0.0048 | **si** |
| MCF7 | ORACULO DE BULTO | 1.0000 | 0.0000 | 4 | +0.5000 | 0.0000 | **si** |
| MCF7 | SIGNOS ALEATORIOS DENSOS | 0.4967 | 0.0413 | 8 | -0.0033 | 0.0146 | no |
| MCF7 | SUELO | 0.5264 | 0.0548 | 8 | +0.0264 | 0.0194 | **si** |
| MCF7 | REAL | 0.5107 | 0.0274 | 8 | +0.0107 | 0.0097 | **si** |
| MCF7 | PERMUTADO | 0.5098 | 0.0472 | 8 | +0.0098 | 0.0167 | no |

---

## 6. PARTE C: las dos cuentas ya medidas, recomprobadas

### TABLA C1: la copia del control contra el `b`, con las 8 semillas

| linea | `s` de 1 semilla (`cable_suelo.md`) | `s` de 8 semillas | `sd` | `SE = sd/sqrt(8)` | `s` > 0 con 8 | \|`s`\| > `SE` | se sostiene |
|---|---|---|---|---|---|---|---|
| A549 | -0.0423 | **+0.0573** | 0.2840 | 0.1004 | **si** | no | no |
| BXPC3 | +0.0870 | **+0.1063** | 0.0892 | 0.0315 | **si** | **si** | **si** |
| HAP1 | +0.2125 | **+0.2838** | 0.1526 | 0.0539 | **si** | **si** | **si** |
| HT29 | -0.0403 | **+0.0088** | 0.1779 | 0.0629 | **si** | no | no |
| K562 | -0.0847 | **-0.0345** | 0.1464 | 0.0518 | no | no | no |
| MCF7 | +0.0197 | **+0.0523** | 0.1267 | 0.0448 | **si** | **si** | **si** |

`s > 0` en **5/6** lineas con 8 semillas; de esas, **3/6** con \|`s`\| por encima del `SE` de la media de 8.

### TABLA C2: con que puntuan las perturbaciones de ED real VACIO

| linea | n ED real vacio | brazo | `pds` medio de los VACIOS | sd | `pds` medio de los LLENOS | sd | dif vacios-llenos | `u` total |
|---|---|---|---|---|---|---|---|---|
| A549 | 15 | SUELO | **0.5156** | 0.0707 | 0.4933 | 0.1080 | +0.0224 | 0.5034 |
| A549 | 15 | gwps REAL | **0.5104** | 0.0267 | 0.5436 | 0.0479 | -0.0332 | 0.5285 |
| A549 | 15 | gwps PERM | **0.5289** | 0.0386 | 0.4249 | 0.0368 | +0.1040 | 0.4722 |
| A549 | 15 | ORACULO COMPLETO (copia) | **1.0000** | - | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| A549 | 15 | ORACULO COMPLETO (razon) | **1.0000** | 0.0000 | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| A549 | 15 | ORACULO DE COLA | **0.4943** | 0.0259 | 0.7613 | 0.0625 | -0.2670 | 0.6399 |
| A549 | 15 | ORACULO DE BULTO | **1.0000** | 0.0000 | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| A549 | 15 | SIGNOS ALEATORIOS DENSOS | **0.5526** | 0.0607 | 0.4264 | 0.0469 | +0.1262 | 0.4838 |
| BXPC3 | 0 | SUELO | **n/d** | - | 0.5272 | 0.0427 | n/d | 0.5272 |
| BXPC3 | 0 | gwps REAL | **n/d** | - | 0.4961 | 0.0157 | n/d | 0.4961 |
| BXPC3 | 0 | gwps PERM | **n/d** | - | 0.4912 | 0.0404 | n/d | 0.4912 |
| BXPC3 | 0 | ORACULO COMPLETO (copia) | **n/d** | - | 1.0000 | 0.0000 | n/d | 1.0000 |
| BXPC3 | 0 | ORACULO COMPLETO (razon) | **n/d** | - | 1.0000 | 0.0000 | n/d | 1.0000 |
| BXPC3 | 0 | ORACULO DE COLA | **n/d** | - | 0.8198 | 0.0360 | n/d | 0.8198 |
| BXPC3 | 0 | ORACULO DE BULTO | **n/d** | - | 1.0000 | 0.0000 | n/d | 1.0000 |
| BXPC3 | 0 | SIGNOS ALEATORIOS DENSOS | **n/d** | - | 0.5123 | 0.0312 | n/d | 0.5123 |
| HAP1 | 14 | SUELO | **0.5689** | 0.0678 | 0.5142 | 0.1004 | +0.0547 | 0.5374 |
| HAP1 | 14 | gwps REAL | **0.5776** | 0.0273 | 0.3945 | 0.0406 | +0.1830 | 0.4722 |
| HAP1 | 14 | gwps PERM | **0.5762** | 0.0831 | 0.4484 | 0.0817 | +0.1278 | 0.5026 |
| HAP1 | 14 | ORACULO COMPLETO (copia) | **1.0000** | - | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| HAP1 | 14 | ORACULO COMPLETO (razon) | **1.0000** | 0.0000 | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| HAP1 | 14 | ORACULO DE COLA | **0.5017** | 0.0692 | 0.7627 | 0.0660 | -0.2611 | 0.6520 |
| HAP1 | 14 | ORACULO DE BULTO | **1.0000** | 0.0000 | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| HAP1 | 14 | SIGNOS ALEATORIOS DENSOS | **0.6180** | 0.0556 | 0.4303 | 0.0572 | +0.1877 | 0.5099 |
| HT29 | 6 | SUELO | **0.4525** | 0.1307 | 0.5159 | 0.0469 | -0.0634 | 0.5044 |
| HT29 | 6 | gwps REAL | **0.4655** | 0.0476 | 0.5025 | 0.0334 | -0.0370 | 0.4957 |
| HT29 | 6 | gwps PERM | **0.5169** | 0.1208 | 0.4910 | 0.0506 | +0.0259 | 0.4957 |
| HT29 | 6 | ORACULO COMPLETO (copia) | **1.0000** | - | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| HT29 | 6 | ORACULO COMPLETO (razon) | **1.0000** | 0.0000 | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| HT29 | 6 | ORACULO DE COLA | **0.4466** | 0.0942 | 0.7610 | 0.0221 | -0.3144 | 0.7038 |
| HT29 | 6 | ORACULO DE BULTO | **1.0000** | 0.0000 | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| HT29 | 6 | SIGNOS ALEATORIOS DENSOS | **0.4603** | 0.0614 | 0.4769 | 0.0267 | -0.0166 | 0.4738 |
| K562 | 8 | SUELO | **0.4302** | 0.1242 | 0.5361 | 0.0633 | -0.1059 | 0.5104 |
| K562 | 8 | gwps REAL | **0.4189** | 0.0399 | 0.5625 | 0.0273 | -0.1436 | 0.5277 |
| K562 | 8 | gwps PERM | **0.4751** | 0.0528 | 0.5073 | 0.0507 | -0.0322 | 0.4995 |
| K562 | 8 | ORACULO COMPLETO (copia) | **1.0000** | - | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| K562 | 8 | ORACULO COMPLETO (razon) | **1.0000** | 0.0000 | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| K562 | 8 | ORACULO DE COLA | **0.4775** | 0.0805 | 0.6541 | 0.0619 | -0.1765 | 0.6113 |
| K562 | 8 | ORACULO DE BULTO | **1.0000** | 0.0000 | 1.0000 | 0.0000 | +0.0000 | 1.0000 |
| K562 | 8 | SIGNOS ALEATORIOS DENSOS | **0.4497** | 0.0749 | 0.5655 | 0.0795 | -0.1158 | 0.5374 |
| MCF7 | 0 | SUELO | **n/d** | - | 0.5264 | 0.0548 | n/d | 0.5264 |
| MCF7 | 0 | gwps REAL | **n/d** | - | 0.5043 | 0.0314 | n/d | 0.5043 |
| MCF7 | 0 | gwps PERM | **n/d** | - | 0.5098 | 0.0472 | n/d | 0.5098 |
| MCF7 | 0 | ORACULO COMPLETO (copia) | **n/d** | - | 1.0000 | 0.0000 | n/d | 1.0000 |
| MCF7 | 0 | ORACULO COMPLETO (razon) | **n/d** | - | 1.0000 | 0.0000 | n/d | 1.0000 |
| MCF7 | 0 | ORACULO DE COLA | **n/d** | - | 0.9605 | 0.0095 | n/d | 0.9605 |
| MCF7 | 0 | ORACULO DE BULTO | **n/d** | - | 1.0000 | 0.0000 | n/d | 1.0000 |
| MCF7 | 0 | SIGNOS ALEATORIOS DENSOS | **n/d** | - | 0.4967 | 0.0413 | n/d | 0.4967 |

---

## 7. Los intermedios crudos, corrida por corrida

### TABLA D: las corridas crudas de los cuatro brazos, una por una

| brazo | linea | `j` | semilla emision | `u` `pds_cosine` | `s` `pds_cosine` | `nnz` pred | max valor | `score_ok` |
|---|---|---|---|---|---|---|---|---|
| ORACULO COMPLETO (copia) | A549 | 0 | 20273900 | **1.000000** | +1.763175 | 29263449 | 2951 | True |
| ORACULO COMPLETO (copia) | BXPC3 | 0 | 20273901 | **1.000000** | +1.093750 | 39342039 | 3589 | True |
| ORACULO COMPLETO (copia) | HAP1 | 0 | 20273902 | **1.000000** | +1.597633 | 44317777 | 2392 | True |
| ORACULO COMPLETO (copia) | HT29 | 0 | 20273903 | **1.000000** | +1.768662 | 37818990 | 4236 | True |
| ORACULO COMPLETO (copia) | K562 | 0 | 20273904 | **1.000000** | +1.587202 | 37135145 | 4888 | True |
| ORACULO COMPLETO (copia) | MCF7 | 0 | 20273905 | **1.000000** | +1.147613 | 50396316 | 4165 | True |
| ORACULO COMPLETO (razon) | A549 | 0 | 20275000 | **1.000000** | +1.763175 | 27485610 | 2054 | True |
| ORACULO COMPLETO (razon) | A549 | 1 | 20275001 | **1.000000** | +1.763175 | 27360404 | 2290 | True |
| ORACULO COMPLETO (razon) | A549 | 2 | 20275002 | **1.000000** | +1.763175 | 27318011 | 2183 | True |
| ORACULO COMPLETO (razon) | A549 | 3 | 20275003 | **1.000000** | +1.763175 | 27331772 | 2280 | True |
| **ORACULO COMPLETO (razon) media** | **A549** | | | **1.000000** (sd 0.000000) | **+1.763175** (sd 0.000000) | | | |
| ORACULO COMPLETO (razon) | BXPC3 | 0 | 20275100 | **1.000000** | +1.093750 | 37385147 | 2465 | True |
| ORACULO COMPLETO (razon) | BXPC3 | 1 | 20275101 | **1.000000** | +1.093750 | 37400998 | 2952 | True |
| ORACULO COMPLETO (razon) | BXPC3 | 2 | 20275102 | **1.000000** | +1.093750 | 37516630 | 2621 | True |
| ORACULO COMPLETO (razon) | BXPC3 | 3 | 20275103 | **1.000000** | +1.093750 | 37242502 | 2974 | True |
| **ORACULO COMPLETO (razon) media** | **BXPC3** | | | **1.000000** (sd 0.000000) | **+1.093750** (sd 0.000000) | | | |
| ORACULO COMPLETO (razon) | HAP1 | 0 | 20275200 | **1.000000** | +1.597633 | 42195318 | 2876 | True |
| ORACULO COMPLETO (razon) | HAP1 | 1 | 20275201 | **1.000000** | +1.597633 | 42058326 | 3064 | True |
| ORACULO COMPLETO (razon) | HAP1 | 2 | 20275202 | **1.000000** | +1.597633 | 42184008 | 2798 | True |
| ORACULO COMPLETO (razon) | HAP1 | 3 | 20275203 | **1.000000** | +1.597633 | 42281063 | 3016 | True |
| **ORACULO COMPLETO (razon) media** | **HAP1** | | | **1.000000** (sd 0.000000) | **+1.597633** (sd 0.000000) | | | |
| ORACULO COMPLETO (razon) | HT29 | 0 | 20275300 | **1.000000** | +1.768662 | 34886411 | 2351 | True |
| ORACULO COMPLETO (razon) | HT29 | 1 | 20275301 | **1.000000** | +1.768662 | 34800340 | 2131 | True |
| ORACULO COMPLETO (razon) | HT29 | 2 | 20275302 | **1.000000** | +1.768662 | 34815862 | 2300 | True |
| ORACULO COMPLETO (razon) | HT29 | 3 | 20275303 | **1.000000** | +1.768662 | 34786383 | 2315 | True |
| **ORACULO COMPLETO (razon) media** | **HT29** | | | **1.000000** (sd 0.000000) | **+1.768662** (sd 0.000000) | | | |
| ORACULO COMPLETO (razon) | K562 | 0 | 20275400 | **1.000000** | +1.587202 | 34540275 | 3220 | True |
| ORACULO COMPLETO (razon) | K562 | 1 | 20275401 | **1.000000** | +1.587202 | 34315954 | 3237 | True |
| ORACULO COMPLETO (razon) | K562 | 2 | 20275402 | **1.000000** | +1.587202 | 34543356 | 3177 | True |
| ORACULO COMPLETO (razon) | K562 | 3 | 20275403 | **1.000000** | +1.587202 | 34451842 | 3244 | True |
| **ORACULO COMPLETO (razon) media** | **K562** | | | **1.000000** (sd 0.000000) | **+1.587202** (sd 0.000000) | | | |
| ORACULO COMPLETO (razon) | MCF7 | 0 | 20275500 | **1.000000** | +1.147613 | 45644574 | 3035 | True |
| ORACULO COMPLETO (razon) | MCF7 | 1 | 20275501 | **1.000000** | +1.147613 | 45307349 | 2971 | True |
| ORACULO COMPLETO (razon) | MCF7 | 2 | 20275502 | **1.000000** | +1.147613 | 45681388 | 3231 | True |
| ORACULO COMPLETO (razon) | MCF7 | 3 | 20275503 | **1.000000** | +1.147613 | 45827964 | 3019 | True |
| **ORACULO COMPLETO (razon) media** | **MCF7** | | | **1.000000** (sd 0.000000) | **+1.147613** (sd 0.000000) | | | |
| ORACULO DE COLA | A549 | 0 | 20276000 | **0.682765** | +0.673390 | 28489291 | 1939 | True |
| ORACULO DE COLA | A549 | 1 | 20276001 | **0.642045** | +0.533507 | 28256703 | 1939 | True |
| ORACULO DE COLA | A549 | 2 | 20276002 | **0.616477** | +0.445673 | 28184304 | 1939 | True |
| ORACULO DE COLA | A549 | 3 | 20276003 | **0.618371** | +0.452180 | 28446185 | 1939 | True |
| **ORACULO DE COLA media** | **A549** | | | **0.639915** (sd 0.030844) | **+0.526187** (sd 0.105959) | | | |
| ORACULO DE COLA | BXPC3 | 0 | 20276100 | **0.860795** | +0.803006 | 38608513 | 2500 | True |
| ORACULO DE COLA | BXPC3 | 1 | 20276101 | **0.839015** | +0.757516 | 38575286 | 2500 | True |
| ORACULO DE COLA | BXPC3 | 2 | 20276102 | **0.794508** | +0.664557 | 38643979 | 2500 | True |
| ORACULO DE COLA | BXPC3 | 3 | 20276103 | **0.785038** | +0.644778 | 38537705 | 2500 | True |
| **ORACULO DE COLA media** | **BXPC3** | | | **0.819839** (sd 0.036046) | **+0.717464** (sd 0.075286) | | | |
| ORACULO DE COLA | HAP1 | 0 | 20276200 | **0.639205** | +0.572889 | 43323664 | 2666 | True |
| ORACULO DE COLA | HAP1 | 1 | 20276201 | **0.589015** | +0.430339 | 43367305 | 2666 | True |
| ORACULO DE COLA | HAP1 | 2 | 20276202 | **0.696970** | +0.736955 | 43424252 | 2666 | True |
| ORACULO DE COLA | HAP1 | 3 | 20276203 | **0.682765** | +0.696611 | 42999078 | 2666 | True |
| **ORACULO DE COLA media** | **HAP1** | | | **0.651989** (sd 0.048647) | **+0.609198** (sd 0.138169) | | | |
| ORACULO DE COLA | HT29 | 0 | 20276300 | **0.728220** | +0.803631 | 35070843 | 2103 | True |
| ORACULO DE COLA | HT29 | 1 | 20276301 | **0.714015** | +0.753194 | 35348628 | 1892 | True |
| ORACULO DE COLA | HT29 | 2 | 20276302 | **0.699811** | +0.702757 | 35379160 | 2139 | True |
| ORACULO DE COLA | HT29 | 3 | 20276303 | **0.673295** | +0.608608 | 35351753 | 2190 | True |
| **ORACULO DE COLA media** | **HT29** | | | **0.703835** (sd 0.023432) | **+0.717048** (sd 0.083200) | | | |
| ORACULO DE COLA | K562 | 0 | 20276400 | **0.692235** | +0.567754 | 35186022 | 2783 | True |
| ORACULO DE COLA | K562 | 1 | 20276401 | **0.538826** | +0.059598 | 34970925 | 2783 | True |
| ORACULO DE COLA | K562 | 2 | 20276402 | **0.589962** | +0.228984 | 35227033 | 2783 | True |
| ORACULO DE COLA | K562 | 3 | 20276403 | **0.624053** | +0.341907 | 35241962 | 2783 | True |
| **ORACULO DE COLA media** | **K562** | | | **0.611269** (sd 0.064345) | **+0.299561** (sd 0.213138) | | | |
| ORACULO DE COLA | MCF7 | 0 | 20276500 | **0.962121** | +1.060009 | 45906689 | 3352 | True |
| ORACULO DE COLA | MCF7 | 1 | 20276501 | **0.972538** | +1.084100 | 45850385 | 3123 | True |
| ORACULO DE COLA | MCF7 | 2 | 20276502 | **0.949811** | +1.031537 | 46115844 | 3286 | True |
| ORACULO DE COLA | MCF7 | 3 | 20276503 | **0.957386** | +1.049058 | 46062782 | 2747 | True |
| **ORACULO DE COLA media** | **MCF7** | | | **0.960464** (sd 0.009513) | **+1.056176** (sd 0.022001) | | | |
| ORACULO DE BULTO | A549 | 0 | 20277000 | **1.000000** | +1.763175 | 27590644 | 2175 | True |
| ORACULO DE BULTO | A549 | 1 | 20277001 | **1.000000** | +1.763175 | 27398775 | 2166 | True |
| ORACULO DE BULTO | A549 | 2 | 20277002 | **1.000000** | +1.763175 | 27603974 | 2165 | True |
| ORACULO DE BULTO | A549 | 3 | 20277003 | **1.000000** | +1.763175 | 27435052 | 2339 | True |
| **ORACULO DE BULTO media** | **A549** | | | **1.000000** (sd 0.000000) | **+1.763175** (sd 0.000000) | | | |
| ORACULO DE BULTO | BXPC3 | 0 | 20277100 | **1.000000** | +1.093750 | 37518281 | 2873 | True |
| ORACULO DE BULTO | BXPC3 | 1 | 20277101 | **1.000000** | +1.093750 | 37490740 | 2866 | True |
| ORACULO DE BULTO | BXPC3 | 2 | 20277102 | **1.000000** | +1.093750 | 37424441 | 3072 | True |
| ORACULO DE BULTO | BXPC3 | 3 | 20277103 | **1.000000** | +1.093750 | 37475933 | 2734 | True |
| **ORACULO DE BULTO media** | **BXPC3** | | | **1.000000** (sd 0.000000) | **+1.093750** (sd 0.000000) | | | |
| ORACULO DE BULTO | HAP1 | 0 | 20277200 | **1.000000** | +1.597633 | 42496850 | 2707 | True |
| ORACULO DE BULTO | HAP1 | 1 | 20277201 | **1.000000** | +1.597633 | 42374696 | 2880 | True |
| ORACULO DE BULTO | HAP1 | 2 | 20277202 | **1.000000** | +1.597633 | 41986891 | 2883 | True |
| ORACULO DE BULTO | HAP1 | 3 | 20277203 | **1.000000** | +1.597633 | 42399890 | 2823 | True |
| **ORACULO DE BULTO media** | **HAP1** | | | **1.000000** (sd 0.000000) | **+1.597633** (sd 0.000000) | | | |
| ORACULO DE BULTO | HT29 | 0 | 20277300 | **1.000000** | +1.768662 | 34904050 | 2258 | True |
| ORACULO DE BULTO | HT29 | 1 | 20277301 | **1.000000** | +1.768662 | 34768416 | 2173 | True |
| ORACULO DE BULTO | HT29 | 2 | 20277302 | **1.000000** | +1.768662 | 34573628 | 2225 | True |
| ORACULO DE BULTO | HT29 | 3 | 20277303 | **1.000000** | +1.768662 | 35111564 | 2108 | True |
| **ORACULO DE BULTO media** | **HT29** | | | **1.000000** (sd 0.000000) | **+1.768662** (sd 0.000000) | | | |
| ORACULO DE BULTO | K562 | 0 | 20277400 | **1.000000** | +1.587202 | 34454477 | 3224 | True |
| ORACULO DE BULTO | K562 | 1 | 20277401 | **1.000000** | +1.587202 | 34346853 | 3255 | True |
| ORACULO DE BULTO | K562 | 2 | 20277402 | **1.000000** | +1.587202 | 34202635 | 3188 | True |
| ORACULO DE BULTO | K562 | 3 | 20277403 | **1.000000** | +1.587202 | 34314594 | 3121 | True |
| **ORACULO DE BULTO media** | **K562** | | | **1.000000** (sd 0.000000) | **+1.587202** (sd 0.000000) | | | |
| ORACULO DE BULTO | MCF7 | 0 | 20277500 | **1.000000** | +1.147613 | 45425898 | 2869 | True |
| ORACULO DE BULTO | MCF7 | 1 | 20277501 | **1.000000** | +1.147613 | 45563590 | 2906 | True |
| ORACULO DE BULTO | MCF7 | 2 | 20277502 | **1.000000** | +1.147613 | 45508944 | 2941 | True |
| ORACULO DE BULTO | MCF7 | 3 | 20277503 | **1.000000** | +1.147613 | 46020301 | 2979 | True |
| **ORACULO DE BULTO media** | **MCF7** | | | **1.000000** (sd 0.000000) | **+1.147613** (sd 0.000000) | | | |
| SIGNOS ALEATORIOS DENSOS | A549 | 0 | 20273000 | **0.454545** | -0.110605 | 25333498 | 2512 | True |
| SIGNOS ALEATORIOS DENSOS | A549 | 1 | 20273001 | **0.448864** | -0.130124 | 25485255 | 2265 | True |
| SIGNOS ALEATORIOS DENSOS | A549 | 2 | 20273002 | **0.475379** | -0.039037 | 25314310 | 2326 | True |
| SIGNOS ALEATORIOS DENSOS | A549 | 3 | 20273003 | **0.499053** | +0.042290 | 25324441 | 2822 | True |
| SIGNOS ALEATORIOS DENSOS | A549 | 4 | 20273004 | **0.516098** | +0.100846 | 25450877 | 2242 | True |
| SIGNOS ALEATORIOS DENSOS | A549 | 5 | 20273005 | **0.470644** | -0.055303 | 25261325 | 2769 | True |
| SIGNOS ALEATORIOS DENSOS | A549 | 6 | 20273006 | **0.476326** | -0.035784 | 25205850 | 2453 | True |
| SIGNOS ALEATORIOS DENSOS | A549 | 7 | 20273007 | **0.529356** | +0.146389 | 25363910 | 2192 | True |
| **SIGNOS ALEATORIOS DENSOS media** | **A549** | | | **0.483783** (sd 0.028604) | **-0.010166** (sd 0.098263) | | | |
| SIGNOS ALEATORIOS DENSOS | BXPC3 | 0 | 20273100 | **0.484848** | +0.017801 | 34426054 | 3771 | True |
| SIGNOS ALEATORIOS DENSOS | BXPC3 | 1 | 20273101 | **0.540720** | +0.134494 | 34608104 | 3217 | True |
| SIGNOS ALEATORIOS DENSOS | BXPC3 | 2 | 20273102 | **0.568182** | +0.191851 | 34186125 | 2500 | True |
| SIGNOS ALEATORIOS DENSOS | BXPC3 | 3 | 20273103 | **0.501894** | +0.053402 | 34634159 | 3023 | True |
| SIGNOS ALEATORIOS DENSOS | BXPC3 | 4 | 20273104 | **0.528409** | +0.108782 | 34763069 | 2431 | True |
| SIGNOS ALEATORIOS DENSOS | BXPC3 | 5 | 20273105 | **0.508523** | +0.067247 | 34627416 | 4286 | True |
| SIGNOS ALEATORIOS DENSOS | BXPC3 | 6 | 20273106 | **0.484848** | +0.017801 | 34245458 | 3370 | True |
| SIGNOS ALEATORIOS DENSOS | BXPC3 | 7 | 20273107 | **0.481061** | +0.009889 | 34475868 | 3483 | True |
| **SIGNOS ALEATORIOS DENSOS media** | **BXPC3** | | | **0.512311** (sd 0.031158) | **+0.075158** (sd 0.065076) | | | |
| SIGNOS ALEATORIOS DENSOS | HAP1 | 0 | 20273200 | **0.520833** | +0.236686 | 38454674 | 2966 | True |
| SIGNOS ALEATORIOS DENSOS | HAP1 | 1 | 20273201 | **0.503788** | +0.188273 | 38393940 | 2666 | True |
| SIGNOS ALEATORIOS DENSOS | HAP1 | 2 | 20273202 | **0.537879** | +0.285100 | 38741000 | 2666 | True |
| SIGNOS ALEATORIOS DENSOS | HAP1 | 3 | 20273203 | **0.502841** | +0.185584 | 38429453 | 2679 | True |
| SIGNOS ALEATORIOS DENSOS | HAP1 | 4 | 20273204 | **0.511364** | +0.209790 | 38630041 | 2666 | True |
| SIGNOS ALEATORIOS DENSOS | HAP1 | 5 | 20273205 | **0.525568** | +0.250134 | 38384977 | 2666 | True |
| SIGNOS ALEATORIOS DENSOS | HAP1 | 6 | 20273206 | **0.477273** | +0.112964 | 38532374 | 2666 | True |
| SIGNOS ALEATORIOS DENSOS | HAP1 | 7 | 20273207 | **0.500000** | +0.177515 | 38862112 | 2855 | True |
| **SIGNOS ALEATORIOS DENSOS media** | **HAP1** | | | **0.509943** (sd 0.018515) | **+0.205756** (sd 0.052588) | | | |
| SIGNOS ALEATORIOS DENSOS | HT29 | 0 | 20273300 | **0.497159** | -0.016812 | 31362464 | 2129 | True |
| SIGNOS ALEATORIOS DENSOS | HT29 | 1 | 20273301 | **0.497159** | -0.016812 | 31445297 | 2108 | True |
| SIGNOS ALEATORIOS DENSOS | HT29 | 2 | 20273302 | **0.497159** | -0.016812 | 31139662 | 2007 | True |
| SIGNOS ALEATORIOS DENSOS | HT29 | 3 | 20273303 | **0.466856** | -0.124412 | 31220725 | 2334 | True |
| SIGNOS ALEATORIOS DENSOS | HT29 | 4 | 20273304 | **0.466856** | -0.124412 | 31362609 | 1892 | True |
| SIGNOS ALEATORIOS DENSOS | HT29 | 5 | 20273305 | **0.501894** | +0.000000 | 31464177 | 2553 | True |
| SIGNOS ALEATORIOS DENSOS | HT29 | 6 | 20273306 | **0.418561** | -0.295898 | 31354031 | 1892 | True |
| SIGNOS ALEATORIOS DENSOS | HT29 | 7 | 20273307 | **0.445076** | -0.201748 | 31486357 | 1892 | True |
| **SIGNOS ALEATORIOS DENSOS media** | **HT29** | | | **0.473840** (sd 0.030235) | **-0.099613** (sd 0.107358) | | | |
| SIGNOS ALEATORIOS DENSOS | K562 | 0 | 20273400 | **0.576705** | +0.185069 | 30859642 | 2783 | True |
| SIGNOS ALEATORIOS DENSOS | K562 | 1 | 20273401 | **0.484848** | -0.119197 | 30933714 | 2783 | True |
| SIGNOS ALEATORIOS DENSOS | K562 | 2 | 20273402 | **0.468750** | -0.172522 | 30851617 | 2783 | True |
| SIGNOS ALEATORIOS DENSOS | K562 | 3 | 20273403 | **0.505682** | -0.050188 | 30873921 | 2783 | True |
| SIGNOS ALEATORIOS DENSOS | K562 | 4 | 20273404 | **0.637311** | +0.385822 | 30916110 | 2783 | True |
| SIGNOS ALEATORIOS DENSOS | K562 | 5 | 20273405 | **0.501894** | -0.062735 | 30902552 | 2783 | True |
| SIGNOS ALEATORIOS DENSOS | K562 | 6 | 20273406 | **0.533144** | +0.040778 | 30867712 | 2783 | True |
| SIGNOS ALEATORIOS DENSOS | K562 | 7 | 20273407 | **0.590909** | +0.232120 | 30895605 | 2783 | True |
| **SIGNOS ALEATORIOS DENSOS media** | **K562** | | | **0.537405** (sd 0.058736) | **+0.054893** (sd 0.194559) | | | |
| SIGNOS ALEATORIOS DENSOS | MCF7 | 0 | 20273500 | **0.492424** | -0.026281 | 41293511 | 3145 | True |
| SIGNOS ALEATORIOS DENSOS | MCF7 | 1 | 20273501 | **0.473485** | -0.070083 | 41587913 | 3651 | True |
| SIGNOS ALEATORIOS DENSOS | MCF7 | 2 | 20273502 | **0.412879** | -0.210250 | 41364874 | 3139 | True |
| SIGNOS ALEATORIOS DENSOS | MCF7 | 3 | 20273503 | **0.508523** | +0.010951 | 41509266 | 3718 | True |
| SIGNOS ALEATORIOS DENSOS | MCF7 | 4 | 20273504 | **0.534091** | +0.070083 | 41188438 | 3317 | True |
| SIGNOS ALEATORIOS DENSOS | MCF7 | 5 | 20273505 | **0.524621** | +0.048182 | 41397112 | 3774 | True |
| SIGNOS ALEATORIOS DENSOS | MCF7 | 6 | 20273506 | **0.540720** | +0.085414 | 41176876 | 3101 | True |
| SIGNOS ALEATORIOS DENSOS | MCF7 | 7 | 20273507 | **0.486742** | -0.039422 | 41089441 | 3003 | True |
| **SIGNOS ALEATORIOS DENSOS media** | **MCF7** | | | **0.496686** (sd 0.041296) | **-0.016426** (sd 0.095507) | | | |

### TABLA E: verificaciones de las corridas

| | |
|---|---|
| corridas sin fallo | **126** |
| `score_ok` | **126/126** |
| `config_hash` de competencia coincide | **126/126** |
| un solo digest de `config_hash` | `a5519017556705abf9a6e553ebb4e517b9612e9546331ace7c5ec6c5c30a1a13` |
| `poblacion_real_identica` | **126/126** |
| huellas reales distintas por linea | **{'A549': 1, 'BXPC3': 1, 'HAP1': 1, 'HT29': 1, 'K562': 1, 'MCF7': 1}** |
| huellas `pred` distintas | **126/126** |
| `todo_entero` | **126/126** |
| corridas con algun valor negativo | **0** |
| semillas distintas declaradas | **134** |
| choques internos de semilla | **0** |
| choques con familias previas | **0** |
| genes con media de control 0 (inalcanzables) | **6361 a 10252** de 34025 (por linea abajo) |
| fraccion de \|delta real\| en esos genes | **0.000684 a 0.005700** |

**Los genes inalcanzables por `_reescala_entero`, por linea** (media de control 0, asi que ninguna razon los mueve):

| linea | genes con media de control 0 | de `n_var` | fraccion de \|delta real\| que llevan |
|---|---|---|---|
| A549 | **10252** | 34025 | 0.005700 |
| BXPC3 | **6361** | 34025 | 0.001486 |
| HAP1 | **6935** | 34025 | 0.001782 |
| HT29 | **8420** | 34025 | 0.001734 |
| K562 | **9340** | 34025 | 0.003794 |
| MCF7 | **6617** | 34025 | 0.000684 |

---

## 8. Defectos de la construccion, como datos sobre la medicion

1. 🔴 **LOS TRES ORACULOS SON TRANSDUCTIVOS Y NO SON ALCANZABLES POR NINGUNA
   SUBMISSION.** Leen las medias del lado real perturbado de la misma poblacion
   que se puntua. Acotan **la metrica**, no un modelo. Ninguna de sus cifras
   puede entrar en una submission ni en un rasgo de un candidato. Esta
   advertencia es el punto 1 de la lista a proposito.
2. **La realizacion por razon no alcanza los genes que el control no expresa.**
   `_reescala_entero` solo toca entradas **no nulas**, asi que un gen con media
   de control 0 queda en 0 por mucho que valga `r`. Medido: ****6361 a 10252****
   de 34025 genes por linea, que llevan el ****0.07 % a 0.57 %**** del
   `|delta real|` total. Es la razon de ser del brazo 1a (copia), que **no**
   tiene esa limitacion, y la distancia entre 1a y 1b la acota medida.
3. 🔴 **El ORACULO DE COLA no emite nada donde el ED real esta vacio.** En esas
   perturbaciones es el SUELO por construccion, y son **15, 0, 14, 6, 8 y 0** de
   33. La media del panel de ese brazo esta **diluida hacia el azar por
   construccion**, no por falta de senal: por eso la fraccion se da dos veces,
   sobre las 33 (la TABLA B3) y solo donde la cola existe (la TABLA B3b). **Las dos son la
   misma medicion contada sobre dos cohortes distintas**, y ninguna de las dos
   es "la" fraccion.
4. **La cola es un conjunto de genes dependiente de la poblacion real**, no una
   particion fija. Sale del DE del propio marcador a `p_adj < 0.05` sobre
   **esta** poblacion real (`seed_real = 20260001 + fold`); otra poblacion daria
   otra cola. **No se midio** con otro umbral, ni con una cola definida por
   magnitud en vez de por significancia.
5. 🔴 **`pds_cosine` no es aditivo sobre genes.** Es un estadistico de **rango**:
   en cuanto el delta propio queda mas cerca que los otros 32, mas senal no
   suma nada. Por eso cola y bulto **particionan los genes pero no el valor de
   la metrica**, y su suma no tiene por que dar el completo. la TABLA B3 reporta la
   suma medida para que se vea cuanto se aparta, no para corregirla.
6. **Los 33 genes blanco del panel estan FUERA del espacio de rasgos** de
   `pds_cosine` (`exclusion_scope="panel"`), asi que cualquier emision sobre
   ellos es **invisible** para esta metrica. Medido: de los pares (blanco, gen)
   del ED real, solo el ****0.44 % a 4.28 %**** cae en genes del panel, asi
   que la exclusion apenas toca la cola -- pero toca, y queda declarado.
7. **`n` desigual entre brazos.** Los tres oraculos llevan **4** semillas y el
   brazo de signos **8**, como el suelo y los dos brazos de gwps. Las
   resoluciones (`SE = sd/sqrt(n)`) **no son comparables entre brazos de `n`
   distinto** y se dan con su `n` en la tabla. La eleccion de 4 y 8 se fijo
   **antes** de la primera corrida y no se cambio al ver ningun valor.
8. **El brazo 1a (copia) tiene `n = 1`** porque es determinista dado
   `seed_real`: las 400 celulas predichas de `h` son las 400 reales de `h`. No
   tiene `sd` y su celda de resolucion queda **n/d**, no 0.
9. **El multiconjunto de `|s|` del brazo 4 es una ELECCION DE ESCALA heredada**,
   no una magnitud medida: es el rango normalizado de `tres_vias_rango.md` §2.2.
   "La misma distribucion de magnitud que gwps" significa *la misma que esa
   eleccion*, no la misma que la biologia. Se verifico que el brazo preserva la
   densidad (**264/264** pares blanco-sorteo con los mismos genes no nulos) y el
   multiconjunto de `|s|` (**264/264** iguales a `1.1e-16`, medio `eps` de
   float64, que es el redondeo de `(1+s)-1` y no de la construccion), y que el
   recorte a `[0, 2]` **nunca actua** porque `|s| <= 1`.
10. **La descomposicion del delta (la TABLA B0) depende del espacio.** Se midio en
    `bulk_lognorm` con `bulk_target_sum = 50000.0`, leido de `run_params.yaml`
    de las corridas ya hechas. Si ese valor cambia, el reparto cambia.
11. **El `b` y el `r` oficiales son UNA realizacion y CINCO particiones**
    respectivamente, no 8 semillas. `sd_b` sale de las 8 realizaciones de
    `dispersion_b`, que es **otra** cantidad que la `sd` de la unica `b` del
    bundle; y `r_sd` es la `sd` entre 5 particiones, no el error de `r`. Las
    tres se dan separadas en la TABLA A1 y no se mezclan.
12. **No se midio**: los oraculos con 8 semillas, ni el brazo de signos con otra
    familia de sorteo, ni ningun brazo con el gen blanco DENTRO del espacio de
    `pds_cosine`, ni una cola definida sobre el ED **predicho** en vez del real,
    ni el efecto de estos brazos en las otras cinco metricas mas alla de lo que
    la tabla de crudas registra.
