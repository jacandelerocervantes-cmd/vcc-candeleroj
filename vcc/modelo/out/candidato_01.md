# Candidato 01: respuesta lineal a lo largo de la covarianza del control

Su trabajo era uno solo: mostrar habilidad **especifica del blanco**, por encima
de un `b` que emite el control escalado a UN perfil de respuesta promedio, el
mismo para las 33 perturbaciones.

**No la muestra.** La refutacion esta medida en las seis lineas y es el resultado
de este encargo. La medicion completa se corrio igual, por tres razones que se
declaran: el punto 5 (como mueve las cohortes una emision con mucha senal) es
independiente de que el candidato tenga habilidad y cierra el hueco que §3.17
dejo abierto; '¿supera al `b`?' merece quedar registrado con su margen; y la
dispersion de un candidato real no se habia medido nunca. **No se itero ni se
ajusto nada mirando la `s`.**

## 1. 🔴 LA REFUTACION, primero porque es el resultado

**La hipotesis.** Para el blanco `g`, los genes que mas covarian con `g` en el
pool de control de ese contexto son los candidatos a moverse cuando `g` se
perturba.

**La medicion.** Spearman del delta de pseudobulk contra el control, sobre los
genes expresados (>= 1 en CP50k en los DOS controles). Cuatro cantidades, y las
dos ultimas son las que impiden una conclusion falsa:

| clave | que es | para que |
|---|---|---|
| **A** | mitad-contra-mitad del REAL | el **techo alcanzable** a 400 celulas. Si fuera ~0, el test no tendria potencia y no se podria concluir NADA |
| **B** | real de un blanco vs real de OTRO | la respuesta **generica** compartida. Es la cantidad que el `b` captura |
| **C** | la prediccion vs el real del MISMO blanco | lo que se quiere medir |
| **D** | la prediccion vs el real de OTRO blanco | **control negativo**. `C - D` es lo especifico del blanco que el candidato captura |

| linea | genes expresados | celulas de control | A techo | B generica | C candidato | D control negativo | C - D | concordancia de signo |
|---|---|---|---|---|---|---|---|---|
| A549 | 6934 | 1379 | +0.1575 | +0.2004 | **+0.0009** | -0.0045 | **+0.0053** | **0.4996** |
| BXPC3 | 7475 | 3228 | +0.1544 | +0.1198 | **+0.0009** | -0.0067 | **+0.0076** | **0.4964** |
| HAP1 | 7857 | 3172 | +0.0820 | +0.1053 | **+0.0046** | +0.0001 | **+0.0044** | **0.5159** |
| HT29 | 6999 | 2736 | +0.1020 | +0.1298 | **-0.0054** | -0.0052 | **-0.0002** | **0.4959** |
| K562 | 7403 | 1398 | +0.1609 | +0.2184 | **+0.0108** | +0.0113 | **-0.0005** | **0.5094** |
| MCF7 | 7266 | 2669 | +0.2012 | +0.2133 | **-0.0042** | -0.0099 | **+0.0057** | **0.5021** |

**`A` corre de +0.0820 a +0.2012 en las seis lineas: el test tiene
potencia.** A 400 celulas por perturbacion el delta real es reproducible, asi que
un modelo con senal tendria donde mostrarla.

**`C` corre de -0.0054 a +0.0108: cero.** Y la concordancia de signo
-- la fraccion de genes que el lado real mueve de verdad (`|dlog| > 0.25`) en que
el signo predicho acierta -- va de 0.4959 a 0.5159. **0.5 es azar
exacto.** El candidato no acierta el signo de la respuesta mejor que una moneda.

**Y `D` descarta la explicacion alternativa.** `C - D` va de -0.0005 a +0.0076: el candidato no captura ni lo
especifico del blanco ni la respuesta generica. Si hubiera capturado lo generico,
`D` habria salido alto y `C - D` bajo; salen los dos en cero.

🔴 **Y el hallazgo que sale de paso y vale mas que el candidato: `B > A` en cinco
de las seis lineas.** La respuesta compartida entre perturbaciones distintas es
mas fuerte que la reproducibilidad de una sola perturbacion consigo misma:

| linea | A (una perturbacion consigo misma) | B (entre perturbaciones distintas) | B - A | B > A |
|---|---|---|---|---|
| A549 | +0.1575 | +0.2004 | +0.0429 | **si** |
| BXPC3 | +0.1544 | +0.1198 | -0.0347 | no |
| HAP1 | +0.0820 | +0.1053 | +0.0233 | **si** |
| HT29 | +0.1020 | +0.1298 | +0.0278 | **si** |
| K562 | +0.1609 | +0.2184 | +0.0575 | **si** |
| MCF7 | +0.2012 | +0.2133 | +0.0120 | **si** |

Medido en contexto, sin tocar el `b` ni el anchor. Explica mecanicamente por que
el `b` es un comparador tan duro: **esta capturando la componente dominante de la
varianza**, y la parte especifica del blanco que queda por encima es la minoritaria
incluso en el lado real. Un candidato que quiera superar al `b` tiene que acertar
la componente pequena, no la grande.

### 1.1 Las dos magnitudes que se midieron y se descartaron ANTES de puntuar

Ninguna de las dos se puntuo nunca: se descartaron por una propiedad medida de la
EMISION, no por su `s`. Queda escrito para que no se repitan.

| version | formula | lo que se midio en A549 | por que se descarto |
|---|---|---|---|
| v1 | `beta_h = cov/var(y_g)`, `Delta_y_h = -beta_h * mean(y_g)` (apagar `g` del todo) | `r_min = 0.0000` con `r_max = 1.0000` en **24 de 33** blancos; cientos de genes por blanco a cero exacto (RAP1GAP 381, ZBP1 335, STAT4 332) | `beta_h = rho * sd_h/sd_g` crece sin cota cuando `sd_h >> sd_g`; la extrapolacion sale del soporte |
| v2 | `Delta_y_h = -rho_h * sd(y_h)` (esperanza condicional a una `sd`) | la **mediana** de los genes que responden va a cero; **18 de 33** blancos saturados; `\|log2FC\|` de pseudobulk hasta **17.4** | para un gen disperso `mean/sd ~ sqrt(p/(1-p))` = 0.23 a `p = 0.05`, asi que cualquier `\|rho\| > 0.23` lo cruza por debajo de cero. log1p de cuentas dispersas esta inflado de ceros, no es una conjunta normal |
| **v3, la que corre** | `r_h = 1 - rho_h`, acotada en `[0, 2]` | **0 de 33** saturados | adimensional y acotada por construccion; la unica cantidad que entra es `rho` en contexto, que es literalmente 'el orden y el signo, nunca el valor' |

La magnitud queda declarada como **hueco**, no como resultado: que la fraccion de
cambio sea `rho` es una eleccion de escala, no una cifra medida de cuanto baja un
gen cuando se perturba su regulador. Eso exigiria el lado real perturbado, que
esta prohibido.

## 2. Lo que vario y lo que no, con la prueba

**Vario solo la semilla de emision.** Familia nueva `20266000 + 100*fold + j`,
`j = 0..7`. Ninguna colision con `2026000x` (verdad), `2026100x` (suelo),
`2026200x`/`2026300x` (techo), `2026400x` (dispersion de `u`), `2026500x`
(dispersion de `b`) ni la `0` del baseline/anchor.

**No vario la poblacion real**, ni **la estructura del candidato**. Los ratios por
blanco son deterministicos dado el pool de control (`construye_respuesta` no usa
`rng`), asi que `sd_s` mide **solo el muestreo de la emision**, que es lo que se
quiere medir.

**No transductivo, verificado en el codigo del arnes y no asumido.**
`build_fold_adatas` llama `predict_fn(control_pool, panel + [CONTROL_LABEL],
rng_pred)` (`suelo_arnes.py:176`), y `control_pool` sale de
`gene_local == CONTROL_LABEL` (147-148). La funcion no recibe la identidad de la
linea ni una sola celula perturbada. `var_names` entra por clausura porque la firma
del arnes no lo lleva y sin el no se puede localizar la columna del gen blanco: son
los nombres de las columnas del dataset, no el lado real perturbado. Queda
declarado. **El arnes no se modifico** -- `predict_fn` es un parametro suyo.

**Ninguna cifra del `b` ni del anchor entra como entrada**: son artefactos de
puntuacion. Nada importado: ni Replogle, ni el benchmark, ni tabla de pesos.

| linea | huellas del REAL | huellas del PREDICHO | huella del real = dispersion_u.md | config_hash de competencia en las 8 | todo entero | hay negativos | nnz | valor max |
|---|---|---|---|---|---|---|---|---|
| A549 | **1** | **8** | **si** | **si** | **si** | **no** | 28,141,912-28,554,810 | 1939 |
| BXPC3 | **1** | **8** | **si** | **si** | **si** | **no** | 38,386,943-38,669,423 | 2500 |
| HAP1 | **1** | **8** | **si** | **si** | **si** | **no** | 43,169,011-43,498,026 | 2874 |
| HT29 | **1** | **8** | **si** | **si** | **si** | **no** | 35,106,403-35,568,304 | 2081 |
| K562 | **1** | **8** | **si** | **si** | **si** | **no** | 34,943,412-35,350,289 | 2783 |
| MCF7 | **1** | **8** | **si** | **si** | **si** | **no** | 45,702,304-46,283,570 | 2310 |

Una sola huella del real por linea y ocho del predicho en las seis. El
`config_hash` recalculado coincide con `a5519017556705ab...` en las **48**.

**El contrato de emision se comprueba en cada corrida, no se asume.** Conteos
enteros muestreados: se parte de celulas de control remuestreadas con reemplazo y
se reescala gen a gen con `r_h <= 1: Binomial(x, r)` y `r_h > 1: x + Poisson(x(r-1))`.
Las dos ramas tienen esperanza exactamente `x*r` y las dos devuelven enteros. No
hay ningun perfil esperado suave, que es lo que `paso1_contrato.md` prohibe por el
tope `HARD_CAP_NNZ`.

Ninguna corrida fallo: 48 de 48 `run` y 48 de 48 `score --real-bundle`.

## 3. El agregado y las 36 celdas

| linea | los 8 avg_score | media | sd (ddof=1) | rango | el suelo (dispersion_u.md) | candidato - suelo |
|---|---|---|---|---|---|---|
| A549 | -0.487 / -0.577 / -0.272 / -0.366 / +0.011 / -0.504 / -0.594 / -1.203 | **-0.4991** | 0.346559 | 1.214119 | -0.1486 | -0.3505 |
| BXPC3 | -0.180 / -0.206 / -0.132 / -0.167 / -0.183 / -0.174 / -0.160 / -0.172 | **-0.1718** | 0.021159 | 0.074054 | -0.1945 | +0.0227 |
| HAP1 | -0.126 / -0.202 / -0.303 / -0.133 / -0.309 / -0.214 / -0.340 / -0.275 | **-0.2377** | 0.081401 | 0.213903 | -0.1163 | -0.1215 |
| HT29 | -0.136 / -0.113 / -0.060 / -0.108 / -0.160 / -0.061 / -0.127 / -0.183 | **-0.1185** | 0.043380 | 0.123189 | -0.1112 | -0.0073 |
| K562 | -0.121 / -0.116 / -0.166 / -0.105 / -0.159 / -0.200 / -0.148 / -0.066 | **-0.1353** | 0.041697 | 0.134117 | -0.1702 | +0.0350 |
| MCF7 | -0.120 / -0.147 / -0.205 / -0.106 / -0.146 / -0.136 / -0.163 / -0.125 | **-0.1435** | 0.030643 | 0.099328 | -0.1591 | +0.0156 |

`s` es la columna `from_replicate` de `scored.csv`, que con un bundle de
competencia es la que toma `avg_score`.

| linea | metrica | sentido | los 8 valores de s | media | sd_s (ddof=1) | sd_s (ddof=0) | rango | min | max | el suelo |
|---|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | alto | +0.322 / +0.133 / +0.098 / +0.016 / +0.224 / +0.394 / -0.234 / -0.153 | **+0.1000** | 0.218700 | 0.204575 | 0.627846 | -0.2342 | +0.3936 | +0.0573 |
| A549 | `mse_capped_norm` | bajo | +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 | **+0.0000** | 0.000000 | 0.000000 | 0.000000 | +0.0000 | +0.0000 | +0.0000 |
| A549 | `lfc_nmae` | bajo | -0.311 / -0.485 / -0.288 / -0.367 / -0.372 / -0.337 / -0.457 / -0.301 | **-0.3648** | 0.072338 | 0.067666 | 0.196414 | -0.4847 | -0.2883 | -0.3954 |
| A549 | `fid` | alto | +0.285 / +0.239 / +0.129 / +0.192 / +0.199 / +0.299 / +0.319 / +0.280 | **+0.2429** | 0.064967 | 0.060771 | 0.189981 | +0.1291 | +0.3191 | +0.0000 |
| A549 | `reach` | alto | +0.023 / -0.113 / +0.051 / -0.419 / +0.012 / -0.138 / +0.046 / -0.568 | **-0.1382** | 0.233685 | 0.218592 | 0.619197 | -0.5678 | +0.0514 | -0.3510 |
| A549 | `jac` | alto | -3.239 / -3.239 / -1.620 / -1.620 / +0.000 / -3.239 / -3.239 / -6.479 | **-2.8346** | 1.886945 | 1.765076 | 6.478978 | -6.4790 | +0.0000 | -0.2025 |
| BXPC3 | `pds_cosine` | alto | -0.008 / -0.045 / +0.243 / +0.099 / -0.047 / +0.018 / +0.220 / +0.020 | **+0.0623** | 0.114144 | 0.106772 | 0.290744 | -0.0475 | +0.2433 | +0.1063 |
| BXPC3 | `mse_capped_norm` | bajo | +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 | **+0.0000** | 0.000000 | 0.000000 | 0.000000 | +0.0000 | +0.0000 | +0.0000 |
| BXPC3 | `lfc_nmae` | bajo | -0.712 / -0.728 / -0.667 / -0.708 / -0.719 / -0.706 / -0.727 / -0.699 | **-0.7083** | 0.019515 | 0.018255 | 0.061013 | -0.7279 | -0.6669 | -0.7241 |
| BXPC3 | `fid` | alto | +0.108 / +0.099 / +0.114 / +0.107 / +0.108 / +0.078 / +0.125 / +0.116 | **+0.1068** | 0.013687 | 0.012803 | 0.046301 | +0.0784 | +0.1247 | -0.0011 |
| BXPC3 | `reach` | alto | -0.471 / -0.562 / -0.486 / -0.502 / -0.443 / -0.432 / -0.585 / -0.468 | **-0.4935** | 0.054354 | 0.050843 | 0.153251 | -0.5849 | -0.4317 | -0.5432 |
| BXPC3 | `jac` | alto | +0.002 / -0.000 / +0.005 / +0.000 / +0.002 / +0.000 / +0.009 / -0.001 | **+0.0021** | 0.003491 | 0.003266 | 0.010572 | -0.0011 | +0.0095 | -0.0049 |
| HAP1 | `pds_cosine` | alto | +0.597 / +0.393 / +0.436 / +0.148 / +0.234 / +0.118 / +0.124 / +0.371 | **+0.3026** | 0.173963 | 0.162727 | 0.478752 | +0.1183 | +0.5971 | +0.2838 |
| HAP1 | `mse_capped_norm` | bajo | +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 | **+0.0000** | 0.000000 | 0.000000 | 0.000000 | +0.0000 | +0.0000 | +0.0000 |
| HAP1 | `lfc_nmae` | bajo | -0.221 / -0.201 / -0.305 / -0.214 / -0.247 / -0.180 / -0.273 / -0.115 | **-0.2193** | 0.058245 | 0.054483 | 0.189666 | -0.3045 | -0.1149 | -0.2079 |
| HAP1 | `fid` | alto | +0.077 / +0.001 / +0.000 / +0.001 / +0.051 / +0.046 / +0.152 / +0.131 | **+0.0574** | 0.059275 | 0.055447 | 0.152389 | +0.0000 | +0.1524 | +0.0154 |
| HAP1 | `reach` | alto | -0.523 / -0.719 / -0.571 / -0.052 / -0.509 / -0.575 / -0.663 / -0.656 | **-0.5336** | 0.207686 | 0.194272 | 0.666873 | -0.7189 | -0.0520 | -0.5301 |
| HAP1 | `jac` | alto | -0.687 / -0.687 / -1.380 / -0.683 / -1.380 / -0.690 / -1.380 / -1.380 | **-1.0334** | 0.370632 | 0.346695 | 0.696691 | -1.3801 | -0.6834 | -0.2588 |
| HT29 | `pds_cosine` | alto | -0.111 / -0.024 / +0.168 / +0.024 / -0.121 / +0.077 / -0.356 / -0.155 | **-0.0622** | 0.161389 | 0.150965 | 0.524546 | -0.3564 | +0.1681 | +0.0088 |
| HT29 | `mse_capped_norm` | bajo | +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 | **+0.0000** | 0.000000 | 0.000000 | 0.000000 | +0.0000 | +0.0000 | +0.0000 |
| HT29 | `lfc_nmae` | bajo | -0.454 / -0.400 / -0.428 / -0.464 / -0.445 / -0.434 / -0.424 / -0.362 | **-0.4264** | 0.032348 | 0.030258 | 0.101089 | -0.4636 | -0.3625 | -0.3973 |
| HT29 | `fid` | alto | -0.042 / +0.089 / +0.081 / -0.046 / -0.045 / -0.103 / -0.053 / -0.040 | **-0.0198** | 0.068094 | 0.063696 | 0.192651 | -0.1032 | +0.0895 | -0.3463 |
| HT29 | `reach` | alto | -0.547 / -0.519 / -0.356 / -0.500 / -0.521 / -0.413 / -0.272 / -0.880 | **-0.5013** | 0.180395 | 0.168744 | 0.607807 | -0.8803 | -0.2725 | -0.4126 |
| HT29 | `jac` | alto | +0.337 / +0.174 / +0.175 / +0.341 / +0.175 / +0.508 / +0.342 / +0.339 | **+0.2988** | 0.117586 | 0.109991 | 0.333687 | +0.1743 | +0.5080 | +0.4803 |
| K562 | `pds_cosine` | alto | -0.003 / +0.144 / -0.191 / -0.003 / +0.006 / -0.411 / +0.207 / +0.038 | **-0.0267** | 0.194544 | 0.181980 | 0.617942 | -0.4109 | +0.2070 | -0.0345 |
| K562 | `mse_capped_norm` | bajo | +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 | **+0.0000** | 0.000000 | 0.000000 | 0.000000 | +0.0000 | +0.0000 | +0.0000 |
| K562 | `lfc_nmae` | bajo | -0.349 / -0.449 / -0.493 / -0.376 / -0.451 / -0.478 / -0.417 / -0.383 | **-0.4245** | 0.051443 | 0.048121 | 0.143318 | -0.4928 | -0.3495 | -0.4803 |
| K562 | `fid` | alto | +0.009 / +0.019 / -0.098 / -0.061 / -0.042 / -0.013 / -0.064 / -0.055 | **-0.0381** | 0.039836 | 0.037263 | 0.116095 | -0.0976 | +0.0185 | -0.1659 |
| K562 | `reach` | alto | -0.387 / -0.420 / -0.330 / -0.198 / -0.579 / -0.308 / -0.620 / -0.111 | **-0.3692** | 0.173582 | 0.162371 | 0.508928 | -0.6200 | -0.1111 | -0.4361 |
| K562 | `jac` | alto | +0.002 / +0.008 / +0.114 / +0.008 / +0.112 / +0.010 / +0.003 / +0.117 | **+0.0469** | 0.056128 | 0.052503 | 0.115539 | +0.0017 | +0.1172 | +0.0954 |
| MCF7 | `pds_cosine` | alto | +0.191 / +0.107 / -0.162 / +0.226 / +0.103 / +0.114 / -0.020 / +0.197 | **+0.0944** | 0.128961 | 0.120632 | 0.387648 | -0.1621 | +0.2256 | +0.0523 |
| MCF7 | `mse_capped_norm` | bajo | +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 / +0.000 | **+0.0000** | 0.000000 | 0.000000 | 0.000000 | +0.0000 | +0.0000 | +0.0000 |
| MCF7 | `lfc_nmae` | bajo | -0.803 / -0.847 / -0.930 / -0.769 / -0.854 / -0.805 / -0.869 / -0.843 | **-0.8398** | 0.049054 | 0.045886 | 0.161002 | -0.9298 | -0.7688 | -0.8778 |
| MCF7 | `fid` | alto | +0.019 / +0.008 / +0.018 / +0.012 / +0.009 / +0.017 / +0.029 / +0.035 | **+0.0183** | 0.009673 | 0.009048 | 0.027527 | +0.0079 | +0.0354 | -0.0003 |
| MCF7 | `reach` | alto | -0.133 / -0.152 / -0.157 / -0.107 / -0.138 / -0.142 / -0.128 / -0.146 | **-0.1378** | 0.015690 | 0.014677 | 0.050258 | -0.1574 | -0.1071 | -0.1281 |
| MCF7 | `jac` | alto | +0.006 / +0.003 / +0.000 / +0.003 / +0.006 / -0.000 / +0.009 / +0.005 | **+0.0039** | 0.003054 | 0.002856 | 0.008981 | -0.0004 | +0.0086 | -0.0004 |

## 4. 🔴 ¿Supera al `b`? Celda por celda, con el margen ABSOLUTO de los tres terminos

Superar al `b` es, por definicion de la escala, `s > 0`: `s = (u - b)/(r - b)` vale
exactamente 0 cuando `u = b`. Esta es la pregunta **absoluta**, no una comparacion
entre dos candidatos, asi que lleva el margen completo de §3.19 con sus tres
terminos -- incluido el de `b`, que no se cancela aqui:

```
margen = sqrt( sd_s^2 + (|s| * SE_r/|r-b|)^2 + (|s-1| * sd_b/|r-b|)^2 )
```

`sd_s` es la **del candidato**, medida con sus 8 semillas, no la del suelo (punto 1
del encargo). `SE_r = sd_r/sqrt(5)` de `anclas_nativas_v2.md`; `sd_b` de las 8
semillas de `dispersion_b.md`. La columna de veredicto exige `s > margen`.

| linea | metrica | s medio | termino de u (sd_s) | termino de r | termino de b | margen absoluto | s / margen | veredicto |
|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | +0.1000 | 0.218700 | 0.004462 | 0.095849 | **0.238823** | +0.42 | indistinguible de `b` |
| A549 | `mse_capped_norm` | +0.0000 | 0.000000 | 0.000000 | 0.074293 | **0.074293** | +0.00 | indistinguible de `b` |
| A549 | `lfc_nmae` | -0.3648 | 0.072338 | 0.009308 | 0.066953 | **0.099006** | -3.68 | **por debajo** |
| A549 | `fid` | +0.2429 | 0.064967 | 0.009612 | 0.067431 | **0.094128** | +2.58 | 🔴 **SUPERA** |
| A549 | `reach` | -0.1382 | 0.233685 | 0.012379 | 0.227034 | **0.326047** | -0.42 | indistinguible de `b` |
| A549 | `jac` | -2.8346 | 1.886945 | 2.372039 | 2.195919 | **3.742886** | -0.76 | indistinguible de `b` |
| BXPC3 | `pds_cosine` | +0.0623 | 0.114144 | 0.000979 | 0.044250 | **0.122425** | +0.51 | indistinguible de `b` |
| BXPC3 | `mse_capped_norm` | +0.0000 | 0.000000 | 0.000000 | 0.051808 | **0.051808** | +0.00 | indistinguible de `b` |
| BXPC3 | `lfc_nmae` | -0.7083 | 0.019515 | 0.060225 | 0.031375 | **0.070656** | -10.03 | **por debajo** |
| BXPC3 | `fid` | +0.1068 | 0.013687 | 0.003658 | 0.003223 | **0.014529** | +7.35 | 🔴 **SUPERA** |
| BXPC3 | `reach` | -0.4935 | 0.054354 | 0.020724 | 0.188457 | **0.197230** | -2.50 | **por debajo** |
| BXPC3 | `jac` | +0.0021 | 0.003491 | 0.000106 | 0.000940 | **0.003617** | +0.57 | indistinguible de `b` |
| HAP1 | `pds_cosine` | +0.3026 | 0.173963 | 0.005861 | 0.068766 | **0.187153** | +1.62 | 🔴 **SUPERA** |
| HAP1 | `mse_capped_norm` | +0.0000 | 0.000000 | 0.000000 | 0.079353 | **0.079353** | +0.00 | indistinguible de `b` |
| HAP1 | `lfc_nmae` | -0.2193 | 0.058245 | 0.007888 | 0.034462 | **0.068134** | -3.22 | **por debajo** |
| HAP1 | `fid` | +0.0574 | 0.059275 | 0.004427 | 0.038128 | **0.070618** | +0.81 | indistinguible de `b` |
| HAP1 | `reach` | -0.5336 | 0.207686 | 0.048901 | 0.282240 | **0.353813** | -1.51 | **por debajo** |
| HAP1 | `jac` | -1.0334 | 0.370632 | 0.888896 | 0.649543 | **1.161641** | -0.89 | indistinguible de `b` |
| HT29 | `pds_cosine` | -0.0622 | 0.161389 | 0.001687 | 0.130991 | **0.207865** | -0.30 | indistinguible de `b` |
| HT29 | `mse_capped_norm` | +0.0000 | 0.000000 | 0.000000 | 0.035282 | **0.035282** | +0.00 | indistinguible de `b` |
| HT29 | `lfc_nmae` | -0.4264 | 0.032348 | 0.006562 | 0.090225 | **0.096073** | -4.44 | **por debajo** |
| HT29 | `fid` | -0.0198 | 0.068094 | 0.001824 | 0.128530 | **0.145465** | -0.14 | indistinguible de `b` |
| HT29 | `reach` | -0.5013 | 0.180395 | 0.031561 | 0.155487 | **0.240238** | -2.09 | **por debajo** |
| HT29 | `jac` | +0.2988 | 0.117586 | 0.043166 | 0.132299 | **0.182189** | +1.64 | 🔴 **SUPERA** |
| K562 | `pds_cosine` | -0.0267 | 0.194544 | 0.000693 | 0.112678 | **0.224820** | -0.12 | indistinguible de `b` |
| K562 | `mse_capped_norm` | +0.0000 | 0.000000 | 0.000000 | 0.040365 | **0.040365** | +0.00 | indistinguible de `b` |
| K562 | `lfc_nmae` | -0.4245 | 0.051443 | 0.004649 | 0.083265 | **0.097985** | -4.33 | **por debajo** |
| K562 | `fid` | -0.0381 | 0.039836 | 0.001755 | 0.000528 | **0.039878** | -0.96 | indistinguible de `b` |
| K562 | `reach` | -0.3692 | 0.173582 | 0.022561 | 0.223848 | **0.284161** | -1.30 | **por debajo** |
| K562 | `jac` | +0.0469 | 0.056128 | 0.005120 | 0.000574 | **0.056364** | +0.83 | indistinguible de `b` |
| MCF7 | `pds_cosine` | +0.0944 | 0.128961 | 0.001596 | 0.068615 | **0.146087** | +0.65 | indistinguible de `b` |
| MCF7 | `mse_capped_norm` | +0.0000 | 0.000000 | 0.000000 | 0.039528 | **0.039528** | +0.00 | indistinguible de `b` |
| MCF7 | `lfc_nmae` | -0.8398 | 0.049054 | 0.020894 | 0.070670 | **0.088527** | -9.49 | **por debajo** |
| MCF7 | `fid` | +0.0183 | 0.009673 | 0.001596 | 0.001822 | **0.009971** | +1.83 | 🔴 **SUPERA** |
| MCF7 | `reach` | -0.1378 | 0.015690 | 0.006820 | 0.040351 | **0.043828** | -3.14 | **por debajo** |
| MCF7 | `jac` | +0.0039 | 0.003054 | 0.000588 | 0.002091 | **0.003748** | +1.05 | 🔴 **SUPERA** |

**Supera al `b` por mas que su margen en 6 de 36 celdas**; queda por
debajo en 11; indistinguible en 19.

Las que superan:

| linea | metrica | s medio | margen | s / margen | u medio | b | r | el suelo en la misma celda |
|---|---|---|---|---|---|---|---|---|
| BXPC3 | `fid` | +0.1068 | 0.014529 | **+7.35** | 0.0653 | 0.0014 | 0.5993 | -0.0011 |
| A549 | `fid` | +0.2429 | 0.094128 | **+2.58** | 0.1015 | 0.0000 | 0.4179 | +0.0000 |
| MCF7 | `fid` | +0.0183 | 0.009971 | **+1.83** | 0.0086 | 0.0002 | 0.4626 | -0.0003 |
| HT29 | `jac` | +0.2988 | 0.182189 | **+1.64** | 0.1452 | 0.0911 | 0.2721 | +0.4803 |
| HAP1 | `pds_cosine` | +0.3026 | 0.187153 | **+1.62** | 0.5440 | 0.4375 | 0.7896 | +0.2838 |
| MCF7 | `jac` | +0.0039 | 0.003748 | **+1.05** | 0.0003 | 0.0000 | 0.0722 | -0.0004 |

🔴 **La comprobacion que hace falta para que esa tabla no engañe: ¿superan
tambien al SUELO?** El suelo es una copia remuestreada del control, sin ninguna
estructura. Si el candidato no le gana a eso en una celda, lo que esa celda
premia no es su estructura. Se compara contra el mismo `sd_s` del candidato,
porque las dos emisiones se puntuaron contra el mismo bundle y ahi `r` y `b` se
cancelan (§7):

| linea | metrica | s del candidato | s del suelo | diferencia | resolucion k=8 | dif / resolucion | b | el b es ~0 | veredicto contra el suelo |
|---|---|---|---|---|---|---|---|---|---|
| BXPC3 | `fid` | +0.1068 | -0.0011 | **+0.1079** | 0.006843 | +15.77 | 0.0014 | 🔴 **b ~ 0** | gana al suelo |
| A549 | `fid` | +0.2429 | +0.0000 | **+0.2429** | 0.032484 | +7.48 | 0.0000 | 🔴 **b ~ 0** | gana al suelo |
| MCF7 | `fid` | +0.0183 | -0.0003 | **+0.0186** | 0.004836 | +3.84 | 0.0002 | 🔴 **b ~ 0** | gana al suelo |
| HT29 | `jac` | +0.2988 | +0.4803 | **-0.1815** | 0.058793 | -3.09 | 0.0911 |  | **PIERDE con el suelo** |
| HAP1 | `pds_cosine` | +0.3026 | +0.2838 | **+0.0188** | 0.086981 | +0.22 | 0.4375 |  | indistinguible del suelo |
| MCF7 | `jac` | +0.0039 | -0.0004 | **+0.0043** | 0.001527 | +2.83 | 0.0000 | 🔴 **b ~ 0** | gana al suelo |

**De las 6 celdas que superan al `b`, 3 son `fid` y 4 tienen
`b` practicamente cero.** Y una, HT29/`jac`, **queda por debajo del suelo**: el
suelo cobra +0.4803 ahi y el candidato +0.2988. Esa celda no premia estructura,
premia la convencion `J(vacio, vacio) = 1` -- y `dispersion_b.md` §3.2 midio que
la mitad de ese +0.4803 la pone el `b` del bundle por haber caido en el minimo
de sus nueve realizaciones. El candidato, al predecir genes, PIERDE parte de ese
regalo.

**El cierre del punto 3, contado.** Celdas que (a) superan al `b` por mas que su
margen, (b) superan tambien al suelo por mas que la resolucion, y (c) lo hacen
contra un `b` que no es ~0: **0 de 36**.

🔴 **Y lo que esas celdas NO demuestran.** §1 midio que la correlacion con el
delta real es cero y la concordancia de signo es 0.5. Una `s` positiva en `fid`
o `pds_cosine` con signo al azar no es habilidad especifica del blanco: es que
esas dos metricas premian propiedades que una emision con dispersion y con
genes significativos satisface sin acertar CUALES. `fid` es
`k/max(n_pred, n_conf)`: el `b` tiene `n_pred = 0` en A549 y HAP1 y cobra 0
exacto, asi que cualquier emision que prediga algo con una fraccion de acierto
no nula lo supera. **Superar a un comparador que vale 0 por no predecir nada no
es evidencia de habilidad**, y es exactamente la lectura que `baseline.py:6-10`
advierte en el otro sentido. La refutacion de §1 manda sobre esta tabla.

## 5. ⏸️ ¿Cruzo el clamp de `mse`?

**No. En ninguna de las seis lineas.** `u` sigue muy por encima de `b` (sentido
bajo = peor), asi que la celda sigue devolviendo exactamente 0 y sigue ciega.

| linea | los 8 u crudos | u medio | sd_u | b | r | cruzo (u < b) | u del suelo | candidato - suelo | s |
|---|---|---|---|---|---|---|---|---|---|
| A549 | 2.537 / 2.544 / 2.617 / 2.585 / 2.470 / 2.428 / 2.695 / 2.606 | 2.560261 | 0.084864 | 0.959928 | 0.138661 | **no** | 2.543222 | +0.0170 | +0.0000 |
| BXPC3 | 1.394 / 1.438 / 1.402 / 1.351 / 1.380 / 1.391 / 1.395 / 1.402 | 1.394165 | 0.024439 | 0.824367 | 0.068734 | **no** | 1.350893 | +0.0433 | +0.0000 |
| HAP1 | 1.780 / 1.825 / 1.734 / 1.891 / 1.862 / 1.881 / 1.764 / 1.819 | 1.819534 | 0.056847 | 1.003465 | 0.185635 | **no** | 1.834121 | -0.0146 | +0.0000 |
| HT29 | 2.027 / 1.999 / 1.939 / 1.981 / 2.026 / 1.915 / 2.013 / 2.014 | 1.989164 | 0.041767 | 1.137395 | 0.132822 | **no** | 1.953467 | +0.0357 | +0.0000 |
| K562 | 3.073 / 3.120 / 3.272 / 3.188 / 3.059 / 3.158 / 3.147 / 2.964 | 3.122613 | 0.092627 | 1.173286 | 0.126447 | **no** | 3.052408 | +0.0702 | +0.0000 |
| MCF7 | 1.820 / 1.866 / 1.878 / 1.806 / 1.797 / 1.802 / 1.842 / 1.819 | 1.828790 | 0.030401 | 0.922312 | 0.104958 | **no** | 1.757116 | +0.0717 | +0.0000 |

La columna `candidato - suelo` dice algo que conviene no pasar por alto: el
candidato esta **peor** que el suelo en `mse` en las seis lineas (`mse` es de
sentido bajo). Mover genes sin acertar cuales empeora el error cuadratico del
pseudobulk respecto a no moverlos. La celda no lo registra -- el clamp devuelve 0
en los dos casos -- pero el `u` crudo si, y es la misma ceguera que
`dispersion_u.md` §5.1 declaro: dos emisiones distintas, las dos por debajo de `b`,
son indistinguibles.

## 6. 🔴 Las dos cohortes moviles, y la tasa a la que mueve `E`

Esto cierra el hueco de §3.17. El suelo predecia genes significativos en 0-3
perturbaciones por corrida, asi que su `E` se movia por accidente, una unidad cada
varias semillas. Este candidato predice mucho mas, y aqui se ve a que tasa.

| linea | E del candidato por semilla | valores distintos de E | cohorte de fid por semilla | valores distintos de fid | perts con algun gen predicho sig | suma n_DE_pred | suma n_DE_real | suma n_DE_interseccion | E del SUELO (dispersion_u.md) | n_DE_pred del suelo |
|---|---|---|---|---|---|---|---|---|---|---|
| A549 | 13/13/14/14/15/13/13/11 | 4 | 20/20/19/19/18/20/20/22 | 4 | 4-7 | 5-16 | 2537-2537 | 0-0 | 15/15/15/15/15/14/15/15 | 0-1 |
| BXPC3 | 0/0/0/0/0/0/0/0 | 1 | 33/33/33/33/33/33/33/33 | 1 | 11-14 | 53-81 | 4642-4642 | 6-16 | 0/0/0/0/0/0/0/0 | 0-3 |
| HAP1 | 13/13/12/13/12/13/12/12 | 2 | 20/20/21/20/21/20/21/21 | 2 | 2-5 | 3-7 | 817-817 | 0-2 | 13/14/14/14/14/13/14/13 | 0-2 |
| HT29 | 5/4/4/5/4/6/5/5 | 3 | 28/29/29/28/29/27/28/28 | 3 | 10-13 | 54-86 | 2505-2505 | 10-23 | 6/5/6/6/6/6/6/6 | 0-2 |
| K562 | 7/7/8/7/8/7/7/8 | 2 | 26/26/25/26/25/26/26/25 | 2 | 3-5 | 9-19 | 1042-1042 | 1-6 | 8/8/8/7/8/8/8/8 | 0-1 |
| MCF7 | 0/0/0/0/0/0/0/0 | 1 | 33/33/33/33/33/33/33/33 | 1 | 10-13 | 45-95 | 12867-12867 | 0-3 | 0/0/0/0/0/0/0/0 | 0-1 |

**La tasa, dicha con los numeros delante.** El suelo sumaba 0-3 genes predichos
significativos en todo el panel; este candidato suma decenas. Y el efecto sobre
`E` es inmediato: donde el suelo tenia `E` casi fija, el candidato la baja y la
mueve. La interseccion con el lado real es la columna que dice si eso compra algo.

🔴 **A549 es el caso que lo paga entero.** Su span en `jac` es 0.0187, menos que una
unidad de convencion (`1/33 = 0.0303`), asi que cada perturbacion que pierde su
`1.0` gratis cuesta 1.62 puntos de `s`. El candidato, prediciendo genes, pierde
varias: su `jac` sale -2.8346 con `sd_s` 1.8869, y arrastra el `avg_score` de la linea a
-0.4991. Es el mecanismo que `dispersion_u.md` §5.2 aislo con un solo gen
(STAT5A, 1 de 8 semillas); aqui se ve a escala. **Predecir señal en una celda con
span pequeno es caro con independencia de si la señal acierta.**

## 7. La resolucion del candidato

Para comparar DOS candidatos en la misma celda del mismo bundle, `r` y `b` se
cancelan (`s1 - s2 = (u1 - u2)/(r - b)`) y solo entra `sd_s`. Con `k` replicas:

```
resolucion_k = sd_s * sqrt(2/k)
```

Con la `sd_s` **del candidato**, que es lo que el punto 1 pedia y que no se asume
igual a la del suelo -- y no lo es:

| linea | metrica | sd_s del CANDIDATO | sd_s del suelo | razon | resolucion k=1 | resolucion k=2 | resolucion k=4 | resolucion k=8 | margen absoluto de §4 |
|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | 0.218700 | 0.284040 | 0.77x | 0.309288 | 0.218700 | 0.154644 | 0.109350 | **0.238823** |
| A549 | `mse_capped_norm` | 0.000000 | 0.000000 | n/d | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.074293** |
| A549 | `lfc_nmae` | 0.072338 | 0.042320 | 1.71x | 0.102302 | 0.072338 | 0.051151 | 0.036169 | **0.099006** |
| A549 | `fid` | 0.064967 | 0.000000 | n/d | 0.091878 | 0.064967 | 0.045939 | 0.032484 | **0.094128** |
| A549 | `reach` | 0.233685 | 0.157531 | 1.48x | 0.330481 | 0.233685 | 0.165240 | 0.116843 | **0.326047** |
| A549 | `jac` | 1.886945 | 0.572666 | 3.30x | 2.668544 | 1.886945 | 1.334272 | 0.943473 | **3.742886** |
| BXPC3 | `pds_cosine` | 0.114144 | 0.089199 | 1.28x | 0.161424 | 0.114144 | 0.080712 | 0.057072 | **0.122425** |
| BXPC3 | `mse_capped_norm` | 0.000000 | 0.000000 | n/d | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.051808** |
| BXPC3 | `lfc_nmae` | 0.019515 | 0.024938 | 0.78x | 0.027599 | 0.019515 | 0.013799 | 0.009758 | **0.070656** |
| BXPC3 | `fid` | 0.013687 | 0.002910 | 4.70x | 0.019356 | 0.013687 | 0.009678 | 0.006843 | **0.014529** |
| BXPC3 | `reach` | 0.054354 | 0.048389 | 1.12x | 0.076868 | 0.054354 | 0.038434 | 0.027177 | **0.197230** |
| BXPC3 | `jac` | 0.003491 | 0.000000 | n/d | 0.004937 | 0.003491 | 0.002469 | 0.001746 | **0.003617** |
| HAP1 | `pds_cosine` | 0.173963 | 0.152589 | 1.14x | 0.246021 | 0.173963 | 0.123010 | 0.086981 | **0.187153** |
| HAP1 | `mse_capped_norm` | 0.000000 | 0.000000 | n/d | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.079353** |
| HAP1 | `lfc_nmae` | 0.058245 | 0.063788 | 0.91x | 0.082371 | 0.058245 | 0.041185 | 0.029122 | **0.068134** |
| HAP1 | `fid` | 0.059275 | 0.032694 | 1.81x | 0.083828 | 0.059275 | 0.041914 | 0.029638 | **0.070618** |
| HAP1 | `reach` | 0.207686 | 0.145648 | 1.43x | 0.293712 | 0.207686 | 0.146856 | 0.103843 | **0.353813** |
| HAP1 | `jac` | 0.370632 | 0.357138 | 1.04x | 0.524153 | 0.370632 | 0.262076 | 0.185316 | **1.161641** |
| HT29 | `pds_cosine` | 0.161389 | 0.177889 | 0.91x | 0.228238 | 0.161389 | 0.114119 | 0.080694 | **0.207865** |
| HT29 | `mse_capped_norm` | 0.000000 | 0.000000 | n/d | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.035282** |
| HT29 | `lfc_nmae` | 0.032348 | 0.079919 | 0.40x | 0.045746 | 0.032348 | 0.022873 | 0.016174 | **0.096073** |
| HT29 | `fid` | 0.068094 | 0.034050 | 2.00x | 0.096299 | 0.068094 | 0.048150 | 0.034047 | **0.145465** |
| HT29 | `reach` | 0.180395 | 0.161456 | 1.12x | 0.255117 | 0.180395 | 0.127558 | 0.090197 | **0.240238** |
| HT29 | `jac` | 0.117586 | 0.059183 | 1.99x | 0.166291 | 0.117586 | 0.083146 | 0.058793 | **0.182189** |
| K562 | `pds_cosine` | 0.194544 | 0.146438 | 1.33x | 0.275127 | 0.194544 | 0.137564 | 0.097272 | **0.224820** |
| K562 | `mse_capped_norm` | 0.000000 | 0.000000 | n/d | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.040365** |
| K562 | `lfc_nmae` | 0.051443 | 0.046277 | 1.11x | 0.072752 | 0.051443 | 0.036376 | 0.025722 | **0.097985** |
| K562 | `fid` | 0.039836 | 0.031276 | 1.27x | 0.056337 | 0.039836 | 0.028168 | 0.019918 | **0.039878** |
| K562 | `reach` | 0.173582 | 0.152003 | 1.14x | 0.245482 | 0.173582 | 0.122741 | 0.086791 | **0.284161** |
| K562 | `jac` | 0.056128 | 0.038556 | 1.46x | 0.079377 | 0.056128 | 0.039688 | 0.028064 | **0.056364** |
| MCF7 | `pds_cosine` | 0.128961 | 0.126699 | 1.02x | 0.182378 | 0.128961 | 0.091189 | 0.064480 | **0.146087** |
| MCF7 | `mse_capped_norm` | 0.000000 | 0.000000 | n/d | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.039528** |
| MCF7 | `lfc_nmae` | 0.049054 | 0.036243 | 1.35x | 0.069373 | 0.049054 | 0.034687 | 0.024527 | **0.088527** |
| MCF7 | `fid` | 0.009673 | 0.000126 | 76.81x | 0.013679 | 0.009673 | 0.006840 | 0.004836 | **0.009971** |
| MCF7 | `reach` | 0.015690 | 0.015129 | 1.04x | 0.022189 | 0.015690 | 0.011095 | 0.007845 | **0.043828** |
| MCF7 | `jac` | 0.003054 | 0.000000 | n/d | 0.004319 | 0.003054 | 0.002159 | 0.001527 | **0.003748** |

**La `sd_s` del candidato NO es la del suelo**: la razon va de 0.40x a 76.81x en las 27 celdas
comparables. Asumir la del suelo habria dado un margen equivocado en los dos
sentidos segun la celda. Es la razon por la que el punto 1 pedia medirla.

## 8. La `sd` acumulada k=2..8, y donde es Bernoulli

En el orden `j = 0..7` de la familia de semillas, fijado antes de medir.

| linea | metrica | k=2 | k=3 | k=4 | k=5 | k=6 | k=7 | k=8 | arranca en 0 y salta |
|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | 0.133416 | 0.120599 | 0.129454 | 0.117976 | 0.142579 | 0.208844 | 0.218700 | no |
| A549 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| A549 | `lfc_nmae` | 0.123048 | 0.107519 | 0.087844 | 0.076181 | 0.069035 | 0.072969 | 0.072338 | no |
| A549 | `fid` | 0.032434 | 0.080219 | 0.066755 | 0.058060 | 0.063643 | 0.068306 | 0.064967 | no |
| A549 | `reach` | 0.095863 | 0.087769 | 0.215175 | 0.194716 | 0.175316 | 0.168999 | 0.233685 | no |
| A549 | `jac` | 0 | 0.935160 | 0.935160 | 1.355176 | 1.322516 | 1.274408 | 1.886945 | 🔴 **si** |
| BXPC3 | `pds_cosine` | 0.026572 | 0.156999 | 0.129419 | 0.124200 | 0.111782 | 0.121885 | 0.114144 | no |
| BXPC3 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| BXPC3 | `lfc_nmae` | 0.011049 | 0.031694 | 0.026034 | 0.023613 | 0.021122 | 0.020697 | 0.019515 | no |
| BXPC3 | `fid` | 0.005900 | 0.007226 | 0.005906 | 0.005139 | 0.012589 | 0.014269 | 0.013687 | no |
| BXPC3 | `reach` | 0.064103 | 0.048471 | 0.039634 | 0.044289 | 0.046788 | 0.057646 | 0.054354 | no |
| BXPC3 | `jac` | 0.001273 | 0.002578 | 0.002318 | 0.002009 | 0.001881 | 0.003505 | 0.003491 | no |
| HAP1 | `pds_cosine` | 0.144540 | 0.107764 | 0.185777 | 0.175966 | 0.186075 | 0.185502 | 0.173963 | no |
| HAP1 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| HAP1 | `lfc_nmae` | 0.014348 | 0.055054 | 0.047132 | 0.041192 | 0.043673 | 0.043336 | 0.058245 | no |
| HAP1 | `fid` | 0.054131 | 0.044471 | 0.038372 | 0.035983 | 0.033164 | 0.055482 | 0.059275 | no |
| HAP1 | `reach` | 0.138288 | 0.102010 | 0.288444 | 0.250528 | 0.227799 | 0.217819 | 0.207686 | no |
| HAP1 | `jac` | 0 | 0.400328 | 0.347248 | 0.380390 | 0.358067 | 0.370633 | 0.370632 | 🔴 **si** |
| HT29 | `pds_cosine` | 0.061818 | 0.142750 | 0.116718 | 0.117816 | 0.111615 | 0.169584 | 0.161389 | no |
| HT29 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| HT29 | `lfc_nmae` | 0.037916 | 0.026813 | 0.028408 | 0.024925 | 0.022369 | 0.021008 | 0.032348 | no |
| HT29 | `fid` | 0.093046 | 0.073740 | 0.074710 | 0.071108 | 0.078040 | 0.073014 | 0.068094 | no |
| HT29 | `reach` | 0.019841 | 0.103568 | 0.085569 | 0.076289 | 0.074885 | 0.102977 | 0.180395 | no |
| HT29 | `jac` | 0.115336 | 0.094040 | 0.094971 | 0.090085 | 0.135765 | 0.125778 | 0.117586 | no |
| K562 | `pds_cosine` | 0.104247 | 0.168229 | 0.137527 | 0.119424 | 0.195646 | 0.208249 | 0.194544 | no |
| K562 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| K562 | `lfc_nmae` | 0.070182 | 0.073409 | 0.065826 | 0.059055 | 0.057243 | 0.052605 | 0.051443 | no |
| K562 | `fid` | 0.006690 | 0.064470 | 0.055924 | 0.048621 | 0.044410 | 0.042390 | 0.039836 | no |
| K562 | `reach` | 0.023386 | 0.045314 | 0.097730 | 0.138473 | 0.127574 | 0.149874 | 0.173582 | no |
| K562 | `jac` | 0.004673 | 0.063092 | 0.054104 | 0.058761 | 0.054935 | 0.052280 | 0.056128 | no |
| MCF7 | `pds_cosine` | 0.058848 | 0.184311 | 0.175430 | 0.152031 | 0.136252 | 0.131891 | 0.128961 | no |
| MCF7 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | no |
| MCF7 | `lfc_nmae` | 0.031117 | 0.064626 | 0.069622 | 0.060750 | 0.056207 | 0.052967 | 0.049054 | no |
| MCF7 | `fid` | 0.007749 | 0.006041 | 0.005206 | 0.005082 | 0.004826 | 0.007288 | 0.009673 | no |
| MCF7 | `reach` | 0.013567 | 0.013030 | 0.022687 | 0.019656 | 0.017673 | 0.016597 | 0.015690 | no |
| MCF7 | `jac` | 0.002171 | 0.002894 | 0.002376 | 0.002391 | 0.002690 | 0.003254 | 0.003054 | no |

**2 celdas arrancan en 0 y saltan**: A549/`jac`, HAP1/`jac`.

🔴 **Y en A549/`jac` la `sd` no describe una nube ni una Bernoulli: describe una
RETICULA, y se puede verificar con aritmetica exacta.** Sus 8 valores de `s` son

| j | semilla | E del candidato | unidades perdidas vs el suelo (E=15) | u de jac | u en 33avos | s medido | -(unidades) x (1/33)/span | \|dif\| |
|---|---|---|---|---|---|---|---|---|
| 0 | 20266000 | 13 | 2 | 0.39394 | 13/33 | -3.2395 | -3.2395 | 0.0e+00 |
| 1 | 20266001 | 13 | 2 | 0.39394 | 13/33 | -3.2395 | -3.2395 | 0.0e+00 |
| 2 | 20266002 | 14 | 1 | 0.42424 | 14/33 | -1.6197 | -1.6197 | 1.6e-15 |
| 3 | 20266003 | 14 | 1 | 0.42424 | 14/33 | -1.6197 | -1.6197 | 1.6e-15 |
| 4 | 20266004 | 15 | 0 | 0.45455 | 15/33 | +0.0000 | +0.0000 | 0.0e+00 |
| 5 | 20266005 | 13 | 2 | 0.39394 | 13/33 | -3.2395 | -3.2395 | 0.0e+00 |
| 6 | 20266006 | 13 | 2 | 0.39394 | 13/33 | -3.2395 | -3.2395 | 0.0e+00 |
| 7 | 20266007 | 11 | 4 | 0.33333 | 11/33 | -6.4790 | -6.4790 | 0.0e+00 |

`(1/33)/span = (1/33)/0.018709 = 1.6197`, y los 8 `s` son multiplos exactos de esa cantidad con residuo < 1e-15. La celda **solo puede
tomar valores en esa reticula**: no hay nada continuo que promediar. Pedir
mas semillas ahi no estrecha un intervalo, estima la distribucion de cuantas
unidades de convencion pierde la emision. Es la misma aritmetica que
`dispersion_u.md` §5.2 verifico con un solo gen; con un candidato que predice
señal se ve en las 8 semillas a la vez.

**6 celdas con `sd_s = 0`**, todas de `mse_capped_norm`: el clamp (§5).
No son estabilidad, son ceguera; el `u` crudo se mueve (`sd_u` 0.02-0.07) y el `s`
no.

## 9. Las 20 filas trabajadas

**Criterio de seleccion, declarado antes de mirar los valores**: las celdas de
`pds_cosine`, `fid` y `jac` en las seis lineas = 18 filas, mas `lfc_nmae` en A549
y MCF7, los dos extremos del `span` de esa metrica = **20**. Se eligen `pds`, `fid`
y `jac` porque son las tres donde §4 encuentra celdas que superan al `b`, que es la
pregunta del encargo; no por su valor.

Cada fila trae los 8 `u` crudos, los 8 `s`, la aritmetica `(u-b)/(r-b)` rehecha a
mano sobre la media, y los tres terminos del margen.

| linea | metrica | los 8 u crudos | u medio | sd_u | b | r | span r-b | s a mano desde u medio | s medio medido | \|dif\| | sd_s | termino de r | termino de b | margen absoluto | resolucion k=8 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | 0.5805 / 0.5256 / 0.5152 / 0.4915 / 0.5521 / 0.6013 / 0.4186 / 0.4422 | 0.515862 | 0.063663 | 0.486742 | 0.777841 | 0.291098 | +0.1000 | +0.1000 | 2.8e-17 | 0.218700 | 0.004462 | 0.095849 | **0.238823** | 0.109350 |
| A549 | `fid` | 0.1192 / 0.1000 / 0.0539 / 0.0803 / 0.0833 / 0.1250 / 0.1333 / 0.1169 | 0.101491 | 0.027147 | 0.000000 | 0.417863 | 0.417863 | +0.2429 | +0.2429 | 2.8e-17 | 0.064967 | 0.009612 | 0.067431 | **0.094128** | 0.032484 |
| A549 | `jac` | 0.3939 / 0.3939 / 0.4242 / 0.4242 / 0.4545 / 0.3939 / 0.3939 / 0.3333 | 0.401515 | 0.035302 | 0.454545 | 0.473254 | 0.018709 | -2.8346 | -2.8346 | 4.4e-16 | 1.886945 | 2.372039 | 2.195919 | **3.742886** | 0.943473 |
| BXPC3 | `pds_cosine` | 0.4725 / 0.4545 / 0.5928 / 0.5237 / 0.4536 / 0.4848 / 0.5814 / 0.4858 | 0.506155 | 0.054651 | 0.476326 | 0.955114 | 0.478788 | +0.0623 | +0.0623 | 9.7e-17 | 0.114144 | 0.000979 | 0.044250 | **0.122425** | 0.057072 |
| BXPC3 | `fid` | 0.0657 / 0.0607 / 0.0693 / 0.0656 / 0.0660 / 0.0483 / 0.0760 / 0.0706 | 0.065280 | 0.008182 | 0.001444 | 0.599253 | 0.597809 | +0.1068 | +0.1068 | 0.0e+00 | 0.013687 | 0.003658 | 0.003223 | **0.014529** | 0.006843 |
| BXPC3 | `jac` | 0.0013 / 0.0009 / 0.0019 / 0.0010 / 0.0013 / 0.0010 / 0.0029 / 0.0008 | 0.001392 | 0.000697 | 0.000978 | 0.200769 | 0.199792 | +0.0021 | +0.0021 | 4.3e-19 | 0.003491 | 0.000106 | 0.000940 | **0.003617** | 0.001746 |
| HAP1 | `pds_cosine` | 0.6477 / 0.5758 / 0.5909 / 0.4896 / 0.5199 / 0.4792 / 0.4811 / 0.5682 | 0.544034 | 0.061249 | 0.437500 | 0.789583 | 0.352083 | +0.3026 | +0.3026 | 0.0e+00 | 0.173963 | 0.005861 | 0.068766 | **0.187153** | 0.086981 |
| HAP1 | `fid` | 0.0424 / 0.0005 / 0.0000 / 0.0007 / 0.0278 / 0.0250 / 0.0833 / 0.0714 | 0.031394 | 0.032414 | 0.000000 | 0.546846 | 0.546846 | +0.0574 | +0.0574 | 0.0e+00 | 0.059275 | 0.004427 | 0.038128 | **0.070618** | 0.029638 |
| HAP1 | `jac` | 0.3941 / 0.3941 / 0.3636 / 0.3942 / 0.3636 / 0.3939 / 0.3636 / 0.3636 | 0.378861 | 0.016276 | 0.424242 | 0.468156 | 0.043914 | -1.0334 | -1.0334 | 2.2e-16 | 0.370632 | 0.888896 | 0.649543 | **1.161641** | 0.185316 |
| HT29 | `pds_cosine` | 0.4706 / 0.4953 / 0.5492 / 0.5085 / 0.4678 / 0.5237 / 0.4015 / 0.4583 | 0.484375 | 0.045452 | 0.501894 | 0.783523 | 0.281629 | -0.0622 | -0.0622 | 4.2e-17 | 0.161389 | 0.001687 | 0.130991 | **0.207865** | 0.080694 |
| HT29 | `fid` | 0.1008 / 0.1426 / 0.1400 / 0.0997 / 0.0998 / 0.0814 / 0.0973 / 0.1015 | 0.107903 | 0.021627 | 0.114189 | 0.431793 | 0.317604 | -0.0198 | -0.0198 | 3.1e-17 | 0.068094 | 0.001824 | 0.128530 | **0.145465** | 0.034047 |
| HT29 | `jac` | 0.1522 / 0.1226 / 0.1227 / 0.1527 / 0.1227 / 0.1830 / 0.1529 / 0.1525 | 0.145181 | 0.021286 | 0.091090 | 0.272116 | 0.181026 | +0.2988 | +0.2988 | 0.0e+00 | 0.117586 | 0.043166 | 0.132299 | **0.182189** | 0.058793 |
| K562 | `pds_cosine` | 0.5199 / 0.5644 / 0.4631 / 0.5199 / 0.5227 / 0.3968 / 0.5833 / 0.5322 | 0.512784 | 0.058732 | 0.520833 | 0.822727 | 0.301894 | -0.0267 | -0.0267 | 4.2e-17 | 0.194544 | 0.000693 | 0.112678 | **0.224820** | 0.097272 |
| K562 | `fid` | 0.0809 / 0.0850 / 0.0345 / 0.0504 / 0.0585 / 0.0714 / 0.0492 / 0.0530 | 0.060351 | 0.017320 | 0.076923 | 0.511705 | 0.434782 | -0.0381 | -0.0381 | 6.9e-18 | 0.039836 | 0.001755 | 0.000528 | **0.039878** | 0.019918 |
| K562 | `jac` | 0.2126 / 0.2144 / 0.2438 / 0.2144 / 0.2433 / 0.2148 / 0.2130 / 0.2447 | 0.225142 | 0.015596 | 0.212121 | 0.489996 | 0.277874 | +0.0469 | +0.0469 | 1.4e-17 | 0.056128 | 0.005120 | 0.000574 | **0.056364** | 0.028064 |
| MCF7 | `pds_cosine` | 0.5862 / 0.5502 / 0.4337 / 0.6013 / 0.5483 / 0.5530 / 0.4953 / 0.5890 | 0.544626 | 0.055761 | 0.503788 | 0.936174 | 0.432386 | +0.0944 | +0.0944 | 1.2e-16 | 0.128961 | 0.001596 | 0.068615 | **0.146087** | 0.064480 |
| MCF7 | `fid` | 0.0089 / 0.0038 / 0.0084 / 0.0055 / 0.0042 / 0.0080 / 0.0136 / 0.0166 | 0.008621 | 0.004473 | 0.000173 | 0.462605 | 0.462432 | +0.0183 | +0.0183 | 3.5e-18 | 0.009673 | 0.001596 | 0.001822 | **0.009971** | 0.004836 |
| MCF7 | `jac` | 0.0005 / 0.0002 / 0.0000 / 0.0003 / 0.0004 / 0.0000 / 0.0006 / 0.0004 | 0.000312 | 0.000221 | 0.000029 | 0.072243 | 0.072214 | +0.0039 | +0.0039 | 8.7e-19 | 0.003054 | 0.000588 | 0.002091 | **0.003748** | 0.001527 |
| A549 | `lfc_nmae` | 0.9812 / 1.0510 / 0.9722 / 1.0040 / 1.0057 / 0.9919 / 1.0400 / 0.9773 | 1.002918 | 0.028997 | 0.856686 | 0.455829 | -0.400857 | -0.3648 | -0.3648 | 1.7e-16 | 0.072338 | 0.009308 | 0.066953 | **0.099006** | 0.036169 |
| MCF7 | `lfc_nmae` | 1.0394 / 1.0486 / 1.0659 / 1.0324 / 1.0500 / 1.0399 / 1.0533 / 1.0478 | 1.047187 | 0.010222 | 0.872196 | 0.663824 | -0.208372 | -0.8398 | -0.8398 | 6.7e-16 | 0.049054 | 0.020894 | 0.070670 | **0.088527** | 0.024527 |

Son **20** filas.

## Lo que no se hizo

- **No se leyo nada del lado real perturbado para construir la prediccion.** El
  arnes lo impide estructuralmente y se verifico en su codigo (§2). No hay fallo
  del arnes que reportar.
- **No se uso ninguna cifra del `b` ni del anchor como entrada.** Entran solo en la
  puntuacion, que es su sitio.
- **No se ajusto el candidato mirando su `s`.** Las dos magnitudes descartadas
  (§1.1) se descartaron por propiedades medidas de la EMISION y nunca se puntuaron.
  La version que corre se midio una vez, con 8 semillas, y no se toco despues.
- **No se itero.** Una iteracion seria un candidato nuevo con su propia medicion.
- **No se proponen umbrales ni se decide adopcion.** §4 aplica la regla ya escrita
  en §3.16 y §3.18; no se añade ninguna.
- No se leyo nada de `~/vcc/benchmark/` ni de
  `/datos/benchmark_VEDADO_techo_metricas/`.
- Sin envios, sin `vcc sample`, sin commits.

