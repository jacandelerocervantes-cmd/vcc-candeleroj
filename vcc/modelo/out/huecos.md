# Los cuatro huecos de §3.14

Cuatro preguntas independientes, medidas, no inferidas. Ningun criterio de
adopcion se redacta aqui: no es de este encargo.

**Fuentes de medicion** (todas de esta sesion, todas sobre la MISMA poblacion de
siempre -- `build_fold_adatas` sin modificar, `seed_real = 20260001 + fold_idx`,
panel comun de 33 de `suelo_panel.json`, etiqueta de control `non-targeting`):

| fuente | que contiene | como se obtuvo |
|---|---|---|
| `descomp_rgen.json` | las 5 particiones del anchor recomputadas, **las 6 lineas** | `cell_eval2.ceiling._disjoint_halves` sin modificar + las semillas que estampo cada bundle |
| `huecos_cohorte.json` | lado real por perturbacion (`n_DE_real`, `n_conf`, `q`) y valores por perturbacion de los 6 metricas, para nuestro suelo y para el `b` oficial | `cell-eval2 run --write-degenes` + `cell_eval2.metrics.direction._reference_stats` sin reimplementar |
| `huecos_pds.json` | distribucion de `pds_cosine` por perturbacion | `results.csv` / `baseline_results.csv` ya en disco |
| `anclas_nativas_comp.json` | `b` y `r` oficiales de las 6 metricas por linea | los 6 bundles de competencia de `anclas_nativas_v2.md` |

**Criterios de seleccion, declarados antes de medir.** Las SEIS lineas. Los 33 del
panel comun. Ninguna celda excluida por su valor ni por su signo. Donde una tabla
muestra menos de 33 filas es porque la metrica no emitio la fila (gate del lado
real), y el conteo de cuantas faltan ES el resultado, no un recorte.

**Una desviacion instrumental, declarada.** Para obtener las tablas DE del `b`
oficial hubo que volver a correr su emision por `cell-eval2 run`, y `run` la rechaza
con `declared input_type='counts' but values are fractional` (`norm.py:482`): el
perfil del baseline es una MEDIA, fraccionaria en cualquier espacio. Se paso
`--set allow_fractional_counts=true`, que es **lo que `cell-eval2 baseline` se pone
a si mismo** (`baseline.py:1320-1321`) y que esta en `DIGEST_EXEMPT_FIELDS`
(`baseline.py:692-703`: *"a validation ALLOWANCE, not a scoring semantic"*). Aplica
al lado PRED solamente (`run.py:397-401`), el lado real sigue validado estricto. No
cambia ninguna semantica de puntuacion; la verificacion de abajo lo comprueba
numero a numero.

## Los cuatro veredictos, en una linea cada uno

| hueco | la pregunta | el veredicto |
|---|---|---|
| **1** | ¿la convencion resta del span, o eso pasa solo donde `E` es grande? | **No se generaliza.** Resta en 2 de 6 (A549, HAP1) y suma en 4. El signo lo decide cuanta senal DE tiene el panel, no el tamano de `E`. |
| **2** | ¿la cohorte del lado real es la misma para los tres objetos? | **Si entre submissions** (medido celda por celda, 6 lineas), **no entre poblaciones**: las particiones del anchor corren sobre media poblacion y su cohorte de `fid` y `reach` se mueve. Y la cohorte que la componente genuina de `jac` necesita -- union no vacia -- NO es la del lado real y si la mueve la submission. |
| **3** | ¿nuestro `b` de `reach` cae sobre nuestro `c/N_conf`? | **No, y por razones distintas en cada punta.** 4 lineas por encima (hasta 9.9x), K562 encima, A549 un 6% por debajo. El `b` esta hecho en un 53-88% de 4 a 7 perturbaciones con `N_conf <= 10` que sacan 1.0 exacto. |
| **4** | ¿explica el empate en 0.5 el `b < u` de HAP1 y BXPC3? | **No. Hipotesis refutada.** Ningun delta tiene norma cero, el `b` tiene mas 0.5 exactos que el suelo, y 0.5 esta en la rejilla asi que ese conteo no mide empates. Lo que explica las seis diferencias es el `SE` propio de una media de 33 rangos casi uniformes, 0.0518: HAP1 queda en 0.97 sigma. |

**Verificacion de que el `b` por perturbacion de este reporte ES el `b` oficial.**
Dos comprobaciones encadenadas, las dos celda por celda:

1. la emision del baseline regenerada con `--seed 0` reproduce el
   `baseline_results.csv` del `BASE_<linea>` que construyo los bundles;
2. `cell-eval2 run` sobre esa emision reproduce los mismos valores por
   perturbacion que `cell-eval2 baseline`.

| linea | celdas (1) baseline regenerado vs oficial | difieren | max \|dif\| | celdas (2) run vs baseline | difieren | max \|dif\| |
|---|---|---|---|---|---|---|
| A549 | 272 | 2 | 3.33e-16 | 110 | 2 | 2.22e-16 |
| BXPC3 | 287 | 13 | 1.22e-15 | 155 | 0 | 0 |
| HAP1 | 272 | 0 | 0 | 112 | 0 | 0 |
| HT29 | 279 | 0 | 0 | 138 | 0 | 0 |
| K562 | 273 | 0 | 0 | 126 | 5 | 3.33e-16 |
| MCF7 | 297 | 10 | 9.99e-16 | 165 | 24 | 8.88e-16 |

Ninguna diferencia pasa de **1.22e-15**. No son diferencias: son los 2-3 ULP
de coma flotante que `anclas_nativas_v2.md` ya habia aislado en el PASO 0, y que
vuelven a salir **en las mismas celdas**:

| linea | perturbacion | metrica | oficial | regenerado | dif |
|---|---|---|---|---|---|
| A549 | IFNAR2 | de_wilcoxon_lfc_nmae | 0.8561795666261268 | 0.8561795666261265 | -3.331e-16 |
| A549 | TYK2 | de_wilcoxon_lfc_nmae | 0.8391460433144295 | 0.8391460433144297 | +1.110e-16 |
| BXPC3 | H1F0 | de_wilcoxon_lfc_nmae | 0.34696050258406996 | 0.34696050258406985 | -1.110e-16 |
| BXPC3 | ID2 | de_wilcoxon_lfc_nmae | 0.8669434558170855 | 0.8669434558170854 | -1.110e-16 |
| BXPC3 | IFNAR2 | de_wilcoxon_lfc_nmae | 0.8675108281990093 | 0.867510828199008 | -1.221e-15 |
| BXPC3 | JAK1 | de_wilcoxon_lfc_nmae | 0.7898088406774896 | 0.7898088406774894 | -2.220e-16 |
| BXPC3 | JUN | de_wilcoxon_lfc_nmae | 0.793852223825327 | 0.7938522238253269 | -1.110e-16 |
| BXPC3 | MAP3K14 | de_wilcoxon_lfc_nmae | 0.8023020887196471 | 0.8023020887196473 | +1.110e-16 |
| BXPC3 | MAPK8 | de_wilcoxon_lfc_nmae | 0.9656100202965551 | 0.9656100202965553 | +2.220e-16 |
| BXPC3 | RNF114 | de_wilcoxon_lfc_nmae | 0.7745536317928219 | 0.7745536317928222 | +3.331e-16 |
| BXPC3 | SOCS1 | de_wilcoxon_lfc_nmae | 0.8079679563964015 | 0.8079679563964016 | +1.110e-16 |
| BXPC3 | STAT1 | de_wilcoxon_lfc_nmae | 0.7966633817146717 | 0.7966633817146722 | +5.551e-16 |
| BXPC3 | STAT3 | de_wilcoxon_lfc_nmae | 0.8517714322883864 | 0.8517714322883863 | -1.110e-16 |
| BXPC3 | STAT4 | de_wilcoxon_lfc_nmae | 0.8041620807040033 | 0.8041620807040031 | -2.220e-16 |
| BXPC3 | ZBP1 | de_wilcoxon_lfc_nmae | 0.8584032614744544 | 0.8584032614744548 | +3.331e-16 |
| MCF7 | H1F0 | de_wilcoxon_lfc_nmae | 0.9949881050028376 | 0.9949881050028379 | +2.220e-16 |
| MCF7 | HERC6 | de_wilcoxon_lfc_nmae | 0.8666969373797873 | 0.8666969373797871 | -1.110e-16 |
| MCF7 | IFNAR2 | de_wilcoxon_lfc_nmae | 0.9927355699371916 | 0.9927355699371915 | -1.110e-16 |
| MCF7 | IRF1 | de_wilcoxon_lfc_nmae | 0.8795594075118879 | 0.8795594075118883 | +4.441e-16 |
| MCF7 | IRF7 | de_wilcoxon_lfc_nmae | 0.9052005115472296 | 0.9052005115472286 | -9.992e-16 |
| MCF7 | MAPK8 | de_wilcoxon_lfc_nmae | 0.9042809187438366 | 0.9042809187438362 | -4.441e-16 |
| MCF7 | STAT5A | de_wilcoxon_lfc_nmae | 0.8461249555691421 | 0.8461249555691422 | +1.110e-16 |
| MCF7 | TYK2 | de_wilcoxon_lfc_nmae | 0.9232235146465162 | 0.923223514646517 | +7.772e-16 |
| MCF7 | VAV1 | de_wilcoxon_lfc_nmae | 0.8656398923715515 | 0.8656398923715509 | -5.551e-16 |
| MCF7 | ZBP1 | de_wilcoxon_lfc_nmae | 0.8964283323723601 | 0.8964283323723595 | -5.551e-16 |

Todas en `de_wilcoxon_lfc_nmae`. Que reaparezcan las mismas perturbaciones en
una corrida nueva confirma lo que el PASO 0 habia concluido: la cadena no es
bit-reproducible corrida a corrida, y el efecto esta localizado. **No depende
de `allow_fractional_counts`**: la comprobacion (1) no lo usa en absoluto y
produce la misma desviacion que la (2).

El `b` por perturbacion de este reporte es, a 2-3 ULP, el `b` oficial, y
`allow_fractional_counts` no movio nada.

---

## Hueco 1 — la descomposicion en las cuatro lineas que faltaban

`descomposicion.md` recomputo las 5 particiones del anchor solo en A549 y HAP1.
Aqui estan BXPC3, HT29, K562 y MCF7 con el mismo metodo, y las seis juntas.

### 1.1 Lo primero: el valor recomputado reproduce el del bundle

Antes de descomponer nada. Si la reproduccion fallara, la descomposicion no
describiria el `r` que puntua.

| linea | particiones | r crudo (media) | max \|recomputado - bundle\| | reproduce |
|---|---|---|---|---|
| A549 | 5 | 0.473254 | 1.11e-16 | si |
| BXPC3 | 5 | 0.200769 | 2.78e-17 | si |
| HAP1 | 5 | 0.468156 | 5.55e-17 | si |
| HT29 | 5 | 0.272116 | 5.55e-17 | si |
| K562 | 5 | 0.489996 | 1.11e-16 | si |
| MCF7 | 5 | 0.072243 | 1.39e-17 | si |

La desviacion maxima en las 30 particiones es **1.11e-16** -- ruido de ultimo bit
de coma flotante, no una diferencia. La descomposicion describe el `r` que puntua.

### 1.2 Las 20 particiones nuevas, fila por fila

Los intermedios crudos de los que sale todo lo demas. `E` = perturbaciones con
`|union| = 0`, que cobran 1.0 cada una por la convencion D-023.

| linea | split | semilla | n half_a | n half_b | E | n genuinos | suma genuinos | media genuinos | E/33 | suma/33 | jac recomputado | jac del bundle | dif |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| BXPC3 | 0 | 2968811710 | 8214 | 8214 | 1 | 32 | 4.3971 | 0.137409 | 0.030303 | 0.133245 | 0.163549 | 0.163549 | 2.8e-17 |
| BXPC3 | 1 | 3677149159 | 8214 | 8214 | 3 | 30 | 3.8424 | 0.128079 | 0.090909 | 0.116436 | 0.207345 | 0.207345 | 0.0e+00 |
| BXPC3 | 2 | 745650761 | 8214 | 8214 | 1 | 32 | 5.3812 | 0.168162 | 0.030303 | 0.163066 | 0.193369 | 0.193369 | 2.8e-17 |
| BXPC3 | 3 | 2884920346 | 8214 | 8214 | 3 | 30 | 4.7114 | 0.157047 | 0.090909 | 0.142770 | 0.233679 | 0.233679 | 2.8e-17 |
| BXPC3 | 4 | 2642120001 | 8214 | 8214 | 2 | 31 | 4.7948 | 0.154671 | 0.060606 | 0.145297 | 0.205903 | 0.205903 | 0.0e+00 |
| HT29 | 0 | 2968811710 | 7968 | 7968 | 8 | 25 | 2.6545 | 0.106181 | 0.242424 | 0.080440 | 0.322864 | 0.322864 | 0.0e+00 |
| HT29 | 1 | 3677149159 | 7968 | 7968 | 6 | 27 | 3.6589 | 0.135514 | 0.181818 | 0.110875 | 0.292693 | 0.292693 | 0.0e+00 |
| HT29 | 2 | 745650761 | 7968 | 7968 | 6 | 27 | 2.3101 | 0.085561 | 0.181818 | 0.070004 | 0.251823 | 0.251823 | 0.0e+00 |
| HT29 | 3 | 2884920346 | 7968 | 7968 | 8 | 25 | 2.7347 | 0.109389 | 0.242424 | 0.082871 | 0.325295 | 0.325295 | 5.6e-17 |
| HT29 | 4 | 2642120001 | 7968 | 7968 | 3 | 30 | 2.5409 | 0.084696 | 0.090909 | 0.076997 | 0.167906 | 0.167906 | 0.0e+00 |
| K562 | 0 | 2968811710 | 7299 | 7299 | 10 | 23 | 4.3033 | 0.187099 | 0.303030 | 0.130402 | 0.433433 | 0.433433 | 5.6e-17 |
| K562 | 1 | 3677149159 | 7299 | 7299 | 14 | 19 | 4.9303 | 0.259492 | 0.424242 | 0.149404 | 0.573647 | 0.573647 | 1.1e-16 |
| K562 | 2 | 745650761 | 7299 | 7299 | 10 | 23 | 4.6354 | 0.201539 | 0.303030 | 0.140467 | 0.443497 | 0.443497 | 0.0e+00 |
| K562 | 3 | 2884920346 | 7299 | 7299 | 11 | 22 | 3.0995 | 0.140886 | 0.333333 | 0.093924 | 0.427257 | 0.427257 | 5.6e-17 |
| K562 | 4 | 2642120001 | 7299 | 7299 | 15 | 18 | 3.8808 | 0.215598 | 0.454545 | 0.117599 | 0.572145 | 0.572145 | 0.0e+00 |
| MCF7 | 0 | 2968811710 | 7934 | 7934 | 0 | 33 | 1.7458 | 0.052904 | 0.000000 | 0.052904 | 0.052904 | 0.052904 | 0.0e+00 |
| MCF7 | 1 | 3677149159 | 7934 | 7934 | 2 | 31 | 1.7842 | 0.057555 | 0.060606 | 0.054067 | 0.114673 | 0.114673 | 1.4e-17 |
| MCF7 | 2 | 745650761 | 7934 | 7934 | 0 | 33 | 1.8273 | 0.055373 | 0.000000 | 0.055373 | 0.055373 | 0.055373 | 6.9e-18 |
| MCF7 | 3 | 2884920346 | 7934 | 7934 | 1 | 32 | 1.7813 | 0.055665 | 0.030303 | 0.053978 | 0.084281 | 0.084281 | 0.0e+00 |
| MCF7 | 4 | 2642120001 | 7934 | 7934 | 0 | 33 | 1.7815 | 0.053985 | 0.000000 | 0.053985 | 0.053985 | 0.053985 | 6.9e-18 |

Son **20** filas: 4 lineas x 5 particiones. La identidad
`jac = (E + suma_genuinos)/33` se lee directamente en cada una: la columna `E/33`
mas la columna `suma/33` da la columna `jac recomputado`.

### 1.3 E y la media de los genuinos, las seis lineas

| linea | E por particion | E medio | rango E | media genuinos por particion | media genuinos | r crudo | sd r crudo | componente convencion de r | componente genuina de r |
|---|---|---|---|---|---|---|---|---|---|
| A549 | 14/13/15/13/14 | 13.8 | 13-15 | 0.1157 / 0.0684 / 0.1128 / 0.0564 / 0.1244 | 0.095526 | 0.473254 | 0.035008 | 0.418182 | 0.055072 |
| BXPC3 | 1/3/1/3/2 | 2.0 | 1-3 | 0.1374 / 0.1281 / 0.1682 / 0.1570 / 0.1547 | 0.149074 | 0.200769 | 0.022772 | 0.060606 | 0.140163 |
| HAP1 | 13/14/9/13/7 | 11.2 | 7-14 | 0.2147 / 0.2480 / 0.1853 / 0.1860 / 0.1566 | 0.198131 | 0.468156 | 0.084462 | 0.339394 | 0.128762 |
| HT29 | 8/6/6/8/3 | 6.2 | 3-8 | 0.1062 / 0.1355 / 0.0856 / 0.1094 / 0.0847 | 0.104268 | 0.272116 | 0.058478 | 0.187879 | 0.084237 |
| K562 | 10/14/10/11/15 | 12.0 | 10-15 | 0.1871 / 0.2595 / 0.2015 / 0.1409 / 0.2156 | 0.200923 | 0.489996 | 0.067887 | 0.363636 | 0.126359 |
| MCF7 | 0/2/0/1/0 | 0.6 | 0-2 | 0.0529 / 0.0576 / 0.0554 / 0.0557 / 0.0540 | 0.055096 | 0.072243 | 0.024237 | 0.018182 | 0.054061 |

### 1.4 El aporte de la convencion al span, contra el aporte de lo genuino

El span crudo es `r - b`. Las dos componentes son aditivas porque la convencion
entra como `E/33` y lo genuino como `suma_genuinos/33`, y las dos suman el crudo.
Asi que:

```
span_crudo      = r_crudo      - b_crudo
aporte_conv     = E_r/33       - E_b/33
aporte_genuino  = suma_gen_r/33 - suma_gen_b/33
span_crudo      = aporte_conv  + aporte_genuino      (identidad, se verifica abajo)
```

`E_b` y `suma_gen_b` del `b` oficial salen de `huecos_cohorte.json`, que trae
`|union|` por perturbacion MEDIDO (no inferido del `jac == 1.0`).

| linea | E del b | E del r (medio) | conv b = E_b/33 | conv r = E_r/33 | genuina b | genuina r | b crudo | \|b - (conv+gen)\| | r crudo | span crudo | aporte CONVENCION | aporte GENUINO | suma de aportes | \|span - suma\| |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A549 | 15 | 13.8 | 0.454545 | 0.418182 | 0.000000 | 0.055072 | 0.454545 | 0.0e+00 | 0.473254 | 0.018709 | **-0.036364** | **+0.055072** | 0.018709 | 5.6e-17 |
| BXPC3 | 0 | 2.0 | 0.000000 | 0.060606 | 0.000978 | 0.140163 | 0.000978 | 0.0e+00 | 0.200769 | 0.199792 | **+0.060606** | **+0.139185** | 0.199792 | 2.8e-17 |
| HAP1 | 14 | 11.2 | 0.424242 | 0.339394 | 0.000000 | 0.128762 | 0.424242 | 0.0e+00 | 0.468156 | 0.043914 | **-0.084848** | **+0.128762** | 0.043914 | 1.1e-16 |
| HT29 | 3 | 6.2 | 0.090909 | 0.187879 | 0.000181 | 0.084237 | 0.091090 | 1.4e-17 | 0.272116 | 0.181026 | **+0.096970** | **+0.084057** | 0.181026 | 2.8e-17 |
| K562 | 7 | 12.0 | 0.212121 | 0.363636 | 0.000000 | 0.126359 | 0.212121 | 0.0e+00 | 0.489996 | 0.277874 | **+0.151515** | **+0.126359** | 0.277874 | 1.1e-16 |
| MCF7 | 0 | 0.6 | 0.000000 | 0.018182 | 0.000029 | 0.054061 | 0.000029 | 0.0e+00 | 0.072243 | 0.072214 | **+0.018182** | **+0.054033** | 0.072214 | 2.8e-17 |

Dos verificaciones en la misma tabla. `|b - (conv+gen)|` comprueba que el `b`
oficial del bundle se reconstruye desde los conteos por perturbacion que medimos
aqui -- es la identidad `jac = (E + suma_genuinos)/33` aplicada al `b`. Y
`|span - suma|` confirma que los dos aportes suman el span crudo en las seis
lineas.

#### Veredicto de §1.4: la afirmacion con las seis a la vista

**La convencion resta del span en 2 de 6 lineas y suma en 4.**

- **A549**: E_b = 15, E_r medio = 13.8 -> la convencion resta 0.036364 del span de 0.018709 (194.4% de su magnitud); lo genuino aporta +0.055072.
- **BXPC3**: E_b = 0, E_r medio = 2.0 -> la convencion suma 0.060606 del span de 0.199792 (30.3% de su magnitud); lo genuino aporta +0.139185.
- **HAP1**: E_b = 14, E_r medio = 11.2 -> la convencion resta 0.084848 del span de 0.043914 (193.2% de su magnitud); lo genuino aporta +0.128762.
- **HT29**: E_b = 3, E_r medio = 6.2 -> la convencion suma 0.096970 del span de 0.181026 (53.6% de su magnitud); lo genuino aporta +0.084057.
- **K562**: E_b = 7, E_r medio = 12.0 -> la convencion suma 0.151515 del span de 0.277874 (54.5% de su magnitud); lo genuino aporta +0.126359.
- **MCF7**: E_b = 0, E_r medio = 0.6 -> la convencion suma 0.018182 del span de 0.072214 (25.2% de su magnitud); lo genuino aporta +0.054033.

**La afirmacion de `descomposicion.md` -- "la convencion resta del span" -- no se
generaliza.** Estaba medida en A549 y HAP1, y en las otras cuatro lineas la
convencion **suma**. Lo que hay no es un mecanismo con un signo fijo sino una
resta de dos cantidades con origenes distintos, y con las seis a la vista se ve
cual es cual:

- **`E_b` es, salvo cuatro celdas, el complemento de la cohorte del lado real**:
  `33 - |n_conf > 0|`. Son 15, 0, 14, 6, 8 y 0 (§2.5 lo mide y nombra las cuatro
  excepciones). Mide **cuantas perturbaciones del panel no tienen senal DE en la
  referencia**. No tiene nada que ver con el baseline: nuestro suelo da el mismo
  numero.
- **`E_r` lo decide la perdida de potencia de media poblacion**: 13.8, 2.0, 11.2,
  6.2, 12.0 y 0.6. Es cuantas quedan sin senal **a los dos lados de la particion**.

Y por eso el signo cambia: donde el panel tiene poca senal (A549 con 18 de 33 en
cohorte, HAP1 con 19) `E_b` ya es alto y la particion no puede subirlo mucho mas,
asi que `E_r < E_b` y la convencion resta. Donde el panel tiene mucha senal (BXPC3
y MCF7, 33 de 33) `E_b = 0` por construccion y cualquier perdida de potencia lo
sube, asi que la convencion suma. La correlacion entre `E_b` y el aporte de la
convencion sobre las seis lineas es **-0.589**.

**No hace falta que `E` sea grande para que reste; hace falta que el panel tenga
poca senal.** Son la misma cosa en estas seis lineas, y por eso los dos casos donde
resta son exactamente los dos de cohorte mas chica.

Y la magnitud esta cuantizada: cada unidad de E mueve el crudo exactamente
`1/33 = 0.030303`. Por eso el `sd` entre particiones del crudo es grande (E
cambia de particion a particion) mientras el de la componente genuina es chico.

---

## Hueco 2 — la cohorte comun de la componente genuina

### 2.1 La hipotesis, y por que esa y no otra

> La cohorte debe definirse por el LADO REAL unicamente: las perturbaciones donde
> la referencia tiene senal DE. Razon: eso la vuelve invariante a la submission.

El codigo dice por que esa es la propiedad que importa, y lo dice dos veces:

- `metrics/de.py:608-612` (gate de `lfc_nmae`): *"the gate rules below drop a
  perturbation whose gate is empty, smaller than `min_gate_size`, or whose
  denominator is zero or non-finite -- those score nothing rather than 1.0. **All
  three are real-side decisions, so the surviving set is the same for every
  submission.**"*
- `metrics/direction.py:296-299` (`_reference_stats`, que produce el `n_conf` de
  `reach`): *"Never touches the join: `_direction_frame` is an INNER join, so
  anything computed from it is prediction-dependent, and `q`/`N_conf` must not be
  -- **that invariance is what stops a model shrinking its own BUDGET by omitting
  genes**."*

Asi que la hipotesis se prueba en dos mitades, y las dos se miden:

1. ¿es la misma entre SUBMISSIONS sobre la misma poblacion? (nuestro suelo contra
   el `b` oficial)
2. ¿es la misma entre OBJETOS, incluidas las particiones del anchor, que corren
   sobre media poblacion?

### 2.2 Mitad 1: el lado real es invariante a la submission. Medido.

Dos emisiones completamente distintas sobre el MISMO `real.h5ad` de cada linea:
nuestro suelo (control remuestreado, sin escalar) y el `b` oficial (control
escalado al perfil de respuesta promedio). Se comparan los conteos del lado real
perturbacion por perturbacion.

| linea | perturbaciones | difieren en n_DE_real | difieren en n_conf | difieren en n_DE_pred (lado predicho) | suma n_DE_pred suelo | suma n_DE_pred b | lado real identico |
|---|---|---|---|---|---|---|---|
| A549 | 33 | 0 | 0 | 0 | 0 | 0 | **si** |
| BXPC3 | 33 | 0 | 0 | 3 | 1 | 2 | **si** |
| HAP1 | 33 | 0 | 0 | 0 | 0 | 0 | **si** |
| HT29 | 33 | 0 | 0 | 12 | 3 | 14 | **si** |
| K562 | 33 | 0 | 0 | 2 | 0 | 2 | **si** |
| MCF7 | 33 | 0 | 0 | 6 | 0 | 7 | **si** |

**El lado real es identico celda por celda en las seis lineas.**

Las tres ultimas columnas estan ahi como control, y hay que leerlas con cuidado:
el lado predicho difiere en 0, 3, 0, 12, 2 y 6 perturbaciones de 33 -- **no en
casi todas**. La razon esta en las dos columnas de suma: las dos emisiones
predicen **casi ningun gen significativo en ningun sitio**. Nuestro suelo suma 0,
1, 0, 3, 0 y 0 genes significativos sobre las 33 perturbaciones de cada linea; el
`b` oficial suma 0, 2, 0, 14, 2 y 7. En A549 y HAP1 las dos suman cero, asi que
ahi el control positivo **no discrimina**: las dos emisiones coinciden en el lado
predicho porque las dos estan vacias, no porque sean la misma. Donde si
discrimina -- HT29, con 12 perturbaciones distintas -- el lado real sigue saliendo
identico celda por celda. Esa es la evidencia que sostiene la invariancia; las
otras cinco lineas la acompanan sin probarla por si solas.

### 2.3 Mitad 2: pero NO es invariante a la POBLACION

El `r` del anchor no se calcula sobre la poblacion completa: cada particion puntua
`half_b` contra `half_a` (`anchor.py:120`), media poblacion a cada lado. Y el gate
es del lado real, luego depende de cuantas celulas tiene ese lado.

`anchor_splits.parquet` estampa `n_perturbations` por metrica y particion. Eso es
el tamano de la cohorte, directamente del artefacto oficial:

| linea | pds (5 particiones) | jac | lfc_nmae | fid | reach | run completo: lfc/fid/reach |
|---|---|---|---|---|---|---|
| A549 | 33/33/33/33/33 | 33/33/33/33/33 | 8/8/8/8/8 | 19/20/18/20/19 | 14/16/15/17/18 | 8 / 18 / 18 |
| BXPC3 | 33/33/33/33/33 | 33/33/33/33/33 | 23/23/23/23/23 | 32/30/32/30/31 | 30/29/31/30/30 | 23 / 33 / 33 |
| HAP1 | 33/33/33/33/33 | 33/33/33/33/33 | 8/8/8/8/8 | 20/19/24/20/26 | 16/17/20/16/20 | 8 / 19 / 19 |
| HT29 | 33/33/33/33/33 | 33/33/33/33/33 | 15/15/15/15/15 | 25/27/27/25/30 | 23/24/25/23/25 | 15 / 30 / 27 |
| K562 | 33/33/33/33/33 | 33/33/33/33/33 | 9/9/9/9/9 | 23/19/23/22/18 | 16/18/23/13/17 | 9 / 26 / 25 |
| MCF7 | 33/33/33/33/33 | 33/33/33/33/33 | 33/33/33/33/33 | 33/31/33/32/33 | 32/30/31/31/32 | 33 / 33 / 33 |

Dos cosas, las dos medidas:

- **`lfc_nmae` no se mueve.** Su cohorte es la misma en las 5 particiones Y en el
  run completo, en las seis lineas. Su gate es `>= min_gate_size = 10` genes
  significativos del lado real (`metrics/de.py:601`), y perder la mitad de las
  celulas no cruza ese umbral en ninguna perturbacion de este panel.
- **`fid` y `reach` si se mueven**, y entre particiones de la misma linea. Su
  cohorte es `n_conf > 0`: basta un gen que pierda la significancia al quedarse
  con media poblacion para que la perturbacion salga.

### 2.4 Las cohortes, por linea

La cohorte del lado real que la hipotesis propone -- `n_conf > 0`, las
perturbaciones donde la referencia tiene senal DE -- medida sobre la poblacion
completa, con `_reference_stats` del marcador:

| linea | panel | n_DE_real > 0 | n_conf > 0 (**la cohorte**) | n_conf >= 10 (gate de lfc) | q definido | n_conf min/mediana/max |
|---|---|---|---|---|---|---|
| A549 | 33 | 18 | 18 | 8 | 18 | 0 / 1 / 1021 |
| BXPC3 | 33 | 33 | 33 | 23 | 33 | 1 / 32 / 1164 |
| HAP1 | 33 | 19 | 19 | 8 | 19 | 0 / 2 / 207 |
| HT29 | 33 | 27 | 27 | 15 | 27 | 0 / 4 / 615 |
| K562 | 33 | 25 | 25 | 9 | 25 | 0 / 2 / 281 |
| MCF7 | 33 | 33 | 33 | 33 | 33 | 11 / 197 / 1438 |

`n_DE_real > 0` y `n_conf > 0` dan el mismo numero: los dos cuentan genes
significativos del lado real con el gen propio excluido, y son el mismo conteo
escrito en dos sitios del paquete (`contar_de_por_perturbacion`, que repite las dos
lineas de filtrado de `de_sig_jaccard`, y `_reference_stats`, direction.py:343-350).
Que coincidan no era obligatorio y se comprueba aqui.

### 2.5 Veredicto de Hueco 2

**La cohorte del lado real SI es invariante a la submission -- medido, celda por
celda, en las seis lineas -- y NO es invariante a la poblacion.** Las dos mitades
de la hipotesis no caen del mismo lado, y la consecuencia es concreta:

- comparar `u` contra `b` sobre la cohorte del lado real es legitimo: es la misma
  cohorte, por construccion y por medicion;
- comparar cualquiera de los dos contra `r` sobre "la" cohorte no lo es sin mas,
  porque `r` se calcula sobre media poblacion y su cohorte de `fid` y `reach` no es
  la del run completo. En `jac` y `pds` no estorba (33/33 siempre); en `lfc_nmae`
  tampoco (no se movio en ninguna linea ni particion); en `fid` y `reach` si.

Lo que la cohorte del lado real NO resuelve: la cohorte de `|union| != 0` que la
componente genuina de `jac` necesita. `|union| = |R u P|` incluye el lado predicho,
asi que **esa** cohorte no es invariante a la submission. Medido:

| linea | E de nuestro suelo | E del b oficial | en los dos | en uno solo | E de las 5 particiones de r | misma cohorte u/b | 33 - cohorte | celdas que la mueven |
|---|---|---|---|---|---|---|---|---|
| A549 | 15 | 15 | 15 | 0 | 14/13/15/13/14 | si | 15 | -- |
| BXPC3 | 0 | 0 | 0 | 0 | 1/3/1/3/2 | si | 0 | -- |
| HAP1 | 14 | 14 | 14 | 0 | 13/14/9/13/7 | si | 14 | -- |
| HT29 | 6 | 3 | 3 | 3 | 8/6/6/8/3 | **no** | 6 | BRD9, IRF7, RAP1GAP |
| K562 | 8 | 7 | 7 | 1 | 10/14/10/11/15 | **no** | 8 | ZBP1 |
| MCF7 | 0 | 0 | 0 | 0 | 0/2/0/1/0 | si | 0 | -- |

**Las dos cohortes estan ligadas por una identidad, y la medicion la fija.** La
columna `E de nuestro suelo` coincide exactamente con `33 - cohorte` en las seis
lineas. No es casualidad: `|union| = |R u P|`, nuestro suelo no predice **ningun**
gen significativo en casi ningun sitio (§2.2: suma 0 en cuatro lineas de seis),
asi que su union es vacia exactamente donde el lado real esta vacio. **La `E` de
nuestro suelo ES el complemento de la cohorte del lado real.**

Y la columna `celdas que la mueven` nombra lo que rompe esa identidad para el `b`:
las perturbaciones con `n_conf = 0` donde el `b` SI predijo algo significativo.
Son cuatro en total en las seis lineas: **BRD9, IRF7 y RAP1GAP en HT29, y ZBP1 en
K562**. En esas cuatro el `b` pierde el 1.0 gratis de la convencion y se queda con
un jaccard genuino de 0 (interseccion 0 sobre union 2), que es por lo que su `E` es
3 y 7 donde la de nuestro suelo es 6 y 8. Cada una cuesta `1/33 = 0.030303`.

Eso es la demostracion concreta, con nombre y apellido, de que la cohorte de union
no vacia **la mueve la submission**: cuatro perturbaciones de 198, y mueven el `b`
de dos lineas.

La cohorte de union no vacia **cambia con la submission** -- es la columna `en uno
solo` -- y cambia tambien de particion a particion del anchor. Con `u` y `b`
genuinos en 0 no estorbaba porque el promedio de un conjunto vacio de jaccards
genuinos es 0 en los dos casos. Con un candidato que SI tenga senal, la media
genuina de `u` y la de `b` se calcularian sobre conjuntos distintos.

**Eso no lo arregla la cohorte del lado real**: una perturbacion con
`n_conf = 0` puede tener union no vacia (si la submission predice algo
significativo) y una con `n_conf > 0` siempre la tiene. La cohorte del lado real
es un SUBCONJUNTO de la de union no vacia, y la diferencia es exactamente lo que la
submission controla. Lo que la medicion dice es cual de las dos elegir si se quiere
una cohorte no manipulable; no dice que la componente genuina de `jac` quede bien
definida sobre ella. Para eso haria falta decidir que recibe una perturbacion con
`n_conf > 0` cuyo `|union|` sea no vacio solo por el lado predicho, y esa decision
no es de este encargo.

---

## Hueco 3 — ¿nuestro `b` de `reach` cae sobre NUESTRO `c/N_conf`?

### 3.1 De donde sale cada pieza

- **`reach = k* / N_conf`**, `metrics/direction.py:1054`.
- **`N_conf`** lo construye `_reference_stats` desde `prepared.real_df` SOLO
  (`direction.py:293`, `343-350`): genes del lado real con `p_adj < alpha`,
  excluido el gen propio del target. No se persiste en ningun artefacto; aqui se
  obtiene llamando a esa funcion, sin reimplementarla.
- **`c`**: se usa **`c = 0.96`**, de `direction.py:967`: *"the fitted `c` is ~0.96
  for a coin-flip predictor against ~1.90 for one that ranks a single
  HABITUAL-DIRECTION gene first"*. Es un ajuste medido por el propio paquete en
  0.13.0 sobre 200 replicas por celda (`direction.py:965-966`), no una constante
  del codigo ejecutable: no hay ningun `c` en el modulo. El punto sin habilidad es
  `~c/N_conf` (`direction.py:946`), **no** una constante.

⚠️ El propio codigo acota donde la ley falla, y hay que decirlo antes de usarla:
`E[reach_raw]` medido *"runs 0.5100 at `N_conf = 1` to 0.0019 at 500"*
(`direction.py:965-966`). En `N_conf = 1` la ley daria 0.96 y lo medido fue 0.51:
la ley sobreestima donde `N_conf` es diminuto. La tabla trae la distribucion de
`N_conf` para que se vea cuanto peso tiene esa zona en nuestras lineas.

### 3.2 El azar esperado sobre NUESTRAS perturbaciones

Para cada linea: la media de `c/N_conf` sobre las perturbaciones que **entran** a
`reach` (las que tienen `reach` no-NaN en el run, que son exactamente las de
`N_conf > 0`, `direction.py:1048-1053`), contra nuestro `b` oficial de `reach`.

| linea | n entran a reach | N_conf min | N_conf mediana | N_conf max | N_conf <= 3 | azar = media(0.96/N_conf) | nuestro b (bundle) | b recomputado (media por pert) | b - azar | b / azar | c implicado por nuestro b | reach de nuestro suelo |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A549 | 18 | 1 | 5 | 1021 | 8 | 0.335352 | 0.314403 | 0.314403 | -0.020949 | 0.938x | 0.9000 | 0.144540 |
| BXPC3 | 33 | 1 | 32 | 1164 | 2 | 0.091684 | 0.344686 | 0.344686 | +0.253002 | 3.759x | 3.6091 | 0.142889 |
| HAP1 | 19 | 2 | 7 | 207 | 6 | 0.208835 | 0.354552 | 0.354552 | +0.145716 | 1.698x | 1.6298 | 0.192408 |
| HT29 | 27 | 1 | 10 | 615 | 10 | 0.224932 | 0.327624 | 0.327624 | +0.102692 | 1.457x | 1.3983 | 0.216033 |
| K562 | 25 | 1 | 4 | 281 | 12 | 0.402085 | 0.406814 | 0.406814 | +0.004730 | 1.012x | 0.9713 | 0.200448 |
| MCF7 | 33 | 11 | 197 | 1438 | 0 | 0.010304 | 0.101892 | 0.101892 | +0.091589 | 9.889x | 9.4935 | 0.013563 |

La columna `b recomputado` esta como control: es la media de los `reach` por
perturbacion del run del `b`, y debe coincidir con el `b` del bundle, que es
`agg = mean` sobre las que entran. La ultima columna es nuestro suelo sobre la
MISMA cohorte -- §2.2 midio que `n_conf` es identico entre los dos objetos, asi
que las dos medias se calculan sobre las mismas perturbaciones.

### 3.3 Veredicto de Hueco 3

**No coinciden, y la discrepancia no tiene un signo unico.** Linea por linea, el
cociente `b / azar`:

- A549: **0.94x** -- por DEBAJO del azar
- K562: **1.01x** -- practicamente ENCIMA del azar
- HT29: **1.46x** -- por ENCIMA del azar
- HAP1: **1.70x** -- por ENCIMA del azar
- BXPC3: **3.76x** -- por ENCIMA del azar
- MCF7: **9.89x** -- por ENCIMA del azar

Asi que **4 por encima, 1 encima y 1 por debajo.** No hay una direccion
unica que atribuir a nada.

Pero lo decisivo no es el signo: es la **amplitud**. El azar esperado varia
**39x** entre lineas (0.0103 en MCF7 a 0.4021 en K562), porque la distribucion de `N_conf` lo es. El `b` varia
solo **4.0x** (0.1019 a 0.4068).

La correlacion entre los dos es **+0.734**, asi que el `b` se mueve en la
MISMA direccion que la ley -- no es independiente de ella. Lo que no hace es
moverse por lo mismo: un objeto que estuviera en el punto sin habilidad tendria que
escalar con `1/N_conf`, y un cambio de 39x en el azar le produce 4x. **La forma de
la dependencia es otra**, y la seccion siguiente mide cual.

### 3.3bis De que esta hecho nuestro `b` de `reach`

La tabla por perturbacion dice de donde sale, y no hace falta inventar nada:

| linea | n entran | perts con reach = 1.0 exacto | sus N_conf | suma de reach | suma de los 1.0 | fraccion del total | b | media sin los 1.0 | perts con N_conf <= 3 | N_conf min |
|---|---|---|---|---|---|---|---|---|---|---|
| A549 | 18 | 5 | [1, 2, 2, 3, 10] | 5.6592 | 5.0000 | 0.884 | 0.314403 | 0.050711 | 8 | 1 |
| BXPC3 | 33 | 6 | [1, 3, 4, 5, 6, 10] | 11.3746 | 6.0000 | 0.527 | 0.344686 | 0.199061 | 2 | 1 |
| HAP1 | 19 | 4 | [2, 2, 2, 3] | 6.7365 | 4.0000 | 0.594 | 0.354552 | 0.182432 | 6 | 2 |
| HT29 | 27 | 5 | [1, 1, 2, 3, 3] | 8.8459 | 5.0000 | 0.565 | 0.327624 | 0.174811 | 10 | 1 |
| K562 | 25 | 7 | [1, 1, 1, 1, 2, 2, 2] | 10.1704 | 7.0000 | 0.688 | 0.406814 | 0.176131 | 12 | 1 |
| MCF7 | 33 | 0 | -- | 3.3624 | 0.0000 | 0.000 | 0.101892 | 0.101892 | 0 | 11 |

**Entre el 53% y el 88% del `b` de `reach` lo ponen entre 4 y 7 perturbaciones que
sacan 1.0 exacto, y todas tienen `N_conf <= 10`** -- la mayoria `<= 3`. Sin ellas,
la media del resto cae a 0.05-0.20.

El codigo describe exactamente esa banda y dice que en la val oficial A **esta
vacia**: *"val A's smallest N_conf is 4, so the `N_conf <= 3` band where the head
flip is worth 0.29-0.51 is empty"* (`direction.py:977-979`). En nuestras lineas no
lo esta: hay 8, 2, 6, 10, 12 y 0 perturbaciones en esa banda. **MCF7 es la unica
con la banda vacia** (`N_conf` minimo 11) y es tambien la unica sin ninguna
perturbacion en 1.0, y la de `b` mas bajo (0.1019).

Los dos regimenes, dichos aparte porque son distintos:

- **Lineas con la banda poblada** (las otras cinco): el `b` lo decide el CONTEO de
  perturbaciones que resuelven del todo, no el presupuesto. Por eso es plano.
- **MCF7, banda vacia**: ninguna resuelve del todo, y aun asi el `b` (0.1019) queda
  casi 10x por encima del azar (0.0103). La caida `1/N_conf` de la ley es mucho mas
  empinada que lo que el `b` consigue de verdad sobre perturbaciones con `N_conf`
  de tres digitos.

En **ninguna** linea el `b` se apoya sobre `c/N_conf`, y en las dos puntas la razon
medida es distinta. No se propone ninguna causa comun: no esta en el codigo.

### 3.3ter Lo que eso dice

- El extremo 0 de la escala de `reach` **no** es el piso de azar sobre nuestra
  poblacion. En 4 lineas el `b` esta por
  encima, asi que un predictor de moneda cobraria `s < 0`; en A549 esta un 6% por
  debajo, asi que ahi cobraria `s > 0`; en K562 cae encima por 1%. Lo que no hay
  es coincidencia sistematica.
- La cifra con que el codigo ilustra la coincidencia (*"the baseline arm reads
  0.042314 against an ESTIMATED no-skill 0.0378, so the baseline already sits at
  the no-skill point"*, `direction.py:971-972`) es de la val oficial A, cuya
  distribucion de `N_conf` es otra: *"val A's smallest N_conf is 4"*
  (`direction.py:978`). En nuestras lineas la mediana de `N_conf` es de un digito
  en varias y el minimo es 1. La ley `c/N_conf` es justamente la que el codigo
  declara sobreestimada en esa zona, asi que el azar que calculamos con ella es si
  acaso un techo del azar, y nuestro `b` queda aun mas arriba.
- La columna `c implicado` resuelve `media(c/N_conf) = b` por `c`. Es la lectura
  no-inventiva de la discrepancia: nuestro `b` se comporta como un predictor con
  `c` entre 0.90 y 9.49, contra el 0.96 de la moneda y el 1.90
  del que *"ranks a single HABITUAL-DIRECTION gene first"* (`direction.py:967`).
  El codigo dice de ese segundo predictor que es *"close to free"*. **No se afirma
  que el `b` sea ese predictor**: se dice donde cae su `c` implicado y de donde
  salen los dos puntos de referencia.

### 3.4 Las filas trabajadas de Hueco 3

Los intermedios crudos: `N_conf` por perturbacion, el azar que implica, y el
`reach` que el `b` saco de verdad.

**Criterio de seleccion, declarado antes de mirar los valores**: la PRIMERA y la
ULTIMA linea en orden alfabetico -- A549 y MCF7 -- y de cada una **todas** las
perturbaciones que entran a `reach`, sin truncar. Es un criterio posicional y da
los dos conjuntos completos, no una muestra. Las filas van ordenadas por `N_conf`
creciente, que es una funcion del lado real y por tanto invariante a la submission
(§2.2).

| linea | perturbacion | N_conf | azar 0.96/N_conf | reach del b | reach - azar | reach = 1.0 | q | n_DE_real | n_DE_pred | \|union\| |
|---|---|---|---|---|---|---|---|---|---|---|
| A549 | BRD9 | 1 | 0.960000 | 1.000000 | +0.040000 | **si** | 1.0000 | 1 | 0 | 1 |
| A549 | MAPK8 | 1 | 0.960000 | 0.000000 | -0.960000 |  | 1.0000 | 1 | 0 | 1 |
| A549 | SP100 | 1 | 0.960000 | 0.000000 | -0.960000 |  | 1.0000 | 1 | 0 | 1 |
| A549 | STAT1 | 1 | 0.960000 | 0.000000 | -0.960000 |  | 1.0000 | 1 | 0 | 1 |
| A549 | RAP1GAP | 2 | 0.480000 | 1.000000 | +0.520000 | **si** | 1.0000 | 2 | 0 | 2 |
| A549 | ZBP1 | 2 | 0.480000 | 1.000000 | +0.520000 | **si** | 1.0000 | 2 | 0 | 2 |
| A549 | DTX3L | 3 | 0.320000 | 1.000000 | +0.680000 | **si** | 0.6667 | 3 | 0 | 3 |
| A549 | TRAFD1 | 3 | 0.320000 | 0.000000 | -0.320000 |  | 0.6667 | 3 | 0 | 3 |
| A549 | CEBPG | 4 | 0.240000 | 0.250000 | +0.010000 |  | 0.7500 | 4 | 0 | 4 |
| A549 | JUN | 6 | 0.160000 | 0.166667 | +0.006667 |  | 0.6667 | 6 | 0 | 6 |
| A549 | IRF1 | 10 | 0.096000 | 1.000000 | +0.904000 | **si** | 0.5000 | 10 | 0 | 10 |
| A549 | ETS2 | 19 | 0.050526 | 0.052632 | +0.002105 |  | 0.5789 | 19 | 0 | 19 |
| A549 | MAP3K14 | 40 | 0.024000 | 0.025000 | +0.001000 |  | 0.6250 | 40 | 0 | 40 |
| A549 | STAT3 | 84 | 0.011429 | 0.047619 | +0.036190 |  | 0.6667 | 84 | 0 | 84 |
| A549 | JAK1 | 94 | 0.010213 | 0.000000 | -0.010213 |  | 0.5638 | 94 | 0 | 94 |
| A549 | IFNAR1 | 494 | 0.001943 | 0.010121 | +0.008178 |  | 0.6781 | 494 | 0 | 494 |
| A549 | IFNAR2 | 751 | 0.001278 | 0.057257 | +0.055979 |  | 0.5366 | 751 | 0 | 751 |
| A549 | TYK2 | 1021 | 0.000940 | 0.049951 | +0.049011 |  | 0.6592 | 1021 | 0 | 1021 |
| MCF7 | PARP12 | 11 | 0.087273 | 0.272727 | +0.185455 |  | 0.8182 | 11 | 0 | 11 |
| MCF7 | STAT4 | 20 | 0.048000 | 0.400000 | +0.352000 |  | 0.9500 | 20 | 0 | 20 |
| MCF7 | RAP1GAP | 31 | 0.030968 | 0.064516 | +0.033548 |  | 0.9677 | 31 | 0 | 31 |
| MCF7 | STAT3 | 44 | 0.021818 | 0.045455 | +0.023636 |  | 1.0000 | 44 | 0 | 44 |
| MCF7 | BRD9 | 56 | 0.017143 | 0.071429 | +0.054286 |  | 0.9821 | 56 | 1 | 57 |
| MCF7 | HERC6 | 67 | 0.014328 | 0.029851 | +0.015522 |  | 0.8657 | 67 | 0 | 67 |
| MCF7 | STAT1 | 67 | 0.014328 | 0.223881 | +0.209552 |  | 0.9403 | 67 | 0 | 67 |
| MCF7 | ETS2 | 72 | 0.013333 | 0.458333 | +0.445000 |  | 0.8056 | 72 | 0 | 72 |
| MCF7 | SP100 | 93 | 0.010323 | 0.000000 | -0.010323 |  | 0.9677 | 93 | 0 | 93 |
| MCF7 | DTX3L | 102 | 0.009412 | 0.156863 | +0.147451 |  | 0.9510 | 102 | 0 | 102 |
| MCF7 | JAK1 | 115 | 0.008348 | 0.008696 | +0.000348 |  | 1.0000 | 115 | 0 | 115 |
| MCF7 | BATF2 | 131 | 0.007328 | 0.122137 | +0.114809 |  | 0.8702 | 131 | 0 | 131 |
| MCF7 | IRF1 | 163 | 0.005890 | 0.061350 | +0.055460 |  | 0.8528 | 163 | 0 | 163 |
| MCF7 | ELK1 | 164 | 0.005854 | 0.073171 | +0.067317 |  | 0.8049 | 164 | 0 | 164 |
| MCF7 | MAP3K14 | 175 | 0.005486 | 0.000000 | -0.005486 |  | 0.9829 | 175 | 0 | 175 |
| MCF7 | STAT5A | 187 | 0.005134 | 0.010695 | +0.005561 |  | 0.9893 | 187 | 0 | 187 |
| MCF7 | AKT1 | 197 | 0.004873 | 0.406091 | +0.401218 |  | 0.9797 | 197 | 0 | 197 |
| MCF7 | RNF114 | 233 | 0.004120 | 0.536481 | +0.532361 |  | 0.9914 | 233 | 0 | 233 |
| MCF7 | ETS1 | 236 | 0.004068 | 0.021186 | +0.017119 |  | 0.9873 | 236 | 0 | 236 |
| MCF7 | CEBPG | 294 | 0.003265 | 0.003401 | +0.000136 |  | 0.9626 | 294 | 0 | 294 |
| MCF7 | VAV1 | 415 | 0.002313 | 0.019277 | +0.016964 |  | 0.9639 | 415 | 0 | 415 |
| MCF7 | HERC5 | 458 | 0.002096 | 0.006550 | +0.004454 |  | 0.9651 | 458 | 0 | 458 |
| MCF7 | SOCS1 | 460 | 0.002087 | 0.065217 | +0.063130 |  | 0.9935 | 460 | 1 | 461 |
| MCF7 | IRF7 | 474 | 0.002025 | 0.002110 | +0.000084 |  | 0.9662 | 474 | 0 | 474 |
| MCF7 | H1F0 | 529 | 0.001815 | 0.137996 | +0.136181 |  | 1.0000 | 529 | 1 | 530 |
| MCF7 | ZBP1 | 578 | 0.001661 | 0.000000 | -0.001661 |  | 0.9792 | 578 | 0 | 578 |
| MCF7 | JUN | 627 | 0.001531 | 0.003190 | +0.001659 |  | 0.9697 | 627 | 0 | 627 |
| MCF7 | IFNAR1 | 842 | 0.001140 | 0.036817 | +0.035677 |  | 0.6805 | 842 | 0 | 842 |
| MCF7 | TYK2 | 1027 | 0.000935 | 0.012658 | +0.011723 |  | 0.7595 | 1027 | 0 | 1027 |
| MCF7 | MAPK8 | 1053 | 0.000912 | 0.088319 | +0.087407 |  | 0.9696 | 1053 | 1 | 1053 |
| MCF7 | TRAFD1 | 1140 | 0.000842 | 0.004386 | +0.003544 |  | 0.9763 | 1140 | 0 | 1140 |
| MCF7 | IFNAR2 | 1368 | 0.000702 | 0.018275 | +0.017573 |  | 0.9064 | 1368 | 1 | 1369 |
| MCF7 | ID2 | 1438 | 0.000668 | 0.001391 | +0.000723 |  | 0.9722 | 1438 | 2 | 1440 |

Son **51** filas: 18 de A549 y 33 de MCF7, los dos conjuntos enteros.

Tres cosas se leen directamente en ellas:

1. La columna `azar` es `0.96/N_conf` celda a celda, y su media sobre cada linea
   es el `azar` de la tabla de §3.2.
2. `n_DE_pred = 0` en todas las filas de A549 -- el `b` no predice ni un gen
   significativo -- y aun asi cinco de ellas sacan `reach = 1.0`. Es la propiedad
   que `descomposicion.md` ya habia explicado: con `universe='adjudicated'` el
   ranking recorre solo los genes significativos de la REFERENCIA
   (`direction.py:980-982`), asi que el conjunto predicho vacio no impide nada.
3. `n_DE_real` y `|union|` coinciden en todas las filas de A549, por lo mismo: la
   union es el lado real solo. Y `n_DE_real` coincide con `N_conf`, que es la
   identidad que §2.4 comprueba.

---

## Hueco 4 — `pds`: la causa en HAP1 y BXPC3

### 4.1 La hipotesis, declarada antes de medir

> El estadistico asigna exactamente 0.5 a una fila empatada
> (`discrimination.py:70`, `tie_policy="midrank"`). Nuestro suelo esta MENOS
> disperso, asi que produciria mas casi-empates y se pegaria mas a 0.5; el `b`, mas
> disperso, se alejaria de 0.5 en una direccion u otra segun la linea.

### 4.2 La distribucion por perturbacion

`pds_cosine` emite una fila por perturbacion (`agg="mean"`, `catalog.py:319-324`).
Las 33 de cada linea, para los dos objetos:

| linea | n perts | media suelo | sd suelo | media b | sd b | b - suelo | corr(suelo, b) por pert |
|---|---|---|---|---|---|---|---|
| A549 | 33 | 0.474432 | 0.312543 | 0.486742 | 0.302396 | +0.012311 | -0.4081 |
| BXPC3 | 33 | 0.517992 | 0.301711 | 0.476326 | 0.296067 | -0.041667 | -0.1079 |
| HAP1 | 33 | 0.512311 | 0.258261 | 0.437500 | 0.303370 | -0.074811 | -0.2432 |
| HT29 | 33 | 0.490530 | 0.317431 | 0.501894 | 0.298050 | +0.011364 | -0.1331 |
| K562 | 33 | 0.495265 | 0.315715 | 0.520833 | 0.301184 | +0.025568 | +0.1115 |
| MCF7 | 33 | 0.512311 | 0.291728 | 0.503788 | 0.305292 | -0.008523 | -0.3448 |

Y la fraccion que cae **exactamente** en 0.5, que es lo que la hipotesis predecia:

| linea | n suelo | exactamente 0.5 (suelo) | fraccion | n b | exactamente 0.5 (b) | fraccion |
|---|---|---|---|---|---|---|
| A549 | 33 | 0 | 0.0000 | 33 | 2 | 0.0606 |
| BXPC3 | 33 | 0 | 0.0000 | 33 | 2 | 0.0606 |
| HAP1 | 33 | 0 | 0.0000 | 33 | 3 | 0.0909 |
| HT29 | 33 | 1 | 0.0303 | 33 | 0 | 0.0000 |
| K562 | 33 | 1 | 0.0303 | 33 | 4 | 0.1212 |
| MCF7 | 33 | 0 | 0.0000 | 33 | 0 | 0.0000 |

**2 de 198 en nuestro suelo, 11 de 198 en el `b`.** El `b` tiene MAS celdas en 0.5 exacto que
nuestro suelo. La hipotesis predecia lo contrario.

### 4.3 Dos razones por las que el conteo de 0.5 no era la prueba que parecia

**Primera: 0.5 es un punto de la rejilla.** `vcc2026` usa
`rank_denominator: n-1` (`configs/vcc2026.yaml:33`), asi que con n = 33
perturbaciones el valor es `1 - rango/32` y vive en la rejilla
`{0, 1/32, ..., 32/32}` -- 33 puntos, y `16/32 = 0.5` es uno de ellos. Verificado:
todos los valores de los dos objetos en las seis lineas son multiplos exactos de
1/32. Una fila SIN ningun empate puede caer en 0.5 por rango 16. Contar 0.5
exactos no cuenta empates.

**Segunda: el unico mecanismo de empate conocido no se activa.** El codigo nombra
el caso (`discrimination.py:63-68`): *"under `cosine` a zero-norm predicted effect
(a submission pasting the reference control cells for a target) ties the entire row
at distance 1.0"*, y `midrank` *"scores every such target exactly 0.5"*
(`discrimination.py:69`). `descomp_pds.json` ya midio `frac_norma_cero` en los tres
objetos y las seis lineas: **0 en todos**. Ningun delta tiene norma cero, asi que
ninguna fila esta empatada entera. Sobre flotantes continuos, un empate exacto
entre dos distancias cosine distintas tiene probabilidad practicamente nula.

### 4.4 Entonces de que esta hecha la variacion: el yardstick correcto

Lo que la tabla de §4.2 muestra es que **ninguno de los dos objetos esta
concentrado en 0.5**. Los dos estan repartidos casi uniformemente sobre la rejilla:

- `sd` teorico de la rejilla uniforme de 33 puntos: **0.297560**
- `sd` medido: suelo 0.3125 / 0.3017 / 0.2583 / 0.3174 / 0.3157 / 0.2917; `b` 0.3024 / 0.2961 / 0.3034 / 0.2981 / 0.3012 / 0.3053

Eso es exactamente lo que el codigo dice que debe pasar con una emision sin
informacion: *"the ranks are a bijection onto `{0..n-1}`, so the panel mean is
0.5"* (`discrimination.py:73-74`), y `catalog.py:313-314`: *"Its comparator sits at
the metric's chance point (a uniform rank is 0.5)"*. La media de panel sale ~0.5
**no porque las filas se peguen a 0.5, sino porque los rangos se reparten**.

Y una media de 33 rangos casi uniformes tiene ruido de muestreo propio:

```
SE(media de panel) = sd_rejilla / sqrt(33) = 0.297560 / 5.744563 = 0.051798
```

**Ese** es el yardstick de `|b - u|`, no el `sd` entre particiones del anchor. Con
el `SE` de la diferencia calculado sobre las diferencias por perturbacion
observadas (que no asume independencia):

| linea | b - u | SE media u | SE media b | SE dif (indep) | SE dif (observado) | \|b-u\| / SE dif | sd entre particiones de r | \|b-u\| / sd del anchor |
|---|---|---|---|---|---|---|---|---|
| A549 | +0.012311 | 0.0544 | 0.0526 | 0.0757 | 0.0898 | **0.14** | 0.0290 | 0.42 |
| BXPC3 | -0.041667 | 0.0525 | 0.0515 | 0.0736 | 0.0775 | **0.54** | 0.0168 | 2.48 |
| HAP1 | -0.074811 | 0.0450 | 0.0528 | 0.0694 | 0.0772 | **0.97** | 0.0153 | 4.91 |
| HT29 | +0.011364 | 0.0553 | 0.0519 | 0.0758 | 0.0807 | **0.14** | 0.0171 | 0.67 |
| K562 | +0.025568 | 0.0550 | 0.0524 | 0.0760 | 0.0716 | **0.36** | 0.0176 | 1.46 |
| MCF7 | -0.008523 | 0.0508 | 0.0531 | 0.0735 | 0.0852 | **0.10** | 0.0163 | 0.52 |

**Ninguna linea pasa de 0.97 SE.** Las dos que "excedian" el `sd` del
anchor caen en HAP1 0.97, BXPC3 0.54 contra el yardstick correcto. Las dos ultimas columnas muestran la diferencia
entre los dos yardsticks: el del anchor declara 2.5 y 4.9 sigma donde el del estadistico declara 0.5 y 1.0.

**Por que el `sd` entre particiones del anchor era el yardstick equivocado**, con
los numeros a la vista: las medias de panel de `pds` de las 5 particiones son

| linea | pds de las 5 particiones de r | media | sd | media u | media b |
|---|---|---|---|---|---|
| A549 | 0.8097 / 0.7803 / 0.7964 / 0.7784 / 0.7244 | 0.7778 | 0.0290 | 0.4744 | 0.4867 |
| BXPC3 | 0.9527 / 0.9489 / 0.9631 / 0.9811 / 0.9299 | 0.9551 | 0.0168 | 0.5180 | 0.4763 |
| HAP1 | 0.7652 / 0.7926 / 0.7803 / 0.8049 / 0.8049 | 0.7896 | 0.0153 | 0.5123 | 0.4375 |
| HT29 | 0.8163 / 0.7822 / 0.7756 / 0.7670 / 0.7765 | 0.7835 | 0.0171 | 0.4905 | 0.5019 |
| K562 | 0.8248 / 0.8542 / 0.8040 / 0.8220 / 0.8087 | 0.8227 | 0.0176 | 0.4953 | 0.5208 |
| MCF7 | 0.9527 / 0.9489 / 0.9384 / 0.9347 / 0.9062 | 0.9362 | 0.0163 | 0.5123 | 0.5038 |

Las particiones del anchor son real-contra-real: sus valores por perturbacion se
apilan cerca de 1.0, donde la rejilla deja poco sitio para moverse, y su `sd` de
panel es chico por eso. Es el `sd` de otro regimen. Usarlo para juzgar la
variacion de un objeto sin informacion, cuyos rangos se reparten por toda la
rejilla, subestima el ruido por un factor de 2 a 3.

### 4.5 Veredicto de Hueco 4

**La hipotesis no se sostiene, en sus dos mitades.**

1. El mecanismo de empate que proponia no se activa: ningun delta tiene norma cero
   (`descomp_pds.json`, seis lineas, tres objetos), que es el unico caso de fila
   entera empatada que el codigo nombra.
2. La prediccion direccional sale al reves: el `b` tiene mas celdas en 0.5 exacto
   (11) que nuestro suelo (2) -- y de todos modos ese conteo no mide
   empates, porque 0.5 esta en la rejilla.
3. Nuestro suelo no esta "menos disperso" en el espacio del estadistico: su `sd`
   por perturbacion es indistinguible del del `b`, y los dos coinciden con el de
   la rejilla uniforme. La menor dispersion que `descomposicion.md` midio estaba en
   el espacio de los DELTAS (distancia cosine entre pares), no en el de los RANGOS,
   y el estadistico puntua rangos.

**Lo que queda explicado**: el signo cambiante y los dos "excesos". Los dos objetos
estan en el punto sin informacion del estadistico; sus medias de panel son medias
de 33 rangos casi uniformes, con un `SE` propio de ~0.052;
las seis diferencias caben en ese ruido. Nada que atribuir al escalado de
`baseline.py:344`.

**Lo que queda sin explicar**: nada de esta pregunta. Lo que no se midio aqui es
si un candidato CON habilidad en `pds` tiene ese mismo `SE` -- no lo tendria: el
regimen del anchor muestra que cuando los valores se apilan arriba el `sd` baja.
Eso no era la pregunta y no se estima.

### 4.6 Las 20 filas trabajadas de Hueco 4

Criterio declarado antes de mirar los valores: las perturbaciones de **HAP1** -- la
linea del `|b - u|` mas grande, que es la que el encargo nombra -- en orden
alfabetico, las 20 primeras.

| linea | perturbacion | pds suelo | rango suelo | pds b | rango b | b - suelo | suelo en 0.5 | b en 0.5 |
|---|---|---|---|---|---|---|---|---|
| HAP1 | AKT1 | 0.750000 | 24/32 | 0.031250 | 1/32 | -0.718750 |  |  |
| HAP1 | BATF2 | 0.625000 | 20/32 | 0.125000 | 4/32 | -0.500000 |  |  |
| HAP1 | BRD9 | 0.406250 | 13/32 | 0.062500 | 2/32 | -0.343750 |  |  |
| HAP1 | CEBPG | 0.812500 | 26/32 | 0.125000 | 4/32 | -0.687500 |  |  |
| HAP1 | DTX3L | 0.843750 | 27/32 | 0.500000 | 16/32 | -0.343750 |  | si |
| HAP1 | ELK1 | 1.000000 | 32/32 | 0.718750 | 23/32 | -0.281250 |  |  |
| HAP1 | ETS1 | 0.468750 | 15/32 | 0.750000 | 24/32 | +0.281250 |  |  |
| HAP1 | ETS2 | 0.468750 | 15/32 | 0.500000 | 16/32 | +0.031250 |  | si |
| HAP1 | H1F0 | 0.625000 | 20/32 | 0.000000 | 0/32 | -0.625000 |  |  |
| HAP1 | HERC5 | 0.218750 | 7/32 | 0.406250 | 13/32 | +0.187500 |  |  |
| HAP1 | HERC6 | 0.812500 | 26/32 | 0.218750 | 7/32 | -0.593750 |  |  |
| HAP1 | ID2 | 0.656250 | 21/32 | 0.343750 | 11/32 | -0.312500 |  |  |
| HAP1 | IFNAR1 | 0.312500 | 10/32 | 1.000000 | 32/32 | +0.687500 |  |  |
| HAP1 | IFNAR2 | 0.000000 | 0/32 | 0.218750 | 7/32 | +0.218750 |  |  |
| HAP1 | IRF1 | 0.406250 | 13/32 | 0.843750 | 27/32 | +0.437500 |  |  |
| HAP1 | IRF7 | 0.125000 | 4/32 | 0.406250 | 13/32 | +0.281250 |  |  |
| HAP1 | JAK1 | 0.593750 | 19/32 | 0.937500 | 30/32 | +0.343750 |  |  |
| HAP1 | JUN | 0.187500 | 6/32 | 0.562500 | 18/32 | +0.375000 |  |  |
| HAP1 | MAP3K14 | 0.156250 | 5/32 | 0.531250 | 17/32 | +0.375000 |  |  |
| HAP1 | MAPK8 | 0.968750 | 31/32 | 0.343750 | 11/32 | -0.625000 |  |  |

Son 20 filas. Las columnas `rango` muestran que los valores son enteros
sobre 32: la rejilla. La media de la columna `pds b` sobre las 33 (no solo estas
20) es 0.437500, y la de `pds suelo` 0.512311; su
diferencia -0.074811 es la celda de HAP1 de la tabla de §4.4.

---

## Lo que no se hizo

- **No se redacto ningun criterio de adopcion.** No es de este encargo.
- **No se construyo ningun candidato.** Las dos emisiones que se corrieron son
  nuestro suelo (`predict_floor` sin modificar) y la del `b` oficial
  (`cell-eval2 baseline --seed 0`), las dos ya conocidas.
- **No se modifico ninguna metrica ni ningun arnes.** `build_fold_adatas`,
  `predict_floor`, `contar_de_por_perturbacion`, `_disjoint_halves` y
  `_reference_stats` se usaron tal cual. `descomp_rgen.py` -- script de analisis de
  esta sesion, no un arnes -- gano un argumento de linea de comandos y una columna
  de conteos por perturbacion; su metodo no cambio, y §1.1 lo verifica contra los
  valores que ya habia estampado el bundle.
- **No se selecciono ninguna celda por su valor ni por su signo.** Los criterios de
  seleccion que aparecen (las filas trabajadas de §3.4 y §4.6, y las lineas de cada
  tabla) estan declarados antes de los numeros y son posicionales o funcion del
  lado real. §3.4 muestra los dos conjuntos ENTEROS, sin truncar.
- **No se derivo ningun objetivo numerico por linea.**
- No se leyo nada de `~/vcc/benchmark/` ni de
  `/datos/benchmark_VEDADO_techo_metricas/`.
- Sin envios, sin `vcc sample`, sin commits.

