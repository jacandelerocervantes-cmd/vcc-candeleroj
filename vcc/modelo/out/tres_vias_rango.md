# Tres vias lado a lado como RANGO: suelo, candidato 01, estructura externa (gwps)

## 🔴 0. Que es este reporte y que no es

**Presenta datos.** Declara la construccion antes de medir, da las tablas con sus
controles y reporta los defectos de la propia construccion. **No lleva veredicto,
ni conclusiones, ni nada dirigido al siguiente candidato.**

🔴 **EL PUNTO DEL FORMATO.** Los pesos de A/B/C no transfieren a D/E/F; lo que
transfiere es un **rango de aproximacion**. Por eso aqui nada se juzga con una
correlacion calculada a mano: cada via se **emite por linea**, se **puntua contra
el bundle DE ESA LINEA**, y lo que sale son los **seis `s` normalizados**. Es la
unica cantidad comparable entre contextos, porque el extremo 0 y el extremo 1 son
de cada linea y difieren enormemente. **El entregable es el RANGO por metrica
(§5), no el promedio.**

**🔴 DOS POBLACIONES DISTINTAS en la via nueva.** gwps = K562 **sin estimular**,
knockdown **CRISPRi**. Jiang = seis lineas **con IFNB**. No son la misma cosa y
nada de esta construccion supone que lo sean.

---

## 1. La prediccion declarada ANTES de correr

Escrita antes de la primera corrida, para que la corrida la confronte:

> **H3 del traspaso**: `jac` en celdas de span pequeno dominara el agregado de
> cualquier candidato que prediga de verdad. A549 tiene span `0.0187 < 1/33`.
> Si sale asi, **NO se arregla prediciendo menos**: eso seria explotar la
> convencion (D-023). Se declara y se leen las otras celdas.

**Lo que salio, en §8.** Se confronta ahi y las otras celdas se leen en §5 y §6.

---

## 2. Construccion de la via de estructura externa

### 2.1 De donde sale la emision

| | |
|---|---|
| archivo | `/datos/ref/replogle/k562_gwps_raw_bulk.h5ad`, matriz `X` (11258 perturbaciones x 8248 genes) |
| perfil | `p_g = log2( (X_g + eps) / (media de X sobre las 514 filas core_control + eps) )`, `eps = 1e-8` -- **igual que `gwps_jiang_datos.md` §1.1** |
| genes | 8248 columnas, 2 nombres duplicados (`HSPA14`, `TBCE`, 4 filas) excluidos enteros; de los 8246 restantes **8124 caen** en los 34025 `var_names` de Jiang |
| codigo | `modelo/src/gwps_emision.py`, `modelo/src/gwps_via.py` |

### 2.2 🔴 LA MAGNITUD NO VIAJA: es una ELECCION DE ESCALA, no una cifra medida

Sigue declarada como hueco desde `candidato_01.md`. De `p_g` se usan **solo el
ORDEN y el SIGNO**:

```
dentro de los genes que BAJAN en gwps (m de ellos), ordenados por |p| descendente,
  el i-esimo recibe   s = -(m - i + 1)/m
dentro de los que SUBEN (n), igual con signo +:   s = +(n - i + 1)/n
r_h = 1 + s_h
```

El gen que mas baja recibe `s = -1` y se apaga del todo; el que menos, `s = -1/m`.
**`r` vive en [0, 2]** por construccion, es **adimensional**, y es funcion
**unicamente del rango**: el valor de `log2FC` de gwps no entra en ninguna cifra,
solo su orden y su signo. Es la misma clase de eleccion que la v3 de
`candidato_01` (`r = 1 - rho`, acotada en [0,2]): acotada, adimensional y sin
ninguna constante elegida aqui. **No es una medicion de cuanto baja un gen.**

Medido en la emision: 8124 genes participan por blanco, con **mediana `|s|` =
0.500**; en `STAT1` bajan 3055 y suben 5069.

### 2.3 Sin puerta de significancia, y la variante con puerta que NO se puntuo

La via se especifica como "el perfil da el orden y el signo", asi que entran los
8124 genes mapeados. Meter una puerta anadiria una constante de seleccion que no
esta en la especificacion.

🔴 **Se midio la puerta analoga a la de `candidato_01`** (BH al 0.05 del preset
`vcc2026` sobre `z = p_g / sd_nulo`, con `sd_nulo` la sd de las 514 filas
`core_control`) y el resultado esta en §9.3: de **0 a 829 genes** por blanco,
**mediana 1**, **7 de 27 blancos en cero**. **Esa variante NO se puntuo**, y queda
declarada como **no medida**: elegir entre las dos mirando el `s` esta prohibido.

### 2.4 🔴 El repliegue al SUELO de los 6 blancos ausentes

Los 6 blancos del panel que **no estan en gwps** -- `DTX3L`, `H1F0`, `HERC5`,
`HERC6`, `MAP3K14`, `ZBP1` -- reciben `r_h = 1` en todos los genes, que es
**exactamente la emision del suelo** (control remuestreado sin tocar). **Es un
repliegue declarado, no una eleccion**: no hay con que emitirlos por esta via.
Cuanto del resultado de cada linea viene de esos 6 esta en **§7**.

### 2.5 🔴 ID2 y RAP1GAP: P1, sin promediar

Los dos tienen filas P1 y P2 en gwps. **No se promedian** (`gwps_fiabilidad.md`
midio que P1 domina en 67.6 % de 546 pares, asi que la media no es neutral).
**Se usa P1** en la corrida principal. La variante con **P2 en esas dos filas**
se corrio aparte con una semilla y esta en **§10**.

### 2.6 El contrato de emision y el arnes

| | |
|---|---|
| arnes | `suelo_arnes.build_fold_adatas` **SIN MODIFICAR**; `predict_fn(control_pool_X, labels, rng)` |
| conteos | **ENTEROS muestreados**, nunca un perfil suave (`paso1_contrato.md`, HARD_CAP_NNZ). Se reusa `_reescala_entero` de `cand01_covar` sin modificar: Binomial para `r<=1`, Poisson aditivo para `r>1`; esperanza `x*r` y enteros en las dos ramas |
| comprobado en las 48 | `todo_entero` **48/48**, `hay_negativos` **0/48**, `min_valor` 1.0 |
| no transductivo | `predict_fn` recibe solo el pool de control. La tabla de gwps entra por clausura y **es la estructura externa que la via trae a proposito** |
| 🔴 propiedad de la via | los ratios **no dependen del pool de control**: son los mismos en las seis lineas. Esta emision **no estima nada en contexto** |

### 2.7 Semillas y verificaciones

| | |
|---|---|
| semillas de emision | familia **nueva** `20271000 + 100*fold + j`, `j = 0..7`. Sin colision con 2026000x (verdad), 2026100x (suelo), 2026200x/2026300x (techo), 2026400x (dispersion de `u`), 2026500x (dispersion de `b`), 2026600x (candidato 01) ni la 0 del baseline/anchor |
| poblacion real FIJA | `seed_real = 20260001 + fold_idx`, probada con **sha256** como en `dispersion_u.md`: **1 sola huella real por linea** en las 8 semillas, `poblacion_real_identica` **48/48** |
| emisiones distintas | **48 huellas `pred` distintas** de 48 corridas |
| `config_hash` de competencia | **verificado 48/48** (`rule_config_hash` == `competition_payload()['config_hash']`) |
| puntuacion | `score --real-bundle /datos/modelo/out/nativas_comp/BUNDLE_<linea>`, **48/48 aceptadas**, `real_bundle_id` correcto en las 6 |
| coste | 60-83 s por corrida; 48 corridas + 6 de la variante P2 |

---

## 3. Las tres vias

| via | que emite | de donde salen sus cifras |
|---|---|---|
| **SUELO** | control remuestreado con reemplazo | `dispersion_u.json`, **ya medido**, 48 corridas |
| **CANDIDATO 01** | respuesta lineal a lo largo de la covarianza del propio control, `r = 1 - rho` | `candidato_01*.json`, **ya medido**, 48 corridas |
| **GWPS** | orden y signo del perfil de k562_gwps, `r = 1 + s` por rango | **esta corrida**, 48 corridas |

Las tres con 8 semillas de emision, la misma poblacion real por linea, y
puntuadas contra el mismo bundle de competencia de cada linea.

---

## 4. Los seis `s` por linea y por via

Media sobre las 8 semillas, con su `sd`. `s` es `from_replicate` de
`score --real-bundle`.

| linea | metrica | **SUELO** | sd | **CAND01** | sd | **GWPS** | sd |
|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | **+0.0573** | 0.2840 | **+0.1000** | 0.2187 | **+0.1281** | 0.0965 |
| A549 | `mse_capped_norm` | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 |
| A549 | `lfc_nmae` | **-0.3954** | 0.0423 | **-0.3648** | 0.0723 | **-1.4362** | 0.0410 |
| A549 | `fid` | **+0.0000** | 0.0000 | **+0.2429** | 0.0650 | **+1.0665** | 0.0017 |
| A549 | `reach` | **-0.3510** | 0.1575 | **-0.1382** | 0.2337 | **-0.2638** | 0.0810 |
| A549 | `jac` | **-0.2025** | 0.5727 | **-2.8346** | 1.8869 | **-19.0726** | 0.0042 |
| BXPC3 | `pds_cosine` | **+0.1063** | 0.0892 | **+0.0623** | 0.1141 | **+0.0299** | 0.0327 |
| BXPC3 | `mse_capped_norm` | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 |
| BXPC3 | `lfc_nmae` | **-0.7241** | 0.0249 | **-0.7083** | 0.0195 | **-1.9442** | 0.1136 |
| BXPC3 | `fid` | **-0.0011** | 0.0029 | **+0.1068** | 0.0137 | **+0.6868** | 0.0011 |
| BXPC3 | `reach` | **-0.5432** | 0.0484 | **-0.4935** | 0.0544 | **-0.5368** | 0.0312 |
| BXPC3 | `jac` | **-0.0049** | 0.0000 | **+0.0021** | 0.0035 | **+0.0549** | 0.0006 |
| HAP1 | `pds_cosine` | **+0.2838** | 0.1526 | **+0.3026** | 0.1740 | **+0.1019** | 0.0821 |
| HAP1 | `mse_capped_norm` | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 |
| HAP1 | `lfc_nmae` | **-0.2079** | 0.0638 | **-0.2193** | 0.0582 | **-0.7978** | 0.0759 |
| HAP1 | `fid` | **+0.0154** | 0.0327 | **+0.0574** | 0.0593 | **+0.8066** | 0.0013 |
| HAP1 | `reach` | **-0.5301** | 0.1456 | **-0.5336** | 0.2077 | **-0.2616** | 0.0943 |
| HAP1 | `jac` | **-0.2588** | 0.3571 | **-1.0334** | 0.3706 | **-8.2329** | 0.0010 |
| HT29 | `pds_cosine` | **+0.0088** | 0.1779 | **-0.0622** | 0.1614 | **-0.0412** | 0.0685 |
| HT29 | `mse_capped_norm` | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 |
| HT29 | `lfc_nmae` | **-0.3973** | 0.0799 | **-0.4264** | 0.0323 | **-1.7634** | 0.1682 |
| HT29 | `fid` | **-0.3463** | 0.0340 | **-0.0198** | 0.0681 | **+1.0090** | 0.0014 |
| HT29 | `reach` | **-0.4126** | 0.1615 | **-0.5013** | 0.1804 | **-0.5250** | 0.1123 |
| HT29 | `jac` | **+0.4803** | 0.0592 | **+0.2988** | 0.1176 | **-0.2967** | 0.0007 |
| K562 | `pds_cosine` | **-0.0345** | 0.1464 | **-0.0267** | 0.1945 | **-0.0027** | 0.0734 |
| K562 | `mse_capped_norm` | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 |
| K562 | `lfc_nmae` | **-0.4803** | 0.0463 | **-0.4245** | 0.0514 | **-2.0643** | 0.1749 |
| K562 | `fid` | **-0.1659** | 0.0313 | **-0.0381** | 0.0398 | **+0.8195** | 0.0010 |
| K562 | `reach` | **-0.4361** | 0.1520 | **-0.3692** | 0.1736 | **-0.4402** | 0.0805 |
| K562 | `jac` | **+0.0954** | 0.0386 | **+0.0469** | 0.0561 | **-0.6435** | 0.0002 |
| MCF7 | `pds_cosine` | **+0.0523** | 0.1267 | **+0.0944** | 0.1290 | **+0.0159** | 0.0634 |
| MCF7 | `mse_capped_norm` | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 | **+0.0000** | 0.0000 |
| MCF7 | `lfc_nmae` | **-0.8778** | 0.0362 | **-0.8398** | 0.0491 | **-5.2308** | 0.2283 |
| MCF7 | `fid` | **-0.0003** | 0.0001 | **+0.0183** | 0.0097 | **+0.8895** | 0.0014 |
| MCF7 | `reach` | **-0.1281** | 0.0151 | **-0.1378** | 0.0157 | **-0.1577** | 0.0046 |
| MCF7 | `jac` | **-0.0004** | 0.0000 | **+0.0039** | 0.0031 | **+0.3338** | 0.0029 |

### 4.1 `avg_score` por linea y por via

| linea | **SUELO** | sd | **CAND01** | sd | **GWPS** | sd |
|---|---|---|---|---|---|---|
| A549 | **-0.1486** | 0.0792 | **-0.4991** | 0.3466 | **-3.2630** | 0.0230 |
| BXPC3 | **-0.1945** | 0.0217 | **-0.1718** | 0.0212 | **-0.2849** | 0.0246 |
| HAP1 | **-0.1163** | 0.0685 | **-0.2377** | 0.0814 | **-1.3973** | 0.0255 |
| HT29 | **-0.1112** | 0.0276 | **-0.1185** | 0.0434 | **-0.2696** | 0.0406 |
| K562 | **-0.1702** | 0.0389 | **-0.1353** | 0.0417 | **-0.3885** | 0.0340 |
| MCF7 | **-0.1591** | 0.0249 | **-0.1435** | 0.0306 | **-0.6915** | 0.0448 |

**`avg_score` es la media aritmetica de los seis `s`**, comprobado: A549/GWPS da
`(0.1281 + 0 - 1.4362 + 1.0665 - 0.2638 - 19.0726)/6 = -3.2630`, el valor
reportado.

---

## 5. 🔴 EL RANGO por metrica sobre las seis lineas: `[min, max]`

**Este es el entregable.** No el promedio.

| metrica | **SUELO [min, max]** | **CAND01 [min, max]** | **GWPS [min, max]** |
|---|---|---|---|
| `pds_cosine` | **[-0.0345, +0.2838]** | **[-0.0622, +0.3026]** | **[-0.0412, +0.1281]** |
| `mse_capped_norm` | **[+0.0000, +0.0000]** | **[+0.0000, +0.0000]** | **[+0.0000, +0.0000]** |
| `lfc_nmae` | **[-0.8778, -0.2079]** | **[-0.8398, -0.2193]** | **[-5.2308, -0.7978]** |
| `fid` | **[-0.3463, +0.0154]** | **[-0.0381, +0.2429]** | **[+0.6868, +1.0665]** |
| `reach` | **[-0.5432, -0.1281]** | **[-0.5336, -0.1378]** | **[-0.5368, -0.1577]** |
| `jac` | **[-0.2588, +0.4803]** | **[-2.8346, +0.2988]** | **[-19.0726, +0.3338]** |
| `avg_score` | **[-0.1945, -0.1112]** | **[-0.4991, -0.1185]** | **[-3.2630, -0.2696]** |

**Anchos de rango** (max - min): `fid` SUELO 0.3617, CAND01 0.2810, GWPS 0.3797;
`jac` SUELO 0.7391, CAND01 3.1334, **GWPS 19.4064**; `avg_score` SUELO 0.0833,
CAND01 0.3806, **GWPS 2.9934**.

`mse_capped_norm` es **+0.0000 exacto** en las tres vias y las seis lineas.

---

## 6. Los conteos contra los dos referentes

### 6.1 Contra el `b` (`s > 0`), con el margen absoluto de tres terminos

El margen de `dispersion_b.md` §4 (`margen COMPLETO`, celda por celda) es el que
corresponde a situar **una** `s` en la escala absoluta.

| metrica | SUELO `s>0` | SUELO `s>`margen | CAND01 `s>0` | CAND01 `s>`margen | **GWPS `s>0`** | **GWPS `s>`margen** |
|---|---|---|---|---|---|---|
| `pds_cosine` | 5/6 | 2/6 | 4/6 | 1/6 | 4/6 | **0/6** |
| `mse_capped_norm` | 0/6 | 0/6 | 0/6 | 0/6 | 0/6 | **0/6** |
| `lfc_nmae` | 0/6 | 0/6 | 0/6 | 0/6 | 0/6 | **0/6** |
| `fid` | 1/6 | 0/6 | 4/6 | 4/6 | **6/6** | **6/6** |
| `reach` | 0/6 | 0/6 | 0/6 | 0/6 | 0/6 | **0/6** |
| `jac` | 2/6 | 2/6 | 4/6 | 4/6 | 2/6 | **2/6** |

### 6.2 Contra el SUELO, con `resolucion_k`

`resolucion_8 = sqrt( (sd_a^2 + sd_b^2) / 8 )`, que es el error estandar de la
**diferencia** de dos medias de 8. Se reduce a `sd*sqrt(2/k)` de
`dispersion_b.md` §5 cuando las dos `sd` son iguales; aqui no lo son, asi que se
usa la forma general y se declara.

| metrica | CAND01 > SUELO | de ellas > `resolucion_8` | **GWPS > SUELO** | **de ellas > `resolucion_8`** |
|---|---|---|---|---|
| `pds_cosine` | 4/6 | 0/6 | 2/6 | **0/6** |
| `mse_capped_norm` | 0/6 | 0/6 | 0/6 | **0/6** |
| `lfc_nmae` | 4/6 | 4/6 | 0/6 | **0/6** |
| `fid` | 6/6 | 6/6 | **6/6** | **6/6** |
| `reach` | 3/6 | 2/6 | 3/6 | **2/6** |
| `jac` | 2/6 | 2/6 | 2/6 | **2/6** |

### 6.3 GWPS - SUELO celda por celda

| linea | metrica | SUELO | GWPS | **GWPS - SUELO** | `resolucion_8` | \|dif\|>res | margen abs. §4 |
|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | +0.0573 | +0.1281 | **+0.0708** | 0.1061 | no | 0.3013 |
| A549 | `mse_capped_norm` | +0.0000 | +0.0000 | **+0.0000** | 0.0000 | no | 0.0743 |
| A549 | `lfc_nmae` | -0.3954 | -1.4362 | **-1.0408** | 0.0208 | si | 0.0811 |
| A549 | `fid` | +0.0000 | +1.0665 | **+1.0665** | 0.0006 | si | 0.0891 |
| A549 | `reach` | -0.3510 | -0.2638 | **+0.0872** | 0.0626 | si | 0.3137 |
| A549 | `jac` | -0.2025 | -19.0726 | **-18.8702** | 0.2025 | si | 0.9115 |
| BXPC3 | `pds_cosine` | +0.1063 | +0.0299 | **-0.0764** | 0.0336 | si | 0.0987 |
| BXPC3 | `mse_capped_norm` | +0.0000 | +0.0000 | **+0.0000** | 0.0000 | no | 0.0518 |
| BXPC3 | `lfc_nmae` | -0.7241 | -1.9442 | **-1.2201** | 0.0411 | si | 0.0736 |
| BXPC3 | `fid` | -0.0011 | +0.6868 | **+0.6879** | 0.0011 | si | 0.0046 |
| BXPC3 | `reach` | -0.5432 | -0.5368 | **+0.0064** | 0.0204 | no | 0.2019 |
| BXPC3 | `jac` | -0.0049 | +0.0549 | **+0.0598** | 0.0002 | si | 0.0010 |
| HAP1 | `pds_cosine` | +0.2838 | +0.1019 | **-0.1819** | 0.0613 | si | 0.1682 |
| HAP1 | `mse_capped_norm` | +0.0000 | +0.0000 | **+0.0000** | 0.0000 | no | 0.0794 |
| HAP1 | `lfc_nmae` | -0.2079 | -0.7978 | **-0.5898** | 0.0350 | si | 0.0727 |
| HAP1 | `fid` | +0.0154 | +0.8066 | **+0.7912** | 0.0116 | si | 0.0515 |
| HAP1 | `reach` | -0.5301 | -0.2616 | **+0.2685** | 0.0613 | si | 0.3207 |
| HAP1 | `jac` | -0.2588 | -8.2329 | **-7.9741** | 0.1263 | si | 0.5820 |
| HT29 | `pds_cosine` | +0.0088 | -0.0412 | **-0.0500** | 0.0674 | no | 0.2158 |
| HT29 | `mse_capped_norm` | +0.0000 | +0.0000 | **+0.0000** | 0.0000 | no | 0.0353 |
| HT29 | `lfc_nmae` | -0.3973 | -1.7634 | **-1.3661** | 0.0658 | si | 0.1193 |
| HT29 | `fid` | -0.3463 | +1.0090 | **+1.3553** | 0.0120 | si | 0.1760 |
| HT29 | `reach` | -0.4126 | -0.5250 | **-0.1125** | 0.0695 | si | 0.2194 |
| HT29 | `jac` | +0.4803 | -0.2967 | **-0.7770** | 0.0209 | si | 0.1339 |
| K562 | `pds_cosine` | -0.0345 | -0.0027 | **+0.0318** | 0.0579 | no | 0.1853 |
| K562 | `mse_capped_norm` | +0.0000 | +0.0000 | **+0.0000** | 0.0000 | no | 0.0404 |
| K562 | `lfc_nmae` | -0.4803 | -2.0643 | **-1.5840** | 0.0640 | si | 0.0983 |
| K562 | `fid` | -0.1659 | +0.8195 | **+0.9854** | 0.0111 | si | 0.0322 |
| K562 | `reach` | -0.4361 | -0.4402 | **-0.0041** | 0.0608 | no | 0.2810 |
| K562 | `jac` | +0.0954 | -0.6435 | **-0.7389** | 0.0136 | si | 0.0399 |
| MCF7 | `pds_cosine` | +0.0523 | +0.0159 | **-0.0364** | 0.0501 | no | 0.1456 |
| MCF7 | `mse_capped_norm` | +0.0000 | +0.0000 | **+0.0000** | 0.0000 | no | 0.0395 |
| MCF7 | `lfc_nmae` | -0.8778 | -5.2308 | **-4.3530** | 0.0817 | si | 0.0836 |
| MCF7 | `fid` | -0.0003 | +0.8895 | **+0.8898** | 0.0005 | si | 0.0019 |
| MCF7 | `reach` | -0.1281 | -0.1577 | **-0.0296** | 0.0056 | si | 0.0432 |
| MCF7 | `jac` | -0.0004 | +0.3338 | **+0.3342** | 0.0010 | si | 0.0021 |

---

## 7. 🔴 Cuanto de cada linea viene de los 6 blancos REPLEGADOS al suelo

El agregado del marcador es la **media sobre perturbaciones**, asi que la media
con y sin esos 6 es una descomposicion **exacta a nivel `u`**. 🔴 **No a nivel
`s`**, que es un reescalado no lineal por celda: se declara y no se convierte.

Media sobre las 8 semillas. `aporte` = `n_replegados_con_valor * u_solo_6 / n_todas`.

| linea | metrica | `u` todas (33) | `u` sin los 6 (27) | `u` solo los 6 | n de los 6 con valor | **aporte de los 6** |
|---|---|---|---|---|---|---|
| A549 | `jac` | 0.0977 | 0.0083 | 0.5000 | 6.0 | **0.0909 (93 % del total)** |
| A549 | `fid` | 0.4457 | 0.4952 | 0.0000 | 3.0 | 0.0000 (0 %) |
| A549 | `lfc_nmae` | 1.4324 | 1.4899 | 1.0302 | 1.0 | 0.1288 (9 %) |
| A549 | `reach` | 0.2054 | 0.2075 | 0.1948 | 3.0 | 0.0325 (16 %) |
| BXPC3 | `jac` | 0.0119 | 0.0146 | 0.0000 | 6.0 | 0.0000 (0 %) |
| BXPC3 | `fid` | 0.4120 | 0.5036 | 0.0000 | 6.0 | 0.0000 (0 %) |
| BXPC3 | `lfc_nmae` | 1.3350 | 1.4063 | 0.9964 | 4.0 | 0.1733 (13 %) |
| BXPC3 | `reach` | 0.0817 | 0.0821 | 0.0800 | 6.0 | 0.0145 (18 %) |
| HAP1 | `jac` | 0.0627 | 0.0026 | 0.3333 | 6.0 | **0.0606 (97 % del total)** |
| HAP1 | `fid` | 0.4411 | 0.5065 | 0.0000 | 4.0 | 0.0000 (0 %) |
| HAP1 | `lfc_nmae` | 1.3007 | 1.3007 | n/d | 0.0 | n/d |
| HAP1 | `reach` | 0.2575 | 0.2317 | 0.3542 | 4.0 | 0.0746 (29 %) |
| HT29 | `jac` | 0.0374 | 0.0086 | 0.1667 | 6.0 | **0.0303 (81 % del total)** |
| HT29 | `fid` | 0.4347 | 0.5152 | 0.0000 | 5.0 | 0.0000 (0 %) |
| HT29 | `lfc_nmae` | 1.5377 | 1.6117 | 1.0568 | 2.0 | 0.1409 (9 %) |
| HT29 | `reach` | 0.1748 | 0.1377 | 0.3381 | 5.0 | 0.0626 (36 %) |
| K562 | `jac` | 0.0333 | 0.0037 | 0.1667 | 6.0 | **0.0303 (91 % del total)** |
| K562 | `fid` | 0.4332 | 0.5134 | 0.0000 | 5.0 | 0.0000 (0 %) |
| K562 | `lfc_nmae` | 1.5674 | 1.5674 | n/d | 0.0 | n/d |
| K562 | `reach` | 0.2296 | 0.1839 | 0.4125 | 5.0 | 0.0825 (36 %) |
| MCF7 | `jac` | 0.0241 | 0.0295 | 0.0000 | 6.0 | 0.0000 (0 %) |
| MCF7 | `fid` | 0.4115 | 0.5030 | 0.0000 | 6.0 | 0.0000 (0 %) |
| MCF7 | `lfc_nmae` | 1.9621 | 2.1625 | 1.0604 | 6.0 | 0.1928 (10 %) |
| MCF7 | `reach` | 0.0113 | 0.0113 | 0.0114 | 6.0 | 0.0021 (18 %) |

En `fid` los 6 replegados valen **0.0000** donde tienen valor, y quitarlos SUBE
el `u` de 0.41-0.45 a 0.50-0.52 en las seis lineas. En `jac` los 6 valen 0.167 a
0.500 contra 0.003-0.030 de los otros 27.

---

## 8. 🔴 La prediccion H3 contra lo que salio

| linea | **span `jac` (r-b)** | `s` `jac` GWPS | `s` `jac` SUELO | `s` `jac` CAND01 | aporte de `jac` al `avg` GWPS | % del `avg` |
|---|---|---|---|---|---|---|
| **A549** | **0.018709** | **-19.0726** | -0.2025 | -2.8346 | **-3.1788** | **97 %** |
| **HAP1** | **0.043914** | **-8.2329** | -0.2588 | -1.0334 | **-1.3721** | **98 %** |
| MCF7 | 0.072214 | +0.3338 | -0.0004 | +0.0039 | +0.0556 | -8 % |
| HT29 | 0.181026 | -0.2967 | +0.4803 | +0.2988 | -0.0495 | 18 % |
| BXPC3 | 0.199792 | +0.0549 | -0.0049 | +0.0021 | +0.0091 | -3 % |
| K562 | 0.277874 | -0.6435 | +0.0954 | +0.0469 | -0.1072 | 28 % |

Ordenado por span. Las **dos lineas de span mas chico** son las dos con
`|s|` de `jac` de dos digitos y las dos donde `jac` aporta **97 % y 98 %** del
`avg_score`. En las otras cuatro el aporte de `jac` va de **-8 % a 28 %**.
A549 tiene span `0.018709`, y `1/33 = 0.030303`.

**Las otras celdas, leidas** (§4, §5, §6.1): `fid` da **6/6** por encima del `b`
y **6/6** por encima del margen absoluto, con rango **[+0.6868, +1.0665]** y dos
lineas por encima de 1.0 (A549 +1.0665, HT29 +1.0090). `lfc_nmae` da **0/6** y es
peor que el suelo en **6/6** con `|dif| > resolucion_8` en las seis.
`mse_capped_norm` es 0.0000 exacto. `pds_cosine` da 4/6 positivas y **0/6** por
encima del margen. `reach` da 0/6 positivas en las tres vias.

🔴 **Y lo que la prediccion mandaba no hacer, no se hizo.** No se bajo la
emision: las 48 corridas llaman **66 180 a 98 172 genes DE** en total (§9.1), y
las unicas uniones vacias son de los **6 blancos replegados al suelo** (§9.1),
que son los que no tienen perfil de gwps. No se toco ninguna constante despues de
ver un `s`.

---

## 9. Las dos cohortes moviles, y los defectos de la construccion

### 9.1 Cohortes moviles por semilla

| linea | `n_fid` media | `n_fid` rango | union no vacia media | rango | `E` union vacia | replegados con union vacia | DEpred suma |
|---|---|---|---|---|---|---|---|
| A549 | 30.00 | [30, 30] | 30.00 | [30, 30] | 3.00 | `H1F0`, `HERC5`, `HERC6` | 66 180 |
| BXPC3 | 33.00 | [33, 33] | 33.00 | [33, 33] | 0.00 | ninguno | 86 139 |
| HAP1 | 31.00 | [31, 31] | 31.00 | [31, 31] | 2.00 | `H1F0`, `MAP3K14` | 98 172 |
| HT29 | 32.00 | [32, 32] | 32.00 | [32, 32] | 1.00 | `HERC5` | 83 639 |
| K562 | 32.00 | [32, 32] | 32.00 | [32, 32] | 1.00 | `ZBP1` | 85 931 |
| MCF7 | 33.00 | [33, 33] | 33.00 | [33, 33] | 0.00 | ninguno | 96 536 |

Las dos cohortes **no se mueven entre semillas**: rango de ancho 0 en las seis
lineas. **Todas las uniones vacias son de blancos replegados** (`E` = 3, 2, 1, 1,
0, 0), y ninguna de los 27 con perfil de gwps.

### 9.2 Los ratios no dependen del contexto

Los `r_h` salen solo de gwps, asi que son **identicos en las seis lineas**. Esta
emision no estima nada del pool de control mas que las celulas que remuestrea.
Es una propiedad de la via, declarada en §2.6, y es la razon de que `sd` de GWPS
sea mucho mas chica que la de las otras dos vias en `fid` y `jac` (0.0002 a
0.0042 contra 0.0035 a 1.8869 de CAND01).

### 9.3 La variante con puerta de significancia, medida y NO puntuada

BH al 0.05 del preset sobre `z = p_g / sd_nulo`. Los 27 blancos con perfil:

| # | blanco | constr. | genes que pasan | # | blanco | constr. | genes que pasan |
|---|---|---|---|---|---|---|---|
| 1 | `AKT1` | P1P2 | 2 | 15 | `MAPK8` | P1P2 | 1 |
| 2 | `BATF2` | P1P2 | **0** | 16 | `PARP12` | P1P2 | 14 |
| 3 | `BRD9` | P1P2 | **127** | 17 | `RAP1GAP` | P1 | **0** |
| 4 | `CEBPG` | P1P2 | **829** | 18 | `RNF114` | P1P2 | 3 |
| 5 | `ELK1` | P1P2 | 3 | 19 | `SOCS1` | P1P2 | 3 |
| 6 | `ETS1` | P1P2 | 1 | 20 | `SP100` | P1P2 | **0** |
| 7 | `ETS2` | P1P2 | **630** | 21 | `STAT1` | P1P2 | 3 |
| 8 | `ID2` | P1 | 2 | 22 | `STAT3` | P1P2 | **0** |
| 9 | `IFNAR1` | P1P2 | 5 | 23 | `STAT4` | P1P2 | 1 |
| 10 | `IFNAR2` | P1P2 | 1 | 24 | `STAT5A` | P1P2 | 1 |
| 11 | `IRF1` | P1P2 | **0** | 25 | `TRAFD1` | P1P2 | 1 |
| 12 | `IRF7` | P1P2 | **0** | 26 | `TYK2` | P1P2 | **96** |
| 13 | `JAK1` | P1P2 | 11 | 27 | `VAV1` | P1P2 | **0** |
| 14 | `JUN` | P1P2 | 1 | | | | |

**De 0 a 829 genes de 8248, mediana 1, 7 de 27 blancos en cero**, y solo 4 por
encima de 90. **Esta variante no se puntuo** y queda declarada como no medida.

### 9.4 Defectos de la construccion, como datos sobre la medicion

1. **La magnitud es una eleccion de escala, no una medicion** (§2.2). El rango
   normalizado reparte `|s|` con **mediana 0.500** sobre 8124 genes por blanco:
   cada gen seleccionado cambia su esperanza en ~50 %. Ninguna cifra de gwps
   fija esa escala; solo fija el orden y el signo.
2. **Los 6 blancos replegados no son una emision de la via** (§2.4, §7). En
   `jac` aportan 81-97 % del `u` en cuatro lineas.
3. **Un solo perfil por construccion en gwps**: no hay replicas dentro de una
   construccion, asi que el ruido de gwps entra sin descontar. Las unicas dos
   excepciones son `ID2` y `RAP1GAP` (§10), y **no se promediaron**.
4. **8124 de 34025 genes de Jiang reciben `r != 1`**: el 76 % del universo de
   Jiang queda en `r = 1` por no estar en gwps, y eso no se distingue de "gwps
   dice que no se mueve".
5. **`mse_capped_norm` es +0.0000 exacto** en las 18 celdas (3 vias x 6 lineas).
   No discrimina entre estas vias.
6. **El margen de §4 y la `resolucion_8` de §6.2 responden preguntas distintas**
   y no son intercambiables: el primero situa una `s` en la escala absoluta, la
   segunda compara dos vias. Se usan separados.
7. **La `sd` de GWPS en `jac` y `fid` es diminuta** (0.0002-0.0042) porque la
   estructura no depende del contexto (§9.2), asi que `resolucion_8` en esas
   celdas es casi solo la `sd` de la otra via.
8. **No se midio** la via con la puerta de §9.3, ni con P2 en las 8 semillas, ni
   con mas de 8 semillas, ni ninguna variante de escala distinta de la de §2.2.

---

## 10. P1 contra P2 en `ID2` y `RAP1GAP`, sin promediar

Una semilla (`j = 0`, la misma en las dos: `20271000 + 100*fold`), cambiando
**solo** la fila de gwps de esos dos blancos.

| linea | `avg` con **P1** | `avg` con **P2** | dif | DEpred P1 | DEpred P2 |
|---|---|---|---|---|---|
| A549 | -3.2814 | -3.2484 | **+0.0330** | 66 224 | 67 388 |
| BXPC3 | -0.2928 | -0.2998 | **-0.0070** | 85 794 | 87 191 |
| HAP1 | -1.3816 | -1.4108 | **-0.0293** | 98 151 | 99 270 |
| HT29 | -0.3361 | -0.2633 | **+0.0729** | 83 397 | 84 247 |
| K562 | -0.4050 | -0.3979 | **+0.0071** | 86 239 | 86 986 |
| MCF7 | -0.6598 | -0.7398 | **-0.0800** | 96 681 | 97 430 |

Por metrica, `s` con P2 menos `s` con P1:

| linea | `pds_cosine` | `mse_capped_norm` | `lfc_nmae` | `fid` | `reach` | `jac` |
|---|---|---|---|---|---|---|
| A549 | +0.1236 | +0.0000 | -0.0664 | -0.0042 | +0.1426 | +0.0025 |
| BXPC3 | +0.0910 | +0.0000 | -0.0839 | -0.0006 | -0.0475 | -0.0008 |
| HAP1 | -0.0269 | +0.0000 | -0.0214 | +0.0014 | -0.1285 | -0.0003 |
| HT29 | +0.1547 | +0.0000 | +0.3083 | -0.0025 | -0.0217 | -0.0016 |
| K562 | +0.2039 | +0.0000 | -0.0654 | +0.0018 | -0.0980 | +0.0002 |
| MCF7 | -0.1183 | +0.0000 | -0.3492 | -0.0012 | -0.0045 | -0.0068 |

El cambio de dos filas de 33 mueve el `avg` entre **-0.0800 y +0.0729**, y en
`pds_cosine` hasta **+0.2039**. **No se promedian**, y la corrida principal (§4,
§5, §6) es la de **P1**.

---

## 11. Lo que no se hizo, y las entradas

| prohibicion | estado |
|---|---|
| usar cifras de `~/vcc/benchmark/` o `/datos/benchmark_VEDADO_*` | **no leidos**. El perfil sale del h5ad crudo en `/datos/ref/replogle/`; los bundles de puntuacion son `/datos/modelo/out/nativas_comp/`, que es la ruta que el protocolo de `candidato_01.md` y `dispersion_u.md` ya usaba. El dato de P1 de `gwps_fiabilidad.md` (67.6 % de 546 pares) se uso **solo** como la instruccion de no promediar, tal como venia en el encargo; no entra en ningun calculo |
| ajustar nada mirando el `s` | **no**. La escala de §2.2 y la decision de §2.3 se fijaron y se declararon antes de la primera puntuacion; la variante con puerta no se puntuo |
| promediar P1 y P2 | **no**. §10 las da separadas |
| modificar el arnes | **no**. `suelo_arnes.py` y `cand01_covar._reescala_entero` se importan sin tocar |
| escribir veredicto o conclusiones | **no** |
| envios, `vcc sample`, commits | **ninguno** |

| | |
|---|---|
| entradas | `/datos/ref/replogle/k562_gwps_raw_bulk.h5ad`; `/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad`; `modelo/out/suelo_panel.json`; bundles `/datos/modelo/out/nativas_comp/BUNDLE_<linea>`; `modelo/out/dispersion_u.json`, `modelo/out/candidato_01{,_p2,_p3}.json` (vias ya medidas); `modelo/out/dispersion_b.md` §4 y §5 (margen y resolucion) |
| codigo | `modelo/src/gwps_emision.py`, `modelo/src/gwps_via.py`, `modelo/src/tres_vias_tablas.py` |
| salidas | `modelo/out/gwps_via.json` (48), `modelo/out/gwps_via_p2.json` (6), `modelo/out/gwps_via_run.log`, `modelo/out/gwps_via_p2_run.log` |
| semillas | emision `20271000 + 100*fold + j`, `j = 0..7`; poblacion real `20260001 + fold` |
| verificaciones | `score_ok` 48/48, `config_hash` 48/48, `poblacion_real_identica` 48/48, 1 huella real por linea, 48 huellas `pred` distintas, `todo_entero` 48/48, 0 negativos |
