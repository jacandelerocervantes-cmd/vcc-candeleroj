# El control permutado de la via de estructura externa

## 🔴 0. Que es este reporte y que no es

**Presenta datos.** Declara la construccion antes de medir, da las tablas con sus
controles y reporta los defectos de la propia construccion. **No lleva veredicto,
ni conclusiones, ni nada dirigido al siguiente candidato.**

🔴 **LA PREGUNTA QUE ESTE BRAZO MIDE.** `fid` pregunta: *de los genes que de
verdad se mueven, ¿acerte el signo?* **No** pregunta si acerte el signo **de ese
blanco**. Si gwps emite una direccion **compartida** entre los 33, `fid` la premia
igual. Este brazo mide cuanto del `fid` sobrevive cuando el perfil se asigna al
blanco **equivocado**. La referencia del otro lado es `pds_cosine`, que es
invariante de escala y pregunta "¿es el blanco correcto entre 33?", y que en
`tres_vias_rango.md` §6.1 dio **0 de 6** por encima del margen.

---

## 1. Construccion: un DESARREGLO, no una emision nueva

### 1.1 Que cambia y que no

| | |
|---|---|
| **lo unico que cambia** | para la perturbacion `h` se emite el perfil de gwps de **OTRO** blanco `g != h` |
| el desarreglo | **solo sobre los 27 blancos CON perfil** de gwps |
| 🔴 los 6 replegados | `DTX3L`, `H1F0`, `HERC5`, `HERC6`, `MAP3K14`, `ZBP1` siguen **replegados al suelo en LOS DOS BRAZOS** (`r = 1`), asi que su aporte es **identico** y la comparacion queda limpia sobre los 27 |
| 🔴 la regla de escala | **SIN TOCAR**. La misma de `tres_vias_rango.md` §2.2: orden y signo, `s` por rango normalizado dentro de cada sentido, `r = 1 + s` en `[0, 2]`. Se llama a `gwps_emision.construye_ratios` sin modificarla y solo se **reasigna** que vector recibe cada etiqueta |
| 🔴 el muestreo entero | **SIN TOCAR**. El mismo `cand01_covar._reescala_entero` (Binomial para `r<=1`, Poisson aditivo para `r>1`) |
| 🔴 el arnes | **SIN TOCAR**. `suelo_arnes.build_fold_adatas`, `predict_fn(control_pool_X, labels, rng)` |
| 🔴 la poblacion real | la misma, `seed_real = 20260001 + fold` |
| 🔴 el bundle | el mismo de competencia, `/datos/modelo/out/nativas_comp/BUNDLE_<linea>` |
| P1/P2 | **P1 en los dos brazos**, igual que la corrida principal. No se promedia |

### 1.2 🔴 Por que este control y no otro

El brazo permutado emite el **mismo numero de genes** con la **misma
distribucion de magnitud**: es el mismo conjunto de 27 vectores `r`, solo
reasignado. **Esta igualado en potencia por construccion**, y eso se comprueba
medido en §8: `DEpred` difiere entre brazos en **-342 a +110 genes** sobre
66 000-98 000 (0.03 % a 0.4 %), y la cohorte `n_fid` es **identica** en las seis
lineas.

### 1.3 Semillas

| | |
|---|---|
| desarreglo | **`20272500 + j`**, `j = 0..7`. 🔴 **Los mismos 8 desarreglos en las seis lineas**, para que las lineas sean comparables |
| remuestreo de la emision | familia **nueva** **`20272000 + 100*fold + j`**. Sin colision con 2026000x-2026600x ni 2027100x, **pero SI con la familia del desarreglo de esta misma corrida**: ver abajo |
| poblacion real | `20260001 + fold` |

Los desarreglos se generan por **rechazo**: se sortea una permutacion de los 27 y
se descarta si deja algun punto fijo.

🔴 **UNA COLISION DE SEMILLA, medida y declarada.** `MCF7` es `fold = 5`, asi que
su semilla de remuestreo es `20272000 + 500 + j = 20272500 + j`, que es
**exactamente la familia de la semilla del desarreglo**. Contado en el JSON:
`seed_emision == seed_desarreglo` en **8 de las 48 corridas, las 8 de MCF7**, y en
**0 de las 40** de las otras cinco lineas. El brazo REAL no la tiene
(`20271000 + 100*fold + j` contra `20272500 + j`: **0 colisiones**), y las dos
familias de emision de los dos brazos **no colisionan entre si** (0).

Lo que eso hace y lo que no: son **dos objetos `default_rng` distintos**, uno que
consume permutaciones de 27 y otro que consume indices de remuestreo y sorteos
Binomial/Poisson, pero **arrancan del mismo estado**. En MCF7, y solo en MCF7, el
sorteo de remuestreo de cada `j` es una funcion determinista del **mismo entero**
que fijo el desarreglo `sigma_j`, asi que **el remuestreo de MCF7 no se sorteo con
independencia de la familia de desarreglos**. No cambia ninguna cifra ya medida
--las 8 corridas de MCF7 son las que son, y su `sd` esta en §2 y §9-- y no afecta
a `sigma_j`, que es **identico en las seis lineas** (verificado: una sola firma
por `j`, §1.4). Queda como **defecto declarado** en §10.8, y la frase "familia
nueva sin colision" **solo vale para las otras cinco lineas**.

### 1.4 Los 8 desarreglos: ninguno tiene punto fijo

| `j` | semilla | intentos hasta el desarreglo | puntos fijos | `STAT1` recibe | `IFNAR2` recibe |
|---|---|---|---|---|---|
| 0 | 20272500 | 6 | **0** | `ETS1` | `STAT3` |
| 1 | 20272501 | 1 | **0** | `JAK1` | `STAT3` |
| 2 | 20272502 | 2 | **0** | `ETS1` | `BATF2` |
| 3 | 20272503 | 1 | **0** | `AKT1` | `ID2` |
| 4 | 20272504 | 4 | **0** | `IFNAR2` | `JUN` |
| 5 | 20272505 | 6 | **0** | `RNF114` | `IRF7` |
| 6 | 20272506 | 1 | **0** | `SOCS1` | `RNF114` |
| 7 | 20272507 | 7 | **0** | `ELK1` | `ID2` |

**`sin_punto_fijo` verificado en las 48 corridas**, y el `sigma` de cada `j` es
**identico en las seis lineas** (una sola firma por `j`).

### 1.5 Verificaciones de las 48 corridas

| | |
|---|---|
| `score_ok` | **48/48** |
| `config_hash` de competencia | **48/48** coincide con `competition_payload()['config_hash']`, y es **un solo digest en las 48**: `a5519017556705abf9a6e553ebb4e517b9612e9546331ace7c5ec6c5c30a1a13`. El mismo digest en las 48 del brazo REAL (`gwps_via.json`): **48/48** |
| `poblacion_real_identica` (sha256) | **48/48**, **1 sola huella real por linea** |
| huellas `pred` distintas | **48/48** |
| huellas `pred` compartidas con el brazo REAL | **0** |
| `todo_entero` | **48/48** |
| valores negativos | **0** |
| 🔴 `seed_emision == seed_desarreglo` | **8 de 48** (las 8 de MCF7), 0 de 40 en las otras cinco lineas, 0 en el brazo REAL. Colision declarada en §1.3 y §10.8 |

---

## 2. Los seis `s` del brazo PERMUTADO

Media sobre los 8 desarreglos, con su `sd`.

| linea | `pds_cosine` | `mse_capped_norm` | `lfc_nmae` | `fid` | `reach` | `jac` | `avg_score` |
|---|---|---|---|---|---|---|---|
| A549 | **-0.0500** ± 0.0915 | +0.0000 ± 0.0000 | **-1.4712** ± 0.2786 | **+1.0672** ± 0.0025 | **-0.1414** ± 0.1482 | **-19.0853** ± 0.0136 | **-3.2801** ± 0.0542 |
| BXPC3 | **+0.0312** ± 0.0843 | +0.0000 ± 0.0000 | **-1.9096** ± 0.1218 | **+0.6857** ± 0.0021 | **-0.5388** ± 0.0572 | **+0.0529** ± 0.0017 | **-0.2798** ± 0.0263 |
| HAP1 | **+0.1849** ± 0.1460 | +0.0000 ± 0.0000 | **-0.8167** ± 0.1742 | **+0.8069** ± 0.0030 | **-0.5092** ± 0.1784 | **-8.2315** ± 0.0042 | **-1.4276** ± 0.0487 |
| HT29 | **-0.0219** ± 0.1475 | +0.0000 ± 0.0000 | **-1.5463** ± 0.1810 | **+1.0087** ± 0.0044 | **-0.6934** ± 0.1790 | **-0.2986** ± 0.0025 | **-0.2586** ± 0.0513 |
| K562 | **-0.0706** ± 0.1592 | +0.0000 ± 0.0000 | **-1.9997** ± 0.4068 | **+0.8129** ± 0.0034 | **-0.5837** ± 0.2064 | **-0.6437** ± 0.0004 | **-0.4141** ± 0.0954 |
| MCF7 | **+0.0140** ± 0.1092 | +0.0000 ± 0.0000 | **-5.6192** ± 0.2228 | **+0.8879** ± 0.0035 | **-0.1651** ± 0.0094 | **+0.3386** ± 0.0087 | **-0.7573** ± 0.0450 |

---

## 3. Los tres brazos lado a lado

SUELO de `dispersion_u.json` (ya medido), REAL de `gwps_via.json` (ya medido),
PERMUTADO de esta corrida. Los tres con 8 semillas y el mismo bundle por linea.

| linea | metrica | SUELO | sd | REAL | sd | PERMUTADO | sd |
|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | +0.0573 | 0.2840 | +0.1281 | 0.0965 | -0.0500 | 0.0915 |
| A549 | `mse_capped_norm` | +0.0000 | 0.0000 | +0.0000 | 0.0000 | +0.0000 | 0.0000 |
| A549 | `lfc_nmae` | -0.3954 | 0.0423 | -1.4362 | 0.0410 | -1.4712 | 0.2786 |
| A549 | `fid` | +0.0000 | 0.0000 | **+1.0665** | 0.0017 | **+1.0672** | 0.0025 |
| A549 | `reach` | -0.3510 | 0.1575 | -0.2638 | 0.0810 | -0.1414 | 0.1482 |
| A549 | `jac` | -0.2025 | 0.5727 | -19.0726 | 0.0042 | -19.0853 | 0.0136 |
| BXPC3 | `pds_cosine` | +0.1063 | 0.0892 | +0.0299 | 0.0327 | +0.0312 | 0.0843 |
| BXPC3 | `mse_capped_norm` | +0.0000 | 0.0000 | +0.0000 | 0.0000 | +0.0000 | 0.0000 |
| BXPC3 | `lfc_nmae` | -0.7241 | 0.0249 | -1.9442 | 0.1136 | -1.9096 | 0.1218 |
| BXPC3 | `fid` | -0.0011 | 0.0029 | **+0.6868** | 0.0011 | **+0.6857** | 0.0021 |
| BXPC3 | `reach` | -0.5432 | 0.0484 | -0.5368 | 0.0312 | -0.5388 | 0.0572 |
| BXPC3 | `jac` | -0.0049 | 0.0000 | +0.0549 | 0.0006 | +0.0529 | 0.0017 |
| HAP1 | `pds_cosine` | +0.2838 | 0.1526 | +0.1019 | 0.0821 | +0.1849 | 0.1460 |
| HAP1 | `mse_capped_norm` | +0.0000 | 0.0000 | +0.0000 | 0.0000 | +0.0000 | 0.0000 |
| HAP1 | `lfc_nmae` | -0.2079 | 0.0638 | -0.7978 | 0.0759 | -0.8167 | 0.1742 |
| HAP1 | `fid` | +0.0154 | 0.0327 | **+0.8066** | 0.0013 | **+0.8069** | 0.0030 |
| HAP1 | `reach` | -0.5301 | 0.1456 | -0.2616 | 0.0943 | -0.5092 | 0.1784 |
| HAP1 | `jac` | -0.2588 | 0.3571 | -8.2329 | 0.0010 | -8.2315 | 0.0042 |
| HT29 | `pds_cosine` | +0.0088 | 0.1779 | -0.0412 | 0.0685 | -0.0219 | 0.1475 |
| HT29 | `mse_capped_norm` | +0.0000 | 0.0000 | +0.0000 | 0.0000 | +0.0000 | 0.0000 |
| HT29 | `lfc_nmae` | -0.3973 | 0.0799 | -1.7634 | 0.1682 | -1.5463 | 0.1810 |
| HT29 | `fid` | -0.3463 | 0.0340 | **+1.0090** | 0.0014 | **+1.0087** | 0.0044 |
| HT29 | `reach` | -0.4126 | 0.1615 | -0.5250 | 0.1123 | -0.6934 | 0.1790 |
| HT29 | `jac` | +0.4803 | 0.0592 | -0.2967 | 0.0007 | -0.2986 | 0.0025 |
| K562 | `pds_cosine` | -0.0345 | 0.1464 | -0.0027 | 0.0734 | -0.0706 | 0.1592 |
| K562 | `mse_capped_norm` | +0.0000 | 0.0000 | +0.0000 | 0.0000 | +0.0000 | 0.0000 |
| K562 | `lfc_nmae` | -0.4803 | 0.0463 | -2.0643 | 0.1749 | -1.9997 | 0.4068 |
| K562 | `fid` | -0.1659 | 0.0313 | **+0.8195** | 0.0010 | **+0.8129** | 0.0034 |
| K562 | `reach` | -0.4361 | 0.1520 | -0.4402 | 0.0805 | -0.5837 | 0.2064 |
| K562 | `jac` | +0.0954 | 0.0386 | -0.6435 | 0.0002 | -0.6437 | 0.0004 |
| MCF7 | `pds_cosine` | +0.0523 | 0.1267 | +0.0159 | 0.0634 | +0.0140 | 0.1092 |
| MCF7 | `mse_capped_norm` | +0.0000 | 0.0000 | +0.0000 | 0.0000 | +0.0000 | 0.0000 |
| MCF7 | `lfc_nmae` | -0.8778 | 0.0362 | -5.2308 | 0.2283 | -5.6192 | 0.2228 |
| MCF7 | `fid` | -0.0003 | 0.0001 | **+0.8895** | 0.0014 | **+0.8879** | 0.0035 |
| MCF7 | `reach` | -0.1281 | 0.0151 | -0.1577 | 0.0046 | -0.1651 | 0.0094 |
| MCF7 | `jac` | -0.0004 | 0.0000 | +0.3338 | 0.0029 | +0.3386 | 0.0087 |
| A549 | `avg_score` | -0.1486 | 0.0792 | -3.2630 | 0.0230 | -3.2801 | 0.0542 |
| BXPC3 | `avg_score` | -0.1945 | 0.0217 | -0.2849 | 0.0246 | -0.2798 | 0.0263 |
| HAP1 | `avg_score` | -0.1163 | 0.0685 | -1.3973 | 0.0255 | -1.4276 | 0.0487 |
| HT29 | `avg_score` | -0.1112 | 0.0276 | -0.2696 | 0.0406 | -0.2586 | 0.0513 |
| K562 | `avg_score` | -0.1702 | 0.0389 | -0.3885 | 0.0340 | -0.4141 | 0.0954 |
| MCF7 | `avg_score` | -0.1591 | 0.0249 | -0.6915 | 0.0448 | -0.7573 | 0.0450 |

---

## 4. 🔴 REAL - PERMUTADO celda por celda

`resolucion_8 = sqrt( (sd_real^2 + sd_perm^2) / 8 )`, el error estandar de la
diferencia de dos medias de 8.

| linea | metrica | REAL | PERMUTADO | **REAL - PERM** | `resolucion_8` | \|dif\|>res |
|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | +0.1281 | -0.0500 | **+0.1781** | 0.0470 | **si** |
| A549 | `mse_capped_norm` | +0.0000 | +0.0000 | +0.0000 | 0.0000 | no |
| A549 | `lfc_nmae` | -1.4362 | -1.4712 | +0.0350 | 0.0996 | no |
| A549 | `fid` | +1.0665 | +1.0672 | **-0.0007** | 0.0011 | no |
| A549 | `reach` | -0.2638 | -0.1414 | -0.1225 | 0.0597 | **si** |
| A549 | `jac` | -19.0726 | -19.0853 | +0.0127 | 0.0050 | **si** |
| BXPC3 | `pds_cosine` | +0.0299 | +0.0312 | -0.0012 | 0.0320 | no |
| BXPC3 | `mse_capped_norm` | +0.0000 | +0.0000 | +0.0000 | 0.0000 | no |
| BXPC3 | `lfc_nmae` | -1.9442 | -1.9096 | -0.0346 | 0.0589 | no |
| BXPC3 | `fid` | +0.6868 | +0.6857 | **+0.0011** | 0.0008 | **si** |
| BXPC3 | `reach` | -0.5368 | -0.5388 | +0.0020 | 0.0231 | no |
| BXPC3 | `jac` | +0.0549 | +0.0529 | +0.0020 | 0.0007 | **si** |
| HAP1 | `pds_cosine` | +0.1019 | +0.1849 | -0.0830 | 0.0592 | **si** |
| HAP1 | `mse_capped_norm` | +0.0000 | +0.0000 | +0.0000 | 0.0000 | no |
| HAP1 | `lfc_nmae` | -0.7978 | -0.8167 | +0.0189 | 0.0672 | no |
| HAP1 | `fid` | +0.8066 | +0.8069 | **-0.0003** | 0.0011 | no |
| HAP1 | `reach` | -0.2616 | -0.5092 | +0.2476 | 0.0713 | **si** |
| HAP1 | `jac` | -8.2329 | -8.2315 | -0.0014 | 0.0015 | no |
| HT29 | `pds_cosine` | -0.0412 | -0.0219 | -0.0193 | 0.0575 | no |
| HT29 | `mse_capped_norm` | +0.0000 | +0.0000 | +0.0000 | 0.0000 | no |
| HT29 | `lfc_nmae` | -1.7634 | -1.5463 | -0.2171 | 0.0874 | **si** |
| HT29 | `fid` | +1.0090 | +1.0087 | **+0.0004** | 0.0016 | no |
| HT29 | `reach` | -0.5250 | -0.6934 | +0.1683 | 0.0747 | **si** |
| HT29 | `jac` | -0.2967 | -0.2986 | +0.0019 | 0.0009 | **si** |
| K562 | `pds_cosine` | -0.0027 | -0.0706 | +0.0678 | 0.0620 | **si** |
| K562 | `mse_capped_norm` | +0.0000 | +0.0000 | +0.0000 | 0.0000 | no |
| K562 | `lfc_nmae` | -2.0643 | -1.9997 | -0.0646 | 0.1565 | no |
| K562 | `fid` | +0.8195 | +0.8129 | **+0.0066** | 0.0013 | **si** |
| K562 | `reach` | -0.4402 | -0.5837 | +0.1434 | 0.0783 | **si** |
| K562 | `jac` | -0.6435 | -0.6437 | +0.0003 | 0.0002 | **si** |
| MCF7 | `pds_cosine` | +0.0159 | +0.0140 | +0.0019 | 0.0447 | no |
| MCF7 | `mse_capped_norm` | +0.0000 | +0.0000 | +0.0000 | 0.0000 | no |
| MCF7 | `lfc_nmae` | -5.2308 | -5.6192 | +0.3884 | 0.1128 | **si** |
| MCF7 | `fid` | +0.8895 | +0.8879 | **+0.0016** | 0.0013 | **si** |
| MCF7 | `reach` | -0.1577 | -0.1651 | +0.0074 | 0.0037 | **si** |
| MCF7 | `jac` | +0.3338 | +0.3386 | -0.0048 | 0.0032 | **si** |

### 4.1 `fid`: la diferencia en fraccion del valor

| linea | `fid` REAL | `fid` PERM | REAL - PERM | **PERM / REAL** |
|---|---|---|---|---|
| A549 | +1.0665 | +1.0672 | -0.0007 | **100.06 %** |
| BXPC3 | +0.6868 | +0.6857 | +0.0011 | **99.84 %** |
| HAP1 | +0.8066 | +0.8069 | -0.0003 | **100.03 %** |
| HT29 | +1.0090 | +1.0087 | +0.0004 | **99.96 %** |
| K562 | +0.8195 | +0.8129 | +0.0066 | **99.19 %** |
| MCF7 | +0.8895 | +0.8879 | +0.0016 | **99.82 %** |

Promediado sobre las 48 corridas de cada brazo: REAL `0.8797`, PERMUTADO
`0.8782`, es decir el permutado es el **99.83 %** del real.

**En `u` (sin normalizar), `fid` y `fid` sin los 6 replegados:**

| linea | `u` `fid` REAL | `u` `fid` PERM | dif | `u` sin los 6 REAL | `u` sin los 6 PERM | dif |
|---|---|---|---|---|---|---|
| A549 | 0.4457 | 0.4459 | -0.0003 | **0.4952** | **0.4955** | -0.0003 |
| BXPC3 | 0.4120 | 0.4113 | +0.0007 | **0.5036** | **0.5028** | +0.0008 |
| HAP1 | 0.4411 | 0.4412 | -0.0001 | **0.5065** | **0.5066** | -0.0002 |
| HT29 | 0.4347 | 0.4345 | +0.0001 | **0.5152** | **0.5150** | +0.0001 |
| K562 | 0.4332 | 0.4303 | +0.0029 | **0.5134** | **0.5100** | +0.0034 |
| MCF7 | 0.4115 | 0.4108 | +0.0008 | **0.5030** | **0.5020** | +0.0009 |

Sobre los 27 blancos con perfil, `u` de `fid` esta entre **0.4952 y 0.5152** en
el brazo REAL y entre **0.5020 y 0.5066** en el PERMUTADO.

---

## 5. 🔴 Los conteos que contestan la pregunta

| metrica | **REAL > PERM mas alla de `resolucion_8`** | REAL `s>0` | **PERM `s>0`** | REAL `s>`margen | **PERM `s>`margen** |
|---|---|---|---|---|---|
| `pds_cosine` | **2/6** | 4/6 | 3/6 | 0/6 | **1/6** |
| `mse_capped_norm` | **0/6** | 0/6 | 0/6 | 0/6 | **0/6** |
| `lfc_nmae` | **1/6** | 0/6 | 0/6 | 0/6 | **0/6** |
| **`fid`** | **3/6** | **6/6** | **6/6** | **6/6** | **6/6** |
| `reach` | **4/6** | 0/6 | 0/6 | 0/6 | **0/6** |
| `jac` | **4/6** | 2/6 | 2/6 | 2/6 | **2/6** |

El margen es el absoluto de tres terminos de `dispersion_b.md` §4
(`margen COMPLETO`), celda por celda.

**En `fid`: el brazo PERMUTADO da `s > 0` en 6/6 y `s >` margen en 6/6**, los
mismos conteos que el brazo REAL. Las 3 celdas donde `REAL > PERM` supera la
`resolucion_8` son BXPC3 (+0.0011), K562 (+0.0066) y MCF7 (+0.0016), sobre
valores de `fid` de +0.69 a +0.89.

---

## 6. 🔴 La tabla 2x2 de Venn por brazo y por linea

Media sobre las 8 semillas de cada brazo. "lleno" = conjunto de DE no vacio.
El `jaccard genuino medio` es `interseccion / union` promediado **solo sobre las
perturbaciones de la celda `ambos llenos`** -- sin la convencion que devuelve 1.0
cuando la union es vacia.

| linea | brazo | ambos vacios | real ∅ / pred lleno (FP) | real lleno / pred ∅ | ambos llenos | **jaccard genuino medio** |
|---|---|---|---|---|---|---|
| A549 | SUELO | 14.88 | 0.12 | 18.00 | 0.00 | n/d (sin `por_pert`) |
| A549 | **REAL** | 3.00 | **12.00** | 3.00 | 15.00 | **0.0150** |
| A549 | **PERMUTADO** | 3.00 | **12.00** | 3.00 | 15.00 | **0.0145** |
| BXPC3 | SUELO | 0.00 | 0.00 | 31.88 | 1.12 | n/d (sin `por_pert`) |
| BXPC3 | **REAL** | 0.00 | 0.00 | 6.00 | 27.00 | **0.0146** |
| BXPC3 | **PERMUTADO** | 0.00 | 0.00 | 6.00 | 27.00 | **0.0141** |
| HAP1 | SUELO | 13.62 | 0.38 | 18.75 | 0.25 | n/d (sin `por_pert`) |
| HAP1 | **REAL** | 2.00 | **12.00** | 3.75 | 15.25 | **0.0045** |
| HAP1 | **PERMUTADO** | 2.00 | **12.00** | 4.00 | 15.00 | **0.0048** |
| HT29 | SUELO | 5.88 | 0.12 | 26.62 | 0.38 | n/d (sin `por_pert`) |
| HT29 | **REAL** | 1.00 | 5.00 | 5.00 | 22.00 | **0.0106** |
| HT29 | **PERMUTADO** | 1.00 | 5.00 | 5.00 | 22.00 | **0.0101** |
| K562 | SUELO | 7.88 | 0.12 | 25.00 | 0.00 | n/d (sin `por_pert`) |
| K562 | **REAL** | 1.00 | 7.00 | 5.00 | 20.00 | **0.0050** |
| K562 | **PERMUTADO** | 1.00 | 7.00 | 4.88 | 20.12 | **0.0048** |
| MCF7 | SUELO | 0.00 | 0.00 | 32.62 | 0.38 | n/d (sin `por_pert`) |
| MCF7 | **REAL** | 0.00 | 0.00 | 6.00 | 27.00 | **0.0295** |
| MCF7 | **PERMUTADO** | 0.00 | 0.00 | 5.88 | 27.12 | **0.0298** |

**Como se obtuvo cada fila:**

* **REAL y PERMUTADO**: contados directamente de `cohortes.por_pert`
  (`n_DE_real`, `n_DE_pred`, `n_DE_interseccion`, `n_DE_union`) de cada corrida.
* 🔴 **SUELO**: `dispersion_u.json` **no guarda `por_pert`**, asi que las cuatro
  celdas se **derivan exactamente** de tres cantidades que si guarda --
  `R` (# real vacio, del lado real fijo), `P = n_perts_con_pred_sig` y
  `E = E_union_vacia`:

  ```
  ambos vacios      = E
  real ∅ / pred +   = R - E
  real + / pred ∅   = 33 - P - E
  ambos llenos      = P - (R - E)
  ```

  **Comprobado en las 48 corridas del suelo**: las cuatro celdas suman 33 y
  ninguna sale negativa, en 48 de 48. El `jaccard genuino` del suelo **no se
  puede derivar** de esas tres cantidades y queda **n/d, no medido**; su celda
  `ambos llenos` promedia 0.00 a 1.12 perturbaciones.

**Dos lecturas de la tabla que son aritmetica, no interpretacion:** en los dos
brazos de gwps la columna `real lleno / pred ∅` mas la de `ambos vacios` suma
**6.00** en A549, BXPC3, HT29 y MCF7 -- exactamente los 6 blancos replegados, que
son los unicos que no predicen nada. Y la celda `real ∅ / pred lleno` (falsos
positivos) vale **12.00** en A549 y HAP1, **7.00** en K562, **5.00** en HT29 y
**0.00** en BXPC3 y MCF7, **identica entre los dos brazos**.

---

## 7. 🔴 Perturbaciones con DE REAL vacio, MEDIDO en las seis lineas

Contado de `n_DE_real == 0` en `cohortes.por_pert`. El lado real depende solo de
la poblacion real, que esta fijada por `seed_real`, asi que es el mismo en los
dos brazos de gwps y en las 8 semillas de cada uno.

| linea | **# perturbaciones con DE real vacio** | de 33 | identico en las 8 semillas | `E_union_vacia` del SUELO (**otra cantidad**) |
|---|---|---|---|---|
| A549 | **15** | 33 | si | [14, 15] |
| BXPC3 | **0** | 33 | si | [0] |
| HAP1 | **14** | 33 | si | [13, 14] |
| HT29 | **6** | 33 | si | [5, 6] |
| K562 | **8** | 33 | si | [7, 8] |
| MCF7 | **0** | 33 | si | [0] |

🔴 **Las dos ultimas columnas son cantidades distintas y no hay que confundirlas.**
`# DE real vacio` es del lado real solo. `E_union_vacia` del suelo es
`# (real vacio Y pred vacio)`, que es `<=` la anterior y **varia con la semilla
del suelo**, porque en 18 de sus 48 corridas el suelo si llama algun gen
significativo. En K562 eso da `E_union_vacia` de 7 u 8 segun la semilla, mientras
que **el # de DE real vacio medido es 8**.

---

## 8. El rango `[min, max]` por metrica y por brazo

| metrica | SUELO | REAL | **PERMUTADO** |
|---|---|---|---|
| `pds_cosine` | [-0.0345, +0.2838] | [-0.0412, +0.1281] | **[-0.0706, +0.1849]** |
| `mse_capped_norm` | [+0.0000, +0.0000] | [+0.0000, +0.0000] | **[+0.0000, +0.0000]** |
| `lfc_nmae` | [-0.8778, -0.2079] | [-5.2308, -0.7978] | **[-5.6192, -0.8167]** |
| `fid` | [-0.3463, +0.0154] | [+0.6868, +1.0665] | **[+0.6857, +1.0672]** |
| `reach` | [-0.5432, -0.1281] | [-0.5368, -0.1577] | **[-0.6934, -0.1414]** |
| `jac` | [-0.2588, +0.4803] | [-19.0726, +0.3338] | **[-19.0853, +0.3386]** |
| `avg_score` | [-0.1945, -0.1112] | [-3.2630, -0.2696] | **[-3.2801, -0.2586]** |

### 8.1 La igualacion de potencia, medida

| linea | REAL `DEpred` | PERM `DEpred` | dif | dif relativa | REAL `n_fid` | PERM `n_fid` |
|---|---|---|---|---|---|---|
| A549 | 66 181 | 66 164 | -17 | -0.03 % | 30.00 | 30.00 |
| BXPC3 | 86 139 | 86 051 | -89 | -0.10 % | 33.00 | 33.00 |
| HAP1 | 98 172 | 98 127 | -45 | -0.05 % | 31.00 | 31.00 |
| HT29 | 83 640 | 83 749 | +110 | +0.13 % | 32.00 | 32.00 |
| K562 | 85 932 | 85 899 | -33 | -0.04 % | 32.00 | 32.00 |
| MCF7 | 96 536 | 96 194 | -342 | -0.35 % | 33.00 | 33.00 |

---

## 🔴 9. Los intermedios crudos: las 48 corridas una por una

Todas las cifras de §2 a §8 son **medias sobre 8**. Aqui estan las filas crudas
que hay detras de cada una, sin agregar, para que cada numero cargado del reporte
se pueda recomputar a mano. Las medias y `sd` al pie de cada bloque son las
mismas que aparecen en §2, §4, §6, §7 y §8.1.

### TABLA A: las 48 corridas crudas del brazo PERMUTADO

| linea | `j` | semilla desarreglo | semilla emision | `pds_cosine` | `mse_capped_norm` | `lfc_nmae` | `fid` | `reach` | `jac` | `avg_score` |
|---|---|---|---|---|---|---|---|---|---|---|
| A549 | 0 | 20272500 | 20272000 | +0.0163 | +0.0000 | -1.5467 | **+1.0668** | +0.0429 | -19.0900 | -3.2518 |
| A549 | 1 | 20272501 | 20272001 | -0.0455 | +0.0000 | -1.1392 | **+1.0661** | -0.4000 | -19.0933 | -3.2687 |
| A549 | 2 | 20272502 | 20272002 | +0.0683 | +0.0000 | -1.7991 | **+1.0658** | -0.0856 | -19.0787 | -3.3049 |
| A549 | 3 | 20272503 | 20272003 | -0.1529 | +0.0000 | -1.1802 | **+1.0720** | -0.1755 | -19.0608 | -3.2496 |
| A549 | 4 | 20272504 | 20272004 | -0.0553 | +0.0000 | -1.3077 | **+1.0634** | -0.0703 | -19.1056 | -3.2459 |
| A549 | 5 | 20272505 | 20272005 | -0.0390 | +0.0000 | -1.4606 | **+1.0668** | +0.0090 | -19.0756 | -3.2499 |
| A549 | 6 | 20272506 | 20272006 | +0.0163 | +0.0000 | -1.4125 | **+1.0688** | -0.1673 | -19.0926 | -3.2646 |
| A549 | 7 | 20272507 | 20272007 | -0.2082 | +0.0000 | -1.9240 | **+1.0680** | -0.2842 | -19.0861 | -3.4057 |
| **A549 media** | | | | -0.0500 | +0.0000 | -1.4712 | **+1.0672** | -0.1414 | -19.0853 | -3.2801 |
| **A549 sd** | | | | 0.0915 | 0.0000 | 0.2786 | **0.0025** | 0.1482 | 0.0136 | 0.0542 |
| BXPC3 | 0 | 20272500 | 20272100 | +0.0593 | +0.0000 | -1.9200 | **+0.6839** | -0.5619 | +0.0523 | -0.2811 |
| BXPC3 | 1 | 20272501 | 20272101 | +0.0218 | +0.0000 | -1.8701 | **+0.6851** | -0.4393 | +0.0539 | -0.2581 |
| BXPC3 | 2 | 20272502 | 20272102 | -0.0237 | +0.0000 | -1.7811 | **+0.6834** | -0.4973 | +0.0517 | -0.2612 |
| BXPC3 | 3 | 20272503 | 20272103 | +0.0613 | +0.0000 | -1.9482 | **+0.6878** | -0.6409 | +0.0552 | -0.2975 |
| BXPC3 | 4 | 20272504 | 20272104 | -0.0752 | +0.0000 | -1.8727 | **+0.6838** | -0.5421 | +0.0522 | -0.2923 |
| BXPC3 | 5 | 20272505 | 20272105 | -0.0396 | +0.0000 | -2.1302 | **+0.6885** | -0.5482 | +0.0504 | -0.3298 |
| BXPC3 | 6 | 20272506 | 20272106 | +0.1998 | +0.0000 | -2.0039 | **+0.6881** | -0.5514 | +0.0553 | -0.2687 |
| BXPC3 | 7 | 20272507 | 20272107 | +0.0455 | +0.0000 | -1.7508 | **+0.6849** | -0.5292 | +0.0519 | -0.2496 |
| **BXPC3 media** | | | | +0.0312 | +0.0000 | -1.9096 | **+0.6857** | -0.5388 | +0.0529 | -0.2798 |
| **BXPC3 sd** | | | | 0.0843 | 0.0000 | 0.1218 | **0.0021** | 0.0572 | 0.0017 | 0.0263 |
| HAP1 | 0 | 20272500 | 20272200 | +0.1022 | +0.0000 | -0.7690 | **+0.8070** | -0.7890 | -8.2362 | -1.4808 |
| HAP1 | 1 | 20272501 | 20272201 | +0.3147 | +0.0000 | -0.7717 | **+0.8111** | -0.3515 | -8.2343 | -1.3719 |
| HAP1 | 2 | 20272502 | 20272202 | +0.2959 | +0.0000 | -0.7894 | **+0.8085** | -0.6959 | -8.2317 | -1.4354 |
| HAP1 | 3 | 20272503 | 20272203 | +0.1291 | +0.0000 | -0.6548 | **+0.8068** | -0.6653 | -8.2266 | -1.4351 |
| HAP1 | 4 | 20272504 | 20272204 | +0.2985 | +0.0000 | -0.7361 | **+0.8044** | -0.3890 | -8.2315 | -1.3756 |
| HAP1 | 5 | 20272505 | 20272205 | -0.1022 | +0.0000 | -1.0075 | **+0.8013** | -0.3663 | -8.2350 | -1.4850 |
| HAP1 | 6 | 20272506 | 20272206 | +0.2985 | +0.0000 | -0.6537 | **+0.8091** | -0.4597 | -8.2239 | -1.3716 |
| HAP1 | 7 | 20272507 | 20272207 | +0.1425 | +0.0000 | -1.1511 | **+0.8069** | -0.3571 | -8.2325 | -1.4652 |
| **HAP1 media** | | | | +0.1849 | +0.0000 | -0.8167 | **+0.8069** | -0.5092 | -8.2315 | -1.4276 |
| **HAP1 sd** | | | | 0.1460 | 0.0000 | 0.1742 | **0.0030** | 0.1784 | 0.0042 | 0.0487 |
| HT29 | 0 | 20272500 | 20272300 | +0.0706 | +0.0000 | -1.2612 | **+1.0146** | -0.6832 | -0.3002 | -0.1932 |
| HT29 | 1 | 20272501 | 20272301 | +0.1042 | +0.0000 | -1.4544 | **+1.0062** | -0.5904 | -0.2983 | -0.2054 |
| HT29 | 2 | 20272502 | 20272302 | -0.0773 | +0.0000 | -1.7939 | **+1.0084** | -0.3825 | -0.2992 | -0.2574 |
| HT29 | 3 | 20272503 | 20272303 | +0.1143 | +0.0000 | -1.3640 | **+1.0149** | -0.6886 | -0.2973 | -0.2034 |
| HT29 | 4 | 20272504 | 20272304 | +0.0235 | +0.0000 | -1.7100 | **+1.0047** | -0.9413 | -0.2954 | -0.3197 |
| HT29 | 5 | 20272505 | 20272305 | -0.2320 | +0.0000 | -1.5551 | **+1.0049** | -0.6178 | -0.3018 | -0.2836 |
| HT29 | 6 | 20272506 | 20272306 | -0.2488 | +0.0000 | -1.5508 | **+1.0115** | -0.7302 | -0.2956 | -0.3023 |
| HT29 | 7 | 20272507 | 20272307 | +0.0706 | +0.0000 | -1.6813 | **+1.0042** | -0.9130 | -0.3014 | -0.3035 |
| **HT29 media** | | | | -0.0219 | +0.0000 | -1.5463 | **+1.0087** | -0.6934 | -0.2986 | -0.2586 |
| **HT29 sd** | | | | 0.1475 | 0.0000 | 0.1810 | **0.0044** | 0.1790 | 0.0025 | 0.0513 |
| K562 | 0 | 20272500 | 20272400 | -0.1474 | +0.0000 | -1.5534 | **+0.8087** | -0.5752 | -0.6439 | -0.3519 |
| K562 | 1 | 20272501 | 20272401 | -0.0878 | +0.0000 | -2.2928 | **+0.8099** | -0.9750 | -0.6433 | -0.5315 |
| K562 | 2 | 20272502 | 20272402 | -0.3074 | +0.0000 | -2.1449 | **+0.8130** | -0.5610 | -0.6433 | -0.4739 |
| K562 | 3 | 20272503 | 20272403 | +0.0220 | +0.0000 | -2.0537 | **+0.8131** | -0.5946 | -0.6437 | -0.4095 |
| K562 | 4 | 20272504 | 20272404 | -0.2666 | +0.0000 | -2.5460 | **+0.8100** | -0.6517 | -0.6434 | -0.5496 |
| K562 | 5 | 20272505 | 20272405 | +0.0345 | +0.0000 | -1.4598 | **+0.8142** | -0.6473 | -0.6440 | -0.3171 |
| K562 | 6 | 20272506 | 20272406 | +0.0565 | +0.0000 | -1.6145 | **+0.8193** | -0.4007 | -0.6437 | -0.2972 |
| K562 | 7 | 20272507 | 20272407 | +0.1317 | +0.0000 | -2.3328 | **+0.8147** | -0.2638 | -0.6445 | -0.3824 |
| **K562 media** | | | | -0.0706 | +0.0000 | -1.9997 | **+0.8129** | -0.5837 | -0.6437 | -0.4141 |
| **K562 sd** | | | | 0.1592 | 0.0000 | 0.4068 | **0.0034** | 0.2064 | 0.0004 | 0.0954 |
| MCF7 | 0 | 20272500 | 20272500 🔴 | +0.1073 | +0.0000 | -5.8025 | **+0.8868** | -0.1668 | +0.3383 | -0.7728 |
| MCF7 | 1 | 20272501 | 20272501 🔴 | +0.0810 | +0.0000 | -5.5447 | **+0.8826** | -0.1451 | +0.3470 | -0.7299 |
| MCF7 | 2 | 20272502 | 20272502 🔴 | +0.1183 | +0.0000 | -5.6214 | **+0.8884** | -0.1724 | +0.3526 | -0.7391 |
| MCF7 | 3 | 20272503 | 20272503 🔴 | -0.0810 | +0.0000 | -5.6518 | **+0.8884** | -0.1668 | +0.3320 | -0.7799 |
| MCF7 | 4 | 20272504 | 20272504 🔴 | +0.0394 | +0.0000 | -5.2124 | **+0.8922** | -0.1624 | +0.3388 | -0.6841 |
| MCF7 | 5 | 20272505 | 20272505 🔴 | +0.0788 | +0.0000 | -5.5206 | **+0.8898** | -0.1606 | +0.3244 | -0.7314 |
| MCF7 | 6 | 20272506 | 20272506 🔴 | -0.1927 | +0.0000 | -5.6177 | **+0.8835** | -0.1738 | +0.3404 | -0.7934 |
| MCF7 | 7 | 20272507 | 20272507 🔴 | -0.0394 | +0.0000 | -5.9824 | **+0.8914** | -0.1727 | +0.3353 | -0.8280 |
| **MCF7 media** | | | | +0.0140 | +0.0000 | -5.6192 | **+0.8879** | -0.1651 | +0.3386 | -0.7573 |
| **MCF7 sd** | | | | 0.1092 | 0.0000 | 0.2228 | **0.0035** | 0.0094 | 0.0087 | 0.0450 |

### TABLA B: los 96 `fid` crudos, los dos brazos

| `j` | linea | `fid` REAL | semilla emision REAL | `fid` PERM | semilla emision PERM | REAL - PERM (mismo indice `j`) |
|---|---|---|---|---|---|---|
| 0 | A549 | +1.0653 | 20271000 | +1.0668 | 20272000 | -0.0015 |
| 1 | A549 | +1.0648 | 20271001 | +1.0661 | 20272001 | -0.0013 |
| 2 | A549 | +1.0658 | 20271002 | +1.0658 | 20272002 | -0.0000 |
| 3 | A549 | +1.0665 | 20271003 | +1.0720 | 20272003 | -0.0055 |
| 4 | A549 | +1.0688 | 20271004 | +1.0634 | 20272004 | +0.0054 |
| 5 | A549 | +1.0664 | 20271005 | +1.0668 | 20272005 | -0.0005 |
| 6 | A549 | +1.0653 | 20271006 | +1.0688 | 20272006 | -0.0035 |
| 7 | A549 | +1.0696 | 20271007 | +1.0680 | 20272007 | +0.0016 |
| **media** | **A549** | **+1.0665** | | **+1.0672** | | **-0.0007** |
| **sd** | **A549** | **0.0017** | | **0.0025** | | |
| 0 | BXPC3 | +0.6878 | 20271100 | +0.6839 | 20272100 | +0.0039 |
| 1 | BXPC3 | +0.6858 | 20271101 | +0.6851 | 20272101 | +0.0007 |
| 2 | BXPC3 | +0.6883 | 20271102 | +0.6834 | 20272102 | +0.0049 |
| 3 | BXPC3 | +0.6873 | 20271103 | +0.6878 | 20272103 | -0.0004 |
| 4 | BXPC3 | +0.6856 | 20271104 | +0.6838 | 20272104 | +0.0018 |
| 5 | BXPC3 | +0.6876 | 20271105 | +0.6885 | 20272105 | -0.0009 |
| 6 | BXPC3 | +0.6863 | 20271106 | +0.6881 | 20272106 | -0.0019 |
| 7 | BXPC3 | +0.6857 | 20271107 | +0.6849 | 20272107 | +0.0008 |
| **media** | **BXPC3** | **+0.6868** | | **+0.6857** | | **+0.0011** |
| **sd** | **BXPC3** | **0.0011** | | **0.0021** | | |
| 0 | HAP1 | +0.8067 | 20271200 | +0.8070 | 20272200 | -0.0003 |
| 1 | HAP1 | +0.8067 | 20271201 | +0.8111 | 20272201 | -0.0044 |
| 2 | HAP1 | +0.8088 | 20271202 | +0.8085 | 20272202 | +0.0003 |
| 3 | HAP1 | +0.8047 | 20271203 | +0.8068 | 20272203 | -0.0020 |
| 4 | HAP1 | +0.8054 | 20271204 | +0.8044 | 20272204 | +0.0010 |
| 5 | HAP1 | +0.8064 | 20271205 | +0.8013 | 20272205 | +0.0051 |
| 6 | HAP1 | +0.8066 | 20271206 | +0.8091 | 20272206 | -0.0025 |
| 7 | HAP1 | +0.8077 | 20271207 | +0.8069 | 20272207 | +0.0008 |
| **media** | **HAP1** | **+0.8066** | | **+0.8069** | | **-0.0003** |
| **sd** | **HAP1** | **0.0013** | | **0.0030** | | |
| 0 | HT29 | +1.0097 | 20271300 | +1.0146 | 20272300 | -0.0049 |
| 1 | HT29 | +1.0080 | 20271301 | +1.0062 | 20272301 | +0.0018 |
| 2 | HT29 | +1.0117 | 20271302 | +1.0084 | 20272302 | +0.0033 |
| 3 | HT29 | +1.0080 | 20271303 | +1.0149 | 20272303 | -0.0069 |
| 4 | HT29 | +1.0081 | 20271304 | +1.0047 | 20272304 | +0.0034 |
| 5 | HT29 | +1.0091 | 20271305 | +1.0049 | 20272305 | +0.0042 |
| 6 | HT29 | +1.0076 | 20271306 | +1.0115 | 20272306 | -0.0039 |
| 7 | HT29 | +1.0100 | 20271307 | +1.0042 | 20272307 | +0.0058 |
| **media** | **HT29** | **+1.0090** | | **+1.0087** | | **+0.0004** |
| **sd** | **HT29** | **0.0014** | | **0.0044** | | |
| 0 | K562 | +0.8192 | 20271400 | +0.8087 | 20272400 | +0.0106 |
| 1 | K562 | +0.8193 | 20271401 | +0.8099 | 20272401 | +0.0094 |
| 2 | K562 | +0.8209 | 20271402 | +0.8130 | 20272402 | +0.0079 |
| 3 | K562 | +0.8202 | 20271403 | +0.8131 | 20272403 | +0.0071 |
| 4 | K562 | +0.8207 | 20271404 | +0.8100 | 20272404 | +0.0106 |
| 5 | K562 | +0.8178 | 20271405 | +0.8142 | 20272405 | +0.0036 |
| 6 | K562 | +0.8191 | 20271406 | +0.8193 | 20272406 | -0.0002 |
| 7 | K562 | +0.8187 | 20271407 | +0.8147 | 20272407 | +0.0040 |
| **media** | **K562** | **+0.8195** | | **+0.8129** | | **+0.0066** |
| **sd** | **K562** | **0.0010** | | **0.0034** | | |
| 0 | MCF7 | +0.8898 | 20271500 | +0.8868 | 20272500 | +0.0030 |
| 1 | MCF7 | +0.8879 | 20271501 | +0.8826 | 20272501 | +0.0053 |
| 2 | MCF7 | +0.8886 | 20271502 | +0.8884 | 20272502 | +0.0002 |
| 3 | MCF7 | +0.8897 | 20271503 | +0.8884 | 20272503 | +0.0013 |
| 4 | MCF7 | +0.8909 | 20271504 | +0.8922 | 20272504 | -0.0013 |
| 5 | MCF7 | +0.8883 | 20271505 | +0.8898 | 20272505 | -0.0015 |
| 6 | MCF7 | +0.8891 | 20271506 | +0.8835 | 20272506 | +0.0055 |
| 7 | MCF7 | +0.8919 | 20271507 | +0.8914 | 20272507 | +0.0004 |
| **media** | **MCF7** | **+0.8895** | | **+0.8879** | | **+0.0016** |
| **sd** | **MCF7** | **0.0014** | | **0.0035** | | |

🔴 **La TABLA B no son pares emparejados.** Los dos brazos usan familias de
semilla de remuestreo distintas (`20271000 + 100*fold + j` el REAL,
`20272000 + 100*fold + j` el PERMUTADO), asi que el `j` de una fila indexa
**sorteos de remuestreo diferentes** en los dos brazos. La columna
`REAL - PERM (mismo indice j)` es una resta por indice, **no** una diferencia
pareada, y no se usa en ningun conteo: los conteos de §5 y las diferencias de §4
salen de las **medias** de los 8, como declara §4.

### TABLA C: las cuatro celdas de Venn CONTADAS corrida por corrida, brazo PERMUTADO

| linea | `j` | ambos vacios | real ∅ / pred lleno | real lleno / pred ∅ | ambos llenos | suma | `jaccard` genuino | # DE real vacio | `DEpred` |
|---|---|---|---|---|---|---|---|---|---|
| A549 | 0 | 3 | 12 | 3 | 15 | 33 | 0.0143 | 15 | 66201 |
| A549 | 1 | 3 | 12 | 3 | 15 | 33 | 0.0141 | 15 | 66050 |
| A549 | 2 | 3 | 12 | 3 | 15 | 33 | 0.0147 | 15 | 65692 |
| A549 | 3 | 3 | 12 | 3 | 15 | 33 | 0.0155 | 15 | 65915 |
| A549 | 4 | 3 | 12 | 3 | 15 | 33 | 0.0136 | 15 | 66043 |
| A549 | 5 | 3 | 12 | 3 | 15 | 33 | 0.0149 | 15 | 66685 |
| A549 | 6 | 3 | 12 | 3 | 15 | 33 | 0.0142 | 15 | 66218 |
| A549 | 7 | 3 | 12 | 3 | 15 | 33 | 0.0144 | 15 | 66505 |
| **A549 media** | | **3.00** | **12.00** | **3.00** | **15.00** | 33.00 | **0.0145** | **15.00** | **66163.6** |
| BXPC3 | 0 | 0 | 0 | 6 | 27 | 33 | 0.0140 | 0 | 86088 |
| BXPC3 | 1 | 0 | 0 | 6 | 27 | 33 | 0.0144 | 0 | 85958 |
| BXPC3 | 2 | 0 | 0 | 6 | 27 | 33 | 0.0138 | 0 | 86044 |
| BXPC3 | 3 | 0 | 0 | 6 | 27 | 33 | 0.0147 | 0 | 86251 |
| BXPC3 | 4 | 0 | 0 | 6 | 27 | 33 | 0.0139 | 0 | 85813 |
| BXPC3 | 5 | 0 | 0 | 6 | 27 | 33 | 0.0135 | 0 | 86033 |
| BXPC3 | 6 | 0 | 0 | 6 | 27 | 33 | 0.0147 | 0 | 85879 |
| BXPC3 | 7 | 0 | 0 | 6 | 27 | 33 | 0.0139 | 0 | 86339 |
| **BXPC3 media** | | **0.00** | **0.00** | **6.00** | **27.00** | 33.00 | **0.0141** | **0.00** | **86050.6** |
| HAP1 | 0 | 2 | 12 | 4 | 15 | 33 | 0.0043 | 14 | 98369 |
| HAP1 | 1 | 2 | 12 | 4 | 15 | 33 | 0.0045 | 14 | 98012 |
| HAP1 | 2 | 2 | 12 | 4 | 15 | 33 | 0.0047 | 14 | 98422 |
| HAP1 | 3 | 2 | 12 | 4 | 15 | 33 | 0.0052 | 14 | 97434 |
| HAP1 | 4 | 2 | 12 | 4 | 15 | 33 | 0.0048 | 14 | 98213 |
| HAP1 | 5 | 2 | 12 | 4 | 15 | 33 | 0.0044 | 14 | 98247 |
| HAP1 | 6 | 2 | 12 | 4 | 15 | 33 | 0.0055 | 14 | 97841 |
| HAP1 | 7 | 2 | 12 | 4 | 15 | 33 | 0.0047 | 14 | 98481 |
| **HAP1 media** | | **2.00** | **12.00** | **4.00** | **15.00** | 33.00 | **0.0048** | **14.00** | **98127.4** |
| HT29 | 0 | 1 | 5 | 5 | 22 | 33 | 0.0097 | 6 | 83439 |
| HT29 | 1 | 1 | 5 | 5 | 22 | 33 | 0.0102 | 6 | 83843 |
| HT29 | 2 | 1 | 5 | 5 | 22 | 33 | 0.0099 | 6 | 83784 |
| HT29 | 3 | 1 | 5 | 5 | 22 | 33 | 0.0105 | 6 | 83714 |
| HT29 | 4 | 1 | 5 | 5 | 22 | 33 | 0.0110 | 6 | 84213 |
| HT29 | 5 | 1 | 5 | 5 | 22 | 33 | 0.0092 | 6 | 83955 |
| HT29 | 6 | 1 | 5 | 5 | 22 | 33 | 0.0109 | 6 | 83550 |
| HT29 | 7 | 1 | 5 | 5 | 22 | 33 | 0.0093 | 6 | 83496 |
| **HT29 media** | | **1.00** | **5.00** | **5.00** | **22.00** | 33.00 | **0.0101** | **6.00** | **83749.2** |
| K562 | 0 | 1 | 7 | 5 | 20 | 33 | 0.0048 | 8 | 86052 |
| K562 | 1 | 1 | 7 | 5 | 20 | 33 | 0.0051 | 8 | 85648 |
| K562 | 2 | 1 | 7 | 5 | 20 | 33 | 0.0051 | 8 | 85830 |
| K562 | 3 | 1 | 7 | 4 | 21 | 33 | 0.0046 | 8 | 85945 |
| K562 | 4 | 1 | 7 | 5 | 20 | 33 | 0.0050 | 8 | 86054 |
| K562 | 5 | 1 | 7 | 5 | 20 | 33 | 0.0047 | 8 | 85949 |
| K562 | 6 | 1 | 7 | 5 | 20 | 33 | 0.0049 | 8 | 85761 |
| K562 | 7 | 1 | 7 | 5 | 20 | 33 | 0.0045 | 8 | 85952 |
| **K562 media** | | **1.00** | **7.00** | **4.88** | **20.12** | 33.00 | **0.0048** | **8.00** | **85898.9** |
| MCF7 | 0 | 0 | 0 | 6 | 27 | 33 | 0.0299 | 0 | 96459 |
| MCF7 | 1 | 0 | 0 | 6 | 27 | 33 | 0.0307 | 0 | 95692 |
| MCF7 | 2 | 0 | 0 | 6 | 27 | 33 | 0.0312 | 0 | 96068 |
| MCF7 | 3 | 0 | 0 | 6 | 27 | 33 | 0.0293 | 0 | 96158 |
| MCF7 | 4 | 0 | 0 | 5 | 28 | 33 | 0.0289 | 0 | 95908 |
| MCF7 | 5 | 0 | 0 | 6 | 27 | 33 | 0.0287 | 0 | 96190 |
| MCF7 | 6 | 0 | 0 | 6 | 27 | 33 | 0.0301 | 0 | 96651 |
| MCF7 | 7 | 0 | 0 | 6 | 27 | 33 | 0.0296 | 0 | 96426 |
| **MCF7 media** | | **0.00** | **0.00** | **5.88** | **27.12** | 33.00 | **0.0298** | **0.00** | **96194.0** |

### TABLA D: las mismas, brazo REAL

| linea | `j` | ambos vacios | real ∅ / pred lleno | real lleno / pred ∅ | ambos llenos | suma | `jaccard` genuino | # DE real vacio | `DEpred` |
|---|---|---|---|---|---|---|---|---|---|
| A549 | 0 | 3 | 12 | 3 | 15 | 33 | 0.0149 | 15 | 66224 |
| A549 | 1 | 3 | 12 | 3 | 15 | 33 | 0.0151 | 15 | 65984 |
| A549 | 2 | 3 | 12 | 3 | 15 | 33 | 0.0153 | 15 | 65800 |
| A549 | 3 | 3 | 12 | 3 | 15 | 33 | 0.0148 | 15 | 66385 |
| A549 | 4 | 3 | 12 | 3 | 15 | 33 | 0.0148 | 15 | 66103 |
| A549 | 5 | 3 | 12 | 3 | 15 | 33 | 0.0152 | 15 | 66013 |
| A549 | 6 | 3 | 12 | 3 | 15 | 33 | 0.0149 | 15 | 66572 |
| A549 | 7 | 3 | 12 | 3 | 15 | 33 | 0.0149 | 15 | 66364 |
| **A549 media** | | **3.00** | **12.00** | **3.00** | **15.00** | 33.00 | **0.0150** | **15.00** | **66180.6** |
| BXPC3 | 0 | 0 | 0 | 6 | 27 | 33 | 0.0147 | 0 | 85794 |
| BXPC3 | 1 | 0 | 0 | 6 | 27 | 33 | 0.0147 | 0 | 86205 |
| BXPC3 | 2 | 0 | 0 | 6 | 27 | 33 | 0.0143 | 0 | 86445 |
| BXPC3 | 3 | 0 | 0 | 6 | 27 | 33 | 0.0147 | 0 | 86515 |
| BXPC3 | 4 | 0 | 0 | 6 | 27 | 33 | 0.0144 | 0 | 85935 |
| BXPC3 | 5 | 0 | 0 | 6 | 27 | 33 | 0.0146 | 0 | 86048 |
| BXPC3 | 6 | 0 | 0 | 6 | 27 | 33 | 0.0147 | 0 | 85913 |
| BXPC3 | 7 | 0 | 0 | 6 | 27 | 33 | 0.0146 | 0 | 86259 |
| **BXPC3 media** | | **0.00** | **0.00** | **6.00** | **27.00** | 33.00 | **0.0146** | **0.00** | **86139.2** |
| HAP1 | 0 | 2 | 12 | 4 | 15 | 33 | 0.0046 | 14 | 98151 |
| HAP1 | 1 | 2 | 12 | 4 | 15 | 33 | 0.0047 | 14 | 98198 |
| HAP1 | 2 | 2 | 12 | 4 | 15 | 33 | 0.0045 | 14 | 98313 |
| HAP1 | 3 | 2 | 12 | 3 | 16 | 33 | 0.0045 | 14 | 98140 |
| HAP1 | 4 | 2 | 12 | 4 | 15 | 33 | 0.0046 | 14 | 97889 |
| HAP1 | 5 | 2 | 12 | 4 | 15 | 33 | 0.0046 | 14 | 98434 |
| HAP1 | 6 | 2 | 12 | 3 | 16 | 33 | 0.0044 | 14 | 97985 |
| HAP1 | 7 | 2 | 12 | 4 | 15 | 33 | 0.0045 | 14 | 98270 |
| **HAP1 media** | | **2.00** | **12.00** | **3.75** | **15.25** | 33.00 | **0.0045** | **14.00** | **98172.5** |
| HT29 | 0 | 1 | 5 | 5 | 22 | 33 | 0.0108 | 6 | 83397 |
| HT29 | 1 | 1 | 5 | 5 | 22 | 33 | 0.0108 | 6 | 83992 |
| HT29 | 2 | 1 | 5 | 5 | 22 | 33 | 0.0108 | 6 | 83670 |
| HT29 | 3 | 1 | 5 | 5 | 22 | 33 | 0.0105 | 6 | 83447 |
| HT29 | 4 | 1 | 5 | 5 | 22 | 33 | 0.0107 | 6 | 83961 |
| HT29 | 5 | 1 | 5 | 5 | 22 | 33 | 0.0103 | 6 | 83575 |
| HT29 | 6 | 1 | 5 | 5 | 22 | 33 | 0.0105 | 6 | 83521 |
| HT29 | 7 | 1 | 5 | 5 | 22 | 33 | 0.0105 | 6 | 83554 |
| **HT29 media** | | **1.00** | **5.00** | **5.00** | **22.00** | 33.00 | **0.0106** | **6.00** | **83639.6** |
| K562 | 0 | 1 | 7 | 5 | 20 | 33 | 0.0048 | 8 | 86239 |
| K562 | 1 | 1 | 7 | 5 | 20 | 33 | 0.0051 | 8 | 85720 |
| K562 | 2 | 1 | 7 | 5 | 20 | 33 | 0.0050 | 8 | 85904 |
| K562 | 3 | 1 | 7 | 5 | 20 | 33 | 0.0049 | 8 | 85865 |
| K562 | 4 | 1 | 7 | 5 | 20 | 33 | 0.0051 | 8 | 85915 |
| K562 | 5 | 1 | 7 | 5 | 20 | 33 | 0.0050 | 8 | 86056 |
| K562 | 6 | 1 | 7 | 5 | 20 | 33 | 0.0051 | 8 | 85730 |
| K562 | 7 | 1 | 7 | 5 | 20 | 33 | 0.0049 | 8 | 86026 |
| **K562 media** | | **1.00** | **7.00** | **5.00** | **20.00** | 33.00 | **0.0050** | **8.00** | **85931.9** |
| MCF7 | 0 | 0 | 0 | 6 | 27 | 33 | 0.0299 | 0 | 96681 |
| MCF7 | 1 | 0 | 0 | 6 | 27 | 33 | 0.0291 | 0 | 96559 |
| MCF7 | 2 | 0 | 0 | 6 | 27 | 33 | 0.0295 | 0 | 96805 |
| MCF7 | 3 | 0 | 0 | 6 | 27 | 33 | 0.0297 | 0 | 96182 |
| MCF7 | 4 | 0 | 0 | 6 | 27 | 33 | 0.0297 | 0 | 96561 |
| MCF7 | 5 | 0 | 0 | 6 | 27 | 33 | 0.0295 | 0 | 96346 |
| MCF7 | 6 | 0 | 0 | 6 | 27 | 33 | 0.0292 | 0 | 96400 |
| MCF7 | 7 | 0 | 0 | 6 | 27 | 33 | 0.0295 | 0 | 96757 |
| **MCF7 media** | | **0.00** | **0.00** | **6.00** | **27.00** | 33.00 | **0.0295** | **0.00** | **96536.4** |

**Lo que las TABLAS C y D muestran corrida por corrida:** las cuatro celdas suman
**33 en 96 de 96 corridas**. En A549, BXPC3 y HT29 las cuatro celdas son
**constantes en las 8 semillas de los dos brazos**; las unicas celdas que se
mueven son HAP1 REAL (`real lleno / pred ∅` 3 o 4), K562 PERMUTADO y MCF7
PERMUTADO (una perturbacion que cruza en 1 de las 8 semillas), que es de donde
salen los 3.75, 4.88, 15.25, 20.12, 5.88 y 27.12 de §6. La columna
**`# DE real vacio` es constante en las 8 semillas de los dos brazos en las seis
lineas**, que es la verificacion que §7 declara.

---

## 10. Defectos de la construccion, como datos sobre la medicion

1. **El desarreglo no controla la direccion compartida con el suelo, solo la
   identidad del blanco.** Los dos brazos emiten vectores `r` del mismo conjunto
   de 27; ninguno de los dos aisla lo que `fid` cobraria con una direccion
   generica construida desde el propio control.
2. **8 desarreglos no cubren el espacio.** Hay 27! permutaciones y se sortearon 8
   desarreglos. La `sd` sobre 8 no es la `sd` sobre todos los desarreglos
   posibles.
3. **Un desarreglo deja a cada blanco con el perfil de otro blanco DEL PANEL**,
   no de un gen cualquiera de gwps. Los 33 del panel son de vias de interferon,
   asi que el brazo permutado sigue emitiendo perfiles de la misma familia. **No
   se midio** la version que asigna perfiles de blancos de gwps fuera del panel.
4. **`jaccard genuino` del SUELO: n/d.** `dispersion_u.json` no guarda `por_pert`
   y la celda `ambos llenos` del suelo promedia 0.00-1.12 perturbaciones. No se
   re-corrio el suelo para obtenerlo.
5. **Las cuatro celdas del SUELO son derivadas, no contadas**, aunque la
   derivacion sea exacta y este comprobada en 48/48 (§6).
6. **`resolucion_8` en `fid` y `jac` es diminuta** (0.0002-0.0016) porque las
   `sd` de los dos brazos de gwps son diminutas; por eso diferencias de 0.0003 a
   0.0066 la superan. **El tamano de la diferencia y su resolucion son dos cosas
   distintas** y las dos estan en §4.
7. **No se cambio la escala de magnitud** en esta corrida, por construccion
   (§1.1). Lo que esa escala vale o no vale sigue sin medirse aqui.
8. 🔴 **COLISION DE SEMILLA EN MCF7, contada: 8 de 48 corridas** (§1.3). La
   familia de remuestreo `20272000 + 100*fold + j` da `20272500 + j` en
   `fold = 5`, que es la familia del desarreglo. En MCF7, y solo en MCF7, el
   remuestreo de cada `j` arranca del mismo estado de RNG que fijo `sigma_j`, asi
   que **no se sorteo con independencia de la familia de desarreglos**. Las otras
   cinco lineas: 0 colisiones. El brazo REAL: 0. `sigma_j` no se altera (identico
   en las seis lineas, §1.4). **No se re-corrio MCF7 con otra familia de
   remuestreo**, asi que no esta medido cuanto -- si algo -- mueve esa colision.
9. **No se midio** el brazo permutado con P2, ni con mas de 8 desarreglos, ni
   ningun brazo con el gen blanco dentro de la comparacion DE.

---

## 11. Lo que no se hizo, y las entradas

| prohibicion | estado |
|---|---|
| cambiar la escala de magnitud | **no**. Se llama `gwps_emision.construye_ratios` sin modificar y solo se reasigna que vector recibe cada etiqueta |
| usar cifras de `~/vcc/benchmark/` o `/datos/benchmark_VEDADO_*` | **no leidos** |
| ajustar algo mirando el `s` | **no**. Los 8 desarreglos salen de `20272500 + j` por rechazo de puntos fijos, fijados antes de puntuar; no se descarto ni se repitio ninguno |
| promediar P1 y P2 | **no**. **P1 en los dos brazos** |
| modificar el arnes | **no**. `suelo_arnes.py`, `cand01_covar._reescala_entero` y `gwps_emision.construye_ratios` se importan sin tocar |
| escribir veredicto o conclusiones | **no** |
| envios, `vcc sample`, commits | **ninguno** |

| | |
|---|---|
| entradas | `/datos/ref/replogle/k562_gwps_raw_bulk.h5ad`; `/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad`; `modelo/out/suelo_panel.json`; bundles `/datos/modelo/out/nativas_comp/BUNDLE_<linea>`; `modelo/out/gwps_via.json` y `modelo/out/dispersion_u.json` (brazos ya medidos); `modelo/out/dispersion_b.md` §4 (margen) |
| codigo | `modelo/src/gwps_permutado.py`, `modelo/src/gwps_perm_via.py`, `modelo/src/control_permutado_tablas.py` |
| salidas | `modelo/out/control_permutado.json` (48), `modelo/out/control_permutado_run.log` |
| semillas | desarreglo `20272500 + j`; emision `20272000 + 100*fold + j`; poblacion real `20260001 + fold` |
| verificaciones | `score_ok` 48/48, `config_hash` 48/48, `poblacion_real_identica` 48/48 (1 huella por linea), 48 huellas `pred` distintas y 0 compartidas con el brazo REAL, `todo_entero` 48/48, 0 negativos, desarreglos sin punto fijo 48/48, `sigma` identico entre lineas para cada `j` |
