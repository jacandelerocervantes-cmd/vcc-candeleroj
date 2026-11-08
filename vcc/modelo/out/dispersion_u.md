# La dispersion de `u`

Cuanto se mueve el `s` de una MISMA emision al cambiar **solo** la semilla de
emision. 8 semillas x 6 lineas = 48 corridas `run` -> `score --real-bundle`
contra los bundles de competencia de `/datos/modelo/out/nativas_comp/`.

La emision es `suelo_floor.predict_floor` sin modificar -- la misma desechable de
siempre, no un candidato. Las 8 corridas de una linea **no se promedian**: se
puntuan por separado y lo que se reporta es su dispersion.

## 1. Lo que vario y lo que no, con la prueba

**Vario solo la semilla de emision.** Familia nueva `20264000 + 100*fold + j`,
`j = 0..7`. Ninguna colision con `2026000x` (verdad), `2026100x` (suelo),
`2026200x`/`2026300x` (techo) ni la `0` del anchor.

**No vario la poblacion real.** `build_fold_adatas` sin modificar con
`seed_real = 20260001 + fold_idx`, etiqueta `non-targeting`, identica a
`cable_suelo.md`. No se asume: de cada construccion se tomo un `sha256` de la
matriz real (`data`, `indices`, `indptr`) y de sus etiquetas. Las 8 de cada linea
tienen que coincidir, y el script aborta si no.

| linea | corridas | semillas de emision | huellas distintas del REAL | huellas distintas del PRED | config_hash de competencia en todas | score corrio en todas |
|---|---|---|---|---|---|---|
| A549 | 8 | 20264000..20264007 | **1** | **8** | **si** | **si** |
| BXPC3 | 8 | 20264100..20264107 | **1** | **8** | **si** | **si** |
| HAP1 | 8 | 20264200..20264207 | **1** | **8** | **si** | **si** |
| HT29 | 8 | 20264300..20264307 | **1** | **8** | **si** | **si** |
| K562 | 8 | 20264400..20264407 | **1** | **8** | **si** | **si** |
| MCF7 | 8 | 20264500..20264507 | **1** | **8** | **si** | **si** |

**Una sola huella del real por linea y ocho del predicho**: la poblacion esta
fija y la emision cambia en las ocho. Es la condicion que hace que la dispersion
medida tenga una sola fuente.

El `config_hash` recalculado desde el `run_params.yaml` de cada corrida coincide
con el de referencia `a5519017556705ab...` en las 48. `--write-degenes` -- que se
pasa para poder MEDIR las dos cohortes moviles en vez de inferirlas -- no lo
desvia.

Ninguna corrida fallo: 48 de 48 `run` y 48 de 48 `score --real-bundle`.

## 2. El margen combinado: la formula, y por que esa

La escala es `s = (u - b)/(r - b)` (`score.py:395-399`). De las tres cantidades,
esta medicion mueve **una**:

| cantidad | que es | su ruido aqui |
|---|---|---|
| `u` | la emision del candidato | **lo que mide este encargo**: 8 semillas, `sd_u` |
| `r` | el ancla replicada, media de 5 particiones | `SE_r = sd_r/sqrt(5)`, de `anclas_nativas_v2.md` |
| `b` | el baseline generico, **una sola realizacion** con `seed 0` | no estimado: ver la reserva de abajo |

Propagacion de primer orden. Con `b` fijo y `u`, `r` independientes:

```
ds/du = 1/(r - b)
ds/dr = -(u - b)/(r - b)^2 = -s/(r - b)

sd(s) de la emision   = sd_u / |r - b|          <- ya viene medida directa en s
sd(s) del ancla       = |s| * SE_r / |r - b|

margen_celda = sqrt( sd_s^2 + ( |s| * SE_r / |r - b| )^2 )
```

Dos detalles que importan y que no se aplican en silencio:

1. **`sd_s` se mide sobre los `s`, no sobre los `u`.** Las 8 corridas se puntuaron
   por separado contra el mismo bundle, asi que la division por `(r - b)` ya esta
   dentro de cada `s`. Dividir otra vez seria contarla dos veces. El `sd_u` crudo
   se reporta igual, en su columna, porque es la cantidad que no depende del
   bundle.
2. **El termino del ancla escala con `|s|`.** En estas 48 corridas los `s` son ~0 o
   negativos y pequenos, asi que ese termino es pequeno **por el valor de `s`, no
   porque `r` sea estable**. Para un candidato con `s` grande pesaria mas.

**La reserva sobre `b`.** `b` es una medicion, no una constante: `cell-eval2
baseline --seed 0` sobre la poblacion evaluada. Su propia dispersion no se estimo
-- habria hecho falta correr el baseline con varias semillas, que no es este
encargo -- asi que el margen de arriba **omite un termino**. Queda declarado: el
margen es una cota INFERIOR del ruido total de `s`.

**Y el caso que no necesita el termino del ancla.** Para comparar DOS candidatos en
la MISMA celda del MISMO bundle, `r` y `b` son comunes y se cancelan:

```
s1 - s2 = (u1 - u2)/(r - b)
```

El ruido del ancla no entra en la diferencia: mueve las dos `s` juntas. Ahi la
cantidad relevante es solo `sd_s`. Las dos columnas estan en la tabla para que se
pueda usar la que corresponda; **no se propone ningun umbral**.

**Convencion de `sd`, declarada.** `sd_r` se toma tal como lo estampa
`anclas_nativas_v2.md`, que uso `ddof=0` (verificado: su 0.035008 para A549/jac se
reproduce con `ddof=0` y no con `ddof=1`, que da 0.039140). Para `u` se reporta
`ddof=1` como principal -- es el estimador no sesgado de la dispersion de UNA
emision, que es la pregunta -- y `ddof=0` al lado. Con 8 draws la razon entre los
dos es `sqrt(8/7) = 1.069`.

## 3. Las 36 celdas

Antes de desglosarlas, el agregado, porque es el numero que decide y sale gratis
de las mismas 48 corridas. `avg_score` es el promedio plano de las seis `s`
(`from_replicate`; `from_baseline` queda null con un bundle de competencia):

| linea | los 8 avg_score | media | sd (ddof=1) | rango | min | max |
|---|---|---|---|---|---|---|
| A549 | -0.0714 / -0.1404 / -0.1831 / -0.0909 / -0.1485 / -0.3190 / -0.0827 / -0.1528 | -0.148595 | 0.079230 | 0.247555 | -0.3190 | -0.0714 |
| BXPC3 | -0.1808 / -0.2156 / -0.1714 / -0.2240 / -0.1988 / -0.2121 / -0.1887 / -0.1646 | -0.194501 | 0.021674 | 0.059347 | -0.2240 | -0.1646 |
| HAP1 | -0.2070 / -0.0598 / -0.1484 / -0.0408 / -0.0921 / -0.1462 / -0.0367 / -0.1991 | -0.116274 | 0.068469 | 0.170262 | -0.2070 | -0.0367 |
| HT29 | -0.0841 / -0.1507 / -0.1058 / -0.1215 / -0.0964 / -0.0993 / -0.1512 / -0.0803 | -0.111165 | 0.027628 | 0.070872 | -0.1512 | -0.0803 |
| K562 | -0.1761 / -0.1942 / -0.1014 / -0.1373 / -0.1989 / -0.1809 / -0.1502 / -0.2231 | -0.170239 | 0.038850 | 0.121690 | -0.2231 | -0.1014 |
| MCF7 | -0.1319 / -0.1511 / -0.1422 / -0.1935 / -0.1297 / -0.1629 / -0.1919 / -0.1693 | -0.159052 | 0.024858 | 0.063783 | -0.1935 | -0.1297 |

**Solo cambiando la semilla de emision, el `avg_score` de una misma emision
se mueve en un rango de 0.0593 a 0.2476** segun la
linea, con `sd` de 0.0217 a 0.0792. Es la cifra que faltaba
para poder separar una diferencia entre candidatos del sorteo.

Y conviene decir de donde NO sale: las seis `s` que promedia no tienen la misma
dispersion -- hay celdas con `sd = 0` exacta por razones distintas (§5) -- asi
que todo el movimiento del agregado lo ponen las celdas que si se mueven,
repartido a 1/6 cada una.

⚠️ **El caso de A549 no es ruido continuo.** Su rango (0.2476) y su `sd` (0.0792)
son los mayores de las seis, y 0.270 de ese rango los pone una sola celda en una
sola semilla: `j=5` perdio el 1.0 de convencion de **STAT5A** en `jac` y cobro
`s = -1.6197` en vez de 0, que a 1/6 son -0.270 de agregado. Las otras siete
semillas caen entre -0.0714 y -0.1831. §5.2 lo desarma entero. Leer la `sd` de
A549 como la anchura de una nube seria leer mal la celda.

Los 8 valores de `s` de cada celda, con su dispersion. `s` es la columna
`from_replicate` de `scored.csv`, que con un bundle de competencia es la que toma
`avg_score`.

| linea | metrica | los 8 valores de s | media | sd (ddof=1) | sd (ddof=0) | rango | min | max |
|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | +0.4684 / -0.2993 / -0.1334 / +0.0455 / -0.1919 / +0.3578 / +0.2993 / -0.0878 | 0.057336 | 0.284040 | 0.265695 | 0.767729 | -0.2993 | 0.4684 |
| A549 | `mse_capped_norm` | +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.0000 | 0.0000 |
| A549 | `lfc_nmae` | -0.3675 / -0.4172 / -0.4066 / -0.3631 / -0.4576 / -0.3594 / -0.4458 / -0.3462 | -0.395435 | 0.042320 | 0.039586 | 0.111401 | -0.4576 | -0.3462 |
| A549 | `fid` | +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.0000 | 0.0000 |
| A549 | `reach` | -0.5295 / -0.1259 / -0.5586 / -0.2275 / -0.2414 / -0.2926 / -0.3497 / -0.4829 | -0.351003 | 0.157531 | 0.147357 | 0.432720 | -0.5586 | -0.1259 |
| A549 | `jac` | +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / -1.6197 / +0.0000 / +0.0000 | -0.202468 | 0.572666 | 0.535680 | 1.619745 | -1.6197 | 0.0000 |
| BXPC3 | `pds_cosine` | +0.0771 / +0.0257 / +0.2176 / +0.0020 / +0.0989 / +0.0653 / +0.1048 / +0.2591 | 0.106309 | 0.089199 | 0.083438 | 0.257120 | 0.0020 | 0.2591 |
| BXPC3 | `mse_capped_norm` | +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.0000 | 0.0000 |
| BXPC3 | `lfc_nmae` | -0.7029 / -0.7520 / -0.7296 / -0.7205 / -0.7403 / -0.7521 / -0.7157 / -0.6798 | -0.724108 | 0.024938 | 0.023327 | 0.072312 | -0.7521 | -0.6798 |
| BXPC3 | `fid` | +0.0060 / -0.0016 / -0.0024 / -0.0024 / -0.0023 / -0.0016 / -0.0024 / -0.0022 | -0.001123 | 0.002910 | 0.002722 | 0.008448 | -0.0024 | 0.0060 |
| BXPC3 | `reach` | -0.4602 / -0.5607 / -0.5093 / -0.6180 / -0.5444 / -0.5794 / -0.5138 / -0.5599 | -0.543194 | 0.048389 | 0.045263 | 0.157782 | -0.6180 | -0.4602 |
| BXPC3 | `jac` | -0.0049 / -0.0049 / -0.0049 / -0.0049 / -0.0049 / -0.0049 / -0.0049 / -0.0049 | -0.004893 | 0.000000 | 0.000000 | 0.000000 | -0.0049 | -0.0049 |
| HAP1 | `pds_cosine` | +0.1775 / +0.4196 / +0.1210 / +0.3228 / +0.0780 / +0.2421 / +0.4922 / +0.4169 | 0.283755 | 0.152589 | 0.142734 | 0.414201 | 0.0780 | 0.4922 |
| HAP1 | `mse_capped_norm` | +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.0000 | 0.0000 |
| HAP1 | `lfc_nmae` | -0.2501 / -0.2565 / -0.2601 / -0.1202 / -0.1454 / -0.2422 / -0.1293 / -0.2595 | -0.207934 | 0.063788 | 0.059669 | 0.139880 | -0.2601 | -0.1202 |
| HAP1 | `fid` | +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0914 / +0.0321 / +0.0000 | 0.015439 | 0.032694 | 0.030583 | 0.091433 | 0.0000 | 0.0914 |
| HAP1 | `reach` | -0.4791 / -0.5221 / -0.7514 / -0.4472 / -0.4854 / -0.2784 / -0.6152 / -0.6621 | -0.530135 | 0.145648 | 0.136241 | 0.472931 | -0.7514 | -0.2784 |
| HAP1 | `jac` | -0.6901 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / -0.6901 / +0.0000 / -0.6901 | -0.258771 | 0.357138 | 0.334072 | 0.690056 | -0.6901 | 0.0000 |
| HT29 | `pds_cosine` | +0.3228 / -0.1244 / +0.0941 / -0.1210 / +0.0471 / +0.1648 / -0.1379 / -0.1748 | 0.008826 | 0.177889 | 0.166400 | 0.497646 | -0.1748 | 0.3228 |
| HT29 | `mse_capped_norm` | +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.0000 | 0.0000 |
| HT29 | `lfc_nmae` | -0.5103 / -0.4754 / -0.3137 / -0.3660 / -0.4625 / -0.4051 / -0.3571 / -0.2880 | -0.397264 | 0.079919 | 0.074758 | 0.222275 | -0.5103 | -0.2880 |
| HT29 | `fid` | -0.3595 / -0.3595 / -0.3595 / -0.3595 / -0.3595 / -0.3506 / -0.2624 / -0.3595 | -0.346264 | 0.034050 | 0.031850 | 0.097178 | -0.3595 | -0.2624 |
| HT29 | `reach` | -0.4588 / -0.2785 / -0.5570 / -0.3837 / -0.3045 / -0.5060 / -0.6511 / -0.1608 | -0.412553 | 0.161456 | 0.151028 | 0.490304 | -0.6511 | -0.1608 |
| HT29 | `jac` | +0.5012 / +0.3338 / +0.5012 / +0.5012 / +0.5012 / +0.5012 / +0.5012 / +0.5012 | 0.480264 | 0.059183 | 0.055361 | 0.167396 | 0.3338 | 0.5012 |
| K562 | `pds_cosine` | +0.0376 / -0.1004 / +0.1757 / +0.1506 / -0.0816 / -0.0439 / -0.1945 / -0.2196 | -0.034504 | 0.146438 | 0.136980 | 0.395232 | -0.2196 | 0.1757 |
| K562 | `mse_capped_norm` | +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.0000 | 0.0000 |
| K562 | `lfc_nmae` | -0.4200 / -0.5418 / -0.4075 / -0.4960 / -0.4774 / -0.4824 / -0.4949 / -0.5227 | -0.480347 | 0.046277 | 0.043289 | 0.134279 | -0.5418 | -0.4075 |
| K562 | `fid` | -0.1769 / -0.1769 / -0.1769 / -0.0885 / -0.1769 / -0.1769 / -0.1769 / -0.1769 | -0.165865 | 0.031276 | 0.029256 | 0.088462 | -0.1769 | -0.0885 |
| K562 | `reach` | -0.6062 / -0.4549 / -0.3085 / -0.3897 / -0.5669 / -0.4909 / -0.1439 / -0.5282 | -0.436141 | 0.152003 | 0.142186 | 0.462362 | -0.6062 | -0.1439 |
| K562 | `jac` | +0.1091 / +0.1091 / +0.1091 / +0.0000 / +0.1091 / +0.1091 / +0.1091 / +0.1091 | 0.095421 | 0.038556 | 0.036066 | 0.109053 | 0.0000 | 0.1091 |
| MCF7 | `pds_cosine` | +0.1971 / +0.0504 / +0.1949 / -0.1248 / +0.1796 / +0.0088 / -0.0832 / -0.0044 | 0.052289 | 0.126699 | 0.118516 | 0.321945 | -0.1248 | 0.1971 |
| MCF7 | `mse_capped_norm` | +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 / +0.0000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.0000 | 0.0000 |
| MCF7 | `lfc_nmae` | -0.8712 / -0.8451 / -0.9058 / -0.9066 / -0.8484 / -0.8371 / -0.9399 / -0.8681 | -0.877766 | 0.036243 | 0.033902 | 0.102721 | -0.9399 | -0.8371 |
| MCF7 | `fid` | -0.0003 / -0.0004 / -0.0004 / -0.0004 / -0.0000 / -0.0002 / -0.0004 / -0.0004 | -0.000303 | 0.000126 | 0.000118 | 0.000350 | -0.0004 | -0.0000 |
| MCF7 | `reach` | -0.1168 / -0.1109 / -0.1418 / -0.1285 / -0.1088 / -0.1483 / -0.1274 / -0.1426 | -0.128132 | 0.015129 | 0.014152 | 0.039435 | -0.1483 | -0.1088 |
| MCF7 | `jac` | -0.0004 / -0.0004 / -0.0004 / -0.0004 / -0.0004 / -0.0004 / -0.0004 / -0.0004 | -0.000399 | 0.000000 | 0.000000 | 0.000000 | -0.0004 | -0.0004 |

Y el margen combinado de cada celda, con sus piezas a la vista:

| linea | metrica | sentido | u medio | sd_u (ddof=1) | b | r | span r-b | sd_r (5 splits) | SE_r | s medio | sd_s | termino de r | margen combinado |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | alto | 0.503433 | 0.082684 | 0.486742 | 0.777841 | 0.291098 | 0.029037 | 0.012986 | 0.057336 | 0.284040 | 0.002558 | **0.284051** |
| A549 | `mse_capped_norm` | bajo | 2.543222 | 0.069206 | 0.959928 | 0.138661 | -0.821267 | 0.061351 | 0.027437 | 0.000000 | 0.000000 | 0.000000 | **0.000000** |
| A549 | `lfc_nmae` | bajo | 1.015199 | 0.016964 | 0.856686 | 0.455829 | -0.400857 | 0.022871 | 0.010228 | -0.395435 | 0.042320 | 0.010090 | **0.043506** |
| A549 | `fid` | alto | 0.000000 | 0.000000 | 0.000000 | 0.417863 | 0.417863 | 0.036978 | 0.016537 | 0.000000 | 0.000000 | 0.000000 | **0.000000** |
| A549 | `reach` | alto | 0.169390 | 0.065082 | 0.314403 | 0.727540 | 0.413137 | 0.082754 | 0.037009 | -0.351003 | 0.157531 | 0.031443 | **0.160639** |
| A549 | `jac` | alto | 0.450758 | 0.010714 | 0.454545 | 0.473254 | 0.018709 | 0.035008 | 0.015656 | -0.202468 | 0.572666 | 0.169431 | **0.597205** |
| BXPC3 | `pds_cosine` | alto | 0.527225 | 0.042707 | 0.476326 | 0.955114 | 0.478788 | 0.016829 | 0.007526 | 0.106309 | 0.089199 | 0.001671 | **0.089215** |
| BXPC3 | `mse_capped_norm` | bajo | 1.350893 | 0.025707 | 0.824367 | 0.068734 | -0.755633 | 0.063862 | 0.028560 | 0.000000 | 0.000000 | 0.000000 | **0.000000** |
| BXPC3 | `lfc_nmae` | bajo | 1.008544 | 0.006673 | 0.814780 | 0.547191 | -0.267589 | 0.050873 | 0.022751 | -0.724108 | 0.024938 | 0.061566 | **0.066425** |
| BXPC3 | `fid` | alto | 0.000773 | 0.001740 | 0.001444 | 0.599253 | 0.597809 | 0.045793 | 0.020479 | -0.001123 | 0.002910 | 0.000038 | **0.002910** |
| BXPC3 | `reach` | alto | 0.078621 | 0.023701 | 0.344686 | 0.834502 | 0.489816 | 0.045994 | 0.020569 | -0.543194 | 0.048389 | 0.022811 | **0.053496** |
| BXPC3 | `jac` | alto | 0.000000 | 0.000000 | 0.000978 | 0.200769 | 0.199792 | 0.022772 | 0.010184 | -0.004893 | 0.000000 | 0.000249 | **0.000249** |
| HAP1 | `pds_cosine` | alto | 0.537405 | 0.053724 | 0.437500 | 0.789583 | 0.352083 | 0.015251 | 0.006820 | 0.283755 | 0.152589 | 0.005497 | **0.152688** |
| HAP1 | `mse_capped_norm` | bajo | 1.834121 | 0.050603 | 1.003465 | 0.185635 | -0.817830 | 0.091747 | 0.041031 | 0.000000 | 0.000000 | 0.000000 | **0.000000** |
| HAP1 | `lfc_nmae` | bajo | 1.002410 | 0.032256 | 0.897264 | 0.391590 | -0.505674 | 0.040663 | 0.018185 | -0.207934 | 0.063788 | 0.007478 | **0.064225** |
| HAP1 | `fid` | alto | 0.008443 | 0.017879 | 0.000000 | 0.546846 | 0.546846 | 0.094299 | 0.042172 | 0.015439 | 0.032694 | 0.001191 | **0.032716** |
| HAP1 | `reach` | alto | 0.157879 | 0.054033 | 0.354552 | 0.725538 | 0.370986 | 0.076025 | 0.034000 | -0.530135 | 0.145648 | 0.048585 | **0.153537** |
| HAP1 | `jac` | alto | 0.412879 | 0.015683 | 0.424242 | 0.468156 | 0.043914 | 0.084462 | 0.037772 | -0.258771 | 0.357138 | 0.222582 | **0.420821** |
| HT29 | `pds_cosine` | alto | 0.504380 | 0.050099 | 0.501894 | 0.783523 | 0.281629 | 0.017083 | 0.007640 | 0.008826 | 0.177889 | 0.000239 | **0.177889** |
| HT29 | `mse_capped_norm` | bajo | 1.953467 | 0.042701 | 1.137395 | 0.132822 | -1.004573 | 0.095598 | 0.042753 | 0.000000 | 0.000000 | 0.000000 | **0.000000** |
| HT29 | `lfc_nmae` | bajo | 1.016465 | 0.030492 | 0.864894 | 0.483357 | -0.381536 | 0.013128 | 0.005871 | -0.397264 | 0.079919 | 0.006113 | **0.080153** |
| HT29 | `fid` | alto | 0.004214 | 0.010814 | 0.114189 | 0.431793 | 0.317604 | 0.065468 | 0.029278 | -0.346264 | 0.034050 | 0.031920 | **0.046672** |
| HT29 | `reach` | alto | 0.207573 | 0.046983 | 0.327624 | 0.618619 | 0.290995 | 0.040970 | 0.018322 | -0.412553 | 0.161456 | 0.025976 | **0.163532** |
| HT29 | `jac` | alto | 0.178030 | 0.010714 | 0.091090 | 0.272116 | 0.181026 | 0.058478 | 0.026152 | 0.480264 | 0.059183 | 0.069382 | **0.091195** |
| K562 | `pds_cosine` | alto | 0.510417 | 0.044209 | 0.520833 | 0.822727 | 0.301894 | 0.017558 | 0.007852 | -0.034504 | 0.146438 | 0.000897 | **0.146440** |
| K562 | `mse_capped_norm` | bajo | 3.052408 | 0.046670 | 1.173286 | 0.126447 | -1.046839 | 0.105315 | 0.047098 | 0.000000 | 0.000000 | 0.000000 | **0.000000** |
| K562 | `lfc_nmae` | bajo | 1.014302 | 0.016161 | 0.846558 | 0.497343 | -0.349215 | 0.008552 | 0.003825 | -0.480347 | 0.046277 | 0.005261 | **0.046575** |
| K562 | `fid` | alto | 0.004808 | 0.013598 | 0.076923 | 0.511705 | 0.434782 | 0.044777 | 0.020025 | -0.165865 | 0.031276 | 0.007639 | **0.032195** |
| K562 | `reach` | alto | 0.231273 | 0.061179 | 0.406814 | 0.809303 | 0.402489 | 0.054991 | 0.024593 | -0.436141 | 0.152003 | 0.026649 | **0.154321** |
| K562 | `jac` | alto | 0.238636 | 0.010714 | 0.212121 | 0.489996 | 0.277874 | 0.067887 | 0.030360 | 0.095421 | 0.038556 | 0.010426 | **0.039941** |
| MCF7 | `pds_cosine` | alto | 0.526397 | 0.054783 | 0.503788 | 0.936174 | 0.432386 | 0.016343 | 0.007309 | 0.052289 | 0.126699 | 0.000884 | **0.126702** |
| MCF7 | `mse_capped_norm` | bajo | 1.757116 | 0.049443 | 0.922312 | 0.104958 | -0.817354 | 0.079845 | 0.035708 | 0.000000 | 0.000000 | 0.000000 | **0.000000** |
| MCF7 | `lfc_nmae` | bajo | 1.055097 | 0.007552 | 0.872196 | 0.663824 | -0.208372 | 0.011592 | 0.005184 | -0.877766 | 0.036243 | 0.021838 | **0.042314** |
| MCF7 | `fid` | alto | 0.000033 | 0.000058 | 0.000173 | 0.462605 | 0.462432 | 0.090324 | 0.040394 | -0.000303 | 0.000126 | 0.000027 | **0.000129** |
| MCF7 | `reach` | alto | 0.028280 | 0.008691 | 0.101892 | 0.676395 | 0.574503 | 0.063574 | 0.028431 | -0.128132 | 0.015129 | 0.006341 | **0.016404** |
| MCF7 | `jac` | alto | 0.000000 | 0.000000 | 0.000029 | 0.072243 | 0.072214 | 0.024237 | 0.010839 | -0.000399 | 0.000000 | 0.000060 | **0.000060** |

## 4. La `sd` acumulada con k = 2..8

La `sd` (ddof=1) de las primeras `k` semillas, en el orden `j = 0, 1, ... 7` --
que es el orden de la familia de semillas, fijado antes de medir y no una
reordenacion por valor. Responde a cuantas semillas hacen falta: si la columna se
aplana en k = 4, pedir 8 es gasto.

| linea | metrica | k=2 | k=3 | k=4 | k=5 | k=6 | k=7 | k=8 | primer k dentro del 20% de k=8 |
|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | 0.542867 | 0.403964 | 0.330263 | 0.301358 | 0.310990 | 0.300185 | 0.284040 | 4 |
| A549 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -- |
| A549 | `lfc_nmae` | 0.035077 | 0.026141 | 0.027287 | 0.038862 | 0.038935 | 0.040345 | 0.042320 | 2 |
| A549 | `fid` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -- |
| A549 | `reach` | 0.285422 | 0.241876 | 0.216434 | 0.194846 | 0.175198 | 0.160120 | 0.157531 | 6 |
| A549 | `jac` | 0 | 0 | 0 | 0 | 0.661258 | 0.612206 | 0.572666 | 6 |
| BXPC3 | `pds_cosine` | 0.036362 | 0.099307 | 0.096549 | 0.084013 | 0.075542 | 0.069541 | 0.089199 | 3 |
| BXPC3 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -- |
| BXPC3 | `lfc_nmae` | 0.034722 | 0.024583 | 0.020438 | 0.018791 | 0.019262 | 0.018748 | 0.024938 | 3 |
| BXPC3 | `fid` | 0.005422 | 0.004669 | 0.004111 | 0.003690 | 0.003330 | 0.003106 | 0.002910 | 6 |
| BXPC3 | `reach` | 0.071039 | 0.050236 | 0.067778 | 0.058789 | 0.055170 | 0.051757 | 0.048389 | 3 |
| BXPC3 | `jac` | 0 | 0 | 0 | 0 | 0.000000 | 0.000000 | 0 | -- |
| HAP1 | `pds_cosine` | 0.171166 | 0.158596 | 0.136039 | 0.143251 | 0.128345 | 0.154232 | 0.152589 | 2 |
| HAP1 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -- |
| HAP1 | `lfc_nmae` | 0.004526 | 0.005061 | 0.067802 | 0.067920 | 0.062474 | 0.065115 | 0.063788 | 4 |
| HAP1 | `fid` | 0 | 0 | 0 | 0 | 0.037328 | 0.034665 | 0.032694 | 6 |
| HAP1 | `reach` | 0.030385 | 0.146353 | 0.137732 | 0.122723 | 0.152301 | 0.146392 | 0.145648 | 3 |
| HAP1 | `jac` | 0.487943 | 0.398404 | 0.345028 | 0.308602 | 0.356343 | 0.336713 | 0.357138 | 3 |
| HT29 | `pds_cosine` | 0.316225 | 0.223624 | 0.212792 | 0.184293 | 0.172085 | 0.174621 | 0.177889 | 4 |
| HT29 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -- |
| HT29 | `lfc_nmae` | 0.024688 | 0.104907 | 0.091987 | 0.082298 | 0.074086 | 0.071961 | 0.079919 | 4 |
| HT29 | `fid` | 0 | 0 | 0 | 0 | 0.003662 | 0.036319 | 0.034050 | 7 |
| HT29 | `reach` | 0.127442 | 0.141236 | 0.117759 | 0.114212 | 0.111498 | 0.135411 | 0.161456 | 3 |
| HT29 | `jac` | 0.118367 | 0.096646 | 0.083698 | 0.074862 | 0.068339 | 0.063270 | 0.059183 | 6 |
| K562 | `pds_cosine` | 0.097593 | 0.138018 | 0.126044 | 0.127524 | 0.118679 | 0.135993 | 0.146438 | 3 |
| K562 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -- |
| K562 | `lfc_nmae` | 0.086153 | 0.074196 | 0.063720 | 0.055404 | 0.049879 | 0.046435 | 0.046277 | 5 |
| K562 | `fid` | 0 | 0 | 0.044231 | 0.039561 | 0.036114 | 0.033435 | 0.031276 | 6 |
| K562 | `reach` | 0.107033 | 0.148887 | 0.126068 | 0.123087 | 0.110591 | 0.159194 | 0.152003 | 3 |
| K562 | `jac` | 0 | 0.000000 | 0.054526 | 0.048770 | 0.044521 | 0.041218 | 0.038556 | 6 |
| MCF7 | `pds_cosine` | 0.103759 | 0.084093 | 0.152485 | 0.139451 | 0.130106 | 0.134597 | 0.126699 | 2 |
| MCF7 | `mse_capped_norm` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -- |
| MCF7 | `lfc_nmae` | 0.018474 | 0.030458 | 0.029748 | 0.029866 | 0.030946 | 0.038920 | 0.036243 | 3 |
| MCF7 | `fid` | 0.000041 | 0.000033 | 0.000029 | 0.000152 | 0.000140 | 0.000132 | 0.000126 | 6 |
| MCF7 | `reach` | 0.004200 | 0.016437 | 0.013682 | 0.013765 | 0.016499 | 0.015072 | 0.015129 | 3 |
| MCF7 | `jac` | 0 | 0 | 0 | 0 | 0 | 0 | 0 | -- |

De las **27 celdas con `sd > 0`**, el primer `k` cuya `sd` ya cae
dentro del 20% de la de `k = 8` es: A549/`pds_cosine` k=4, A549/`lfc_nmae` k=2, A549/`reach` k=6, A549/`jac` k=6, BXPC3/`pds_cosine` k=3, BXPC3/`lfc_nmae` k=3, BXPC3/`fid` k=6, BXPC3/`reach` k=3, HAP1/`pds_cosine` k=2, HAP1/`lfc_nmae` k=4, HAP1/`fid` k=6, HAP1/`reach` k=3, HAP1/`jac` k=3, HT29/`pds_cosine` k=4, HT29/`lfc_nmae` k=4, HT29/`fid` k=7, HT29/`reach` k=3, HT29/`jac` k=6, K562/`pds_cosine` k=3, K562/`lfc_nmae` k=5, K562/`fid` k=6, K562/`reach` k=3, K562/`jac` k=6, MCF7/`pds_cosine` k=2, MCF7/`lfc_nmae` k=3, MCF7/`fid` k=6, MCF7/`reach` k=3
.

Mediana de ese `k`: **4**; maximo **7**. Es una
descripcion de estas 48 corridas, no una recomendacion: el encargo mide la
dispersion y no propone cuantas semillas pedir.

⚠️ **Y una celda donde esa columna no significa lo que parece.** Donde la `sd`
arranca en 0 y salta de golpe (A549/`jac` 0 hasta k=5 y 0.661 en k=6; HAP1/`fid`
0 hasta k=5; HT29/`fid` 0 hasta k=5; K562/`fid` y K562/`jac` 0 hasta k=3) la
dispersion no es una nube que se estreche con mas draws: es un evento discreto
que ocurre o no (§5.2). Ahi `k` no dice cuando se estabiliza la `sd`, dice en que
draw cayo el evento. La `sd` acumulada de una Bernoulli con 8 draws no converge
en 4; estima una probabilidad con 8 observaciones.

## 5. Las `sd = 0`, separadas por causa

Son las celdas donde la emision no mueve nada. Tienen causas distintas y una de
ellas **no es estabilidad**.

**9 de 36 celdas tienen `sd = 0` exacta**, y 27 tienen dispersion.

**La particion es limpia y la decide una sola columna: `sd_u`.**

- `sd_u = 0`: la emision no mueve el `u` crudo. La `sd_s = 0` es **estructural**
  y es informacion real -- esa celda no tiene ruido de semilla.
- `sd_u > 0` con `sd_s = 0`: el `u` crudo SI se mueve y el `s` no. Algo entre los
  dos borra la variacion. Eso **no es estabilidad**.

Antes de la tabla, el aviso que cambia como se lee: el encargo esperaba que `jac`
diera `sd = 0` por razon estructural **en las seis lineas**, y solo la da en dos.
Las otras cuatro estan entre las celdas con dispersion, por un mecanismo discreto
que §5.2 nombra perturbacion por perturbacion.

| linea | metrica | u medio | **sd_u** | b | s | causa |
|---|---|---|---|---|---|---|
| A549 | `mse_capped_norm` | 2.543222 | 0.069206 | 0.959928 | 0.000000 | **clamp a [0,1]** (`catalog.py:415-416`) -- ceguera, no estabilidad |
| A549 | `fid` | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **estructural** -- `u` constante en 0.0000 |
| BXPC3 | `mse_capped_norm` | 1.350893 | 0.025707 | 0.824367 | 0.000000 | **clamp a [0,1]** (`catalog.py:415-416`) -- ceguera, no estabilidad |
| BXPC3 | `jac` | 0.000000 | 0.000000 | 0.000978 | -0.004893 | **estructural** -- `E = 0` en las 8: el lado real tiene senal en las 33, no hay 1.0 de convencion que perder (§5.2) |
| HAP1 | `mse_capped_norm` | 1.834121 | 0.050603 | 1.003465 | 0.000000 | **clamp a [0,1]** (`catalog.py:415-416`) -- ceguera, no estabilidad |
| HT29 | `mse_capped_norm` | 1.953467 | 0.042701 | 1.137395 | 0.000000 | **clamp a [0,1]** (`catalog.py:415-416`) -- ceguera, no estabilidad |
| K562 | `mse_capped_norm` | 3.052408 | 0.046670 | 1.173286 | 0.000000 | **clamp a [0,1]** (`catalog.py:415-416`) -- ceguera, no estabilidad |
| MCF7 | `mse_capped_norm` | 1.757116 | 0.049443 | 0.922312 | 0.000000 | **clamp a [0,1]** (`catalog.py:415-416`) -- ceguera, no estabilidad |
| MCF7 | `jac` | 0.000000 | 0.000000 | 0.000029 | -0.000399 | **estructural** -- `E = 0` en las 8: el lado real tiene senal en las 33, no hay 1.0 de convencion que perder (§5.2) |

Reparto: **3 estructurales**, **6 por el clamp** y **0 sin identificar**.

### 5.1 `mse_capped_norm`: ceguera, no estabilidad

**Las seis lineas dan `s = 0.000000` exacto en las 48 corridas.**

| linea | u medio | b | r | u peor que b | los 8 valores de u (CRUDOS) | sd de u | s | sd de s |
|---|---|---|---|---|---|---|---|---|
| A549 | 2.543222 | 0.959928 | 0.138661 | si | 2.442968 / 2.654257 / 2.619887 / 2.501768 / 2.528126 / 2.513136 / 2.508057 / 2.577581 | 0.069206 | 0.000000 | 0.000000 |
| BXPC3 | 1.350893 | 0.824367 | 0.068734 | si | 1.330246 / 1.375076 / 1.313458 / 1.368791 / 1.360340 / 1.329817 / 1.386820 / 1.342596 | 0.025707 | 0.000000 | 0.000000 |
| HAP1 | 1.834121 | 1.003465 | 0.185635 | si | 1.889359 / 1.778160 / 1.902936 / 1.871182 / 1.846043 / 1.800036 / 1.813666 / 1.771588 | 0.050603 | 0.000000 | 0.000000 |
| HT29 | 1.953467 | 1.137395 | 0.132822 | si | 1.933478 / 1.990297 / 1.900939 / 2.019131 / 1.913665 / 1.958209 / 1.921974 / 1.990046 | 0.042701 | 0.000000 | 0.000000 |
| K562 | 3.052408 | 1.173286 | 0.126447 | si | 3.067567 / 3.007006 / 3.027178 / 3.043949 / 3.049387 / 3.062562 / 3.008273 / 3.153339 | 0.046670 | 0.000000 | 0.000000 |
| MCF7 | 1.757116 | 0.922312 | 0.104958 | si | 1.727755 / 1.774452 / 1.721978 / 1.769669 / 1.662584 / 1.801134 / 1.805166 / 1.794193 | 0.049443 | 0.000000 | 0.000000 |

La columna que lo dice todo es `sd de u` contra `sd de s`: **el `u` crudo SI se
mueve de semilla a semilla**, y el `s` no se mueve nada. El clamp lo borra.

`catalog.py:415-416` fija `[0, 1]` sobre esta metrica y solo sobre esta. El suelo
esta por debajo de `b` en las seis lineas (columna `u peor que b`), asi que `s`
sale negativo antes del clamp y el clamp lo devuelve exactamente 0. **Una celda
que marca 0 con `sd = 0` aqui no informa de que un candidato este mejorando: dos
candidatos distintos, los dos por debajo de `b`, son indistinguibles.** Es la
propiedad que `descomposicion.md` ya habia declarado, ahora medida con 8 semillas:
la dispersion del `u` crudo es la prueba de que hay senal que el clamp tira.

### 5.2 `jac`: la prediccion del encargo sale REFUTADA en 4 de 6 lineas

El encargo predijo que la `E` del suelo no depende de la semilla de emision,
porque es el complemento de la cohorte del lado real (§3.15), y que por eso la
`sd` de `jac` deberia ser 0 por razon estructural, **y pidio comprobarlo y, si NO
era 0, decir que de la emision mueve la union.** No es 0:

- `sd = 0` solo en **BXPC3 y MCF7**, y ahi es estructural por una razon mas fuerte
  que la predicha: `E = 0` en las 8 semillas. El lado real tiene senal DE en las 33
  perturbaciones, asi que no hay ningun 1.0 de convencion que perder.
- En **A549, HAP1, HT29 y K562 la `E` SI se mueve con la semilla de emision**, y con
  ella `jac`. `sd_u` = 0.010714, 0.015683, 0.010714, 0.010714. La emision mueve la
  union.

**Que la mueve, nombrado.** La emision predice algun gen significativo en 18 de las
48 corridas. Eso por si solo no cambia `E`: solo la cambia cuando la perturbacion
predicha es una cuyo lado real esta VACIO. Entonces esa perturbacion pierde su 1.0
gratis y cobra un jaccard genuino de 0 (interseccion vacia con un predicho no
vacio). Las ocurrencias, una por una:

| linea | semilla | perturbacion que salio de la union vacia | E antes -> despues | mecanismo |
|---|---|---|---|---|
| A549 | j=5 | **STAT5A** | 15 -> 14 | predicho significativo sobre real vacio |
| HAP1 | j=0 | **RAP1GAP** | 14 -> 13 | predicho significativo sobre real vacio |
| HAP1 | j=5 | **BRD9** | 14 -> 13 | predicho significativo sobre real vacio |
| HAP1 | j=7 | **SP100** | 14 -> 13 | predicho significativo sobre real vacio |
| HT29 | j=1 | **RAP1GAP** | 6 -> 5 | predicho significativo sobre real vacio |
| K562 | j=3 | **SOCS1** | 8 -> 7 | predicho significativo sobre real vacio |

Son los mismos nombres que `huecos.md` §2.5 habia aislado para el `b` (BRD9, IRF7,
RAP1GAP, ZBP1): el conjunto de perturbaciones donde esto puede pasar es el
complemento de la cohorte del lado real, y es el lado real el que lo fija. Lo que
cambia con la semilla es CUAL de ellas cae, no cuales son candidatas.

**Y la amplificacion, que es el hallazgo caro.** Un cambio de `E` de uno mueve el
`jac` crudo en exactamente `1/33 = 0.030303`. Dividido por el span de la celda, eso
es el salto en `s`:

| linea | span `r - b` | `(1/33)/span` | salto de `s` observado | coincide |
|---|---|---|---|---|
| A549 | 0.018709 | **1.619745** | 0.0000 -> -1.619745 | exacto |
| HAP1 | 0.043914 | **0.690056** | 0.0000 -> -0.690056 | exacto |
| HT29 | 0.181026 | **0.167396** | 0.501189 -> 0.333793 | exacto |
| K562 | 0.277874 | **0.109053** | 0.109053 -> 0.0000 | exacto |

En A549 el span entero (0.0187) es **menos que una unidad de convencion** (0.0303):
una sola perturbacion de 33 que cruce el umbral de significancia mueve `s` en 1.62
puntos. Repartido a 1/6 en el agregado son 0.270 de `avg_score` -- y el rango de
`avg_score` de A549 es 0.2476, el mayor de las seis lineas. **Casi toda la
dispersion del agregado de A549 es ese unico evento discreto**, no ruido continuo de
la emision: `j=5` da -0.3190 y las otras siete caen entre -0.0714 y -0.1831.

Consecuencia para un candidato, dicha sin proponer nada: en las celdas de `jac` con
span pequeno la `sd` no describe una nube, describe una moneda. Promediar 8 semillas
ahi no estrecha un intervalo gaussiano; estima la probabilidad de un salto.

Conviene separar dos cosas que no son la misma: que `sd = 0` **no** quiere decir
`s = 0`. En BXPC3 y MCF7, las dos con `sd = 0`, el `s` sale -0.004893 y -0.000399
-- practicamente 0 porque `u = 0` y el `b` tambien es casi 0. En HT29 y K562 el `u`
esta POR ENCIMA del `b` -- las cuatro celdas de `huecos.md` §2.5 -- y `s` sale
positivo (+0.4803 y +0.0954 de media): un puntaje que una copia del control cobra
por convencion, aunque ahi la semilla si lo mueva.

| linea | E por semilla (j=0..7) | valores distintos de E | union no vacia por semilla | CONJUNTOS distintos de union vacia | u medio de jac | b de jac | sd de u | sd de s |
|---|---|---|---|---|---|---|---|---|
| A549 | 15/15/15/15/15/14/15/15 | 2 | 18/18/18/18/18/19/18/18 | 2 | 0.450758 | 0.454545 | 0.010714 | 0.572666 |
| BXPC3 | 0/0/0/0/0/0/0/0 | 1 | 33/33/33/33/33/33/33/33 | 1 | 0.000000 | 0.000978 | 0.000000 | 0.000000 |
| HAP1 | 13/14/14/14/14/13/14/13 | 2 | 20/19/19/19/19/20/19/20 | 4 | 0.412879 | 0.424242 | 0.015683 | 0.357138 |
| HT29 | 6/5/6/6/6/6/6/6 | 2 | 27/28/27/27/27/27/27/27 | 2 | 0.178030 | 0.091090 | 0.010714 | 0.059183 |
| K562 | 8/8/8/7/8/8/8/8 | 2 | 25/25/25/26/25/25/25/25 | 2 | 0.238636 | 0.212121 | 0.010714 | 0.038556 |
| MCF7 | 0/0/0/0/0/0/0/0 | 1 | 33/33/33/33/33/33/33/33 | 1 | 0.000000 | 0.000029 | 0.000000 | 0.000000 |

La columna `CONJUNTOS distintos` es la que hace la comprobacion: donde vale 1
(BXPC3, MCF7) son exactamente las mismas perturbaciones en las 8 semillas; donde
vale 2 o 4 (A549, HT29, K562 y HAP1) **no lo son**, y el conteo coincide con el
numero de configuraciones distintas de la tabla de ocurrencias de arriba.

La parte de la prediccion que SI se sostiene es la mas de fondo: el CONJUNTO de
perturbaciones donde `E` puede moverse es el complemento de la cohorte del lado
real, y eso no depende de la emision. `huecos.md` §2.5 midio para nuestro suelo la
igualdad `E_u = 33 - (tamano de la cohorte del lado real)`; lo que aqui se corrige
es que esa igualdad se cumple **por semilla pero no identicamente en todas**: el
suelo deja de ser una copia exacta del `b` en 6 de las 48 corridas.

## 6. Las dos cohortes moviles, por semilla

`fid` es la unica de las seis cuyo gate puede moverlo la submission: su
denominador es `max(n_pred, N_conf)` (`direction.py:776`) y `_safe_div` devuelve
NaN cuando vale 0 (`direction.py:689-699`), asi que su cohorte es
`max(n_pred, N_conf) > 0`. Las otras cinco dependen del lado real solo.

| linea | cohorte de fid por semilla | rango fid | union no vacia (jac) por semilla | cohorte de reach | cohorte de lfc_nmae | perts con algun gen predicho significativo | suma de n_DE_pred |
|---|---|---|---|---|---|---|---|
| A549 | 18/18/18/18/18/19/18/18 | 18-19 | 18/18/18/18/18/19/18/18 | 18/18/18/18/18/18/18/18 | 8/8/8/8/8/8/8/8 | 0/0/0/0/0/1/0/0 | 0/0/0/0/0/1/0/0 |
| BXPC3 | 33/33/33/33/33/33/33/33 | 33-33 | 33/33/33/33/33/33/33/33 | 33/33/33/33/33/33/33/33 | 23/23/23/23/23/23/23/23 | 1/1/0/0/2/3/1/1 | 1/1/0/0/2/3/1/1 |
| HAP1 | 20/19/19/19/19/20/19/20 | 19-20 | 20/19/19/19/19/20/19/20 | 19/19/19/19/19/19/19/19 | 8/8/8/8/8/8/8/8 | 2/0/0/0/0/1/1/1 | 2/0/0/0/0/1/1/1 |
| HT29 | 27/28/27/27/27/27/27/27 | 27-28 | 27/28/27/27/27/27/27/27 | 27/27/27/27/27/27/27/27 | 15/15/15/15/15/15/15/15 | 0/1/0/0/0/1/2/0 | 0/1/0/0/0/1/2/0 |
| K562 | 25/25/25/26/25/25/25/25 | 25-26 | 25/25/25/26/25/25/25/25 | 25/25/25/25/25/25/25/25 | 9/9/9/9/9/9/9/9 | 0/0/0/1/0/0/0/0 | 0/0/0/1/0/0/0/0 |
| MCF7 | 33/33/33/33/33/33/33/33 | 33-33 | 33/33/33/33/33/33/33/33 | 33/33/33/33/33/33/33/33 | 33/33/33/33/33/33/33/33 | 1/0/0/0/1/1/0/0 | 1/0/0/0/1/1/0/0 |

**La cohorte de `fid` SI se mueve sola, solo por la semilla de emision**, en A549 (18-19), HAP1 (19-20), HT29 (27-28), K562 (25-26).

Las de `reach` y `lfc_nmae` no se mueven en ninguna linea, como tiene que ser:
son del lado real y el lado real esta fijo (§1).

Esto es ruido del candidato, no del ancla: la misma emision con otra semilla
mete o saca perturbaciones del promedio de `fid`. Hay que verlo en cada
corrida, que es lo que §3.16 pide.

## 7. Las 20 filas trabajadas

**Criterio de seleccion, declarado antes de mirar los valores**: las celdas de
`pds_cosine`, `lfc_nmae` y `reach` en las seis lineas = 18 filas, mas las dos de
`jac` de HT29 y K562 -- las dos lineas donde `huecos.md` §2.5 ya habia identificado,
ANTES de este encargo, que la `E` del `b` no coincide con la de nuestro suelo
(BRD9/IRF7/RAP1GAP y ZBP1) = **20**. Es un criterio por metrica y por hallazgo
previo, no por el valor medido aqui.

Cada fila trae los 8 `u` crudos, los 8 `s`, y la aritmetica `(u - b)/(r - b)`
rehecha a mano sobre la media para que se vea que el `s` reportado es ese.

| linea | metrica | los 8 u crudos | u medio | sd_u | b | r | span | s a mano desde u medio | s medio medido | \|dif\| | sd_s | termino de r | margen |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | 0.6231 / 0.3996 / 0.4479 / 0.5000 / 0.4309 / 0.5909 / 0.5739 / 0.4612 | 0.503433 | 0.082684 | 0.486742 | 0.777841 | 0.291098 | 0.057336 | 0.057336 | 2.8e-17 | 0.284040 | 0.002558 | 0.284051 |
| A549 | `lfc_nmae` | 1.0040 / 1.0239 / 1.0197 / 1.0022 / 1.0401 / 1.0008 / 1.0354 / 0.9955 | 1.015199 | 0.016964 | 0.856686 | 0.455829 | -0.400857 | -0.395435 | -0.395435 | 1.1e-16 | 0.042320 | 0.010090 | 0.043506 |
| A549 | `reach` | 0.0956 / 0.2624 / 0.0836 / 0.2204 / 0.2147 / 0.1935 / 0.1699 / 0.1149 | 0.169390 | 0.065082 | 0.314403 | 0.727540 | 0.413137 | -0.351003 | -0.351003 | 5.6e-17 | 0.157531 | 0.031443 | 0.160639 |
| BXPC3 | `pds_cosine` | 0.5133 / 0.4886 / 0.5805 / 0.4773 / 0.5237 / 0.5076 / 0.5265 / 0.6004 | 0.527225 | 0.042707 | 0.476326 | 0.955114 | 0.478788 | 0.106309 | 0.106309 | 1.4e-17 | 0.089199 | 0.001671 | 0.089215 |
| BXPC3 | `lfc_nmae` | 1.0029 / 1.0160 / 1.0100 / 1.0076 / 1.0129 / 1.0160 / 1.0063 / 0.9967 | 1.008544 | 0.006673 | 0.814780 | 0.547191 | -0.267589 | -0.724108 | -0.724108 | 4.4e-16 | 0.024938 | 0.061566 | 0.066425 |
| BXPC3 | `reach` | 0.1193 / 0.0701 / 0.0952 / 0.0420 / 0.0780 / 0.0609 / 0.0930 / 0.0705 | 0.078621 | 0.023701 | 0.344686 | 0.834502 | 0.489816 | -0.543194 | -0.543194 | 1.1e-16 | 0.048389 | 0.022811 | 0.053496 |
| HAP1 | `pds_cosine` | 0.5000 / 0.5852 / 0.4801 / 0.5511 / 0.4650 / 0.5227 / 0.6108 / 0.5843 | 0.537405 | 0.053724 | 0.437500 | 0.789583 | 0.352083 | 0.283755 | 0.283755 | 1.7e-16 | 0.152589 | 0.005497 | 0.152688 |
| HAP1 | `lfc_nmae` | 1.0237 / 1.0270 / 1.0288 / 0.9581 / 0.9708 / 1.0197 / 0.9626 / 1.0285 | 1.002410 | 0.032256 | 0.897264 | 0.391590 | -0.505674 | -0.207934 | -0.207934 | 1.7e-16 | 0.063788 | 0.007478 | 0.064225 |
| HAP1 | `reach` | 0.1768 / 0.1609 / 0.0758 / 0.1886 / 0.1745 / 0.2513 / 0.1263 / 0.1089 | 0.157879 | 0.054033 | 0.354552 | 0.725538 | 0.370986 | -0.530135 | -0.530135 | 1.1e-16 | 0.145648 | 0.048585 | 0.153537 |
| HT29 | `pds_cosine` | 0.5928 / 0.4669 / 0.5284 / 0.4678 / 0.5152 / 0.5483 / 0.4631 / 0.4527 | 0.504380 | 0.050099 | 0.501894 | 0.783523 | 0.281629 | 0.008826 | 0.008826 | 5.4e-17 | 0.177889 | 0.000239 | 0.177889 |
| HT29 | `lfc_nmae` | 1.0596 / 1.0463 / 0.9846 / 1.0046 / 1.0414 / 1.0194 / 1.0011 / 0.9748 | 1.016465 | 0.030492 | 0.864894 | 0.483357 | -0.381536 | -0.397264 | -0.397264 | 2.8e-16 | 0.079919 | 0.006113 | 0.080153 |
| HT29 | `reach` | 0.1941 / 0.2466 / 0.1655 / 0.2160 / 0.2390 / 0.1804 / 0.1382 / 0.2808 | 0.207573 | 0.046983 | 0.327624 | 0.618619 | 0.290995 | -0.412553 | -0.412553 | 5.6e-17 | 0.161456 | 0.025976 | 0.163532 |
| K562 | `pds_cosine` | 0.5322 / 0.4905 / 0.5739 / 0.5663 / 0.4962 / 0.5076 / 0.4621 / 0.4545 | 0.510417 | 0.044209 | 0.520833 | 0.822727 | 0.301894 | -0.034504 | -0.034504 | 2.6e-16 | 0.146438 | 0.000897 | 0.146440 |
| K562 | `lfc_nmae` | 0.9932 / 1.0358 / 0.9889 / 1.0198 / 1.0133 / 1.0150 / 1.0194 / 1.0291 | 1.014302 | 0.016161 | 0.846558 | 0.497343 | -0.349215 | -0.480347 | -0.480347 | 0.0e+00 | 0.046277 | 0.005261 | 0.046575 |
| K562 | `reach` | 0.1628 / 0.2237 / 0.2827 / 0.2500 / 0.1786 / 0.2092 / 0.3489 / 0.1942 | 0.231273 | 0.061179 | 0.406814 | 0.809303 | 0.402489 | -0.436141 | -0.436141 | 0.0e+00 | 0.152003 | 0.026649 | 0.154321 |
| MCF7 | `pds_cosine` | 0.5890 / 0.5256 / 0.5881 / 0.4498 / 0.5814 / 0.5076 / 0.4678 / 0.5019 | 0.526397 | 0.054783 | 0.503788 | 0.936174 | 0.432386 | 0.052289 | 0.052289 | 3.5e-17 | 0.126699 | 0.000884 | 0.126702 |
| MCF7 | `lfc_nmae` | 1.0537 / 1.0483 / 1.0609 / 1.0611 / 1.0490 / 1.0466 / 1.0680 / 1.0531 | 1.055097 | 0.007552 | 0.872196 | 0.663824 | -0.208372 | -0.877766 | -0.877766 | 2.2e-16 | 0.036243 | 0.021838 | 0.042314 |
| MCF7 | `reach` | 0.0348 / 0.0382 / 0.0204 / 0.0281 / 0.0394 / 0.0167 / 0.0287 / 0.0200 | 0.028280 | 0.008691 | 0.101892 | 0.676395 | 0.574503 | -0.128132 | -0.128132 | 2.8e-17 | 0.015129 | 0.006341 | 0.016404 |
| HT29 | `jac` | 0.1818 / 0.1515 / 0.1818 / 0.1818 / 0.1818 / 0.1818 / 0.1818 / 0.1818 | 0.178030 | 0.010714 | 0.091090 | 0.272116 | 0.181026 | 0.480264 | 0.480264 | 0.0e+00 | 0.059183 | 0.069382 | 0.091195 |
| K562 | `jac` | 0.2424 / 0.2424 / 0.2424 / 0.2121 / 0.2424 / 0.2424 / 0.2424 / 0.2424 | 0.238636 | 0.010714 | 0.212121 | 0.489996 | 0.277874 | 0.095421 | 0.095421 | 1.4e-17 | 0.038556 | 0.010426 | 0.039941 |

Son **20** filas. La columna `|dif|` compara el `s` que sale de rehacer
`(u_medio - b)/(r - b)` contra la media de los 8 `s` que escribio `score`. Donde
no es cero del todo es porque la media de una funcion lineal es la funcion de la
media solo si no hay clamp ni gates por el camino -- y donde el clamp actua
(`mse`, que no esta en esta seleccion) la identidad se rompe del todo.

## Lo que no se hizo

- **No se construyo ningun candidato.** La emision es `predict_floor` sin
  modificar, la misma desechable de `cable_suelo.md`.
- **No se tocaron los bundles ni ningun arnes.** `build_fold_adatas`,
  `predict_floor` y `contar_de_por_perturbacion` se usaron tal cual.
- **No se promediaron las 8 corridas en una submission.** Se puntuaron por
  separado contra el mismo bundle; lo que se reporta es su dispersion.
- **No se propone ningun margen ni umbral.** La seccion 2 da la formula y la
  seccion 3 la aplica celda por celda; la regla de decision ya esta escrita en
  §3.16 y no se toca aqui.
- No se leyo nada de `~/vcc/benchmark/` ni de
  `/datos/benchmark_VEDADO_techo_metricas/`.
- Sin envios, sin `vcc sample`, sin commits.

