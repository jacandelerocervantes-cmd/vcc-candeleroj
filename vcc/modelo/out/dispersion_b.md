# La dispersion de `b`

El tercer termino del margen, que hasta hoy faltaba. `dispersion_u.md` §2 dejo
escrito que su margen OMITE el termino de `b` porque `b` era UNA sola
realizacion con `--seed 0`, y que mientras lo fuera todo margen seria cota
inferior. Aqui se estima: `cell-eval2 baseline` con **8 semillas** por linea
sobre la MISMA poblacion real. 48 corridas, ninguna fallida.

Con este termino, la resolucion de cada celda deja de ser una cota y pasa a ser
un numero.

## 1. Lo que vario y lo que no, con la prueba

**Vario solo la semilla del baseline.** Familia nueva `20265000 + 100*fold + j`,
`j = 0..7`. Ninguna colision con `2026000x` (verdad), `2026100x` (suelo),
`2026200x`/`2026300x` (techo), `2026400x` (la dispersion de `u`) ni la `0` del
baseline/anchor.

**No vario la poblacion real.** `build_fold_adatas` sin modificar con
`seed_real = 20260001 + fold_idx`, etiqueta `non-targeting`. Se comprueba con el
mismo `sha256` de `dispersion_u.md` (matriz `data`/`indices`/`indptr` mas las
etiquetas), y el script aborta si no coincide.

🔴 **Y la huella es la MISMA que la de `dispersion_u.md`**, no solo consistente
dentro de esta tanda. Eso es lo que permite componer `sd_b` con `sd_s` en un solo
margen: los dos terminos se midieron sobre la misma poblacion, celda por celda.

| linea | corridas | semillas del baseline | huellas distintas del REAL | huellas distintas del EMITIDO | huella del real = la de dispersion_u.md | config_hash distintos entre las 8 | source_fingerprint_strict |
|---|---|---|---|---|---|---|---|
| A549 | 8 | 20265000..20265007 | **1** | **8** | **si** | **1** | **si** |
| BXPC3 | 8 | 20265100..20265107 | **1** | **8** | **si** | **1** | **si** |
| HAP1 | 8 | 20265200..20265207 | **1** | **8** | **si** | **1** | **si** |
| HT29 | 8 | 20265300..20265307 | **1** | **8** | **si** | **1** | **si** |
| K562 | 8 | 20265400..20265407 | **1** | **8** | **si** | **1** | **si** |
| MCF7 | 8 | 20265500..20265507 | **1** | **8** | **si** | **1** | **si** |

**Una sola huella del real por linea y ocho del h5ad que el baseline emite**
(`--save-pred`): la poblacion esta fija y la emision del baseline cambia en las
ocho. Es la condicion que hace que la dispersion medida tenga una sola fuente.

**La semilla no entra en el `config_hash`.** `run_params.yaml` del baseline no
contiene `seed` -- es parametro de la emision, no de `EvalConfig` -- y el hash
de regla sale identico en las 8 semillas de cada linea: `9d36eb547b5c650d8684cc7e...`.
Ese valor **no** es el de competencia (`a5519017556705ab...`) y no tiene por que
serlo: es el del brazo baseline, y coincide **bit a bit con el del `BASE_<linea>`
oficial** que se uso para construir los bundles. Verificado recalculandolo desde
los dos `run_params.yaml`. Lo que esto descarta es que cambiar la semilla cambie
la configuracion; lo unico que cambia es el sorteo.

`source_fingerprint_strict` sale `true` en las 48.

## 2. El margen con sus TRES terminos, derivado

La escala es `s = (u - b)/(r - b)` (`score.py:395-399`). Las tres cantidades son
mediciones, no constantes. Las derivadas de primer orden, escritas:

```
s = (u - b) / (r - b)

ds/du = 1/(r - b)

ds/dr = -(u - b)/(r - b)^2            = -s/(r - b)

ds/db = -1/(r - b) + (u - b)/(r - b)^2
      = [-(r - b) + (u - b)] / (r - b)^2
      = (u - r)/(r - b)^2             = (s - 1)/(r - b)
```

La forma compacta `ds/db = (s - 1)/(r - b)` dice dos cosas que la forma larga
esconde:

- En `s = 0` (un candidato pegado al baseline) la derivada vale `-1/(r - b)`, el
  mismo tamano que `ds/du` y de signo opuesto. Ahi `b` pesa tanto como la propia
  emision.
- En `s = 1` (un candidato pegado al ancla) la derivada vale **0**. Un candidato
  que alcanza el ancla es insensible a donde cayo el baseline.

Las 48 corridas de `dispersion_u.md` dan `s` entre **-1.6197 y +0.5012** (el
minimo es el salto de convencion de A549/`jac`, §5.2 de ese reporte); por celda,
las medias corren de -0.878 a +0.480. En los dos casos `|s - 1| >= 0.50`, y llega
a **2.62**. Como `|ds/du| = 1/|r - b|` y `|ds/db| = |s - 1|/|r - b|`, el factor que
los separa es exactamente `|s - 1|`: el termino de `b` pesa **mas** que el de `u`
en toda celda con `s < 0`, que son la mayoria de las medidas.

Componiendo en cuadratura, con `u`, `r` y `b` tratados como independientes:

```
margen_celda = sqrt(  sd_s^2
                    + ( |s|     * SE_r / |r - b| )^2
                    + ( |s - 1| * sd_b / |r - b| )^2  )
```

🔴 **Y la asimetria entre `SE_r` y `sd_b`, que no es un descuido.** Para `r` entra
`SE_r = sd_r/sqrt(5)` porque el artefacto del bundle guarda la MEDIA de 5
particiones: el numero que usa el marcador ya es un promedio. Para `b` entra la
`sd_b` ENTERA, sin dividir por `sqrt(8)`, porque el bundle guarda **un solo**
sorteo con `--seed 0` -- las 8 semillas estiman cuanto se mueve ese sorteo, no lo
promedian. Dividir por `sqrt(8)` seria describir un `b` promediado que no existe
en ningun bundle.

**Lo que sigue sin estimarse, declarado.** La poblacion real es una sola
realizacion (`seed_real = 20260001 + fold_idx`) y los tres terminos se miden
sobre ella. El margen de arriba es el ruido de los sorteos DENTRO de esa
poblacion; no cubre cambiar la poblacion. Eso sigue siendo cota inferior respecto
al ruido total, y queda dicho.

## 3. Las 36 celdas: los 8 valores de `b`

`b` es la columna `mean` de `baseline_agg.csv` de cada corrida. El `b` de semilla
0 -- el que esta DENTRO de los 6 bundles de competencia y que no se sustituye ni
se promedia -- se localiza dentro de la distribucion de las 8, marcado como
`sem0`.

| linea | metrica | sentido | los 8 valores de b | media | sd (ddof=1) | sd (ddof=0) | rango | min | max |
|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | alto | 0.4924 / 0.5265 / 0.4962 / 0.4280 / 0.4839 / 0.4527 / 0.5028 / 0.4669 | 0.481179 | 0.031003 | 0.029000 | 0.098485 | 0.4280 | 0.5265 |
| A549 | `mse_capped_norm` | bajo | 1.0641 / 0.9842 / 0.9997 / 1.1385 / 1.0078 / 1.0394 / 0.9508 / 1.0882 | 1.034074 | 0.061015 | 0.057074 | 0.187705 | 0.9508 | 1.1385 |
| A549 | `lfc_nmae` | bajo | 0.8938 / 0.8526 / 0.8392 / 0.8493 / 0.8718 / 0.8870 / 0.8556 / 0.8491 | 0.862296 | 0.019665 | 0.018395 | 0.054624 | 0.8392 | 0.8938 |
| A549 | `fid` | alto | 0.0000 / 0.0000 / 0.0000 / 0.0000 / 0.0000 / 0.1053 / 0.0000 / 0.0000 | 0.013158 | 0.037216 | 0.034813 | 0.105263 | 0.0000 | 0.1053 |
| A549 | `reach` | alto | 0.3515 / 0.2708 / 0.4463 / 0.3916 / 0.3938 / 0.4223 / 0.5443 / 0.4784 | 0.412372 | 0.082408 | 0.077086 | 0.273459 | 0.2708 | 0.5443 |
| A549 | `jac` | alto | 0.4545 / 0.4545 / 0.4545 / 0.4545 / 0.4545 / 0.4242 / 0.4545 / 0.4545 | 0.450758 | 0.010714 | 0.010022 | 0.030303 | 0.4242 | 0.4545 |
| BXPC3 | `pds_cosine` | alto | 0.5350 / 0.5540 / 0.5246 / 0.4991 / 0.4962 / 0.5009 / 0.5218 / 0.5492 | 0.522609 | 0.022594 | 0.021135 | 0.057765 | 0.4962 | 0.5540 |
| BXPC3 | `mse_capped_norm` | bajo | 0.7557 / 0.7577 / 0.7817 / 0.8009 / 0.8338 / 0.8566 / 0.7856 / 0.7459 | 0.789725 | 0.039148 | 0.036620 | 0.110702 | 0.7459 | 0.8566 |
| BXPC3 | `lfc_nmae` | bajo | 0.8095 / 0.8092 / 0.8114 / 0.8159 / 0.8237 / 0.8107 / 0.8156 / 0.8104 | 0.813297 | 0.004915 | 0.004597 | 0.014434 | 0.8092 | 0.8237 |
| BXPC3 | `fid` | alto | 0.0000 / 0.0011 / 0.0065 / 0.0001 / 0.0000 / 0.0004 / 0.0018 / 0.0013 | 0.001392 | 0.002157 | 0.002018 | 0.006470 | 0.0000 | 0.0065 |
| BXPC3 | `reach` | alto | 0.2721 / 0.2794 / 0.2711 / 0.1740 / 0.1878 / 0.3034 / 0.3186 / 0.3549 | 0.270184 | 0.061807 | 0.057815 | 0.180861 | 0.1740 | 0.3549 |
| BXPC3 | `jac` | alto | 0.0000 / 0.0003 / 0.0000 / 0.0000 / 0.0000 / 0.0004 / 0.0000 / 0.0004 | 0.000139 | 0.000188 | 0.000176 | 0.000421 | 0.0000 | 0.0004 |
| HAP1 | `pds_cosine` | alto | 0.5445 / 0.5133 / 0.5407 / 0.4782 / 0.4432 / 0.4782 / 0.5208 / 0.4943 | 0.501657 | 0.034716 | 0.032473 | 0.101326 | 0.4432 | 0.5445 |
| HAP1 | `mse_capped_norm` | bajo | 1.0012 / 0.8817 / 0.8252 / 0.9503 / 0.9745 / 0.9560 / 0.9487 / 1.0283 | 0.945743 | 0.064897 | 0.060706 | 0.203096 | 0.8252 | 1.0283 |
| HAP1 | `lfc_nmae` | bajo | 0.8904 / 0.9133 / 0.8867 / 0.8942 / 0.9110 / 0.8900 / 0.9151 / 0.9240 | 0.903099 | 0.014292 | 0.013369 | 0.037304 | 0.8867 | 0.9240 |
| HAP1 | `fid` | alto | 0.0000 / 0.0000 / 0.0000 / 0.0474 / 0.0000 / 0.0075 / 0.0000 / 0.0500 | 0.013111 | 0.022120 | 0.020692 | 0.050000 | 0.0000 | 0.0500 |
| HAP1 | `reach` | alto | 0.2947 / 0.2741 / 0.1857 / 0.3164 / 0.2322 / 0.4221 / 0.2819 / 0.2875 | 0.286822 | 0.068276 | 0.063866 | 0.236402 | 0.1857 | 0.4221 |
| HAP1 | `jac` | alto | 0.4242 / 0.3939 / 0.4242 / 0.4242 / 0.4242 / 0.4242 / 0.4242 / 0.3939 | 0.416667 | 0.014028 | 0.013122 | 0.030303 | 0.3939 | 0.4242 |
| HT29 | `pds_cosine` | alto | 0.5540 / 0.4451 / 0.5426 / 0.5350 / 0.5095 / 0.5133 / 0.4953 / 0.4905 | 0.510653 | 0.034730 | 0.032487 | 0.108902 | 0.4451 | 0.5540 |
| HT29 | `mse_capped_norm` | bajo | 1.0021 / 1.1076 / 1.0386 / 1.0895 / 1.0926 / 1.0630 / 1.0928 / 1.0887 | 1.071865 | 0.035443 | 0.033154 | 0.105559 | 1.0021 | 1.1076 |
| HT29 | `lfc_nmae` | bajo | 0.9106 / 0.9180 / 0.8572 / 0.8656 / 0.8560 / 0.8832 / 0.8622 / 0.8882 | 0.880140 | 0.024133 | 0.022574 | 0.062018 | 0.8560 | 0.9180 |
| HT29 | `fid` | alto | 0.0279 / 0.0788 / 0.1244 / 0.0604 / 0.0327 / 0.0632 / 0.1166 / 0.0161 | 0.065012 | 0.040029 | 0.037444 | 0.108300 | 0.0161 | 0.1244 |
| HT29 | `reach` | alto | 0.3187 / 0.3424 / 0.3394 / 0.3067 / 0.2499 / 0.3340 / 0.3281 / 0.3026 | 0.315222 | 0.030139 | 0.028192 | 0.092529 | 0.2499 | 0.3424 |
| HT29 | `jac` | alto | 0.1819 / 0.1517 / 0.0911 / 0.1821 / 0.1517 / 0.1212 / 0.1213 / 0.1820 | 0.147889 | 0.034155 | 0.031949 | 0.091006 | 0.0911 | 0.1821 |
| K562 | `pds_cosine` | alto | 0.4612 / 0.4915 / 0.5047 / 0.4451 / 0.5294 / 0.5417 / 0.5199 / 0.4896 | 0.497869 | 0.033133 | 0.030993 | 0.096591 | 0.4451 | 0.5417 |
| K562 | `mse_capped_norm` | bajo | 1.1884 / 1.0876 / 1.1818 / 1.1553 / 1.1504 / 1.2112 / 1.1024 / 1.1684 | 1.155684 | 0.042256 | 0.039527 | 0.123533 | 1.0876 | 1.2112 |
| K562 | `lfc_nmae` | bajo | 0.8511 / 0.8783 / 0.8571 / 0.8477 / 0.8575 / 0.8689 / 0.8092 / 0.8590 | 0.853596 | 0.020413 | 0.019094 | 0.069144 | 0.8092 | 0.8783 |
| K562 | `fid` | alto | 0.0000 / 0.0000 / 0.0000 / 0.0000 / 0.0000 / 0.0000 / 0.0006 / 0.0000 | 0.000078 | 0.000221 | 0.000207 | 0.000625 | 0.0000 | 0.0006 |
| K562 | `reach` | alto | 0.5890 / 0.4616 / 0.5093 / 0.3914 / 0.5303 / 0.4072 / 0.4393 / 0.4818 | 0.476232 | 0.065801 | 0.061551 | 0.197583 | 0.3914 | 0.5890 |
| K562 | `jac` | alto | 0.2424 / 0.2424 / 0.2424 / 0.2424 / 0.2424 / 0.2424 / 0.2429 / 0.2424 | 0.242483 | 0.000167 | 0.000157 | 0.000473 | 0.2424 | 0.2429 |
| MCF7 | `pds_cosine` | alto | 0.4640 / 0.5237 / 0.5009 / 0.4820 / 0.5123 / 0.4261 / 0.5142 / 0.5095 | 0.491596 | 0.032763 | 0.030647 | 0.097538 | 0.4261 | 0.5237 |
| MCF7 | `mse_capped_norm` | bajo | 0.9332 / 0.8819 / 0.8747 / 0.8830 / 0.8998 / 0.9607 / 0.8771 / 0.8722 | 0.897818 | 0.032309 | 0.030222 | 0.088506 | 0.8722 | 0.9607 |
| MCF7 | `lfc_nmae` | bajo | 0.8688 / 0.8864 / 0.8763 / 0.8615 / 0.8751 / 0.8745 / 0.8771 / 0.8638 | 0.872926 | 0.008004 | 0.007487 | 0.024910 | 0.8615 | 0.8864 |
| MCF7 | `fid` | alto | 0.0026 / 0.0020 / 0.0022 / 0.0001 / 0.0014 / 0.0017 / 0.0022 / 0.0006 | 0.001586 | 0.000858 | 0.000803 | 0.002493 | 0.0001 | 0.0026 |
| MCF7 | `reach` | alto | 0.1280 / 0.1049 / 0.0778 / 0.1136 / 0.1418 / 0.1130 / 0.1350 / 0.1036 | 0.114708 | 0.020374 | 0.019058 | 0.063961 | 0.0778 | 0.1418 |
| MCF7 | `jac` | alto | 0.0000 / 0.0000 / 0.0000 / 0.0000 / 0.0005 / 0.0001 / 0.0001 / 0.0001 | 0.000085 | 0.000152 | 0.000142 | 0.000452 | 0.0000 | 0.0005 |

### 3.1 Donde cae el `b` de semilla 0 -- el que esta en los bundles

`z` son desviaciones desde la media de las 8; las columnas de conteo dicen
cuantas de las 8 quedan por debajo, iguales y por encima de `sem0`. Se separan
tres situaciones distintas que el percentil solo confundiria:

- **fuera**: `sem0` es estrictamente mas extremo que las ocho. El bundle se
  construyo con una realizacion que ningun otro sorteo reprodujo.
- **empate**: `sem0` coincide con el valor de varias de las 8. Suele ser un `b`
  exactamente 0, no una punta.
- **dentro**: hay valores a los dos lados.

🔴 **Y una columna mas, porque sin ella la `z` enganaria.** Donde las 8 semillas
toman muy pocos valores distintos, la `sd` es casi 0 y la `z` se dispara sin que
eso signifique una cola: significa que `sem0` cayo en **otro valor discreto**. La
columna `valores distintos de las 8` lo marca, y en esas celdas la `z` no se debe
leer como desviaciones de una nube.

| linea | metrica | b de semilla 0 | media de las 8 | sd de las 8 | valores distintos de las 8 | z de sem0 | debajo/igual/encima | sd del pool de 9 (marcado) | situacion de sem0 | sd degenerada |
|---|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | 0.486742 | 0.481179 | 0.031003 | 8 | +0.18 | 4/0/4 | 0.029060 | dentro |  |
| A549 | `mse_capped_norm` | 0.959928 | 1.034074 | 0.061015 | 8 | -1.22 | 1/0/7 | 0.062196 | dentro |  |
| A549 | `lfc_nmae` | 0.856686 | 0.862296 | 0.019665 | 8 | -0.29 | 5/0/3 | 0.018490 | dentro |  |
| A549 | `fid` | 0.000000 | 0.013158 | 0.037216 | 2 | -0.35 | 0/7/1 | 0.035088 | empate | 🔴 **si** |
| A549 | `reach` | 0.314403 | 0.412372 | 0.082408 | 8 | -1.19 | 1/0/7 | 0.083718 | dentro |  |
| A549 | `jac` | 0.454545 | 0.450758 | 0.010714 | 2 | +0.35 | 1/7/0 | 0.010101 | empate | 🔴 **si** |
| BXPC3 | `pds_cosine` | 0.476326 | 0.522609 | 0.022594 | 8 | -2.05 | 0/0/8 | 0.026167 | 🔴 **fuera** |  |
| BXPC3 | `mse_capped_norm` | 0.824367 | 0.789725 | 0.039148 | 8 | +0.88 | 6/0/2 | 0.038397 | dentro |  |
| BXPC3 | `lfc_nmae` | 0.814780 | 0.813297 | 0.004915 | 8 | +0.30 | 5/0/3 | 0.004624 | dentro |  |
| BXPC3 | `fid` | 0.001444 | 0.001392 | 0.002157 | 8 | +0.02 | 6/0/2 | 0.002018 | dentro |  |
| BXPC3 | `reach` | 0.344686 | 0.270184 | 0.061807 | 8 | +1.21 | 7/0/1 | 0.062923 | dentro |  |
| BXPC3 | `jac` | 0.000978 | 0.000139 | 0.000188 | 5 | +4.46 | 8/0/0 | 0.000330 | 🔴 **fuera** |  |
| HAP1 | `pds_cosine` | 0.437500 | 0.501657 | 0.034716 | 7 | -1.85 | 0/0/8 | 0.038883 | 🔴 **fuera** |  |
| HAP1 | `mse_capped_norm` | 1.003465 | 0.945743 | 0.064897 | 8 | +0.89 | 7/0/1 | 0.063682 | dentro |  |
| HAP1 | `lfc_nmae` | 0.897264 | 0.903099 | 0.014292 | 8 | -0.41 | 4/0/4 | 0.013509 | dentro |  |
| HAP1 | `fid` | 0.000000 | 0.013111 | 0.022120 | 4 | -0.59 | 0/5/3 | 0.021148 | empate |  |
| HAP1 | `reach` | 0.354552 | 0.286822 | 0.068276 | 8 | +0.99 | 7/0/1 | 0.067739 | dentro |  |
| HAP1 | `jac` | 0.424242 | 0.416667 | 0.014028 | 2 | +0.54 | 2/6/0 | 0.013362 | empate | 🔴 **si** |
| HT29 | `pds_cosine` | 0.501894 | 0.510653 | 0.034730 | 8 | -0.25 | 3/0/5 | 0.032618 | dentro |  |
| HT29 | `mse_capped_norm` | 1.137395 | 1.071865 | 0.035443 | 8 | +1.85 | 8/0/0 | 0.039703 | 🔴 **fuera** |  |
| HT29 | `lfc_nmae` | 0.864894 | 0.880140 | 0.024133 | 8 | -0.63 | 3/0/5 | 0.023139 | dentro |  |
| HT29 | `fid` | 0.114189 | 0.065012 | 0.040029 | 8 | +1.23 | 6/0/2 | 0.040875 | dentro |  |
| HT29 | `reach` | 0.327624 | 0.315222 | 0.030139 | 8 | +0.41 | 4/0/4 | 0.028494 | dentro |  |
| HT29 | `jac` | 0.091090 | 0.147889 | 0.034155 | 8 | -1.66 | 0/0/8 | 0.037138 | 🔴 **fuera** |  |
| K562 | `pds_cosine` | 0.520833 | 0.497869 | 0.033133 | 8 | +0.69 | 6/0/2 | 0.031925 | dentro |  |
| K562 | `mse_capped_norm` | 1.173286 | 1.155684 | 0.042256 | 8 | +0.42 | 5/0/3 | 0.039960 | dentro |  |
| K562 | `lfc_nmae` | 0.846558 | 0.853596 | 0.020413 | 8 | -0.34 | 1/0/7 | 0.019238 | dentro |  |
| K562 | `fid` | 0.076923 | 0.000078 | 0.000221 | 2 | +347.76 | 8/0/0 | 0.025616 | 🔴 **fuera** | 🔴 **si** |
| K562 | `reach` | 0.406814 | 0.476232 | 0.065801 | 8 | -1.05 | 1/0/7 | 0.065757 | dentro |  |
| K562 | `jac` | 0.212121 | 0.242483 | 0.000167 | 2 | -181.37 | 0/0/8 | 0.010122 | 🔴 **fuera** | 🔴 **si** |
| MCF7 | `pds_cosine` | 0.503788 | 0.491596 | 0.032763 | 8 | +0.37 | 4/0/4 | 0.030915 | dentro |  |
| MCF7 | `mse_capped_norm` | 0.922312 | 0.897818 | 0.032309 | 8 | +0.76 | 6/0/2 | 0.031305 | dentro |  |
| MCF7 | `lfc_nmae` | 0.872196 | 0.872926 | 0.008004 | 8 | -0.09 | 3/0/5 | 0.007491 | dentro |  |
| MCF7 | `fid` | 0.000173 | 0.001586 | 0.000858 | 8 | -1.65 | 1/0/7 | 0.000931 | dentro |  |
| MCF7 | `reach` | 0.101892 | 0.114708 | 0.020374 | 8 | -0.63 | 1/0/7 | 0.019531 | dentro |  |
| MCF7 | `jac` | 0.000029 | 0.000085 | 0.000152 | 6 | -0.37 | 3/1/4 | 0.000143 | dentro |  |

**7 de 36 celdas tienen el `b` del bundle estrictamente fuera** del
rango de las 8 semillas: el bundle se construyo con una realizacion que ningun
otro sorteo reprodujo.

En las 31 celdas con `sd` no degenerada, `|z|` mediano 0.69 y maximo 4.46 -- una dispersion normal. Las 5 degeneradas son otra cosa y van aparte: K562/`fid` (`z` = +348, 2 valor(es) distinto(s) en las 8), K562/`jac` (`z` = +181, 2 valor(es) distinto(s) en las 8), HAP1/`jac` (`z` = +1, 2 valor(es) distinto(s) en las 8), A549/`fid` (`z` = +0, 2 valor(es) distinto(s) en las 8). Ahi las 8 coinciden en un valor y `sem0` cayo en otro: la `z` es grande
porque el denominador es ~0, no porque `sem0` este lejos en unidades de la
metrica. §6 lo explica: son saltos de la convencion `1/33`.

La columna `sd del pool de 9` cumple con que el `b` de semilla 0 **entre como una
observacion mas**: es la `sd` de las 8 nuevas junto con esa. Esta al lado, no en
lugar de la de 8, y no se usa para el margen: para el margen entra la `sd` de las
8, porque `sem0` es el valor cuyo error se esta describiendo y no puede ayudar a
estimar su propia desviacion. 🔴 **Ningun bundle se toco y ningun `b` se
sustituyo por una media.**

### 3.2 Cuanto vale esa punta, en `s` -- DIAGNOSTICO, no cifra de competencia

🔴 **Marcado arriba y no en una nota: esta tabla NO es la puntuacion.** La `s` de
competencia es la de la columna `s con el b del bundle`, que sale de los 6 bundles
tal como estan y es la que `dispersion_u.md` midio. La segunda columna es un
CONTRAFACTUAL: que `s` habria salido si el baseline hubiera caido en la media de
sus 9 realizaciones en vez de en la que cayo. Sirve para una sola cosa -- poner en
unidades de `s` lo que §3.1 dice en unidades de `b` -- y para nada mas. **Ningun
bundle se modifico y ningun `b` se sustituyo por una media.**

Se calcula exacto, `s = (u - b)/(r - b)` con las dos `b`, no con la derivada de
§2: donde `|b - b0|` es grande respecto del span la aproximacion de primer orden
se rompe, y precisamente esas son las celdas que interesan. `u` es la media de las
8 emisiones de `dispersion_u.md` y `r` no se toca. El centro es la media de las **8
semillas nuevas**, sin `sem0`, por el mismo motivo que la `sd` de §3.1: `sem0` es
el valor cuyo desplazamiento se esta midiendo y no puede entrar en su propio
centro.

🔴 **Y el clamp, que si no se aplica da una columna falsa.** `catalog.py:415-416`
acota `expr_mse_unbiased_capped_norm` a `[0, 1]`, y solo esa. Las dos columnas de
`s` llevan el clamp en esas seis celdas, igual que lo lleva el marcador. La columna
`s medida por score` es la que escribio `score --real-bundle` en las 48 corridas de
`dispersion_u.md` (media de las 8), y esta al lado como comprobacion: la `s`
rehecha a mano tiene que reproducirla, y donde no lo hace es porque el clamp
intervino.

| linea | metrica | u medio | b del bundle (sem0) | media de las 8 | r | s medida por score (COMPETENCIA) | s rehecha con el b del bundle | \|rehecha - medida\| | s con la b media (DIAGNOSTICO) | cuanto se movia | signo |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | 0.5034 | 0.4867 | 0.4812 | 0.7778 | +0.0573 | +0.0573 | 2.8e-17 | +0.0750 | **+0.0177** |  |
| A549 | `mse_capped_norm` | 2.5432 | 0.9599 | 1.0341 | 0.1387 | +0.0000 | +0.0000 | 0.0e+00 | +0.0000 | **+0.0000** |  |
| A549 | `lfc_nmae` | 1.0152 | 0.8567 | 0.8623 | 0.4558 | -0.3954 | -0.3954 | 1.1e-16 | -0.3762 | **+0.0193** |  |
| A549 | `fid` | 0.0000 | 0.0000 | 0.0132 | 0.4179 | +0.0000 | +0.0000 | 0.0e+00 | -0.0325 | **-0.0325** |  |
| A549 | `reach` | 0.1694 | 0.3144 | 0.4124 | 0.7275 | -0.3510 | -0.3510 | 5.6e-17 | -0.7710 | **-0.4200** |  |
| A549 | `jac` | 0.4508 | 0.4545 | 0.4508 | 0.4733 | -0.2025 | -0.2025 | 0.0e+00 | +0.0000 | **+0.2025** |  |
| BXPC3 | `pds_cosine` | 0.5272 | 0.4763 | 0.5226 | 0.9551 | +0.1063 | +0.1063 | 1.4e-17 | +0.0107 | **-0.0956** |  |
| BXPC3 | `mse_capped_norm` | 1.3509 | 0.8244 | 0.7897 | 0.0687 | +0.0000 | +0.0000 | 0.0e+00 | +0.0000 | **+0.0000** |  |
| BXPC3 | `lfc_nmae` | 1.0085 | 0.8148 | 0.8133 | 0.5472 | -0.7241 | -0.7241 | 4.4e-16 | -0.7337 | **-0.0096** |  |
| BXPC3 | `fid` | 0.0008 | 0.0014 | 0.0014 | 0.5993 | -0.0011 | -0.0011 | 0.0e+00 | -0.0010 | **+0.0001** |  |
| BXPC3 | `reach` | 0.0786 | 0.3447 | 0.2702 | 0.8345 | -0.5432 | -0.5432 | 1.1e-16 | -0.3395 | **+0.2037** |  |
| BXPC3 | `jac` | 0.0000 | 0.0010 | 0.0001 | 0.2008 | -0.0049 | -0.0049 | 0.0e+00 | -0.0007 | **+0.0042** |  |
| HAP1 | `pds_cosine` | 0.5374 | 0.4375 | 0.5017 | 0.7896 | +0.2838 | +0.2838 | 1.7e-16 | +0.1242 | **-0.1596** |  |
| HAP1 | `mse_capped_norm` | 1.8341 | 1.0035 | 0.9457 | 0.1856 | +0.0000 | +0.0000 | 0.0e+00 | +0.0000 | **+0.0000** |  |
| HAP1 | `lfc_nmae` | 1.0024 | 0.8973 | 0.9031 | 0.3916 | -0.2079 | -0.2079 | 1.7e-16 | -0.1942 | **+0.0138** |  |
| HAP1 | `fid` | 0.0084 | 0.0000 | 0.0131 | 0.5468 | +0.0154 | +0.0154 | 1.7e-18 | -0.0087 | **-0.0242** | 🔴 **cambia de signo** |
| HAP1 | `reach` | 0.1579 | 0.3546 | 0.2868 | 0.7255 | -0.5301 | -0.5301 | 1.1e-16 | -0.2939 | **+0.2362** |  |
| HAP1 | `jac` | 0.4129 | 0.4242 | 0.4167 | 0.4682 | -0.2588 | -0.2588 | 8.3e-16 | -0.0736 | **+0.1852** |  |
| HT29 | `pds_cosine` | 0.5044 | 0.5019 | 0.5107 | 0.7835 | +0.0088 | +0.0088 | 5.4e-17 | -0.0230 | **-0.0318** | 🔴 **cambia de signo** |
| HT29 | `mse_capped_norm` | 1.9535 | 1.1374 | 1.0719 | 0.1328 | +0.0000 | +0.0000 | 0.0e+00 | +0.0000 | **+0.0000** |  |
| HT29 | `lfc_nmae` | 1.0165 | 0.8649 | 0.8801 | 0.4834 | -0.3973 | -0.3973 | 2.8e-16 | -0.3436 | **+0.0537** |  |
| HT29 | `fid` | 0.0042 | 0.1142 | 0.0650 | 0.4318 | -0.3463 | -0.3463 | 0.0e+00 | -0.1658 | **+0.1805** |  |
| HT29 | `reach` | 0.2076 | 0.3276 | 0.3152 | 0.6186 | -0.4126 | -0.4126 | 5.6e-17 | -0.3548 | **+0.0577** |  |
| HT29 | `jac` | 0.1780 | 0.0911 | 0.1479 | 0.2721 | +0.4803 | +0.4803 | 0.0e+00 | +0.2426 | **-0.2376** |  |
| K562 | `pds_cosine` | 0.5104 | 0.5208 | 0.4979 | 0.8227 | -0.0345 | -0.0345 | 2.6e-16 | +0.0386 | **+0.0731** | 🔴 **cambia de signo** |
| K562 | `mse_capped_norm` | 3.0524 | 1.1733 | 1.1557 | 0.1264 | +0.0000 | +0.0000 | 0.0e+00 | +0.0000 | **+0.0000** |  |
| K562 | `lfc_nmae` | 1.0143 | 0.8466 | 0.8536 | 0.4973 | -0.4803 | -0.4803 | 0.0e+00 | -0.4511 | **+0.0292** |  |
| K562 | `fid` | 0.0048 | 0.0769 | 0.0001 | 0.5117 | -0.1659 | -0.1659 | 0.0e+00 | +0.0092 | **+0.1751** | 🔴 **cambia de signo** |
| K562 | `reach` | 0.2313 | 0.4068 | 0.4762 | 0.8093 | -0.4361 | -0.4361 | 0.0e+00 | -0.7355 | **-0.2993** |  |
| K562 | `jac` | 0.2386 | 0.2121 | 0.2425 | 0.4900 | +0.0954 | +0.0954 | 1.4e-17 | -0.0155 | **-0.1110** | 🔴 **cambia de signo** |
| MCF7 | `pds_cosine` | 0.5264 | 0.5038 | 0.4916 | 0.9362 | +0.0523 | +0.0523 | 3.5e-17 | +0.0783 | **+0.0260** |  |
| MCF7 | `mse_capped_norm` | 1.7571 | 0.9223 | 0.8978 | 0.1050 | +0.0000 | +0.0000 | 0.0e+00 | +0.0000 | **+0.0000** |  |
| MCF7 | `lfc_nmae` | 1.0551 | 0.8722 | 0.8729 | 0.6638 | -0.8778 | -0.8778 | 2.2e-16 | -0.8712 | **+0.0066** |  |
| MCF7 | `fid` | 0.0000 | 0.0002 | 0.0016 | 0.4626 | -0.0003 | -0.0003 | 0.0e+00 | -0.0034 | **-0.0031** |  |
| MCF7 | `reach` | 0.0283 | 0.1019 | 0.1147 | 0.6764 | -0.1281 | -0.1281 | 2.8e-17 | -0.1539 | **-0.0257** |  |
| MCF7 | `jac` | 0.0000 | 0.0000 | 0.0001 | 0.0722 | -0.0004 | -0.0004 | 0.0e+00 | -0.0012 | **-0.0008** |  |

**Lo que esto mide.** En **15 de 36 celdas el sorteo del baseline
mueve la `s` mas que toda la dispersion de la emision junta** (`|cuanto se movia| > sd_s` de §4). Las cinco mayores: A549/`reach` (-0.4200 contra sd_s 0.1575), K562/`reach` (-0.2993 contra sd_s 0.1520), HT29/`jac` (-0.2376 contra sd_s 0.0592), HAP1/`reach` (+0.2362 contra sd_s 0.1456), BXPC3/`reach` (+0.2037 contra sd_s 0.0484).

🔴 **Y en 5 celdas la `s` cambia de SIGNO**: HAP1/`fid` (+0.0154 -> -0.0087), HT29/`pds_cosine` (+0.0088 -> -0.0230), K562/`pds_cosine` (-0.0345 -> +0.0386), K562/`fid` (-0.1659 -> +0.0092), K562/`jac` (+0.0954 -> -0.0155).

El cambio va en las dos direcciones, no en una: en 3 celdas (HAP1/`fid`, HT29/`pds_cosine`, K562/`jac`) el `b` del bundle le da a la emision un puntaje positivo que con el `b` medio
seria negativo, y en 2 (K562/`pds_cosine`, K562/`fid`) ocurre lo contrario. No hay sesgo de direccion: hay sensibilidad. **El signo
de esas celdas lo decide el sorteo del baseline, no la emision**, y en las dos
de `pds_cosine` la `s` de competencia esta a menos de 0.04 de cero, que es
menos que lo que el sorteo la mueve.

**El caso que ya estaba senalado y ahora tiene causa.** `cable_suelo.md`
marco HT29/`jac` como la celda donde una copia del control cobra `s` = +0.5012
por convencion, y `descomposicion.md` y `huecos.md` la persiguieron por el lado
de la `E`. Aqui se ve de donde sale la mitad: el `b` del bundle (0.0911) es **el minimo de las nueve realizaciones** (media 0.1479, maximo 0.1821), y con la media la misma emision cobraria +0.2426 en vez de +0.4803. **49%
de ese puntaje es la realizacion concreta del baseline**, no la emision ni la
convencion. Las dos explicaciones no compiten: la convencion pone el `u` arriba
y el sorteo del baseline pone el `b` abajo.

## 4. El margen completo, y cuanto crecia

Los tres terminos de §2 celda por celda, contra el margen de `dispersion_u.md`
que omitia el tercero. `sd_s` y `s` vienen de `dispersion_u.json` -- misma
poblacion, verificada en §1.

| linea | metrica | s medio | span r-b | \|ds/db\| | sd_s (termino u) | termino de r | sd_b | termino de b | margen de dispersion_u.md | margen COMPLETO | crece |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | 0.057336 | 0.291098 | 3.2383 | 0.284040 | 0.002558 | 0.031003 | 0.100396 | 0.284051 | **0.301271** | 1.06x |
| A549 | `mse_capped_norm` | 0.000000 | -0.821267 | 1.2176 | 0.000000 | 0.000000 | 0.061015 | 0.074293 | 0.000000 | **0.074293** | n/d |
| A549 | `lfc_nmae` | -0.395435 | -0.400857 | 3.4811 | 0.042320 | 0.010090 | 0.019665 | 0.068456 | 0.043506 | **0.081111** | 1.86x |
| A549 | `fid` | 0.000000 | 0.417863 | 2.3931 | 0.000000 | 0.000000 | 0.037216 | 0.089063 | 0.000000 | **0.089063** | n/d |
| A549 | `reach` | -0.351003 | 0.413137 | 3.2701 | 0.157531 | 0.031443 | 0.082408 | 0.269484 | 0.160639 | **0.313729** | 1.95x |
| A549 | `jac` | -0.202468 | 0.018709 | 64.2738 | 0.572666 | 0.169431 | 0.010714 | 0.688613 | 0.597205 | **0.911505** | 1.53x |
| BXPC3 | `pds_cosine` | 0.106309 | 0.478788 | 1.8666 | 0.089199 | 0.001671 | 0.022594 | 0.042173 | 0.089215 | **0.098681** | 1.11x |
| BXPC3 | `mse_capped_norm` | 0.000000 | -0.755633 | 1.3234 | 0.000000 | 0.000000 | 0.039148 | 0.051808 | 0.000000 | **0.051808** | n/d |
| BXPC3 | `lfc_nmae` | -0.724108 | -0.267589 | 6.4431 | 0.024938 | 0.061566 | 0.004915 | 0.031665 | 0.066425 | **0.073586** | 1.11x |
| BXPC3 | `fid` | -0.001123 | 0.597809 | 1.6747 | 0.002910 | 0.000038 | 0.002157 | 0.003613 | 0.002910 | **0.004639** | 1.59x |
| BXPC3 | `reach` | -0.543194 | 0.489816 | 3.1506 | 0.048389 | 0.022811 | 0.061807 | 0.194727 | 0.053496 | **0.201941** | 3.77x |
| BXPC3 | `jac` | -0.004893 | 0.199792 | 5.0297 | 0.000000 | 0.000249 | 0.000188 | 0.000946 | 0.000249 | **0.000979** | 3.92x |
| HAP1 | `pds_cosine` | 0.283755 | 0.352083 | 2.0343 | 0.152589 | 0.005497 | 0.034716 | 0.070622 | 0.152688 | **0.168229** | 1.10x |
| HAP1 | `mse_capped_norm` | 0.000000 | -0.817830 | 1.2227 | 0.000000 | 0.000000 | 0.064897 | 0.079353 | 0.000000 | **0.079353** | n/d |
| HAP1 | `lfc_nmae` | -0.207934 | -0.505674 | 2.3888 | 0.063788 | 0.007478 | 0.014292 | 0.034139 | 0.064225 | **0.072735** | 1.13x |
| HAP1 | `fid` | 0.015439 | 0.546846 | 1.8004 | 0.032694 | 0.001191 | 0.022120 | 0.039826 | 0.032716 | **0.051541** | 1.58x |
| HAP1 | `reach` | -0.530135 | 0.370986 | 4.1245 | 0.145648 | 0.048585 | 0.068276 | 0.281605 | 0.153537 | **0.320741** | 2.09x |
| HAP1 | `jac` | -0.258771 | 0.043914 | 28.6645 | 0.357138 | 0.222582 | 0.014028 | 0.402094 | 0.420821 | **0.582039** | 1.38x |
| HT29 | `pds_cosine` | 0.008826 | 0.281629 | 3.5194 | 0.177889 | 0.000239 | 0.034730 | 0.122231 | 0.177889 | **0.215835** | 1.21x |
| HT29 | `mse_capped_norm` | 0.000000 | -1.004573 | 0.9954 | 0.000000 | 0.000000 | 0.035443 | 0.035282 | 0.000000 | **0.035282** | n/d |
| HT29 | `lfc_nmae` | -0.397264 | -0.381536 | 3.6622 | 0.079919 | 0.006113 | 0.024133 | 0.088379 | 0.080153 | **0.119312** | 1.49x |
| HT29 | `fid` | -0.346264 | 0.317604 | 4.2388 | 0.034050 | 0.031920 | 0.040029 | 0.169676 | 0.046672 | **0.175978** | 3.77x |
| HT29 | `reach` | -0.412553 | 0.290995 | 4.8542 | 0.161456 | 0.025976 | 0.030139 | 0.146300 | 0.163532 | **0.219423** | 1.34x |
| HT29 | `jac` | 0.480264 | 0.181026 | 2.8710 | 0.059183 | 0.069382 | 0.034155 | 0.098061 | 0.091195 | **0.133912** | 1.47x |
| K562 | `pds_cosine` | -0.034504 | 0.301894 | 3.4267 | 0.146438 | 0.000897 | 0.033133 | 0.113538 | 0.146440 | **0.185299** | 1.27x |
| K562 | `mse_capped_norm` | 0.000000 | -1.046839 | 0.9553 | 0.000000 | 0.000000 | 0.042256 | 0.040365 | 0.000000 | **0.040365** | n/d |
| K562 | `lfc_nmae` | -0.480347 | -0.349215 | 4.2391 | 0.046277 | 0.005261 | 0.020413 | 0.086530 | 0.046575 | **0.098269** | 2.11x |
| K562 | `fid` | -0.165865 | 0.434782 | 2.6815 | 0.031276 | 0.007639 | 0.000221 | 0.000593 | 0.032195 | **0.032201** | 1.00x |
| K562 | `reach` | -0.436141 | 0.402489 | 3.5682 | 0.152003 | 0.026649 | 0.065801 | 0.234787 | 0.154321 | **0.280962** | 1.82x |
| K562 | `jac` | 0.095421 | 0.277874 | 3.2554 | 0.038556 | 0.010426 | 0.000167 | 0.000545 | 0.039941 | **0.039944** | 1.00x |
| MCF7 | `pds_cosine` | 0.052289 | 0.432386 | 2.1918 | 0.126699 | 0.000884 | 0.032763 | 0.071810 | 0.126702 | **0.145637** | 1.15x |
| MCF7 | `mse_capped_norm` | 0.000000 | -0.817354 | 1.2235 | 0.000000 | 0.000000 | 0.032309 | 0.039528 | 0.000000 | **0.039528** | n/d |
| MCF7 | `lfc_nmae` | -0.877766 | -0.208372 | 9.0116 | 0.036243 | 0.021838 | 0.008004 | 0.072128 | 0.042314 | **0.083624** | 1.98x |
| MCF7 | `fid` | -0.000303 | 0.462432 | 2.1631 | 0.000126 | 0.000027 | 0.000858 | 0.001857 | 0.000129 | **0.001861** | 14.46x |
| MCF7 | `reach` | -0.128132 | 0.574503 | 1.9637 | 0.015129 | 0.006341 | 0.020374 | 0.040008 | 0.016404 | **0.043240** | 2.64x |
| MCF7 | `jac` | -0.000399 | 0.072214 | 13.8532 | 0.000000 | 0.000060 | 0.000152 | 0.002100 | 0.000060 | **0.002101** | 35.13x |

**El margen crece en 27 de 29 celdas comparables**, por factores de 1.00x a 35.13x. El maximo es MCF7/`jac`.

Las otras 7 no llevan factor porque su margen anterior era **exactamente 0** y no se puede dividir por el. Son las `sd_s = 0` de `dispersion_u.md` §5, y ahi el cambio es mayor que en ninguna otra: pasan de 0 a un margen real de 0.0353 a 0.0891. Una celda que parecia sin ruido tiene el ruido entero en `b`.

🔴 **Y en 19 celdas el termino de `b` es el MAYOR de los tres** -- mayor que la dispersion de la propia emision:

| linea | metrica | sd_s (u) | termino de r | termino de b | b / u |
|---|---|---|---|---|---|
| MCF7 | `jac` | 0.000000 | 0.000060 | **0.002100** | inf |
| BXPC3 | `jac` | 0.000000 | 0.000249 | **0.000946** | inf |
| MCF7 | `fid` | 0.000126 | 0.000027 | **0.001857** | 14.7x |
| HT29 | `fid` | 0.034050 | 0.031920 | **0.169676** | 5.0x |
| BXPC3 | `reach` | 0.048389 | 0.022811 | **0.194727** | 4.0x |
| MCF7 | `reach` | 0.015129 | 0.006341 | **0.040008** | 2.6x |
| MCF7 | `lfc_nmae` | 0.036243 | 0.021838 | **0.072128** | 2.0x |
| HAP1 | `reach` | 0.145648 | 0.048585 | **0.281605** | 1.9x |
| K562 | `lfc_nmae` | 0.046277 | 0.005261 | **0.086530** | 1.9x |
| A549 | `reach` | 0.157531 | 0.031443 | **0.269484** | 1.7x |
| HT29 | `jac` | 0.059183 | 0.069382 | **0.098061** | 1.7x |
| A549 | `lfc_nmae` | 0.042320 | 0.010090 | **0.068456** | 1.6x |

Es la consecuencia directa de §2: `|ds/db| = |s - 1|/|r - b|` y los `s`
medidos son ~0 o negativos, asi que `|s - 1| >= 1` en casi todas las celdas
mientras `|ds/du| = 1/|r - b|`. Con `sd_b` del orden de `sd_u`, el termino de
`b` sale mayor que el de `u`. **El baseline era la fuente de ruido mas grande
de la escala y estaba sin medir.**

## 5. La resolucion de cada celda

El tamano minimo de diferencia entre dos candidatos que la celda distingue. Aqui
`r` y `b` **se cancelan**, y eso cambia el resultado por completo:

```
s1 - s2 = (u1 - b)/(r - b) - (u2 - b)/(r - b) = (u1 - u2)/(r - b)
```

Los dos candidatos se puntuan contra el MISMO bundle, asi que `r` y `b` son
comunes: su ruido mueve las dos `s` juntas y no entra en la diferencia. Con `k`
replicas (semillas de emision) por candidato, promediadas:

```
SE( s1_media - s2_media ) = sd_s * sqrt(2/k)

resolucion_k = sd_s * sqrt(2/k)
```

Es el margen de la DIFERENCIA, que es la cantidad que §3.16 compara. No se
multiplica por ningun factor de confianza: el encargo pide la resolucion, no un
corte, y la regla de decision ya esta escrita.

🔴 **Y la diferencia con §4 que conviene no confundir.** El margen completo de §4
es para situar UNA `s` en la escala absoluta (por ejemplo, decir si un candidato
supera el ancla). La resolucion de aqui es para comparar DOS candidatos. La
segunda es mucho mas pequena que la primera porque no carga ni `r` ni `b`. Las dos
columnas estan para que se use la que corresponda a la pregunta.

| linea | metrica | sd_s | resolucion k=1 | resolucion k=2 | resolucion k=4 | resolucion k=8 | margen absoluto de §4 (otra pregunta) |
|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | 0.284040 | 0.401693 | 0.284040 | 0.200846 | 0.142020 | **0.301271** |
| A549 | `mse_capped_norm` | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.074293** |
| A549 | `lfc_nmae` | 0.042320 | 0.059849 | 0.042320 | 0.029925 | 0.021160 | **0.081111** |
| A549 | `fid` | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.089063** |
| A549 | `reach` | 0.157531 | 0.222783 | 0.157531 | 0.111391 | 0.078766 | **0.313729** |
| A549 | `jac` | 0.572666 | 0.809872 | 0.572666 | 0.404936 | 0.286333 | **0.911505** |
| BXPC3 | `pds_cosine` | 0.089199 | 0.126147 | 0.089199 | 0.063073 | 0.044600 | **0.098681** |
| BXPC3 | `mse_capped_norm` | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.051808** |
| BXPC3 | `lfc_nmae` | 0.024938 | 0.035267 | 0.024938 | 0.017634 | 0.012469 | **0.073586** |
| BXPC3 | `fid` | 0.002910 | 0.004116 | 0.002910 | 0.002058 | 0.001455 | **0.004639** |
| BXPC3 | `reach` | 0.048389 | 0.068432 | 0.048389 | 0.034216 | 0.024194 | **0.201941** |
| BXPC3 | `jac` | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.000979** |
| HAP1 | `pds_cosine` | 0.152589 | 0.215793 | 0.152589 | 0.107896 | 0.076294 | **0.168229** |
| HAP1 | `mse_capped_norm` | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.079353** |
| HAP1 | `lfc_nmae` | 0.063788 | 0.090210 | 0.063788 | 0.045105 | 0.031894 | **0.072735** |
| HAP1 | `fid` | 0.032694 | 0.046236 | 0.032694 | 0.023118 | 0.016347 | **0.051541** |
| HAP1 | `reach` | 0.145648 | 0.205977 | 0.145648 | 0.102988 | 0.072824 | **0.320741** |
| HAP1 | `jac` | 0.357138 | 0.505069 | 0.357138 | 0.252535 | 0.178569 | **0.582039** |
| HT29 | `pds_cosine` | 0.177889 | 0.251572 | 0.177889 | 0.125786 | 0.088944 | **0.215835** |
| HT29 | `mse_capped_norm` | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.035282** |
| HT29 | `lfc_nmae` | 0.079919 | 0.113023 | 0.079919 | 0.056511 | 0.039960 | **0.119312** |
| HT29 | `fid` | 0.034050 | 0.048153 | 0.034050 | 0.024077 | 0.017025 | **0.175978** |
| HT29 | `reach` | 0.161456 | 0.228333 | 0.161456 | 0.114166 | 0.080728 | **0.219423** |
| HT29 | `jac` | 0.059183 | 0.083698 | 0.059183 | 0.041849 | 0.029592 | **0.133912** |
| K562 | `pds_cosine` | 0.146438 | 0.207094 | 0.146438 | 0.103547 | 0.073219 | **0.185299** |
| K562 | `mse_capped_norm` | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.040365** |
| K562 | `lfc_nmae` | 0.046277 | 0.065446 | 0.046277 | 0.032723 | 0.023139 | **0.098269** |
| K562 | `fid` | 0.031276 | 0.044231 | 0.031276 | 0.022115 | 0.015638 | **0.032201** |
| K562 | `reach` | 0.152003 | 0.214965 | 0.152003 | 0.107482 | 0.076002 | **0.280962** |
| K562 | `jac` | 0.038556 | 0.054526 | 0.038556 | 0.027263 | 0.019278 | **0.039944** |
| MCF7 | `pds_cosine` | 0.126699 | 0.179179 | 0.126699 | 0.089590 | 0.063349 | **0.145637** |
| MCF7 | `mse_capped_norm` | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.039528** |
| MCF7 | `lfc_nmae` | 0.036243 | 0.051256 | 0.036243 | 0.025628 | 0.018122 | **0.083624** |
| MCF7 | `fid` | 0.000126 | 0.000178 | 0.000126 | 0.000089 | 0.000063 | **0.001861** |
| MCF7 | `reach` | 0.015129 | 0.021395 | 0.015129 | 0.010698 | 0.007564 | **0.043240** |
| MCF7 | `jac` | 0.000000 | 0.000000 | 0.000000 | 0.000000 | 0.000000 | **0.002101** |

**9 celdas dan resolucion 0 en todo `k`**, y eso NO quiere decir que
distingan diferencias arbitrariamente pequenas. Son exactamente las de
`dispersion_u.md` §5: las 6 de `mse_capped_norm`, donde el clamp a
`[0,1]` (`catalog.py:415-416`) devuelve 0 para los dos candidatos si los dos estan
por debajo de `b`, y la celda no distingue **nada**; y A549/`fid`, BXPC3/`jac`, MCF7/`jac`, donde la emision no mueve el `u` crudo (`sd_u = 0`) y el cero si es
estructural. Una resolucion de 0 por ceguera y una por estabilidad se escriben
igual y significan lo contrario; esta columna no las separa y §5 de
`dispersion_u.md` si.

🔴 **Y lo que `b` le hace a esas mismas celdas.** Su resolucion sigue siendo 0 --
comparar dos candidatos ahi no depende de `b`, que se cancela -- pero su margen
ABSOLUTO deja de ser 0: §4 les da 0.0353 a 0.0891. `dispersion_u.md` escribio 0
en las dos columnas porque solo tenia la primera. Con el tercer termino las dos
dejan de coincidir, y es justo lo que debe pasar: una celda puede ser ciega a la
diferencia entre candidatos y seguir estando mal situada en la escala absoluta.

## 6. La `E` del baseline, sus cohortes y si alguna `sd` es Bernoulli

`dispersion_u.md` §5.2 destapo que la `E` del SUELO si se mueve con su semilla, y
que un solo gen predicho sobre una referencia vacia cuesta `1/33` de `jac` crudo,
amplificado por `1/span`. La misma pregunta, ahora para el baseline.

**Como se mide, y la validacion del proxy.** `baseline_results.csv` da el valor
por perturbacion, asi que las cohortes de `fid`, `reach` y `lfc_nmae` son el
conteo de filas no-NaN, y la `E` se lee como `jac == 1.0`. Que `jac == 1.0`
equivalga a `|union| = 0` es un PROXY, y esta validado antes de usarse: contra los
conteos de union de `huecos_cohorte.json` -- que salieron de `de_real.parquet` y
`de_pred.parquet` con `--write-degenes` -- los dos conjuntos coinciden exactamente
en las seis lineas con semilla 0 (E = 15/0/14/3/7/0). Ademas cada semilla lleva su
guarda: union vacia implica `n_real = n_pred = 0`, luego `max(n_pred, n_conf) = 0`
y `fid` sale NaN (`direction.py:776` con `_safe_div`, 689-699); asi que toda
perturbacion con `jac == 1.0` tiene que tener `fid` NaN.

| linea | E por semilla (j=0..7) | valores distintos de E | CONJUNTOS distintos de union vacia | cohorte de fid por semilla | cohorte de reach | cohorte de lfc_nmae | violaciones de la guarda del proxy |
|---|---|---|---|---|---|---|---|
| A549 | 15/15/15/15/15/14/15/15 | 2 | 2 | 18/18/18/18/18/19/18/18 | 18/18/18/18/18/18/18/18 | 8/8/8/8/8/8/8/8 | **0** |
| BXPC3 | 0/0/0/0/0/0/0/0 | 1 | 1 | 33/33/33/33/33/33/33/33 | 33/33/33/33/33/33/33/33 | 23/23/23/23/23/23/23/23 | **0** |
| HAP1 | 14/13/14/14/14/14/14/13 | 2 | 3 | 19/20/19/19/19/19/19/20 | 19/19/19/19/19/19/19/19 | 8/8/8/8/8/8/8/8 | **0** |
| HT29 | 6/5/3/6/5/4/4/6 | 4 | 6 | 27/28/30/27/28/29/29/27 | 27/27/27/27/27/27/27/27 | 15/15/15/15/15/15/15/15 | **0** |
| K562 | 8/8/8/8/8/8/8/8 | 1 | 1 | 25/25/25/25/25/25/25/25 | 25/25/25/25/25/25/25/25 | 9/9/9/9/9/9/9/9 | **0** |
| MCF7 | 0/0/0/0/0/0/0/0 | 1 | 1 | 33/33/33/33/33/33/33/33 | 33/33/33/33/33/33/33/33 | 33/33/33/33/33/33/33/33 | **0** |

Violaciones de la guarda en las 48 corridas: **0**. El proxy se sostiene
en todas, no solo en la semilla 0.

### 6.1 Que perturbaciones mueve la semilla del baseline

`huecos.md` §2.5 nombro BRD9, IRF7 y RAP1GAP (HT29) y ZBP1 (K562) como las que el
`b` de semilla 0 saco de la union vacia. Con 8 semillas se ve si son siempre esas.
El conjunto de referencia es el de `dispersion_u.md`: la union vacia de NUESTRO
suelo, que es `33 - |cohorte del lado real|` y no depende de ninguna emision.

| linea | E del suelo (referencia) | perturbacion que el baseline saca de la union vacia | en cuantas de las 8 semillas nuevas | la saco tambien la semilla 0 del bundle | en cuales (j) |
|---|---|---|---|---|---|
| A549 | 15 | **STAT5A** | 1/8 | no | 5 |
| BXPC3 | 0 | -- ninguna -- | 0/8 | no | -- |
| HAP1 | 14 | **BRD9** | 1/8 | no | 1 |
| HAP1 | 14 | **MAP3K14** | 1/8 | no | 7 |
| HT29 | 6 | **RAP1GAP** | 3/8 | **si** | 1/2/6 |
| HT29 | 6 | **BRD9** | 2/8 | **si** | 2/5 |
| HT29 | 6 | **STAT3** | 2/8 | no | 4/5 |
| HT29 | 6 | **HERC5** | 1/8 | no | 6 |
| HT29 | 6 | **IRF7** | 1/8 | **si** | 2 |
| K562 | 8 | **ZBP1** | 0/8 | **si** | -- |
| MCF7 | 0 | -- ninguna -- | 0/8 | no | -- |

🔴 **La respuesta: no son siempre esas.** El encargo preguntaba si las cuatro de
`huecos.md` §2.5 son estables o cambian. Cambian, y de dos maneras distintas:

- **HT29.** La semilla 0 saco **las tres a la vez** (BRD9, IRF7 y RAP1GAP), que es
  `E = 3`, el **minimo** de las nueve realizaciones -- igualado por una sola de las
  8 nuevas, `j=2`, que saca exactamente esas mismas tres. Las otras siete sacan
  menos: `E` = 6/5/6/5/4/4/6. Entre las 8 esos nombres aparecen sueltos --
  RAP1GAP en 3/8, BRD9 en 2/8, IRF7 en 1/8 -- y aparecen ademas **dos que la
  semilla 0 no saco: STAT3 (2/8) y HERC5 (1/8)**. El
  conjunto de candidatas es mayor que lo que revelaba una sola realizacion; lo que
  el sorteo decide es cuales caen.
- **K562.** La semilla 0 saco ZBP1. **Ninguna de las 8 nuevas saca nada**, ZBP1
  incluido: las ocho dan `E = 8`, que es exactamente la `E` del suelo. Ese evento
  es 1 de 9 sorteos.
- **A549 y HAP1.** La semilla 0 no saco ninguna, pero las nuevas si: STAT5A (1/8)
  en A549, BRD9 y MAP3K14 (1/8 cada una) en HAP1.

**Y una coincidencia que vale anotar sin explicarla.** STAT5A es la unica que el
baseline saca en A549, y es la misma que el SUELO saco en su `j=5`
(`dispersion_u.md` §5.2). Las dos emisiones son distintas y los sorteos
independientes; que la marginal de A549 sea la misma para las dos apunta a una
propiedad de esa perturbacion en esa linea, no de la emision. Queda dicho como
observacion, no como mecanismo: no se midio.

**Lo que esto le hace a `huecos.md` §2.5.** Esa seccion midio bien lo que midio --
con el `b` del bundle, esas cuatro celdas son las que rompen `E_u = E_b` -- pero
la lectura de que BRD9/IRF7/RAP1GAP/ZBP1 fueran *las* perturbaciones donde el
baseline se despega no se sostiene con nueve realizaciones a la vista. Son las que
se despegaron en una.

### 6.2 La `sd` acumulada de `b` con k = 2..8: nube o Bernoulli

`dispersion_u.md` §4 encontro celdas donde la `sd` arranca en 0 y salta: ahi la
dispersion no es una nube que se estreche con mas draws, es un evento discreto. La
misma comprobacion sobre `b`, en el orden `j = 0..7` fijado antes de medir.

| linea | metrica | k=2 | k=3 | k=4 | k=5 | k=6 | k=7 | k=8 | arranca en 0 y salta |
|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | 0.024106 | 0.018685 | 0.041422 | 0.035883 | 0.034770 | 0.032898 | 0.031003 | no |
| A549 | `mse_capped_norm` | 0.056496 | 0.042361 | 0.070327 | 0.063335 | 0.056649 | 0.061523 | 0.061015 | no |
| A549 | `lfc_nmae` | 0.029148 | 0.028468 | 0.024078 | 0.021656 | 0.022014 | 0.020450 | 0.019665 | no |
| A549 | `fid` | 0 | 0 | 0 | 0 | 0.042974 | 0.039786 | 0.037216 | 🔴 **si** |
| A549 | `reach` | 0.057062 | 0.087857 | 0.073883 | 0.065264 | 0.062047 | 0.084211 | 0.082408 | no |
| A549 | `jac` | 0 | 0 | 0 | 0.000000 | 0.012371 | 0.011453 | 0.010714 | 🔴 **si** |
| BXPC3 | `pds_cosine` | 0.013392 | 0.014883 | 0.022903 | 0.024448 | 0.023462 | 0.021458 | 0.022594 | no |
| BXPC3 | `mse_capped_norm` | 0.001468 | 0.014459 | 0.021492 | 0.032593 | 0.040991 | 0.037701 | 0.039148 | no |
| BXPC3 | `lfc_nmae` | 0.000157 | 0.001206 | 0.003079 | 0.006061 | 0.005575 | 0.005157 | 0.004915 | no |
| BXPC3 | `fid` | 0.000746 | 0.003457 | 0.003080 | 0.002801 | 0.002546 | 0.002330 | 0.002157 | no |
| BXPC3 | `reach` | 0.005169 | 0.004540 | 0.050249 | 0.051440 | 0.053423 | 0.055593 | 0.061807 | no |
| BXPC3 | `jac` | 0.000182 | 0.000157 | 0.000138 | 0.000124 | 0.000184 | 0.000174 | 0.000188 | no |
| HAP1 | `pds_cosine` | 0.022097 | 0.017054 | 0.030650 | 0.043123 | 0.039978 | 0.037360 | 0.034716 | no |
| HAP1 | `mse_capped_norm` | 0.084513 | 0.089868 | 0.077127 | 0.071961 | 0.065471 | 0.060119 | 0.064897 | no |
| HAP1 | `lfc_nmae` | 0.016207 | 0.014404 | 0.011835 | 0.012209 | 0.011532 | 0.012440 | 0.014292 | no |
| HAP1 | `fid` | 0 | 0 | 0.023684 | 0.021184 | 0.018964 | 0.017654 | 0.022120 | 🔴 **si** |
| HAP1 | `reach` | 0.014592 | 0.057933 | 0.057364 | 0.052156 | 0.080751 | 0.073746 | 0.068276 | no |
| HAP1 | `jac` | 0.021427 | 0.017495 | 0.015152 | 0.013552 | 0.012371 | 0.011453 | 0.014028 | no |
| HT29 | `pds_cosine` | 0.077005 | 0.059864 | 0.050010 | 0.043527 | 0.038965 | 0.036470 | 0.034730 | no |
| HT29 | `mse_capped_norm` | 0.074641 | 0.053606 | 0.048142 | 0.044244 | 0.039593 | 0.037574 | 0.035443 | no |
| HT29 | `lfc_nmae` | 0.005252 | 0.033160 | 0.030865 | 0.030293 | 0.027103 | 0.025827 | 0.024133 | no |
| HT29 | `fid` | 0.035954 | 0.048258 | 0.040270 | 0.039223 | 0.035088 | 0.037594 | 0.040029 | no |
| HT29 | `reach` | 0.016759 | 0.012911 | 0.017054 | 0.037442 | 0.034735 | 0.032086 | 0.030139 | no |
| HT29 | `jac` | 0.021393 | 0.046240 | 0.042855 | 0.037113 | 0.035453 | 0.033748 | 0.034155 | no |
| K562 | `pds_cosine` | 0.021427 | 0.022329 | 0.027330 | 0.033731 | 0.037683 | 0.035605 | 0.033133 | no |
| K562 | `mse_capped_norm` | 0.071269 | 0.056390 | 0.046061 | 0.039911 | 0.042938 | 0.045304 | 0.042256 | no |
| K562 | `lfc_nmae` | 0.019217 | 0.014288 | 0.013727 | 0.011897 | 0.011485 | 0.021924 | 0.020413 | no |
| K562 | `fid` | 0 | 0 | 0 | 0 | 0 | 0.000236 | 0.000221 | 🔴 **si** |
| K562 | `reach` | 0.090098 | 0.064372 | 0.083031 | 0.074372 | 0.075823 | 0.071031 | 0.065801 | no |
| K562 | `jac` | 0 | 0 | 0 | 0 | 0 | 0.000179 | 0.000167 | 🔴 **si** |
| MCF7 | `pds_cosine` | 0.042185 | 0.030110 | 0.025590 | 0.023840 | 0.035805 | 0.034517 | 0.032763 | no |
| MCF7 | `mse_capped_norm` | 0.036246 | 0.031883 | 0.026907 | 0.023490 | 0.034245 | 0.033063 | 0.032309 | no |
| MCF7 | `lfc_nmae` | 0.012450 | 0.008836 | 0.010650 | 0.009261 | 0.008291 | 0.007674 | 0.008004 | no |
| MCF7 | `fid` | 0.000456 | 0.000327 | 0.001102 | 0.000970 | 0.000867 | 0.000815 | 0.000858 | no |
| MCF7 | `reach` | 0.016330 | 0.025125 | 0.021122 | 0.024274 | 0.021711 | 0.021464 | 0.020374 | no |
| MCF7 | `jac` | 0.000020 | 0.000017 | 0.000014 | 0.000199 | 0.000179 | 0.000163 | 0.000152 | no |

**5 celdas arrancan en 0 y saltan**: A549/`fid`, A549/`jac`, HAP1/`fid`, K562/`fid`, K562/`jac`. Ahi la `sd_b` de 8 draws estima la probabilidad de un salto, no la
anchura de una nube.

Donde esto se tiene que leer con cuidado es en el **margen de §4**, que es el
que carga `sd_b` -- no en la resolucion de §5, que solo carga `sd_s` porque `b`
se cancela entre dos candidatos. Y las dos de K562 son el ejemplo mas claro del
problema: su `sd_b` sale ~0.0002 porque las 8 semillas coinciden, asi que su
termino de `b` en §4 es despreciable (0.0006 y 0.0005) -- y sin embargo §3.2
mide que el `b` de semilla 0 mueve su `s` en 0.1751 y 0.1110, con cambio de
signo en las dos. **El margen de primer orden no ve ese salto**: lo estima con
una `sd` que describe ocho draws de acuerdo, mientras el noveno cayo en otro
valor discreto. En esas celdas el contrafactual exacto de §3.2 es la cifra
informativa y el termino de `b` de §4 no lo es.

## 7. Las 20 filas trabajadas

**Criterio de seleccion, declarado antes de mirar los valores**: las celdas de
`pds_cosine`, `lfc_nmae` y `reach` en las seis lineas = 18 filas, mas las dos de
`jac` de HT29 y K562 -- las dos lineas donde `huecos.md` §2.5 habia identificado,
ANTES de este encargo, que la `E` del `b` no coincide con la del suelo
(BRD9/IRF7/RAP1GAP y ZBP1) = **20**. Es el mismo criterio de `dispersion_u.md`
§7, para que las filas se puedan poner lado a lado.

Cada fila trae los 8 `b` crudos, el `b` del bundle, y los tres terminos del margen
con su aritmetica a la vista.

| linea | metrica | los 8 b crudos | media de b | sd_b | b del bundle (sem0) | z de sem0 | r | span r-b | s medio | \|ds/db\| | sd_s | termino de r | termino de b | margen COMPLETO | resolucion k=8 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A549 | `pds_cosine` | 0.4924 / 0.5265 / 0.4962 / 0.4280 / 0.4839 / 0.4527 / 0.5028 / 0.4669 | 0.481179 | 0.031003 | 0.486742 | +0.18 | 0.777841 | 0.291098 | 0.057336 | 3.2383 | 0.284040 | 0.002558 | 0.100396 | **0.301271** | 0.142020 |
| A549 | `lfc_nmae` | 0.8938 / 0.8526 / 0.8392 / 0.8493 / 0.8718 / 0.8870 / 0.8556 / 0.8491 | 0.862296 | 0.019665 | 0.856686 | -0.29 | 0.455829 | -0.400857 | -0.395435 | 3.4811 | 0.042320 | 0.010090 | 0.068456 | **0.081111** | 0.021160 |
| A549 | `reach` | 0.3515 / 0.2708 / 0.4463 / 0.3916 / 0.3938 / 0.4223 / 0.5443 / 0.4784 | 0.412372 | 0.082408 | 0.314403 | -1.19 | 0.727540 | 0.413137 | -0.351003 | 3.2701 | 0.157531 | 0.031443 | 0.269484 | **0.313729** | 0.078766 |
| BXPC3 | `pds_cosine` | 0.5350 / 0.5540 / 0.5246 / 0.4991 / 0.4962 / 0.5009 / 0.5218 / 0.5492 | 0.522609 | 0.022594 | 0.476326 | -2.05 | 0.955114 | 0.478788 | 0.106309 | 1.8666 | 0.089199 | 0.001671 | 0.042173 | **0.098681** | 0.044600 |
| BXPC3 | `lfc_nmae` | 0.8095 / 0.8092 / 0.8114 / 0.8159 / 0.8237 / 0.8107 / 0.8156 / 0.8104 | 0.813297 | 0.004915 | 0.814780 | +0.30 | 0.547191 | -0.267589 | -0.724108 | 6.4431 | 0.024938 | 0.061566 | 0.031665 | **0.073586** | 0.012469 |
| BXPC3 | `reach` | 0.2721 / 0.2794 / 0.2711 / 0.1740 / 0.1878 / 0.3034 / 0.3186 / 0.3549 | 0.270184 | 0.061807 | 0.344686 | +1.21 | 0.834502 | 0.489816 | -0.543194 | 3.1506 | 0.048389 | 0.022811 | 0.194727 | **0.201941** | 0.024194 |
| HAP1 | `pds_cosine` | 0.5445 / 0.5133 / 0.5407 / 0.4782 / 0.4432 / 0.4782 / 0.5208 / 0.4943 | 0.501657 | 0.034716 | 0.437500 | -1.85 | 0.789583 | 0.352083 | 0.283755 | 2.0343 | 0.152589 | 0.005497 | 0.070622 | **0.168229** | 0.076294 |
| HAP1 | `lfc_nmae` | 0.8904 / 0.9133 / 0.8867 / 0.8942 / 0.9110 / 0.8900 / 0.9151 / 0.9240 | 0.903099 | 0.014292 | 0.897264 | -0.41 | 0.391590 | -0.505674 | -0.207934 | 2.3888 | 0.063788 | 0.007478 | 0.034139 | **0.072735** | 0.031894 |
| HAP1 | `reach` | 0.2947 / 0.2741 / 0.1857 / 0.3164 / 0.2322 / 0.4221 / 0.2819 / 0.2875 | 0.286822 | 0.068276 | 0.354552 | +0.99 | 0.725538 | 0.370986 | -0.530135 | 4.1245 | 0.145648 | 0.048585 | 0.281605 | **0.320741** | 0.072824 |
| HT29 | `pds_cosine` | 0.5540 / 0.4451 / 0.5426 / 0.5350 / 0.5095 / 0.5133 / 0.4953 / 0.4905 | 0.510653 | 0.034730 | 0.501894 | -0.25 | 0.783523 | 0.281629 | 0.008826 | 3.5194 | 0.177889 | 0.000239 | 0.122231 | **0.215835** | 0.088944 |
| HT29 | `lfc_nmae` | 0.9106 / 0.9180 / 0.8572 / 0.8656 / 0.8560 / 0.8832 / 0.8622 / 0.8882 | 0.880140 | 0.024133 | 0.864894 | -0.63 | 0.483357 | -0.381536 | -0.397264 | 3.6622 | 0.079919 | 0.006113 | 0.088379 | **0.119312** | 0.039960 |
| HT29 | `reach` | 0.3187 / 0.3424 / 0.3394 / 0.3067 / 0.2499 / 0.3340 / 0.3281 / 0.3026 | 0.315222 | 0.030139 | 0.327624 | +0.41 | 0.618619 | 0.290995 | -0.412553 | 4.8542 | 0.161456 | 0.025976 | 0.146300 | **0.219423** | 0.080728 |
| K562 | `pds_cosine` | 0.4612 / 0.4915 / 0.5047 / 0.4451 / 0.5294 / 0.5417 / 0.5199 / 0.4896 | 0.497869 | 0.033133 | 0.520833 | +0.69 | 0.822727 | 0.301894 | -0.034504 | 3.4267 | 0.146438 | 0.000897 | 0.113538 | **0.185299** | 0.073219 |
| K562 | `lfc_nmae` | 0.8511 / 0.8783 / 0.8571 / 0.8477 / 0.8575 / 0.8689 / 0.8092 / 0.8590 | 0.853596 | 0.020413 | 0.846558 | -0.34 | 0.497343 | -0.349215 | -0.480347 | 4.2391 | 0.046277 | 0.005261 | 0.086530 | **0.098269** | 0.023139 |
| K562 | `reach` | 0.5890 / 0.4616 / 0.5093 / 0.3914 / 0.5303 / 0.4072 / 0.4393 / 0.4818 | 0.476232 | 0.065801 | 0.406814 | -1.05 | 0.809303 | 0.402489 | -0.436141 | 3.5682 | 0.152003 | 0.026649 | 0.234787 | **0.280962** | 0.076002 |
| MCF7 | `pds_cosine` | 0.4640 / 0.5237 / 0.5009 / 0.4820 / 0.5123 / 0.4261 / 0.5142 / 0.5095 | 0.491596 | 0.032763 | 0.503788 | +0.37 | 0.936174 | 0.432386 | 0.052289 | 2.1918 | 0.126699 | 0.000884 | 0.071810 | **0.145637** | 0.063349 |
| MCF7 | `lfc_nmae` | 0.8688 / 0.8864 / 0.8763 / 0.8615 / 0.8751 / 0.8745 / 0.8771 / 0.8638 | 0.872926 | 0.008004 | 0.872196 | -0.09 | 0.663824 | -0.208372 | -0.877766 | 9.0116 | 0.036243 | 0.021838 | 0.072128 | **0.083624** | 0.018122 |
| MCF7 | `reach` | 0.1280 / 0.1049 / 0.0778 / 0.1136 / 0.1418 / 0.1130 / 0.1350 / 0.1036 | 0.114708 | 0.020374 | 0.101892 | -0.63 | 0.676395 | 0.574503 | -0.128132 | 1.9637 | 0.015129 | 0.006341 | 0.040008 | **0.043240** | 0.007564 |
| HT29 | `jac` | 0.1819 / 0.1517 / 0.0911 / 0.1821 / 0.1517 / 0.1212 / 0.1213 / 0.1820 | 0.147889 | 0.034155 | 0.091090 | -1.66 | 0.272116 | 0.181026 | 0.480264 | 2.8710 | 0.059183 | 0.069382 | 0.098061 | **0.133912** | 0.029592 |
| K562 | `jac` | 0.2424 / 0.2424 / 0.2424 / 0.2424 / 0.2424 / 0.2424 / 0.2429 / 0.2424 | 0.242483 | 0.000167 | 0.212121 | -181.37 | 0.489996 | 0.277874 | 0.095421 | 3.2554 | 0.038556 | 0.010426 | 0.000545 | **0.039944** | 0.019278 |

Son **20** filas.

## Lo que no se hizo

- **No se reconstruyo ningun bundle.** Los 6 de competencia de
  `/datos/modelo/out/nativas_comp/` estan exactamente como estaban.
- **No se sustituyo el `b` de los bundles por una media.** El `b` de semilla 0
  entra como una observacion mas y marcado (§3.1); el margen usa la `sd` de las 8.
- **No se construyo ningun candidato.** Aqui no hay emision de candidato: lo que
  se emite es el baseline del propio marcador.
- **No se propone ningun umbral.** §4 da el margen y §5 la resolucion con su
  formula; la regla de decision ya esta escrita en §3.16 y no se toca.
- No se leyo nada de `~/vcc/benchmark/` ni de
  `/datos/benchmark_VEDADO_techo_metricas/`.
- Sin envios, sin `vcc sample`, sin commits.

