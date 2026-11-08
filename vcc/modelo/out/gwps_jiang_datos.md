# gwps (K562, sin estimular) contra el delta especifico de Jiang (seis lineas, con IFNB): los datos

## 🔴 0. Que es este reporte y que no es

**Presenta datos.** Declara la construccion, da las tablas con sus controles, y
reporta los defectos de la propia construccion. **No lleva veredicto, ni seccion
de conclusiones, ni «lo que esto significa», ni nada dirigido a un candidato.**
Donde un numero sale sorprendente, esta el numero.

**🔴 DOS POBLACIONES DISTINTAS, y no son la misma cosa:**

| | poblacion | estimulo | perturbacion | linea |
|---|---|---|---|---|
| **gwps** | Replogle genome-wide perturb-seq | **SIN estimular** | knockdown **CRISPRi** | **K562** |
| **Jiang** | Jiang perturb-seq | **CON IFNB** | (segun el dataset) | **seis lineas** |

Lo que se mide es si una estructura medida en la primera se relaciona con la
segunda. **No se mide si son equivalentes, y nada de aqui supone que lo sean.**

Esto es transductivo como todo lo de esta serie: `E_g` de Jiang se construye con
la verdad de Jiang. Ninguna cifra de aqui entra a una submission.

---

## 1. Construccion, declarada antes de medir

### 1.1 El archivo de gwps y la columna que es el perfil

| | |
|---|---|
| archivo | `/datos/ref/replogle/k562_gwps_raw_bulk.h5ad` |
| matriz | **`X`**, forma **(11258 perturbaciones, 8248 genes)**, pseudobulk normalizado por perturbacion. **Es la unica matriz del archivo**: las claves de primer nivel son `X`, `obs`, `var`; no hay `layers` ni `obsm` |
| identidad de la perturbacion | `obs/gene_transcript`, formato `{idx}_{GEN}_{P1|P2|P1P2}_{ENSG}` |
| controles | `obs/core_control == True` -> **514 filas**, todas `non-targeting` (comprobado: 0 filas `core_control` que no sean `non-targeting`). Hay **71 filas `non-targeting` que NO son `core_control`** y **no se usan** |
| perturbaciones dirigidas | **10744 filas** (`core_control == False`) |

**El perfil de respuesta, declarado como formula:**

```
p_g = log2( (X_g + eps) / (media de X sobre las 514 filas core_control + eps) )
eps = 1e-8
```

Es log-fold-change contra el pseudobulk de los no-dirigidos. **No es la columna
`fold_expr` de `obs`**, que es un escalar por perturbacion (el knockdown del
propio gen diana) y no un vector sobre genes.

### 1.2 Como se alinean los universos de genes, y cuantos quedan

Por **simbolo**: `var/gene_name` de gwps (categorico, 8248 entradas) contra el
array `genes` de `rasgos_<linea>.npz`, que es el universo de Jiang de
`control_o_no_rasgos.py`.

🔴 **Dos nombres de gen estan duplicados en `var` de gwps (`HSPA14`, `TBCE`):
se excluyen ENTEROS**, no se elige una de las dos filas.

| linea | genes de Jiang | **genes en comun** | fraccion |
|---|---|---|---|
| A549 | 7169 | **4825** | 67.3 % |
| BXPC3 | 7673 | **5213** | 67.9 % |
| HAP1 | 8045 | **5776** | 71.8 % |
| HT29 | 7224 | **5186** | 71.8 % |
| K562 | 7616 | **5844** | 76.7 % |
| MCF7 | 7444 | **5266** | 70.7 % |

### 1.3 `E_g` de Jiang

`E1` y `E2` de `rasgos_<linea>.npz`, que son **exactamente** el delta ESPECIFICO
con el control **PARTIDO** en dos mitades disjuntas de `control_o_no_rasgos.py`.
No se recalcula nada: se leen y se restringen a los genes en comun.

### 1.4 La metrica: Spearman, con el producto cruzado de mitades

🔴 **Spearman y no Pearson**, porque lo que puede transferir es el ORDEN y nunca
la magnitud. Y con el producto cruzado de mitades de Jiang para que el ruido de
Jiang no entre, como `control_o_no.md` §2.2 -- pero **en espacio de rangos**:

```
p = rango(perfil gwps), centrado
a = rango(E1), centrado      b = rango(E2), centrado
G = <p, (a+b)/2> / ( ||p|| * sqrt(<a,b>) )
<a,b> <= 0  ->  n/d
```

Los rangos de `E1` y `E2` son funciones de mitades de control **disjuntas**, asi
que su ruido sigue siendo independiente y el producto cruzado sigue valiendo.
Empates promediados.

**En columnas aparte** (§3.3): el **Pearson** con el mismo producto cruzado, y el
**Spearman ingenuo** contra `(E1+E2)/2`, que esta atenuado por el ruido de Jiang.

### 1.5 Las cuatro cantidades

| | que es |
|---|---|
| **`G_mismo`** | perfil gwps de `g` contra `E_g` de Jiang -- lo que se mide |
| **`G_otro`** | perfil gwps de `g` contra `E_q`, `q != g` -- 🔴 **CONTROL NEGATIVO**. Se dan media, mediana, sd, min y max sobre los `q` |
| **`A_esp`** | el techo, **ya medido** en `control_o_no.md` §8.1, leido del npz. No se recalcula |
| **`G_prom`** | perfiles **promedio** de gwps contra `E_g`, para ver si lo que correlaciona es la respuesta generica de gwps y no la del blanco. **Dos versiones**: promedio de los perfiles de los blancos cubiertos, y promedio de las 10744 perturbaciones dirigidas |

### 1.6 Las dos construcciones, y el gen blanco

🔴 **Los blancos con P1 y P2 por separado se reportan POR SEPARADO y NO se
promedian.** `gwps_fiabilidad.md` midio que P1 domina en 67.6 % de 546 pares, asi
que la media no es neutral. Son `ID2` y `RAP1GAP` (§4).

**El gen blanco se reporta con y sin el** (§3.1 y §3.2), porque es el punto de
mayor palanca en los dos lados (`control_o_no.md` §4.3). Las tablas por fila son
**sin** el gen blanco.

### 1.7 Semillas

Familia nueva **`20270000 + 100*fold`**, `fold` = indice de linea (0-5). Gobierna
el **unico** uso de azar: el nulo por permutacion de las etiquetas de gen del
perfil de gwps, **200 permutaciones** por fila. `G_mismo`, `G_otro`, `A_esp` y
`G_prom` son deterministas dadas las mitades de control, que ya venian fijadas
en el npz por la semilla `20267999 + fold` del encargo anterior.

---

## 2. PARTE 1 -- Cobertura del panel en k562_gwps

**De los 33 blancos del panel comun de Jiang (`modelo/out/suelo_panel.json`),
aparecen como gen perturbado en k562_gwps: `N = 27`.**

### 2.1 Los 27 que SI aparecen, con su construccion

| # | blanco | construccion(es) | # | blanco | construccion(es) |
|---|---|---|---|---|---|
| 1 | `AKT1` | P1P2 | 15 | `MAPK8` | P1P2 |
| 2 | `BATF2` | P1P2 | 16 | `PARP12` | P1P2 |
| 3 | `BRD9` | P1P2 | 17 | `RAP1GAP` | **P1 y P2** |
| 4 | `CEBPG` | P1P2 | 18 | `RNF114` | P1P2 |
| 5 | `ELK1` | P1P2 | 19 | `SOCS1` | P1P2 |
| 6 | `ETS1` | P1P2 | 20 | `SP100` | P1P2 |
| 7 | `ETS2` | P1P2 | 21 | `STAT1` | P1P2 |
| 8 | `ID2` | **P1 y P2** | 22 | `STAT3` | P1P2 |
| 9 | `IFNAR1` | P1P2 | 23 | `STAT4` | P1P2 |
| 10 | `IFNAR2` | P1P2 | 24 | `STAT5A` | P1P2 |
| 11 | `IRF1` | P1P2 | 25 | `TRAFD1` | P1P2 |
| 12 | `IRF7` | P1P2 | 26 | `TYK2` | P1P2 |
| 13 | `JAK1` | P1P2 | 27 | `VAV1` | P1P2 |
| 14 | `JUN` | P1P2 | | | |

### 2.2 Los 6 que NO aparecen

`DTX3L`, `H1F0`, `HERC5`, `HERC6`, `MAP3K14`, `ZBP1`

### 2.3 Los que aparecen con mas de una construccion

**2 de los 27**: `ID2` (filas P1 y P2) y `RAP1GAP` (filas P1 y P2). Los otros
**25** tienen una sola fila, de construccion `P1P2`.

Reparto de construcciones en todo gwps, para contexto: `P1P2` 8866 filas, `P2`
745, `P1` 742, mas 585 `non-targeting` y 17 filas con sufijo `nan`.

**Asi que las filas medidas son 29** (27 blancos, con `ID2` y `RAP1GAP` dos veces)
x 6 lineas = 174, de las que 165 dan `G_mismo` y 9 salen `n/d` (§5.2).

---

## 3. PARTE 2 -- La relacion, con su control negativo

### 3.1 Por linea, sin el gen blanco

`G_mismo` y `G_otro` son **medias sobre las filas** de esa linea. `G_otro sd` es
la media de las sd de cada fila sobre sus `q`. `nulo sd` es la sd del nulo por
permutacion de etiquetas de gen (200 por fila).

| linea | genes comunes | filas | **`G_mismo` media** | mediana | **`G_otro` media** | `G_otro` sd | **`G_mismo` - `G_otro`** | `G_mismo` > 0 | nulo sd | **`G_prom` cubiertos** | **`G_prom` todas** |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A549 | 4825 | 28 | **+0.0303** | +0.0274 | **-0.0013** | 0.0740 | **+0.0316** | 21/28 | 0.0697 | -0.0437 | **+0.0655** |
| BXPC3 | 5213 | 29 | **+0.0022** | +0.0127 | **+0.0006** | 0.0484 | **+0.0016** | 17/29 | 0.0351 | -0.0078 | **+0.0734** |
| HAP1 | 5776 | 26 | **+0.0108** | +0.0094 | **-0.0032** | 0.0892 | **+0.0140** | 16/26 | 0.0724 | -0.0474 | **+0.0876** |
| HT29 | 5186 | 27 | **-0.0041** | -0.0145 | **-0.0019** | 0.0740 | **-0.0022** | 13/27 | 0.0532 | -0.0221 | **+0.0134** |
| **K562** 🔴 | 5844 | 27 | **+0.0200** | +0.0116 | **+0.0008** | 0.1423 | **+0.0191** | 15/27 | 0.0517 | -0.0130 | **+0.0246** |
| MCF7 | 5266 | 28 | **+0.0017** | -0.0064 | **-0.0015** | 0.0520 | **+0.0032** | 11/28 | 0.0435 | -0.0069 | **+0.0044** |

### 3.2 La misma tabla CON el gen blanco dentro

| linea | `G_mismo` media | `G_otro` media | `G_mismo` - `G_otro` | `G_mismo` > 0 |
|---|---|---|---|---|
| A549 | +0.0315 | -0.0016 | +0.0330 | 21/28 |
| BXPC3 | +0.0032 | +0.0003 | +0.0028 | 17/29 |
| HAP1 | +0.0123 | -0.0037 | +0.0160 | 16/26 |
| HT29 | -0.0025 | -0.0023 | -0.0002 | 13/27 |
| **K562** 🔴 | **-0.0148** | -0.0030 | **-0.0119** | 15/28 |
| MCF7 | +0.0030 | -0.0017 | +0.0047 | 11/28 |

**En K562 el gen blanco cambia el signo de la media**: `-0.0148` con el,
`+0.0200` sin el. En las otras cinco lineas la media se mueve menos de 0.002.

### 3.3 Spearman contra Pearson contra el Spearman ingenuo

| linea | **`G_mismo` Spearman, cruzado** (el declarado) | `G_mismo` Pearson, cruzado | Spearman ingenuo (contra `(E1+E2)/2`) |
|---|---|---|---|
| A549 | **+0.0303** | +0.0307 | +0.0056 |
| BXPC3 | **+0.0022** | +0.0020 | +0.0021 |
| HAP1 | **+0.0108** | +0.0081 | -0.0000 |
| HT29 | **-0.0041** | +0.0002 | -0.0031 |
| **K562** 🔴 | **+0.0200** | +0.0240 | +0.0057 |
| MCF7 | **+0.0017** | -0.0015 | -0.0006 |

### 3.4 Conteos contra los controles

| linea | `G_mismo` > su `G_otro` | \|`G_mismo`\| > 2 sd del nulo | filas con z > +2 | filas con z < -2 |
|---|---|---|---|---|
| A549 | 21/28 | 5/28 | 2 | 3 |
| BXPC3 | 16/29 | 7/29 | 4 | 3 |
| HAP1 | 14/26 | 1/26 | 0 | 1 |
| HT29 | 13/27 | 2/27 | 0 | 3 |
| **K562** 🔴 | 15/27 | 4/27 | 4 | 1 |
| MCF7 | 12/28 | 3/28 | 1 | 2 |

### 3.5 El techo `A_esp`, ya medido, para escala

De `control_o_no.md` §8.1, leido del npz, sobre las filas medidas:

| linea | `A_esp` min | mediana | max |
|---|---|---|---|
| A549 | -0.0493 | +0.0399 | +0.4998 |
| BXPC3 | +0.0347 | +0.1274 | +0.6180 |
| HAP1 | +0.0152 | +0.0603 | +0.4844 |
| HT29 | +0.0010 | +0.0585 | +0.5398 |
| K562 | +0.0023 | +0.0631 | +0.4216 |
| MCF7 | +0.0231 | +0.1253 | +0.4782 |

---

## 4. Las dos construcciones, por separado y sin promediar

`G_mismo` sin gen blanco. 🔴 Estas filas **no** se promedian entre si.

| blanco | constr. | A549 | BXPC3 | HAP1 | HT29 | **K562** 🔴 | MCF7 |
|---|---|---|---|---|---|---|---|
| `ID2` | **P1** | -0.0965 | -0.0104 | n/d | -0.0708 | -0.1100 | -0.0657 |
| `ID2` | **P2** | **+0.3156** | +0.0127 | n/d | +0.0077 | -0.0102 | +0.0339 |
| `RAP1GAP` | **P1** | +0.0600 | -0.0663 | +0.1967 | -0.0210 | +0.0460 | +0.0976 |
| `RAP1GAP` | **P2** | +0.0412 | +0.0445 | +0.2082 | -0.0749 | -0.0468 | -0.0108 |

Diferencia **P1 - P2** por linea:

| blanco | A549 | BXPC3 | HAP1 | HT29 | **K562** 🔴 | MCF7 |
|---|---|---|---|---|---|---|
| `ID2` | **-0.4121** | -0.0231 | n/d | -0.0785 | -0.0997 | -0.0996 |
| `RAP1GAP` | +0.0188 | -0.1108 | -0.0116 | +0.0540 | +0.0928 | +0.1084 |

`ID2` da P1 negativo y P2 positivo en 5 de las 5 lineas con dato; `RAP1GAP`
cambia de signo la diferencia segun la linea. El `+0.3156` de `ID2` P2 en A549 es
el `G_mismo` mayor de las 165 filas medidas, y su `A_esp` es **-0.0067** (§5.1).

---

## 5. PARTE 3 -- Por linea, con K562 senalada

### 5.1 Las 29 filas de K562, ordenadas por `G_mismo`, sin gen blanco

🔴 **K562 es la linea de gwps.** La tabla, sin comentario del patron.

| # | blanco | constr. | n genes | **`G_mismo`** | **`G_otro` media** | `G_otro` min | `G_otro` max | z vs nulo | `A_esp` | `G_prom` todas |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `STAT5A` | P1P2 | 5843 | **+0.1581** | -0.0251 | -0.8643 | +0.1134 | +4.02 | +0.1055 | -0.1052 |
| 2 | `RNF114` | P1P2 | 5843 | **+0.1407** | -0.0090 | -0.5452 | +0.1620 | +1.65 | +0.0342 | +0.1711 |
| 3 | `SOCS1` | P1P2 | 5844 | **+0.1374** | +0.0037 | -0.4135 | +0.1975 | +1.57 | +0.0224 | -0.0609 |
| 4 | `BRD9` | P1P2 | 5843 | **+0.1321** | -0.0406 | -1.0959 | +0.2649 | +4.68 | +0.1678 | +0.4030 |
| 5 | `CEBPG` | P1P2 | 5843 | **+0.1149** | +0.0108 | -0.2952 | +0.4035 | +3.92 | +0.1306 | +0.1983 |
| 6 | `ETS2` | P1P2 | 5843 | **+0.0992** | -0.0004 | -0.1956 | +0.1527 | +2.61 | +0.1332 | +0.0529 |
| 7 | `VAV1` | P1P2 | 5844 | **+0.0860** | +0.0134 | -0.1644 | +0.4116 | +1.04 | +0.0049 | +0.4677 |
| 8 | `RAP1GAP` | **P1** | 5844 | **+0.0460** | +0.0088 | -0.1779 | +0.1301 | +0.91 | +0.0264 | -0.0922 |
| 9 | `JAK1` | P1P2 | 5843 | **+0.0458** | -0.0480 | -0.6143 | +0.1423 | +1.48 | +0.2666 | -0.2466 |
| 10 | `JUN` | P1P2 | 5843 | **+0.0312** | -0.0093 | -0.5162 | +0.2248 | +0.55 | +0.0607 | -0.1247 |
| 11 | `STAT1` | P1P2 | 5843 | **+0.0300** | +0.0197 | -0.1299 | +0.1617 | +0.56 | +0.0440 | -0.2956 |
| 12 | `AKT1` | P1P2 | 5843 | **+0.0299** | -0.0043 | -0.5386 | +0.3980 | +0.55 | +0.0600 | +0.1298 |
| 13 | `TYK2` | P1P2 | 5843 | **+0.0240** | -0.0213 | -0.4133 | +0.2170 | +0.82 | +0.3058 | -0.3101 |
| 14 | `IFNAR1` | P1P2 | 5843 | **+0.0116** | -0.0118 | -0.2273 | +0.1079 | +0.37 | +0.3218 | -0.2112 |
| 15 | `ELK1` | P1P2 | 5843 | **+0.0115** | +0.0784 | -0.3195 | **+2.7142** | +0.52 | +0.2409 | +0.3805 |
| 16 | `BATF2` | P1P2 | 5844 | **-0.0064** | -0.0128 | -0.6139 | +0.2037 | -0.26 | +0.1562 | +0.0576 |
| 17 | `SP100` | P1P2 | 5843 | **-0.0086** | +0.0576 | -0.2379 | **+1.7837** | -0.14 | +0.0708 | +0.2494 |
| 18 | `TRAFD1` | P1P2 | 5843 | **-0.0090** | +0.0135 | -0.0968 | +0.7173 | -0.18 | +0.0631 | -0.1742 |
| 19 | `ID2` | **P2** | 5844 | **-0.0102** | -0.0006 | -0.0983 | +0.4031 | -0.29 | +0.0237 | -0.1008 |
| 20 | `STAT3` | P1P2 | 5843 | **-0.0188** | +0.0556 | -0.1791 | **+1.4214** | -0.25 | +0.0423 | -0.3656 |
| 21 | `MAPK8` | P1P2 | 5843 | **-0.0311** | -0.0092 | -0.1133 | +0.0746 | -0.75 | +0.0766 | +0.3713 |
| 22 | `IFNAR2` | P1P2 | 5843 | **-0.0316** | -0.0200 | -0.6575 | +0.2994 | -1.48 | +0.4216 | -0.1969 |
| 23 | `RAP1GAP` | **P2** | 5844 | **-0.0468** | -0.0328 | -0.9588 | +0.1184 | -0.79 | +0.0264 | -0.0922 |
| 24 | `ETS1` | P1P2 | 5844 | **-0.0742** | +0.0239 | -0.1372 | +0.5067 | -1.90 | +0.0854 | +0.0487 |
| 25 | `STAT4` | P1P2 | 5844 | **-0.0782** | -0.0061 | -0.2645 | +0.1926 | -1.04 | +0.0023 | +0.0797 |
| 26 | `ID2` | **P1** | 5844 | **-0.1100** | -0.0147 | -0.3262 | +0.1266 | -2.11 | +0.0237 | -0.1008 |
| 27 | `IRF1` | P1P2 | 5844 | **-0.1343** | +0.0033 | -0.1122 | +0.1733 | -1.76 | +0.0668 | +0.5319 |
| 28 | `IRF7` | P1P2 | 5844 | **n/d** | -0.0300 | -0.8391 | +0.1285 | n/d | +0.0149 | n/d |
| 29 | `PARP12` | P1P2 | 5843 | **n/d** | -0.0118 | -0.2021 | +0.1433 | n/d | +0.0173 | n/d |

### 5.2 Las 29 filas de A549, ordenadas por `G_mismo`, sin gen blanco

| # | blanco | constr. | **`G_mismo`** | **`G_otro` media** | z vs nulo | `A_esp` |
|---|---|---|---|---|---|---|
| 1 | `ID2` | **P2** | **+0.3156** | -0.0212 | +1.21 | -0.0067 |
| 2 | `CEBPG` | P1P2 | **+0.1525** | +0.0097 | +2.53 | +0.0515 |
| 3 | `SP100` | P1P2 | **+0.1077** | +0.0049 | +1.94 | +0.0399 |
| 4 | `IRF1` | P1P2 | **+0.0945** | -0.0060 | +2.15 | +0.1091 |
| 5 | `ETS2` | P1P2 | **+0.0914** | +0.0084 | +1.62 | +0.1015 |
| 6 | `BRD9` | P1P2 | **+0.0902** | +0.0212 | +1.12 | +0.0322 |
| 7 | `RAP1GAP` | **P1** | **+0.0600** | -0.0057 | +1.31 | +0.0290 |
| 8 | `MAPK8` | P1P2 | **+0.0595** | -0.0242 | +1.02 | +0.0380 |
| 9 | `AKT1` | P1P2 | **+0.0582** | +0.0103 | +0.96 | +0.0351 |
| 10 | `RNF114` | P1P2 | **+0.0574** | +0.0240 | +0.76 | +0.0252 |
| 11 | `STAT1` | P1P2 | **+0.0565** | +0.0115 | +1.61 | +0.0912 |
| 12 | `RAP1GAP` | **P2** | **+0.0412** | -0.0110 | +0.76 | +0.0290 |
| 13 | `JUN` | P1P2 | **+0.0388** | +0.0069 | +0.96 | +0.1476 |
| 14 | `ELK1` | P1P2 | **+0.0311** | -0.0121 | +0.67 | +0.0978 |
| 15 | `JAK1` | P1P2 | **+0.0237** | -0.0107 | +0.67 | +0.1832 |
| 16 | `STAT3` | P1P2 | **+0.0222** | +0.0116 | +0.81 | +0.2238 |
| 17 | `STAT4` | P1P2 | **+0.0173** | -0.0250 | +0.15 | +0.0218 |
| 18 | `IRF7` | P1P2 | **+0.0168** | -0.0041 | +0.35 | +0.0683 |
| 19 | `IFNAR1` | P1P2 | **+0.0129** | -0.0040 | +0.73 | +0.3308 |
| 20 | `SOCS1` | P1P2 | **+0.0067** | -0.0082 | +0.05 | +0.0346 |
| 21 | `TRAFD1` | P1P2 | **+0.0038** | -0.0103 | +0.17 | +0.1061 |
| 22 | `PARP12` | P1P2 | **-0.0295** | -0.0170 | -0.40 | +0.0013 |
| 23 | `TYK2` | P1P2 | **-0.0400** | +0.0012 | -2.03 | +0.4975 |
| 24 | `IFNAR2` | P1P2 | **-0.0495** | -0.0028 | -2.41 | +0.4998 |
| 25 | `STAT5A` | P1P2 | **-0.0511** | +0.0128 | -0.86 | +0.0368 |
| 26 | `BATF2` | P1P2 | **-0.0567** | -0.0140 | -1.72 | +0.1011 |
| 27 | `ID2` | **P1** | **-0.0965** | -0.0097 | -0.31 | -0.0067 |
| 28 | `ETS1` | P1P2 | **-0.1869** | +0.0267 | -2.37 | +0.0087 |
| 29 | `VAV1` | P1P2 | **n/d** | -0.0217 | n/d | -0.0493 |

---

## 6. Defectos de la construccion, reportados como datos sobre la medicion

### 6.1 `G_otro` sale fuera de [-1, 1] en 4 filas de 165

Valores observados: `ELK1` en K562 da `G_otro` max **+2.7142**; `SP100`
**+1.7837**; `STAT3` **+1.4214**; `BRD9` min **-1.0959**. Las cuatro en K562.

Es la misma propiedad del estimador que `control_o_no.md` §8.1 punto 1: `G`
divide por `sqrt(<a,b>)`, que es una **estimacion** de `||senal||`, y cuando
`<a,b>` es chico el cociente no esta acotado por 1. **No es una Spearman
acotada.** Por eso §3.1 da tambien la mediana y no solo la media.

El maximo \|`G_otro`\| por linea: A549 0.6240, BXPC3 0.1677, HAP1 0.4713,
HT29 0.4876, **K562 2.7142**, MCF7 0.2910.

### 6.2 `n/d`: 9 filas de 174, por `<a,b> <= 0` en espacio de rangos

| linea | filas `n/d` |
|---|---|
| A549 | `VAV1` |
| BXPC3 | ninguna |
| HAP1 | `ID2` (P1), `ID2` (P2), `IRF7` |
| HT29 | `IRF7`, `SOCS1` |
| K562 | `IRF7`, `PARP12` |
| MCF7 | `ELK1` |

El producto cruzado de las dos mitades de Jiang sale <= 0 para ese blanco, asi
que el denominador no existe. `VAV1` en A549 tiene `cruz` en rangos
**-2.13e+08** y su `A_esp` es **-0.0493**, el unico negativo de A549 en
`control_o_no.md` §8.1.

### 6.3 Las dos cantidades `G_prom` tienen signos opuestos

`G_prom` **cubiertos** (promedio de los perfiles de los blancos cubiertos) es
**negativo en las 6 lineas**: -0.0069 a -0.0474.
`G_prom` **todas** (promedio de las 10744 perturbaciones dirigidas) es
**positivo en las 6 lineas**: +0.0044 a +0.0876.

Son dos promedios distintos del mismo objeto y dan signo contrario. Los dos estan
en §3.1 y en las tablas por fila; **no se elige uno**.

`G_prom` todas (+0.0655, +0.0734, +0.0876, +0.0134, +0.0246, +0.0044) es mayor
que `G_mismo` media (+0.0303, +0.0022, +0.0108, -0.0041, +0.0200, +0.0017) en las
6 lineas.

### 6.4 La fraccion de genes que se pierde al alinear

Se mide sobre **4825 a 5844 genes**, que es el **67.3 % al 76.7 %** del universo
de Jiang. El 23-33 % restante no esta en `var/gene_name` de gwps y **no entra en
ninguna cifra de este reporte**.

### 6.5 El numero de filas por linea no es el mismo

26 a 29 de las 29 posibles, por los `n/d` de §6.2. Las medias de §3.1 son sobre
distintos conjuntos de filas segun la linea, y **no estan pareadas entre lineas**.

### 6.6 Los controles `non-targeting` que no se usan

71 filas son `non-targeting` pero no `core_control`, y quedan fuera del
denominador del perfil. **No se midio** el perfil con las 585.

### 6.7 El perfil de gwps es de una sola replica por construccion

`X` tiene una fila por perturbacion; no hay replicas dentro de una construccion,
asi que **no se puede aplicar al lado de gwps el mismo producto cruzado que se
aplica al de Jiang**. El ruido de gwps entra sin descontar en todas las cifras
de `G_mismo`, `G_otro` y `G_prom`. La unica excepcion son `ID2` y `RAP1GAP`,
donde P1 y P2 son dos mediciones del mismo gen (§4), y **no se promediaron**.

---

## 7. Lo que no se hizo, y las entradas

| prohibicion | estado |
|---|---|
| construir un candidato, tocar arneses o bundles | **no**. El script lee dos archivos y escribe un JSON |
| usar cifras de `~/vcc/benchmark/` o `/datos/benchmark_VEDADO_techo_metricas/` | **no leidos**. 🔴 Los productos de gwps del proyecto viven en `~/vcc/benchmark/out/` y **no se abrio ninguno**: el perfil se construyo desde el h5ad crudo en `/datos/ref/replogle/`. El dato de `gwps_fiabilidad.md` sobre P1 (67.6 % de 546 pares) se usa **solo** como la instruccion de no promediar, tal como venia en el encargo; no entra en ningun calculo |
| escribir conclusiones | **no**. No hay seccion de veredicto ni de interpretacion |
| envios, `vcc sample`, commits | **ninguno** |

| | |
|---|---|
| entradas | `/datos/ref/replogle/k562_gwps_raw_bulk.h5ad`; `modelo/out/control_o_no/rasgos_{A549,BXPC3,HAP1,HT29,K562,MCF7}.npz`; `modelo/out/suelo_panel.json` |
| codigo | `modelo/src/gwps_jiang_datos.py` |
| salidas | `modelo/out/gwps_jiang_datos.json`, `modelo/out/control_o_no/gwps_jiang_run.log` |
| semillas | `20270000 + 100*fold`, solo para el nulo por permutacion (200 por fila) |

**No se midio:**

1. El perfil de gwps con las 585 filas `non-targeting` en vez de las 514
   `core_control`.
2. Ninguna version que promedie P1 y P2.
3. Los 6 blancos del panel que no estan en gwps (§2.2).
4. Los 23-33 % de genes de Jiang que no estan en el universo de gwps.
5. `k562_essential` ni `rpe1`, que estan en el mismo directorio.
6. Cualquier cantidad sobre el lado perturbado de Jiang distinta de `E1`/`E2`
   ya calculadas.
