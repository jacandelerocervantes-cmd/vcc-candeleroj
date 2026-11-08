# El Venn de los conjuntos ED entre contextos (Jiang, 6 lineas, 33 perturbaciones)

## 🔴 0. Encabezado obligatorio: esto es una COTA DE JIANG

Todo lo que sigue se mide sobre las **seis lineas celulares de Jiang**, condicion IFNB, con el panel comun de 33 perturbaciones. **No se extrapola a D/E/F en ningun sentido, ni como cota superior ni como cota inferior.**

La razon esta MEDIDA, no supuesta. La heterogeneidad que produce esta cota -- que las seis lineas tengan conjuntos ED de tamanos y contenidos muy distintos -- esta **medida como AUSENTE en el reto** (`ntc_nulo.md` §3):

| #  | cantidad | reto 2026: max/min entre contextos | Jiang: max/min entre lineas |
|---|---|---|---|
| 1 | celulas de control por unidad | **x1.000** | **x2.341** |
| 2 | profundidad mediana (UMIs/celula) | **x1.028** | **x1.869** |
| 3 | `nnz` mediana por celula | **x1.092** | **x1.602** |
| 4 | panel probable (el `m` de BH) | **x1.062** | **x1.077** |

Y el nulo de `ntc_id` en contexto **da cero** conjunto ED falso positivo (`ntc_nulo.md` §2.1). Las dos cosas juntas dicen que el mecanismo que genera esta cota en Jiang **no esta medido como presente en el reto**. Lo que esta cota ES: una medicion de Jiang. Lo que NO es: una prediccion sobre D/E/F.

## 1. La construccion, declarada antes de medir

| #  | cantidad | valor |
|---|---|---|
| 1 | conjunto ED real | `cell-eval2 run --write-degenes` sobre la poblacion real + `prepare_de` + `p_adj < 0.05` + `_exclude_own_gene` -- la maquinaria del propio marcador, igual que `pds_objetivo.md` §3.2 |
| 2 | umbral | `p_adj < 0.05` (el del preset `vcc2026`) |
| 3 | poblacion real | `seed_real = 20260001 + fold`, la MISMA de siempre, verificada por sha256 |
| 4 | sondas | 6, una por linea (el lado real depende solo de `seed_real`) |
| 5 | genes del espacio | 34025 |
| 6 | panel | 33 perturbaciones, el comun de las 6 lineas |
| 7 | semilla de desempate del recorte | `20280101` (familia nueva 20280000+) |
| 8 | `jac` | interseccion / union, crudo. **Dos convenciones declaradas** para la union vacia: la del marcador (`= 1.0`) y excluyendo esas celdas |

### 1.1 Verificacion de la poblacion real y recomprobacion del cache

| #  | linea | huella real (sha256, 16) | coincide con la de disco | # json que la comparten | perts identicas al cache | perts distintas | genes solo ahora | genes solo en cache | duracion |
|---|---|---|---|---|---|---|---|---|---|
| 1 | **A549** | 07194569d8167727... | si | 45 | 33 | 0 | 0 | 0 | 71.7 s |
| 2 | **BXPC3** | 41383fe61b7f5bf9... | si | 45 | 33 | 0 | 0 | 0 | 1.8 min (106 s) |
| 3 | **HAP1** | 3a28f744c9fc9c83... | si | 45 | 33 | 0 | 0 | 0 | 1.8 min (108 s) |
| 4 | **HT29** | f65c2ef7f2d9c97a... | si | 45 | 33 | 0 | 0 | 0 | 1.5 min (91 s) |
| 5 | **K562** | c7d34fa6b06f7f0d... | si | 45 | 33 | 0 | 0 | 0 | 73.8 s |
| 6 | **MCF7** | d1c8674360476df0... | si | 45 | 33 | 0 | 0 | 0 | 1.7 min (105 s) |

El cache `/datos/modelo/out/pds_obj/ed_real.json` se recomprobo **gen a gen**, no por conteo. Perturbaciones identicas: **198 de 198**.

## 2. Tamano del conjunto ED real, por linea y por perturbacion

🔴 **Rango, nunca media.** Con 15 de 33 perturbaciones vacias en A549, el techo ahi es **0 por construccion en casi medio panel**, y una media lo esconderia.

| #  | linea | **# perts con ED vacio** | # con ED no vacio | rango [min, max] de |ED| | mediana |ED| | suma |ED| |
|---|---|---|---|---|---|---|
| 1 | **A549** | 15 | 18 | [0, 1021] | 1 | 2537 |
| 2 | **BXPC3** | 0 | 33 | [1, 1164] | 32 | 4642 |
| 3 | **HAP1** | 14 | 19 | [0, 207] | 2 | 817 |
| 4 | **HT29** | 6 | 27 | [0, 615] | 4 | 2505 |
| 5 | **K562** | 8 | 25 | [0, 281] | 2 | 1042 |
| 6 | **MCF7** | 0 | 33 | [11, 1438] | 197 | 12867 |

La columna de vacios reproduce `control_permutado.md` §7 (15, 0, 14, 6, 8, 0): **coincide en las 6**.

### 🔴 TABLA 2.1 -- los intermedios crudos: |ED| por perturbacion y linea, las 33 filas

| #  | pert | **A549** | **BXPC3** | **HAP1** | **HT29** | **K562** | **MCF7** | # lineas vacias | rango |ED| | genes en 6 | genes en 1 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `AKT1` | 0 | 3 | 0 | 3 | 1 | 197 | 2 | [0, 197] | 0 | 204 |
| 2 | `BATF2` | 0 | 74 | 7 | 2 | 1 | 131 | 1 | [0, 131] | 0 | 211 |
| 3 | `BRD9` | 1 | 107 | 0 | 0 | 32 | 56 | 2 | [0, 107] | 0 | 194 |
| 4 | `CEBPG` | 4 | 33 | 18 | 10 | 2 | 294 | 0 | [2, 294] | 0 | 353 |
| 5 | `DTX3L` | 3 | 9 | 3 | 14 | 1 | 102 | 0 | [1, 102] | 1 | 126 |
| 6 | `ELK1` | 0 | 17 | 15 | 1 | 35 | 164 | 1 | [0, 164] | 0 | 228 |
| 7 | `ETS1` | 0 | 385 | 2 | 5 | 8 | 236 | 1 | [0, 385] | 0 | 622 |
| 8 | `ETS2` | 19 | 9 | 0 | 192 | 31 | 72 | 1 | [0, 192] | 0 | 313 |
| 9 | `H1F0` | 0 | 10 | 0 | 22 | 2 | 529 | 2 | [0, 529] | 0 | 555 |
| 10 | `HERC5` | 0 | 8 | 4 | 0 | 2 | 458 | 2 | [0, 458] | 0 | 472 |
| 11 | `HERC6` | 0 | 31 | 4 | 1 | 1 | 67 | 1 | [0, 67] | 0 | 99 |
| 12 | `ID2` | 0 | 80 | 0 | 15 | 1 | 1438 | 2 | [0, 1438] | 0 | 1518 |
| 13 | `IFNAR1` | 494 | 531 | 184 | 275 | 189 | 842 | 0 | [184, 842] | 36 | 1435 |
| 14 | `IFNAR2` | 751 | 1164 | 177 | 615 | 281 | 1368 | 0 | [177, 1368] | 44 | 2540 |
| 15 | `IRF1` | 10 | 72 | 0 | 2 | 0 | 163 | 2 | [0, 163] | 0 | 231 |
| 16 | `IRF7` | 0 | 5 | 10 | 0 | 4 | 474 | 2 | [0, 474] | 0 | 488 |
| 17 | `JAK1` | 94 | 363 | 60 | 423 | 139 | 115 | 0 | [60, 423] | 0 | 810 |
| 18 | `JUN` | 6 | 39 | 0 | 42 | 5 | 627 | 1 | [0, 627] | 0 | 705 |
| 19 | `MAP3K14` | 40 | 50 | 0 | 2 | 1 | 175 | 1 | [0, 175] | 0 | 266 |
| 20 | `MAPK8` | 1 | 72 | 0 | 3 | 0 | 1053 | 2 | [0, 1053] | 0 | 1121 |
| 21 | `PARP12` | 0 | 4 | 2 | 51 | 0 | 11 | 2 | [0, 51] | 0 | 66 |
| 22 | `RAP1GAP` | 2 | 8 | 0 | 0 | 0 | 31 | 3 | [0, 31] | 0 | 41 |
| 23 | `RNF114` | 0 | 38 | 4 | 49 | 4 | 233 | 1 | [0, 233] | 0 | 322 |
| 24 | `SOCS1` | 0 | 32 | 7 | 3 | 0 | 460 | 2 | [0, 460] | 0 | 500 |
| 25 | `SP100` | 1 | 65 | 0 | 13 | 1 | 93 | 1 | [0, 93] | 0 | 165 |
| 26 | `STAT1` | 1 | 264 | 207 | 465 | 64 | 67 | 0 | [1, 465] | 0 | 801 |
| 27 | `STAT3` | 84 | 15 | 2 | 0 | 2 | 44 | 1 | [0, 84] | 0 | 147 |
| 28 | `STAT4` | 0 | 13 | 0 | 0 | 0 | 20 | 4 | [0, 20] | 0 | 33 |
| 29 | `STAT5A` | 0 | 1 | 2 | 3 | 61 | 187 | 1 | [0, 187] | 0 | 252 |
| 30 | `TRAFD1` | 3 | 9 | 0 | 4 | 0 | 1140 | 2 | [0, 1140] | 0 | 1150 |
| 31 | `TYK2` | 1021 | 1108 | 107 | 278 | 172 | 1027 | 0 | [107, 1108] | 45 | 2068 |
| 32 | `VAV1` | 0 | 6 | 0 | 10 | 2 | 415 | 2 | [0, 415] | 0 | 431 |
| 33 | `ZBP1` | 2 | 17 | 2 | 2 | 0 | 578 | 1 | [0, 578] | 0 | 601 |

## 3. Cuantos genes se llaman en 6, 5, 4, 3, 2 y 1 sola linea

🔴 **Dos definiciones distintas, las dos reportadas, porque no son la misma cantidad.**

* **A -- pares (perturbacion, gen)**: para cada perturbacion se cuenta en cuantas lineas se llama cada gen, y se SUMAN las 33 distribuciones. Es la definicion que cita el encargo.
* **B -- genes distintos**: un gen cuenta las lineas en que aparece en **alguna** perturbacion. Colapsa el panel, asi que da otro numero.

**(a) CRUDA**

| # lineas en que se llama | **A: pares (pert, gen)** | frac A | **B: genes distintos** | frac B |
|---|---|---|---|---|
| 6 | 126 | 0.60 % | 72 | 0.86 % |
| 5 | 115 | 0.55 % | 71 | 0.85 % |
| 4 | 133 | 0.63 % | 172 | 2.05 % |
| 3 | 355 | 1.69 % | 598 | 7.14 % |
| 2 | 1207 | 5.75 % | 1871 | 22.33 % |
| 1 | 19068 | 90.78 % | 5593 | 66.77 % |

Total A = **21004** pares, total B = **8377** genes distintos.

**(b) IGUALADA EN POTENCIA**

| # lineas en que se llama | **A: pares (pert, gen)** | frac A | **B: genes distintos** | frac B |
|---|---|---|---|---|
| 6 | 73 | 3.39 % | 34 | 2.61 % |
| 5 | 83 | 3.85 % | 35 | 2.68 % |
| 4 | 70 | 3.25 % | 30 | 2.30 % |
| 3 | 108 | 5.01 % | 42 | 3.22 % |
| 2 | 192 | 8.91 % | 126 | 9.66 % |
| 1 | 1629 | 75.59 % | 1038 | 79.54 % |

Total A = **2155** pares, total B = **1305** genes distintos.

🔴 **La recomprobacion del encargo.** El cache citado decia 19,068 de 21,004 en una sola linea y 126 en las seis. Medido ahora, definicion A: **19068 de 21004** en una sola y **126** en las seis. Coincide.

## 4. 🔴 El techo de `jac` para una emision IDENTICA en las seis lineas

**La interseccion NO es el optimo.** Con `R_1={a,b}` y `R_2={a,c}`, `S={a}` da media `0.5000` y `S={a,b}` da `0.6667`. Por eso no se toma la interseccion: se **barre** el corte de frecuencia. `S_k` = los genes llamados en **k o mas** lineas, con `k = 6, 5, 4, 3, 2, 1`.

🔴 El maximo de la curva es una **COTA INFERIOR del techo, restringida a esta familia de `S`** (los seis cortes de frecuencia). **No es el techo exacto**: el optimo sobre TODOS los subconjuntos posibles no se calcula aqui y puede ser mayor.

### 4.1 (a) CRUDA -- la curva entera, agregada sobre las 33 perturbaciones

| corte | rango |S| sobre las 33 perts | mediana |S| | rango de la media de `jac` (conv. marcador) | media de medias | rango (sin union vacia) | media de medias | # celdas union vacia |
|---|---|---|---|---|---|---|---|
| `>= 6` | [0, 45] | 0 | [+0.0000, +0.6667] | 0.2381 | [+0.0000, +0.3098] | 0.0209 | 43 |
| `>= 5` | [0, 87] | 0 | [+0.0000, +0.6667] | 0.2461 | [+0.0000, +0.3098] | 0.0289 | 43 |
| `>= 4` | [0, 117] | 0 | [+0.0000, +0.6667] | 0.2525 | [+0.0000, +0.3098] | 0.0454 | 41 |
| `>= 3` | [0, 244] | 0 | [+0.0000, +0.6667] | 0.2482 | [+0.0000, +0.3098] | 0.0562 | 38 |
| `>= 2` | [0, 664] | 2 | [+0.0036, +0.6667] | 0.1431 | [+0.0000, +0.3098] | 0.0775 | 13 |
| `>= 1` | [33, 3204] | 357 | [+0.1667, +0.2331] | 0.1760 | [+0.1667, +0.2331] | 0.1760 | 0 |

**El maximo de la curva, como COTA INFERIOR del techo (familia de cortes de frecuencia):** con la convencion del marcador, el rango sobre las 33 perturbaciones es [+0.1667, +0.6667] y el maximo lo da `STAT4` con **0.6667**; sin las celdas de union vacia, rango [+0.1667, +0.3098], maximo `DTX3L` con **0.3098**.

#### 🔴 TABLA 4.1b -- los intermedios crudos: la curva por perturbacion, las 33 filas

| #  | pert | genes llamados en alguna | \|S\| >=6 | \|S\| >=5 | \|S\| >=4 | \|S\| >=3 | \|S\| >=2 | \|S\| >=1 | jac >=6 | jac >=5 | jac >=4 | jac >=3 | jac >=2 | jac >=1 | argmax | **max (cota inf.)** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `AKT1` | 204 | 0 | 0 | 0 | 0 | 0 | 204 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 2 | `BATF2` | 213 | 0 | 0 | 0 | 0 | 2 | 213 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.0266 | 0.1682 | `>= 1` | 0.1682 |
| 3 | `BRD9` | 195 | 0 | 0 | 0 | 0 | 1 | 195 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.0045 | 0.1675 | `>= 6` | 0.3333 |
| 4 | `CEBPG` | 357 | 0 | 0 | 0 | 0 | 4 | 357 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0858 | 0.1685 | `>= 1` | 0.1685 |
| 5 | `DTX3L` | 127 | 1 | 1 | 1 | 1 | 1 | 127 | 0.3098 | 0.3098 | 0.3098 | 0.3098 | 0.3098 | 0.1732 | `>= 6` | 0.3098 |
| 6 | `ELK1` | 230 | 0 | 0 | 0 | 0 | 2 | 230 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.0252 | 0.1681 | `>= 1` | 0.1681 |
| 7 | `ETS1` | 629 | 0 | 0 | 0 | 0 | 7 | 629 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.0395 | 0.1685 | `>= 1` | 0.1685 |
| 8 | `ETS2` | 318 | 0 | 0 | 0 | 0 | 5 | 318 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.0512 | 0.1693 | `>= 1` | 0.1693 |
| 9 | `H1F0` | 559 | 0 | 0 | 0 | 0 | 4 | 559 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.0429 | 0.1679 | `>= 6` | 0.3333 |
| 10 | `HERC5` | 472 | 0 | 0 | 0 | 0 | 0 | 472 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 11 | `HERC6` | 101 | 0 | 0 | 0 | 1 | 2 | 101 | 0.1667 | 0.1667 | 0.1667 | 0.0495 | 0.0491 | 0.1716 | `>= 1` | 0.1716 |
| 12 | `ID2` | 1526 | 0 | 0 | 0 | 0 | 8 | 1526 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.0290 | 0.1675 | `>= 6` | 0.3333 |
| 13 | `IFNAR1` | 1798 | 36 | 68 | 97 | 153 | 363 | 1798 | 0.1167 | 0.1977 | 0.2443 | 0.2849 | 0.3064 | 0.2331 | `>= 2` | 0.3064 |
| 14 | `IFNAR2` | 3204 | 44 | 87 | 113 | 244 | 664 | 3204 | 0.1009 | 0.1699 | 0.1972 | 0.2460 | 0.2786 | 0.2266 | `>= 2` | 0.2786 |
| 15 | `IRF1` | 239 | 0 | 0 | 0 | 0 | 8 | 239 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.0529 | 0.1722 | `>= 6` | 0.3333 |
| 16 | `IRF7` | 490 | 0 | 0 | 0 | 1 | 2 | 490 | 0.3333 | 0.3333 | 0.3333 | 0.0504 | 0.0618 | 0.1677 | `>= 6` | 0.3333 |
| 17 | `JAK1` | 953 | 0 | 7 | 30 | 61 | 143 | 953 | 0.0000 | 0.0462 | 0.1560 | 0.2160 | 0.2371 | 0.2088 | `>= 2` | 0.2371 |
| 18 | `JUN` | 712 | 0 | 0 | 0 | 0 | 7 | 712 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.0365 | 0.1683 | `>= 1` | 0.1683 |
| 19 | `MAP3K14` | 267 | 0 | 0 | 0 | 0 | 1 | 267 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.0075 | 0.1673 | `>= 1` | 0.1673 |
| 20 | `MAPK8` | 1125 | 0 | 0 | 0 | 0 | 4 | 1125 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.0099 | 0.1673 | `>= 6` | 0.3333 |
| 21 | `PARP12` | 67 | 0 | 0 | 0 | 0 | 1 | 67 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.0985 | 0.1692 | `>= 6` | 0.3333 |
| 22 | `RAP1GAP` | 41 | 0 | 0 | 0 | 0 | 0 | 41 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.1667 | `>= 6` | 0.5000 |
| 23 | `RNF114` | 324 | 0 | 0 | 1 | 1 | 2 | 324 | 0.1667 | 0.1667 | 0.0911 | 0.0911 | 0.0822 | 0.1687 | `>= 1` | 0.1687 |
| 24 | `SOCS1` | 501 | 0 | 0 | 0 | 0 | 1 | 501 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.0056 | 0.1670 | `>= 6` | 0.3333 |
| 25 | `SP100` | 168 | 0 | 0 | 1 | 1 | 3 | 168 | 0.1667 | 0.1667 | 0.1838 | 0.1838 | 0.0898 | 0.1716 | `>= 4` | 0.1838 |
| 26 | `STAT1` | 906 | 0 | 1 | 14 | 42 | 105 | 906 | 0.0000 | 0.0069 | 0.0637 | 0.1541 | 0.1828 | 0.1965 | `>= 1` | 0.1965 |
| 27 | `STAT3` | 147 | 0 | 0 | 0 | 0 | 0 | 147 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | `>= 6` | 0.1667 |
| 28 | `STAT4` | 33 | 0 | 0 | 0 | 0 | 0 | 33 | 0.6667 | 0.6667 | 0.6667 | 0.6667 | 0.6667 | 0.1667 | `>= 6` | 0.6667 |
| 29 | `STAT5A` | 253 | 0 | 0 | 0 | 0 | 1 | 253 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.0036 | 0.1673 | `>= 1` | 0.1673 |
| 30 | `TRAFD1` | 1153 | 0 | 0 | 0 | 0 | 3 | 1153 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.0823 | 0.1671 | `>= 6` | 0.3333 |
| 31 | `TYK2` | 2659 | 45 | 77 | 117 | 224 | 591 | 2659 | 0.1621 | 0.2224 | 0.2525 | 0.2701 | 0.2732 | 0.2327 | `>= 2` | 0.2732 |
| 32 | `VAV1` | 432 | 0 | 0 | 0 | 0 | 1 | 432 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.0837 | 0.1671 | `>= 6` | 0.3333 |
| 33 | `ZBP1` | 601 | 0 | 0 | 0 | 0 | 0 | 601 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | `>= 6` | 0.1667 |

### 4.2 (b) IGUALADA EN POTENCIA -- la curva entera, agregada sobre las 33 perturbaciones

| corte | rango |S| sobre las 33 perts | mediana |S| | rango de la media de `jac` (conv. marcador) | media de medias | rango (sin union vacia) | media de medias | # celdas union vacia |
|---|---|---|---|---|---|---|---|
| `>= 6` | [0, 26] | 0 | [+0.0000, +1.0000] | 0.2620 | [+0.0000, +1.0000] | 0.0448 | 43 |
| `>= 5` | [0, 54] | 0 | [+0.0000, +1.0000] | 0.2765 | [+0.0000, +1.0000] | 0.0593 | 43 |
| `>= 4` | [0, 76] | 0 | [+0.0000, +1.0000] | 0.2860 | [+0.0000, +1.0000] | 0.0739 | 42 |
| `>= 3` | [0, 105] | 0 | [+0.0000, +1.0000] | 0.2969 | [+0.0000, +1.0000] | 0.0848 | 42 |
| `>= 2` | [0, 173] | 0 | [+0.0000, +1.0000] | 0.2945 | [+0.0000, +1.0000] | 0.1228 | 34 |
| `>= 1` | [1, 670] | 8 | [+0.1667, +1.0000] | 0.2115 | [+0.1667, +1.0000] | 0.2115 | 0 |

**El maximo de la curva, como COTA INFERIOR del techo (familia de cortes de frecuencia):** con la convencion del marcador, el rango sobre las 33 perturbaciones es [+0.1667, +1.0000] y el maximo lo da `DTX3L` con **1.0000**; sin las celdas de union vacia, rango [+0.1667, +1.0000], maximo `DTX3L` con **1.0000**.

#### 🔴 TABLA 4.2b -- los intermedios crudos: la curva por perturbacion, las 33 filas

| #  | pert | genes llamados en alguna | \|S\| >=6 | \|S\| >=5 | \|S\| >=4 | \|S\| >=3 | \|S\| >=2 | \|S\| >=1 | jac >=6 | jac >=5 | jac >=4 | jac >=3 | jac >=2 | jac >=1 | argmax | **max (cota inf.)** |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `AKT1` | 4 | 0 | 0 | 0 | 0 | 0 | 4 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 2 | `BATF2` | 5 | 0 | 0 | 0 | 0 | 0 | 5 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | `>= 6` | 0.1667 |
| 3 | `BRD9` | 4 | 0 | 0 | 0 | 0 | 0 | 4 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 4 | `CEBPG` | 11 | 0 | 0 | 0 | 0 | 1 | 11 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1667 | 0.1818 | `>= 1` | 0.1818 |
| 5 | `DTX3L` | 1 | 1 | 1 | 1 | 1 | 1 | 1 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | `>= 6` | 1.0000 |
| 6 | `ELK1` | 4 | 0 | 0 | 0 | 0 | 1 | 4 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.3333 | 0.2083 | `>= 2` | 0.3333 |
| 7 | `ETS1` | 10 | 0 | 0 | 0 | 0 | 0 | 10 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | `>= 6` | 0.1667 |
| 8 | `ETS2` | 44 | 0 | 0 | 0 | 0 | 1 | 44 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.0370 | 0.1705 | `>= 1` | 0.1705 |
| 9 | `H1F0` | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 10 | `HERC5` | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 11 | `HERC6` | 5 | 0 | 0 | 0 | 0 | 0 | 5 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | `>= 6` | 0.1667 |
| 12 | `ID2` | 4 | 0 | 0 | 0 | 0 | 0 | 4 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 13 | `IFNAR1` | 670 | 26 | 54 | 76 | 105 | 173 | 670 | 0.1413 | 0.2626 | 0.3303 | 0.3771 | 0.3993 | 0.2746 | `>= 2` | 0.3993 |
| 14 | `IFNAR2` | 641 | 25 | 50 | 74 | 103 | 169 | 641 | 0.1412 | 0.2534 | 0.3292 | 0.3812 | 0.4053 | 0.2761 | `>= 2` | 0.4053 |
| 15 | `IRF1` | 7 | 0 | 0 | 0 | 0 | 1 | 7 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | 0.1905 | `>= 6` | 0.3333 |
| 16 | `IRF7` | 15 | 0 | 0 | 0 | 0 | 1 | 15 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.0833 | 0.1778 | `>= 6` | 0.3333 |
| 17 | `JAK1` | 238 | 0 | 4 | 14 | 45 | 59 | 238 | 0.0000 | 0.0556 | 0.1636 | 0.3785 | 0.3944 | 0.2521 | `>= 2` | 0.3944 |
| 18 | `JUN` | 24 | 0 | 0 | 0 | 0 | 1 | 24 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.0667 | 0.1736 | `>= 1` | 0.1736 |
| 19 | `MAP3K14` | 5 | 0 | 0 | 0 | 0 | 0 | 5 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | `>= 6` | 0.1667 |
| 20 | `MAPK8` | 4 | 0 | 0 | 0 | 0 | 0 | 4 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 21 | `PARP12` | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 22 | `RAP1GAP` | 6 | 0 | 0 | 0 | 0 | 0 | 6 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.5000 | 0.1667 | `>= 6` | 0.5000 |
| 23 | `RNF114` | 17 | 0 | 0 | 1 | 1 | 1 | 17 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1961 | `>= 1` | 0.1961 |
| 24 | `SOCS1` | 12 | 0 | 0 | 0 | 0 | 0 | 12 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 25 | `SP100` | 4 | 0 | 0 | 0 | 0 | 1 | 4 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.3333 | 0.2083 | `>= 2` | 0.3333 |
| 26 | `STAT1` | 6 | 0 | 0 | 0 | 0 | 0 | 6 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.0000 | 0.1667 | `>= 1` | 0.1667 |
| 27 | `STAT3` | 10 | 0 | 0 | 0 | 0 | 0 | 10 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | `>= 6` | 0.1667 |
| 28 | `STAT4` | 26 | 0 | 0 | 0 | 0 | 0 | 26 | 0.6667 | 0.6667 | 0.6667 | 0.6667 | 0.6667 | 0.1667 | `>= 6` | 0.6667 |
| 29 | `STAT5A` | 5 | 0 | 0 | 0 | 0 | 0 | 5 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | `>= 6` | 0.1667 |
| 30 | `TRAFD1` | 12 | 0 | 0 | 0 | 0 | 0 | 12 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 31 | `TYK2` | 319 | 21 | 47 | 60 | 79 | 116 | 319 | 0.1963 | 0.3851 | 0.4480 | 0.4946 | 0.4988 | 0.3354 | `>= 2` | 0.4988 |
| 32 | `VAV1` | 8 | 0 | 0 | 0 | 0 | 0 | 8 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.3333 | 0.1667 | `>= 6` | 0.3333 |
| 33 | `ZBP1` | 10 | 0 | 0 | 0 | 0 | 0 | 10 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | 0.1667 | `>= 6` | 0.1667 |

### 4.3 Rango POR LINEA en el corte de cada maximo, nunca media

| #  | linea | # perts con ED real vacio | # celdas con jac definido | rango de `jac` en el corte del maximo | mediana |
|---|---|---|---|---|---|
| 1 | **A549** | 15 | 24 | [+0.0000, +0.3580] | 0.0000 |
| 2 | **BXPC3** | 0 | 33 | [+0.0000, +0.6121] | 0.0040 |
| 3 | **HAP1** | 14 | 23 | [+0.0000, +0.3624] | 0.0079 |
| 4 | **HT29** | 6 | 27 | [+0.0000, +0.6038] | 0.0094 |
| 5 | **K562** | 8 | 25 | [+0.0000, +1.0000] | 0.0099 |
| 6 | **MCF7** | 0 | 33 | [+0.0000, +0.8806] | 0.0098 |

🔴 En A549 y HAP1 el techo es **0 por construccion** en 15 y 14 de las 33 perturbaciones: no hay conjunto real con el que solapar. Esa es la razon de no dar media de panel.

## 5. 🔴 LA CIFRA QUE HABLA DE CONTEXTO: la diferencia entre (a) y (b)

Las dos versiones son **construcciones distintas**, declaradas como distintas. (a) usa los conjuntos tal como salen. (b) recorta las seis al tamano del conjunto ED mas pequeno **no vacio** de esa perturbacion, tomando los de `p_adj` mas baja. Si (b) sube mucho respecto de (a), lo que limitaba el solape era el **tamano** (la potencia); si no sube, lo que limitaba era **que genes** (el contexto).

| #  | pert | m (tamano del mas chico no vacio) | # lineas no vacias | max jac (a) cruda | max jac (b) igualada | **(b) - (a)** | genes llamados (a) | genes llamados (b) |
|---|---|---|---|---|---|---|---|---|
| 1 | `AKT1` | 1 | 4 | 0.3333 | 0.3333 | +0.0000 | 204 | 4 |
| 2 | `BATF2` | 1 | 5 | 0.1682 | 0.1667 | -0.0016 | 213 | 5 |
| 3 | `BRD9` | 1 | 4 | 0.3333 | 0.3333 | +0.0000 | 195 | 4 |
| 4 | `CEBPG` | 2 | 6 | 0.1685 | 0.1818 | +0.0133 | 357 | 11 |
| 5 | `DTX3L` | 1 | 6 | 0.3098 | 1.0000 | +0.6902 | 127 | 1 |
| 6 | `ELK1` | 1 | 5 | 0.1681 | 0.3333 | +0.1652 | 230 | 4 |
| 7 | `ETS1` | 2 | 5 | 0.1685 | 0.1667 | -0.0019 | 629 | 10 |
| 8 | `ETS2` | 9 | 5 | 0.1693 | 0.1705 | +0.0012 | 318 | 44 |
| 9 | `H1F0` | 2 | 4 | 0.3333 | 0.3333 | +0.0000 | 559 | 8 |
| 10 | `HERC5` | 2 | 4 | 0.3333 | 0.3333 | +0.0000 | 472 | 8 |
| 11 | `HERC6` | 1 | 5 | 0.1716 | 0.1667 | -0.0050 | 101 | 5 |
| 12 | `ID2` | 1 | 4 | 0.3333 | 0.3333 | +0.0000 | 1526 | 4 |
| 13 | `IFNAR1` | 184 | 6 | 0.3064 | 0.3993 | +0.0929 | 1798 | 670 |
| 14 | `IFNAR2` | 177 | 6 | 0.2786 | 0.4053 | +0.1266 | 3204 | 641 |
| 15 | `IRF1` | 2 | 4 | 0.3333 | 0.3333 | +0.0000 | 239 | 7 |
| 16 | `IRF7` | 4 | 4 | 0.3333 | 0.3333 | +0.0000 | 490 | 15 |
| 17 | `JAK1` | 60 | 6 | 0.2371 | 0.3944 | +0.1574 | 953 | 238 |
| 18 | `JUN` | 5 | 5 | 0.1683 | 0.1736 | +0.0053 | 712 | 24 |
| 19 | `MAP3K14` | 1 | 5 | 0.1673 | 0.1667 | -0.0006 | 267 | 5 |
| 20 | `MAPK8` | 1 | 4 | 0.3333 | 0.3333 | +0.0000 | 1125 | 4 |
| 21 | `PARP12` | 2 | 4 | 0.3333 | 0.3333 | +0.0000 | 67 | 8 |
| 22 | `RAP1GAP` | 2 | 3 | 0.5000 | 0.5000 | +0.0000 | 41 | 6 |
| 23 | `RNF114` | 4 | 5 | 0.1687 | 0.1961 | +0.0274 | 324 | 17 |
| 24 | `SOCS1` | 3 | 4 | 0.3333 | 0.3333 | +0.0000 | 501 | 12 |
| 25 | `SP100` | 1 | 5 | 0.1838 | 0.3333 | +0.1495 | 168 | 4 |
| 26 | `STAT1` | 1 | 6 | 0.1965 | 0.1667 | -0.0298 | 906 | 6 |
| 27 | `STAT3` | 2 | 5 | 0.1667 | 0.1667 | +0.0000 | 147 | 10 |
| 28 | `STAT4` | 13 | 2 | 0.6667 | 0.6667 | +0.0000 | 33 | 26 |
| 29 | `STAT5A` | 1 | 5 | 0.1673 | 0.1667 | -0.0007 | 253 | 5 |
| 30 | `TRAFD1` | 3 | 4 | 0.3333 | 0.3333 | +0.0000 | 1153 | 12 |
| 31 | `TYK2` | 107 | 6 | 0.2732 | 0.4988 | +0.2256 | 2659 | 319 |
| 32 | `VAV1` | 2 | 4 | 0.3333 | 0.3333 | +0.0000 | 432 | 8 |
| 33 | `ZBP1` | 2 | 5 | 0.1667 | 0.1667 | +0.0000 | 601 | 10 |

Rango de `(b) - (a)` sobre las 33 perturbaciones donde las dos estan definidas: **[-0.0298, +0.6902]**, mediana **+0.0000**. Se presenta la diferencia; no se concluye de ella.

### 5.1 Los intermedios del recorte: tamanos crudos, `m`, y empates en el corte

| #  | pert | \|ED\| A549 | \|ED\| BXPC3 | \|ED\| HAP1 | \|ED\| HT29 | \|ED\| K562 | \|ED\| MCF7 | **m** | empates en el corte, por linea |
|---|---|---|---|---|---|---|---|---|---|
| 1 | `AKT1` | 0 | 3 | 0 | 3 | 1 | 197 | 1 | BXPC3:1, HT29:1, MCF7:1 |
| 2 | `BATF2` | 0 | 74 | 7 | 2 | 1 | 131 | 1 | BXPC3:1, HAP1:1, HT29:2, MCF7:3 |
| 3 | `BRD9` | 1 | 107 | 0 | 0 | 32 | 56 | 1 | BXPC3:1, K562:1, MCF7:2 |
| 4 | `CEBPG` | 4 | 33 | 18 | 10 | 2 | 294 | 2 | A549:1, BXPC3:2, HAP1:1, HT29:1, MCF7:1 |
| 5 | `DTX3L` | 3 | 9 | 3 | 14 | 1 | 102 | 1 | A549:1, BXPC3:1, HAP1:1, HT29:1, MCF7:1 |
| 6 | `ELK1` | 0 | 17 | 15 | 1 | 35 | 164 | 1 | BXPC3:1, HAP1:1, K562:1, MCF7:1 |
| 7 | `ETS1` | 0 | 385 | 2 | 5 | 8 | 236 | 2 | BXPC3:1, HT29:5, K562:1, MCF7:5 |
| 8 | `ETS2` | 19 | 9 | 0 | 192 | 31 | 72 | 9 | A549:4, HT29:2, K562:1, MCF7:1 |
| 9 | `H1F0` | 0 | 10 | 0 | 22 | 2 | 529 | 2 | BXPC3:1, HT29:3, MCF7:6 |
| 10 | `HERC5` | 0 | 8 | 4 | 0 | 2 | 458 | 2 | BXPC3:1, HAP1:3, MCF7:1 |
| 11 | `HERC6` | 0 | 31 | 4 | 1 | 1 | 67 | 1 | BXPC3:1, HAP1:1, MCF7:1 |
| 12 | `ID2` | 0 | 80 | 0 | 15 | 1 | 1438 | 1 | BXPC3:1, HT29:1, MCF7:1 |
| 13 | `IFNAR1` | 494 | 531 | 184 | 275 | 189 | 842 | 184 | A549:1, BXPC3:1, HT29:1, K562:3, MCF7:1 |
| 14 | `IFNAR2` | 751 | 1164 | 177 | 615 | 281 | 1368 | 177 | A549:1, BXPC3:1, HT29:1, K562:1, MCF7:1 |
| 15 | `IRF1` | 10 | 72 | 0 | 2 | 0 | 163 | 2 | A549:1, BXPC3:1, MCF7:1 |
| 16 | `IRF7` | 0 | 5 | 10 | 0 | 4 | 474 | 4 | BXPC3:1, HAP1:5, MCF7:1 |
| 17 | `JAK1` | 94 | 363 | 60 | 423 | 139 | 115 | 60 | A549:2, BXPC3:1, HT29:1, K562:1, MCF7:3 |
| 18 | `JUN` | 6 | 39 | 0 | 42 | 5 | 627 | 5 | A549:2, BXPC3:2, HT29:3, MCF7:3 |
| 19 | `MAP3K14` | 40 | 50 | 0 | 2 | 1 | 175 | 1 | A549:1, BXPC3:1, HT29:1, MCF7:1 |
| 20 | `MAPK8` | 1 | 72 | 0 | 3 | 0 | 1053 | 1 | BXPC3:1, HT29:1, MCF7:1 |
| 21 | `PARP12` | 0 | 4 | 2 | 51 | 0 | 11 | 2 | BXPC3:1, HT29:5, MCF7:2 |
| 22 | `RAP1GAP` | 2 | 8 | 0 | 0 | 0 | 31 | 2 | BXPC3:1, MCF7:1 |
| 23 | `RNF114` | 0 | 38 | 4 | 49 | 4 | 233 | 4 | BXPC3:3, HT29:2, MCF7:1 |
| 24 | `SOCS1` | 0 | 32 | 7 | 3 | 0 | 460 | 3 | BXPC3:1, HAP1:2, MCF7:2 |
| 25 | `SP100` | 1 | 65 | 0 | 13 | 1 | 93 | 1 | BXPC3:1, HT29:1, MCF7:1 |
| 26 | `STAT1` | 1 | 264 | 207 | 465 | 64 | 67 | 1 | BXPC3:1, HAP1:1, HT29:1, K562:1, MCF7:1 |
| 27 | `STAT3` | 84 | 15 | 2 | 0 | 2 | 44 | 2 | A549:1, BXPC3:1, MCF7:2 |
| 28 | `STAT4` | 0 | 13 | 0 | 0 | 0 | 20 | 13 | MCF7:13 |
| 29 | `STAT5A` | 0 | 1 | 2 | 3 | 61 | 187 | 1 | HAP1:1, HT29:1, K562:1, MCF7:1 |
| 30 | `TRAFD1` | 3 | 9 | 0 | 4 | 0 | 1140 | 3 | BXPC3:3, HT29:2, MCF7:2 |
| 31 | `TYK2` | 1021 | 1108 | 107 | 278 | 172 | 1027 | 107 | A549:1, BXPC3:1, HT29:1, K562:1, MCF7:1 |
| 32 | `VAV1` | 0 | 6 | 0 | 10 | 2 | 415 | 2 | BXPC3:1, HT29:3, MCF7:1 |
| 33 | `ZBP1` | 2 | 17 | 2 | 2 | 0 | 578 | 2 | BXPC3:1, MCF7:1 |

## 6. Defectos de la construccion, como datos sobre la medicion

| #  | defecto | que implica |
|---|---|---|
| 1 | El maximo del barrido es una **cota inferior** del techo, restringida a la familia de los seis cortes de frecuencia. | El optimo sobre los `2^n` subconjuntos no se calcula. El techo verdadero de esta familia de emisiones es **>=** lo reportado. |
| 2 | La convencion de union vacia mueve la media. | Por eso se dan las DOS columnas. Con la del marcador, una perturbacion vacia en las seis lineas y `S` vacio aporta `1.0`. |
| 3 | Una sola realizacion de la poblacion real (`seed_real` fijo). | El lado real no se remuestrea, asi que no hay barra de error sobre `R(L,h)`. Es el mismo lado real de todos los reportes previos, a proposito. |
| 4 | El recorte (b) usa `p_adj` del propio ajuste BH de cada linea. | Los `p_adj` no son comparables entre lineas (distinto `m` de BH). El recorte ordena DENTRO de cada linea, nunca entre lineas. |
| 5 | `m` es el tamano del mas chico NO VACIO. | Las lineas con ED vacio siguen vacias en (b): el recorte no las rellena. En A549, 15 de 33 siguen aportando `jac = 0`. |
| 6 | El espacio es el de Jiang (34,025 genes), no el del reto (18,533). | El panel probable de BH difiere (`ntc_nulo.md` §3: x0.79 a x0.91). Otra razon para no extrapolar. |
| 7 | `_exclude_own_gene` quita el propio blanco del conjunto. | Es la convencion del marcador, no una eleccion de este reporte. |

## 7. Las entradas

| entrada | valor |
|---|---|
| IFNB | `/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad` |
| panel | `modelo/out/suelo_panel.json` (`panel_comun`, 33) |
| cache recomprobado | `/datos/modelo/out/pds_obj/ed_real.json` |
| marcador | `cell-eval2` 0.20.0, preset/profile `vcc2026` |
| salidas | `modelo/out/venn_ed.md`, `modelo/out/venn_ed.json` |

*Corrida: 2026-11-05 05:01:05 -> 2026-11-05 05:10:24 UTC, duracion 9.3 min (559 s).*
