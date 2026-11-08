# Bucle del techo sobre 6 pliegues -- condicion IFNB

El techo es el MISMO arnes del suelo con la funcion de emision reemplazada: en vez de emitir una prediccion degenerada (400 celulas del control remuestreadas con reemplazo), emite LA OTRA MITAD DE LAS CELULAS REALES de cada perturbacion. Las seis metricas salen de `cell-eval2 run --preset vcc2026 --profile vcc2026` igual que en el suelo, asi que la tabla resta columna a columna contra `bucle_suelo.md`. No hay ninguna formula nueva.

**24 corridas**: 6 pliegues x 2 variantes de control (S compartido / I partido) x 2 direcciones (A como prediccion / B como prediccion). Ninguna se promedia con otra.

## Punto 1 -- lo que se reuso y lo que no, con el motivo

| pieza de suelo_arnes.py | reusada sin modificar | motivo si no |
|---|---|---|
| `load_global_obs`, `load_var_names`, `load_line_block`, `read_rows_csr` | si | - |
| `score_vs_self`, `read_agg_metrics`, `read_avg_score` | si | - |
| `build_fold_adatas` | **no** | Su contrato de emision es `predict_fn(control_pool_X, labels, rng)`: la funcion de emision recibe SOLO el pool de control, nunca las celulas perturbadas, asi que por construccion no puede emitir la otra mitad de las celulas reales. Ademas fija el lado VERDAD (400 celulas sorteadas con `seed_real`) y el techo necesita gobernar LOS DOS lados para correr las dos direcciones. `techo_arnes.build_techo_pair` ensambla los AnnData con las MISMAS reglas (columna obs `target_gene`, mismo var index, `obs_names_make_unique`). |
| `run_cell_eval2` | **no** | Su argv es fijo y no admite banderas extra. El techo necesita dos: `--write-degenes`, para contar n_DE y |union| con la definicion del propio marcador (punto 5), y `--set control_source=pred` para la variante I (punto 3). `techo_arnes.run_cell_eval2_techo` reproduce ese argv exacto (`argv_suelo`) y agrega solo esas banderas. |

Las tres piezas del suelo **no fueron modificadas**. No hay `git` en este arbol, asi que la evidencia es mtime + md5, leidos al escribir este reporte: siguen siendo los de la corrida del suelo (04:16-04:19 UTC), anteriores a la primera linea de `techo_arnes.py`. Las piezas nuevas viven solo en `modelo/src/techo_arnes.py`, `modelo/src/bucle_techo.py` y `modelo/src/techo_reporte.py`.

| archivo | mtime | md5 |
|---|---|---|
| `suelo_arnes.py` | 2026-10-29 04:19:48 UTC | `26ec515ad24885753fe069c3220946ee` |
| `suelo_floor.py` | 2026-10-29 04:16:50 UTC | `9138f651ff0a0f1432ad2363f707bce4` |
| `bucle_suelo.py` | 2026-10-29 04:17:18 UTC | `92c09ac4981d7f27f9443c4805ee2491` |
| `techo_arnes.py` | 2026-10-29 14:04:49 UTC | `f1ee8376895ee9fdd9af1248833ac1e0` |
| `bucle_techo.py` | 2026-10-29 14:22:26 UTC | `479435924fed7f1e76cccde112420152` |
## Punto 3 -- las dos variantes se deciden por `control_source`, no por el h5ad

El preset `vcc2026` trae `control_source: real` (`cell-eval2/src/cell_eval2/configs/vcc2026.yaml`). Bajo ese valor el DE del lado PREDICHO toma su referencia del bloque de control del lado REAL -- de modo que poner un control distinto en el h5ad de prediccion no lo independiza. Por eso:

| variante | control en el h5ad | `control_source` | que comparte | techo |
|---|---|---|---|---|
| S (compartido) | VERDAD = pool de control completo de la linea; PRED = 400 celulas de ese mismo pool con reemplazo (la regla exacta del suelo) | `real` (preset tal cual) | las dos mitades hacen su DE contra el MISMO control: errores correlacionados | **optimista** |
| I (partido) | el control tambien se parte en mitades; cada mitad usa la suya | `pred` (`--set control_source=pred`) | nada: cada mitad trae su propio control | **pesimista** |

Consecuencia que conviene tener presente: **el suelo ya estaba medido con control compartido** (`control_source: real`), asi que la variante S es la unica que comparte convencion de puntuacion con `bucle_suelo.md`. La variante I cambia esa convencion, de modo que su resta contra el suelo mezcla dos efectos (mitad de poblacion + control independiente) y no solo uno.

Esto no es una hipotesis de esta sesion: `cell_eval2/ceiling.py` implementa un techo por self-split disjunto y **fuerza** `control_source="pred"`, documentando que bajo `"real"` ambas mitades calculan su log2FC contra el control de una sola de ellas, lo que comparte el ruido de muestreo entre las dos cantidades cuya concordancia se esta midiendo y **sesga el techo hacia arriba**. Trae la medicion: `lfc_spearman` 0.54 -> 0.74, y `nsig_spearman` paso de NaN (sin techo defendible) a un 0.94 de apariencia confiable. En sus propias palabras, el override "no es opcional".

Ese `ceiling.py` NO se uso como sustituto, por dos razones: aplica la correccion Spearman-Brown `2r/(1+r)` (que el encargo pide no hacer), y su lista verificada `SB_METRICS` contiene solo UNA de las seis metricas vcc2026 (`pds_cosine`): las otras cinco saldrian `NaN`. El arnes del suelo con la emision reemplazada si entrega las seis.

## Punto 4 -- el marcador NO exige un numero fijo de celulas por perturbacion

`cell_eval2.io.validate_pair` exige exactamente tres cosas: eje de genes identico en nombre y orden, conjuntos de etiquetas de perturbacion iguales entre PRED y REAL, y que la etiqueta de control este presente en REAL. **Nada sobre el numero de celulas por etiqueta.**

Lo que hizo el arnes del suelo en ese punto, para copiarlo:

| lado | regla del suelo | profundidad |
|---|---|---|
| VERDAD, perturbaciones | 400 celulas sin reemplazo, con `assert pool >= 400` | 400 |
| VERDAD, control | TODO el pool de control, sin submuestrear | 1379-3228 |
| PREDICCION, toda etiqueta | 400 celulas con reemplazo | 400 |

Es decir: el suelo **ya enviaba profundidades desiguales** entre PRED y REAL (su control real tenia 1379-3228 celulas contra 400 en pred; los totales fueron 14579 vs 13600 filas en A549). Las 400 eran convencion del arnes y cifra del reto, no requisito del marcador.

**Decision, declarada y no silenciosa:** las mitades se emiten a su tamano natural, sin remuestrear ni rellenar a 400. Emitir la mitad tal cual es lo que mide la reproducibilidad de la propia verdad; rellenar a 400 con reemplazo habria inventado celulas en el lado VERDAD y habria sido, de hecho, una correccion por la perdida de potencia -- justo lo que el encargo prohibe.

**Consecuencia, declarada y NO corregida:** el techo se mide a ~la mitad de la profundidad por perturbacion del suelo. La resta columna a columna contra `bucle_suelo.md` mezcla el efecto techo con ese cambio de profundidad. Queda como **cota inferior**.

## Punto 2 -- n de celulas por perturbacion y sus mitades (no son 400)

400 es la cifra del reto, no de Jiang. Los n reales de IFNB y las mitades efectivas que de ahi salen, por linea. `extra` marca a que mitad fue la celula sobrante de un n impar (desempate alternado por ordinal de grupo, no dirigido).

| perturbacion | A549 n (A/B) | BXPC3 n (A/B) | HAP1 n (A/B) | HT29 n (A/B) | K562 n (A/B) | MCF7 n (A/B) |
|---|---|---|---|---|---|---|
| AKT1 | 606 (303/303) | 2499 (1250/1249) | 1014 (507/507) | 1165 (583/582) | 481 (241/240) | 682 (341/341) |
| BATF2 | 617 (308/309) | 712 (356/356) | 1326 (663/663) | 1002 (501/501) | 483 (241/242) | 1203 (601/602) |
| BRD9 | 699 (350/349) | 1676 (838/838) | 975 (488/487) | 996 (498/498) | 692 (346/346) | 796 (398/398) |
| CEBPG | 701 (350/351) | 1814 (907/907) | 1485 (742/743) | 1118 (559/559) | 443 (221/222) | 827 (413/414) |
| DTX3L | 544 (272/272) | 2255 (1128/1127) | 1321 (661/660) | 879 (440/439) | 525 (263/262) | 1786 (893/893) |
| ELK1 | 856 (428/428) | 2174 (1087/1087) | 641 (320/321) | 1015 (507/508) | 695 (347/348) | 2212 (1106/1106) |
| ETS1 | 622 (311/311) | 984 (492/492) | 1044 (522/522) | 1083 (542/541) | 485 (243/242) | 1045 (523/522) |
| ETS2 | 576 (288/288) | 2518 (1259/1259) | 1587 (793/794) | 968 (484/484) | 700 (350/350) | 1368 (684/684) |
| H1F0 | 671 (336/335) | 1170 (585/585) | 1604 (802/802) | 1281 (641/640) | 629 (315/314) | 1128 (564/564) |
| HERC5 | 468 (234/234) | 1066 (533/533) | 995 (497/498) | 1040 (520/520) | 462 (231/231) | 971 (485/486) |
| HERC6 | 440 (220/220) | 1640 (820/820) | 1169 (585/584) | 826 (413/413) | 408 (204/204) | 1281 (641/640) |
| ID2 | 576 (288/288) | 2031 (1015/1016) | 1394 (697/697) | 796 (398/398) | 531 (265/266) | 716 (358/358) |
| IFNAR1 | 634 (317/317) | 1104 (552/552) | 1152 (576/576) | 1148 (574/574) | 460 (230/230) | 644 (322/322) |
| IFNAR2 | 723 (361/362) | 2901 (1450/1451) | 1261 (630/631) | 1117 (558/559) | 439 (219/220) | 450 (225/225) |
| IRF1 | 950 (475/475) | 1447 (724/723) | 1698 (849/849) | 1351 (676/675) | 568 (284/284) | 1803 (902/901) |
| IRF7 | 436 (218/218) | 3860 (1930/1930) | 1328 (664/664) | 1086 (543/543) | 528 (264/264) | 557 (278/279) |
| JAK1 | 693 (347/346) | 1045 (523/522) | 1471 (736/735) | 1692 (846/846) | 587 (294/293) | 1289 (645/644) |
| JUN | 525 (262/263) | 1709 (854/855) | 1604 (802/802) | 787 (393/394) | 503 (251/252) | 1181 (590/591) |
| MAP3K14 | 407 (204/203) | 1734 (867/867) | 982 (491/491) | 980 (490/490) | 429 (215/214) | 1041 (521/520) |
| MAPK8 | 581 (290/291) | 1029 (514/515) | 1385 (692/693) | 999 (499/500) | 509 (254/255) | 1119 (559/560) |
| PARP12 | 682 (341/341) | 1178 (589/589) | 1481 (741/740) | 1062 (531/531) | 647 (324/323) | 1376 (688/688) |
| RAP1GAP | 675 (337/338) | 2832 (1416/1416) | 1733 (866/867) | 1332 (666/666) | 582 (291/291) | 1287 (643/644) |
| RNF114 | 688 (344/344) | 1007 (504/503) | 1738 (869/869) | 1439 (720/719) | 602 (301/301) | 984 (492/492) |
| SOCS1 | 597 (298/299) | 1907 (953/954) | 1526 (763/763) | 1188 (594/594) | 583 (291/292) | 863 (431/432) |
| SP100 | 496 (248/248) | 2057 (1029/1028) | 1234 (617/617) | 1534 (767/767) | 579 (290/289) | 724 (362/362) |
| STAT1 | 630 (315/315) | 1702 (851/851) | 2025 (1012/1013) | 1317 (658/659) | 565 (282/283) | 1277 (638/639) |
| STAT3 | 505 (253/252) | 1354 (677/677) | 1165 (583/582) | 961 (481/480) | 478 (239/239) | 1141 (571/570) |
| STAT4 | 510 (255/255) | 1520 (760/760) | 1369 (684/685) | 855 (427/428) | 466 (233/233) | 855 (427/428) |
| STAT5A | 842 (421/421) | 3004 (1502/1502) | 1793 (897/896) | 1278 (639/639) | 869 (435/434) | 1195 (598/597) |
| TRAFD1 | 478 (239/239) | 1231 (615/616) | 796 (398/398) | 1230 (615/615) | 620 (310/310) | 984 (492/492) |
| TYK2 | 737 (369/368) | 2687 (1344/1343) | 736 (368/368) | 1108 (554/554) | 532 (266/266) | 822 (411/411) |
| VAV1 | 784 (392/392) | 2172 (1086/1086) | 2049 (1024/1025) | 1096 (548/548) | 595 (297/298) | 1137 (568/569) |
| ZBP1 | 509 (255/254) | 1297 (649/648) | 1218 (609/609) | 1236 (618/618) | 482 (241/241) | 1598 (799/799) |
| **NT (control)** | 1379 (689/690) | 3228 (1614/1614) | 3172 (1586/1586) | 2736 (1368/1368) | 1398 (699/699) | 2669 (1334/1335) |

Rango de n por (linea, perturbacion): **407 a 3860** celulas, mediana 1007. Ninguna es 400. La mitad mas delgada del estudio tiene **203 celulas**, frente a las 400 por perturbacion del suelo.

## Las tablas -- una por variante y por direccion (el entregable)

Mismas seis columnas y mismo orden de lineas que `bucle_suelo.md`. `avg_score` se puntua contra si mismo igual que en el suelo, asi que vale 0.0 por construccion y no es informativo (ver nota al final).


### Variante S (control COMPARTIDO) -- direccion: mitad A como PREDICCION

| linea excluida | n control (real/pred) | pds_cosine | mse_capped_norm | lfc_nmae | fid | reach | jac | avg_score |
|---|---|---|---|---|---|---|---|---|
| A549 | 1379/400 | 0.865530 | 0.022814 | 0.444261 | 0.636183 | 0.910526 | 0.403355 | 0.000000 |
| BXPC3 | 3228/400 | 0.996212 | 0.014616 | 0.613665 | 0.822708 | 0.983263 | 0.247759 | 0.000000 |
| HAP1 | 3172/400 | 0.892045 | -0.000847 | 0.334892 | 0.727757 | 0.952498 | 0.338312 | NaN |
| HT29 | 2736/400 | 0.914773 | 0.047808 | 0.594210 | 0.606566 | 0.746267 | 0.105407 | 0.000000 |
| K562 | 1398/400 | 0.855114 | 0.136006 | 0.426234 | 0.671101 | 0.882125 | 0.282639 | 0.000000 |
| MCF7 | 2669/400 | 0.987689 | 0.041080 | 0.601472 | 0.619716 | 0.857840 | 0.136259 | 0.000000 |

### Variante S (control COMPARTIDO) -- direccion: mitad B como PREDICCION

| linea excluida | n control (real/pred) | pds_cosine | mse_capped_norm | lfc_nmae | fid | reach | jac | avg_score |
|---|---|---|---|---|---|---|---|---|
| A549 | 1379/400 | 0.835227 | 0.086362 | 0.396085 | 0.639951 | 0.962963 | 0.403355 | 0.000000 |
| BXPC3 | 3228/400 | 0.997159 | 0.021209 | 0.521412 | 0.842348 | 0.991702 | 0.247759 | 0.000000 |
| HAP1 | 3172/400 | 0.903409 | 0.024103 | 0.317906 | 0.679394 | 0.921664 | 0.338312 | 0.000000 |
| HT29 | 2736/400 | 0.920455 | 0.014729 | 0.588089 | 0.538807 | 0.820564 | 0.105407 | 0.000000 |
| K562 | 1398/400 | 0.831439 | 0.145790 | 0.417848 | 0.548864 | 0.920355 | 0.282639 | 0.000000 |
| MCF7 | 2669/400 | 0.984848 | 0.053387 | 0.605358 | 0.569063 | 0.861717 | 0.136259 | 0.000000 |

### Variante I (control PARTIDO) -- direccion: mitad A como PREDICCION

| linea excluida | n control (real/pred) | pds_cosine | mse_capped_norm | lfc_nmae | fid | reach | jac | avg_score |
|---|---|---|---|---|---|---|---|---|
| A549 | 690/689 | 0.841856 | 0.021735 | 0.492830 | 0.681106 | 0.949497 | 0.562245 | 0.000000 |
| BXPC3 | 1614/1614 | 0.997159 | 0.014231 | 0.707090 | 0.790451 | 0.968083 | 0.211960 | 0.000000 |
| HAP1 | 1586/1586 | 0.928977 | -0.000803 | 0.414945 | 0.621897 | 0.890658 | 0.318423 | NaN |
| HT29 | 1368/1368 | 0.908144 | 0.046899 | 0.553540 | 0.530503 | 0.540714 | 0.133187 | 0.000000 |
| K562 | 699/699 | 0.856061 | 0.129427 | 0.345512 | 0.636055 | 0.823956 | 0.403280 | 0.000000 |
| MCF7 | 1335/1334 | 0.988636 | 0.041203 | 0.613280 | 0.684198 | 0.861142 | 0.095610 | 0.000000 |

### Variante I (control PARTIDO) -- direccion: mitad B como PREDICCION

| linea excluida | n control (real/pred) | pds_cosine | mse_capped_norm | lfc_nmae | fid | reach | jac | avg_score |
|---|---|---|---|---|---|---|---|---|
| A549 | 689/690 | 0.782197 | 0.094039 | 0.473392 | 0.748485 | 0.922787 | 0.562245 | 0.000000 |
| BXPC3 | 1614/1614 | 0.991477 | 0.021559 | 0.628513 | 0.778034 | 0.943964 | 0.211960 | 0.000000 |
| HAP1 | 1586/1586 | 0.910985 | 0.025049 | 0.362917 | 0.646667 | 0.845436 | 0.318423 | 0.000000 |
| HT29 | 1368/1368 | 0.840909 | 0.015348 | 0.642498 | 0.407195 | 0.704372 | 0.133187 | 0.000000 |
| K562 | 699/699 | 0.809659 | 0.153676 | 0.518946 | 0.519768 | 0.743649 | 0.403280 | 0.000000 |
| MCF7 | 1334/1335 | 0.985795 | 0.053492 | 0.776169 | 0.350310 | 0.673735 | 0.095610 | 0.000000 |

## La resta columna a columna contra el suelo

`techo - suelo` por metrica y por pliegue. El signo util depende del sentido de cada metrica: en `mse_capped_norm` y `lfc_nmae` mas bajo es mejor, en las otras cuatro mas alto es mejor. La columna `mejor` dice si el techo quedo del lado bueno del suelo.


### Variante S, direccion A: techo - suelo

| linea | pds_cosine | mse_capped_norm | lfc_nmae | fid | reach | jac | mejor en |
|---|---|---|---|---|---|---|---|
| A549 | +0.391098 | -2.539358 | -0.591757 | +0.636183 | +0.765987 | -0.051190 | 5/6 |
| BXPC3 | +0.478220 | -1.318144 | -0.383300 | +0.822682 | +0.840374 | +0.247759 | 6/6 |
| HAP1 | +0.379735 | -1.869574 | -0.665929 | +0.727757 | +0.760090 | -0.085931 | 5/6 |
| HT29 | +0.424242 | -1.883758 | -0.450104 | +0.597307 | +0.530234 | -0.076411 | 5/6 |
| K562 | +0.359848 | -2.777221 | -0.602612 | +0.671101 | +0.681677 | +0.040215 | 6/6 |
| MCF7 | +0.475379 | -1.733794 | -0.455336 | +0.619716 | +0.844278 | +0.136259 | 6/6 |

### Variante S, direccion B: techo - suelo

| linea | pds_cosine | mse_capped_norm | lfc_nmae | fid | reach | jac | mejor en |
|---|---|---|---|---|---|---|---|
| A549 | +0.360795 | -2.475810 | -0.639933 | +0.639951 | +0.818423 | -0.051190 | 5/6 |
| BXPC3 | +0.479167 | -1.311550 | -0.475554 | +0.842322 | +0.848813 | +0.247759 | 6/6 |
| HAP1 | +0.391098 | -1.844625 | -0.682915 | +0.679394 | +0.729256 | -0.085931 | 5/6 |
| HT29 | +0.429924 | -1.916837 | -0.456225 | +0.529548 | +0.604531 | -0.076411 | 5/6 |
| K562 | +0.336174 | -2.767437 | -0.610999 | +0.548864 | +0.719906 | +0.040215 | 6/6 |
| MCF7 | +0.472538 | -1.721487 | -0.451450 | +0.569063 | +0.848155 | +0.136259 | 6/6 |

### Variante I, direccion A: techo - suelo

| linea | pds_cosine | mse_capped_norm | lfc_nmae | fid | reach | jac | mejor en |
|---|---|---|---|---|---|---|---|
| A549 | +0.367424 | -2.540438 | -0.543189 | +0.681106 | +0.804957 | +0.107700 | 6/6 |
| BXPC3 | +0.479167 | -1.318529 | -0.289876 | +0.790425 | +0.825194 | +0.211960 | 6/6 |
| HAP1 | +0.416667 | -1.869531 | -0.585876 | +0.621897 | +0.698250 | -0.105819 | 5/6 |
| HT29 | +0.417614 | -1.884666 | -0.490774 | +0.521243 | +0.324682 | -0.048631 | 5/6 |
| K562 | +0.360795 | -2.783800 | -0.683334 | +0.636055 | +0.623508 | +0.160856 | 6/6 |
| MCF7 | +0.476326 | -1.733670 | -0.443528 | +0.684198 | +0.847579 | +0.095610 | 6/6 |

### Variante I, direccion B: techo - suelo

| linea | pds_cosine | mse_capped_norm | lfc_nmae | fid | reach | jac | mejor en |
|---|---|---|---|---|---|---|---|
| A549 | +0.307765 | -2.468133 | -0.562626 | +0.748485 | +0.778248 | +0.107700 | 6/6 |
| BXPC3 | +0.473485 | -1.311201 | -0.368452 | +0.778008 | +0.801075 | +0.211960 | 6/6 |
| HAP1 | +0.398674 | -1.843678 | -0.637904 | +0.646667 | +0.653028 | -0.105819 | 5/6 |
| HT29 | +0.350379 | -1.916217 | -0.401816 | +0.397936 | +0.488340 | -0.048631 | 5/6 |
| K562 | +0.314394 | -2.759551 | -0.509900 | +0.519768 | +0.543201 | +0.160856 | 6/6 |
| MCF7 | +0.473485 | -1.721382 | -0.280639 | +0.350310 | +0.660172 | +0.095610 | 6/6 |

## El denominador `(r - b)`: donde la escala funciona y donde NO

La escala `(u - b) / (r - b)` solo tiene sentido si el techo `r` esta del lado BUENO del suelo `b`. Aqui el `span` esta orientado por el sentido de cada metrica: `span = r - b` donde mas alto es mejor, `span = b - r` donde mas bajo es mejor. **Un `span` negativo significa que el techo quedo peor que el suelo y que esa metrica no admite esta escala.**

| metrica | var | dir | spans por linea (A549,BXPC3,HAP1,HT29,K562,MCF7) | n negativos |
|---|---|---|---|---|
| pds_cosine | S | A | +0.3911, +0.4782, +0.3797, +0.4242, +0.3598, +0.4754 | 0/6 |
| pds_cosine | S | B | +0.3608, +0.4792, +0.3911, +0.4299, +0.3362, +0.4725 | 0/6 |
| pds_cosine | I | A | +0.3674, +0.4792, +0.4167, +0.4176, +0.3608, +0.4763 | 0/6 |
| pds_cosine | I | B | +0.3078, +0.4735, +0.3987, +0.3504, +0.3144, +0.4735 | 0/6 |
| mse_capped_norm | S | A | +2.5394, +1.3181, +1.8696, +1.8838, +2.7772, +1.7338 | 0/6 |
| mse_capped_norm | S | B | +2.4758, +1.3116, +1.8446, +1.9168, +2.7674, +1.7215 | 0/6 |
| mse_capped_norm | I | A | +2.5404, +1.3185, +1.8695, +1.8847, +2.7838, +1.7337 | 0/6 |
| mse_capped_norm | I | B | +2.4681, +1.3112, +1.8437, +1.9162, +2.7596, +1.7214 | 0/6 |
| lfc_nmae | S | A | +0.5918, +0.3833, +0.6659, +0.4501, +0.6026, +0.4553 | 0/6 |
| lfc_nmae | S | B | +0.6399, +0.4756, +0.6829, +0.4562, +0.6110, +0.4515 | 0/6 |
| lfc_nmae | I | A | +0.5432, +0.2899, +0.5859, +0.4908, +0.6833, +0.4435 | 0/6 |
| lfc_nmae | I | B | +0.5626, +0.3685, +0.6379, +0.4018, +0.5099, +0.2806 | 0/6 |
| fid | S | A | +0.6362, +0.8227, +0.7278, +0.5973, +0.6711, +0.6197 | 0/6 |
| fid | S | B | +0.6400, +0.8423, +0.6794, +0.5295, +0.5489, +0.5691 | 0/6 |
| fid | I | A | +0.6811, +0.7904, +0.6219, +0.5212, +0.6361, +0.6842 | 0/6 |
| fid | I | B | +0.7485, +0.7780, +0.6467, +0.3979, +0.5198, +0.3503 | 0/6 |
| reach | S | A | +0.7660, +0.8404, +0.7601, +0.5302, +0.6817, +0.8443 | 0/6 |
| reach | S | B | +0.8184, +0.8488, +0.7293, +0.6045, +0.7199, +0.8482 | 0/6 |
| reach | I | A | +0.8050, +0.8252, +0.6982, +0.3247, +0.6235, +0.8476 | 0/6 |
| reach | I | B | +0.7782, +0.8011, +0.6530, +0.4883, +0.5432, +0.6602 | 0/6 |
| jac | S | A | -0.0512, +0.2478, -0.0859, -0.0764, +0.0402, +0.1363 | 3/6 |
| jac | S | B | -0.0512, +0.2478, -0.0859, -0.0764, +0.0402, +0.1363 | 3/6 |
| jac | I | A | +0.1077, +0.2120, -0.1058, -0.0486, +0.1609, +0.0956 | 2/6 |
| jac | I | B | +0.1077, +0.2120, -0.1058, -0.0486, +0.1609, +0.0956 | 2/6 |

**Lectura.**

- Metricas donde el techo quedo del lado bueno del suelo en TODAS las lineas, variantes y direcciones (la escala funciona): `pds_cosine`, `mse_capped_norm`, `lfc_nmae`, `fid`, `reach`.

- Metricas con al menos un `span` negativo o nulo (la escala NO funciona ahi): `jac`.

- El caso de `jac` no es un fallo de medicion sino el artefacto D-023 operando en el SUELO: la prediccion degenerada del suelo no llama nada significativo, de modo que toda perturbacion cuya mitad real tampoco llama nada queda con union vacia y cobra un 1.0 regalado. Ese 1.0 inflaba el suelo. El techo, con senal real en los dos lados, tiene uniones NO vacias y jaccards genuinos menores que 1 -- y sale por debajo. **El suelo de `jac` no es un suelo: es un artefacto de union vacia.** Mientras eso no se arregle, `jac` no admite la escala `(u - b)/(r - b)`, y la tabla de uniones vacias de mas abajo es el lugar donde se ve.


## El intervalo del techo: entre S y I (no se elige una)

El techo verdadero esta ENTRE la variante optimista (S, control compartido) y la pesimista (I, control partido). Aqui las dos, con el ancho del intervalo. Se promedia sobre las dos DIRECCIONES solo para dar el ancho; las direcciones crudas estan arriba, sin promediar.

| linea | metrica | S dir A | S dir B | I dir A | I dir B | intervalo S..I | ancho |
|---|---|---|---|---|---|---|---|
| A549 | pds_cosine | 0.865530 | 0.835227 | 0.841856 | 0.782197 | 0.812027 .. 0.850379 | 0.038352 |
| A549 | mse_capped_norm | 0.022814 | 0.086362 | 0.021735 | 0.094039 | 0.054588 .. 0.057887 | 0.003299 |
| A549 | lfc_nmae | 0.444261 | 0.396085 | 0.492830 | 0.473392 | 0.420173 .. 0.483111 | 0.062938 |
| A549 | fid | 0.636183 | 0.639951 | 0.681106 | 0.748485 | 0.638067 .. 0.714795 | 0.076729 |
| A549 | reach | 0.910526 | 0.962963 | 0.949497 | 0.922787 | 0.936142 .. 0.936745 | 0.000603 |
| A549 | jac | 0.403355 | 0.403355 | 0.562245 | 0.562245 | 0.403355 .. 0.562245 | 0.158891 |
| BXPC3 | pds_cosine | 0.996212 | 0.997159 | 0.997159 | 0.991477 | 0.994318 .. 0.996686 | 0.002367 |
| BXPC3 | mse_capped_norm | 0.014616 | 0.021209 | 0.014231 | 0.021559 | 0.017895 .. 0.017912 | 0.000018 |
| BXPC3 | lfc_nmae | 0.613665 | 0.521412 | 0.707090 | 0.628513 | 0.567538 .. 0.667801 | 0.100263 |
| BXPC3 | fid | 0.822708 | 0.842348 | 0.790451 | 0.778034 | 0.784243 .. 0.832528 | 0.048285 |
| BXPC3 | reach | 0.983263 | 0.991702 | 0.968083 | 0.943964 | 0.956023 .. 0.987482 | 0.031459 |
| BXPC3 | jac | 0.247759 | 0.247759 | 0.211960 | 0.211960 | 0.211960 .. 0.247759 | 0.035799 |
| HAP1 | pds_cosine | 0.892045 | 0.903409 | 0.928977 | 0.910985 | 0.897727 .. 0.919981 | 0.022254 |
| HAP1 | mse_capped_norm | -0.000847 | 0.024103 | -0.000803 | 0.025049 | 0.011628 .. 0.012123 | 0.000495 |
| HAP1 | lfc_nmae | 0.334892 | 0.317906 | 0.414945 | 0.362917 | 0.326399 .. 0.388931 | 0.062532 |
| HAP1 | fid | 0.727757 | 0.679394 | 0.621897 | 0.646667 | 0.634282 .. 0.703575 | 0.069293 |
| HAP1 | reach | 0.952498 | 0.921664 | 0.890658 | 0.845436 | 0.868047 .. 0.937081 | 0.069034 |
| HAP1 | jac | 0.338312 | 0.338312 | 0.318423 | 0.318423 | 0.318423 .. 0.338312 | 0.019888 |
| HT29 | pds_cosine | 0.914773 | 0.920455 | 0.908144 | 0.840909 | 0.874527 .. 0.917614 | 0.043087 |
| HT29 | mse_capped_norm | 0.047808 | 0.014729 | 0.046899 | 0.015348 | 0.031124 .. 0.031268 | 0.000145 |
| HT29 | lfc_nmae | 0.594210 | 0.588089 | 0.553540 | 0.642498 | 0.591150 .. 0.598019 | 0.006869 |
| HT29 | fid | 0.606566 | 0.538807 | 0.530503 | 0.407195 | 0.468849 .. 0.572687 | 0.103838 |
| HT29 | reach | 0.746267 | 0.820564 | 0.540714 | 0.704372 | 0.622543 .. 0.783415 | 0.160872 |
| HT29 | jac | 0.105407 | 0.105407 | 0.133187 | 0.133187 | 0.105407 .. 0.133187 | 0.027781 |
| K562 | pds_cosine | 0.855114 | 0.831439 | 0.856061 | 0.809659 | 0.832860 .. 0.843277 | 0.010417 |
| K562 | mse_capped_norm | 0.136006 | 0.145790 | 0.129427 | 0.153676 | 0.140898 .. 0.141551 | 0.000654 |
| K562 | lfc_nmae | 0.426234 | 0.417848 | 0.345512 | 0.518946 | 0.422041 .. 0.432229 | 0.010189 |
| K562 | fid | 0.671101 | 0.548864 | 0.636055 | 0.519768 | 0.577911 .. 0.609983 | 0.032071 |
| K562 | reach | 0.882125 | 0.920355 | 0.823956 | 0.743649 | 0.783803 .. 0.901240 | 0.117437 |
| K562 | jac | 0.282639 | 0.282639 | 0.403280 | 0.403280 | 0.282639 .. 0.403280 | 0.120641 |
| MCF7 | pds_cosine | 0.987689 | 0.984848 | 0.988636 | 0.985795 | 0.986269 .. 0.987216 | 0.000947 |
| MCF7 | mse_capped_norm | 0.041080 | 0.053387 | 0.041203 | 0.053492 | 0.047233 .. 0.047347 | 0.000114 |
| MCF7 | lfc_nmae | 0.601472 | 0.605358 | 0.613280 | 0.776169 | 0.603415 .. 0.694724 | 0.091309 |
| MCF7 | fid | 0.619716 | 0.569063 | 0.684198 | 0.350310 | 0.517254 .. 0.594389 | 0.077135 |
| MCF7 | reach | 0.857840 | 0.861717 | 0.861142 | 0.673735 | 0.767438 .. 0.859779 | 0.092340 |
| MCF7 | jac | 0.136259 | 0.136259 | 0.095610 | 0.095610 | 0.095610 .. 0.136259 | 0.040649 |

### El suelo de `jac` medido por dentro: es la tasa de union vacia

La afirmacion de arriba no se infiere, se cuenta. En `results.csv` de cada pliegue del suelo, el jaccard por perturbacion toma **solo dos valores: exactamente 0.0 o exactamente 1.0**, nunca uno intermedio (se verifico en los seis pliegues, 33 de 33 perturbaciones cada uno). Con una prediccion degenerada el conjunto predicho esta vacio, asi que un 1.0 solo puede venir de `|R u P| = 0` y un 0.0 significa que el lado real si tenia senal. De modo que la columna `jac` del suelo **es la fraccion de perturbaciones sin ninguna senal de DE**, y no mide nada de la prediccion.

| linea | suelo: perts jac=1.0 (union vacia) | suelo: perts jac=0.0 | suelo: otros valores | suelo `jac` = 1.0s/33 | techo S: perts union vacia | techo S `jac` | span (techo-suelo) |
|---|---|---|---|---|---|---|---|
| A549 | 15 | 18 | 0 | 15/33 = 0.454545 | 10 | 0.403355 | -0.051190 |
| BXPC3 | 0 | 33 | 0 | 0/33 = 0.000000 | 0 | 0.247759 | +0.247759 |
| HAP1 | 14 | 19 | 0 | 14/33 = 0.424242 | 4 | 0.338312 | -0.085931 |
| HT29 | 6 | 27 | 0 | 6/33 = 0.181818 | 1 | 0.105407 | -0.076411 |
| K562 | 8 | 25 | 0 | 8/33 = 0.242424 | 6 | 0.282639 | +0.040215 |
| MCF7 | 0 | 33 | 0 | 0/33 = 0.000000 | 0 | 0.136259 | +0.136259 |

La quinta columna reproduce exactamente la columna `jac` de `bucle_suelo.md` como un simple conteo de uniones vacias dividido por 33. Donde el suelo tiene 0 uniones vacias (BXPC3 y MCF7) el techo queda limpiamente por encima; donde el suelo cobraba mas 1.0 regalados que el techo (A549 15 vs 10, HAP1 14 vs 4, HT29 6 vs 1) el span sale negativo. **Eso es todo lo que significa el span negativo de `jac`: no que la verdad sea menos reproducible que una prediccion degenerada, sino que el suelo de esa metrica estaba midiendo la tasa de perturbaciones sin senal.**


### ¿Es S realmente la optimista? (se comprueba, no se asume)

Por metrica, en cuantas de las 6 lineas la variante S quedo del lado BUENO de la variante I (promediando las dos direcciones solo para esta comparacion). El encargo esperaba S optimista en las seis.

| metrica | S mejor que I en | S peor que I en | veredicto |
|---|---|---|---|
| pds_cosine | 4/6 | 2/6 | mezclado, no domina ninguna |
| mse_capped_norm | 4/6 | 2/6 | mezclado, no domina ninguna |
| lfc_nmae | 6/6 | 0/6 | S optimista, como se esperaba |
| fid | 5/6 | 1/6 | mezclado, no domina ninguna |
| reach | 6/6 | 0/6 | S optimista, como se esperaba |
| jac | 3/6 | 3/6 | mezclado, no domina ninguna |

Donde S no sale optimista, el sospechoso es otra vez la union vacia: la variante I parte el control, de modo que cada mitad hace su DE con la mitad de las celulas de referencia, detecta MENOS genes significativos, y produce MAS uniones vacias -- cada una un 1.0 regalado en `jac`. El efecto es un techo mas alto por perdida de potencia, no por mejor concordancia. La tabla de uniones vacias por variante, mas abajo, es la que lo distingue.


## Punto 5 -- por perturbacion: n_DE mitad A, n_DE mitad B, |union|

`de_wilcoxon_sig_jaccard` devuelve **1.0 cuando la union es vacia** (convencion `J(vacio, vacio) = 1`, D-023). Sin la columna `|union|` a la vista, una perturbacion sin ninguna senal aparece con techo perfecto y es invisible. Las cuentas se calculan con la definicion del propio marcador: `p_adj < 0.05`, `.unique()` sobre (target, feature) y el anti-join que excluye el gen propio de cada target en LOS DOS lados (issue #172), llamando a `cell_eval2.de.prepare_de` sobre las tablas que escribio `--write-degenes`.

`n_DE_mitadA` y `n_DE_mitadB` son invariantes a la direccion dentro de una variante (cada mitad hace su DE contra el mismo control en las dos direcciones), asi que se tabulan una vez por variante.

### Resumen: perturbaciones con union vacia (jaccard=1.0 regalado)

La ultima columna es un **diagnostico, no una metrica**: es la media de los jaccards por perturbacion dejando fuera las de union vacia, calculada aqui para medir cuanto del numero reportado es artefacto. NO es comparable contra las columnas de `bucle_suelo.md` y no sustituye a lo que agrega el marcador.

| variante | linea | n perts con \|union\|=0 | de 33 | jac reportado (marcador) | jac sin esas perts (diagnostico) |
|---|---|---|---|---|---|
| S | A549 | 10 | 33 | 0.403355 | 0.143944 |
| S | BXPC3 | 0 | 33 | 0.247759 | 0.247759 |
| S | HAP1 | 4 | 33 | 0.338312 | 0.247044 |
| S | HT29 | 1 | 33 | 0.105407 | 0.077451 |
| S | K562 | 6 | 33 | 0.282639 | 0.123225 |
| S | MCF7 | 0 | 33 | 0.136259 | 0.136259 |
| I | A549 | 16 | 33 | 0.562245 | 0.150241 |
| I | BXPC3 | 0 | 33 | 0.211960 | 0.211960 |
| I | HAP1 | 4 | 33 | 0.318423 | 0.224413 |
| I | HT29 | 2 | 33 | 0.133187 | 0.077264 |
| I | K562 | 9 | 33 | 0.403280 | 0.179510 |
| I | MCF7 | 0 | 33 | 0.095610 | 0.095610 |

### Variante S -- detalle por (linea, perturbacion)

| linea | perturbacion | n celulas | n mitad A | n mitad B | n_DE_mitadA | n_DE_mitadB | interseccion | **\|union\|** | jac marcador |
|---|---|---|---|---|---|---|---|---|---|
| A549 | AKT1 | 606 | 303 | 303 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | BATF2 | 617 | 308 | 309 | 0 | 1 | 0 | 1 | 0.0000 |
| A549 | BRD9 | 699 | 350 | 349 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | CEBPG | 701 | 350 | 351 | 1 | 2 | 0 | 3 | 0.0000 |
| A549 | DTX3L | 544 | 272 | 272 | 2 | 1 | 1 | 2 | 0.5000 |
| A549 | ELK1 | 856 | 428 | 428 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | ETS1 | 622 | 311 | 311 | 0 | 1 | 0 | 1 | 0.0000 |
| A549 | ETS2 | 576 | 288 | 288 | 9 | 19 | 2 | 26 | 0.0769 |
| A549 | H1F0 | 671 | 336 | 335 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | HERC5 | 468 | 234 | 234 | 0 | 1 | 0 | 1 | 0.0000 |
| A549 | HERC6 | 440 | 220 | 220 | 3 | 0 | 0 | 3 | 0.0000 |
| A549 | ID2 | 576 | 288 | 288 | 0 | 1 | 0 | 1 | 0.0000 |
| A549 | IFNAR1 | 634 | 317 | 317 | 348 | 392 | 157 | 583 | 0.2693 |
| A549 | IFNAR2 | 723 | 361 | 362 | 589 | 697 | 303 | 983 | 0.3082 |
| A549 | IRF1 | 950 | 475 | 475 | 7 | 8 | 4 | 11 | 0.3636 |
| A549 | IRF7 | 436 | 218 | 218 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | JAK1 | 693 | 347 | 346 | 84 | 70 | 24 | 130 | 0.1846 |
| A549 | JUN | 525 | 262 | 263 | 2 | 7 | 2 | 7 | 0.2857 |
| A549 | MAP3K14 | 407 | 204 | 203 | 5 | 5 | 2 | 8 | 0.2500 |
| A549 | MAPK8 | 581 | 290 | 291 | 3 | 1 | 0 | 4 | 0.0000 |
| A549 | PARP12 | 682 | 341 | 341 | 1 | 0 | 0 | 1 | 0.0000 |
| A549 | RAP1GAP | 675 | 337 | 338 | 1 | 0 | 0 | 1 | 0.0000 |
| A549 | RNF114 | 688 | 344 | 344 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | SOCS1 | 597 | 298 | 299 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | SP100 | 496 | 248 | 248 | 0 | 2 | 0 | 2 | 0.0000 |
| A549 | STAT1 | 630 | 315 | 315 | 3 | 2 | 1 | 4 | 0.2500 |
| A549 | STAT3 | 505 | 253 | 252 | 34 | 54 | 11 | 77 | 0.1429 |
| A549 | STAT4 | 510 | 255 | 255 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | STAT5A | 842 | 421 | 421 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | TRAFD1 | 478 | 239 | 239 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | TYK2 | 737 | 369 | 368 | 882 | 946 | 470 | 1358 | 0.3461 |
| A549 | VAV1 | 784 | 392 | 392 | 1 | 0 | 0 | 1 | 0.0000 |
| A549 | ZBP1 | 509 | 255 | 254 | 3 | 1 | 1 | 3 | 0.3333 |
| BXPC3 | AKT1 | 2499 | 1250 | 1249 | 31 | 26 | 12 | 45 | 0.2667 |
| BXPC3 | BATF2 | 712 | 356 | 356 | 36 | 43 | 11 | 68 | 0.1618 |
| BXPC3 | BRD9 | 1676 | 838 | 838 | 262 | 133 | 78 | 317 | 0.2461 |
| BXPC3 | CEBPG | 1814 | 907 | 907 | 157 | 148 | 49 | 256 | 0.1914 |
| BXPC3 | DTX3L | 2255 | 1128 | 1127 | 83 | 94 | 38 | 139 | 0.2734 |
| BXPC3 | ELK1 | 2174 | 1087 | 1087 | 39 | 66 | 19 | 86 | 0.2209 |
| BXPC3 | ETS1 | 984 | 492 | 492 | 454 | 540 | 277 | 717 | 0.3863 |
| BXPC3 | ETS2 | 2518 | 1259 | 1259 | 158 | 112 | 40 | 230 | 0.1739 |
| BXPC3 | H1F0 | 1170 | 585 | 585 | 16 | 52 | 6 | 62 | 0.0968 |
| BXPC3 | HERC5 | 1066 | 533 | 533 | 7 | 15 | 4 | 18 | 0.2222 |
| BXPC3 | HERC6 | 1640 | 820 | 820 | 251 | 176 | 92 | 335 | 0.2746 |
| BXPC3 | ID2 | 2031 | 1015 | 1016 | 498 | 573 | 205 | 866 | 0.2367 |
| BXPC3 | IFNAR1 | 1104 | 552 | 552 | 815 | 558 | 324 | 1049 | 0.3089 |
| BXPC3 | IFNAR2 | 2901 | 1450 | 1451 | 2992 | 3050 | 2159 | 3883 | 0.5560 |
| BXPC3 | IRF1 | 1447 | 724 | 723 | 131 | 138 | 52 | 217 | 0.2396 |
| BXPC3 | IRF7 | 3860 | 1930 | 1930 | 153 | 87 | 51 | 189 | 0.2698 |
| BXPC3 | JAK1 | 1045 | 523 | 522 | 714 | 579 | 252 | 1041 | 0.2421 |
| BXPC3 | JUN | 1709 | 854 | 855 | 126 | 124 | 56 | 194 | 0.2887 |
| BXPC3 | MAP3K14 | 1734 | 867 | 867 | 136 | 137 | 54 | 219 | 0.2466 |
| BXPC3 | MAPK8 | 1029 | 514 | 515 | 99 | 112 | 25 | 186 | 0.1344 |
| BXPC3 | PARP12 | 1178 | 589 | 589 | 14 | 13 | 4 | 23 | 0.1739 |
| BXPC3 | RAP1GAP | 2832 | 1416 | 1416 | 40 | 40 | 17 | 63 | 0.2698 |
| BXPC3 | RNF114 | 1007 | 504 | 503 | 39 | 56 | 8 | 87 | 0.0920 |
| BXPC3 | SOCS1 | 1907 | 953 | 954 | 190 | 157 | 73 | 274 | 0.2664 |
| BXPC3 | SP100 | 2057 | 1029 | 1028 | 363 | 212 | 112 | 463 | 0.2419 |
| BXPC3 | STAT1 | 1702 | 851 | 851 | 566 | 599 | 298 | 867 | 0.3437 |
| BXPC3 | STAT3 | 1354 | 677 | 677 | 41 | 31 | 12 | 60 | 0.2000 |
| BXPC3 | STAT4 | 1520 | 760 | 760 | 32 | 87 | 13 | 106 | 0.1226 |
| BXPC3 | STAT5A | 3004 | 1502 | 1502 | 13 | 10 | 4 | 19 | 0.2105 |
| BXPC3 | TRAFD1 | 1231 | 615 | 616 | 19 | 37 | 3 | 53 | 0.0566 |
| BXPC3 | TYK2 | 2687 | 1344 | 1343 | 2945 | 2630 | 1920 | 3655 | 0.5253 |
| BXPC3 | VAV1 | 2172 | 1086 | 1086 | 29 | 39 | 13 | 55 | 0.2364 |
| BXPC3 | ZBP1 | 1297 | 649 | 648 | 15 | 20 | 10 | 25 | 0.4000 |
| HAP1 | AKT1 | 1014 | 507 | 507 | 1 | 0 | 0 | 1 | 0.0000 |
| HAP1 | BATF2 | 1326 | 663 | 663 | 18 | 8 | 7 | 19 | 0.3684 |
| HAP1 | BRD9 | 975 | 488 | 487 | 2 | 0 | 0 | 2 | 0.0000 |
| HAP1 | CEBPG | 1485 | 742 | 743 | 85 | 113 | 25 | 173 | 0.1445 |
| HAP1 | DTX3L | 1321 | 661 | 660 | 2 | 8 | 2 | 8 | 0.2500 |
| HAP1 | ELK1 | 641 | 320 | 321 | 11 | 10 | 3 | 18 | 0.1667 |
| HAP1 | ETS1 | 1044 | 522 | 522 | 7 | 8 | 3 | 12 | 0.2500 |
| HAP1 | ETS2 | 1587 | 793 | 794 | 3 | 1 | 0 | 4 | 0.0000 |
| HAP1 | H1F0 | 1604 | 802 | 802 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| HAP1 | HERC5 | 995 | 497 | 498 | 1 | 4 | 0 | 5 | 0.0000 |
| HAP1 | HERC6 | 1169 | 585 | 584 | 13 | 6 | 4 | 15 | 0.2667 |
| HAP1 | ID2 | 1394 | 697 | 697 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| HAP1 | IFNAR1 | 1152 | 576 | 576 | 231 | 238 | 160 | 309 | 0.5178 |
| HAP1 | IFNAR2 | 1261 | 630 | 631 | 291 | 285 | 178 | 398 | 0.4472 |
| HAP1 | IRF1 | 1698 | 849 | 849 | 18 | 17 | 6 | 29 | 0.2069 |
| HAP1 | IRF7 | 1328 | 664 | 664 | 4 | 5 | 2 | 7 | 0.2857 |
| HAP1 | JAK1 | 1471 | 736 | 735 | 98 | 115 | 75 | 138 | 0.5435 |
| HAP1 | JUN | 1604 | 802 | 802 | 0 | 1 | 0 | 1 | 0.0000 |
| HAP1 | MAP3K14 | 982 | 491 | 491 | 1 | 0 | 0 | 1 | 0.0000 |
| HAP1 | MAPK8 | 1385 | 692 | 693 | 0 | 2 | 0 | 2 | 0.0000 |
| HAP1 | PARP12 | 1481 | 741 | 740 | 6 | 8 | 4 | 10 | 0.4000 |
| HAP1 | RAP1GAP | 1733 | 866 | 867 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| HAP1 | RNF114 | 1738 | 869 | 869 | 4 | 5 | 4 | 5 | 0.8000 |
| HAP1 | SOCS1 | 1526 | 763 | 763 | 16 | 33 | 8 | 41 | 0.1951 |
| HAP1 | SP100 | 1234 | 617 | 617 | 1 | 1 | 0 | 2 | 0.0000 |
| HAP1 | STAT1 | 2025 | 1012 | 1013 | 415 | 490 | 257 | 648 | 0.3966 |
| HAP1 | STAT3 | 1165 | 583 | 582 | 3 | 1 | 1 | 3 | 0.3333 |
| HAP1 | STAT4 | 1369 | 684 | 685 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| HAP1 | STAT5A | 1793 | 897 | 896 | 9 | 6 | 5 | 10 | 0.5000 |
| HAP1 | TRAFD1 | 796 | 398 | 398 | 1 | 0 | 0 | 1 | 0.0000 |
| HAP1 | TYK2 | 736 | 368 | 368 | 115 | 119 | 87 | 147 | 0.5918 |
| HAP1 | VAV1 | 2049 | 1024 | 1025 | 1 | 0 | 0 | 1 | 0.0000 |
| HAP1 | ZBP1 | 1218 | 609 | 609 | 1 | 2 | 1 | 2 | 0.5000 |
| HT29 | AKT1 | 1165 | 583 | 582 | 4 | 7 | 0 | 11 | 0.0000 |
| HT29 | BATF2 | 1002 | 501 | 501 | 3 | 4 | 0 | 7 | 0.0000 |
| HT29 | BRD9 | 996 | 498 | 498 | 3 | 24 | 0 | 27 | 0.0000 |
| HT29 | CEBPG | 1118 | 559 | 559 | 25 | 22 | 3 | 44 | 0.0682 |
| HT29 | DTX3L | 879 | 440 | 439 | 15 | 12 | 1 | 26 | 0.0385 |
| HT29 | ELK1 | 1015 | 507 | 508 | 2 | 10 | 0 | 12 | 0.0000 |
| HT29 | ETS1 | 1083 | 542 | 541 | 4 | 3 | 0 | 7 | 0.0000 |
| HT29 | ETS2 | 968 | 484 | 484 | 179 | 203 | 51 | 331 | 0.1541 |
| HT29 | H1F0 | 1281 | 641 | 640 | 0 | 1 | 0 | 1 | 0.0000 |
| HT29 | HERC5 | 1040 | 520 | 520 | 2 | 0 | 0 | 2 | 0.0000 |
| HT29 | HERC6 | 826 | 413 | 413 | 35 | 23 | 1 | 57 | 0.0175 |
| HT29 | ID2 | 796 | 398 | 398 | 21 | 34 | 10 | 45 | 0.2222 |
| HT29 | IFNAR1 | 1148 | 574 | 574 | 637 | 463 | 287 | 813 | 0.3530 |
| HT29 | IFNAR2 | 1117 | 558 | 559 | 750 | 733 | 429 | 1054 | 0.4070 |
| HT29 | IRF1 | 1351 | 676 | 675 | 28 | 2 | 0 | 30 | 0.0000 |
| HT29 | IRF7 | 1086 | 543 | 543 | 14 | 2 | 0 | 16 | 0.0000 |
| HT29 | JAK1 | 1692 | 846 | 846 | 782 | 925 | 505 | 1202 | 0.4201 |
| HT29 | JUN | 787 | 393 | 394 | 39 | 69 | 9 | 99 | 0.0909 |
| HT29 | MAP3K14 | 980 | 490 | 490 | 3 | 24 | 1 | 26 | 0.0385 |
| HT29 | MAPK8 | 999 | 499 | 500 | 42 | 23 | 0 | 65 | 0.0000 |
| HT29 | PARP12 | 1062 | 531 | 531 | 88 | 46 | 2 | 132 | 0.0152 |
| HT29 | RAP1GAP | 1332 | 666 | 666 | 1 | 31 | 0 | 32 | 0.0000 |
| HT29 | RNF114 | 1439 | 720 | 719 | 80 | 6 | 2 | 84 | 0.0238 |
| HT29 | SOCS1 | 1188 | 594 | 594 | 25 | 0 | 0 | 25 | 0.0000 |
| HT29 | SP100 | 1534 | 767 | 767 | 12 | 23 | 1 | 34 | 0.0294 |
| HT29 | STAT1 | 1317 | 658 | 659 | 352 | 439 | 145 | 646 | 0.2245 |
| HT29 | STAT3 | 961 | 481 | 480 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| HT29 | STAT4 | 855 | 427 | 428 | 2 | 4 | 0 | 6 | 0.0000 |
| HT29 | STAT5A | 1278 | 639 | 639 | 3 | 0 | 0 | 3 | 0.0000 |
| HT29 | TRAFD1 | 1230 | 615 | 615 | 37 | 2 | 0 | 39 | 0.0000 |
| HT29 | TYK2 | 1108 | 554 | 554 | 399 | 513 | 249 | 663 | 0.3756 |
| HT29 | VAV1 | 1096 | 548 | 548 | 10 | 9 | 0 | 19 | 0.0000 |
| HT29 | ZBP1 | 1236 | 618 | 618 | 8 | 10 | 0 | 18 | 0.0000 |
| K562 | AKT1 | 481 | 241 | 240 | 0 | 3 | 0 | 3 | 0.0000 |
| K562 | BATF2 | 483 | 241 | 242 | 2 | 33 | 1 | 34 | 0.0294 |
| K562 | BRD9 | 692 | 346 | 346 | 1 | 26 | 0 | 27 | 0.0000 |
| K562 | CEBPG | 443 | 221 | 222 | 1 | 1 | 0 | 2 | 0.0000 |
| K562 | DTX3L | 525 | 263 | 262 | 1 | 2 | 1 | 2 | 0.5000 |
| K562 | ELK1 | 695 | 347 | 348 | 36 | 31 | 6 | 61 | 0.0984 |
| K562 | ETS1 | 485 | 243 | 242 | 4 | 3 | 2 | 5 | 0.4000 |
| K562 | ETS2 | 700 | 350 | 350 | 50 | 20 | 6 | 64 | 0.0938 |
| K562 | H1F0 | 629 | 315 | 314 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | HERC5 | 462 | 231 | 231 | 0 | 1 | 0 | 1 | 0.0000 |
| K562 | HERC6 | 408 | 204 | 204 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | ID2 | 531 | 265 | 266 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | IFNAR1 | 460 | 230 | 230 | 99 | 99 | 63 | 135 | 0.4667 |
| K562 | IFNAR2 | 439 | 219 | 220 | 236 | 103 | 81 | 258 | 0.3140 |
| K562 | IRF1 | 568 | 284 | 284 | 2 | 0 | 0 | 2 | 0.0000 |
| K562 | IRF7 | 528 | 264 | 264 | 4 | 0 | 0 | 4 | 0.0000 |
| K562 | JAK1 | 587 | 294 | 293 | 125 | 107 | 60 | 172 | 0.3488 |
| K562 | JUN | 503 | 251 | 252 | 1 | 4 | 0 | 5 | 0.0000 |
| K562 | MAP3K14 | 429 | 215 | 214 | 2 | 2 | 0 | 4 | 0.0000 |
| K562 | MAPK8 | 509 | 254 | 255 | 0 | 2 | 0 | 2 | 0.0000 |
| K562 | PARP12 | 647 | 324 | 323 | 8 | 2 | 0 | 10 | 0.0000 |
| K562 | RAP1GAP | 582 | 291 | 291 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | RNF114 | 602 | 301 | 301 | 5 | 1 | 0 | 6 | 0.0000 |
| K562 | SOCS1 | 583 | 291 | 292 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | SP100 | 579 | 290 | 289 | 4 | 1 | 1 | 4 | 0.2500 |
| K562 | STAT1 | 565 | 282 | 283 | 59 | 54 | 30 | 83 | 0.3614 |
| K562 | STAT3 | 478 | 239 | 239 | 2 | 0 | 0 | 2 | 0.0000 |
| K562 | STAT4 | 466 | 233 | 233 | 3 | 0 | 0 | 3 | 0.0000 |
| K562 | STAT5A | 869 | 435 | 434 | 76 | 26 | 12 | 90 | 0.1333 |
| K562 | TRAFD1 | 620 | 310 | 310 | 1 | 2 | 0 | 3 | 0.0000 |
| K562 | TYK2 | 532 | 266 | 266 | 116 | 105 | 55 | 166 | 0.3313 |
| K562 | VAV1 | 595 | 297 | 298 | 1 | 0 | 0 | 1 | 0.0000 |
| K562 | ZBP1 | 482 | 241 | 241 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| MCF7 | AKT1 | 682 | 341 | 341 | 414 | 33 | 8 | 439 | 0.0182 |
| MCF7 | BATF2 | 1203 | 601 | 602 | 166 | 181 | 57 | 290 | 0.1966 |
| MCF7 | BRD9 | 796 | 398 | 398 | 7 | 67 | 1 | 73 | 0.0137 |
| MCF7 | CEBPG | 827 | 413 | 414 | 182 | 757 | 75 | 864 | 0.0868 |
| MCF7 | DTX3L | 1786 | 893 | 893 | 204 | 117 | 56 | 265 | 0.2113 |
| MCF7 | ELK1 | 2212 | 1106 | 1106 | 406 | 791 | 259 | 938 | 0.2761 |
| MCF7 | ETS1 | 1045 | 523 | 522 | 90 | 81 | 12 | 159 | 0.0755 |
| MCF7 | ETS2 | 1368 | 684 | 684 | 221 | 287 | 105 | 403 | 0.2605 |
| MCF7 | H1F0 | 1128 | 564 | 564 | 458 | 123 | 32 | 549 | 0.0583 |
| MCF7 | HERC5 | 971 | 485 | 486 | 165 | 1093 | 101 | 1157 | 0.0873 |
| MCF7 | HERC6 | 1281 | 641 | 640 | 406 | 227 | 110 | 523 | 0.2103 |
| MCF7 | ID2 | 716 | 358 | 358 | 1078 | 475 | 154 | 1399 | 0.1101 |
| MCF7 | IFNAR1 | 644 | 322 | 322 | 593 | 710 | 347 | 956 | 0.3630 |
| MCF7 | IFNAR2 | 450 | 225 | 225 | 454 | 942 | 134 | 1262 | 0.1062 |
| MCF7 | IRF1 | 1803 | 902 | 901 | 863 | 1272 | 394 | 1741 | 0.2263 |
| MCF7 | IRF7 | 557 | 278 | 279 | 940 | 46 | 21 | 965 | 0.0218 |
| MCF7 | JAK1 | 1289 | 645 | 644 | 1045 | 215 | 76 | 1184 | 0.0642 |
| MCF7 | JUN | 1181 | 590 | 591 | 272 | 421 | 86 | 607 | 0.1417 |
| MCF7 | MAP3K14 | 1041 | 521 | 520 | 106 | 54 | 22 | 138 | 0.1594 |
| MCF7 | MAPK8 | 1119 | 559 | 560 | 533 | 986 | 257 | 1262 | 0.2036 |
| MCF7 | PARP12 | 1376 | 688 | 688 | 218 | 22 | 14 | 226 | 0.0619 |
| MCF7 | RAP1GAP | 1287 | 643 | 644 | 326 | 53 | 22 | 357 | 0.0616 |
| MCF7 | RNF114 | 984 | 492 | 492 | 122 | 325 | 34 | 413 | 0.0823 |
| MCF7 | SOCS1 | 863 | 431 | 432 | 587 | 144 | 41 | 690 | 0.0594 |
| MCF7 | SP100 | 724 | 362 | 362 | 83 | 470 | 36 | 517 | 0.0696 |
| MCF7 | STAT1 | 1277 | 638 | 639 | 456 | 224 | 79 | 601 | 0.1314 |
| MCF7 | STAT3 | 1141 | 571 | 570 | 132 | 97 | 18 | 211 | 0.0853 |
| MCF7 | STAT4 | 855 | 427 | 428 | 86 | 22 | 9 | 99 | 0.0909 |
| MCF7 | STAT5A | 1195 | 598 | 597 | 259 | 967 | 131 | 1095 | 0.1196 |
| MCF7 | TRAFD1 | 984 | 492 | 492 | 867 | 587 | 189 | 1265 | 0.1494 |
| MCF7 | TYK2 | 822 | 411 | 411 | 1868 | 1256 | 666 | 2458 | 0.2710 |
| MCF7 | VAV1 | 1137 | 568 | 569 | 1683 | 1684 | 645 | 2722 | 0.2370 |
| MCF7 | ZBP1 | 1598 | 799 | 799 | 720 | 1090 | 284 | 1526 | 0.1861 |

### Variante I -- detalle por (linea, perturbacion)

| linea | perturbacion | n celulas | n mitad A | n mitad B | n_DE_mitadA | n_DE_mitadB | interseccion | **\|union\|** | jac marcador |
|---|---|---|---|---|---|---|---|---|---|
| A549 | AKT1 | 606 | 303 | 303 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | BATF2 | 617 | 308 | 309 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | BRD9 | 699 | 350 | 349 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | CEBPG | 701 | 350 | 351 | 1 | 2 | 0 | 3 | 0.0000 |
| A549 | DTX3L | 544 | 272 | 272 | 1 | 1 | 1 | 1 | 1.0000 |
| A549 | ELK1 | 856 | 428 | 428 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | ETS1 | 622 | 311 | 311 | 0 | 2 | 0 | 2 | 0.0000 |
| A549 | ETS2 | 576 | 288 | 288 | 5 | 4 | 1 | 8 | 0.1250 |
| A549 | H1F0 | 671 | 336 | 335 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | HERC5 | 468 | 234 | 234 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | HERC6 | 440 | 220 | 220 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | ID2 | 576 | 288 | 288 | 0 | 1 | 0 | 1 | 0.0000 |
| A549 | IFNAR1 | 634 | 317 | 317 | 263 | 265 | 91 | 437 | 0.2082 |
| A549 | IFNAR2 | 723 | 361 | 362 | 411 | 494 | 193 | 712 | 0.2711 |
| A549 | IRF1 | 950 | 475 | 475 | 4 | 4 | 1 | 7 | 0.1429 |
| A549 | IRF7 | 436 | 218 | 218 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | JAK1 | 693 | 347 | 346 | 32 | 53 | 10 | 75 | 0.1333 |
| A549 | JUN | 525 | 262 | 263 | 4 | 6 | 2 | 8 | 0.2500 |
| A549 | MAP3K14 | 407 | 204 | 203 | 7 | 0 | 0 | 7 | 0.0000 |
| A549 | MAPK8 | 581 | 290 | 291 | 1 | 0 | 0 | 1 | 0.0000 |
| A549 | PARP12 | 682 | 341 | 341 | 1 | 1 | 0 | 2 | 0.0000 |
| A549 | RAP1GAP | 675 | 337 | 338 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | RNF114 | 688 | 344 | 344 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | SOCS1 | 597 | 298 | 299 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | SP100 | 496 | 248 | 248 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | STAT1 | 630 | 315 | 315 | 1 | 2 | 0 | 3 | 0.0000 |
| A549 | STAT3 | 505 | 253 | 252 | 13 | 33 | 6 | 40 | 0.1500 |
| A549 | STAT4 | 510 | 255 | 255 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | STAT5A | 842 | 421 | 421 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | TRAFD1 | 478 | 239 | 239 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | TYK2 | 737 | 369 | 368 | 523 | 664 | 255 | 932 | 0.2736 |
| A549 | VAV1 | 784 | 392 | 392 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| A549 | ZBP1 | 509 | 255 | 254 | 1 | 0 | 0 | 1 | 0.0000 |
| BXPC3 | AKT1 | 2499 | 1250 | 1249 | 32 | 16 | 8 | 40 | 0.2000 |
| BXPC3 | BATF2 | 712 | 356 | 356 | 37 | 35 | 10 | 62 | 0.1613 |
| BXPC3 | BRD9 | 1676 | 838 | 838 | 213 | 84 | 48 | 249 | 0.1928 |
| BXPC3 | CEBPG | 1814 | 907 | 907 | 112 | 89 | 22 | 179 | 0.1229 |
| BXPC3 | DTX3L | 2255 | 1128 | 1127 | 42 | 49 | 16 | 75 | 0.2133 |
| BXPC3 | ELK1 | 2174 | 1087 | 1087 | 23 | 47 | 12 | 58 | 0.2069 |
| BXPC3 | ETS1 | 984 | 492 | 492 | 369 | 452 | 217 | 604 | 0.3593 |
| BXPC3 | ETS2 | 2518 | 1259 | 1259 | 67 | 65 | 15 | 117 | 0.1282 |
| BXPC3 | H1F0 | 1170 | 585 | 585 | 9 | 21 | 3 | 27 | 0.1111 |
| BXPC3 | HERC5 | 1066 | 533 | 533 | 9 | 13 | 3 | 19 | 0.1579 |
| BXPC3 | HERC6 | 1640 | 820 | 820 | 169 | 134 | 56 | 247 | 0.2267 |
| BXPC3 | ID2 | 2031 | 1015 | 1016 | 298 | 311 | 92 | 517 | 0.1779 |
| BXPC3 | IFNAR1 | 1104 | 552 | 552 | 674 | 486 | 254 | 906 | 0.2804 |
| BXPC3 | IFNAR2 | 2901 | 1450 | 1451 | 2242 | 2500 | 1452 | 3290 | 0.4413 |
| BXPC3 | IRF1 | 1447 | 724 | 723 | 76 | 87 | 30 | 133 | 0.2256 |
| BXPC3 | IRF7 | 3860 | 1930 | 1930 | 67 | 54 | 23 | 98 | 0.2347 |
| BXPC3 | JAK1 | 1045 | 523 | 522 | 521 | 437 | 175 | 783 | 0.2235 |
| BXPC3 | JUN | 1709 | 854 | 855 | 82 | 94 | 33 | 143 | 0.2308 |
| BXPC3 | MAP3K14 | 1734 | 867 | 867 | 117 | 92 | 39 | 170 | 0.2294 |
| BXPC3 | MAPK8 | 1029 | 514 | 515 | 70 | 70 | 16 | 124 | 0.1290 |
| BXPC3 | PARP12 | 1178 | 589 | 589 | 16 | 9 | 3 | 22 | 0.1364 |
| BXPC3 | RAP1GAP | 2832 | 1416 | 1416 | 17 | 28 | 9 | 36 | 0.2500 |
| BXPC3 | RNF114 | 1007 | 504 | 503 | 29 | 17 | 4 | 42 | 0.0952 |
| BXPC3 | SOCS1 | 1907 | 953 | 954 | 133 | 90 | 44 | 179 | 0.2458 |
| BXPC3 | SP100 | 2057 | 1029 | 1028 | 194 | 150 | 55 | 289 | 0.1903 |
| BXPC3 | STAT1 | 1702 | 851 | 851 | 376 | 467 | 188 | 655 | 0.2870 |
| BXPC3 | STAT3 | 1354 | 677 | 677 | 20 | 24 | 6 | 38 | 0.1579 |
| BXPC3 | STAT4 | 1520 | 760 | 760 | 15 | 39 | 5 | 49 | 0.1020 |
| BXPC3 | STAT5A | 3004 | 1502 | 1502 | 7 | 3 | 2 | 8 | 0.2500 |
| BXPC3 | TRAFD1 | 1231 | 615 | 616 | 17 | 4 | 1 | 20 | 0.0500 |
| BXPC3 | TYK2 | 2687 | 1344 | 1343 | 2317 | 2258 | 1428 | 3147 | 0.4538 |
| BXPC3 | VAV1 | 2172 | 1086 | 1086 | 19 | 19 | 3 | 35 | 0.0857 |
| BXPC3 | ZBP1 | 1297 | 649 | 648 | 9 | 14 | 7 | 16 | 0.4375 |
| HAP1 | AKT1 | 1014 | 507 | 507 | 1 | 0 | 0 | 1 | 0.0000 |
| HAP1 | BATF2 | 1326 | 663 | 663 | 10 | 12 | 5 | 17 | 0.2941 |
| HAP1 | BRD9 | 975 | 488 | 487 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| HAP1 | CEBPG | 1485 | 742 | 743 | 39 | 83 | 13 | 109 | 0.1193 |
| HAP1 | DTX3L | 1321 | 661 | 660 | 2 | 5 | 2 | 5 | 0.4000 |
| HAP1 | ELK1 | 641 | 320 | 321 | 7 | 5 | 2 | 10 | 0.2000 |
| HAP1 | ETS1 | 1044 | 522 | 522 | 8 | 9 | 3 | 14 | 0.2143 |
| HAP1 | ETS2 | 1587 | 793 | 794 | 1 | 1 | 0 | 2 | 0.0000 |
| HAP1 | H1F0 | 1604 | 802 | 802 | 2 | 0 | 0 | 2 | 0.0000 |
| HAP1 | HERC5 | 995 | 497 | 498 | 0 | 1 | 0 | 1 | 0.0000 |
| HAP1 | HERC6 | 1169 | 585 | 584 | 7 | 5 | 4 | 8 | 0.5000 |
| HAP1 | ID2 | 1394 | 697 | 697 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| HAP1 | IFNAR1 | 1152 | 576 | 576 | 207 | 226 | 141 | 292 | 0.4829 |
| HAP1 | IFNAR2 | 1261 | 630 | 631 | 241 | 245 | 151 | 335 | 0.4507 |
| HAP1 | IRF1 | 1698 | 849 | 849 | 11 | 12 | 5 | 18 | 0.2778 |
| HAP1 | IRF7 | 1328 | 664 | 664 | 2 | 6 | 2 | 6 | 0.3333 |
| HAP1 | JAK1 | 1471 | 736 | 735 | 85 | 100 | 65 | 120 | 0.5417 |
| HAP1 | JUN | 1604 | 802 | 802 | 0 | 1 | 0 | 1 | 0.0000 |
| HAP1 | MAP3K14 | 982 | 491 | 491 | 1 | 0 | 0 | 1 | 0.0000 |
| HAP1 | MAPK8 | 1385 | 692 | 693 | 2 | 1 | 0 | 3 | 0.0000 |
| HAP1 | PARP12 | 1481 | 741 | 740 | 5 | 5 | 2 | 8 | 0.2500 |
| HAP1 | RAP1GAP | 1733 | 866 | 867 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| HAP1 | RNF114 | 1738 | 869 | 869 | 3 | 4 | 3 | 4 | 0.7500 |
| HAP1 | SOCS1 | 1526 | 763 | 763 | 16 | 16 | 5 | 27 | 0.1852 |
| HAP1 | SP100 | 1234 | 617 | 617 | 0 | 2 | 0 | 2 | 0.0000 |
| HAP1 | STAT1 | 2025 | 1012 | 1013 | 321 | 390 | 187 | 524 | 0.3569 |
| HAP1 | STAT3 | 1165 | 583 | 582 | 4 | 0 | 0 | 4 | 0.0000 |
| HAP1 | STAT4 | 1369 | 684 | 685 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| HAP1 | STAT5A | 1793 | 897 | 896 | 7 | 4 | 4 | 7 | 0.5714 |
| HAP1 | TRAFD1 | 796 | 398 | 398 | 2 | 0 | 0 | 2 | 0.0000 |
| HAP1 | TYK2 | 736 | 368 | 368 | 108 | 118 | 83 | 143 | 0.5804 |
| HAP1 | VAV1 | 2049 | 1024 | 1025 | 2 | 0 | 0 | 2 | 0.0000 |
| HAP1 | ZBP1 | 1218 | 609 | 609 | 0 | 1 | 0 | 1 | 0.0000 |
| HT29 | AKT1 | 1165 | 583 | 582 | 1 | 4 | 0 | 5 | 0.0000 |
| HT29 | BATF2 | 1002 | 501 | 501 | 1 | 0 | 0 | 1 | 0.0000 |
| HT29 | BRD9 | 996 | 498 | 498 | 0 | 3 | 0 | 3 | 0.0000 |
| HT29 | CEBPG | 1118 | 559 | 559 | 15 | 3 | 1 | 17 | 0.0588 |
| HT29 | DTX3L | 879 | 440 | 439 | 7 | 10 | 1 | 16 | 0.0625 |
| HT29 | ELK1 | 1015 | 507 | 508 | 3 | 0 | 0 | 3 | 0.0000 |
| HT29 | ETS1 | 1083 | 542 | 541 | 2 | 2 | 0 | 4 | 0.0000 |
| HT29 | ETS2 | 968 | 484 | 484 | 130 | 116 | 32 | 214 | 0.1495 |
| HT29 | H1F0 | 1281 | 641 | 640 | 0 | 1 | 0 | 1 | 0.0000 |
| HT29 | HERC5 | 1040 | 520 | 520 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| HT29 | HERC6 | 826 | 413 | 413 | 0 | 8 | 0 | 8 | 0.0000 |
| HT29 | ID2 | 796 | 398 | 398 | 13 | 16 | 6 | 23 | 0.2609 |
| HT29 | IFNAR1 | 1148 | 574 | 574 | 475 | 392 | 226 | 641 | 0.3526 |
| HT29 | IFNAR2 | 1117 | 558 | 559 | 589 | 629 | 320 | 898 | 0.3563 |
| HT29 | IRF1 | 1351 | 676 | 675 | 10 | 3 | 0 | 13 | 0.0000 |
| HT29 | IRF7 | 1086 | 543 | 543 | 15 | 0 | 0 | 15 | 0.0000 |
| HT29 | JAK1 | 1692 | 846 | 846 | 602 | 704 | 368 | 938 | 0.3923 |
| HT29 | JUN | 787 | 393 | 394 | 21 | 33 | 4 | 50 | 0.0800 |
| HT29 | MAP3K14 | 980 | 490 | 490 | 1 | 0 | 0 | 1 | 0.0000 |
| HT29 | MAPK8 | 999 | 499 | 500 | 24 | 5 | 0 | 29 | 0.0000 |
| HT29 | PARP12 | 1062 | 531 | 531 | 33 | 8 | 0 | 41 | 0.0000 |
| HT29 | RAP1GAP | 1332 | 666 | 666 | 0 | 1 | 0 | 1 | 0.0000 |
| HT29 | RNF114 | 1439 | 720 | 719 | 47 | 2 | 1 | 48 | 0.0208 |
| HT29 | SOCS1 | 1188 | 594 | 594 | 37 | 2 | 0 | 39 | 0.0000 |
| HT29 | SP100 | 1534 | 767 | 767 | 6 | 13 | 1 | 18 | 0.0556 |
| HT29 | STAT1 | 1317 | 658 | 659 | 253 | 269 | 94 | 428 | 0.2196 |
| HT29 | STAT3 | 961 | 481 | 480 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| HT29 | STAT4 | 855 | 427 | 428 | 0 | 1 | 0 | 1 | 0.0000 |
| HT29 | STAT5A | 1278 | 639 | 639 | 1 | 1 | 0 | 2 | 0.0000 |
| HT29 | TRAFD1 | 1230 | 615 | 615 | 5 | 1 | 0 | 6 | 0.0000 |
| HT29 | TYK2 | 1108 | 554 | 554 | 313 | 430 | 207 | 536 | 0.3862 |
| HT29 | VAV1 | 1096 | 548 | 548 | 2 | 2 | 0 | 4 | 0.0000 |
| HT29 | ZBP1 | 1236 | 618 | 618 | 5 | 16 | 0 | 21 | 0.0000 |
| K562 | AKT1 | 481 | 241 | 240 | 0 | 1 | 0 | 1 | 0.0000 |
| K562 | BATF2 | 483 | 241 | 242 | 4 | 8 | 1 | 11 | 0.0909 |
| K562 | BRD9 | 692 | 346 | 346 | 2 | 3 | 0 | 5 | 0.0000 |
| K562 | CEBPG | 443 | 221 | 222 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | DTX3L | 525 | 263 | 262 | 1 | 1 | 1 | 1 | 1.0000 |
| K562 | ELK1 | 695 | 347 | 348 | 8 | 9 | 1 | 16 | 0.0625 |
| K562 | ETS1 | 485 | 243 | 242 | 0 | 1 | 0 | 1 | 0.0000 |
| K562 | ETS2 | 700 | 350 | 350 | 32 | 7 | 5 | 34 | 0.1471 |
| K562 | H1F0 | 629 | 315 | 314 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | HERC5 | 462 | 231 | 231 | 3 | 0 | 0 | 3 | 0.0000 |
| K562 | HERC6 | 408 | 204 | 204 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | ID2 | 531 | 265 | 266 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | IFNAR1 | 460 | 230 | 230 | 89 | 89 | 58 | 120 | 0.4833 |
| K562 | IFNAR2 | 439 | 219 | 220 | 187 | 93 | 64 | 216 | 0.2963 |
| K562 | IRF1 | 568 | 284 | 284 | 1 | 1 | 0 | 2 | 0.0000 |
| K562 | IRF7 | 528 | 264 | 264 | 6 | 0 | 0 | 6 | 0.0000 |
| K562 | JAK1 | 587 | 294 | 293 | 94 | 89 | 51 | 132 | 0.3864 |
| K562 | JUN | 503 | 251 | 252 | 2 | 2 | 0 | 4 | 0.0000 |
| K562 | MAP3K14 | 429 | 215 | 214 | 2 | 0 | 0 | 2 | 0.0000 |
| K562 | MAPK8 | 509 | 254 | 255 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | PARP12 | 647 | 324 | 323 | 8 | 0 | 0 | 8 | 0.0000 |
| K562 | RAP1GAP | 582 | 291 | 291 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | RNF114 | 602 | 301 | 301 | 3 | 1 | 0 | 4 | 0.0000 |
| K562 | SOCS1 | 583 | 291 | 292 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | SP100 | 579 | 290 | 289 | 1 | 1 | 1 | 1 | 1.0000 |
| K562 | STAT1 | 565 | 282 | 283 | 58 | 42 | 26 | 74 | 0.3514 |
| K562 | STAT3 | 478 | 239 | 239 | 4 | 0 | 0 | 4 | 0.0000 |
| K562 | STAT4 | 466 | 233 | 233 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| K562 | STAT5A | 869 | 435 | 434 | 37 | 9 | 6 | 40 | 0.1500 |
| K562 | TRAFD1 | 620 | 310 | 310 | 0 | 1 | 0 | 1 | 0.0000 |
| K562 | TYK2 | 532 | 266 | 266 | 102 | 87 | 48 | 141 | 0.3404 |
| K562 | VAV1 | 595 | 297 | 298 | 0 | 1 | 0 | 1 | 0.0000 |
| K562 | ZBP1 | 482 | 241 | 241 | 0 | 0 | 0 | 0 | 1.0000 **(union vacia)** |
| MCF7 | AKT1 | 682 | 341 | 341 | 605 | 12 | 2 | 615 | 0.0033 |
| MCF7 | BATF2 | 1203 | 601 | 602 | 181 | 102 | 41 | 242 | 0.1694 |
| MCF7 | BRD9 | 796 | 398 | 398 | 6 | 4 | 0 | 10 | 0.0000 |
| MCF7 | CEBPG | 827 | 413 | 414 | 186 | 224 | 34 | 376 | 0.0904 |
| MCF7 | DTX3L | 1786 | 893 | 893 | 192 | 54 | 19 | 227 | 0.0837 |
| MCF7 | ELK1 | 2212 | 1106 | 1106 | 353 | 356 | 145 | 564 | 0.2571 |
| MCF7 | ETS1 | 1045 | 523 | 522 | 102 | 27 | 6 | 123 | 0.0488 |
| MCF7 | ETS2 | 1368 | 684 | 684 | 208 | 163 | 75 | 296 | 0.2534 |
| MCF7 | H1F0 | 1128 | 564 | 564 | 495 | 9 | 2 | 502 | 0.0040 |
| MCF7 | HERC5 | 971 | 485 | 486 | 190 | 358 | 65 | 483 | 0.1346 |
| MCF7 | HERC6 | 1281 | 641 | 640 | 419 | 134 | 62 | 491 | 0.1263 |
| MCF7 | ID2 | 716 | 358 | 358 | 1214 | 153 | 66 | 1301 | 0.0507 |
| MCF7 | IFNAR1 | 644 | 322 | 322 | 567 | 547 | 286 | 828 | 0.3454 |
| MCF7 | IFNAR2 | 450 | 225 | 225 | 522 | 518 | 98 | 942 | 0.1040 |
| MCF7 | IRF1 | 1803 | 902 | 901 | 866 | 416 | 147 | 1135 | 0.1295 |
| MCF7 | IRF7 | 557 | 278 | 279 | 1121 | 11 | 9 | 1123 | 0.0080 |
| MCF7 | JAK1 | 1289 | 645 | 644 | 1173 | 36 | 14 | 1195 | 0.0117 |
| MCF7 | JUN | 1181 | 590 | 591 | 315 | 130 | 41 | 404 | 0.1015 |
| MCF7 | MAP3K14 | 1041 | 521 | 520 | 109 | 46 | 16 | 139 | 0.1151 |
| MCF7 | MAPK8 | 1119 | 559 | 560 | 601 | 426 | 148 | 879 | 0.1684 |
| MCF7 | PARP12 | 1376 | 688 | 688 | 229 | 14 | 8 | 235 | 0.0340 |
| MCF7 | RAP1GAP | 1287 | 643 | 644 | 387 | 5 | 2 | 390 | 0.0051 |
| MCF7 | RNF114 | 984 | 492 | 492 | 146 | 49 | 10 | 185 | 0.0541 |
| MCF7 | SOCS1 | 863 | 431 | 432 | 727 | 26 | 10 | 743 | 0.0135 |
| MCF7 | SP100 | 724 | 362 | 362 | 111 | 155 | 25 | 241 | 0.1037 |
| MCF7 | STAT1 | 1277 | 638 | 639 | 544 | 68 | 34 | 578 | 0.0588 |
| MCF7 | STAT3 | 1141 | 571 | 570 | 203 | 30 | 5 | 228 | 0.0219 |
| MCF7 | STAT4 | 855 | 427 | 428 | 124 | 11 | 7 | 128 | 0.0547 |
| MCF7 | STAT5A | 1195 | 598 | 597 | 266 | 289 | 58 | 497 | 0.1167 |
| MCF7 | TRAFD1 | 984 | 492 | 492 | 955 | 238 | 91 | 1102 | 0.0826 |
| MCF7 | TYK2 | 822 | 411 | 411 | 1942 | 803 | 439 | 2306 | 0.1904 |
| MCF7 | VAV1 | 1137 | 568 | 569 | 1824 | 608 | 262 | 2170 | 0.1207 |
| MCF7 | ZBP1 | 1598 | 799 | 799 | 871 | 309 | 101 | 1079 | 0.0936 |

## 20 filas trabajadas, con los intermedios completos

Veinte filas (linea, perturbacion) elegidas para cubrir los dos extremos de n, las dos variantes y los casos de union vacia. Cada fila trae TODO lo que entra en sus numeros, de la celula cruda al jaccard: no hay ningun paso intermedio oculto.

| # | var | linea | pert | n celulas | mitad A | mitad B | n_DE_A | n_DE_B | inter | \|union\| | jac recomputado | jac marcador | dif | pds_cosine | reach | fid |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | S | A549 | MAP3K14 | 407 | 204 | 203 | 5 | 5 | 2 | 8 | 0.250000 | 0.250000 | 0.00e+00 | 1.000000 | 0.800000 | 1.000000 |
| 2 | I | A549 | MAP3K14 | 407 | 204 | 203 | 7 | 0 | 0 | 7 | 0.000000 | 0.000000 | 0.00e+00 | 1.000000 | NaN | 1.000000 |
| 3 | S | K562 | HERC6 | 408 | 204 | 204 | 0 | 0 | 0 | 0 | 1.000000 | 1.000000 | 0.00e+00 | 0.687500 | NaN | NaN |
| 4 | I | K562 | HERC6 | 408 | 204 | 204 | 0 | 0 | 0 | 0 | 1.000000 | 1.000000 | 0.00e+00 | 0.968750 | NaN | NaN |
| 5 | S | BXPC3 | STAT5A | 3004 | 1502 | 1502 | 13 | 10 | 4 | 19 | 0.210526 | 0.210526 | 0.00e+00 | 1.000000 | 1.000000 | 0.923077 |
| 6 | I | BXPC3 | STAT5A | 3004 | 1502 | 1502 | 7 | 3 | 2 | 8 | 0.250000 | 0.250000 | 0.00e+00 | 1.000000 | 1.000000 | 0.857143 |
| 7 | S | BXPC3 | IRF7 | 3860 | 1930 | 1930 | 153 | 87 | 51 | 189 | 0.269841 | 0.269841 | 0.00e+00 | 1.000000 | 1.000000 | 0.986928 |
| 8 | I | BXPC3 | IRF7 | 3860 | 1930 | 1930 | 67 | 54 | 23 | 98 | 0.234694 | 0.234694 | 0.00e+00 | 1.000000 | 0.981481 | 0.939394 |
| 9 | S | A549 | AKT1 | 606 | 303 | 303 | 0 | 0 | 0 | 0 | 1.000000 | 1.000000 | 0.00e+00 | 0.968750 | NaN | NaN |
| 10 | S | A549 | BRD9 | 699 | 350 | 349 | 0 | 0 | 0 | 0 | 1.000000 | 1.000000 | 0.00e+00 | 0.906250 | NaN | NaN |
| 11 | S | A549 | ELK1 | 856 | 428 | 428 | 0 | 0 | 0 | 0 | 1.000000 | 1.000000 | 0.00e+00 | 1.000000 | NaN | NaN |
| 12 | S | A549 | H1F0 | 671 | 336 | 335 | 0 | 0 | 0 | 0 | 1.000000 | 1.000000 | 0.00e+00 | 0.656250 | NaN | NaN |
| 13 | S | BXPC3 | IFNAR2 | 2901 | 1450 | 1451 | 2992 | 3050 | 2159 | 3883 | 0.556013 | 0.556013 | 0.00e+00 | 1.000000 | 1.000000 | 0.967213 |
| 14 | S | BXPC3 | TYK2 | 2687 | 1344 | 1343 | 2945 | 2630 | 1920 | 3655 | 0.525308 | 0.525308 | 0.00e+00 | 1.000000 | 1.000000 | 0.980645 |
| 15 | I | BXPC3 | IFNAR2 | 2901 | 1450 | 1451 | 2242 | 2500 | 1452 | 3290 | 0.441337 | 0.441337 | 0.00e+00 | 1.000000 | 0.992000 | 0.856400 |
| 16 | I | BXPC3 | TYK2 | 2687 | 1344 | 1343 | 2317 | 2258 | 1428 | 3147 | 0.453765 | 0.453765 | 0.00e+00 | 1.000000 | 0.989371 | 0.964023 |
| 17 | S | A549 | BATF2 | 617 | 308 | 309 | 0 | 1 | 0 | 1 | 0.000000 | 0.000000 | 0.00e+00 | 0.937500 | 1.000000 | 0.000000 |
| 18 | S | A549 | CEBPG | 701 | 350 | 351 | 1 | 2 | 0 | 3 | 0.000000 | 0.000000 | 0.00e+00 | 1.000000 | 1.000000 | 0.500000 |
| 19 | S | A549 | DTX3L | 544 | 272 | 272 | 2 | 1 | 1 | 2 | 0.500000 | 0.500000 | 0.00e+00 | 0.500000 | 1.000000 | 1.000000 |
| 20 | S | A549 | ETS1 | 622 | 311 | 311 | 0 | 1 | 0 | 1 | 0.000000 | 0.000000 | 0.00e+00 | 0.968750 | 1.000000 | 0.000000 |

La columna `dif` es |jac recomputado - jac marcador|: cierra el circuito entre las cuentas de DE de esta tabla y el numero que el marcador puntuo. Si fuera distinto de ~0, las columnas n_DE/|union| no describirian lo que se puntuo.

## Verificacion: las filas por perturbacion componen el agregado

`agg_results.csv` reporta la fila `mean` sobre las 33 perturbaciones. Aqui se recalcula esa media desde `results.csv` (por perturbacion) y se compara con la que se reporto en las tablas de arriba, para cada corrida. `jac cruce` es el maximo |jac recomputado - jac marcador| sobre las 33 perturbaciones de esa corrida.

| linea | var | dir | jac agregado | media de las 33 | dif | pds agregado | media de las 33 | dif | jac cruce max |
|---|---|---|---|---|---|---|---|---|---|
| A549 | S | A | 0.403355 | 0.403355 | 5.55e-17 | 0.865530 | 0.865530 | 0.00e+00 | 0.00e+00 |
| BXPC3 | S | A | 0.247759 | 0.247759 | 0.00e+00 | 0.996212 | 0.996212 | 0.00e+00 | 0.00e+00 |
| HAP1 | S | A | 0.338312 | 0.338312 | 5.55e-17 | 0.892045 | 0.892045 | 0.00e+00 | 0.00e+00 |
| HT29 | S | A | 0.105407 | 0.105407 | 1.39e-17 | 0.914773 | 0.914773 | 0.00e+00 | 0.00e+00 |
| K562 | S | A | 0.282639 | 0.282639 | 0.00e+00 | 0.855114 | 0.855114 | 0.00e+00 | 0.00e+00 |
| MCF7 | S | A | 0.136259 | 0.136259 | 5.55e-17 | 0.987689 | 0.987689 | 0.00e+00 | 0.00e+00 |
| A549 | S | B | 0.403355 | 0.403355 | 5.55e-17 | 0.835227 | 0.835227 | 0.00e+00 | 0.00e+00 |
| BXPC3 | S | B | 0.247759 | 0.247759 | 0.00e+00 | 0.997159 | 0.997159 | 0.00e+00 | 0.00e+00 |
| HAP1 | S | B | 0.338312 | 0.338312 | 5.55e-17 | 0.903409 | 0.903409 | 0.00e+00 | 0.00e+00 |
| HT29 | S | B | 0.105407 | 0.105407 | 1.39e-17 | 0.920455 | 0.920455 | 0.00e+00 | 0.00e+00 |
| K562 | S | B | 0.282639 | 0.282639 | 0.00e+00 | 0.831439 | 0.831439 | 0.00e+00 | 0.00e+00 |
| MCF7 | S | B | 0.136259 | 0.136259 | 5.55e-17 | 0.984848 | 0.984848 | 0.00e+00 | 0.00e+00 |
| A549 | I | A | 0.562245 | 0.562245 | 0.00e+00 | 0.841856 | 0.841856 | 0.00e+00 | 0.00e+00 |
| BXPC3 | I | A | 0.211960 | 0.211960 | 5.55e-17 | 0.997159 | 0.997159 | 0.00e+00 | 0.00e+00 |
| HAP1 | I | A | 0.318423 | 0.318423 | 5.55e-17 | 0.928977 | 0.928977 | 0.00e+00 | 0.00e+00 |
| HT29 | I | A | 0.133187 | 0.133187 | 2.78e-17 | 0.908144 | 0.908144 | 0.00e+00 | 0.00e+00 |
| K562 | I | A | 0.403280 | 0.403280 | 5.55e-17 | 0.856061 | 0.856061 | 0.00e+00 | 0.00e+00 |
| MCF7 | I | A | 0.095610 | 0.095610 | 1.39e-17 | 0.988636 | 0.988636 | 0.00e+00 | 0.00e+00 |
| A549 | I | B | 0.562245 | 0.562245 | 0.00e+00 | 0.782197 | 0.782197 | 0.00e+00 | 0.00e+00 |
| BXPC3 | I | B | 0.211960 | 0.211960 | 5.55e-17 | 0.991477 | 0.991477 | 0.00e+00 | 0.00e+00 |
| HAP1 | I | B | 0.318423 | 0.318423 | 5.55e-17 | 0.910985 | 0.910985 | 0.00e+00 | 0.00e+00 |
| HT29 | I | B | 0.133187 | 0.133187 | 2.78e-17 | 0.840909 | 0.840909 | 0.00e+00 | 0.00e+00 |
| K562 | I | B | 0.403280 | 0.403280 | 5.55e-17 | 0.809659 | 0.809659 | 0.00e+00 | 0.00e+00 |
| MCF7 | I | B | 0.095610 | 0.095610 | 1.39e-17 | 0.985795 | 0.985795 | 0.00e+00 | 0.00e+00 |

## Intermedios de corrida (tiempo, disco, formas, semillas)

| linea | var | dir | control_source | real_shape | pred_shape | run_s | total_s | free_gib | seed_split | seed_ctrl_boot | score |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A549 | S | A | real | [11608, 34025] | [10629, 34025] | 55.59 | 61.33 | 36.13 | 20262000 | 20263000 | ok |
| BXPC3 | S | A | real | [32885, 34025] | [30059, 34025] | 128.32 | 134.29 | 36.13 | 20262001 | 20263001 | ok |
| HAP1 | S | A | real | [25323, 34025] | [22548, 34025] | 111.32 | 118.25 | 36.11 | 20262002 | 20263002 | nego (ver nota) |
| HT29 | S | A | real | [21218, 34025] | [18883, 34025] | 89.38 | 94.26 | 36.13 | 20262003 | 20263003 | ok |
| K562 | S | A | real | [10477, 34025] | [9478, 34025] | 53.84 | 57.94 | 36.13 | 20262004 | 20263004 | ok |
| MCF7 | S | A | real | [20842, 34025] | [18569, 34025] | 98.19 | 103.51 | 36.13 | 20262005 | 20263005 | ok |
| A549 | S | B | real | [11608, 34025] | [10629, 34025] | 54.01 | 58.2 | 36.13 | 20262000 | 20263000 | ok |
| BXPC3 | S | B | real | [32887, 34025] | [30057, 34025] | 127.91 | 133.87 | 36.13 | 20262001 | 20263001 | ok |
| HAP1 | S | B | real | [25320, 34025] | [22551, 34025] | 111.11 | 116.64 | 36.13 | 20262002 | 20263002 | ok |
| HT29 | S | B | real | [21219, 34025] | [18882, 34025] | 89.35 | 94.25 | 36.13 | 20262003 | 20263003 | ok |
| K562 | S | B | real | [10476, 34025] | [9479, 34025] | 53.64 | 57.75 | 36.13 | 20262004 | 20263004 | ok |
| MCF7 | S | B | real | [20838, 34025] | [18573, 34025] | 98.01 | 103.36 | 36.13 | 20262005 | 20263005 | ok |
| A549 | I | A | pred | [10919, 34025] | [10918, 34025] | 50.31 | 54.38 | 36.13 | 20262000 | None | ok |
| BXPC3 | I | A | pred | [31271, 34025] | [31273, 34025] | 115.44 | 121.45 | 36.13 | 20262001 | None | ok |
| HAP1 | I | A | pred | [23737, 34025] | [23734, 34025] | 97.94 | 103.39 | 36.13 | 20262002 | None | nego (ver nota) |
| HT29 | I | A | pred | [19850, 34025] | [19851, 34025] | 80.06 | 84.94 | 36.13 | 20262003 | None | ok |
| K562 | I | A | pred | [9778, 34025] | [9777, 34025] | 49.17 | 53.26 | 36.13 | 20262004 | None | ok |
| MCF7 | I | A | pred | [19508, 34025] | [19503, 34025] | 86.21 | 91.53 | 36.13 | 20262005 | None | ok |
| A549 | I | B | pred | [10918, 34025] | [10919, 34025] | 50.07 | 54.16 | 36.13 | 20262000 | None | ok |
| BXPC3 | I | B | pred | [31273, 34025] | [31271, 34025] | 115.43 | 121.43 | 36.13 | 20262001 | None | ok |
| HAP1 | I | B | pred | [23734, 34025] | [23737, 34025] | 97.7 | 103.18 | 36.13 | 20262002 | None | ok |
| HT29 | I | B | pred | [19851, 34025] | [19850, 34025] | 79.96 | 84.82 | 36.13 | 20262003 | None | ok |
| K562 | I | B | pred | [9777, 34025] | [9778, 34025] | 49.17 | 53.25 | 36.13 | 20262004 | None | ok |
| MCF7 | I | B | pred | [19503, 34025] | [19508, 34025] | 86.25 | 91.59 | 36.13 | 20262005 | None | ok |

## Semillas

| uso | familia | valor |
|---|---|---|
| VERDAD del suelo | 2026000x | 20260001..20260006 (no reusada) |
| PREDICCION del suelo | 2026100x | 20261001..20261006 (no reusada) |
| particion A/B del techo | 2026200x | 20262000..20262005 |
| bootstrap del control de PRED, variante S | 2026300x | 20263000..20263005 |

Nunca 42.

## Lo que NO se hizo, a proposito

- **No se promedio ninguna metrica a mano.** Las seis columnas son la fila `mean` de `agg_results.csv` que agrega el marcador, igual que en el suelo. La seccion de verificacion de arriba recalcula esa media solo para comprobar que coincide, no para sustituirla.
- **No se promediaron las dos direcciones** en las tablas del entregable: se reportan las dos. El unico promedio de direcciones esta en la tabla del intervalo S..I, y solo para dar su ancho.
- **No se eligio entre S e I.** El techo es el intervalo.
- **No se corrigio el techo por la perdida de potencia** de usar media poblacion (ni con Spearman-Brown ni de otro modo). Queda declarado como cota inferior.
- **No se explica por que `reach` da no-nulo con conjunto predicho vacio.** Sigue declarado como no entendido, pendiente de medirse aparte.
- No se envio nada ni se corrio `vcc sample`. No se hizo ningun commit.
- No se leyo NADA bajo `~/vcc/benchmark/` ni los `.h5ad` de `/datos/benchmark_VEDADO_techo_metricas`. Ver la nota al final.

## `avg_score`: 0.0 cuando corre, y donde `score` se niega a correr

`score_vs_self` (reusado sin modificar) pasa el mismo `agg_results.csv` como `--user-agg` y `--baseline-agg`, asi que `from_baseline` cae en 0.0 por construccion. La columna se conserva solo por paridad de forma con `bucle_suelo.md`; lo informativo son las seis columnas de metricas crudas.

En **2 de 24 corridas** `score` se nego a puntuar, y el motivo es un resultado del techo, no un fallo del arnes:

| linea | var | dir | mse_capped_norm del techo | motivo |
|---|---|---|---|---|
| HAP1 | I | A | -0.000803 | cell-eval2 score fallo (exit=1). Ver /datos/modelo/out/techo/run_HAP1_I_A/cell_eval2_score.log |
| HAP1 | S | A | -0.000847 | cell-eval2 score fallo (exit=1). Ver /datos/modelo/out/techo/run_HAP1_S_A/cell_eval2_score.log |

`expr_mse_unbiased_capped_norm` es un estimador **INSESGADO** y **con signo**. Cuando las dos mitades concuerdan -- que es exactamente lo que mide el techo -- el MSE verdadero tiende a 0, y un estimador insesgado de algo que vale ~0 sale negativo la mitad de las veces. `score` entonces no puede definir una escala (`denominator ... is not a finite positive number`) y se niega, incluso con `--allow-config-mismatch`. Consecuencia practica, que conviene no descubrir mas tarde: **el techo no sirve como `--baseline-agg` para esa metrica**; sirve como el `r` de la escala `(u - b)/(r - b)` calculada a mano, con `r` cerca de 0. Es tambien la razon por la que `cell_eval2/ceiling.py` excluye esa metrica de su lista `SB_METRICS`, citando que "is additionally signed, so it is not a bounded reliability in the first place".

## La escala que esto habilita

Con los dos anclajes medidos en la misma definicion del marcador, un candidato `u` se lee por metrica como `(u - b) / (r - b)`, con `b` el suelo de `bucle_suelo.md` y `r` el techo de aqui. Dos advertencias que vienen de este reporte y no deben perderse al usar esa escala: `r` esta medido a la mitad de la profundidad de `b` (cota inferior, sin corregir), y solo la variante S comparte `control_source` con `b`.

## Nota sobre el limite duro de `~/vcc/benchmark/`

No se abrio ningun archivo bajo `~/vcc/benchmark/`, ni `techo_metricas.md`, ni `techo_estado.md`, ni `b_techo_i6_ruido.md`, ni nada en `benchmark/work/techo_metricas_work/`, ni los `.h5ad` de `/datos/benchmark_VEDADO_techo_metricas`. Lo unico que esta sesion vio de ahi fue el **nombre** del directorio `benchmark/` en un `ls` del raiz del proyecto, al ubicar `modelo/`; no se listo su interior ni se leyo nada. Ese trabajo responde a otra pregunta (como responde el marcador a un envio transformado) sobre otra poblacion (170,846 celulas del split test 2025) y no es un techo, asi que no habria servido aqui ni para comparar.

Las unicas fuentes de datos de este reporte son: `/datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad` (IFNB), `modelo/out/suelo_panel.json` (el panel comun de 33), `modelo/out/bucle_suelo.json` y los `results.csv` de `/datos/modelo/out/run_*` (el suelo), y el codigo de `cell_eval2` 0.20.0 en `~/cell-eval2/src/cell_eval2`.

