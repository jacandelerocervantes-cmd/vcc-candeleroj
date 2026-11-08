# PASO 5a — La seleccion: lo que se computa SOLO del control

## 🔴 CÓMO SE LEE

* **La Parte 1 es sobre la superficie de ENTREGA (D/E/F) y la de validacion (A/B/C).** Sus cifras son de ellas, computadas de sus 18,400 controles. El filtro del preset lo decide la referencia; en la puntuacion real la referencia es ese pool, asi que el filtro de cada contexto que se da aqui no es una aproximacion: es la misma cuenta (comprobada gen a gen contra `compute_de`, §1.1).
* **La Parte 2 es sobre OTRAS superficies** (Jiang, 6 lineas; VCC 2025, un contexto), las que tienen verdad. De ella solo viaja el **MECANISMO** — que el filtro del contexto excluye genes de R(g), que la potencia selecciona — **nunca las fracciones**.
* Las dos partes **no se mezclan en ninguna tabla**.
* Formato: datos. Sin veredicto, sin conclusiones.

## 0. Construccion, declarada antes de medir

| #  | cantidad | valor |
|---|---|---|
| 1 | filtro | `filter_gene_min_cpm_cell = 5.0` del preset `vcc2026`, leido de `EvalConfig`; pasa si CPM medio en la referencia `>` umbral |
| 2 | referencia en Parte 1 | los 18,400 controles de cada contexto, todos |
| 3 | comprobacion del filtro | `compute_de` del preset, referencia = 18,400 controles, blanco = 400 controles; genes probados = genes de esta cuenta, gen a gen |
| 4 | cantidades por gen (Parte 1) | CPM medio, tasa de deteccion, Fano de cuentas crudas |
| 5 | control de 1.4 | dos mitades de 9,200 del MISMO contexto: aleatoria (semilla) y por guia `ntc_id` (23/23, determinista) |
| 6 | R(g,L) en Jiang | `compute_de` + `prep_de_side` del preset, `p_adj < 0.05`, sin el propio gen; comprobado gen a gen contra `venn_ed.json` |
| 7 | unidad de la Parte 2 | (g, gen, L1, L2), L1 ≠ L2, gen ∈ R(g,L1); y = 1 si gen ∉ R(g,L2) |
| 8 | predictores (todos del control de L2), orden | p1 pasa filtro → p2 CPM medio (deciles) → p3 deteccion (deciles) → p4 potencia (log n, log UMI mediana); fijado antes de medir |
| 9 | medida A | fraccion de las y = 1 cuyo gen no pasa el filtro de L2 (determinista) |
| 10 | medida B | D² de una logistica (devianza explicada), SOLO y ACUMULADO |
| 11 | nulo | 20 permutaciones uniformes de las 6 lineas: las estadisticas de control de L2 se toman de π(L2); semillas 20281100 + k |
| 12 | anclaje | (obs − media nulo) / (1 − media nulo) |
| 13 | igualacion en potencia (Jiang) | control submuestreado a n_min y TODAS las celulas adelgazadas binomialmente a la UMI mediana minima; se re-corre la puerta |
| 14 | VCC 2025 | la poblacion real de `vcc2025_suelo` (semillas 20280201/20280202 reusadas a proposito), verificada por sha256 |

### 0.1 Semillas: familia nueva 20281000+, comprobadas EN CODIGO

| uso | semillas |
|---|---|
| p1_verif_filtro | 20281001, 20281002, 20281003, 20281004, 20281005, 20281006 |
| p1_mitades | 20281021, 20281022, 20281023, 20281024, 20281025, 20281026 |
| p2_nulo | 20281100, 20281101, 20281102, 20281103, 20281104, 20281105, 20281106, 20281107, 20281108, 20281109, 20281110, 20281111, 20281112, 20281113, 20281114, 20281115, 20281116, 20281117, 20281118, 20281119 |
| p2_ig_control | 20281201, 20281202, 20281203, 20281204, 20281205, 20281206 |
| p2_ig_adelgaza | 20281211, 20281212, 20281213, 20281214, 20281215, 20281216 |

Choques internos: 0; con rangos gastados: 0; con la enumeracion de `noche/comun.py`: 0; con el barrido de 740 enteros 2026xxxx–2028xxxx del codigo y las salidas (maximo previo de la familia 2028: 20280456): 0. **Sin choques: si.** Guardadas en `/home/ubuntu/vcc/modelo/out/seleccion_control.json` y en `modelo/src/seleccion/semillas.py`.

---

## PARTE 1 — SUPERFICIE DE ENTREGA (D/E/F) Y DE VALIDACION (A/B/C), SOLO DE SUS CONTROLES

Toda cifra de esta parte sale de los 18,400 controles de cada contexto. No hay verdad en ella y no la necesita: el filtro del preset lo decide la referencia, y la referencia en la puntuacion real ES este pool.

### 1.1 Que genes pasan el filtro, por contexto

Filtro: CPM medio en el control `> 5.0` (estricto), de `EvalConfig.from_preset('vcc2026')`. Columna *verif*: el conjunto de genes que `compute_de` del preset PRUEBA con referencia = los 18,400 controles y blanco = 400 controles (semilla en la columna), comparado gen a gen con esta cuenta.

| #  | contexto | particion | n control | UMI mediana/celula | nnz mediana/celula | **# pasan** | de 18,533 | # probados por compute_de | identicos gen a gen | semilla verif |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **A** | val | 18,400 | 20,109 | 6,147 | **9,929** | 53.57 % | 9,929 | **si** | 20281001 |
| 2 | **B** | val | 18,400 | 19,946 | 5,756 | **9,626** | 51.94 % | 9,626 | **si** | 20281002 |
| 3 | **C** | val | 18,400 | 20,034 | 6,006 | **10,124** | 54.63 % | 10,124 | **si** | 20281003 |
| 4 | **D** | final | 18,400 | 20,236 | 5,892 | **9,633** | 51.98 % | 9,633 | **si** | 20281004 |
| 5 | **E** | final | 18,400 | 19,681 | 5,634 | **9,743** | 52.57 % | 9,743 | **si** | 20281005 |
| 6 | **F** | final | 18,400 | 19,845 | 5,630 | **9,537** | 51.46 % | 9,537 | **si** | 20281006 |

Las listas completas: `/home/ubuntu/vcc/modelo/out/seleccion_filtrables_listas.json`. Orden de genes identico en los 6 contextos y a `gene_names.csv`: **si**.

### 1.2 🔴 El Venn de los filtrables, dentro de cada particion

**D/E/F**

| #  | region | # genes | frac de la union | genes (si <= 40) |
|---|---|---|---|---|
| 1 | solo D | 1,183 | 10.20 % | 1183 genes (lista en `seleccion_filtrables_listas.json`) |
| 2 | solo E | 483 | 4.17 % | 483 genes (lista en `seleccion_filtrables_listas.json`) |
| 3 | solo F | 379 | 3.27 % | 379 genes (lista en `seleccion_filtrables_listas.json`) |
| 4 | D∩E sin F | 390 | 3.36 % | 390 genes (lista en `seleccion_filtrables_listas.json`) |
| 5 | D∩F sin E | 288 | 2.48 % | 288 genes (lista en `seleccion_filtrables_listas.json`) |
| 6 | E∩F sin D | 1,098 | 9.47 % | 1098 genes (lista en `seleccion_filtrables_listas.json`) |
| 7 | **D∩E∩F** | 7,772 | 67.04 % | 7772 genes (lista en `seleccion_filtrables_listas.json`) |

| pasan en | # genes | frac de 18,533 |
|---|---|---|
| los 3 | 7,772 | 41.94 % |
| exactamente 2 | 1,776 | 9.58 % |
| exactamente 1 | 2,045 | 11.03 % |
| ninguno | 6,940 | 37.45 % |
| union | 11,593 | 62.55 % |

**A/B/C**

| #  | region | # genes | frac de la union | genes (si <= 40) |
|---|---|---|---|---|
| 1 | solo A | 989 | 8.27 % | 989 genes (lista en `seleccion_filtrables_listas.json`) |
| 2 | solo B | 459 | 3.84 % | 459 genes (lista en `seleccion_filtrables_listas.json`) |
| 3 | solo C | 786 | 6.57 % | 786 genes (lista en `seleccion_filtrables_listas.json`) |
| 4 | A∩B sin C | 389 | 3.25 % | 389 genes (lista en `seleccion_filtrables_listas.json`) |
| 5 | A∩C sin B | 560 | 4.68 % | 560 genes (lista en `seleccion_filtrables_listas.json`) |
| 6 | B∩C sin A | 787 | 6.58 % | 787 genes (lista en `seleccion_filtrables_listas.json`) |
| 7 | **A∩B∩C** | 7,991 | 66.81 % | 7991 genes (lista en `seleccion_filtrables_listas.json`) |

| pasan en | # genes | frac de 18,533 |
|---|---|---|
| los 3 | 7,991 | 43.12 % |
| exactamente 2 | 1,736 | 9.37 % |
| exactamente 1 | 2,234 | 12.05 % |
| ninguno | 6,572 | 35.46 % |
| union | 11,961 | 64.54 % |

### 1.3 Cantidades por gen y contexto, todas del control

Cuantiles sobre genes. CPM medio = media de x/L·1e6; deteccion = fraccion de celulas con cuenta > 0; Fano = var/media de las CUENTAS CRUDAS (depende de la profundidad: no es comparable entre contextos de profundidad distinta sin esa salvedad).

**todos los genes con media > 0**

| ctx | # genes | CPM p10 | CPM p50 | CPM p90 | det p10 | det p50 | det p90 | Fano p10 | Fano p50 | Fano p90 |
|---|---|---|---|---|---|---|---|---|---|---|
| **A** | 16,135 | 0.00928 | 19.3 | 161 | 0.000163 | 0.299 | 0.894 | 1 | 1.31 | 2.27 |
| **B** | 15,680 | 0.0189 | 16.1 | 163 | 0.000217 | 0.252 | 0.892 | 1 | 1.18 | 2 |
| **C** | 16,216 | 0.018 | 17.1 | 160 | 0.00038 | 0.279 | 0.896 | 1 | 1.22 | 2.28 |
| **D** | 15,558 | 0.0142 | 20.7 | 152 | 0.000109 | 0.312 | 0.87 | 1 | 1.29 | 2.5 |
| **E** | 15,824 | 0.0129 | 15.4 | 155 | 0.000217 | 0.252 | 0.886 | 1 | 1.23 | 2.39 |
| **F** | 16,371 | 0.0126 | 12.3 | 159 | 0.000217 | 0.206 | 0.887 | 1 | 1.23 | 2.45 |

**solo los que pasan el filtro del contexto**

| ctx | # genes | CPM p10 | CPM p50 | CPM p90 | det p10 | det p50 | det p90 | Fano p10 | Fano p50 | Fano p90 |
|---|---|---|---|---|---|---|---|---|---|---|
| **A** | 9,929 | 11.6 | 57.4 | 220 | 0.195 | 0.624 | 0.937 | 1.18 | 1.5 | 2.7 |
| **B** | 9,626 | 10.2 | 50.2 | 228 | 0.17 | 0.576 | 0.935 | 1.09 | 1.32 | 2.33 |
| **C** | 10,124 | 10.1 | 49.6 | 218 | 0.178 | 0.58 | 0.943 | 1.12 | 1.41 | 2.75 |
| **D** | 9,633 | 12.3 | 56.1 | 203 | 0.199 | 0.608 | 0.915 | 1.16 | 1.5 | 3 |
| **E** | 9,743 | 9.85 | 47.2 | 213 | 0.173 | 0.561 | 0.936 | 1.12 | 1.41 | 2.77 |
| **F** | 9,537 | 9.78 | 49.8 | 230 | 0.17 | 0.576 | 0.938 | 1.12 | 1.45 | 2.95 |

**En que se diferencian los tres contextos de cada particion** (pares; sobre los genes que pasan en al menos uno del par)

| particion | par | # genes (union) | jaccard filtrables | pasan en 1º no 2º | pasan en 2º no 1º | Spearman log2 CPM | mediana |log2 ratio CPM| | p90 |log2 ratio| | mediana dif. deteccion (1º−2º) | ratio UMI mediana (1º/2º) |
|---|---|---|---|---|---|---|---|---|---|---|
| D/E/F | D–E | 11,214 | 0.7278 | 1,471 | 1,581 | 0.4341 | 1.278 | 8.947 | +0.0214 | 1.028 |
| D/E/F | D–F | 11,110 | 0.7255 | 1,573 | 1,477 | 0.4544 | 1.264 | 8.212 | +0.0147 | 1.020 |
| D/E/F | E–F | 10,410 | 0.8521 | 873 | 667 | 0.7932 | 0.650 | 2.895 | +0.0003 | 0.992 |
| A/B/C | A–B | 11,175 | 0.7499 | 1,549 | 1,246 | 0.5365 | 1.130 | 7.489 | +0.0495 | 1.008 |
| A/B/C | A–C | 11,502 | 0.7434 | 1,378 | 1,573 | 0.4544 | 1.332 | 8.042 | +0.0279 | 1.004 |
| A/B/C | B–C | 10,972 | 0.8000 | 848 | 1,346 | 0.6129 | 0.945 | 4.878 | -0.0114 | 0.996 |

🔴 **Intermedios D/E/F: 20 genes que pasan en algun contexto pero no en los 3**, equiespaciados por CPM en D (3,821 candidatos)

| #  | gen | CPM D | det D | Fano D | pasa D | CPM E | det E | Fano E | pasa E | CPM F | det F | Fano F | pasa F |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `MMP24` | 0.000 | 0.0000 | n/d | no | 3.974 | 0.0747 | 1.129 | no | 13.545 | 0.2156 | 1.490 | si |
| 2 | `RIPOR3` | 0.000 | 0.0000 | n/d | no | 13.385 | 0.1837 | 1.720 | si | 0.892 | 0.0106 | 1.807 | no |
| 3 | `SAMD14` | 0.002 | 0.0001 | 1.000 | no | 9.013 | 0.1621 | 1.175 | si | 1.180 | 0.0222 | 1.075 | no |
| 4 | `RBPMS` | 0.009 | 0.0002 | 1.000 | no | 265.176 | 0.9582 | 2.976 | si | 293.825 | 0.9579 | 3.522 | si |
| 5 | `LRP4` | 0.034 | 0.0008 | 1.234 | no | 9.975 | 0.1858 | 1.143 | si | 7.102 | 0.1276 | 1.114 | si |
| 6 | `FMN1` | 0.097 | 0.0014 | 1.343 | no | 16.939 | 0.2802 | 1.188 | si | 23.390 | 0.3497 | 1.404 | si |
| 7 | `BEX1` | 0.252 | 0.0040 | 1.301 | no | 9.399 | 0.1109 | 3.204 | si | 135.358 | 0.7029 | 4.588 | si |
| 8 | `OLFML3` | 0.803 | 0.0100 | 2.071 | no | 0.415 | 0.0080 | 1.218 | no | 8.365 | 0.1457 | 1.282 | si |
| 9 | `TFAP2A` | 1.840 | 0.0298 | 1.715 | no | 4.939 | 0.0866 | 1.346 | no | 49.551 | 0.5593 | 1.651 | si |
| 10 | `CASTOR3` | 3.248 | 0.0515 | 1.080 | no | 1.980 | 0.0403 | 1.031 | no | 7.487 | 0.1308 | 1.088 | si |
| 11 | `APOC1` | 5.509 | 0.0951 | 1.139 | si | 0.001 | 0.0001 | 1.000 | no | 0.016 | 0.0003 | 1.000 | no |
| 12 | `MROH6` | 7.550 | 0.1277 | 1.148 | si | 0.416 | 0.0090 | 0.991 | no | 0.413 | 0.0074 | 0.993 | no |
| 13 | `ARSG` | 10.561 | 0.1633 | 1.209 | si | 0.675 | 0.0127 | 1.093 | no | 3.572 | 0.0667 | 1.087 | no |
| 14 | `ISL2` | 14.702 | 0.2414 | 1.157 | si | 4.831 | 0.0922 | 1.057 | no | 4.071 | 0.0781 | 1.073 | no |
| 15 | `QPRT` | 20.157 | 0.2630 | 1.628 | si | 0.512 | 0.0109 | 1.101 | no | 2.227 | 0.0418 | 1.062 | no |
| 16 | `SEMA6B` | 28.546 | 0.1979 | 3.761 | si | 2.592 | 0.0426 | 2.467 | no | 5.063 | 0.0739 | 1.488 | si |
| 17 | `CD244` | 41.177 | 0.4691 | 1.417 | si | 0.000 | 0.0000 | n/d | no | 0.010 | 0.0001 | 1.000 | no |
| 18 | `FBXO41` | 68.018 | 0.6653 | 1.521 | si | 6.808 | 0.1259 | 1.087 | si | 0.179 | 0.0033 | 1.061 | no |
| 19 | `EMILIN2` | 141.626 | 0.8570 | 1.986 | si | 4.795 | 0.0844 | 1.574 | no | 92.758 | 0.7009 | 2.775 | si |
| 20 | `PRTN3` | 20408.962 | 1.0000 | 169.177 | si | 0.000 | 0.0000 | n/d | no | 0.188 | 0.0031 | 1.065 | no |

🔴 **Intermedios A/B/C: 20 genes que pasan en algun contexto pero no en los 3**, equiespaciados por CPM en A (3,970 candidatos)

| #  | gen | CPM A | det A | Fano A | pasa A | CPM B | det B | Fano B | pasa B | CPM C | det C | Fano C | pasa C |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `WNT3A` | 0.000 | 0.0000 | n/d | no | 0.000 | 0.0000 | n/d | no | 10.297 | 0.1798 | 1.145 | si |
| 2 | `VSNL1` | 0.000 | 0.0000 | n/d | no | 0.000 | 0.0000 | n/d | no | 123.444 | 0.8156 | 2.598 | si |
| 3 | `SPTY2D1OS` | 0.003 | 0.0001 | 1.000 | no | 4.999 | 0.0910 | 1.029 | no | 5.834 | 0.1108 | 1.060 | si |
| 4 | `SELENBP1` | 0.012 | 0.0002 | 1.000 | no | 92.147 | 0.7447 | 1.704 | si | 0.042 | 0.0009 | 1.110 | no |
| 5 | `LCN2` | 0.046 | 0.0010 | 1.199 | no | 0.172 | 0.0037 | 1.025 | no | 169.039 | 0.8461 | 8.241 | si |
| 6 | `SYT12` | 0.136 | 0.0024 | 1.122 | no | 97.577 | 0.7180 | 2.165 | si | 6.237 | 0.1158 | 1.167 | si |
| 7 | `ETV4` | 0.381 | 0.0067 | 1.264 | no | 5.351 | 0.0877 | 1.436 | si | 40.526 | 0.5136 | 1.464 | si |
| 8 | `PCDHB16` | 0.953 | 0.0167 | 1.310 | no | 11.970 | 0.1870 | 1.372 | si | 0.053 | 0.0008 | 1.133 | no |
| 9 | `ACTRT3` | 2.033 | 0.0344 | 1.332 | no | 10.359 | 0.1686 | 1.265 | si | 6.272 | 0.1114 | 1.187 | si |
| 10 | `TJP3` | 3.649 | 0.0683 | 1.240 | no | 0.067 | 0.0009 | 0.999 | no | 6.829 | 0.1235 | 1.163 | si |
| 11 | `ZNF674` | 5.340 | 0.1041 | 1.095 | si | 6.436 | 0.1184 | 1.075 | si | 3.688 | 0.0750 | 1.060 | no |
| 12 | `RPH3AL` | 7.108 | 0.1335 | 1.102 | si | 1.490 | 0.0293 | 1.021 | no | 1.256 | 0.0249 | 1.038 | no |
| 13 | `FRRS1` | 9.450 | 0.1710 | 1.167 | si | 4.152 | 0.0794 | 1.056 | no | 41.138 | 0.5329 | 1.296 | si |
| 14 | `NTNG2` | 12.243 | 0.1899 | 1.416 | si | 0.024 | 0.0006 | 0.999 | no | 3.643 | 0.0666 | 1.228 | no |
| 15 | `FCMR` | 16.421 | 0.2577 | 1.380 | si | 0.548 | 0.0105 | 1.059 | no | 0.798 | 0.0171 | 0.995 | no |
| 16 | `XAF1` | 22.095 | 0.3393 | 1.410 | si | 0.010 | 0.0002 | 1.000 | no | 20.667 | 0.2963 | 1.515 | si |
| 17 | `NSD3` | 30.049 | 0.4359 | 1.156 | si | 2.299 | 0.0426 | 1.009 | no | 2.483 | 0.0511 | 1.023 | no |
| 18 | `RGS14` | 44.463 | 0.5369 | 1.509 | si | 1.637 | 0.0275 | 1.048 | no | 6.306 | 0.1177 | 1.339 | si |
| 19 | `CMPK2` | 90.243 | 0.7202 | 2.166 | si | 0.002 | 0.0001 | 1.000 | no | 159.239 | 0.5821 | 14.311 | si |
| 20 | `HIST1H1B` | 5731.277 | 0.9949 | 208.102 | si | 906.116 | 0.7751 | 117.604 | si | 0.082 | 0.0010 | 9.921 | no |

### 1.4 🔴 Por cuanto fallan: genes que pasan en X y NO en Y, su CPM en Y respecto de 5.0

Columnas de distribucion: fraccion de esos genes con CPM_Y/5.0 en cada intervalo. Si se agolpan junto a 1, el filtro es fragil; si estan lejos, es estructural.

🔴 **CONTROL** (filas `mitad`): el MISMO conteo entre dos mitades del MISMO contexto (9,200 vs 9,200), que no difieren en nada salvo muestreo — (a) aleatoria de celulas, (b) por guia `ntc_id` (23 vs 23). Defecto: cada mitad tiene la mitad de celulas que un contexto, asi que su ruido de muestreo es ~√2 mayor que el de un contexto entero; el control ACOTA POR ARRIBA la fragilidad por muestreo.

| comparacion | # genes | CPM_Y/5 en (0.95,1] | CPM_Y/5 en (0.9,0.95] | CPM_Y/5 en (0.8,0.9] | CPM_Y/5 en (0.5,0.8] | CPM_Y/5 en (0.2,0.5] | CPM_Y/5 en ≤0.2 | mediana CPM_Y/5 | mediana CPM_X/5 |
|---|---|---|---|---|---|---|---|---|---|
| D/E/F: pasa **D**, no **E** | 1,471 | 3.1 % | 2.9 % | 4.9 % | 15.3 % | 17.5 % | 56.4 % | 0.124 | 4.426 |
| D/E/F: pasa **D**, no **F** | 1,573 | 2.9 % | 2.0 % | 4.6 % | 16.9 % | 19.8 % | 53.7 % | 0.165 | 4.381 |
| D/E/F: pasa **E**, no **D** | 1,581 | 1.0 % | 0.9 % | 2.3 % | 9.0 % | 13.7 % | 73.1 % | 0.020 | 5.061 |
| D/E/F: pasa **E**, no **F** | 873 | 4.9 % | 3.6 % | 7.9 % | 23.1 % | 22.1 % | 38.4 % | 0.362 | 2.117 |
| D/E/F: pasa **F**, no **D** | 1,477 | 0.9 % | 1.0 % | 2.3 % | 9.6 % | 13.5 % | 72.6 % | 0.022 | 4.509 |
| D/E/F: pasa **F**, no **E** | 667 | 4.3 % | 4.5 % | 7.0 % | 14.1 % | 18.0 % | 52.0 % | 0.178 | 2.150 |
| A/B/C: pasa **A**, no **B** | 1,549 | 3.2 % | 2.5 % | 5.1 % | 16.3 % | 21.9 % | 51.1 % | 0.187 | 3.364 |
| A/B/C: pasa **A**, no **C** | 1,378 | 2.2 % | 2.5 % | 5.7 % | 17.5 % | 20.5 % | 51.6 % | 0.178 | 3.733 |
| A/B/C: pasa **B**, no **A** | 1,246 | 0.9 % | 1.8 % | 2.4 % | 8.6 % | 13.2 % | 73.1 % | 0.019 | 5.396 |
| A/B/C: pasa **B**, no **C** | 848 | 2.7 % | 3.2 % | 6.7 % | 16.9 % | 18.8 % | 51.8 % | 0.173 | 2.908 |
| A/B/C: pasa **C**, no **A** | 1,573 | 1.3 % | 1.8 % | 3.0 % | 10.2 % | 14.0 % | 69.7 % | 0.033 | 5.297 |
| A/B/C: pasa **C**, no **B** | 1,346 | 3.4 % | 2.8 % | 5.8 % | 16.7 % | 20.4 % | 50.8 % | 0.192 | 3.347 |
| mitad A (aleatoria), ambos sentidos | 46 | 78.3 % | 15.2 % | 6.5 % | 0.0 % | 0.0 % | 0.0 % | 0.973 | 1.046 |
| mitad A (por_guia), ambos sentidos | 40 | 77.5 % | 20.0 % | 2.5 % | 0.0 % | 0.0 % | 0.0 % | 0.974 | 1.041 |
| mitad B (aleatoria), ambos sentidos | 59 | 62.7 % | 30.5 % | 6.8 % | 0.0 % | 0.0 % | 0.0 % | 0.966 | 1.027 |
| mitad B (por_guia), ambos sentidos | 68 | 66.2 % | 26.5 % | 7.4 % | 0.0 % | 0.0 % | 0.0 % | 0.967 | 1.036 |
| mitad C (aleatoria), ambos sentidos | 46 | 69.6 % | 21.7 % | 4.3 % | 4.3 % | 0.0 % | 0.0 % | 0.972 | 1.028 |
| mitad C (por_guia), ambos sentidos | 44 | 70.5 % | 25.0 % | 2.3 % | 2.3 % | 0.0 % | 0.0 % | 0.974 | 1.022 |
| mitad D (aleatoria), ambos sentidos | 29 | 82.8 % | 10.3 % | 3.4 % | 3.4 % | 0.0 % | 0.0 % | 0.979 | 1.031 |
| mitad D (por_guia), ambos sentidos | 57 | 70.2 % | 15.8 % | 10.5 % | 3.5 % | 0.0 % | 0.0 % | 0.969 | 1.037 |
| mitad E (aleatoria), ambos sentidos | 53 | 69.8 % | 24.5 % | 3.8 % | 1.9 % | 0.0 % | 0.0 % | 0.964 | 1.028 |
| mitad E (por_guia), ambos sentidos | 55 | 70.9 % | 29.1 % | 0.0 % | 0.0 % | 0.0 % | 0.0 % | 0.977 | 1.027 |
| mitad F (aleatoria), ambos sentidos | 53 | 73.6 % | 22.6 % | 3.8 % | 0.0 % | 0.0 % | 0.0 % | 0.974 | 1.026 |
| mitad F (por_guia), ambos sentidos | 58 | 74.1 % | 15.5 % | 10.3 % | 0.0 % | 0.0 % | 0.0 % | 0.968 | 1.024 |

🔴 **Intermedios: 20 de los 1,471 genes que pasan en D y no en E**, equiespaciados por CPM en E

| #  | gen | CPM D | CPM E | CPM_E/5 | CPM F | det D | det E | pasa F |
|---|---|---|---|---|---|---|---|---|
| 1 | `PIK3R6` | 12.897 | 0.000 | 0.000 | 0.052 | 0.1642 | 0.0000 | no |
| 2 | `NPIPB11` | 7.461 | 0.000 | 0.000 | 0.000 | 0.1399 | 0.0000 | no |
| 3 | `CNR2` | 26.210 | 0.002 | 0.000 | 0.003 | 0.2500 | 0.0001 | no |
| 4 | `SUSD3` | 87.259 | 0.007 | 0.001 | 0.334 | 0.7515 | 0.0001 | no |
| 5 | `IRAK3` | 148.561 | 0.018 | 0.004 | 2.888 | 0.8846 | 0.0003 | no |
| 6 | `OSM` | 26.797 | 0.040 | 0.008 | 0.016 | 0.2824 | 0.0007 | no |
| 7 | `TAS2R30` | 18.453 | 0.080 | 0.016 | 0.008 | 0.2960 | 0.0018 | no |
| 8 | `CENPV` | 169.827 | 0.159 | 0.032 | 131.031 | 0.8611 | 0.0026 | si |
| 9 | `ZNF512` | 44.054 | 0.295 | 0.059 | 76.046 | 0.5436 | 0.0046 | si |
| 10 | `CCNJ` | 12.996 | 0.502 | 0.100 | 0.102 | 0.2232 | 0.0099 | no |
| 11 | `JMJD8` | 51.271 | 0.767 | 0.153 | 22.275 | 0.5875 | 0.0160 | si |
| 12 | `IL32` | 36.630 | 1.081 | 0.216 | 5.032 | 0.2067 | 0.0191 | si |
| 13 | `PDXDC1` | 48.638 | 1.532 | 0.306 | 0.317 | 0.5996 | 0.0310 | no |
| 14 | `GNRH1` | 19.404 | 1.974 | 0.395 | 2.838 | 0.2835 | 0.0384 | no |
| 15 | `ZNF34` | 18.910 | 2.474 | 0.495 | 5.670 | 0.2950 | 0.0507 | si |
| 16 | `NEK3` | 66.175 | 2.916 | 0.583 | 20.796 | 0.6760 | 0.0565 | si |
| 17 | `ARID5A` | 32.196 | 3.485 | 0.697 | 8.341 | 0.4360 | 0.0674 | si |
| 18 | `MXD3` | 5.899 | 4.015 | 0.803 | 4.666 | 0.1098 | 0.0748 | no |
| 19 | `ZNF81` | 6.921 | 4.566 | 0.913 | 3.122 | 0.1257 | 0.0903 | no |
| 20 | `PAFAH2` | 9.675 | 4.994 | 0.999 | 2.190 | 0.1727 | 0.0983 | no |

### 1.5 El Venn entre las dos particiones: ¿los filtrables de A/B/C son los de D/E/F?

Tabla cruzada: en cuantos de A/B/C pasa el gen (filas) contra en cuantos de D/E/F (columnas).

|  | 0 de D/E/F | 1 de D/E/F | 2 de D/E/F | 3 de D/E/F | total |
|---|---|---|---|---|---|
| **0** de A/B/C | 5,717 | 656 | 176 | 23 | 6,572 |
| **1** de A/B/C | 974 | 791 | 374 | 95 | 2,234 |
| **2** de A/B/C | 202 | 473 | 612 | 449 | 1,736 |
| **3** de A/B/C | 47 | 125 | 614 | 7,205 | 7,991 |
| total | 6,940 | 2,045 | 1,776 | 7,772 | 18,533 |

| conjunto | A/B/C | D/E/F | comun | solo A/B/C | solo D/E/F | jaccard | genes solo A/B/C (<= 25) | genes solo D/E/F (<= 25) |
|---|---|---|---|---|---|---|---|---|
| nucleo (pasa en los 3) | 7,991 | 7,772 | 7,205 | 786 | 567 | 0.8419 | 786 genes (lista en `seleccion_filtrables_listas.json`) | 567 genes (lista en `seleccion_filtrables_listas.json`) |
| union (pasa en >= 1) | 11,961 | 11,593 | 10,738 | 1,223 | 855 | 0.8379 | 1223 genes (lista en `seleccion_filtrables_listas.json`) | 855 genes (lista en `seleccion_filtrables_listas.json`) |

**Jaccard de filtrables, contexto a contexto (A/B/C × D/E/F)**

|  | D | E | F |
|---|---|---|---|
| **A** | 0.7787 | 0.7458 | 0.7454 |
| **B** | 0.7338 | 0.8300 | 0.8483 |
| **C** | 0.7523 | 0.7794 | 0.7807 |

Referencia dentro de particion: D–E 0.7278, D–F 0.7255, E–F 0.8521, A–B 0.7499, A–C 0.7434, B–C 0.8000

---

## PARTE 2 — JIANG Y VCC 2025 (CON VERDAD): EL MECANISMO

🔴 De esta parte viaja solo el MECANISMO. Ninguna fraccion de aqui es una cifra de D/E/F.

### 2.0 Las poblaciones, verificadas

| #  | linea | huella real = venn_ed | R(g,L) = venn_ed gen a gen | n control | UMI mediana | # pasan filtro | Σ|R| crudo | p adelgaza | n control ig. | UMI mediana ig. | # pasan ig. | Σ|R| igualado | semillas ig. (control / adelgaza) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **A549** | si | 33/33 | 1,379 | 3,319 | 11,733 | 2537 | 1.0000 | 1,379 | 3,319 | 11,733 | 2537 | 20281201 / 20281211 |
| 2 | **BXPC3** | si | 33/33 | 3,228 | 4,754 | 12,005 | 4642 | 0.6981 | 1,379 | 3,301 | 11,980 | 2519 | 20281202 / 20281212 |
| 3 | **HAP1** | si | 33/33 | 3,172 | 5,231 | 12,041 | 817 | 0.6345 | 1,379 | 3,377 | 12,039 | 559 | 20281203 / 20281213 |
| 4 | **HT29** | si | 33/33 | 2,736 | 4,396 | 11,181 | 2505 | 0.7550 | 1,379 | 3,271 | 11,203 | 1824 | 20281204 / 20281214 |
| 5 | **K562** | si | 33/33 | 1,398 | 4,302 | 11,750 | 1042 | 0.7715 | 1,379 | 3,335 | 11,727 | 832 | 20281205 / 20281215 |
| 6 | **MCF7** | si | 33/33 | 2,669 | 6,203 | 11,534 | 12867 | 0.5351 | 1,379 | 3,315 | 11,523 | 6779 | 20281206 / 20281216 |

Igualacion: n_min = 1,379 controles, D_min = 3,319 UMI mediana/celula. Genes del espacio Jiang: 34,025.

VCC 2025 (`test`): huella del real = la de `real_2025_test_meta.json`: **si**; control 18,400 celulas, UMI mediana 54,637; filtrables 10,776 de 18,080; Σ|R(g)| = 173,473 sobre 100 perturbaciones.

Chequeo de construccion: genes de R(g,L1) que NO pasan el filtro de L1 — crudo **0**, igualado **0** (tiene que ser 0).

Unidades (g, gen, L1, L2) — CRUDO: **122,050** filas, de ellas **109,830** especificas (gen en R(g,L1), no en R(g,L2)) = 89.99 %.

Unidades (g, gen, L1, L2) — IGUALADO EN POTENCIA: **75,250** filas, de ellas **66,868** especificas (gen en R(g,L1), no en R(g,L2)) = 88.86 %.

### 2.1 De los genes en R(g,L1) y no en R(g,L2): que fraccion NO PASA el filtro de L2

Sin biologia: un gen que no pasa el filtro de L2 no se PRUEBA en L2, asi que no puede estar en R(g,L2).

**CRUDO** — agregado: **14.63 %** de 109,830; nulo (20 permutaciones del contexto de las estadisticas de control) 12.01 % [11.06, 12.99]; anclado **0.0299**.

| #  | L1 | L2 | # en R(g,L1) | # no en R(g,L2) | # de esos que no pasan filtro L2 | fraccion |
|---|---|---|---|---|---|---|
| 1 | A549 | BXPC3 | 2,537 | 1,822 | 218 | 12.0 % |
| 2 | A549 | HAP1 | 2,537 | 2,309 | 347 | 15.0 % |
| 3 | A549 | HT29 | 2,537 | 2,115 | 343 | 16.2 % |
| 4 | A549 | K562 | 2,537 | 2,303 | 420 | 18.2 % |
| 5 | A549 | MCF7 | 2,537 | 2,036 | 351 | 17.2 % |
| 6 | BXPC3 | A549 | 4,642 | 3,927 | 641 | 16.3 % |
| 7 | BXPC3 | HAP1 | 4,642 | 4,298 | 936 | 21.8 % |
| 8 | BXPC3 | HT29 | 4,642 | 3,968 | 771 | 19.4 % |
| 9 | BXPC3 | K562 | 4,642 | 4,288 | 950 | 22.2 % |
| 10 | BXPC3 | MCF7 | 4,642 | 3,995 | 834 | 20.9 % |
| 11 | HAP1 | A549 | 817 | 589 | 63 | 10.7 % |
| 12 | HAP1 | BXPC3 | 817 | 473 | 50 | 10.6 % |
| 13 | HAP1 | HT29 | 817 | 440 | 81 | 18.4 % |
| 14 | HAP1 | K562 | 817 | 532 | 81 | 15.2 % |
| 15 | HAP1 | MCF7 | 817 | 570 | 76 | 13.3 % |
| 16 | HT29 | A549 | 2,505 | 2,083 | 196 | 9.4 % |
| 17 | HT29 | BXPC3 | 2,505 | 1,831 | 124 | 6.8 % |
| 18 | HT29 | HAP1 | 2,505 | 2,128 | 358 | 16.8 % |
| 19 | HT29 | K562 | 2,505 | 2,102 | 340 | 16.2 % |
| 20 | HT29 | MCF7 | 2,505 | 2,084 | 223 | 10.7 % |
| 21 | K562 | A549 | 1,042 | 808 | 118 | 14.6 % |
| 22 | K562 | BXPC3 | 1,042 | 688 | 104 | 15.1 % |
| 23 | K562 | HAP1 | 1,042 | 757 | 150 | 19.8 % |
| 24 | K562 | HT29 | 1,042 | 639 | 149 | 23.3 % |
| 25 | K562 | MCF7 | 1,042 | 784 | 145 | 18.5 % |
| 26 | MCF7 | A549 | 12,867 | 12,366 | 1,512 | 12.2 % |
| 27 | MCF7 | BXPC3 | 12,867 | 12,220 | 1,354 | 11.1 % |
| 28 | MCF7 | HAP1 | 12,867 | 12,620 | 1,523 | 12.1 % |
| 29 | MCF7 | HT29 | 12,867 | 12,446 | 1,643 | 13.2 % |
| 30 | MCF7 | K562 | 12,867 | 12,609 | 1,972 | 15.6 % |

**IGUALADO** — agregado: **13.86 %** de 66,868; nulo (20 permutaciones del contexto de las estadisticas de control) 11.38 % [10.55, 12.24]; anclado **0.0280**.

| #  | L1 | L2 | # en R(g,L1) | # no en R(g,L2) | # de esos que no pasan filtro L2 | fraccion |
|---|---|---|---|---|---|---|
| 1 | A549 | BXPC3 | 2,537 | 2,002 | 213 | 10.6 % |
| 2 | A549 | HAP1 | 2,537 | 2,343 | 358 | 15.3 % |
| 3 | A549 | HT29 | 2,537 | 2,189 | 343 | 15.7 % |
| 4 | A549 | K562 | 2,537 | 2,343 | 416 | 17.8 % |
| 5 | A549 | MCF7 | 2,537 | 2,232 | 359 | 16.1 % |
| 6 | BXPC3 | A549 | 2,519 | 1,984 | 338 | 17.0 % |
| 7 | BXPC3 | HAP1 | 2,519 | 2,270 | 534 | 23.5 % |
| 8 | BXPC3 | HT29 | 2,519 | 2,089 | 395 | 18.9 % |
| 9 | BXPC3 | K562 | 2,519 | 2,261 | 492 | 21.8 % |
| 10 | BXPC3 | MCF7 | 2,519 | 2,233 | 434 | 19.4 % |
| 11 | HAP1 | A549 | 559 | 365 | 34 | 9.3 % |
| 12 | HAP1 | BXPC3 | 559 | 310 | 27 | 8.7 % |
| 13 | HAP1 | HT29 | 559 | 258 | 44 | 17.1 % |
| 14 | HAP1 | K562 | 559 | 338 | 42 | 12.4 % |
| 15 | HAP1 | MCF7 | 559 | 400 | 42 | 10.5 % |
| 16 | HT29 | A549 | 1,824 | 1,476 | 113 | 7.7 % |
| 17 | HT29 | BXPC3 | 1,824 | 1,394 | 75 | 5.4 % |
| 18 | HT29 | HAP1 | 1,824 | 1,523 | 228 | 15.0 % |
| 19 | HT29 | K562 | 1,824 | 1,517 | 226 | 14.9 % |
| 20 | HT29 | MCF7 | 1,824 | 1,588 | 132 | 8.3 % |
| 21 | K562 | A549 | 832 | 638 | 85 | 13.3 % |
| 22 | K562 | BXPC3 | 832 | 574 | 79 | 13.8 % |
| 23 | K562 | HAP1 | 832 | 611 | 111 | 18.2 % |
| 24 | K562 | HT29 | 832 | 525 | 114 | 21.7 % |
| 25 | K562 | MCF7 | 832 | 664 | 107 | 16.1 % |
| 26 | MCF7 | A549 | 6,779 | 6,474 | 736 | 11.4 % |
| 27 | MCF7 | BXPC3 | 6,779 | 6,493 | 683 | 10.5 % |
| 28 | MCF7 | HAP1 | 6,779 | 6,620 | 732 | 11.1 % |
| 29 | MCF7 | HT29 | 6,779 | 6,543 | 852 | 13.0 % |
| 30 | MCF7 | K562 | 6,779 | 6,611 | 923 | 14.0 % |

🔴 **Intermedios (CRUDO): 20 filas especificas**, equiespaciadas por CPM del control de L2

| #  | pert | gen | L1 | L2 | CPM ctrl L1 | CPM ctrl L2 | pasa L2 | det ctrl L2 | n ctrl L2 | UMI med L2 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `AKT1` | `LINC00376` | BXPC3 | A549 | 25.575 | 0.000 | no | 0.0000 | 1379 | 3319 |
| 2 | `ID2` | `DSG2-AS1` | BXPC3 | K562 | 5.963 | 1.103 | no | 0.0072 | 1398 | 4302 |
| 3 | `TYK2` | `SLC22A3` | A549 | K562 | 209.790 | 2.879 | no | 0.0100 | 1398 | 4302 |
| 4 | `MAPK8` | `MXD3` | MCF7 | K562 | 8.715 | 5.585 | si | 0.0215 | 1398 | 4302 |
| 5 | `VAV1` | `PPP1R11` | MCF7 | HT29 | 8.824 | 8.949 | si | 0.0450 | 2736 | 4396 |
| 6 | `TYK2` | `METTL5` | MCF7 | HT29 | 12.634 | 12.816 | si | 0.0680 | 2736 | 4396 |
| 7 | `VAV1` | `PLXDC2` | MCF7 | K562 | 450.558 | 17.162 | si | 0.0608 | 1398 | 4302 |
| 8 | `IRF7` | `E4F1` | MCF7 | HAP1 | 13.446 | 22.271 | si | 0.1277 | 3172 | 5231 |
| 9 | `IFNAR1` | `SAMD12` | MCF7 | A549 | 500.147 | 28.076 | si | 0.0885 | 1379 | 3319 |
| 10 | `IFNAR2` | `SNRPD1` | A549 | HT29 | 25.168 | 34.678 | si | 0.1564 | 2736 | 4396 |
| 11 | `IFNAR2` | `PRDM16` | A549 | MCF7 | 201.961 | 42.528 | si | 0.2596 | 2669 | 6203 |
| 12 | `TYK2` | `NOP2` | A549 | K562 | 53.370 | 52.314 | si | 0.2146 | 1398 | 4302 |
| 13 | `JAK1` | `USP18` | HT29 | A549 | 56.763 | 64.082 | si | 0.2103 | 1379 | 3319 |
| 14 | `BATF2` | `MINDY2` | MCF7 | HAP1 | 53.716 | 78.910 | si | 0.3439 | 3172 | 5231 |
| 15 | `STAT1` | `MAML3` | HT29 | BXPC3 | 646.449 | 98.317 | si | 0.3163 | 3228 | 4754 |
| 16 | `ID2` | `RASEF` | MCF7 | BXPC3 | 85.577 | 125.085 | si | 0.4315 | 3228 | 4754 |
| 17 | `JAK1` | `TAPBP` | K562 | A549 | 138.671 | 165.361 | si | 0.4170 | 1379 | 3319 |
| 18 | `ETS2` | `SCMH1` | MCF7 | HT29 | 136.120 | 232.291 | si | 0.6009 | 2736 | 4396 |
| 19 | `IFNAR2` | `SEMA3C` | MCF7 | BXPC3 | 1936.670 | 375.483 | si | 0.7159 | 3228 | 4754 |
| 20 | `VAV1` | `MT-RNR1` | MCF7 | K562 | 17783.186 | 55998.997 | si | 1.0000 | 1398 | 4302 |

### 2.2 La atribucion: predictores del control de L2, orden DECLARADO p1 → p2 → p3 → p4

D² = 1 − devianza(modelo)/devianza(solo intercepto), logistica. Las dos columnas: cada predictor SOLO, y el ACUMULADO en el orden declarado. Nulo: media [min, max] de 20 sorteos. Anclado = (obs − media nulo)/(1 − media nulo).

**CRUDO**

| predictor | SOLO D² | SOLO nulo | SOLO anclado | ACUM D² | ACUM nulo | ACUM anclado | incremento ACUM |
|---|---|---|---|---|---|---|---|
| p1 pasa el filtro de L2 | **0.0460** | 0.0183 [0.0114, 0.0269] | 0.0282 | **0.0460** | 0.0183 [0.0114, 0.0269] | 0.0282 | +0.0460 |
| p2 CPM medio (deciles) | **0.1545** | 0.1007 [0.0862, 0.1245] | 0.0599 | **0.1552** | 0.1007 [0.0863, 0.1245] | 0.0607 | +0.1093 |
| p3 tasa de deteccion (deciles) | **0.1538** | 0.0950 [0.0783, 0.1119] | 0.0650 | **0.1570** | 0.1039 [0.0885, 0.1279] | 0.0592 | +0.0017 |
| p4 potencia (log n, log UMI mediana) | **0.0031** | 0.0091 [0.0000, 0.0191] | -0.0061 | **0.1575** | 0.1117 [0.0936, 0.1324] | 0.0516 | +0.0006 |

**IGUALADO EN POTENCIA**

| predictor | SOLO D² | SOLO nulo | SOLO anclado | ACUM D² | ACUM nulo | ACUM anclado | incremento ACUM |
|---|---|---|---|---|---|---|---|
| p1 pasa el filtro de L2 | **0.0446** | 0.0196 [0.0135, 0.0311] | 0.0255 | **0.0446** | 0.0196 [0.0135, 0.0311] | 0.0255 | +0.0446 |
| p2 CPM medio (deciles) | **0.1811** | 0.1223 [0.1066, 0.1532] | 0.0670 | **0.1813** | 0.1223 [0.1068, 0.1533] | 0.0672 | +0.1367 |
| p3 tasa de deteccion (deciles) | **0.1777** | 0.1202 [0.1065, 0.1487] | 0.0654 | **0.1835** | 0.1248 [0.1093, 0.1551] | 0.0671 | +0.0022 |
| p4 potencia (log n, log UMI mediana) | **0.0056** | 0.0021 [0.0000, 0.0065] | 0.0035 | **0.1878** | 0.1269 [0.1095, 0.1579] | 0.0697 | +0.0043 |

🔴 **Intermedios del nulo (CRUDO): los 20 sorteos**

| k | semilla | permutacion (L2 → π(L2), indices en A549/BXPC3/HAP1/HT29/K562/MCF7) | # puntos fijos | A det p1 | SOLO p1 | SOLO p2 | SOLO p3 | SOLO p4 | ACUM p1..p4 |
|---|---|---|---|---|---|---|---|---|---|
| 0 | 20281100 | [0, 1, 5, 2, 4, 3] | 3 | 0.1254 | 0.0223 | 0.1162 | 0.1072 | 0.0191 | 0.1324 |
| 1 | 20281101 | [3, 5, 0, 4, 1, 2] | 0 | 0.1138 | 0.0136 | 0.0896 | 0.0952 | 0.0119 | 0.1039 |
| 2 | 20281102 | [5, 0, 1, 2, 3, 4] | 0 | 0.1106 | 0.0114 | 0.0862 | 0.0783 | 0.0175 | 0.1073 |
| 3 | 20281103 | [0, 1, 5, 4, 3, 2] | 2 | 0.1250 | 0.0208 | 0.1160 | 0.1086 | 0.0022 | 0.1186 |
| 4 | 20281104 | [4, 1, 2, 5, 3, 0] | 2 | 0.1236 | 0.0231 | 0.1103 | 0.1026 | 0.0043 | 0.1168 |
| 5 | 20281105 | [5, 0, 1, 3, 4, 2] | 2 | 0.1211 | 0.0218 | 0.1122 | 0.1057 | 0.0012 | 0.1155 |
| 6 | 20281106 | [5, 4, 1, 0, 3, 2] | 0 | 0.1135 | 0.0128 | 0.0908 | 0.0868 | 0.0060 | 0.1028 |
| 7 | 20281107 | [5, 2, 4, 3, 0, 1] | 1 | 0.1174 | 0.0154 | 0.0929 | 0.0948 | 0.0191 | 0.1082 |
| 8 | 20281108 | [2, 0, 3, 4, 1, 5] | 1 | 0.1299 | 0.0203 | 0.0988 | 0.0946 | 0.0104 | 0.1145 |
| 9 | 20281109 | [4, 1, 5, 3, 0, 2] | 2 | 0.1245 | 0.0203 | 0.1094 | 0.1070 | 0.0094 | 0.1182 |
| 10 | 20281110 | [5, 1, 3, 0, 2, 4] | 1 | 0.1187 | 0.0180 | 0.1026 | 0.0955 | 0.0053 | 0.1095 |
| 11 | 20281111 | [1, 4, 2, 5, 0, 3] | 1 | 0.1152 | 0.0174 | 0.0920 | 0.0886 | 0.0010 | 0.0969 |
| 12 | 20281112 | [0, 1, 2, 3, 5, 4] | 4 | 0.1258 | 0.0269 | 0.1245 | 0.1119 | 0.0039 | 0.1310 |
| 13 | 20281113 | [4, 3, 5, 1, 2, 0] | 0 | 0.1159 | 0.0168 | 0.0985 | 0.0856 | 0.0190 | 0.1183 |
| 14 | 20281114 | [4, 0, 3, 2, 1, 5] | 1 | 0.1285 | 0.0177 | 0.0924 | 0.0900 | 0.0080 | 0.1049 |
| 15 | 20281115 | [2, 0, 4, 3, 5, 1] | 1 | 0.1209 | 0.0209 | 0.1044 | 0.0933 | 0.0161 | 0.1196 |
| 16 | 20281116 | [1, 2, 5, 3, 0, 4] | 1 | 0.1127 | 0.0131 | 0.0914 | 0.0875 | 0.0000 | 0.0936 |
| 17 | 20281117 | [1, 4, 2, 0, 5, 3] | 1 | 0.1153 | 0.0180 | 0.0956 | 0.0834 | 0.0066 | 0.1061 |
| 18 | 20281118 | [4, 3, 1, 5, 2, 0] | 0 | 0.1149 | 0.0153 | 0.0907 | 0.0827 | 0.0082 | 0.1023 |
| 19 | 20281119 | [3, 4, 0, 1, 2, 5] | 1 | 0.1286 | 0.0190 | 0.0986 | 0.1006 | 0.0132 | 0.1138 |

Bordes de los deciles (CRUDO): p2 log10(CPM+1e-3) = [0.495, 0.987, 1.267, 1.478, 1.653, 1.826, 2.004, 2.2, 2.465]; p3 deteccion = [0.0149, 0.0457, 0.0863, 0.1352, 0.1946, 0.2698, 0.3668, 0.495, 0.6736].

### 2.3 🔴 EL RESIDUO: la especificidad de contexto que NO explica ninguna cantidad del control

| R(g,L) | # especificas | **residuo B = 1 − D² acum (crudo)** | **residuo B anclado = 1 − D² acum anclado** | residuo A = 1 − frac. que no pasa filtro L2 | residuo A anclado |
|---|---|---|---|---|---|
| CRUDO | 109,830 | **0.8425** | **0.9484** | 0.8537 | 0.9701 |
| IGUALADO EN POTENCIA | 66,868 | **0.8122** | **0.9303** | 0.8614 | 0.9720 |

### 2.4 El universo U(g)

Jiang: U(g) = ∪_L R(g,L) sobre las 6 lineas. VCC 2025: un solo contexto, asi que se da |R(g)| / |filtrables|.

**Jiang CRUDO**: |U(g)| rango **[33, 3204]**, mediana 357; filtrables en al menos una linea: 15,928 (|U(g)|/eso: [0.0021, 0.2012]).

**Jiang IGUALADO**: |U(g)| rango **[7, 2296]**, mediana 192; filtrables en al menos una linea: 15,988 (|U(g)|/eso: [0.0004, 0.1436]).

🔴 **Intermedios: las 33 perturbaciones**

| #  | pert | |U| crudo | |R| A549 | |R| BXPC3 | |R| HAP1 | |R| HT29 | |R| K562 | |R| MCF7 | en las 6 | # lineas no vacias | |U| igualado |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | `AKT1` | **204** | 0 | 3 | 0 | 3 | 1 | 197 | 0 | 4 | 57 |
| 2 | `BATF2` | **213** | 0 | 74 | 7 | 2 | 1 | 131 | 0 | 5 | 81 |
| 3 | `BRD9` | **195** | 1 | 107 | 0 | 0 | 32 | 56 | 0 | 4 | 80 |
| 4 | `CEBPG` | **357** | 4 | 33 | 18 | 10 | 2 | 294 | 0 | 6 | 153 |
| 5 | `DTX3L` | **127** | 3 | 9 | 3 | 14 | 1 | 102 | 1 | 6 | 56 |
| 6 | `ELK1` | **230** | 0 | 17 | 15 | 1 | 35 | 164 | 0 | 5 | 102 |
| 7 | `ETS1` | **629** | 0 | 385 | 2 | 5 | 8 | 236 | 0 | 5 | 336 |
| 8 | `ETS2` | **318** | 19 | 9 | 0 | 192 | 31 | 72 | 0 | 5 | 214 |
| 9 | `H1F0` | **559** | 0 | 10 | 0 | 22 | 2 | 529 | 0 | 4 | 295 |
| 10 | `HERC5` | **472** | 0 | 8 | 4 | 0 | 2 | 458 | 0 | 4 | 271 |
| 11 | `HERC6` | **101** | 0 | 31 | 4 | 1 | 1 | 67 | 0 | 5 | 42 |
| 12 | `ID2` | **1526** | 0 | 80 | 0 | 15 | 1 | 1438 | 0 | 4 | 933 |
| 13 | `IFNAR1` | **1798** | 494 | 531 | 184 | 275 | 189 | 842 | 36 | 6 | 1174 |
| 14 | `IFNAR2` | **3204** | 751 | 1164 | 177 | 615 | 281 | 1368 | 44 | 6 | 2296 |
| 15 | `IRF1` | **239** | 10 | 72 | 0 | 2 | 0 | 163 | 0 | 4 | 110 |
| 16 | `IRF7` | **490** | 0 | 5 | 10 | 0 | 4 | 474 | 0 | 4 | 259 |
| 17 | `JAK1` | **953** | 94 | 363 | 60 | 423 | 139 | 115 | 0 | 6 | 581 |
| 18 | `JUN` | **712** | 6 | 39 | 0 | 42 | 5 | 627 | 0 | 5 | 306 |
| 19 | `MAP3K14` | **267** | 40 | 50 | 0 | 2 | 1 | 175 | 0 | 5 | 135 |
| 20 | `MAPK8` | **1125** | 1 | 72 | 0 | 3 | 0 | 1053 | 0 | 4 | 689 |
| 21 | `PARP12` | **67** | 0 | 4 | 2 | 51 | 0 | 11 | 0 | 4 | 20 |
| 22 | `RAP1GAP` | **41** | 2 | 8 | 0 | 0 | 0 | 31 | 0 | 3 | 7 |
| 23 | `RNF114` | **324** | 0 | 38 | 4 | 49 | 4 | 233 | 0 | 5 | 173 |
| 24 | `SOCS1` | **501** | 0 | 32 | 7 | 3 | 0 | 460 | 0 | 4 | 192 |
| 25 | `SP100` | **168** | 1 | 65 | 0 | 13 | 1 | 93 | 0 | 5 | 68 |
| 26 | `STAT1` | **906** | 1 | 264 | 207 | 465 | 64 | 67 | 0 | 6 | 559 |
| 27 | `STAT3` | **147** | 84 | 15 | 2 | 0 | 2 | 44 | 0 | 5 | 109 |
| 28 | `STAT4` | **33** | 0 | 13 | 0 | 0 | 0 | 20 | 0 | 2 | 10 |
| 29 | `STAT5A` | **253** | 0 | 1 | 2 | 3 | 61 | 187 | 0 | 5 | 118 |
| 30 | `TRAFD1` | **1153** | 3 | 9 | 0 | 4 | 0 | 1140 | 0 | 4 | 873 |
| 31 | `TYK2` | **2659** | 1021 | 1108 | 107 | 278 | 172 | 1027 | 45 | 6 | 1921 |
| 32 | `VAV1` | **432** | 0 | 6 | 0 | 10 | 2 | 415 | 0 | 4 | 271 |
| 33 | `ZBP1` | **601** | 2 | 17 | 2 | 2 | 0 | 578 | 0 | 5 | 340 |

**VCC 2025**: filtrables 10,776. |R(g)| rango [0, 7201], mediana 888; fraccion de los filtrables en R(g): rango **[0.0000, 0.6682]**, mediana 0.0824. Perts con R vacio: 3. Genes de R(g) fuera de los filtrables (tiene que ser 0): 0. ∪_g R(g) sobre las 100 perts: 10,767 = 0.9992 de los filtrables.

🔴 **Intermedios 2025: 20 perturbaciones equiespaciadas por |R(g)|**

| #  | pert | celulas | |R(g)| | fraccion de filtrables |
|---|---|---|---|---|
| 1 | `CNN3` | 400 | 0 | 0.0000 |
| 2 | `APOE` | 400 | 3 | 0.0003 |
| 3 | `MPC2` | 400 | 4 | 0.0004 |
| 4 | `ZNF664` | 400 | 10 | 0.0009 |
| 5 | `GULP1` | 400 | 15 | 0.0014 |
| 6 | `STXBP1` | 400 | 23 | 0.0021 |
| 7 | `HDGF` | 400 | 27 | 0.0025 |
| 8 | `NFX1` | 400 | 50 | 0.0046 |
| 9 | `PLAGL2` | 400 | 303 | 0.0281 |
| 10 | `LTBR` | 400 | 776 | 0.0720 |
| 11 | `MLLT1` | 400 | 998 | 0.0926 |
| 12 | `PIAS1` | 400 | 1229 | 0.1140 |
| 13 | `BECN1` | 349 | 2017 | 0.1872 |
| 14 | `ZRSR2` | 76 | 2472 | 0.2294 |
| 15 | `SAFB2` | 400 | 3130 | 0.2905 |
| 16 | `CDC25B` | 400 | 3708 | 0.3441 |
| 17 | `BRPF1` | 366 | 4316 | 0.4005 |
| 18 | `ADNP` | 291 | 4741 | 0.4400 |
| 19 | `PGAM1` | 400 | 5449 | 0.5057 |
| 20 | `MED19` | 400 | 7201 | 0.6682 |

---

## Defectos propios

1. **El residuo de 2.3 incluye el ruido de muestreo de las DOS llamadas** (R(g,L1) y R(g,L2)), no solo biologia. En Jiang no hay replica dentro de linea para medir ese suelo, asi que el residuo es una cota superior de lo que el control no explica por razones biologicas.
2. **El nulo usa permutaciones uniformes**, que dejan puntos fijos (columna `# puntos fijos`): en esos sorteos parte de las estadisticas son las verdaderas, lo que sube la media del nulo y hace conservador el anclaje.
3. **p4 es de linea**: con 6 lineas, 2 coeficientes sobre 6 valores capturan cualquier diferencia de linea correlacionada con n y profundidad, no solo potencia.
4. **D² depende de la codificacion** (deciles, logistica con L2 debil C = 1e4). Otra codificacion da otros numeros; esta se fijo antes de medir.
5. **Las filas no son independientes** (el mismo gen entra en varios pares y perturbaciones): no se dan intervalos; la unica referencia de dispersion es el rango del nulo.
6. **La igualacion solo iguala n de control y profundidad**. El lado perturbado ya tiene 400 celulas en todas las lineas; la eficiencia de la guia es del lado perturbado y queda fuera por diseño.
7. **Jiang tiene pools de control de 1,379–3,228 celulas** y un espacio de 34,025 genes; el reto, 18,400 y 18,533. El regimen de potencia es otro: por eso de la Parte 2 no viaja ninguna fraccion.
8. **VCC 2025 es un solo contexto**: no da especificidad de contexto, solo 2.4.
9. **El Fano es de cuentas crudas** y depende de la profundidad.
10. **El control de 1.4 usa mitades** (9,200 celulas): su ruido es ~√2 mayor que el de un contexto entero; acota por arriba la fragilidad por muestreo, no la mide al tamaño real.
11. **La medida A da a p1 todo el credito determinista** porque p1 va primero en el orden declarado; la columna SOLO de 2.2 es la que permite ver cuanto toma cada predictor por su cuenta.
12. **Las cifras agregadas de 2.1–2.3 estan dominadas por MCF7 como L1** (Σ|R| = 12,867 de 24,410 en crudo): cada fila pesa igual, asi que las lineas con R(g,L) grande pesan mas. La tabla por par ordenado de 2.1 es la que lo deja ver.
13. **El nulo de 2.1/2.2 ya es alto** (p. ej. ~12 % en la medida A) porque la expresion de un gen se parece entre lineas: las estadisticas de otra linea ya predicen parte de y. El anclaje descuenta eso; la cifra cruda no.

_Generado 2026-11-05 14:30:08 UTC. Fuentes: `/datos/modelo/out/seleccion/p1_filtrables.json`, `p2_datos.json`, `p2_mecanismo.json`; codigo en `modelo/src/seleccion/`._
