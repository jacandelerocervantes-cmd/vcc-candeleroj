# Bucle del suelo sobre 6 pliegues -- condicion IFNB

Panel comun (perturbaciones con >=400 celulas en las SEIS lineas de IFNB): **33** de 61 vistas en total.

## Descartes por linea (perturbaciones con <400 celulas EN ESA linea)

| linea | n descartadas | perturbaciones |
|---|---|---|
| A549 | 17 | ADAR, CEBPB, CRKL, HES4, ID1, IFI16, IRF3, IRF9, NFKB1, POU2F1, RAPGEF1, SMARCA5, SP110, TRIM21, TRIM22, USP18, ZNFX1 |
| HAP1 | 6 | ADAR, CRKL, IFI16, NFKB1, POU2F1, STAT6 |
| HT29 | 6 | ADAR, ID1, IFI16, NFKB1, POU2F1, STAT6 |
| K562 | 22 | CEBPB, DRAP1, ID1, ID3, IFI16, IRF3, IRF9, MAPK14, MYC, NFE2L3, NFKB1, POU2F1, RAPGEF1, SMARCA5, SP110, STAT2, STAT6, TRIM21, TRIM22, UBA7, UBE2L6, ZNFX1 |
| BXPC3 | 8 | FOS, ID1, IFI16, IRF3, MEF2A, NFKB1, STAT6, TRIM21 |
| MCF7 | 4 | IFI16, IRF9, NFKB1, POU2F1 |

## Panel comun -- 33 perturbaciones, celulas por linea (intermedio)

| perturbacion | A549 | BXPC3 | HAP1 | HT29 | K562 | MCF7 |
|---|---|---|---|---|---|---|
| AKT1 | 606 | 2499 | 1014 | 1165 | 481 | 682 |
| BATF2 | 617 | 712 | 1326 | 1002 | 483 | 1203 |
| BRD9 | 699 | 1676 | 975 | 996 | 692 | 796 |
| CEBPG | 701 | 1814 | 1485 | 1118 | 443 | 827 |
| DTX3L | 544 | 2255 | 1321 | 879 | 525 | 1786 |
| ELK1 | 856 | 2174 | 641 | 1015 | 695 | 2212 |
| ETS1 | 622 | 984 | 1044 | 1083 | 485 | 1045 |
| ETS2 | 576 | 2518 | 1587 | 968 | 700 | 1368 |
| H1F0 | 671 | 1170 | 1604 | 1281 | 629 | 1128 |
| HERC5 | 468 | 1066 | 995 | 1040 | 462 | 971 |
| HERC6 | 440 | 1640 | 1169 | 826 | 408 | 1281 |
| ID2 | 576 | 2031 | 1394 | 796 | 531 | 716 |
| IFNAR1 | 634 | 1104 | 1152 | 1148 | 460 | 644 |
| IFNAR2 | 723 | 2901 | 1261 | 1117 | 439 | 450 |
| IRF1 | 950 | 1447 | 1698 | 1351 | 568 | 1803 |
| IRF7 | 436 | 3860 | 1328 | 1086 | 528 | 557 |
| JAK1 | 693 | 1045 | 1471 | 1692 | 587 | 1289 |
| JUN | 525 | 1709 | 1604 | 787 | 503 | 1181 |
| MAP3K14 | 407 | 1734 | 982 | 980 | 429 | 1041 |
| MAPK8 | 581 | 1029 | 1385 | 999 | 509 | 1119 |
| PARP12 | 682 | 1178 | 1481 | 1062 | 647 | 1376 |
| RAP1GAP | 675 | 2832 | 1733 | 1332 | 582 | 1287 |
| RNF114 | 688 | 1007 | 1738 | 1439 | 602 | 984 |
| SOCS1 | 597 | 1907 | 1526 | 1188 | 583 | 863 |
| SP100 | 496 | 2057 | 1234 | 1534 | 579 | 724 |
| STAT1 | 630 | 1702 | 2025 | 1317 | 565 | 1277 |
| STAT3 | 505 | 1354 | 1165 | 961 | 478 | 1141 |
| STAT4 | 510 | 1520 | 1369 | 855 | 466 | 855 |
| STAT5A | 842 | 3004 | 1793 | 1278 | 869 | 1195 |
| TRAFD1 | 478 | 1231 | 796 | 1230 | 620 | 984 |
| TYK2 | 737 | 2687 | 736 | 1108 | 532 | 822 |
| VAV1 | 784 | 2172 | 2049 | 1096 | 595 | 1137 |
| ZBP1 | 509 | 1297 | 1218 | 1236 | 482 | 1598 |

## La tabla -- 6 pliegues x metricas (el entregable)

| linea excluida | n control cells | n control guides | pds_cosine | expr_mse_unbiased_capped_norm | de_wilcoxon_lfc_nmae | de_wilcoxon_direction_fidelity_yield_raw | de_wilcoxon_direction_reach_raw | de_wilcoxon_sig_jaccard | avg_score |
|---|---|---|---|---|---|---|---|---|---|
| A549 | 1379 | 14 | 0.474432 | 2.562172 | 1.036018 | 0.000000 | 0.144540 | 0.454545 | 0.000000 |
| BXPC3 | 3228 | 14 | 0.517992 | 1.332759 | 0.996965 | 0.000026 | 0.142889 | 0.000000 | 0.000000 |
| HAP1 | 3172 | 14 | 0.512311 | 1.868728 | 1.000821 | 0.000000 | 0.192408 | 0.424242 | 0.000000 |
| HT29 | 2736 | 14 | 0.490530 | 1.931565 | 1.044314 | 0.009259 | 0.216033 | 0.181818 | 0.000000 |
| K562 | 1398 | 14 | 0.495265 | 2.913227 | 1.028846 | 0.000000 | 0.200448 | 0.242424 | 0.000000 |
| MCF7 | 2669 | 14 | 0.512311 | 1.774874 | 1.056808 | 0.000000 | 0.013563 | 0.000000 | 0.000000 |

## Intermedios de tiempo y limpieza, por pliegue

| linea | fold_duration_s | run_duration_s | free_gib_before | real_shape | pred_shape | seed_real | seed_pred |
|---|---|---|---|---|---|---|---|
| A549 | 71.49 | 61.23 | 36.13 | [14579, 34025] | [13600, 34025] | 20260001 | 20261001 |
| BXPC3 | 107.92 | 81.3 | 36.13 | [16428, 34025] | [13600, 34025] | 20260002 | 20261002 |
| HAP1 | 109.01 | 84.49 | 36.13 | [16372, 34025] | [13600, 34025] | 20260003 | 20261003 |
| HT29 | 92.63 | 74.46 | 36.13 | [15936, 34025] | [13600, 34025] | 20260004 | 20261004 |
| K562 | 69.08 | 63.97 | 36.13 | [14598, 34025] | [13600, 34025] | 20260005 | 20261005 |
| MCF7 | 105.82 | 81.93 | 36.13 | [15869, 34025] | [13600, 34025] | 20260006 | 20261006 |

## Correlacion con el grosor del control (n=6, solo descriptiva)

| metrica | pearson r vs n_control_cells |
|---|---|
| pds_cosine | 0.792 |
| expr_mse_unbiased_capped_norm | -0.932 |
| de_wilcoxon_lfc_nmae | -0.423 |
| de_wilcoxon_direction_fidelity_yield_raw | 0.180 |
| de_wilcoxon_direction_reach_raw | -0.106 |
| de_wilcoxon_sig_jaccard | -0.432 |

`expr_mse_unbiased_capped_norm` correlaciona fuerte y negativamente (r=-0.93): con un control mas fino (menos celulas), el estimador insesgado de MSE normalizado sale mas alto (peor) -- consistente con mas ruido de muestreo en el perfil de control. `pds_cosine` correlaciona positivo (r=0.79): mas celulas de control, mejor separabilidad. Las columnas de DE (`lfc_nmae`, `sig_jaccard`) correlacionan negativo, mas debil. Con solo 6 puntos esto es descriptivo, no una prueba estadistica.

## avg_score = 0.0 en los seis pliegues -- no es un error

El suelo se puntua CONTRA SI MISMO (`--user-agg` y `--baseline-agg` son el mismo `agg_results.csv`), por el camino diagnostico de `score` que reporta `from_baseline` sin exigir ancla (`--allow-config-mismatch`, sin `--anchor` ni `--real-bundle`). En la escala vcc2026, 0 es por definicion el suelo de respuesta media -- asi que un suelo puntuado contra si mismo cae en 0.0 por construccion, en los seis pliegues. Lo informativo no es esta columna (siempre 0 aqui): son las seis columnas de metricas crudas de arriba, que SI varian por pliegue y que quedan fijadas como el `agg_results.csv` de referencia (`--baseline-agg`) para puntuar futuros candidatos contra el suelo.

## Decisiones de construccion, no triviales

- `cell_eval2.io.validate_pair` exige que PRED y REAL compartan exactamente el mismo conjunto de etiquetas de perturbacion, y que la etiqueta de control este presente en REAL. El encargo solo describe la VERDAD como "las celulas perturbadas" -- se le agrego el bloque de control de la linea (TODO el pool real de esa linea, sin submuestrear, sin modificar) para cumplir ese contrato de la herramienta. El bloque de control de SUELO se construye con la MISMA regla que cualquier otra etiqueta: 400 celulas remuestreadas con reemplazo del propio pool de control (semilla distinta a la de la VERDAD), no una copia de las celulas reales -- asi la funcion de prediccion no necesita un caso especial para el control.
- El arnes (`modelo/src/suelo_arnes.py`) nunca pasa el nombre ni la identidad de la linea a la funcion de prediccion: `predict_fn(control_pool_X, labels, rng)` solo recibe la matriz de control y la lista de etiquetas a predecir. `suelo_floor.predict_floor` (modelo/src/suelo_floor.py) es una implementacion de esa firma entre las que acepta `suelo_arnes.build_fold_adatas` -- un candidato real se conecta sin tocar el bucle (`modelo/src/bucle_suelo.py`).
- `--preset vcc2026 --profile vcc2026 --pert-col target_gene --control NT` en `cell-eval2 run`, usando `~/celleval2_venv` (0.20.0, instalado con `pdex` para DE en CPU -- sin GPU en esta VM, confirmado por `de.backend='auto'` cayendo a `pdex` en el log de cada pliegue).
- Lectura de `/X` sin densificar nunca la matriz completa: para cada linea se agrupan sus filas (dispersas en ~32 bloques contiguos dentro del h5ad de 328,542x34,025) y se leen esos bloques con `h5py`; el submuestreo (VERDAD sin reemplazo, SUELO con reemplazo) se hace en memoria sobre ese bloque ya leido, no celda por celda contra disco.

## Criterio de adopcion (recordatorio, no aplica aqui)

Esta corrida establece el suelo, no lo compara contra nada: no hay ningun candidato todavia. El criterio declarado para cuando exista uno sigue siendo el mismo -- un candidato sustituye al suelo solo si gana en los SEIS pliegues; si gana en 5, se reporta en cual perdio; si ninguno domina en los seis, no hay metodo.

---

# ADDENDUM 2026-10-29 04:51 UTC -- cierre de sesiones duplicadas/completadas en esta VM

Antes de continuar, se verifico con `tmux list-panes -a` que `cc:0` (ventana `claude`, PID 1378, 11h11m activa) y `cc:1` (ventana `estado_actual`, PID 31794, 2h41m activa) seguian siendo exactamente esas -- mismo indice, mismo nombre, mismo PID -- y se cerraron con `tmux kill-window -t cc:0` / `-t cc:1` (nunca `kill` directo al PID). Confirmado despues con `tmux ls` y `tmux list-panes -a`: solo queda esta ventana (`cc:2`, `bucle_suelo`). `cc:0` tenia su tarea de reestructuracion de `benchmark/`/`out/`/S3 ya completada y verificada (sin pendientes); `cc:1` era la tarea duplicada del suelo IFNB, ya confirmada y detenida por su propia sesion antes del cierre.

**Sobre el "escritor no identificado" de `s3://vc-challenge-data/out/`:** `syncout.service` (systemd, activo, no tocado) solo sincroniza `~/vcc/` hacia el prefijo `s3://vc-challenge-data/vcc/` y `/datos/_scratch_techo_metricas/` -- **nunca** hacia el prefijo raiz `s3://vc-challenge-data/out/`. Sin embargo los objetos en `out/` (raiz) tienen nombres identicos a los que `cc:0` escribio localmente en `~/vcc/out/` (`mapa_s3.txt`, `mapa.md`, `reestructura.md`, `estado_s3.txt`, `estado_local.txt`, `estado_actual_tabla.md`, `estado_local_after.txt`, `estado_actual.md`), con timestamps 2-3 segundos despues de cada escritura local -- demasiado rapido para el ciclo de 120s de `syncout.service`, y a un prefijo que ese servicio ni cubre. Esto apunta a que `cc:0` empujaba esos archivos a S3 por su cuenta (parte de su propia tarea), no un demonio aparte. El ultimo objeto nuevo en `out/` es de las 02:11:06 UTC, 2.5h antes de que esta sesion cerrara la ventana -- `cc:0` ya estaba inactiva (esperando respuesta del usuario) desde entonces, asi que "deja de escribir tras cerrarla" no es una prueba disponible aqui; la correlacion de nombres y timestamps es la evidencia real, y es mas fuerte que esa prueba hubiera sido.

No se toco ningun otro proceso ni `syncout.service`.

