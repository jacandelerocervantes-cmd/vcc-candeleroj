# El mapa de sensibilidad de la superficie de entrega

## 🔴 ENCABEZADO OBLIGATORIO

El **brazo B** y el **brazo C** miden **LA SENSIBILIDAD DEL MARCADOR** sobre **los controles de la entrega**, no la biologia de D/E/F, **que no existe**: de D/E/F solo hay 18,400 celulas de control por contexto y ninguna perturbada. La estructura inyectada viene del **brazo A** y es una **TRANSFERENCIA DE ESTRUCTURA declarada, dada como rango** — tres numeros de forma (amplitud, fraccion de genes movidos, fraccion comun), nunca valores de gen ni de perturbacion. **Ninguna cifra de aqui es una prediccion del score real.**

El **brazo A** se mide en Jiang (6 contextos, somera) y VCC 2025 (1 contexto, profunda), que **no son** la superficie de entrega.

Dos referencias externas que fijan la escala, y que no se recalculan aqui:

* el **nulo** de `ntc_nulo.md`: a `n = 400` y en estos mismos contextos, partiendo controles contra controles, el conjunto ED sale **vacio** y las metricas de DE no registran nada. Cero medido, no supuesto.
* el **azar** de `pds_cosine`: **0.5 exacto**, derivado del codigo en `pds_objetivo.md` §1.2 (`1 - ((n-1)/2)/(n-1)`), para todo `n`.

---

## 0. La construccion, declarada antes de medir

### 0.1 Las tres superficies y que mide cada brazo

| brazo | superficie | n contextos | que mide | que NO mide |
|---|---|---|---|---|
| A | Jiang IFNB + VCC 2025 | 6 + 1 | que estructura tiene el delta REAL | nada de la entrega |
| B | controles de D/E/F (y A/B/C aparte) | 3 + 3 | que registra el marcador con verdad SINTETICA de rejilla declarada | biologia; la rejilla la elegi yo |
| C | controles de D/E/F | 3 | lo mismo con la estructura MEDIDA en A, el techo y el mapa calidad->puntuacion | no transfiere valores, solo forma |

### 0.2 Semillas: familia NUEVA, comprobadas EN CODIGO y guardadas enteras

Familia: `20283000 + 100*fold + s, fold desde 1 (D=1,E=2,F=3,A=4,B=5,C=6)`. **300 semillas**, guardadas enteras en `/datos/modelo/out/sensib/semillas.json` y en el json de este reporte (`lista_plana`), no solo su conteo. Minimo 20283100, maximo 20283951.

| comprobacion | resultado |
|---|---|
| choques internos entre familias | **0** |
| dentro de los rangos vetados [[20260000, 20266009], [20270000, 20282999]] | **0** |
| contra la enumeracion de `noche/comun.py` | **0** |
| contra PASO 5a (`seleccion/semillas.py`) | **0** |
| contra el suelo de llamada (`sl_semillas.py`) | **0** |
| contra el barrido de 861 enteros 202[6-8]xxxx en `modelo/src/**/*.py` y los json de `modelo/out/` | **0** |
| **sin choques** | **si** |

### 0.3 La superficie de entrega como banco de sensibilidad

| que | valor |
|---|---|
| contextos con verdad en disco | **ninguno** (D/E/F son solo control) |
| celulas de control por contexto | 18,400 (46 `ntc_id` x 400 exactas) |
| genes | 18,533 |
| nnz del contexto leido | 104,487,981 |
| X entero | si |
| particion POR GUIA (nunca por celula) | 9,200 celulas de control (23 guias) / 23 pseudo-perturbaciones de 400 |
| guia partida | **no** |
| regimen de la llamada | 400 contra 9,200 — el brazo `n400_ref_completa` de `ntc_nulo.md`, el declarado ahi como el que se parece a la entrega |
| etiquetas | 23 nombres de gen del propio reto, sorteados con semilla `20283950` entre las 300 del panel (el marcador exige que resuelvan a gen: con etiquetas inventadas `baseline` se niega) |
| genes filtrables (CPM control > 5, menos los 23 de etiqueta) | 9,608 |
| lado real | [18400, 18533] filas x genes |
| lado pred | [9600, 18533] |

### 0.4 El dispositivo de inyeccion y la rejilla del brazo B

La inyeccion es `cand01_covar._reescala_entero` **sin modificar**: Binomial si la razon es <= 1, Poisson aditivo si es > 1; la esperanza es `x * r` en las dos ramas y la salida sigue siendo cuentas enteras. Es el mismo dispositivo con el que `pds_objetivo.md` §3.1 emitio sus oraculos.

🔴 La rejilla del brazo B esta **declarada antes de correr** y es mia — por eso es un peso inventado, y por eso existe el brazo C:

| eje | valores | que es |
|---|---|---|
| amplitud | 0.25, 0.5, 1.0, 2.0 | `|log2FC|` inyectado |
| fraccion de genes | 0.002, 0.02, 0.2 | de los filtrables, por pseudo-perturbacion |
| fraccion comun | 0.0, 0.5 | del movimiento, comun a las 23 (mismo gen, mismo signo) |
| celdas | 24 (producto) | centro declarado: celda 9 = {'k': 9, 'amplitud': 0.5, 'frac_genes': 0.02, 'frac_comun': 0.5} |

Y el barrido de calidad del brazo C: `f` en {0.0, 0.25, 0.5, 0.75, 1.0}, en dos variantes declaradas distintas (**magnitud**: direccion correcta y magnitud encogida; **identidad**: `floor(f*n)` genes correctos y el resto en genes equivocados con la misma amplitud). En las dos el eje comun se da por acertado, para que `f` module solo la componente especifica.

---

## 1. BRAZO A — la biologia: que estructura tiene el delta real

🔴 Jiang y VCC 2025. **Ninguna cifra de esta seccion es de la entrega.**

_(brazo A sin correr: no hay `a_de.json`)_

---

## 2. BRAZO B — que registra el marcador a n=400

🔴 Verdad SINTETICA, rejilla mia. Cada celda tiene su propio lado real, su propio `baseline` y su propio `prep-real-bundle` (el bundle no se puede compartir entre celdas porque el lado real cambia), y dos brazos de prediccion contra el mismo bundle: **suelo** (sin estructura) y **oraculo** (la misma estructura, celulas distintas, ruido nuevo).

### 2.1 Las celdas corridas

| celda | ctx | amplitud | frac genes | frac comun | genes movidos | comunes | especificos | config_hash de competencia | duracion |
|---|---|---|---|---|---|---|---|---|---|
| `D_B09` | D | 0.5 | 0.02 | 0.5 | 192 | 96 | 96 | si | 17.0 min (1022 s) |

### 2.2 🔴 La superficie: las seis metricas, `u` y `s`, celda por celda y brazo por brazo

`u` es la fila `mean` de `agg_results.csv`; `s` es `from_replicate` de `scored.csv` contra el bundle de esa misma celda. El **recorrido** es `oraculo − suelo`: lo que el marcador puede registrar para esa estructura.

**`D_B09`** — ctx D, amplitud 0.5, frac genes 0.02, frac comun 0.5 (192 genes movidos)

| cantidad | pds | mse | lfc | fid | reach | jac | avg |
|---|---|---|---|---|---|---|---|
| `u` suelo | 0.4545 | 1.0929 | 0.9978 | 0.0000 | 0.0080 | 0.0000 |  |
| `u` oraculo | 1.0000 | 0.1655 | 0.2400 | 0.9350 | 0.9950 | 0.6992 |  |
| `s` suelo | -0.1058 | +0.0000 | -1.3968 | -0.8549 | -2.3808 | -1.3539 | -1.0154 |
| `s` oraculo | +1.0171 | +0.7690 | +1.1414 | +1.0159 | +0.9968 | +1.4930 | +1.0722 |
| **recorrido `u`** | +0.5455 | -0.9275 | -0.7578 | +0.9350 | +0.9870 | +0.6992 |  |
| **recorrido `s`** | +1.1229 | +0.7690 | +2.5382 | +1.8708 | +3.3776 | +2.8469 | +2.0876 |

### 2.3 Como se mueve cada metrica al mover UN eje

Del subconjunto en cruz: el centro y, desde el, un eje a la vez. Cada fila es el `u` del ORACULO (lo mejor alcanzable para esa estructura).

| celda | amplitud | frac genes | frac comun | pds | mse | lfc | fid | reach | jac | s avg |
|---|---|---|---|---|---|---|---|---|---|---|
| `D_B09` | 0.5 | 0.02 | 0.5 | 1.0000 | 0.1655 | 0.2400 | 0.9350 | 0.9950 | 0.6992 | +1.0722 |

La misma tabla para el SUELO (sin estructura): es el control de que lo que se mueve arriba es la estructura y no la celda.

| celda | amplitud | frac genes | frac comun | pds | mse | lfc | fid | reach | jac | s avg |
|---|---|---|---|---|---|---|---|---|---|---|
| `D_B09` | 0.5 | 0.02 | 0.5 | 0.4545 | 1.0929 | 0.9978 | 0.0000 | 0.0080 | 0.0000 | -1.0154 |

---

## 3. BRAZO C — la union, de donde saldrian los pesos

_(sin celdas del brazo C en disco)_

---

## 4. Defectos de la construccion, como datos sobre la medicion

1. **El banco tiene 23 pseudo-perturbaciones, la entrega tiene 300.** Con 18,400 celulas de control por contexto, repartirlas en pool de control y pseudo-perturbadas es un trueque: 23/23 da el regimen 400 contra 9,200 y 23 competidores para `pds_cosine` (granularidad 1/22 = 0.0455), no 300 (1/299).
2. **El lado pred se remuestrea del MISMO pool de control que el lado real usa de referencia.** Es el contrato del arnes (`predict_floor`), pero significa que las celulas del pred son filas del control del real, con reemplazo.
3. **La rejilla del brazo B la elegi yo.** Es el peso inventado que el encargo nombra; el brazo C existe para sustituirla por estructura medida, y las dos cosas no se mezclan en ninguna tabla.
4. **La transferencia del brazo C es de FORMA, y la forma se mide en otras superficies.** Amplitud, fraccion de genes y fraccion comun vienen de Jiang y 2025; que esos tres numeros sean los de D/E/F no se puede comprobar.
5. **La amplitud inyectada es la MISMA en todos los genes movidos.** El delta real tiene una distribucion de amplitudes (seccion 1.1 la da): el banco inyecta su mediana, no su forma.
6. **El eje comun del brazo C se da por acertado** en las dos variantes de `f`. Si un modelo no acierta el eje comun, su `s` es menor que el de este mapa.
7. **El gen propio de cada etiqueta no se mueve** y los 23 genes-etiqueta salen del sorteo de genes movidos: la exclusion del panel los quita del espacio. No se modela eficiencia de knockdown.
8. **Un `s` de aqui no es un `s` de la entrega.** El bundle es el de esta superficie sintetica, con sus 23 etiquetas; el `b` y el `r` de la entrega real no existen porque no hay verdad.
9. **El ruido de la inyeccion es Binomial/Poisson sobre una celula de control.** `sobredispersion.md` ya midio que ese dispositivo no reproduce la sobredispersion del delta real; lo que el banco mide es sensibilidad a un desplazamiento de medias, no a un cambio de dispersion.

---

## 5. Las entradas

| que | de donde |
|---|---|
| controles de la entrega | `/datos/reto/final/context_{D,E,F}.h5ad`, `/datos/reto/val/context_{A,B,C}.h5ad` (solo lectura) |
| etiquetas | `pert_counts.csv` de la particion de cada contexto |
| Jiang | `comun.arma_linea` (= `suelo_arnes.build_fold_adatas`), `seed_real = 20260001 + fold` |
| VCC 2025 | `v25.selecciona` (20280201 / 20280202) |
| deltas de Jiang para A2b | `/datos/modelo/out/magnitud_centro/deltas_*.npz` |
| reparto entre contextos | `modelo/out/magnitud_centro/{calc,ctrl}.json` |
| |ED| de referencia | `/datos/modelo/out/seleccion/p2_datos.json` |
| codigo | `modelo/src/sensib/` |
| celdas crudas | `/datos/modelo/out/sensib/celdas/*.json` |

_Generado 2026-11-05 18:43:13 UTC. Celdas en disco: 1._
