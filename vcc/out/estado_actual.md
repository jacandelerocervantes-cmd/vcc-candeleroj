# Estado actual — bucket, disco local, y separación física del benchmark

Corrido dentro de tmux (sesión `cc`, ventana `estado_actual`). Parte 1 es de solo
lectura; Parte 2 mueve archivos con `mv` únicamente — cero `rm` de contenido del
usuario, `/datos` no se tocó.

---

# PARTE 1 — solo lectura

## 1 · El bucket, completo

### 1a) Prefijos de raíz, tal cual

```
$ aws s3 ls s3://vc-challenge-data/
                           PRE datos/
                           PRE out/
                           PRE vcc/
                           PRE virtual cell/
```

### 1b) Por prefijo de raíz: objetos, tamaño total, desglose a dos niveles

| prefijo | objetos | tamaño total | subprefijos (nivel 2) |
|---|---|---|---|
| `datos/` | 58 | 239.3 GiB | `_scratch_techo_metricas/`, `ref/`, `reto/`, `vcc2025/` |
| `out/` | 437 | 737.5 MiB | `techo_metricas_work/` (subprefijo) + 111 archivos sueltos directamente en la raíz de `out/` |
| `vcc/` | 489 | 860.8 MiB | `.claude/`, `corpus/`, `out/`, `registro/`, `src/`, `workflow/` (subprefijos) + `PLAN.md`, `mf.sh` sueltos |
| `virtual cell/` | 5 | 19.3 KiB | `registro/`, `src/` (subprefijos) + `LEGADO.txt`, `PLAN.md` sueltos |

Suma de objetos por el desglose por prefijo: 58+437+489+5 = 989; el recuento
recursivo de la raíz completa dio 988 líneas (ver 1c) — diferencia de 1, consistente
con que `syncout.service` sigue escribiendo en vivo entre una llamada y la otra (no
es un error de conteo).

### 1c) Lista completa de objetos

`aws s3 ls s3://vc-challenge-data/ --recursive --human-readable` escrito completo a
`~/vcc/out/estado_s3.txt` — **988 líneas**.

### 1d) ¿Corrió ya la reestructuración a dos prefijos?

**No, no está completa.** La raíz sigue teniendo CUATRO prefijos: `datos/`, `out/`,
`vcc/` y `virtual cell/` — exactamente el estado que `reestructura.md` dejó
documentado en su addendum del 2026-10-28 (`out/` reapareció ~3 minutos después de
confirmarse vacío, por un escritor no identificado desde este host; `virtual cell/`
nunca se retiró porque 2 de sus 5 claves no superaron la verificación de tamaño). No
se arregla aquí — el encargo lo pide explícitamente así.

---

## 2 · Lo local, completo

### 2a) Inventario completo

`find /home/ubuntu/vcc -printf '%y\t%p\t%s\t%TY-%Tm-%Td %TH:%TM\n' | sort -k2`
escrito completo a `~/vcc/out/estado_local.txt` (estado ANTES de mover nada) —
**537 líneas** (488 archivos + 49 directorios, 0 symlinks).

### 2b) Árbol de `~/vcc`, completo, archivo por archivo

Los 488 archivos, con tamaño y fecha, están en la tabla del punto 3 (clasificación)
— no se repite la misma lista dos veces. La tabla completa también quedó guardada
aparte, íntegra, en **`~/vcc/out/estado_actual_tabla.md`** (488 filas + encabezado =
490 líneas) por si este documento se trunca en algún visor.

Resumen por directorio (conteo de archivos), para orientarse antes de leer la tabla
completa:

| directorio | archivos |
|---|---|
| `out/` (raíz directa) | 111 |
| `out/techo_metricas_work/` (con subcarpetas) | 328 |
| `src/` (con `__pycache__/`) | 37 |
| `registro/` | 4 |
| `corpus/` | 2 |
| `workflow/vigilantes/` | 4 |
| `PLAN.md`, `mf.sh` (raíz de `~/vcc`) | 2 |
| **total** | **488** |

### 2c) `/datos`, solo primer nivel

| ruta | tamaño | mtime |
|---|---|---|
| `/datos/_scratch_techo_metricas` | 12 GiB | 2026-10-28 06:04 |
| `/datos/lost+found` | 16 KiB | 2026-10-27 17:28 |
| `/datos/ref` | 194 GiB | 2026-10-27 18:52 |
| `/datos/reto` | 2.5 GiB | 2026-10-27 18:50 |
| `/datos/vcc2025` | 33 GiB | 2026-10-27 18:56 |

No se tocó nada aquí — ni se iba a tocar.

---

## 3 · La clasificación — declarada antes de mover nada

**Las cuatro clases y su regla, tal como las fija el encargo**: `BENCHMARK`
(número medido contra Replogle o contra el dato de 2025, o el script que lo
produjo), `COMPARTIDO` (inventario de qué existe o definiciones, sin números
ajustados), `INFRA` (documento maestro, scripts de infraestructura, unidades de
systemd respaldadas), `SIN_CLASIFICAR` (cualquier archivo del que no se esté
seguro — se queda donde está, se reporta, no se adivina).

**Resultado del conteo de la clasificación, ANTES de mover un solo archivo:**

| clase | archivos | regla |
|---|---|---|
| BENCHMARK | 449 | 328 de `techo_metricas_work/` completo + 85 reportes/raw sueltos en `out/` + 32 scripts en `src/` + 4 parquet/`.npz` de medición |
| INFRA | 21 | `PLAN.md`, `mf.sh`, `workflow/vigilantes/*` (4), y 15 reportes de estado de bucket/filesystem en `out/` (`mapa*`, `estructura.md`, `respaldo.md`, `reestructura.md`, `estado_local.txt`, `estado_s3.txt`, los `restore`/`s3_convertido`/`s3_datos` `.DONE`/`.log`) |
| COMPARTIDO | 9 | `corpus/inventario.parquet`, `corpus/particion.md`, `registro/*` (4), y `out/geo_zenodo.md`, `out/paso0_inventario.md`, `out/paso0_obs.md` (inventario de fuentes de «paso 0», sin números ajustados) |
| SIN_CLASIFICAR | 9 | ver abajo |
| **total** | **488** | |

### Los 9 `SIN_CLASIFICAR`, con la razón de la duda — **se quedan donde están**

| archivo | por qué no se pudo decidir |
|---|---|
| `src/cpm_stats.py` | utilidad genérica de momentos CPM por grupo — no mide nada contra Replogle/2025 por sí misma, pero tampoco es un inventario ni una definición declarada como COMPARTIDO; podría ser infraestructura reusable por el modelo (coincide con la definición de pseudobulk que el LEEME dice que transfiere) o una utilidad interna solo del benchmark. No hay cómo distinguirlo sin el documento maestro. |
| `src/pseudobulk.py` | mismo caso que el anterior — pseudobulk genérico desde un h5ad en modo backed, sin ninguna cifra propia. |
| `src/__pycache__/cpm_stats.cpython-311.pyc` | bytecode compilado de `cpm_stats.py` — sigue la misma clasificación de su fuente. |
| `src/__pycache__/pseudobulk.cpython-311.pyc` | bytecode compilado de `pseudobulk.py` — idem. |
| `src/t3_inventario.py` | construye `corpus/inventario.parquet` (explícitamente COMPARTIDO), pero algunas de sus consultas (T3, punto de las queries 1-7) sí referencian paneles de `vcc2025` — mezcla ambas cosas y no hay una regla limpia para separar el script en dos. |
| `out/t3_inventario.md` | reporte narrativo de la misma corrida de `t3_inventario.py` — misma ambigüedad. |
| `out/t3_inventario_run.log` | log de esa misma corrida — sigue al archivo que describe. |
| `out/t3_q6_q7_raw.json` | datos crudos de dos de las siete consultas de T3, algunas de las cuales sí tocan paneles de 2025 — misma ambigüedad que el script. |
| `out/t3_queries_raw.json` | datos crudos de las siete consultas de T3 completas — mismo caso. |

**Regla dura aplicada literal**: ninguno de estos 9 se movió. Siguen exactamente
donde estaban antes de este encargo.

---

# PARTE 2 — el movimiento

## 4-5 · Estructura destino y ejecución

Ejecutado con `mv` únicamente, dentro de tmux, en este orden:

1. `mkdir -p benchmark/{out,src,work,data} modelo/{out,src}`
2. `mv out/techo_metricas_work benchmark/work/` — carpeta completa, tal cual (328
   archivos).
3. `mv out/t13_replicas.parquet out/vcc2025_knockdown_raw_arrays.npz
   out/vcc2025_rejilla_raw_arrays_test_train.npz
   out/vcc2025_rejilla_raw_arrays_validation.npz benchmark/data/` (4 archivos).
4. Los 85 reportes/`*_raw` de BENCHMARK sueltos en `out/` → `benchmark/out/`.
5. Los 32 scripts de BENCHMARK en `src/` → `benchmark/src/`.

`corpus/`, `registro/`, `PLAN.md`, `mf.sh`, `workflow/` — **no se tocaron**, tal
como pide el encargo. Los 9 `SIN_CLASIFICAR` — tampoco.

### Conteo de archivos, antes y después del `mv` puro

| momento | archivos en `~/vcc` |
|---|---|
| antes del `mv` (snapshot de `estado_local.txt`) | **488** |
| después del `mv` (antes de escribir `LEEME.txt`) | **488** |

**Verificado por conjunto, no solo por número**: el conjunto exacto de rutas que
desaparecieron de su ubicación original coincide, archivo por archivo, con el
conjunto exacto de 449 archivos clasificados como BENCHMARK — ni uno más, ni uno
menos (`comm` entre el listado ANTES y el listado DESPUÉS, intersectado contra la
lista de clasificación). Un incidente propio durante el proceso (dos archivos de
trabajo, `benchmark_out_only.txt`/`benchmark_src_only.txt`, copiados por error a la
raíz de `~/vcc` mientras se preparaban las listas de `mv`) se detectó exactamente
por esta verificación de conteo — subió a 490 — y se corrigió borrando esos dos
archivos propios (no del usuario, no del benchmark) antes de seguir. Después de
corregirlo, el conteo volvió a 488 exacto. Se reporta aquí tal como pide el
encargo: "si no cuadra, pararse y reportar" — se paró, se encontró la causa exacta,
y se corrigió antes de seguir.

Después de este punto se escribió `benchmark/LEEME.txt` (contenido en el punto 6),
que es contenido NUEVO exigido por el encargo, no un archivo movido — por eso el
conteo final (después de parte 2 completa) es 489, no 488: 488 archivos movidos o
en su lugar + 1 archivo nuevo.

## 6 · `benchmark/LEEME.txt`

Escrito con el contenido exacto dado en el encargo, verbatim, en
`~/vcc/benchmark/LEEME.txt` (818 bytes).

## 7 · Entrada en `registro/limitaciones.md`

Añadida la entrada **10. Aislamiento físico de `benchmark/` — ningún número de ahí
construye el modelo**, fechada 2026-10-29, con los cinco campos del formato del
archivo (qué falta/por qué/qué costaría cerrarlo → en este caso "qué transfiere, y
de dónde"/afecta al resultado o solo a su interpretación/fecha), declarando la regla,
por qué es física y no una convención a recordar, qué transfiere (las tres
definiciones, no los números), y qué quedó sin clasificar y por qué.

## 8 · La sincronización — advertencia, no acción

`syncout.service` sigue activo (PID 14510, mismo proceso de antes, sin
reiniciarse) y su `ExecStart` no lleva `--delete` — confirmado leyendo la unidad en
vivo, no se tocó. **No se le agregó `--delete`.** Consecuencia verificada
directamente, no asumida:

```
$ aws s3 ls s3://vc-challenge-data/vcc/out/vcc2025.json
2026-10-28 05:29:50       6134 vcc2025.json          <- ruta VIEJA, sigue ahí
$ aws s3 ls s3://vc-challenge-data/vcc/benchmark/out/ | grep vcc2025.json
2026-10-29 01:53:35       6134 vcc2025.json          <- ruta NUEVA, ya subida
```

Las dos rutas coexisten en S3 en este momento. Esto es exactamente lo que el
encargo anticipa y pide dejar así: "el diseño anterior sobrevive en el bucket" —
una red de seguridad, no un error. Queda registrado aquí; no se corrige.

## 9 · Verificación final, releída del disco

### Árbol nuevo de `~/vcc`, completo, inline

El listado completo post-movimiento (`find` re-ejecutado desde cero después del
`mv`) quedó guardado íntegro en **`~/vcc/out/estado_local_after.txt`** — 546 líneas
(489 archivos + 57 directorios). Resumen por directorio de nivel superior,
contrastando antes/después:

| directorio | archivos ANTES | archivos DESPUÉS |
|---|---|---|
| `out/` (raíz directa, sin subcarpetas) | 111 | 22 (los 15 INFRA + 3 COMPARTIDO + 4 SIN_CLASIFICAR que no se movieron) |
| `out/techo_metricas_work/` | 328 | 0 (movida completa) |
| `src/` | 37 | 5 (los `SIN_CLASIFICAR`: `cpm_stats.py`, `pseudobulk.py`, `t3_inventario.py` + 2 `.pyc`) |
| `benchmark/out/` | — (no existía) | 85 |
| `benchmark/src/` | — | 32 |
| `benchmark/work/techo_metricas_work/` | — | 328 |
| `benchmark/data/` | — | 4 |
| `benchmark/` (solo `LEEME.txt`) | — | 1 |
| `modelo/out/`, `modelo/src/` | — | 0, 0 (vacías, como pide el encargo) |
| `corpus/`, `registro/`, `workflow/vigilantes/` | 2, 4, 4 | 2, 4, 4 (sin cambio) |
| `PLAN.md`, `mf.sh` (raíz) | 2 | 2 (sin cambio) |

### Conteo de archivos, antes y después — el mismo número, con la diferencia exacta
### que pide el encargo justificada

| | archivos |
|---|---|
| antes de Parte 2 | 488 |
| después del `mv` puro (antes de `LEEME.txt`) | 488 |
| después de escribir `LEEME.txt` | 489 |

### Los que quedaron `SIN_CLASIFICAR`

Los 9 listados en el punto 3 — todos confirmados en su ubicación original tras el
movimiento (`ls` directo sobre cada ruta, ver comando en la sesión de tmux):
`src/cpm_stats.py`, `src/pseudobulk.py`, `src/__pycache__/cpm_stats.cpython-311.pyc`,
`src/__pycache__/pseudobulk.cpython-311.pyc`, `src/t3_inventario.py`,
`out/t3_inventario.md`, `out/t3_inventario_run.log`, `out/t3_q6_q7_raw.json`,
`out/t3_queries_raw.json`.

### Un ciclo de sync después, y `aws s3 ls` con las rutas nuevas ya subidas

Esperado un ciclo completo del `syncout.service` (el bucle duerme 120s entre
corridas). Confirmado:

```
$ aws s3 ls s3://vc-challenge-data/vcc/benchmark/
                           PRE data/
                           PRE out/
                           PRE src/
                           PRE work/
2026-10-29 01:55:38        818 LEEME.txt
```

`benchmark/` ya subió completo (sus cuatro subcarpetas más `LEEME.txt`) bajo
`s3://vc-challenge-data/vcc/benchmark/`, con fecha de sync 2026-10-29 01:55:38 —
minutos después del `mv` local. Las rutas viejas bajo `vcc/out/` y `vcc/src/`
siguen presentes también (punto 8), tal como se esperaba.

---

## PARA

**Establecido:**
- La reestructuración de S3 a dos prefijos NO está completa — la raíz sigue
  teniendo 4 prefijos (`datos/`, `out/`, `vcc/`, `virtual cell/`), sin cambios
  desde el estado que dejó documentado `reestructura.md`. No se tocó, como pedía
  el encargo.
- Los 488 archivos de `~/vcc` quedaron clasificados sin ambigüedad en 479 casos
  (449 BENCHMARK, 21 INFRA, 9 COMPARTIDO) y con ambigüedad declarada, no adivinada,
  en 9 (`SIN_CLASIFICAR`, todos ligados a `cpm_stats.py`/`pseudobulk.py` o a la
  familia `t3_inventario`).
- Los 449 archivos BENCHMARK se movieron a `benchmark/{out,src,work,data}/` con
  `mv` puro — verificado por conjunto exacto de rutas, no solo por conteo: ni un
  archivo de más, ni uno de menos. `corpus/`, `registro/`, `PLAN.md`, `mf.sh`,
  `workflow/` y los 9 `SIN_CLASIFICAR` no se tocaron.
- `benchmark/LEEME.txt` y la entrada 10 de `registro/limitaciones.md` quedaron
  escritos con la regla de aislamiento físico, declarando qué transfiere (solo las
  tres definiciones) y qué no (ningún número).
- `syncout.service` sigue activo, sin `--delete`, sin modificarse — y el bucket ya
  muestra, verificado en vivo, tanto las rutas nuevas (`vcc/benchmark/...`) como
  las viejas (`vcc/out/...`, `vcc/src/...`) coexistiendo, exactamente como el
  encargo anticipaba.
- Un incidente propio (dos archivos de trabajo copiados por error a la raíz de
  `~/vcc`) se detectó por la propia verificación de conteo exigida por el encargo
  y se corrigió antes de continuar — exactamente el mecanismo de seguridad que esa
  verificación estaba diseñada para dar.

**Abierto, con dependencia:**
- Los 9 `SIN_CLASIFICAR` necesitan una decisión humana (o el documento maestro que
  el LEEME.txt menciona) para saber si `cpm_stats.py`/`pseudobulk.py` son
  infraestructura reusable por el modelo o utilidad interna del benchmark, y si la
  familia `t3_inventario` pertenece a `corpus/` (inventario) o a `benchmark/`
  (toca paneles de 2025 en algunas de sus consultas). Quedan donde están hasta esa
  decisión.
- El escritor no identificado de `out/` en la raíz del bucket (documentado en
  `reestructura.md`) sigue sin identificarse — no se investigó de nuevo aquí,
  fuera del alcance de este encargo.
- `virtual cell/` sigue sin retirarse (2 de 5 claves no verifican) — mismo estado
  que `limitaciones.md` entrada 9 ya documenta, sin cambios.

---

# ADDENDUM 2026-10-29 — cierre de los 9 `SIN_CLASIFICAR`

Decisión tomada y ejecutada con `mv` únicamente. Nada se borró. No se tocó
`/datos` ni la unidad de systemd. No se investigó el escritor no identificado de
`out/` (sigue aparcado). No se reintentó la reestructuración del bucket a dos
prefijos.

## Disposición final de los 9, archivo por archivo

| # | archivo | acción | destino | clase final |
|---|---|---|---|---|
| 1 | `src/cpm_stats.py` | `mv` | `benchmark/src/cpm_stats.py` | BENCHMARK |
| 2 | `src/pseudobulk.py` | `mv` | `benchmark/src/pseudobulk.py` | BENCHMARK |
| 3 | `src/__pycache__/cpm_stats.cpython-311.pyc` | `mv` | `benchmark/src/cpm_stats.cpython-311.pyc` | BENCHMARK |
| 4 | `src/__pycache__/pseudobulk.cpython-311.pyc` | `mv` | `benchmark/src/pseudobulk.cpython-311.pyc` | BENCHMARK |
| 5 | `src/t3_inventario.py` | sin mover | `src/t3_inventario.py` (sin cambio) | COMPARTIDO |
| 6 | `out/t3_inventario.md` | sin mover | `out/t3_inventario.md` (sin cambio) | COMPARTIDO |
| 7 | `out/t3_inventario_run.log` | sin mover | `out/t3_inventario_run.log` (sin cambio) | COMPARTIDO |
| 8 | `out/t3_q6_q7_raw.json` | sin mover | `out/t3_q6_q7_raw.json` (sin cambio) | COMPARTIDO |
| 9 | `out/t3_queries_raw.json` | sin mover | `out/t3_queries_raw.json` (sin cambio) | COMPARTIDO |

**Razón de los 4 (`BENCHMARK`, movidos):** se escribieron contra el objeto de
Replogle y pueden llevar una convención suya incrustada. Lo que transfiere es la
DEFINICIÓN de pseudobulk, que está en el documento maestro — no este código.

**Razón de los 5 (`COMPARTIDO`, sin mover):** una cifra de cobertura es un conteo
de qué contiene qué archivo, no una medición contra Replogle. Es inventario, y es
el primitivo que hace imposible el error de «panel equivocado».

**Regla que generó ambas decisiones**, ya añadida como entrada 10 de
`registro/limitaciones.md`: un script transfiere solo si lee sus convenciones del
objetivo en tiempo de ejecución, no si las lleva incrustadas; ante la duda se mueve
a `benchmark/`, porque reescribir es más barato que auditar.

## Verificación, releída del disco después del `mv`

| chequeo | antes | después | ¿coincide lo esperado? |
|---|---|---|---|
| archivos totales en `~/vcc` (`find . -type f \| wc -l`) | 492 | 492 | sí — mismo número |
| archivos en `src/` | 3 (`cpm_stats.py`, `pseudobulk.py`, `t3_inventario.py`) | 1 (`t3_inventario.py`) | sí |
| archivos en `src/__pycache__/` | 2 (`cpm_stats.cpython-311.pyc`, `pseudobulk.cpython-311.pyc`) | 0 — directorio vacío, **no borrado** | sí |
| archivos en `benchmark/src/` con esos nombres | 0 | 4 (`cpm_stats.py`, `pseudobulk.py`, `cpm_stats.cpython-311.pyc`, `pseudobulk.cpython-311.pyc`) | sí |
| los 5 `COMPARTIDO` en su ruta original | 5 | 5 (sin cambio) | sí |
| `SIN_CLASIFICAR` restantes | 9 | 0 | sí |

`src/__pycache__/` quedó vacío tras el `mv` de los dos `.pyc` — se reporta aquí tal
como pide el encargo, no se eliminó el directorio.

Entrada 10 de `registro/limitaciones.md` actualizada con esta resolución, las dos
razones y la regla, bajo un nuevo apartado «Resolución de los 9 `SIN_CLASIFICAR` —
2026-10-29».

**Estado tras este addendum:** 0 `SIN_CLASIFICAR` pendientes. Los otros tres puntos
abiertos de la sección PARA (reestructuración de S3 a dos prefijos, escritor no
identificado de `out/`, `virtual cell/` sin retirarse) quedan exactamente como
estaban — no formaban parte de este encargo.
