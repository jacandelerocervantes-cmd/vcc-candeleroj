# Limitaciones registradas

Cinco campos por entrada: **qué falta** · **por qué** · **qué costaría cerrarlo** ·
**afecta al resultado o solo a su interpretación** · **fecha**.

---

## 1. Desfase de versión en `cell-eval2` (0.15.0 vs 0.20.0)

- **Qué falta:** no podemos reproducir, bit a bit, con el clon actual, los `b` y `r`
  publicados en `docs/vcc2026_metrics/vcc2026-metrics-brief.md:341-346` — esos números
  se midieron con `cell-eval2 0.15.0`, y el clon que tenemos instalado trae `0.20.0`.
- **Por qué:** entre ambas versiones puede haber cambiado cualquier cosa que afecte la
  medición empírica de `r` (la partición en mitades, el `rule_version`, el pipeline de
  cómputo) — no hay garantía de que `0.20.0` reproduzca exactamente esos rangos.
- **Qué costaría cerrarlo:** instalar o pinear `cell-eval2==0.15.0` (o el commit
  correspondiente) y re-ejecutar la medición de `b`/`r` sobre los tres bundles oficiales
  para contrastar contra la tabla publicada.
- **Afecta al resultado o solo a su interpretación:** **NO afecta al benchmark de
  literatura** — los números publicados en el brief son los que se citan, y están
  identificados con su versión exacta. Solo afecta a si *nosotros* podemos
  re-derivarlos de forma independiente con el clon actual.
- **Fecha:** 2026-10-27.

---

## 2. Potencia 80 % del MDE

- **Qué falta:** el umbral de potencia (80 %) usado para fijar el MDE (mínimo efecto
  detectable) no aparece en ninguna especificación de `cell-eval2` ni de la
  competencia.
- **Por qué:** es una convención estadística estándar nuestra (la típica en diseño de
  experimentos), no un valor que la especificación del VCC2026 fije o mencione.
- **Qué costaría cerrarlo:** no se puede "cerrar" contra una fuente primaria porque no
  existe tal fuente — solo se podría documentar explícitamente como decisión propia,
  o sustituir por un barrido de sensibilidad sobre varios niveles de potencia.
- **Afecta al resultado o solo a su interpretación:** afecta **solo a la
  interpretación** del MDE resultante (qué tan conservador o agresivo es el umbral),
  no a ningún cálculo oficial del benchmark.
- **Fecha:** 2026-10-27.

---

## 3. α del intervalo conforme

- **Qué falta:** el nivel de significancia `α` para el intervalo de predicción
  conforme no está fijado todavía.
- **Por qué:** es un parámetro de nuestro propio método de cuantificación de
  incertidumbre, no algo que la especificación de `cell-eval2` determine (esa `α` es
  distinta de la `α = 0.05` de BH de la fila 3 de `CONSTANTES_A_VERIFICAR.md`, que es
  interna al scorer).
- **Qué costaría cerrarlo:** una decisión de diseño — fijar `α` (p. ej. 0.10 o 0.05)
  según la tolerancia a intervalos más anchos vs. cobertura, y documentarla.
- **Afecta al resultado o solo a su interpretación:** afecta **solo a la
  interpretación** de los intervalos (su ancho y cobertura nominal), no a las métricas
  del benchmark en sí.
- **Fecha:** 2026-10-27.

---

## 4. Datasets en Seurat/RDS y tar, no en `h5ad`

- **Qué falta:** `GSE225775`, `GSE139944` y `Zenodo 14518762` se distribuyen en
  formato Seurat/RDS (y empaquetados en `tar`), no en `h5ad` — el formato que el
  pipeline de `cell-eval2` consume.
- **Por qué:** son los formatos en los que los autores originales depositaron los
  datos; no hay control sobre eso.
- **Qué costaría cerrarlo:** una conversión Seurat → AnnData/`h5ad` (vía
  `SeuratDisk`/`sceasy` u otro puente R↔Python), además de desempaquetar los `tar`.
  Es trabajo mecánico pero no trivial: hay que preservar metadata de células/genes y
  verificar que la conversión no introduzca discrepancias de normalización.
- **Afecta al resultado o solo a su interpretación:** potencialmente **afecta al
  resultado** si la conversión se hace mal (p. ej. pierde capas de conteos crudos o
  normalización ya aplicada) — es una limitación de disponibilidad de datos, no de
  interpretación.
- **Fecha:** 2026-10-27.

---

## 5. `Bulk_RNAseq_Seurat_object_IFNG_and_TGFB_stim.rds` — ~~no confirmado~~ RESUELTO 2026-10-27

- **Qué falta:** nada — **confirmado como descargado**. `dl.log` registra
  `05:53 otros/Bulk_RNAseq_Seurat_object_IFNG_and_TGFB_stim.rds rc=0 3.3M`, el archivo
  está en `/datos/ref/otros/`, y aparece listado como ítem 10 en
  `/datos/ref/otros/A_readme.txt` (el depósito es el Zenodo 14518762 del manuscrito
  «Systematic reconstruction of molecular pathway signatures using scalable
  single-cell perturbation screens»). 3.3 MB registrado vs 3.4 MB esperado: dentro del
  margen de redondeo.
- **Por qué (del hallazgo original):** no se había verificado explícitamente contra el
  listado de descargas en el momento de escribir esta entrada.
- **Qué costaría cerrarlo:** cerrado — verificado por inspección directa de `dl.log` y
  del directorio (paso 0, inventario).
- **Afecta al resultado o solo a su interpretación:** ninguno de los dos — la
  limitación no existía; era una brecha de verificación nuestra, no de los datos.
- **Fecha:** apertura 2026-10-27 · cierre 2026-10-27 (mismo día, en el paso 0 de
  inventario).

---

## 6. `GSE225775`, muestras GSM7056149 y GSM7056151: línea celular no determinable

- **Qué falta:** la línea celular de estas dos muestras (de 4 en total) no se puede
  determinar con los archivos de metadatos extraídos (`cell.annotations`,
  `hashTable.out`, `hash_sample_sheet`) — ninguno contiene un token de nombre de
  línea reconocible, a diferencia de las otras dos muestras del mismo GSE (A172 y
  GSC0827, sí confirmadas). Además, el recuento de llamadas `rank==1` de su
  `hashTable.out` da 133,528,266 y 3,334,085 respectivamente — demasiado alto para
  ser "una fila por célula final", así que tampoco se puede dar un nº de células
  real para estas dos.
- **Por qué:** el diseño de hash de estas dos muestras («sciPlexGxE_2»,
  «sciPlex_4») no incluye el nombre de línea en la etiqueta de pozo (solo
  `plate9_H3`, `plate6_0_vehicle_0_rep1`, etc., sin línea) — a diferencia de las
  otras dos muestras del mismo lote de descarga, que sí la incluyen.
- **Qué costaría cerrarlo:** conseguir la tabla de diseño experimental de GEO
  (serie `GSE225775`, metadata SOFT/MINiML) o el material suplementario del
  manuscrito — no está en los archivos que tenemos en `/datos`, requeriría bajar
  algo nuevo.
- **Afecta al resultado o solo a su interpretación:** no afecta al resultado
  principal (ya se confirmó que `GSE225775` no aporta ninguna línea de
  nadig/replogle, vía las otras dos muestras) — solo deja sin resolver si estas dos
  muestras *también* no coinciden, o si una de ellas sí trae alguna línea de interés
  sin poder confirmarlo.
- **Fecha:** 2026-10-27.

---

## 7. `GSE139944`, GSM4150378: alta tasa de ambigüedad en la asignación célula→línea

- **Qué falta:** de 24,751,262 barcodes únicos con al menos una llamada `rank==1` en
  el `hashTable.out`, 7,755,491 (31 %) tienen llamadas `rank==1` empatadas a ≥2
  placas de **distinto** `cell_type` — no se pueden asignar a una sola línea sin una
  regla de desempate adicional.
- **Por qué:** el archivo de hash deconvolution de esta muestra trae empates
  genuinos (múltiples etiquetas con el mismo conteo máximo y `rank=1`) — no es un
  error del script, es una propiedad del archivo tal como se nos dio.
- **Qué costaría cerrarlo:** una regla de desempate explícita (p. ej. descartar o
  usar un segundo criterio como el conteo UMI total por barcode) — no se improvisó
  aquí porque inventar la regla sin el pipeline original del paper sería "estimar",
  justo lo que se pidió evitar.
- **Afecta al resultado o solo a su interpretación:** **no afecta la conclusión
  principal** — las 728,069 células K562 sin ambigüedad ya bastan para confirmar la
  coincidencia de línea a nivel de célula. Afecta solo a qué fracción del archivo se
  podría usar si se quisiera aprovechar más células de esta fuente.
- **Fecha:** 2026-10-27.

---

## 8. `rds2py` no pudo leer los 5 objetos Seurat grandes de Zenodo

- **Qué falta:** `rds2py` (lectura de `.rds` sin R) falló en los 5 objetos
  `Seurat_object_*_Perturb_seq.rds`, con
  `RdsParserError: ... cannot read unknown (or unsupported) SEXP type 3` dentro del
  slot `@commands` (historial de llamadas a funciones de Seurat, no los datos).
  Sí funcionó en el objeto bulk pequeño (`Bulk_RNAseq_Seurat_object_IFNG_and_TGFB_stim.rds`).
- **Por qué:** `rds2py` parsea el archivo `.rds` completo de un tirón y no tiene
  forma de omitir un slot concreto que no sepa interpretar — un SEXP (closure o
  entorno de R, típico de `@commands`) que no está en su lista de tipos soportados
  tira todo el parseo.
- **Qué costaría cerrarlo:** no se cerró porque no hizo falta — el camino de
  respaldo (R base + `readRDS()` + `methods::slot()`, sin el paquete `Seurat`)
  funcionó en los 5 objetos. Se deja registrado porque es información útil para la
  próxima vez que se intente leer un `.rds` de Seurat sin R.
- **Afecta al resultado o solo a su interpretación:** ninguno — el camino de
  respaldo ya dio el resultado completo. Es una nota sobre la herramienta, no sobre
  los datos.
- **Fecha:** 2026-10-27.

---

## 9. Regla fija de estructura de S3 — un prefijo por sistema de archivos

- **Regla:** un prefijo en la raíz de `s3://vc-challenge-data` por cada sistema de
  archivos local, y el nombre del prefijo es el punto de montaje. `datos/` = `/datos`,
  `vcc/` = `/home/ubuntu/vcc`. Nada en la raíz del bucket que no sea un punto de
  montaje.
- **Por qué:** la raíz llegó a tener cinco prefijos (`datos/`, `out/`, `scratch_techo/`,
  `vcc/`, `virtual cell/`) por dos encargos distintos escribiendo a destinos distintos
  sin esta regla — `out/` duplicaba `vcc/out/` completo, y `scratch_techo/` era un
  subdirectorio de `/datos` puesto en la raíz por error. La próxima vez que algo
  necesite un sitio en S3, la respuesta sale de dónde vive el archivo en disco, no de
  una decisión nueva.
- **Qué queda abierto:** `virtual cell/` (5 claves, 19.3 KiB) no se pudo retirar en el
  cierre del 2026-10-28 — 2 de sus 5 claves no superaron la verificación de
  "mismo tamaño que su equivalente en `vcc/`" (`registro/CONSTANTES_A_VERIFICAR.md`:
  5,333B en S3 vs. 7,699B local/`vcc/`; `LEGADO.txt`: sin ningún equivalente local).
  Por la regla de seguridad del encargo, no se borró ninguna de las 5. Ver
  `~/vcc/out/reestructura.md` para el detalle completo.
- **Afecta al resultado o solo a su interpretación:** ninguno — es estructura de
  respaldo, no datos de análisis.
- **Fecha:** 2026-10-28.

---

## 10. Aislamiento físico de `benchmark/` — ningún número de ahí construye el modelo

- **Regla:** `~/vcc/benchmark/` (subcarpetas `out/`, `src/`, `work/`, `data/`) queda
  aislada por diseño. Contiene todo el trabajo de validación de los tres vectores del
  benchmark VCC 2026 contra Replogle 2022 sobre genes esenciales — reportes, sus
  `*_raw.json`/`*_raw.csv`, los scripts que los produjeron, `techo_metricas_work/`
  completo, y los parquet/`.npz` de medición (`t13_replicas.parquet`,
  `vcc2025_knockdown_raw_arrays.npz`,
  `vcc2025_rejilla_raw_arrays_{test_train,validation}.npz`). Su propio `LEEME.txt`
  repite la regla.
- **Por qué:** todo lo que hay en esa carpeta son números MEDIDOS contra un dato
  externo (Replogle 2022 o el dato de 2025) — baselines, anclas, umbrales,
  medianas, `n_real` por gen. Ninguno de esos números debe entrar en la
  construcción del modelo (A/B/C ni D/E/F). El aislamiento es FÍSICO (carpeta
  separada) en vez de depender de que cada paso futuro recuerde la regla por su
  cuenta — una carpeta no se olvida, una convención sí.
- **Qué transfiere, y de dónde:** únicamente las tres DEFINICIONES (no sus
  números), escritas en el documento maestro, no en `benchmark/`: el perfil de
  conteos como pseudobulk de las células de la perturbación contra el pool de
  control; el conjunto DE por un umbral de significancia sobre un test declarado,
  con población y convención de multiplicidad leídas del objetivo; la magnitud
  como cociente lineal contra el control, resumida por la mediana.
- **Qué queda sin clasificar, y por qué:** ~~`src/cpm_stats.py`, `src/pseudobulk.py`
  (y sus `.pyc` en `__pycache__/`), y la familia `t3_inventario.{py,md}` +
  `t3_inventario_run.log` + `t3_q6_q7_raw.json` + `t3_queries_raw.json` se quedan
  donde estaban — no se pudo decidir con certeza si son infraestructura reusable o
  utilidad interna del propio benchmark, y un archivo mal movido es peor que uno
  sin mover.~~ **RESUELTO 2026-10-29** — ver detalle abajo. Ver
  `~/vcc/out/estado_actual.md` para el detalle completo de la clasificación
  original.
- **Afecta al resultado o solo a su interpretación:** ninguno — es organización de
  carpetas, no un cambio de ningún número ya publicado en un reporte.
- **Fecha:** 2026-10-29.

### Resolución de los 9 `SIN_CLASIFICAR` — 2026-10-29

Decisión tomada y ejecutada con `mv` únicamente (nada borrado):

1. **`src/cpm_stats.py`, `src/pseudobulk.py` y sus `.pyc` (4 archivos) → movidos a
   `benchmark/src/`.** Razón: se escribieron contra el objeto de Replogle y pueden
   llevar una convención suya incrustada. Lo que transfiere es la DEFINICIÓN de
   pseudobulk (ver arriba, «Qué transfiere, y de dónde»), que está en el documento
   maestro — no este código.
2. **`src/t3_inventario.py`, `out/t3_inventario.md`, `out/t3_inventario_run.log`,
   `out/t3_q6_q7_raw.json`, `out/t3_queries_raw.json` (5 archivos) → se quedan donde
   estaban, reclasificados como `COMPARTIDO`.** Razón: una cifra de cobertura es un
   conteo de qué contiene qué archivo, no una medición contra Replogle. Es
   inventario, y es el primitivo que hace imposible el error de «panel equivocado».

**Regla que generó ambas decisiones**: un script transfiere solo si lee sus
convenciones del objetivo en tiempo de ejecución, no si las lleva incrustadas; ante
la duda se mueve a `benchmark/`, porque reescribir es más barato que auditar.

Verificado releído del disco: conteo de archivos de `~/vcc` antes y después del
`mv`, el mismo número (492 = 492); los 4 confirmados en `benchmark/src/` y ya no en
`src/`; los 5 confirmados en su ruta original; `src/` quedó con 1 archivo
(`t3_inventario.py`) y `src/__pycache__/` quedó vacío (no se borró, solo se
reporta); ningún `SIN_CLASIFICAR` queda pendiente. Detalle completo del addendum en
`~/vcc/out/estado_actual.md`.
