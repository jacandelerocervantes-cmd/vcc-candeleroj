# Paso 1 de A/B/C — lectura de la entrada y validación del contrato de salida

Todo este encargo corrió fuera de `benchmark/` — no se leyó nada de ahí. No se
envió nada: ni `vcc sample`, ni `vcc submit`. Solo `vcc prep --dry-run`, que
valida en local sin tocar la red. No se modificó ningún archivo de
`/datos/reto`. Scripts en `modelo/src/`, salidas en `modelo/out/`, ejecutados
fuera de `/datos`, con el intérprete de `~/obs_venv` en modo aislado (`python
-I`), paths siempre pasados como argumento.

## Condición de parada, declarada antes de correr nada

Si cualquiera de estos resultados de la Parte 1 no se cumplía, el plan era
pararse ahí mismo — sin pasar a la Parte 2 — porque habría significado que
leemos mal el objeto, no que el contrato tiene un matiz:

- `n_obs != 18400` en cualquiera de los seis `context_{A..F}.h5ad`.
- `ntc_id` con un número de categorías distinto de 46, o alguna categoría con
  un conteo distinto de 400.
- `target_gene` con más de una categoría, o con una categoría distinta de
  `"non-targeting"`.
- El valor de `context` dentro de un archivo no coincide con su letra, o no
  coincide con `manifest.json`.
- El orden de `var` no es **idéntico, elemento por elemento**, al de
  `gene_names.csv` (no solo el conjunto o el conteo).
- `X` contiene valores no enteros, negativos, o no finitos.

Ninguno de estos seis disparadores se activó. Se procedió a la Parte 2.

---

## Parte 1 — ¿leemos bien la entrada?

Método: `h5py.File(path, "r")`, leyendo `/obs`, `/var` y los arrays crudos de
`/X` (`data`, `indices`, `indptr`) de un CSR ya existente en disco — sin
densificar nunca la matriz completa. Script:
`modelo/src/paso1_leer_entrada.py`. Log crudo (JSON completo, los seis
archivos): `modelo/out/paso1_raw.json`.

### a)–e) — forma, categorías y orden, los seis archivos

| # | contexto | archivo | n_obs medido | n_obs esperado | ¿OK? |
|---|---|---|---|---|---|
| 1 | A | `val/context_A.h5ad` | 18400 | 18400 | sí |
| 2 | B | `val/context_B.h5ad` | 18400 | 18400 | sí |
| 3 | C | `val/context_C.h5ad` | 18400 | 18400 | sí |
| 4 | D | `final/context_D.h5ad` | 18400 | 18400 | sí |
| 5 | E | `final/context_E.h5ad` | 18400 | 18400 | sí |
| 6 | F | `final/context_F.h5ad` | 18400 | 18400 | sí |

| # | chequeo | medido | esperado | ¿OK? |
|---|---|---|---|---|
| 7 | `ntc_id` nº categorías, los 6 archivos | 46 en cada uno | 46 | sí |
| 8 | `ntc_id` células por categoría, los 6 archivos | mín=400, máx=400 en cada uno | exactamente 400 | sí |
| 9 | `target_gene` nº categorías, los 6 archivos | 1 en cada uno | 1 | sí |
| 10 | `target_gene` valor de la categoría, los 6 archivos | `"non-targeting"` en cada uno | `"non-targeting"` | sí |
| 11 | `context` valor dentro de cada archivo | A→"A", B→"B", C→"C", D→"D", E→"E", F→"F" | la letra del propio archivo | sí |
| 12 | `context` vs `manifest.json` (`contexts` por partición) | val=[A,B,C], final=[D,E,F] | igual | sí |
| 13 | `n_vars`, los 6 archivos | 18533 en cada uno | 18533 | sí |
| 14 | orden de `var._index` vs `gene_names.csv`, elemento por elemento, los 6 archivos | `identical_order=true`, `first_mismatch=null` en cada uno | idéntico | sí |

### f) — mediana de UMIs por célula (la página dice «about 20,000»)

| # | contexto | mediana medida | media medida | ¿consistente con «~20,000»? |
|---|---|---|---|---|
| 15 | A | 20,109.0 | 21,133.5 | sí |
| 16 | B | 19,946.0 | 19,989.9 | sí |
| 17 | C | 20,034.0 | 21,156.5 | sí |
| 18 | D | 20,235.5 | 21,249.6 | sí |
| 19 | E | 19,681.0 | 20,921.2 | sí |
| 20 | F | 19,845.0 | 21,003.7 | sí |

### g) — `X`: enteros, min/max, dtype real

| # | contexto | `X` dtype real | min valor almacenado | max valor almacenado | ¿todos enteros? (comprobado sobre el array completo, no una muestra) |
|---|---|---|---|---|---|
| 21 | A | `float32` | 1.0 | 977.0 | sí |
| 22 | B | `float32` | 1.0 | 1753.0 | sí |
| 23 | C | `float32` | 1.0 | 2462.0 | sí |
| 24 | D | `float32` | 1.0 | 3956.0 | sí |
| 25 | E | `float32` | 1.0 | 4104.0 | sí |
| 26 | F | `float32` | 1.0 | 3692.0 | sí |

**Nota:** el dtype real es `float32`, no un entero nativo — pero cada valor
almacenado es entero (`data == round(data)` para los ~100–110 millones de
entradas de cada archivo, no solo una muestra). El `min` de 1.0 en los seis es
esperado: es un CSR que no guarda ceros explícitos, así que el mínimo
almacenado nunca puede ser 0.

### h) — no-ceros por célula (decide si la carga trivial cabe en el tope de dispersión)

| # | contexto | nnz mediana | nnz p95 | nnz máx | nnz mín |
|---|---|---|---|---|---|
| 27 | A | 6147 | 7619.0 | 8453 | 2205 |
| 28 | B | 5756 | 6839.0 | 7649 | 580 |
| 29 | C | 6006 | 7549.0 | 8482 | 1937 |
| 30 | D | 5891.5 | 7316.0 | 8146 | 787 |
| 31 | E | 5634 | 7190.0 | 8107 | 1607 |
| 32 | F | 5630 | 7044.0 | 7914 | 1571 |

Con una mediana de ~5600–6150 no-ceros por célula, muy por debajo del
presupuesto de ~13,200/célula que impone el tope de 4,750,000,000 (ver Parte 2),
una carga que reutiliza filas de control tal cual (sin densificar) tiene
margen amplio para caber.

### i) — `final/context_{D,E,F}.h5ad` y su `pert_counts.csv`, contra `val`

| # | chequeo | medido | ¿OK? |
|---|---|---|---|
| 33 | misma forma `[18400, 18533]` en D/E/F | sí, igual que A/B/C | sí |
| 34 | mismos 18,533 genes, mismo orden, en D/E/F | `identical_order=true` en los tres | sí |
| 35 | `gene_names.csv` de `val/` vs `final/`, byte a byte | idénticos (`diff` vacío) | sí |
| 36 | nº de perturbaciones en `pert_counts.csv` | val=300, final=300, sin duplicados en ninguno | sí |
| 37 | ¿`pert_counts.csv` de val y final son el mismo conjunto de 300? | **no** — `val_eq_final_as_sets=false`, solapan 129 de 300 (171 exclusivas de cada lado) | — (son paneles distintos, no un error) |
| 38 | ¿los 46 `ntc_id` coinciden entre val y final, o son guías distintas? | **coinciden exactamente** — mismo conjunto de 46, 0 exclusivos de cada lado | — (mismas guías de control reutilizadas en ambas particiones) |

**Resultado de la Parte 1: los 38 puntos anteriores pasaron. No se disparó la
condición de parada. Se procedió a la Parte 2.**

---

## Parte 2 — el contrato de salida, probado con una carga trivial

**Carga trivial construida** (script `modelo/src/paso1_contrato_trivial.py`):
para cada una de las 300 perturbaciones de `val/pert_counts.csv` y cada uno de
los contextos A, B, C, se muestrearon 400 índices de fila **con reemplazo**
sobre las 18,400 células de control de ESE contexto (`val/context_{A,B,C}.h5ad`),
con semillas fijas y documentadas (`A=1029, B=1030, C=1031`), y se copió su `X`
tal cual — ni un conteo modificado. `obs.context` reutiliza literalmente la
letra que trae el archivo de control de origen; `obs.target_gene` es el
símbolo de `pert_counts.csv`, nunca un id de guía (verificado: 0 entradas de
`pert_counts.csv` con forma `SIMBOLO-número`, el patrón de `ntc_id`).

### Cálculo previo a escribir nada — entradas almacenadas previstas

| # | chequeo | valor | tope / referencia | ¿cabe? |
|---|---|---|---|---|
| 39 | nnz previsto, contexto A (de `indptr`, sin tocar `data`) | 717,022,642 | — | — |
| 40 | nnz previsto, contexto B | 661,574,281 | — | — |
| 41 | nnz previsto, contexto C | 705,732,356 | — | — |
| 42 | **nnz previsto total** | **2,084,329,279** | tope duro **4,750,000,000** (`HARD_CAP_NNZ` en `vcc/sizing.py` del propio CLI instalado — confirmado leyendo el paquete, no supuesto) | sí, con margen ×2.28 |
| 43 | nnz previsto por célula (promedio) | 5,789.8 | presupuesto por célula = 4.75e9 / 360,000 = 13,194.4 | sí, 44 % del presupuesto |
| 44 | espacio libre en disco en `/` antes de escribir | 14 GiB | tamaño estimado por compresión gzip-4 igual que los archivos de origen (~2.05 bytes/entrada comprimida) ≈ 4.3 GB | sí, con margen |

El cálculo se hizo y se reportó **antes** de construir la matriz de
360,000×18,533 — el guion se detiene con `sys.exit(2)` si `fits_under_cap` sale
`false`; no hizo falta, salió `true`.

### Construcción y escritura — verificado, releído tras escribir

| # | chequeo | medido | esperado | ¿OK? |
|---|---|---|---|---|
| 45 | forma de `X` final | `(360000, 18533)` | `300×400×3 = 360000` filas, 18533 columnas | sí |
| 46 | nnz real de `X` tras construir, vs lo previsto en el punto 42 | 2,084,329,279 | igual (son las mismas filas, copiadas, no una estimación) | sí, exacto |
| 47 | UMI máximo observado en las filas muestreadas | 52,420 | < 1,000,000 | sí |
| 48 | tamaño del archivo escrito, `modelo/out/prueba.h5ad` | 4,087,112,150 bytes (3.81 GiB), gzip nivel 4 | cupo en los 14 GiB libres | sí, quedaron ~11 GiB libres |
| 49 | filas `non-targeting` en el `obs` de salida | 0 | 0 | sí |
| 50 | perturbaciones únicas en `obs.target_gene` | 300 | 300 (las de `pert_counts.csv`) | sí |
| 51 | células por perturbación, por contexto | 400 en los 300×3 = 900 cruces perturbación×contexto, sin excepción | exactamente 400 | sí |

### Validación con `vcc prep --dry-run` — salida textual

Comando corrido (desde `~/vcc`, `-g`/`--perts` apuntando a los CSV reales de
`/datos/reto/val/`, sin copiarlos):

```
vcc prep modelo/out/prueba.h5ad -g /datos/reto/val/gene_names.csv \
  --perts /datos/reto/val/pert_counts.csv -o modelo/out/prueba.vcc --dry-run
```

Salida textual completa, verbatim (guardada en
`modelo/out/paso1_prep_dryrun_1.log`):

```
✓ dry run — no .vcc written
  input:         modelo/out/prueba.h5ad
  output:        (dry run)
  cells:         360000
  genes:         18533  (encoding: float32)
  perturbations: column 'target_gene' → 'target_gene'
  contexts:      column 'context' → A=120000, B=120000, C=120000
  targets:       verified against the official list
  normalization: counts-preserved
```

**Pasó, al primer intento, sin ajustes.** No hubo mensaje de error que
reportar — no se violó ninguna cláusula.

La misma corrida con `--json` (`modelo/out/paso1_prep_dryrun_1.json`) confirma
los números exactos, cruzados con la Parte 2:

| # | campo del JSON de `prep` | valor | cruce con lo ya calculado |
|---|---|---|---|
| 52 | `nnz` | 2,084,329,279 | igual al punto 46 — el propio `prep` recontó el nnz real y coincide |
| 53 | `n_genes` | 18533 | igual al punto 13/45 |
| 54 | `reordered_genes` | `false` | confirma que el orden de `var` ya coincidía con `gene_names.csv` — no hubo que reordenar |
| 55 | `dropped` | `[]` | 0 filas descartadas — confirma el punto 49 (cero `non-targeting`) desde el lado del validador, no solo del script propio |
| 56 | `verified_targets` | `true` | las 300 perturbaciones, por contexto, verificadas contra el panel oficial |
| 57 | `notes` | `[]` | sin advertencias de ningún tipo |
| 58 | `encoding` | `32` | `X` se empaqueta como `float32` en el `.vcc` — el contrato pide «enteros» como *valor*, no un dtype entero nativo; el propio validador oficial lo acepta así (coincide con cómo vienen los archivos de origen, ver punto 21-26) |

### Las 8 cláusulas del contrato, una por una

| # | cláusula del contrato | chequeo | resultado |
|---|---|---|---|
| 59 | un solo `.h5ad`, tres contextos, columna `context` reutilizada exactamente | `obs.context` toma literalmente "A"/"B"/"C" del archivo de control de origen, sin re-etiquetar | cumple |
| 60 | 360,000 filas × 18,533 columnas | medido en el punto 45 | cumple |
| 61 | `obs` con `target_gene` (símbolo) y `context` | columnas presentes; `target_gene` verificado sin ids de guía (0 coincidencias con el patrón `SIMBOLO-nº`) | cumple |
| 62 | cero filas `non-targeting` | punto 49 (script propio) + punto 55 (`prep --json`, `dropped=[]`) | cumple, doble verificado |
| 63 | `var` en el orden de `gene_names.csv` | punto 54 (`reordered_genes=false`) | cumple |
| 64 | `X`: conteos crudos, no negativos, enteros, finitos | `X_total.data.min() >= 0` (assert), todos los valores heredados de archivos ya verificados enteros (puntos 21-26); `prep` no reportó ninguna nota de valores inválidos | cumple |
| 65 | ninguna célula > 1,000,000 conteos totales | punto 47, máximo observado 52,420 (5.2 % del tope) | cumple |
| 66 | disperso, máx. 4,750,000,000 entradas almacenadas | punto 46, 2,084,329,279 (43.9 % del tope, margen ×2.28) | cumple |
| 67 | exactamente las 300 perturbaciones de `pert_counts.csv`, 400 células c/u, por contexto | punto 51 (script propio) + punto 56 (`prep --json`, `verified_targets=true`) | cumple, doble verificado |

No hubo ninguna cláusula violada, así que no hay mensaje de error que
diagnosticar.

---

## PARA

**Establecido:**
- Los 38 puntos de la Parte 1 confirman que se lee bien la entrada en los seis
  `context_{A..F}.h5ad`: forma exacta (18,400×18,533), categorías exactas
  (`ntc_id`=46×400, `target_gene`=1 categoría), orden de genes idéntico
  elemento por elemento a `gene_names.csv`, mediana de UMIs consistente con la
  cifra «~20,000» de la página del reto, y `X` con valores enteros (aunque en
  dtype `float32`) sin negativos.
- `val` y `final` comparten exactamente los mismos 18,533 genes en el mismo
  orden, y los mismos 46 `ntc_id` (guías de control idénticas en ambas
  particiones) — pero **no** comparten el mismo panel de 300 perturbaciones
  (solo 129 de 300 se solapan). No es un error: son paneles distintos por
  diseño.
- El contrato de salida se probó entero, con una carga deliberadamente trivial
  (resample con reemplazo de los propios controles, contratos A/B/C), y
  `vcc prep --dry-run` lo validó al primer intento: cero cláusulas violadas,
  cero advertencias, `verified_targets=true`, `dropped=[]`.
- El cálculo de entradas almacenadas se hizo y reportó **antes** de escribir
  la matriz completa (2,084,329,279 previstas, confirmado exacto tras
  construir) — con margen ×2.28 contra el tope real de 4,750,000,000
  (`HARD_CAP_NNZ`, confirmado leyendo el código del CLI instalado, no
  asumido).
- `modelo/out/prueba.h5ad` (3.81 GiB) queda en disco como plantilla verificada
  — no se borró, tal como pide el encargo. No se escribió ningún `.vcc`
  (`--dry-run` lo mantuvo así).

**Abierto:**
- La carga es deliberadamente trivial — no predice nada. El siguiente paso de
  A/B/C (no parte de este encargo) es sustituir el muestreo degenerado por una
  predicción real, manteniendo el mismo contrato ya verificado.
- No se tocó `vcc sample` ni se hizo ningún envío real — el encargo lo reserva
  para una orden explícita del autor.
- El hallazgo del punto 58 (contrato pide «enteros», el `.vcc` empaqueta
  `float32`) no bloquea nada — el propio validador oficial lo acepta — pero
  queda anotado por si una versión futura de `prep` endurece esa cláusula.

---

# ADDENDUM 2026-10-29 — mover las emisiones grandes al volumen correcto

**Problema.** `~/vcc` vive en el volumen raíz (~30 GiB totales, ~11 GiB libres
antes de este addendum). `prueba.h5ad` pesa 3.81 GiB y vivía ahí; cada emisión
futura pesará lo mismo, y con dos o tres más se llena la raíz — lo cual rompe
la instancia entera, no solo el modelo. Además `modelo/out/` caía dentro de
`syncout.service`, así que cada emisión se subiría completa a S3 cada vez que
cambiara.

**Solución, sin tocar systemd.** `syncout.service` sincroniza únicamente
`/home/ubuntu/vcc/` y `/datos/_scratch_techo_metricas/` — confirmado leyendo
la unidad en vivo (`systemctl cat syncout.service`), no modificada. Cualquier
otra ruta de `/datos` queda fuera del sync por construcción, así que basta con
poner los artefactos grandes en `/datos/modelo/out/`.

## 1 · Destino creado

```
mkdir -p /datos/modelo/out
```

No hizo falta `sudo`: `/datos` ya es `ubuntu:ubuntu` (`drwxr-xr-x ubuntu
ubuntu`), y `/datos/modelo/out/` quedó con el mismo dueño sin pasos extra.

## 2 · Espacio comprobado antes de mover

| filesystem | libre ANTES | libre necesario |
|---|---|---|
| `/` (raíz, donde vivía `prueba.h5ad`) | 11 GiB | liberará 3.81 GiB |
| `/datos` | 40 GiB | necesita 3.81 GiB |

Sobraba espacio de los dos lados. El `mv` cruza sistemas de archivos (copia +
borrado) — tardó 2.8 segundos reales, no los «un par de minutos» anticipados
(el archivo ya estaba en el *page cache* del sistema por las corridas previas
de este mismo paso).

## 3 · Movido

```
mv ~/vcc/modelo/out/prueba.h5ad /datos/modelo/out/
```

Único archivo de `modelo/out/` por encima de 100 MB — no había otro que mover.
Los reportes/logs/JSON pequeños se quedaron en `~/vcc/modelo/out/`.

## 4 · Script actualizado — ruta de salida como argumento, no constante

En `modelo/src/paso1_contrato_trivial.py`: la constante `OUT_H5AD` incrustada
se reemplazó por un argumento `--out` con `argparse`, cuyo valor por defecto es
`DEFAULT_OUT_H5AD = "/datos/modelo/out/prueba.h5ad"`. `main()` ahora lee
`args.out` en vez de la constante. Verificado con `--help` y con un chequeo de
sintaxis (`ast.parse`) tras el cambio — no se volvió a correr la construcción
completa (ya está escrita y verificada; rehacerla no aporta nada nuevo).

## 5 · Nota dejada en `/datos/modelo/out/LEEME.txt`

Escrita con el contenido exacto dado en el encargo, verbatim (602 bytes) —
incluye las tres semillas documentadas de la carga trivial (`A=1029, B=1030,
C=1031`) y la excepción del `.vcc` real, que sí se respalda antes de enviarse.

## 6 · Verificación, releída del disco

| # | chequeo | medido | ¿OK? |
|---|---|---|---|
| 68 | `df -h /` después del `mv` | 14 GiB libres (antes: 11 GiB) | sí, ronda los «~15 GiB» esperados |
| 69 | `df -h /datos` después del `mv` | 37 GiB libres (antes: 40 GiB) | sí, bajó ~3 GiB por el archivo nuevo |
| 70 | `/datos/modelo/out/prueba.h5ad`, tamaño exacto | 4,087,112,150 bytes | sí, idéntico al que tenía antes del `mv` |
| 71 | `~/vcc/modelo/out/` tras el `mv` | ya no tiene el `.h5ad`; sí tiene los 6 reportes/logs/JSON chicos | sí |
| 72 | dueño de `/datos/modelo/`, `/datos/modelo/out/`, y sus dos archivos | `ubuntu:ubuntu` en los cuatro | sí |
| 73 | `.h5ad` movido, releído con `h5py` sin densificar | `encoding-type=csr_matrix`, `shape=[360000, 18533]`, `n_obs=360000` (vía `obs/_index`), `indptr` de longitud 360,001 | sí, idéntico a antes del `mv` |
| 74 | unidad `syncout.service`, leída en vivo (no tocada) | `ExecStart` sincroniza solo `/home/ubuntu/vcc/` y `/datos/_scratch_techo_metricas/` | confirma la premisa de la solución |
| 75 | un ciclo de sync después, `aws s3 ls s3://vc-challenge-data/datos/` | `_scratch_techo_metricas/`, `ref/`, `reto/`, `vcc2025/` — **sin** `modelo/` | sí, `modelo/` queda fuera del sync, que es justo lo buscado |

## PARA (actualización)

**Establecido:**
- Los artefactos grandes del modelo (`.h5ad`/`.vcc`) ahora van a
  `/datos/modelo/out/`, fuera del volumen raíz y fuera de `syncout.service` —
  verificado en vivo, sin tocar la unidad de systemd.
- `modelo/src/paso1_contrato_trivial.py` ya no incrusta la ruta de salida:
  es un argumento (`--out`) con ese destino como valor por defecto.
- `prueba.h5ad` se movió, no se borró, y sigue siendo el mismo archivo byte a
  byte (mismo tamaño, mismo `shape`, legible sin densificar).
- `/datos/modelo/out/LEEME.txt` documenta por qué viven ahí, que no están
  respaldados a propósito (son regenerables desde el script con sus semillas),
  y la excepción del `.vcc` que sí se envíe de verdad.
- No se tocó `/datos/ref`, `/datos/reto`, `/datos/vcc2025`, ni
  `/datos/_scratch_techo_metricas`. No se leyó nada de `benchmark/`. No se
  tocó la unidad de systemd.

**Abierto:**
- Ninguno nuevo introducido por este addendum. Las emisiones futuras del paso
  1 (y de A/B/C en general) deben escribirse con el nuevo default de `--out`
  —o pasando explícitamente una ruta bajo `/datos/modelo/out/`— para no volver
  a acumular en la raíz.
