# GEO / Zenodo — ¿alguna fuente trae Jurkat, HepG2, K562 o RPE1 a nivel de célula?

**Objetivo.** El tercer nivel de T12 (varianza entre preparaciones) solo se puede medir
si una de las líneas de `nadig`/`replogle` (Jurkat, HepG2, K562, RPE1) aparece también
en una fuente independiente, a nivel de célula. Esta es la única ruta que queda: los
dos `tar` de GEO (`GSE225775`, `GSE139944`) y los 9 `.rds` de Zenodo 14518762.

Tarea larga, corriendo sin pausas entre pasos. Este archivo se actualiza conforme
avanza cada paso, no solo al final.

---

## 0 · Disco, antes de tocar nada

```
$ df -h /datos
Filesystem      Size  Used Avail Use% Mounted on
/dev/nvme1n1    274G  145G  116G  56% /datos
```

116 G libres. No se borra nada de `/datos` en ninguna parte de esta tarea.

**Estado:** ⏳ en curso.

---

## 1 · GEO tars — en curso

Se extrajeron **solo** los miembros de anotación/metadatos de ambos `tar` (nunca
`UMI.count.matrix.gz` ni los `.RDS`/`.rds` embebidos de CDS/Seurat, que cuentan como
"matriz" por el criterio del encargo) a `/home/ubuntu/geo_extract/{GSE225775,GSE139944}/`
— directorio nuevo y aislado, fuera de `/datos` y de donde viven los scripts. Los
archivos siguen comprimidos; se leen por streaming (`gzip` de Python), sin
descomprimir a disco, para no llenar los 19 G libres de `/` con los `hashTable` de
cientos de MB comprimidos.

**Miembros extraídos — `GSE225775`:** `cell.annotations`, `hashTable.out`,
`hash_sample_sheet` de las 4 muestras (`GSM7056148-151`).
**Miembros extraídos — `GSE139944`:** `cell.annotations`, `hashTable.out`,
`hashSampleSheet`, `pData` de las 4 muestras (`GSM4150376-379`), más
`hashTable_metadata` de `GSM4150378` (la única con "K562" en el nombre de archivo).

### 1.1 · Líneas encontradas por inspección directa (nombre de archivo + sample sheet)

| GSM | fuente | línea(s) visibles | evidencia |
|---|---|---|---|
| GSM7056148 (sciPlexGxE_1) | GSE225775 | **A172** (glioma) | `hash_sample_sheet`: `plate10_A1_A172_CRISPRi_MMR_100_temozolomide_96...`; `hashTable` label: `plate16_H2_A172_CRISPRi_HPRT1_0_dmso_96` |
| GSM7056149 (sciPlexGxE_2) | GSE225775 | **sin línea en las etiquetas** | `hash_sample_sheet`: `plate3_A1`, `plate3_B1`, … — sin nombre de línea. `hashTable` label: `plate9_H3` — igual, sin línea. Sample tag en `cell.annotations`: `sci3_CRISPRi_kinome` — tampoco la nombra |
| GSM7056150 (sciPlex_3) | GSE225775 | **GSC0827** (línea de células madre de glioma, paciente-derivada) | `hash_sample_sheet`: `plate15_A1_GSC0827_trametinib_10_rep_1`. Sample tag: `sci3_Plex_GSC` |
| GSM7056151 (sciPlex_4) | GSE225775 | **sin línea en las etiquetas** | `hash_sample_sheet`: `plate6_0_vehicle_0_rep1` — sin línea. Sample tag: `sci3_sciPlex8_FollowUp` |

**Ninguna de las cuatro muestras de `GSE225775` nombra Jurkat, HepG2, K562 ni RPE1** en
ningún archivo de metadatos extraído. Dos de las cuatro (149, 151) ni siquiera
declaran su línea en estos archivos — queda como limitación abierta, no como "no
coincide" (no se puede afirmar lo segundo sin saber la línea real).

### 1.2 · `GSE225775`, conteos reales por muestra (no solo nombre de archivo)

`cell.annotations.txt.gz` cuenta **todas las combinaciones de barcode vistas**, no
células reales (en `GSM7056149` da 43,209,765 — imposible como nº de células de un
solo experimento). La cifra correcta es el **recuento de llamadas `rank==1`** del
`hashTable.out.txt.gz` (la asignación de hash ya resuelta):

| GSM | células `rank==1` | línea | % que mapea a la línea declarada |
|---|---|---|---|
| GSM7056148 (sciPlexGxE_1) | 740,468 | A172 | 100 % (las 740,468 mapean a A172) |
| GSM7056149 (sciPlexGxE_2) | 133,528,266 | **sin determinar** | — (sin token de línea en ninguna etiqueta) |
| GSM7056150 (sciPlex_3) | 368 | GSC0827 | 100 % (las 368) |
| GSM7056151 (sciPlex_4) | 3,334,085 | **sin determinar** | — (sin token de línea en ninguna etiqueta) |

⚠️ 133.5 millones y 3.3 millones de llamadas `rank==1` son cifras **demasiado altas
para ser células reales** incluso sin conocer la línea — son indicio de que, para
estas dos muestras, el `hashTable.out` que se nos dio está contando algo más fino que
"una fila por célula final" (p. ej. combinaciones RT×Lig sin colapsar a célula, o un
umbral de filtrado que no se aplicó). No lo estimo: lo que puedo afirmar con certeza
es que **ninguna fila de esas dos muestras contiene un token de línea reconocible**
(ver §1), así que aun si se resolviera el recuento real, la línea seguiría sin
determinarse con los archivos que tenemos. → **va a `limitaciones.md`.**

**Ninguna de las cuatro muestras de `GSE225775` coincide con Jurkat/HepG2/K562/RPE1.**

### 1.3 · `GSE139944` — resultado completo

| GSM | líneas declaradas (nombre de archivo / etiquetas) | ¿coincide con nadig/replogle? |
|---|---|---|
| GSM4150376 (sciPlex1) | **HEK293T**, **NIH3T3** — confirmado en el propio `hashTable` (`HEK293T_A06`, `NIH3T3_F12`, …) | no |
| GSM4150377 (sciPlex2) | A549 (por nombre de archivo; etiquetas de `hashTable` son solo dosis de fármaco, sin token de línea) | no |
| GSM4150378 (sciPlex3) | **A549, MCF7, K562** — confirmado a nivel de célula (ver abajo) | **sí, K562** |
| GSM4150379 (sciPlex4) | A549, MCF7 — confirmado en el propio `hashTable` (etiquetas tipo `..._A549_plate11_E3`, `..._MCF7_plate8_G12`) | no |

**GSM4150378 — cómo se resolvió a nivel de célula.** Su `hashTable.out.txt.gz` no
nombra la línea directamente (solo un id de placa, p. ej. `plate17`, `plate2`); pero
el archivo separado `..._hashTable_metadata.txt.gz` (4,993 filas = pozos, 17 columnas)
sí tiene una columna `cell_type` explícita, y se verificó que **`cell_type` es
constante dentro de cada `plate_oligo`** (52 placas, 52 pares (placa, tipo) únicos —
ninguna placa mezcla líneas). Eso permite asignar línea a cada célula por su llamada
de placa `rank==1`, vía un join exacto (`awk`, placa → tipo):

| resultado del join (células con llamada `rank==1` sin ambigüedad de tipo) | nº de células |
|---|---|
| **K562** | **728,069** |
| MCF7 | 1,057,709 |
| A549 | 397,317 |
| (sin tipo — placa fuera del mapa de 52) | 2,781,514 |
| **ambiguas** (la misma célula con llamadas `rank==1` empatadas a ≥2 placas de tipos distintos) | 7,755,491 |

**K562 queda confirmado a nivel de célula en `GSE139944`, con 728,069 células sin
ambigüedad** — una segunda fuente independiente de K562, además de los 5 objetos de
Zenodo (ver §2). La alta tasa de ambigüedad (7.75 M de 24.75 M barcodes únicos
clasificados) es real y no se resuelve con los archivos que tenemos — también
registrada en `limitaciones.md`, pero no bloquea la conclusión: 728,069 células K562
sin ambigüedad es, de por sí, suficiente para el nivel 3 de T12.

**Estado §1:** ✅ completo.

## 2 · `.rds` de Zenodo — en curso

**Camino 1, `rds2py` (sin R): funciona para el objeto chico, falla en los grandes.**

- `Bulk_RNAseq_Seurat_object_IFNG_and_TGFB_stim.rds` (3.3 M): **`rds2py` lo leyó
  completo, sin R.** Es un objeto `Seurat` genérico (rds2py no tiene una clase
  `Seurat` nativa, lo devuelve como `dict` anidado fiel a los slots S4 — perfectamente
  usable). **9 "células"** (en realidad 9 muestras de RNA-seq bulk, confirmado por
  los barcodes y por el nombre del archivo). 6 columnas de `meta.data`: `orig.ident`
  (constante, `"SeuratProject"`, no informativo), `nCount_RNA`, `nFeature_RNA`,
  `sample` (ej. `TGFB1_rep2`, `IFNG_rep1`, `control_rep3`), `stim` (`TGFB1`/`IFNG`/
  `control`), `replicate`. **Ninguna columna nombra una línea celular.**
- `Seurat_object_TGFB_Perturb_seq.rds` (2.5 G, el primero de los 5 grandes de
  Perturb-seq): **`rds2py` FALLA.** Error exacto:
  `RdsParserError: failed to parse an S4 object's body — failed to parse pairlist
  element 'commands' → failed to parse list element 2 → failed to parse pairlist
  element 'params' → failed to parse list element 5 → cannot read unknown (or
  unsupported) SEXP type 3`. Es el slot `@commands` de Seurat (el historial de
  llamadas a funciones, no los datos) — `rds2py` parsea el archivo entero de un tirón
  y no tiene forma de saltarse solo ese slot.

**Camino 2, R por apt: funcionó, para los 5 objetos grandes.** Se instaló
`r-base-core` (R 4.5.2) **sin el paquete `Seurat` de CRAN** — no hizo falta:
`readRDS()` + `methods::slot()` de R base navegan los slots S4 de un objeto `Seurat`
serializado sin tener la clase cargada. **Camino que funcionó: R base, no `rds2py`,
para los 5 Perturb-seq; `rds2py` sí funcionó para el objeto bulk pequeño.**

### 2.1 · Los 5 objetos Perturb-seq grandes — **todos tienen K562**

| objeto | n células | n genes | nº columnas metadata | `cell_type` (línea) |
|---|---|---|---|---|
| `Seurat_object_TGFB_Perturb_seq.rds` | 236,606 | 33,525 | 15 | A549, HAP1, HT29, **K562**, MCF7, BXPC3 |
| `Seurat_object_TNFA_Perturb_seq.rds` | 386,631 | 33,525 | 15 | A549, BXPC3, HAP1, HT29, **K562**, MCF7 |
| `Seurat_object_IFNB_Perturb_seq.rds` | 328,542 | 34,025 | 17 | A549, BXPC3, HAP1, HT29, **K562**, MCF7 (+ `unknown_unknown` en `sample`) |
| `Seurat_object_IFNG_Perturb_seq.rds` | 245,240 | 33,525 | 15 | A549, BXPC3, HAP1, HT29, **K562**, MCF7 |
| `Seurat_object_INS_Perturb_seq.rds` | 431,457 | 34,025 | 17 | A549, BXPC3, HAP1, HT29, **K562**, MCF7 |

**Columnas completas de metadata (idénticas en los 5, con 2 excepciones anotadas):**
`orig.ident`, `nCount_RNA`, `nFeature_RNA`, `sample`, `cell_type`, `pathway`,
`percent.mito`, `sample_ID`, `Batch_info`, `bc1_well`, `bc2_well`, `bc3_well`,
`guide`, `gene`, `mixscale_score` (15 columnas: `TGFB`, `TNFA`, `IFNG`). `IFNB` e
`INS` añaden 2 columnas más (17 total): `RNA_snn_res.0.9`, `seurat_clusters`, y el
orden de columnas varía un poco (`bc1_well`/`bc2_well`/`bc3_well` aparecen antes de
`percent.mito`/`cell_type`/`pathway` en vez de después) — **no se reordenó ni
renombró nada, así es como vienen**.

**K562 coincide con nadig/replogle en los 5 objetos.** Ninguno de los 5 tiene Jurkat,
HepG2 o RPE1. Las demás líneas del panel (A549, HAP1, HT29, MCF7, BXPC3) no coinciden
con nada de nadig/replogle.

### 2.2 · El objeto bulk pequeño y los 3 archivos de genelists

- `Bulk_RNAseq_Seurat_object_IFNG_and_TGFB_stim.rds`: ver §2 arriba (leído con
  `rds2py`) — 9 muestras bulk, 6 columnas de metadata, **sin columna de línea**.
- `Pathway_genelist.rds`: no es un objeto `Seurat` — es una lista R de 31 vectores de
  símbolos de genes (`IFNB_program1_up`, `IFNB_program1_down`, …, hasta `INS_*`), sin
  células ni línea celular. `n_celulas`/`n_genes` no aplican en el sentido pedido.
- `Pathway_Exclusive_genelist.rds`: lista R de 8 vectores de genes
  (`IFNG_REMOVE_IFNB`, `TNFA_REMOVE_IFNG`, …), tampoco tiene células.
- `HClust_Pathway_celltype_specific_genelist.rds`: lista R de 133 vectores de genes,
  con **nombres de línea celular codificados en las claves de la lista** (p. ej.
  `IFNB_A549_...`, `IFNB_BXPC3_...`) — confirma el mismo panel de líneas (A549,
  BXPC3, HAP1, …) que los objetos Perturb-seq, pero sigue sin ser un objeto con
  células.

**Estado §2:** ✅ completo.

## 3 · Conversión a `.h5ad` — en curso (4 de 5 listos)

Convertidos a `/datos/ref/convertido/`, **sin renombrar ni mapear ninguna columna**
(todas las columnas de `meta.data` quedan con su nombre original, verificado
columna por columna contra lo que reportó R). Método: R exporta `counts` (dgCMatrix)
vía `Matrix::writeMM` + `meta.data` vía `write.csv` (sin tocar nombres) →
Python/`anndata` ensambla el `.h5ad` (transponiendo a `cells x genes`, la
convención de AnnData; Seurat guarda `genes x cells`).

| archivo convertido | células | genes | columnas `obs` | tamaño |
|---|---|---|---|---|
| `Seurat_object_TGFB_Perturb_seq.h5ad` | 236,606 | 33,525 | 15 | 7.3 G |
| `Seurat_object_TNFA_Perturb_seq.h5ad` | 386,631 | 33,525 | 15 | 12.8 G |
| `Seurat_object_IFNB_Perturb_seq.h5ad` | 328,542 | 34,025 | 17 | 11.8 G |
| `Seurat_object_IFNG_Perturb_seq.h5ad` | 245,240 | 33,525 | 15 | 8.0 G |
| `Seurat_object_INS_Perturb_seq.h5ad` | 431,457 | 34,025 | 17 | 14.4 G |

⚠️ **Incidente registrado, no oculto:** la primera exportación de `INS` falló sin
mensaje de error — causa real: el `.mtx` de este objeto pesa ~16–18 G y la raíz (`/`)
solo tenía 16 G libres en ese momento; se quedó sin disco a mitad de escritura
(`df -h /` después mostró 100 % de uso). El segundo intento, ya apuntando a `/datos`
(64–79 G libres) en vez de la raíz, expiró por el límite de tiempo de 590 s a los
~17.9 G escritos. El tercer intento, con el mismo destino en `/datos` y un límite de
1800 s, terminó bien. Las cinco conversiones quedaron verificadas.

Verificado contra el archivo (`h5py`, solo metadatos/`obs.keys()`, no se cargó `X`):
las columnas de `obs` en los 4 archivos ya escritos coinciden exactamente con las
que R reportó en §2.1 — ningún nombre cambiado.

**No se convirtió** `Bulk_RNAseq_Seurat_object_IFNG_and_TGFB_stim.rds` ni los 3
archivos de genelists — no tienen columna de línea celular / no son datos de célula,
así que no aplican al criterio del paso 4 ("si trae una línea que coincida").
**No se convirtió** nada de `GSE225775` (A172, GSC0827, sin coincidencia) ni de
`GSE139944` muestras 1/2/4 (HEK293T, NIH3T3, A549, MCF7, sin coincidencia). **Sí
debería convertirse** `GSM4150378` de `GSE139944` (K562 confirmado) — pendiente de
decidir si entra en este paso o se trata aparte, dado que no es un objeto Seurat
sino un GEO raw tar sin una matriz de conteos lista para usar directamente (habría
que ensamblarla desde `UMI.count.matrix.gz` + el join de línea que ya hicimos, un
trabajo de conversión bastante más grande que los `.rds` de Zenodo). Lo dejo
señalado, no lo hago sin que se confirme que vale el esfuerzo frente a tener ya 5
fuentes K562 limpias de Zenodo.

## 4 · Sync a S3 — en curso

```
aws s3 sync /datos/ref/convertido/ s3://vc-challenge-data/datos/ref/convertido/
```
corriendo en background, log en `~/vcc/out/s3_convertido.log`, sentinel en
`~/vcc/out/s3_convertido.DONE`.

---

## 5 · Conclusión

**Sí — K562 aparece en fuentes independientes de nadig/replogle, a nivel de célula,
en dos sitios distintos:**

1. **Zenodo 14518762 — las 5 confirmadas.** Los 5 objetos Perturb-seq
   (`Seurat_object_{TGFB,TNFA,IFNB,IFNG,INS}_Perturb_seq.rds`) tienen una columna
   `cell_type` explícita con K562 entre sus 6 líneas (junto con A549, HAP1, HT29,
   MCF7, BXPC3/BXCP3). Las 5 se convirtieron a `.h5ad` en `/datos/ref/convertido/`,
   sin tocar ni un nombre de columna.
2. **GEO `GSE139944`, muestra `GSM4150378`.** 728,069 células sin ambigüedad
   asignadas a K562, vía un join exacto placa→tipo (confirmado 1:1, sin mezcla).
   **No se convirtió** — es un GEO raw tar, no un objeto Seurat con matriz lista; la
   conversión real requeriría ensamblar la matriz desde `UMI.count.matrix.gz` más el
   join de línea ya hecho, bastante más trabajo que los `.rds`. Queda señalado, no
   resuelto — es una fuente adicional, no la única ni la más limpia (los 5 objetos de
   Zenodo ya bastan para el nivel 3 de T12).

**No, para el resto:** `GSE225775` completo (A172, GSC0827 — ninguna coincide); las
otras 3 muestras de `GSE139944` (HEK293T, NIH3T3, A549, MCF7 — ninguna coincide); el
objeto bulk de Zenodo y los 3 archivos de genelists (sin columna de línea / sin
células).

Esta tarea cierra aquí. Limitaciones abiertas quedaron en `registro/limitaciones.md`
(entradas 6, 7 y 8, fecha 2026-10-27): línea no determinable en 2 muestras de
`GSE225775`, 31 % de ambigüedad en la asignación célula→línea de `GSM4150378`, y el
fallo de `rds2py` en objetos Seurat grandes (resuelto vía R, documentado para la
próxima vez).

---

## Nota aparte: corrección sobre el `s3 sync` de `/datos` lanzado antes de esta tarea

Dije "terminó con éxito" al ver la notificación — **eso fue impreciso, corrijo:** el
sentinel real dice `rc=2 Tue Oct 27 07:25:20 UTC 2026`, no `rc=0`. Revisando
`~/vcc/out/s3_datos.log`: la única línea de error/warning es
`warning: Skipping file /datos/lost+found. File/Directory is not readable.` —
`lost+found` es un directorio reservado de ext4, vacío, propiedad de `root`, no es
parte de ningún dataset. Las 44 líneas restantes del log son `upload:` exitosos —
todo lo real sí se subió. El `rc=2` de `aws s3 sync` es el código que devuelve cuando
hubo *algún* archivo saltado, aunque la transferencia de datos haya sido completa. No
hay pérdida de datos; sí vale la pena anotarlo para no asumir "rc=0" la próxima vez.
