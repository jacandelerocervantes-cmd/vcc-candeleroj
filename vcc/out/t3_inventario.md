# T3 — La matriz de incidencia, y el cribado de admisibilidad sobre ella

Por qué: el error de "cobertura contra el panel equivocado" ya ocurrió dos veces (el
eje de 18,533 genes del envío tratado como filtro; val/final resultando no compartir
panel, 129 de 300). El arreglo no es revisar cifras una por una: es que el **primitivo**
sea la matriz de incidencia `(fuente, archivo, línea_celular, contexto, perturbación,
construcción)` al grano fino, sin colapsar ni normalizar nada, y que toda cifra de
cobertura sea una **consulta** sobre ella.

## Declaración previa (escrita antes de correr las consultas 1–7)

- **Si `target_crudo` de una fuente es una cadena compuesta** (caso de los bulk de
  Replogle, donde no existe una columna "gen blanco" separada de `gene_transcript`):
  las consultas de cobertura por símbolo de gen contra esa fuente **darán 0 por
  diseño**, no porque no haya coincidencia biológica, sino porque el primitivo no
  normaliza. Eso no es un fallo de la consulta — es la prueba de que el primitivo hace
  lo que se le pidió (no inferir ni colapsar).
- **Si dos fuentes nombran la misma línea biológica con cadenas distintas**
  (`"k562_essential"`/`"k562_gwps"` en Replogle vs. `"K562"` en Zenodo): la consulta 3
  (líneas en ≥2 fuentes independientes) **dará 0** bajo cadena cruda, aunque sepamos
  por trabajo previo (`geo_zenodo.md`) que K562 sí es compartido. Eso separa
  explícitamente lo que el primitivo establece (nada, por diseño) de lo que un mapa de
  normalización **externo y declarado aparte** establecería.
- **Si una consulta pide un grano distinto al de la tabla** (p. ej. CPM por gen y por
  contexto, que es grano `(contexto, gen)`, no `(fuente, línea, contexto, perturbación,
  construcción)`): esa consulta **no sale de `inventario.parquet`**, y eso es un defecto
  de diseño a nombrar, con la columna/tabla que faltaría — no una consulta mal escrita.

## Entregable 1 — `corpus/inventario.parquet`

**72,863 filas**, una por `(fuente, archivo, línea_celular_cruda, contexto,
target_crudo, construcción, es_control)`, sobre las 30 fuentes de `/datos/ref`
(Replogle ×3 líneas ×2 granos, Nadig ×2, Orion ×2, Feng ×3 archivos, los 5 convertidos
de Zenodo), `/datos/reto` (val y final: 3+3 archivos de contexto + 2 `pert_counts.csv`),
y `/datos/vcc2025` (train/validation/test: 3 `.h5ad` + 3 `pert_counts_*.csv`).

**Reglas aplicadas sin excepción:**
- No se colapsó a símbolo de gen en ninguna fuente: `construccion` conserva
  `gene_transcript` (Replogle/Nadig) o `guide`/`guide_id` (Zenodo/vcc2025) tal cual.
- No se normalizó ningún nombre de línea: `K562` (Zenodo, confirmado idéntico en los 5
  archivos) y `k562_essential`/`k562_gwps` (Replogle, tomado del nombre de archivo
  porque **no existe como columna** en esos archivos) se guardaron como cadenas
  distintas, a propósito.
- `tipo_valor_en_X` se leyó de una muestra real de `X` (no se infirió de la extensión
  ni del nombre): dense → primeras 20,000 filas; CSR → primeros 5,000,000 de `data`.
- El contexto A–F de `/datos/reto` se tomó de la columna `context` **dentro** del
  propio archivo de control (y se verificó contra `manifest.json` — coinciden los
  tres: nombre de archivo, columna interna, y manifiesto).

**Hallazgo que no se buscaba y que cambia cómo se debe usar Replogle de aquí en
adelante**: el filtro `core_control==True` que reportes anteriores de esta serie usaron
para `k562_essential`/`rpe1` **no es válido para `k562_gwps`** — tiene 71 filas de
pseudo-grupo de control con `core_control==False` (ver `gwps_fiabilidad.md`, ya
corregido ahí). El primitivo usa `"non-targeting" in gene_transcript` (subcadena) en
las tres líneas, consistente.

**Hallazgo sobre `/datos/reto`**: los 6 archivos `context_{A..F}.h5ad` **no contienen
ninguna célula perturbada** — `target_gene` tiene una sola categoría, `"non-targeting"`,
en los 6. Son el pool de control, nada más. Las 138,400 células de verdad por contexto
que promete `manifest.json` (`ground_truth_cells`) no están en ningún archivo que
tengamos.

## Entregable 2 — `corpus/particion.md`, el cribado de admisibilidad

33 celdas `(fuente agrupada, línea)`, clasificadas A/B/C con restricciones y decisiones
libres nombradas, ordenadas por R/D. Resumen: **3 celdas en B** (`orion_*`:
faltan las células, solo queda el resumen ya contrastado; `feng_{genomewide,targeted}`:
falta línea celular y cualquier recuento de células; `reto_*_pert_counts`: faltan las
células perturbadas reales, el csv ni siquiera trae `n_cells`), **1 celda en C**
(`feng_targeted_perline`: 16 códigos de línea crípticos — `eipl_1`, `fiaj_1`, … — sin
mapa a un nombre estándar en ningún archivo del corpus), **el resto en A** (con
limitaciones nombradas caso por caso, ninguna oculta). Detalle completo, con la rúbrica
declarada antes de clasificar, en `corpus/particion.md`.

## Entregable 3 — las siete consultas (test del primitivo)

### 1 · Cobertura del panel FINAL (300 genes) por fuente y por línea

Consulta pura sobre `inventario.parquet` (el panel de 300 sale de las propias filas de
`reto_final_pert_counts`, no se releyó el csv aparte). Top filas por cobertura:

| fuente | línea | n targets distintos | cobertura panel final (de 300) |
|---|---|---|---|
| `reto_final_pert_counts` | — | 300 | 300 (por construcción) |
| `replogle_k562_gwps_sc` | k562_gwps | 9,866 | **277** |
| `feng_genomewide` | — | 6,673 | 185 |
| `orion_hct116` | hct116 | 8,987 | 148 |
| `reto_val_pert_counts` | — | 300 | 129 |
| `orion_hek293t` | hek293t | 4,297 | 14 |
| `vcc2025_train` | — | 150 | 12 |
| `replogle_k562_essential_sc` | k562_essential | 2,057 | **0** |
| `replogle_k562_gwps_bulk` | k562_gwps | 10,673 | **0** |
| `replogle_rpe1_sc` | rpe1 | 2,393 | **0** |

**Los ceros de Replogle no son un error de la consulta.** `k562_essential_sc`/`rpe1_sc`
dan 0 porque son pantallas de genes **esenciales** y el panel del reto es no-esencial
(0 solapamiento real, ya confirmado en `nreal_predictores.md` — el primitivo reproduce
el mismo cero de forma independiente). `k562_gwps_bulk` da 0 **por la cadena
compuesta** (`target_crudo` ahí es el `gene_transcript` completo, no el símbolo) — ésta
sí es una limitación del primitivo para cruces por símbolo contra fuentes pseudobulk de
Replogle, nombrada en la declaración previa, no descubierta a mitad de consulta.
**El 277 de `replogle_k562_gwps_sc` coincide exacto con el número ya establecido a mano
en `nreal_predictores.md`** — el primitivo valida contra un resultado independiente ya
conocido.

### 2 · Cobertura del panel de VALIDACIÓN (300) + intersección

Misma consulta, panel de validación. `replogle_k562_gwps_sc`: 272 (coincide con
`gwps_resolucion.md`). **Intersección panel val ∩ panel final = 129** — confirmado
exacto, leído de las propias filas de `reto_val_pert_counts` y `reto_final_pert_counts`
dentro del parquet, sin volver a los csv originales.

### 3 · Líneas celulares en ≥2 fuentes independientes

Fuentes independientes agrupadas: `replogle` (3 líneas, sc+bulk cuentan como 1 fuente),
`nadig`, `zenodo` (las 5 rutas de pathway son el mismo estudio/depósito, 1 fuente),
`orion`, `feng`.

**Bajo cadena cruda, sin normalizar: 0 líneas en ≥2 fuentes independientes.**
Ni K562 — `"k562_essential"`/`"k562_gwps"` (Replogle) y `"K562"` (Zenodo) son cadenas
distintas, y el primitivo, por diseño, no las junta. Esto **no contradice**
`geo_zenodo.md` (que sí concluyó K562 compartido): ese reporte usó una identificación
humana explícita de la línea, no una igualdad de cadenas. La consulta 3 aísla
exactamente esa diferencia — es el resultado, no un defecto. **Con un mapa de
normalización externo y declarado** (no inferido de la forma de la cadena), K562
aparecería en 2 fuentes independientes (`replogle` + `zenodo`) — esa es la columna que
T12 nivel 3 necesita, y el primitivo la deja lista para construirse aparte, sin
mezclarla con los datos crudos.

### 4 · Genes con >5 CPM por contexto — **no sale de la tabla, defecto de diseño nombrado**

Esta consulta tiene grano `(contexto, gen)` — una cifra de expresión por gen medido,
por contexto. `inventario.parquet` tiene grano `(fuente, línea, contexto, perturbación,
construcción)`: una fila por perturbación, no por gen del eje. **No se puede contestar
desde la tabla tal como está diseñada**, y no se intentó forzarla leyendo `X` de los 6
archivos de contexto (eso excede lo que la tarea autorizó tocar — solo justifica tocar
`X` el punto 6, integridad de Zenodo, y la lectura ya hecha de `tipo_valor_en_X`).
**Columna/tabla que falta**: una tabla aparte, `corpus/cpm_por_gen_contexto.parquet`,
grano `(contexto, gen, cpm_medio_sobre_el_pool_de_control)`, construida leyendo `X` de
`context_{A..F}.h5ad` una vez — most fuera del alcance de esta entrega.

### 5 · Genes del panel final con datos a nivel de CÉLULA en alguna fuente

Filtrando `inventario.parquet` a `grano == "celula_individual"` y `es_control==False`:

**277 de los 300 genes del panel final tienen datos a nivel de célula en alguna
fuente** — y los 277 vienen enteramente de `replogle_k562_gwps_sc` (coincide con la
cobertura de la consulta 1, como debía). El resto de fuentes con grano celular aportan
muy poco de forma independiente: `vcc2025_train` 12, `vcc2025_test` 11,
`vcc2025_validation` 5, `zenodo_TNFA` 5, `zenodo_IFNB` 2, `zenodo_IFNG` 1 — sin ninguno
de estos sumar genes que `k562_gwps_sc` no tuviera ya (no se verificó la unión exacta
de conjuntos más allá del total, queda como chequeo menor pendiente, no bloqueante).

### 6 · Integridad de X en los 5 `.h5ad` convertidos de Zenodo

| archivo | forma (`.h5ad`) | forma (R, `geo_zenodo.md`) | coincide | nnz (`data`) | `indptr` consistente |
|---|---|---|---|---|---|
| IFNB | 328,542 × 34,025 | 328,542 × 34,025 | ✓ | 980,979,021 | ✓ |
| IFNG | 245,240 × 33,525 | 245,240 × 33,525 | ✓ | 667,070,362 | ✓ |
| **INS** | **431,457 × 34,025** | **431,457 × 34,025** | **✓** | **1,280,241,271** | **✓** |
| TGFB | 236,606 × 33,525 | 236,606 × 33,525 | ✓ | 605,806,920 | ✓ |
| TNFA | 386,631 × 33,525 | 386,631 × 33,525 | ✓ | 1,066,935,493 | ✓ |

Las cinco forman coinciden exacto con lo que R reportó, y la consistencia interna
(`len(indptr) == n_filas+1`, `indptr[-1] == len(data) == len(indices)`) se cumple en
las cinco. **`INS` — el que falló dos veces (disco lleno, luego timeout) antes de la
tercera corrida — pasa íntegro**: no quedó un `.mtx` truncado disfrazado de archivo
completo. Esto cierra, con evidencia directa sobre el archivo (no solo sobre el log de
la corrida), la preocupación del LEDGER §9 sobre `INS`.

### 7 · BXPC3 vs BXCP3, y la cadena exacta de K562

Leído directamente de `obs/cell_type` en los 5 `.h5ad` (no del reporte de conversión):
**las 5 fuentes dicen `BXPC3`** (nunca `BXCP3`) — categorías idénticas en los 5:
`['A549', 'BXPC3', 'HAP1', 'HT29', 'K562', 'MCF7']`. **`K562` también es idéntico en
las 5**, sin variante de mayúsculas/espaciado. **Conclusión: `BXCP3` fue un error de
transcripción en `geo_zenodo.md` (línea 243 del propio reporte), no una segunda
categoría real en el dato.** No hay nada que corregir en los archivos — sí en la
redacción de `geo_zenodo.md`, que se deja señalada aquí en vez de editada ahí (no se
reescribe un reporte cerrado; se corrige la cadena en el registro de limitaciones si se
decide después).

## 20 filas trabajadas, cubriendo fuentes y granos distintos

| fuente | línea | contexto | target_crudo | construcción | control | n_células | n_gem_groups | grano | tipo_valor_X | dtype X | espacio símbolo | n genes eje |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| replogle_k562_essential_sc | k562_essential | — | AAAS | 2_AAAS_P1P2_ENSG00000094914 | False | 212 | 47 | celula_individual | conteo_entero | float32 | ambos | 8563 |
| replogle_k562_essential_sc | k562_essential | — | non-targeting | 10753_non-targeting_... | True | 61 | 35 | celula_individual | conteo_entero | float32 | ambos | 8563 |
| replogle_k562_gwps_bulk | k562_gwps | — | 0_A1BG_P1_ENSG00000121410 | (idéntico) | False | 237 | — | pseudobulk | normalizado/log | float32 | ambos | 8248 |
| nadig_hepg2 | hepg2 | — | AAAS | 2_AAAS_P1P2_ENSG00000094914 | False | 39 | 27 | celula_individual | conteo_entero | float32 | ambos | 9624 |
| zenodo_TNFA | A549 | — | NT | NTg1 | True | 152 | 2 | celula_individual | conteo_entero | int64 | simbolo | 33525 |
| zenodo_IFNG | **K562** | — | ATF3 | ATF3g1 | False | 484 | 2 | celula_individual | conteo_entero | int64 | simbolo | 33525 |
| orion_hct116 | hct116 | — | A3GALT2 | — | False | 70 | — | resumen_de_pares(target,feature) | — | — | **ambos (mixto en la misma columna)** | 20348 |
| orion_hek293t | hek293t | — | A3GALT2 | — | False | 67 | — | resumen_de_pares(target,feature) | — | — | ambos | 18518 |
| feng_genomewide | — | — | A1BG | — | False | — | — | resumen_LFC_por_gen | — | — | ambos | 6472 |
| feng_targeted_perline | **eipl_1** | — | AAMP | — | False | — | — | resumen_LFC_por_gen | — | — | ambos | 6520 |
| reto_val_context_A | — | **A** | non-targeting | — | True | 18400 | 46 | celula_individual | conteo_entero | float32 | simbolo | 18533 |
| reto_final_context_D | — | **D** | non-targeting | — | True | 18400 | 46 | celula_individual | conteo_entero | float32 | simbolo | 18533 |
| reto_val_pert_counts | — | aplica a A/B/C | ABCD1 | — | False | — (sin n_cells) | — | — | — | — | simbolo | — |
| vcc2025_train | — | — | ACAT2 | ACAT2_P1P2_A\|ACAT2_P1P2_B | False | 1765 | 48 | celula_individual | conteo_entero | float32 | simbolo | 18080 |
| vcc2025_train | — | — | non-targeting | non-targeting_00005\|...01689 | True | 1269 | 48 | celula_individual | conteo_entero | float32 | simbolo | 18080 |
| vcc2025_test_pert_counts | — | — | GRB10 | — | False | 4074 (del csv) | — | — | — | — | simbolo | — |
| replogle_rpe1_sc | rpe1 | — | AAAS | 2_AAAS_P1P2_ENSG00000094914 | False | 86 | 45 | celula_individual | conteo_entero | float32 | ambos | 8749 |
| replogle_rpe1_bulk | rpe1 | — | 10755_non-targeting_... | (idéntico) | True | — | — | pseudobulk | normalizado/log | float32 | ambos | 8749 |
| nadig_jurkat | jurkat | — | non-targeting | 10755_non-targeting_... | True | 57 | 37 | celula_individual | conteo_entero | float32 | ambos | 8882 |
| vcc2025_validation_pert_counts | — | — | SH3BP4 | — | False | 2925 (del csv) | — | — | — | — | simbolo | — |

Nota forzada por esta tabla, no buscada: `orion_hct116`/`hek293t` tienen `feature` con
**espacio de símbolo mezclado dentro de la misma columna** — la fila de `A3GALT2`
mostrada arriba convive, para el mismo `target`, con filas `feature=ENSG00000291215`
(Ensembl) junto a `feature=LINC01409` (símbolo). No se decidió cuál forma corresponde a
cuál fila — **D5** de la rúbrica de `particion.md`, visible aquí con un ejemplo real,
no solo como categoría abstracta.

## PARA

**Qué decide:**

1. **El error de "panel equivocado" deja de ser posible estructuralmente**, no solo en
   esta entrega: toda cifra de cobertura reportada arriba es trazable a una consulta
   reproducible sobre `inventario.parquet`, y las dos cifras ya conocidas de forma
   independiente (277 en `k562_gwps_sc`, intersección 129 val∩final) salieron
   **idénticas** a lo ya establecido a mano — el primitivo se valida contra sí mismo.
2. **La consulta 3 deja una distinción explícita y permanente**: "mismo string" y
   "misma línea biológica" no son la misma pregunta. El mapa que junta
   `k562_essential`/`k562_gwps`/`K562` tiene que vivir en una tabla aparte, declarada,
   nunca en el primitivo — exactamente lo que TRASPASO §4 deja abierto para la regla de
   colapso.
3. **El cribado de `particion.md` ordena el trabajo**: 3 celdas B (`orion_*`,
   `feng_{genomewide,targeted}`, `reto_*_pert_counts`) y 1 celda C
   (`feng_targeted_perline`) no son utilizables sin conseguir un dato que no tenemos o
   sin aceptar que no se sabrá — el resto (29 celdas) es A, con limitaciones nombradas,
   no ocultas.
4. **Zenodo queda íntegro** (query 6) y **BXPC3/K562 quedan sin ambigüedad** (query 7)
   — dos preocupaciones del LEDGER §9 cerradas con evidencia directa del archivo.

**Qué queda abierto, y con qué depende:**

- **Query 4** (CPM>5 por contexto) necesita una tabla aparte de grano `(contexto, gen)`
  — `corpus/cpm_por_gen_contexto.parquet` — que implica leer `X` de los 6 archivos de
  contexto completos, fuera del alcance autorizado aquí.
- **El mapa de normalización de línea celular** (K562 y cualquier otra equivalencia
  entre fuentes) queda señalado como necesario pero no construido — es la pieza que
  falta para que la consulta 3 refleje la biología y no solo la cadena.
- **La reproducibilidad de los bulk de Replogle** sigue dependiendo de qué tan bien se
  reconstruya la convención exacta de agregación de los autores (trabajo parcial en
  `replogle_controlexpr*.md`) — el primitivo registra que son pseudobulk ya agregado,
  no que esa agregación esté auditada al 100%.
- Con el primitivo en su sitio, **el paso 5 del PLAN (A completo con el procedimiento
  validado) puede empezar sobre una base donde la cobertura ya no es una fuente de
  error** — pero el procedimiento mismo (qué convención de X3 se valida contra el
  benchmark) sigue siendo una pregunta aparte, no resuelta por este inventario.
