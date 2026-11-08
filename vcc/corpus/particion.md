# T3 — Cribado de admisibilidad sobre `corpus/inventario.parquet`

Clasificación de cada celda **(fuente, línea)** de la matriz de incidencia:
**(A)** reproducible/usable tal cual · **(B)** falta el dato (se nombra qué) ·
**(C)** estimando desconocido (se nombra qué no se sabe).

**Rúbrica declarada antes de clasificar** (para que el orden no sea arbitrario):

**Restricciones** (cosas que el propio dato fija, suman punto cada una):
`R1` grano leído, no inferido · `R2` tipo de valor en X leído, no inferido ·
`R3` regla de control exacta y verificada (no ambigua) · `R4` clave de unión limpia e
inyectiva con otra fuente de la misma línea, ya establecida en trabajo previo ·
`R5` línea celular declarada como columna propia del archivo (no tomada del nombre de
archivo) · `R6` eje de genes sin mezcla de espacio de símbolo dentro de la misma
columna · `R7` contexto declarado explícitamente en el archivo de control (no inferido
del nombre de archivo).

**Decisiones libres** (cosas que hay que fijar para poder usar el dato, restan orden):
`D1` normalización de X por decidir · `D2` mapeo de nombre de línea a una llave
normalizada, si se quiere cruzar con otra fuente · `D3` elegir construcción
representativa cuando hay más de una por gen (confirmado no neutral en
`gwps_fiabilidad.md`) · `D4` sin células perturbadas reales, solo metadata/promesa ·
`D5` espacio de símbolo mixto dentro de una misma columna, hay que decidir cuál
interpretar fila por fila · `D6` línea celular no declarada en absoluto en el archivo.

Orden: por **R/D** descendente (más restricción por decisión libre = más barato de usar).

| fuente(s) | línea | clase | R (lista) | D (lista) | #R | #D | R/D |
|---|---|---|---|---|---|---|---|
| `vcc2025_{train,validation,test}_pert_counts` | (no declarada) | A | R1*,R2*,R3 | D6 | 3 | 1 | 3.00 |
| `replogle_{k562_essential,rpe1}_sc` | k562_essential, rpe1 | A | R1,R2,R3,R4,R6 | D1,D2 | 5 | 2 | 2.50 |
| `replogle_{k562_essential,rpe1}_bulk` | k562_essential, rpe1 | A | R1,R2,R3,R4,R6 | D2,D_agregación-ajena¹ | 5 | 2 | 2.50 |
| `nadig_{hepg2,jurkat}` | hepg2, jurkat | A | R1,R2,R3,R6 | D1,D2 | 4 | 2 | 2.00 |
| `zenodo_{IFNB,IFNG,INS,TGFB,TNFA}` × K562 | K562 | A | R1,R2,R3,R6 | D1,D2 | 4 | 2 | 2.00 |
| `vcc2025_{train,validation,test}` | (no declarada) | A | R1,R2,R3 | D1,D6 | 3 | 2 | 1.50 |
| `zenodo_{IFNB,IFNG,INS,TGFB,TNFA}` × {A549,BXPC3,HAP1,HT29,MCF7} | esas 5 | A | R1,R2,R3,R6 | D1,D2² | 4 | 2 | 2.00 |
| `replogle_{k562_gwps}_sc` | k562_gwps | A | R1,R2,R3,R4,R6 | D1,D2,D3 | 5 | 3 | 1.67 |
| `replogle_{k562_gwps}_bulk` | k562_gwps | A | R1,R2,R3,R4,R6 | D2,D3,D_agregación-ajena¹ | 5 | 3 | 1.67 |
| `reto_{val,final}_context_{A..F}` | (no declarada) | A³ | R1,R2,R3,R7 | D4,D2(no aplica⁴) | 4 | 1 | 4.00 |
| `orion_{hct116,hek293t}` | hct116, hek293t | **B** | R3,R_clave-ya-pareada | D5,D_sin-células | 2 | 2 | 1.00 |
| `feng_{genomewide,targeted}` | (no declarada) | **B** | R_clave-ya-pareada | D6,D_sin-n_celulas | 1 | 2 | 0.50 |
| `reto_{val,final}_pert_counts` | (no declarada) | **B** | R6 | D4,D_sin-n_cells⁵ | 1 | 2 | 0.50 |
| `feng_targeted_perline` | 16 códigos crípticos | **C** | R6 | D6(identidad irresoluble) | 1 | 1 | 1.00⁶ |

¹ El pseudobulk de los bulk de Replogle llegó ya agregado por los autores (media sobre
una población de control que ya se intentó reconstruir en `replogle_controlexpr*.md` —
parcialmente, no al 100%). Es una decisión que no tomamos nosotros y que no se puede
auditar completo sin el `_sc`.
² Sin fuente pareja conocida en *este* inventario para cruzar estas 5 líneas — la
decisión de normalización no tiene destino todavía, se nombra igual porque podría
aparecer uno.
³ Clase A **acotada al alcance real del archivo**: el pool de control en sí es completo
y reproducible. Lo que falta (las células perturbadas) es un defecto de dato distinto,
clasificado aparte en la fila de `pert_counts`.
⁴ El contexto SÍ viene declarado (`R7`), así que no hace falta mapearlo — `D2` no aplica
a esta fila, se deja explícito en vez de omitirlo.
⁵ A diferencia de `vcc2025_*_pert_counts`, el `pert_counts.csv` del reto **no tiene
columna `n_cells`** — ni siquiera el recuento prometido.
⁶ Con solo 1 restricción y 1 decisión el cociente es igual al de `orion`, pero la clase
es **C**, no B: no es que falte un dato que se pueda ir a buscar, es que la identidad de
la línea no es recuperable de ningún archivo que tengamos — ninguna cantidad de trabajo
adicional sobre estos mismos archivos lo resuelve, se necesitaría una fuente externa (el
mapa código→línea del estudio de Feng), que no tenemos.

## Qué falta, nombrado (clase B)

- **`orion_{hct116,hek293t}`**: faltan las células originales. Solo tenemos el resumen
  `(target, feature)` ya contrastado (LFC, p-valor, Ueffect) — útil para comparar
  magnitudes publicadas, no para recomputar ningún test propio sobre estos datos.
- **`feng_{genomewide,targeted}`**: falta la línea celular (no declarada en ningún
  archivo) y falta cualquier recuento de células, a cualquier nivel.
- **`reto_{val,final}_pert_counts`**: faltan las células perturbadas reales — el .h5ad
  de control (`context_{A..F}`) solo trae el pool NTC; las 138,400 células de verdad
  por contexto (`ground_truth_cells` del `manifest.json`) no están en ningún archivo
  que tengamos. El csv ni siquiera promete un `n_cells` — solo lista los 300 genes.

## Qué no se sabe, nombrado (clase C)

- **`feng_targeted_perline`**: no se sabe a qué línea celular estándar corresponde cada
  uno de los 16 códigos (`eipl_1`, `fiaj_1`, `iudw_1`, …). No es una normalización
  pendiente (eso sería una decisión, clase A/B) — es una identidad que no está en
  ningún archivo de este corpus. Sin una tabla externa del estudio de Feng, esto se
  queda cerrado.

## Nota transversal (no es una celda, es un hallazgo del cribado completo)

Bajo cadena cruda, **ninguna línea aparece en ≥2 fuentes independientes** (ver query 3
de `t3_inventario.md`) — ni siquiera K562, porque `"k562_essential"`/`"k562_gwps"`
(Replogle) y `"K562"` (Zenodo) son cadenas distintas y el primitivo, por diseño, no las
junta. La fila de `zenodo_*` × K562 en esta tabla se marca (A) porque el dato en sí es
usable; el cruce con Replogle que de verdad importa para T12 nivel 3 requiere el mapa
de normalización **externo** que el primitivo deliberadamente no aplica — ver T3 §3 del
reporte principal.
