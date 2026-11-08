# VCC 2025 como superficie: inventario

## 🔴 0. Que es esta superficie y que no

`/datos/vcc2025/` tiene **train, validation y test CON VERDAD**. Es la unica superficie del proyecto donde se puede ajustar contra verdad **sin ser transductivo**: la verdad de 2026 no existe, y la de Jiang se usa como contexto, no como panel.

Es el eje **PANEL** de los dos ortogonales. Jiang es el eje **CONTEXTO**.

🔴 **Lo que esta superficie NO puede dar** (§5, declarado aqui arriba tambien): es **UN SOLO contexto** -- una sola linea celular, un solo experimento. No mide dominancia entre contextos y no da rango entre contextos. Cualquier cifra de 2025 es de un punto, no de un rango.

## 1. Los tres splits: formas, conteos y sha256

| #  | split | GiB | celulas | genes | perturbaciones | celulas `non-targeting` | nnz total | nnz mediana/celula | UMI mediana/celula (del csv) | sha256 (16) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | **train** | 14.42 | 221,273 | 18,080 | 150 | 38,176 | 1932.6 M | 8761 | 53,994 | a09977104fefb622... |
| 2 | **validation** | 6.45 | 98,927 | 18,080 | 50 | 38,176 | 864.8 M | 8770 | 54,447 | 376f0bab27d9f22e... |
| 3 | **test** | 11.13 | 170,846 | 18,080 | 100 | 38,176 | 1491.6 M | 8761 | 54,285 | ba5ce66729e4c8e5... |

**sha256 completos:**

| split | bytes | sha256 |
|---|---|---|
| train | 15,482,497,461 | `a09977104fefb622368ca74b50c9d3c1e891733e6c83db07acfca49b0219c02b` |
| validation | 6,928,967,541 | `376f0bab27d9f22e00dfd888668d3602006d035edc182721c378f4802c77bfe1` |
| test | 11,950,739,168 | `ba5ce66729e4c8e50b6f147c7d73d3f9f28f76c91fa43ff33e3f6bc2ed831f03` |

### 1.1 Control de consistencia: el csv contra el `obs`

| #  | split | perts en el csv | perts en el obs | `n_cells` del csv coincide con el obs | # lotes | # `guide_id` | columnas de obs | duracion |
|---|---|---|---|---|---|---|---|---|
| 1 | **train** | 150 | 150 | si | 48 | 189 | _index, batch, guide_id, target_gene | 2.0 min (118 s) |
| 2 | **validation** | 50 | 50 | si | 48 | 85 | _index, batch, guide_id, target_gene | 37.7 s |
| 3 | **test** | 100 | 100 | si | 48 | 136 | _index, batch, guide_id, target_gene | 1.5 min (91 s) |

### 1.2 🔴 Potencia por perturbacion: no todas llegan a 400 celulas

| #  | split | min `n_cells` | mediana | max | **perts con < 400 celulas** | rango de UMI mediana |
|---|---|---|---|---|---|---|
| 1 | **train** | 33 | 1045 | 4760 | **24** de 150 | 38,150 - 60,202 |
| 2 | **validation** | 161 | 1090 | 2925 | **6** de 50 | 38,650 - 63,176 |
| 3 | **test** | 67 | 1142 | 4074 | **15** de 100 | 39,764 - 66,959 |

El reto 2026 fija **400 celulas por perturbacion** (`manifest.json: cells_per_pert`). En 2025 eso **no se cumple para todas**, asi que cualquier construccion que pida 400 tiene que declarar que hace con las que no llegan.

#### Las 24 perturbaciones con < 400 celulas en `train`

| #  | pert | n_cells | UMI mediana/celula |
|---|---|---|---|
| 1 | `ATP6V0C` | 33 | 55,142 |
| 2 | `RNF20` | 47 | 60,202 |
| 3 | `OXA1L` | 63 | 46,607 |
| 4 | `DNAJA3` | 65 | 48,578 |
| 5 | `BRD9` | 67 | 50,153 |
| 6 | `TAF13` | 138 | 47,104 |
| 7 | `FDPS` | 154 | 54,408 |
| 8 | `HSBP1` | 166 | 52,200 |
| 9 | `SLC39A6` | 178 | 44,946 |
| 10 | `SALL4` | 180 | 48,238 |
| 11 | `KIF20A` | 183 | 56,345 |
| 12 | `SDC1` | 215 | 51,225 |
| 13 | `EPB41L2` | 219 | 56,419 |
| 14 | `SHPRH` | 241 | 52,685 |
| 15 | `EPHB4` | 241 | 49,905 |
| 16 | `ATP6V0B` | 245 | 55,734 |
| 17 | `ETV4` | 263 | 55,141 |
| 18 | `SEC62` | 286 | 50,149 |
| 19 | `KDM1A` | 352 | 40,262 |
| 20 | `NREP` | 352 | 59,080 |
| 21 | `PMEL` | 356 | 54,419 |
| 22 | `PRDM14` | 359 | 38,150 |
| 23 | `SUPT4H1` | 364 | 42,893 |
| 24 | `SNCA` | 386 | 53,942 |

#### Las 6 perturbaciones con < 400 celulas en `validation`

| #  | pert | n_cells | UMI mediana/celula |
|---|---|---|---|
| 1 | `MAX` | 161 | 38,650 |
| 2 | `ZNF598` | 255 | 54,046 |
| 3 | `WAC` | 267 | 60,363 |
| 4 | `FUBP1` | 288 | 55,334 |
| 5 | `EIF3H` | 323 | 53,504 |
| 6 | `SUPV3L1` | 339 | 51,024 |

#### Las 15 perturbaciones con < 400 celulas en `test`

| #  | pert | n_cells | UMI mediana/celula |
|---|---|---|---|
| 1 | `NELFE` | 67 | 45,739 |
| 2 | `ZRSR2` | 76 | 51,568 |
| 3 | `NCAPH2` | 98 | 52,106 |
| 4 | `PSMG1` | 139 | 42,147 |
| 5 | `RNF40` | 150 | 58,534 |
| 6 | `OSTC` | 159 | 49,117 |
| 7 | `INO80B` | 170 | 44,988 |
| 8 | `ELAC2` | 184 | 47,454 |
| 9 | `ZFP62` | 204 | 47,560 |
| 10 | `TAF11` | 257 | 46,913 |
| 11 | `DHX9` | 273 | 54,599 |
| 12 | `ADNP` | 291 | 46,586 |
| 13 | `BECN1` | 349 | 55,246 |
| 14 | `BRPF1` | 366 | 53,649 |
| 15 | `NSDHL` | 374 | 57,456 |

## 2. 🔴 Los genes pareados con el panel 2026

| #  | cantidad | valor |
|---|---|---|
| 1 | genes medidos en 2025 | 18,080 |
| 2 | genes medidos en 2026 | 18,533 |
| 3 | **genes medidos en los dos** | **18,077** |
| 4 | genes de 2026 que faltan en 2025 | 456: `ABCC6`, `ABHD17A`, `ABITRAM`, `AC005726.1`, `AC007244.1`, `AC023490.5`, `AC048338.2`, `AC084121.12`, `AC084121.13`, `AC084121.9` |
| 5 | perturbaciones 2025 (union de los 3 splits) | 300 |
| 6 | perturbaciones 2026 `final` / `val` / union / interseccion | 300 / 300 / 471 / 129 |

🔴 **El encargo cita '97 genes pareados'. Se calcularon las siete definiciones de 'pareado' que esta superficie admite y se reporta lo medido en cada una.**

| #  | definicion | que es | **n medido** | ¿97? |
|---|---|---|---|---|
| 1 | `d1_perts2025_x_perts2026_final` | perts 2025 ∩ perts 2026 `final` | 28 |  |
| 2 | `d2_perts2025_x_perts2026_val` | perts 2025 ∩ perts 2026 `val` | 25 |  |
| 3 | `d3_perts2025_x_union_perts2026` | perts 2025 ∩ (union de perts 2026) | 37 |  |
| 4 | `d4_perts2025_x_interseccion_perts2026` | perts 2025 ∩ (interseccion de perts 2026) | 16 |  |
| 5 | `d5_perts2025_en_el_espacio_de_genes_2026` | perts 2025 que 2026 MIDE como gen | 300 |  |
| 6 | `d6_perts2026_union_en_el_espacio_de_genes_2025` | perts 2026 que 2025 MIDE como gen | 471 |  |
| 7 | `d7_genes_medidos_en_los_dos` | genes medidos en los dos paneles | 18077 |  |

🔴 **La cifra 97 NO sale de ninguna de las siete definiciones.** Se deja asi, medido y declarado, sin forzar ninguna definicion para que de 97. La que mas se le parece en espiritu -- perturbaciones de 2025 que tambien se perturban en 2026 -- da **37** contra la union de los dos paneles 2026.

### 2.1 Las listas exactas

**`d3_perts2025_x_union_perts2026`** -- perts 2025 ∩ (union de perts 2026), **37** genes:

`ACLY`, `ADNP`, `AKT2`, `ANKZF1`, `ASF1A`, `BRPF1`, `DCK`, `FUBP1`, `HDAC8`, `HMGB2`, `HMGXB4`, `HSBP1`, `LTBR`, `MED13`, `MED15`, `MED25`, `MLLT1`, `MTA1`, `NFE2L1`, `NIT1`, `PHF14`, `PLAGL2`, `PRCP`, `RNF2`, `SHPRH`, `SIN3B`, `SLIRP`, `SMARCA5`, `STAT6`, `STUB1`, `TARBP2`, `TRAM2`, `TRAPPC6A`, `TWF2`, `ZFP62`, `ZNF32`, `ZNF714`

**`d1_perts2025_x_perts2026_final`** -- perts 2025 ∩ perts 2026 `final`, **28** genes:

`ACLY`, `ADNP`, `AKT2`, `ANKZF1`, `ASF1A`, `BRPF1`, `DCK`, `FUBP1`, `HDAC8`, `HMGB2`, `HMGXB4`, `LTBR`, `MED13`, `MED25`, `MLLT1`, `MTA1`, `NIT1`, `PHF14`, `PLAGL2`, `PRCP`, `SHPRH`, `SIN3B`, `SLIRP`, `STAT6`, `STUB1`, `TARBP2`, `TRAM2`, `ZFP62`

**`d2_perts2025_x_perts2026_val`** -- perts 2025 ∩ perts 2026 `val`, **25** genes:

`ACLY`, `ADNP`, `AKT2`, `ANKZF1`, `BRPF1`, `HDAC8`, `HSBP1`, `MED13`, `MED15`, `MED25`, `MTA1`, `NFE2L1`, `PLAGL2`, `RNF2`, `SHPRH`, `SIN3B`, `SLIRP`, `SMARCA5`, `STAT6`, `TARBP2`, `TRAPPC6A`, `TWF2`, `ZFP62`, `ZNF32`, `ZNF714`

**`d4_perts2025_x_interseccion_perts2026`** -- perts 2025 ∩ (interseccion de perts 2026), **16** genes:

`ACLY`, `ADNP`, `AKT2`, `ANKZF1`, `BRPF1`, `HDAC8`, `MED13`, `MED25`, `MTA1`, `PLAGL2`, `SHPRH`, `SIN3B`, `SLIRP`, `STAT6`, `TARBP2`, `ZFP62`

### 2.2 En que split de 2025 esta cada pareado

| #  | gen | split(s) de 2025 | panel(es) de 2026 |
|---|---|---|---|
| 1 | `ACLY` | validation | final, val |
| 2 | `ADNP` | test | final, val |
| 3 | `AKT2` | train | final, val |
| 4 | `ANKZF1` | test | final, val |
| 5 | `ASF1A` | test | final |
| 6 | `BRPF1` | test | final, val |
| 7 | `DCK` | test | final |
| 8 | `FUBP1` | validation | final |
| 9 | `HDAC8` | validation | final, val |
| 10 | `HMGB2` | train | final |
| 11 | `HMGXB4` | validation | final |
| 12 | `HSBP1` | train | val |
| 13 | `LTBR` | test | final |
| 14 | `MED13` | train | final, val |
| 15 | `MED15` | test | val |
| 16 | `MED25` | test | final, val |
| 17 | `MLLT1` | test | final |
| 18 | `MTA1` | train | final, val |
| 19 | `NFE2L1` | validation | val |
| 20 | `NIT1` | train | final |
| 21 | `PHF14` | train | final |
| 22 | `PLAGL2` | test | final, val |
| 23 | `PRCP` | train | final |
| 24 | `RNF2` | train | val |
| 25 | `SHPRH` | train | final, val |
| 26 | `SIN3B` | train | final, val |
| 27 | `SLIRP` | validation | final, val |
| 28 | `SMARCA5` | train | val |
| 29 | `STAT6` | train | final, val |
| 30 | `STUB1` | test | final |
| 31 | `TARBP2` | train | final, val |
| 32 | `TRAM2` | train | final |
| 33 | `TRAPPC6A` | train | val |
| 34 | `TWF2` | train | val |
| 35 | `ZFP62` | test | final, val |
| 36 | `ZNF32` | test | val |
| 37 | `ZNF714` | train | val |

## 3. La superficie que usan 4b y 4c, declarada antes de medir

| #  | cantidad | valor |
|---|---|---|
| 1 | split | `test` (CON verdad) |
| 2 | perturbaciones | 100 |
| 3 | celulas por perturbacion | `min(400, disponibles)` |
| 4 | 🔴 perturbaciones truncadas (< 400 celulas) | **15** de 100 |
| 5 | celulas perturbadas | 37,157 |
| 6 | pool de control `non-targeting` | 38,176 |
| 7 | control usado | 18,400 (el mismo numero que el reto 2026 lleva por contexto) |
| 8 | forma del h5ad real | 55,557 x 18,080 |
| 9 | `nnz` | 0 (0.0 M) |
| 10 | semillas | `20280201` (celulas/pert), `20280202` (control) -- familia NUEVA 20280000+ |
| 11 | huella sha256 del real | (se calcula en 4b/4c)... |
| 12 | dtype | el X de 2025 es `float32` con valores enteros; se convierte a `int64`, como el de Jiang |

## 4. Defectos de la construccion, como datos sobre la medicion

| #  | defecto | que implica |
|---|---|---|
| 1 | El panel de genes no es el mismo: 18,080 contra 18,533. | 18,077 genes en los dos. Cualquier cifra por gen cruzada entre las dos superficies vive en esa interseccion, no en el panel entero. |
| 2 | Las perturbaciones casi no se solapan. | 37 de 300 perturbaciones de 2025 aparecen en algun panel de 2026. El pareado por blanco es chico. |
| 3 | 🔴 15 de las 100 perturbaciones del split `test` tienen < 400 celulas. | La construccion de 4b/4c toma `min(400, disponibles)` y lo declara por perturbacion. No es el regimen del reto, que lleva 400 fijas. |
| 4 | El control se submuestrea de 38,176 a 18,400. | Es para igualar el `n` de control del reto 2026. El submuestreo tiene semilla declarada y descarta informacion a proposito. |
| 5 | La profundidad es muy distinta a la de Jiang. | ~54,000 UMI/celula en 2025 contra 3,319-6,203 de mediana en Jiang (`ntc_nulo.md` §3). Las dos superficies no son intercambiables. |
| 6 | El `X` de 2025 es `float32`. | Los valores son enteros; se convierte a `int64` para que el camino sea el mismo que en Jiang. Si hubiera no enteros, la conversion los truncaria -- no se detecto ninguno al construir. |

## 5. 🔴 Lo que esta superficie NO puede dar

| #  | lo que no da | por que |
|---|---|---|
| 1 | **Dominancia entre contextos** | Es UN solo contexto. No hay dos lineas ni dos experimentos que comparar dentro de 2025. |
| 2 | **Rango entre contextos** | Sin segundo contexto no hay rango. Toda cifra de aqui es un punto. |
| 3 | Heterogeneidad de respuesta entre contextos | La cantidad que `magnitud_centro.md` reparte entre CENTRO y contexto no esta definida con un solo contexto. |
| 4 | Transferencia a 2026 | El panel y la profundidad difieren (§4). Que algo ajuste aqui no dice que ajuste alla. |

## 6. Las entradas

| entrada | valor |
|---|---|
| datos 2025 | `/datos/vcc2025/{train,validation,test}` |
| genes 2025 | `/datos/vcc2025/gene_names.csv` |
| panel 2026 | `/datos/reto/final/` y `/datos/reto/val/` (`gene_names.csv`, `pert_counts.csv`) |
| salidas | `modelo/out/vcc2025_inventario.md`, `modelo/out/vcc2025_inventario.json` |

*Corrida: 2026-11-05 05:15:30 -> 2026-11-05 05:19:37 UTC, duracion 4.1 min (247 s).*
