# La entrega: la funcion, la emision y el contrato

## 🔴 Que es esto y que no

Un solo `.vcc` empaquetado de un solo `.h5ad`: 300 perturbaciones x 400 celulas x 3 contextos (D, E, F), **sin ninguna fila de control** — las celulas de control son la ENTRADA del modelo, no su salida.

🔴 **La verdad del gen no viaja.** De VCC 2025 sale **solo la funcion** (los coeficientes de la seccion 1): ningun delta de 2025 se copia a la emision, tampoco el de los 28 blancos que 2025 y el panel 2026 comparten. Son lineas distintas; esos 28 pasan por la misma funcion que los otros 272.

🔴 **No se ha enviado nada.** No se ha usado `vcc sample`, ni `vcc submit`, ni `~/vcc/benchmark/`, ni `/datos/benchmark_VEDADO_*`. No hay commits.

---

## 1. PASO 1 — la funcion aprendida en VCC 2025

### 1.1 El objetivo y el modelo, declarados antes de ajustar

| que | valor |
|---|---|
| objetivo | `y(g,h) = log2((cpm_pert + 1)/(cpm_ctrl + 1))`, cpm = media sobre celulas de `x/L*1e6` |
| modelo | ridge lineal, 10 rasgos estandarizados, intercepto, sin interacciones; alpha por CV partida por blanco |
| rasgos | `gwps_perfil`, `gwps_rango`, `gwps_signo`, `tiene_gwps_par`, `log10_cpm_h`, `det_h`, `log_fano_h`, `pasa_filtro_h`, `g_es_TF`, `log10_cpm_g` |
| rejilla de `alpha` | 0.1, 1.0, 10.0, 100.0, 1000.0 |
| blancos de 2025 | 300 |
| genes de 2025 | 18080 |
| genes que pasan el filtro CPM>5 del control de 2025 | 10780 |
| forma del diseno | [3234000, 10] |
| blancos de 2025 con perfil gwps | 280 |
| genes de 2025 con perfil gwps | 7581 |
| celulas de 2025 leidas | 491,046 |
| celulas de control de 2025 | 114,528 |

### 1.2 🔴 La validacion cruzada PARTIDA POR BLANCO

5 pliegues, permutacion de los blancos con la semilla 20284950. Ningun blanco esta a la vez en ajuste y en prueba. Se elige `alpha` por **la mediana de la correlacion de Pearson POR BLANCO retenido**.

| alpha | Pearson (todas las filas) | R² | Pearson por blanco (mediana) |  |
|---|---|---|---|---|
| `0.1` | +0.0735 | +0.0052 | **+0.0518** |  |
| `1.0` | +0.0735 | +0.0052 | **+0.0518** |  |
| `10.0` | +0.0735 | +0.0052 | **+0.0518** |  |
| `100.0` | +0.0735 | +0.0052 | **+0.0518** |  |
| `1000.0` | +0.0735 | +0.0052 | **+0.0518** | **elegida** |

**Los 5 pliegues de la `alpha` elegida (1000.0):**

| pliegue | n blancos | Pearson | Spearman | R² | RMSE | Pearson por blanco: mediana | [min, max] |
|---|---|---|---|---|---|---|---|
| 0 | 60 | +0.0799 | +0.0627 | +0.0064 | 0.1697 | +0.0579 | [-0.0548, +0.2485] |
| 1 | 60 | +0.0817 | +0.0605 | +0.0066 | 0.1748 | +0.0488 | [-0.0609, +0.2582] |
| 2 | 60 | +0.0608 | +0.0673 | +0.0030 | 0.1957 | +0.0503 | [-0.0368, +0.3103] |
| 3 | 60 | +0.0712 | +0.0518 | +0.0048 | 0.1583 | +0.0424 | [-0.0972, +0.2517] |
| 4 | 60 | +0.0741 | +0.0774 | +0.0054 | 0.1830 | +0.0595 | [-0.0032, +0.2000] |

🔴 **La distribucion de la correlacion por blanco retenido** (los 300 blancos, cada uno predicho por un modelo que no lo vio):

| cantidad | valor |
|---|---|
| mediana | +0.0507 |
| media | +0.0628 |
| [p10, p90] | [+0.0013, +0.1456] |
| [min, max] | [-0.0972, +0.3103] |
| fraccion > 0 | 0.9000 |
| fraccion > 0.1 | 0.2200 |

**20 blancos retenidos, equiespaciados por su correlacion** (el intermedio crudo de la cifra de arriba):

| # | blanco | Pearson del blanco retenido |
|---|---|---|
| 1 | `NDUFB4` | -0.0972 |
| 2 | `CHMP3` | -0.0140 |
| 3 | `HMGN1` | +0.0016 |
| 4 | `PHF14` | +0.0086 |
| 5 | `MED13L` | +0.0180 |
| 6 | `NEK2` | +0.0242 |
| 7 | `ASF1A` | +0.0284 |
| 8 | `UQCRQ` | +0.0343 |
| 9 | `ZBTB39` | +0.0406 |
| 10 | `NCAPH2` | +0.0483 |
| 11 | `IRF3` | +0.0534 |
| 12 | `LBR` | +0.0604 |
| 13 | `TMBIM6` | +0.0676 |
| 14 | `ARID3A` | +0.0789 |
| 15 | `JAZF1` | +0.0896 |
| 16 | `HMGB2` | +0.1022 |
| 17 | `MAU2` | +0.1165 |
| 18 | `EIF3H` | +0.1416 |
| 19 | `KLHDC2` | +0.1853 |
| 20 | `SOX2` | +0.3103 |

### 1.3 Los dos controles declarados

| control | Pearson | R² |
|---|---|---|
| (a) predecir 0 en todo gen | n/d (constante) | -0.0000 |
| (b) solo el perfil de gwps, una pendiente | +0.0430 | +0.0000 |
| el modelo de 10 rasgos | +0.0735 | +0.0052 |

### 1.4 Los coeficientes (rasgos estandarizados)

Intercepto `b0` = +0.000365, `alpha` = 1000.0.

| rasgo | peso (sobre rasgo estandarizado) | media | sd |
|---|---|---|---|
| `gwps_perfil` | +0.0024 | -3.844e-02 | 5.485e-01 |
| `gwps_rango` | +0.0130 | -5.334e-04 | 2.251e-01 |
| `gwps_signo` | -0.0043 | -5.021e-02 | 7.825e-01 |
| `tiene_gwps_par` | +0.0023 | 6.148e-01 | 4.866e-01 |
| `log10_cpm_h` | -0.0238 | 1.670e+00 | 5.015e-01 |
| `det_h` | +0.0085 | 7.848e-01 | 2.342e-01 |
| `log_fano_h` | +0.0097 | 2.539e-01 | 1.868e-01 |
| `pasa_filtro_h` | +0.0000 | 1.000e+00 | 1.000e+00 |
| `g_es_TF` | +0.0017 | 3.433e-01 | 4.748e-01 |
| `log10_cpm_g` | -0.0008 | 2.006e+00 | 6.523e-01 |

**La escala de lo predicho contra lo observado** (en muestra, sobre las filas del ajuste):

| cantidad | y observada | y predicha |
|---|---|---|
| sd | +0.1772 | +0.0138 |
| p1 | -0.5259 | -0.0302 |
| p50 | +0.0000 | +0.0010 |
| p99 | +0.5305 | +0.0328 |

🔴 La sd de lo predicho es mucho menor que la de lo observado: una ridge encoge, y el encogimiento es parte del modelo declarado, no un ajuste posterior.

---

## 2. PASO 2 — la emision sobre D/E/F

| que | valor |
|---|---|
| filas escritas | 360,000 de 360,000 |
| nnz | 1,996,620,401 |
| recorte de `y` en log2 | ±2.0 |
| blancos sin perfil gwps (repliegue declarado) | 23: `ANG`, `ARSA`, `ATP10D`, `CAPRIN2`, `HCFC2`, `HEATR6`, `ITPRIPL1`, `LYSMD4`, `MAPK7`, `MLKL`, `NEK3`, `PARP3`, `PDK3`, `PSMB9`, `RAB26`, `SEMA4F`, `SIK2`, `STK10`, `TBC1D19`, `TESK2`, `TMEM79`, `UBTD2`, `ZDHHC13` |
| tamano del h5ad | 14.93 GiB |

**Por contexto** — los rasgos del gen se recalculan con el control de ESE contexto:

| ctx | n control | profundidad mediana | genes que pasan CPM>5 | sd de y | y p1 | y p50 | y p99 | razon [min, max] | recortados | semilla control | semilla ruido |
|---|---|---|---|---|---|---|---|---|---|---|---|
| D | 18,400 | 20,236 | 9,633 | +0.0520 | -0.0372 | +0.0115 | +0.1579 | [0.894, 1.145] | 0 (0.00e+00) | 20284100 | 20284101 |
| E | 18,400 | 19,681 | 9,743 | +0.0519 | -0.0371 | +0.0105 | +0.1579 | [0.896, 1.148] | 0 (0.00e+00) | 20284200 | 20284201 |
| F | 18,400 | 19,845 | 9,537 | +0.0533 | -0.0367 | +0.0116 | +0.1565 | [0.896, 1.146] | 0 (0.00e+00) | 20284300 | 20284301 |

---

## 3. PASO 3 — el contrato, comprobado por mi antes de empaquetar

**17/17 comprobaciones OK.** Archivo 14.934 GiB, nnz 1,996,620,401 (margen x2.379 contra el tope de 4,750,000,000), densidad 0.29926, 5546.2 genes no nulos por celula, suma media por celula 20872.2.

| # | comprobacion | medido | esperado | ok |
|---|---|---|---|---|
| 1 | encoding de X | csr_matrix | csr_matrix | **si** |
| 2 | forma | [360000, 18533] | [360000, 18533] | **si** |
| 3 | var = gene_names.csv, mismo orden elemento a elemento | 18533 nombres, primer desajuste: None | 18533, None | **si** |
| 4 | filas de obs | [360000, 360000] | [360000, 360000] | **si** |
| 5 | CERO filas non-targeting | 0 | 0 | **si** |
| 6 | blancos = las 300 de pert_counts.csv | 300 distintos | 300 | **si** |
| 7 | 1,200 filas por blanco (400 x 3) | min 1200 max 1200 | 1200 | **si** |
| 8 | simbolos, no identificadores de constructo | ninguno con sufijo -N | ninguno | **si** |
| 9 | contextos | {'D': 120000, 'E': 120000, 'F': 120000} | {'D': 120000, 'E': 120000, 'F': 120000} | **si** |
| 10 | 120,000 filas por contexto | [120000, 120000, 120000] | [120000, 120000, 120000] | **si** |
| 11 | nnz <= 4,750,000,000 | 1,996,620,401 | <= 4,750,000,000 | **si** |
| 12 | indptr consistente | [0, 1996620401, 360001] | [0, 1996620401, 360001] | **si** |
| 13 | X entera | 0 valores no enteros | 0 | **si** |
| 14 | X no negativa | min 1.0 | >= 0 | **si** |
| 15 | X finita | 0 no finitos | 0 | **si** |
| 16 | indices dentro de rango | [0, 18532] | [0, 18532] | **si** |
| 17 | ninguna celula suma > 1,000,000 | maxima 53,536 | <= 1,000,000 | **si** |

### `vcc prep --dry-run`

```
✓ dry run — no .vcc written
  input:         /datos/modelo/out/entrega/prediction.h5ad
  output:        (dry run)
  cells:         360000
  genes:         18533  (encoding: float32)
  perturbations: column 'target_gene' → 'target_gene'
  contexts:      column 'context' → D=120000, E=120000, F=120000
  targets:       verified against the official list
  normalization: counts-preserved
```


### `vcc prep`

```
✓ prep complete
  input:         /datos/modelo/out/entrega/prediction.h5ad
  output:        /datos/modelo/out/entrega/prediction.vcc
  cells:         360000
  genes:         18533  (encoding: float32)
  perturbations: column 'target_gene' → 'target_gene'
  contexts:      column 'context' → D=120000, E=120000, F=120000
  targets:       verified against the official list
  normalization: counts-preserved
```


---

## 4. Defectos de la construccion

1. **Una sola linea de entrenamiento.** 2025 es un contexto; los rasgos de CONTEXTO (profundidad, n de control) no son identificables ahi y por eso no entran en el ajuste. Lo que se transfiere a D/E/F es una funcion ajustada en un solo contexto.
2. **El perfil de gwps es de K562 y cubre 7,681 de los 18,533 genes** y 277 de las 300 perturbaciones. Para el resto, el rasgo entra en su valor de no-informacion y la prediccion se queda con la parte comun.
3. **La ridge encoge.** La sd de `y` predicha es una fraccion de la observada: la emision es conservadora por construccion.
4. **El emisor desplaza medias, no dispersiones.** `_reescala_entero` aplica Binomial/Poisson sobre una celula de control; `sobredispersion.md` ya midio que ese dispositivo no reproduce la sobredispersion del delta real.
5. **El objetivo es una razon de CPM medio, no un estadistico de DE.** Es lo que el emisor puede aplicar; las metricas de DE del marcador miran otra cosa (Wilcoxon por gen), y la relacion entre ambas no se mide aqui.
6. **Las 400 celulas salen del propio pool de control con reemplazo**, asi que la variabilidad celula a celula de la emision es la del control mas el ruido de la inyeccion.

_Generado 2026-11-05 19:12:07 UTC._
