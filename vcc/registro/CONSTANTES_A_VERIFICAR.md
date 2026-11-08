# Constantes a verificar contra la especificación

**Por qué existe esta lista.** El build anterior sacó **0.135** de Overall — contra un techo de ~1.51,
eso es ~9 % de lo alcanzable. Así que **ninguna constante calibrada dentro de él tiene credibilidad**:
su única evidencia de corrección es un pipeline que casi no funcionó.

Y hay una trampa más fina: `spec_metrics.py` **afirma en su docstring** que sus constantes vienen de
la especificación (`docs/76 §1`). **Ese docstring lo escribió el mismo build de 0.135.** Sus
afirmaciones sobre su propia procedencia merecen el mismo escrutinio que sus números.

→ **Nada entra por lo que diga el archivo. Cada número se vuelve a derivar de la fuente primaria.**

**Fuente primaria:** el clon de `cell-eval2`, en `docs/vcc2026_metrics/`
(`vcc2026-metrics.md`, 46 KB · `vcc2026-metrics-brief.md`, 18 KB). Perfil `vcc2026`, `rule_version` 3.

---

## Lo que hay que confirmar, uno por uno

| # | constante | valor que afirmaba el build viejo | qué confirmar |
|---|---|---|---|
| 1 | piso de pureza de `reach` | `0.9` | el valor **y si la comparación es `>=` o `>`** |
| 2 | compuerta de `nmae` | `10` genes | el número **y si es `< min_gate` o `<= min_gate`** |
| 3 | alfa de BH | `0.05` | el valor, **que sea BH dentro de la perturbación**, y que sea `<` estricto |
| 4 | recorte de `mse` | `[0.0, 1.0]` | 🔴 **lo decisivo: que el recorte sea sobre el valor NORMALIZADO y no sobre el cociente crudo.** De eso depende el techo (1.513 con recorte vs ≈1.52 sin él) |
| 5 | piso de `nmae` | `−6.0` | **no estaba registrado en el documento.** Si se confirma, acota el daño máximo de la magnitud sin calibrar: −6 en la métrica = **−0.33 en el Overall** |
| 6 | `b` y `r` por métrica | puntos medios: pds .455 · mse .95 · fid .305 · reach .895 · jac .365 · nmae .60 | **tomar los RANGOS, no los puntos medios** (colapsarlos es el patrón de error de §18.8). Verificado: los rangos **no se traslapan**, así que el orden fid > jac > pds > nmae se sostiene.<br><br>✅ **Span de `pds_cosine` cerrado — 2026-10-27.** La *especificación* (`vcc2026-metrics.md`) no da el valor numérico; solo el companion `vcc2026-metrics-brief.md` lo mide. `b = 0.500` tiene dos fuentes independientes que coinciden: (a) constante estructural en `src/cell_eval2/scales.py:428` (`"pds_cosine": _higher(0.5, -1.0)`), derivada en el comentario de las líneas 425-427 («base 0.5: PDS is 1 - rank/D, so a uniform rank lands at 0.5» — no calibrada, es geometría del rango uniforme); (b) medición empírica en `docs/vcc2026_metrics/vcc2026-metrics-brief.md:346` (tabla «§8 Reference points»), confirmada como exactamente 0.500000 en los tres bundles oficiales (línea 353). `r = 0.927 – 0.984` **no existe como constante de código** — `r` se define (líneas 58-72 del brief) como la referencia partida a la mitad y comparada contra sí misma, una cantidad medida, no fijada — y su único valor registrado es esa misma fila de la tabla (`vcc2026-metrics-brief.md:346`), medida «on the three official reference bundles at `cell-eval2 0.15.0`, `rule_version` 3» (línea 341). `r − b = 0.427 – 0.484`, que redondea a **0.43 – 0.48** — coincide exactamente con la cuarta columna de esa misma fila. ⚠️ Esa medición es de `cell-eval2` 0.15.0; el clon que tenemos es 0.20.0 (ver `limitaciones.md`) — no afecta al número citado (es el publicado), sí a si *nosotros* podemos reproducirlo bit a bit con el clon actual. |
| 7 | universo de genes de `b` y `r` | 18,533 genes con > 5 CPM | 🔴 **b y r tienen que medirse en el MISMO universo que u.** Los oficiales no sirven para un u calculado sobre otro conjunto |
| 8 | «adjudicable» en `fid` | `ref_lfc` finito y ≠ 0 | ⚠️ **el propio archivo anota que el oficial acepta ±inf como dirección.** Es una discrepancia que él mismo marcó: hay que resolverla |
| 9 | casos borde de `fid` | n_real=0 y n_pred=0 → NaN · n_real=0 → k/n_pred · n_pred=0 → 0 | los tres contra el spec |
| 10 | `jac` con unión vacía | `1.0` | confirmar. Es la trampa de **D-023**: la implementación es fiel al scorer, **usarla no** |
| 11 | cadena de desempate de `reach` | llamadas primero → padj asc → p asc → \|lfc\| desc → **nombre del gen lexicográfico** («G10» < «G2») | ✅ **Confirmado contra el código — 2026-10-27.** `src/cell_eval2/metrics/direction.py:508-511`, dentro de `_purity_curve()` (la curva que `de_direction_reach`, línea 862, consume vía `_k_star`): `.sort(["target", "_sig_pred", "rank_p_adj", "rank_p_value", "abs_lfc_pred", "feature"], descending=[False, True, False, False, True, False])`. Los cinco criterios, en orden: `_sig_pred` desc (llamada primero, línea 506: `p_adj_pred < alpha`) → `rank_p_adj` asc → `rank_p_value` asc → `abs_lfc_pred` desc → `feature` asc. `feature` es el nombre/símbolo del gen (columna `Utf8`, confirmado en `de.py:293` y en el join de `direction.py:184`). `descending=False` sobre un `Utf8` en polars usa orden lexicográfico de bytes/código de punto — no orden numérico ni "natural sort" — así que comparando carácter a carácter «G10» vs «G2»: `G`=`G`, luego `1` < `2`, por lo tanto **«G10» < «G2»** queda confirmado tal cual. |
| 12 | predicción no finita en `nmae` | → `0` | confirmar |
| 13 | filtro de referencia | ≥ 5 CPM en el control | confirmar. **T6 ya retiró el 10 y el 20 CPM** como invenciones nuestras |
| 14 | normalización por célula | `target_sum = 1e6` | confirmar. **Es lo que valida `cpm_stats.py` (CPM sin pseudoconteo) y descarta `stats_pb.py` (log1p CP10k)** |

---

## Lo que NO entra, decidido

- **`leverage_weights()`** — pesos de pérdida ∝ 1/(r−b). Es una decisión de modelo, y **su propio
  docstring lista cuatro razones por las que está mal**: derivadas locales sobre la métrica cruda y no
  sobre la pérdida sustituta · el peso del mse recortado es un máximo, no un gradiente efectivo · el
  mse es razón de sumas, no media por perturbación · fid, jac y reach están acoplados por el conjunto
  de llamadas. El build anterior ya lo auditó y lo dejó marcado.
- **`ivw()`** del arnés — media ponderada por inversa de la varianza. Matemáticamente limpia, **sin
  parámetros libres**, y aun así **metodológicamente superada**: §21.8 estableció que bajo L1 el
  estimador correcto es la **mediana ponderada**, no la media.
- **`pipe_20_harness.py` completo** — no es el antepasado de X3. Es un **LOSO** que contesta «qué tan
  bien transfiere mi modelo entre mis fuentes», no «reproduce el resultado publicado uno a uno».
  Y está saturado: `K0` · `TEMP` · `N_PROGRAMS` · **`KD_TARGET` hardcodeado** · `clip(−3, 3)` ·
  `0.3 * n_genes` · `ns >= 2` · piso de `tau` en `1e-6` · `MIN_CPM` · `hill_dose()` con parámetros
  supuestos · `SEED` · `fit_distortion` · `strength_class` · `activity_gate` · `residual_prior`.

---

## Lo que sí se copió, y por qué

- **`pseudobulk.py`** — suma conteos crudos por perturbación + n células, por bloques. **Cero valores.**
- **`cpm_stats.py`** — Σ CPM y Σ CPM² por grupo, por bloques. **Cero valores calibrados.** El único
  float es `rs[~(rs > 0)] = 1.0`, una guarda de división por cero: si la biblioteca es 0, todos los
  conteos de esa célula son 0, así que el CPM sale 0 igual. Inocuo.
- El `1e6` de los dos es **la definición de CPM**, no una elección.

---

## La regla, para los 428 archivos restantes

> **Si puedes señalar la página de la especificación, la constante se queda.
> Si vino de un ajuste, un barrido o el criterio de alguien, se va.
> Y si solo el archivo lo afirma, eso no es señalar la página.**
