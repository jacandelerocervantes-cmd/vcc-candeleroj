# La convencion de `cell_eval2/ceiling.py` -- lectura de codigo, sin ejecutar

Fuente: `~/cell-eval2/src/cell_eval2/`, version 0.20.0 (la misma que puntuo
`bucle_suelo.md` y `bucle_techo.md`, confirmado por
`cell_eval2.__version__ == "0.20.0"`). Nada se ejecuto, nada se modifico, nada
se reimplemento. Todas las citas llevan archivo y linea.

---

## Respuesta corta al punto 7, porque cambia todo lo demas

**`ceiling.py` NO produce el `r` de `s = (u - b)/(r - b)`. Es una herramienta
aparte, terminal, que nadie consume.** El `r` lo produce **`anchor.py`**
(`compute_replicate_anchor`), y `anchor.py` **rechaza explicitamente** la
correccion de Spearman-Brown que aplica `ceiling.py`.

Asi que si: **la convencion de `r` esta fijada** -- pero por `anchor.py`, no por
`ceiling.py`. Y la convencion oficial resulta ser **5 particiones, crudas, sin
correccion de profundidad, con el control PARTIDO y una sola direccion**. El
detalle y las citas estan en el punto 7.

---

## 1. `SB_METRICS` completa y literal

`ceiling.py:114-135`. Dieciocho metricas, `frozenset`:

```
SB_METRICS: frozenset[str] = frozenset(          # ceiling.py:114
    {
        "delta_pearson",                          # 116
        "pds_l1",                                 # 117
        "pds_l2",                                 # 118
        "pds_cosine",                             # 119
        "de_wilcoxon_overlap",                    # 120
        "de_wilcoxon_overlap_top50",              # 121
        "de_wilcoxon_overlap_top100",             # 122
        "de_wilcoxon_overlap_top200",             # 123
        "de_wilcoxon_overlap_top500",             # 124
        "de_wilcoxon_precision",                  # 125
        "de_wilcoxon_precision_top50",            # 126
        "de_wilcoxon_precision_top100",           # 127
        "de_wilcoxon_precision_top200",           # 128
        "de_wilcoxon_precision_top500",           # 129
        "de_wilcoxon_nsig_spearman",              # 130
        "de_wilcoxon_lfc_spearman",               # 131
        "de_wilcoxon_direction_match",            # 132
        "de_wilcoxon_sig_recall",                 # 133
    }
)                                                 # 135
```

La regla general de pertenencia, `ceiling.py:91-93`:

> `# Reliability metrics that receive the Spearman-Brown ceiling correction. Names are`
> `# the wilcoxon-family canonical v2 names. Verified 1:1 against the validated`
> `# cell-eval implementations (same computation, not just the same name).`

Y el alcance, `ceiling.py:33-40`:

> "Spearman-Brown is applied ONLY to the metrics in `SB_METRICS` - a
> hand-maintained list (deliberately NOT derived from `MetricSpec.scoring`) ...
> Every other metric the run EMITS is reported as `NaN` (no defensible ceiling)"

### Cruce contra las SEIS metricas vcc2026

| metrica vcc2026 | en `SB_METRICS` | linea | razon de exclusion que da el codigo |
|---|---|---|---|
| `pds_cosine` | **SI** | `ceiling.py:119` | - |
| `expr_mse_unbiased_capped_norm` | no | - | `ceiling.py:99-101`: "error metrics (expr_mae/mse, expr_mse_unbiased_capped_norm, delta_mae/mse) and counts (de_*_nsig_counts_*). expr_mse_unbiased_capped_norm is additionally signed, so it is not a bounded reliability in the first place." |
| `de_wilcoxon_lfc_nmae` | no | - | `ceiling.py:99`: cae bajo "error metrics". No se la nombra individualmente. |
| `de_wilcoxon_direction_fidelity_yield_raw` | no | - | `ceiling.py:107-111`: "the #187 direction metrics (direction_precision, direction_sensitivity, direction_sensitivity_universe) - v2-native, no cell-eval equivalent, never validated under doubling". **El nombre exacto `direction_fidelity_yield_raw` NO aparece en el comentario**; queda cubierto por la regla general de la linea 92-93, no por una mencion propia. |
| `de_wilcoxon_direction_reach_raw` | no | - | Igual que la anterior: **no esta nombrada en el codigo**. `metrics/direction.py:28` la menciona en otro contexto ("replaced produced the metric's own ceiling"), no como exclusion. |
| `de_wilcoxon_sig_jaccard` | no | - | **No esta nombrada en el comentario de exclusiones.** La unica base en el codigo es la regla general (`ceiling.py:92-93`, "verified 1:1 against the validated cell-eval implementations") mas `metrics/de.py:195-196`, donde el propio jaccard se declara "This v2-native metric has no upstream cell-eval counterpart". La inferencia es directa pero **el codigo no la escribe**. |

**Lectura para este proyecto:** de las seis metricas que decide `vcc2026`, el
modulo de techo del marcador solo cubre UNA. Las otras cinco saldrian `NaN`
(`ceiling.py:216`). Eso ya se habia anotado en `bucle_techo.md`; aqui queda con
la linea.

Las demas exclusiones que el codigo si nombra, para completar el punto:
`de_wilcoxon_pr_auc` / `de_wilcoxon_roc_auc` (`ceiling.py:96-98`, el piso de
p-valor difiere del de cell-eval); las chance-corrected v2-nativas
`overlap_adjusted`, `precision_adjusted`, `sig_recall_adjusted`, `sig_mcc`
(`102-103`, "no cell-eval equivalent, never validated");
`model_direction_match` y el par `lfc_spearman_pos/neg` (`104-106`, tienen alias
v1 pero no estaban en el conjunto validado); y los nombres del backend deseq2
`de_deseq2_*` (`112-113`).

---

## 2. Como construye la particion

Funcion `_disjoint_halves`, `ceiling.py:138-179`.

**Mitades, no otra proporcion.** `ceiling.py:157`:

```
h = perm.size // 2
```

y `ceiling.py:161-162`:

```
a_idx.append(perm[:h])
b_idx.append(perm[h : 2 * h])
```

Las dos mitades tienen exactamente `floor(n/2)` celulas cada una. El
`perm[h : 2*h]` es lo que hace que la celula sobrante de un `n` impar **se
descarte**, no se asigne. Declarado en `ceiling.py:9-10` y `146-147`:

> "the leftover cell of an odd `n` is dropped" ... "When `n` is odd the one
> leftover cell is discarded (both halves must be the same depth for the doubling
> to hold)."

**Estratifica por perturbacion**, incluido el control. `ceiling.py:155`:

```
for _pert, idx in real.obs.groupby(pert_col, observed=True).indices.items():
```

`pert_col` es el del config del llamador (`ceiling.py:261`). No estratifica por
lote ni por nada mas: **no hay ninguna otra clave de estratificacion en el
codigo**. `ceiling.py:143`: "Each perturbation's cells (including the control's)".

**Reglas de descarte.** Un grupo con `n < 2` se cae de las dos mitades
(`ceiling.py:158-160`), y se avisa por `logger.warning` que la media del techo
queda sobre un conjunto distinto al del run principal (`165-172`). Si NINGUN
grupo se puede partir, `ValueError` (`163-164`). Si el **control** no se puede
partir, es error fatal (`175-178`):

```
if control not in set(half_a.obs[pert_col].astype(str)):
    raise ValueError(
        f"control {control!r} has < 2 cells; cannot compute a disjoint-split ceiling"
    )
```

**La semilla.** Un unico generador para todos los grupos, `ceiling.py:151`:

```
rng = np.random.default_rng(seed)
```

`seed` es parametro de `compute_ceiling` con **default 0** (`ceiling.py:227`:
`seed: int = 0`), y por CLI viene de `--ceiling-seed`, tambien **default 0**
(`cli.py:342-343`). Como es un solo `rng` consumido dentro del bucle del
`groupby`, el reparto de cada grupo depende del ORDEN del `groupby` -- eso el
codigo no lo comenta.

---

## 3. ¿Una particion o varias replicas?

**UNA sola.** `ceiling.py:261`, una unica llamada:

```
half_a, half_b = _disjoint_halves(real_ad, cfg.pert_col, cfg.control, seed)
```

No hay bucle de semillas en `ceiling.py`. Lo dice el docstring del modulo,
`ceiling.py:4-5`:

> "It is an empirical reproducibility estimate rather than a proved upper bound -
> it is measured on **ONE random split**, and a model that denoises the reference
> better than a half of it does can in principle score above it."

Por lo tanto **no hay agregacion entre replicas**: la pregunta de media /
mediana / percentil no aplica a este modulo. (El que si hace replicas es
`anchor.py`, con k=5 y media -- punto 7.)

---

## 4. Spearman-Brown: formula exacta, a que se aplica, con que factor

**Si la aplica.** Funcion `_spearman_brown`, `ceiling.py:182-220`. La formula,
literal, `ceiling.py:214`:

```
ceilings.append(2.0 * v / (1.0 + v))
```

es decir `r' = 2r/(1+r)`, tambien escrita en prosa en `ceiling.py:18`.

**A que cantidad se le aplica.** A `v`, que es la **media por metrica** que
devuelve `aggregate_metrics` sobre el self-split. `ceiling.py:208`:

```
means = dict(zip(agg_mean["metric"].to_list(), agg_mean["mean"].to_list()))
```

alimentada en `ceiling.py:308-309`:

```
agg_ceiling = _spearman_brown(
    aggregate_metrics(results, metrics=sb_run), metric_output_names(cfg)
)
```

O sea: no se aplica por perturbacion, se aplica **una vez a la media agregada de
cada metrica**. El frame por perturbacion (`results`) sale SIN corregir
(`ceiling.py:236-237`).

**Con que factor.** Factor 2, el de duplicar la profundidad: de la profundidad de
media poblacion a la profundidad combinada `2 * floor(n/2)`. `ceiling.py:15-18`:

> "each metric's per-context mean is mapped from half depth to the COMBINED split
> depth `2 * floor(n/2)` - the run's full depth for even `n`, one cell short for
> odd `n` - with the analytical Spearman-Brown correction `r' = 2r/(1+r)`"

**La guarda.** Solo se corrige si la metrica esta en `SB_METRICS` por su nombre
CANONICO y el valor es finito y **estrictamente positivo** (`ceiling.py:213`):

```
if canon in SB_METRICS and v is not None and math.isfinite(v) and v > 0.0:
```

Si no, `NaN` (`ceiling.py:216`). El motivo, `ceiling.py:198-206`: en `r <= 0` la
formula deja de ser una correccion de fiabilidad (su polo esta en `r == -1`, y en
`r = -0.9` devuelve `-18`), asi que una fiabilidad no positiva se reporta como
"no hay techo defendible" en vez de un techo negativo.

**Advertencia de version que el propio modulo trae**, `ceiling.py:52-64`: la
correccion fue validada contra numeros de cell-eval con convenciones
*v1-equivalentes* (`control_source="pred"`, PDS `rank_denominator="n"`,
`nan_lfc_policy="keep"`, `min_abs_log2fc=0`). Bajo los defaults v2 el algoritmo
es el mismo pero las convenciones difieren, y la extrapolacion se aplica "on the
same empirical footing rather than a separately validated one".

---

## 5. La direccion

**Corre UNA sola direccion: `half_b` como prediccion, `half_a` como verdad.**
`ceiling.py:287-289`:

```
results = compute_metrics(
    half_b,
    half_a,
```

El orden de `compute_metrics` es `(pred, real)` -- se confirma en
`ceiling.py:45-46` ("scoring `pred=half_b` against `real=half_a`") y en
`anchor.py:109` ("Score half_b against half_a for ONE seed"), cuyo nucleo repite
la misma llamada en `anchor.py:120`: `tidy = compute_metrics(half_b, half_a, config=inner)`.

**No corre la otra direccion, y el codigo no justifica por que esa y no la
otra.** El docstring solo dice, `ceiling.py:11`: "one half is treated as 'real'
and the other as 'prediction'", sin decir cual ni por que. **No esta en el
codigo**: no hay comentario sobre la asimetria de elegir una direccion, ni sobre
promediarlas, ni sobre si alguna metrica es asimetrica bajo el intercambio.

---

## 6. ¿Comparte el control o lo parte?

**Lo PARTE, y ademas fuerza que cada mitad use el suyo.** Son dos mecanismos
distintos y los dos estan presentes:

1. **El control se parte como cualquier otro grupo.** El `groupby` de
   `ceiling.py:155` no excluye la etiqueta de control; `ceiling.py:143` lo dice:
   "Each perturbation's cells (**including the control's**) are shuffled and
   split".

2. **`control_source` se fuerza a `"pred"`.** `ceiling.py:296`, dentro del
   `replace` del config:

   ```
   control_source="pred",
   ```

   Esto es lo que hace que el DE de cada lado tome su referencia de su PROPIO
   bloque de control. El motivo, `ceiling.py:276-286`:

   > "`control_source`: MUST be "pred" here, whatever the caller uses. Under
   > "real" the pred side's DE takes its reference from the real side, so with
   > pred=half_b and real=half_a BOTH halves' log2FCs are computed against
   > half_a's control. Their sampling noise is then shared rather than independent,
   > which correlates the two quantities whose agreement is being measured and
   > biases split-half reliability upward. **Measured on a small self-split:
   > lfc_spearman 0.54 -> 0.74, and nsig_spearman went from NaN (no defensible
   > ceiling) to a confident-looking 0.94**, i.e. the shared control can
   > manufacture reliability where there is none and defeat the `r > 0` guard."

   Y `ceiling.py:48-49`: "Independent controls per half are what make them
   replicates, so **the override is not optional**".

**Traduccion a nuestras variantes:** `ceiling.py` implementa **nuestra variante I**
(control partido, `control_source=pred`). **No** implementa la variante S, y la
documenta como el error que sesga el techo hacia arriba. `ceiling.py:60-64`
reconoce el costo de esa decision: contra un run principal que usa
`control_source="real"` -- que es exactamente el preset `vcc2026`
(`configs/vcc2026.yaml:15`) -- "the ceiling is not measured under identical
conventions", y lo acepta: "a biased estimator is worse than a convention
mismatch, and an upward-biased ceiling understates every score measured against
it."

---

## 7. 🔴 ¿Este modulo produce el `r` de `s = (u - b)/(r - b)`?

**No. Es una herramienta aparte cuya salida no la consume nadie.**

### Quien importa `ceiling.py`

Busqueda sobre todo el paquete (`grep -rn "ceiling" --include=*.py`). Hay
exactamente cuatro importaciones:

| importador | linea | que importa | para que |
|---|---|---|---|
| `cli.py` | `cli.py:14` | `from .ceiling import compute_ceiling` | la bandera `--ceiling` |
| `__init__.py` | `__init__.py:5` | `from .ceiling import compute_ceiling` | lo expone en la API publica (`__init__.py:37`) |
| `anchor.py` | `anchor.py:35` | `from .ceiling import _disjoint_halves` | **solo la particion**, no la correccion |
| `lfc_nmae_ref.py` | `lfc_nmae_ref.py:60` | `from .ceiling import _disjoint_halves` | **solo la particion** |

**`score.py` NO importa `ceiling.py`.** Sus imports estan en `score.py:28-43`; el
relevante es `score.py:40`:

```
from .anchor import AnchorExpect, anchor_digest, resolve_anchor
```

### Donde muere la salida de `compute_ceiling`

Un unico punto de consumo, `cli.py:747-753`:

```
if args.ceiling:                                                    # cli.py:747
    cres, cagg = compute_ceiling(args.adata_real, config=cfg,       # 748
                                 seed=args.ceiling_seed)           # 749
    ...
    cres.write_csv(os.path.join(cfg.outdir, "ceiling_results.csv")) # 751
    cagg.write_csv(os.path.join(cfg.outdir, "ceiling_agg.csv"))     # 752
    print(os.path.join(cfg.outdir, "ceiling_agg.csv"))              # 753
```

Se escriben dos CSV y se imprime la ruta. **Ningun modulo del paquete lee
`ceiling_agg.csv` ni `ceiling_results.csv`** -- las unicas apariciones de esos
nombres son esas tres lineas de escritura/impresion y dos comentarios en
`ceiling.py:303,307`. Es salida terminal. El propio `compute_ceiling` lo dice,
`ceiling.py:246`: "Neither writes anything: the two frames are returned, not
persisted."

### Quien SI produce el `r`

`anchor.py`, funcion `compute_replicate_anchor` (`anchor.py:178-305`). Se consume en
`score.py`, funcion `_replicate_entries` (`score.py:394`), cuyo docstring es la
cita que cierra la pregunta -- `score.py:395-399`:

> "The replicate scale as `(base, Scoring)` pairs: 0 = baseline, 1 = replicate.
>
> #276 comment 1 fixes the scale; part C makes it POLICY-APPLIED. `score_one`
> with `policy.anchor` set to the measured replicate has **`(u - b) / (r - b)`**
> as its unclamped linear core for both directions"

Y la asignacion concreta, `score.py:455`:

```
policy = replace(spec.scoring, anchor=float(rep), allow_negative_baseline=False)
```

`rep` es el `replicate` del artefacto de `anchor.py`. Es decir: el `r` de la
escala **es `Scoring.anchor`**, y se rellena con la replica medida por
`anchor.py`.

Dos matices de nomenclatura que conviene no confundir:

- `Scoring.anchor` (`scoring.py:45`) es "the value that SCORES 1.0"
  (`scoring.py:37`). En las politicas por defecto es una **constante**:
  `BOUNDED = Scoring(..., anchor=1.0)` (`scoring.py:189`),
  `ERROR = Scoring(..., anchor=0.0)` (`scoring.py:161`). Con `--anchor` esa
  constante se **sustituye** por la replica medida (`score.py:455`).
- El modulo `anchor.py` es el "replicate anchor" de #276, otra cosa que el campo
  `Scoring.anchor`, aunque uno acabe llenando al otro.

### Y lo decisivo: `anchor.py` RECHAZA la correccion de `ceiling.py`

`anchor.py:1-7`, docstring del modulo:

> "Replicate anchor: the 1.0 end of #276's 0=baseline/1=replicate scale.
>
> A RAW split-half replicate, measured on the real data alone. **Deliberately NOT
> Spearman-Brown corrected** and not depth-corrected in any other way -- the
> correction `ceiling.py` applies is a bounded-reliability extrapolation that
> **does not hold for five of the six scored `vcc2026` members**, and as of
> 2026-08-12 **no depth correction is applied anywhere in this scheme**."

Y reusa de `ceiling.py` solo la particion, con la misma razon sobre el control
(`anchor.py:8-18`):

> "ONE shared split-and-score core, so the split semantics never get a second
> implementation. `control_source` is forced to "pred" inside that core and is
> not optional ... **and under this scheme an upward-biased anchor biases every
> submission's score DOWNWARD, uniformly, with nothing in the output that looks
> wrong.** Measured for de_lfc_nmae on a real panel: shared controls read 0.5-2.3%
> optimistic, consistently, on all 12 comparisons. `ceiling.py` measured the same
> mistake turning a 0.54 split-half reliability into 0.74."

### La convencion oficial de `r`, entonces

| aspecto | `ceiling.py` (NO es el `r`) | `anchor.py` (**SI** es el `r`) | cita |
|---|---|---|---|
| particion | `_disjoint_halves`, `floor(n/2)`, impar descartado | **el mismo** `_disjoint_halves` | `anchor.py:35`, `117` |
| replicas | 1 | **5** (`n_splits: int = 5`) | `anchor.py:183`; fijado `N_SPLITS = 5` en `competition.py:41` |
| semillas | `seed=0` | `numpy.random.SeedSequence(base_seed).generate_state(n_splits)` desde `base_seed=0` | `anchor.py:74` (`SEED_DERIVATION`), `135-141`, `182`; `competition.py:38` |
| agregacion entre replicas | no aplica | **media de los AGREGADOS** (no pool de filas) | `anchor.py:285` (`replicate=pl.col("value").mean()`); razon en `189`, `191-196` |
| Spearman-Brown | si, `2r/(1+r)` | **NO**, cruda | `ceiling.py:214` vs `anchor.py:3-7` |
| correccion de profundidad | si (la SB) | **ninguna, en ningun lugar del esquema** | `anchor.py:6-7` |
| direccion | una, `half_b` vs `half_a` | **una, la misma** | `ceiling.py:288-289`; `anchor.py:109`, `117`, `120` |
| control | partido + `control_source="pred"` | **partido + `control_source="pred"`**, mas una asercion de disjuncion | `anchor.py:103` (el `replace`), `91-92` (el motivo), `118` (`_assert_disjoint_controls`) |
| cobertura de metricas | solo `SB_METRICS` (1 de las 6) | **todas las seleccionadas**; si una no sale, `ValueError` | `anchor.py:256-264` |
| `lfc_nmae` | dentro del nucleo uniforme | **estimador DISTINTO** (`FULL_GATE_RAW`), gate y denominador de la tabla real COMPLETA | `anchor.py:57-66`, `154-176` |

**Consecuencia para este proyecto.** Nuestro techo a mano se recalibra contra
`anchor.py`, no compite con `ceiling.py`. En la convencion oficial:

- nuestra **variante I es la correcta** (control partido, `control_source=pred`)
  y la **variante S es justamente el sesgo documentado**;
- nuestra **dirección B-como-prediccion es la oficial**; la A no tiene
  contraparte en el marcador;
- nuestras mitades a tamano natural coinciden con "sin correccion de
  profundidad" (`anchor.py:6-7`), asi que **no haber aplicado Spearman-Brown fue
  la convencion oficial, no una omision**;
- nos faltan **4 particiones** (k=5) y sobra **una direccion**;
- y `de_wilcoxon_lfc_nmae` necesita **otro estimador**, no el nucleo uniforme.

---

## 8. Constantes ancladas (listadas, sin proponer usarlas)

En `ceiling.py`:

| constante | linea | que es |
|---|---|---|
| `seed: int = 0` | `ceiling.py:227` | semilla por defecto de la particion |
| `2.0` y `1.0` en `2.0 * v / (1.0 + v)` | `ceiling.py:214` | el factor de duplicacion de Spearman-Brown |
| `v > 0.0` | `ceiling.py:213` | umbral de fiabilidad no positiva -> `NaN` |
| `h < 1` | `ceiling.py:158` | umbral de celulas para poder partir un grupo (`n >= 2`) |
| `SB_METRICS` | `ceiling.py:114-135` | la lista de 18 nombres, "hand-maintained" |
| `float("nan")` | `ceiling.py:216` | el valor que recibe toda metrica fuera de la lista |

**No hay ningun `r` de referencia constante en `ceiling.py`**: el techo siempre
se mide, nunca se lee de una constante.

Constantes del esquema que SI fija el `r` oficial (`anchor.py` / `competition.py`),
listadas porque son las que mandan:

| constante | linea | valor |
|---|---|---|
| `BASE_SEED` | `competition.py:38` | `0` |
| `N_SPLITS` | `competition.py:41` | `5` ("k = 5, ruled 2026-08-13") |
| `BULK_TARGET_SUM` | `competition.py:45` | `50_000.0` |
| `SEED_DERIVATION` | `anchor.py:74` | `"numpy.random.SeedSequence(base_seed).generate_state(n_splits)"` |
| `SPLIT_HALF_RAW` | `anchor.py:55` | `"split_half_raw"` (estimador estampado por metrica) |
| `--ceiling-seed` default | `cli.py:342` | `0` |
| `--anchor-base-seed` default | `cli.py` (`run`) | `0` |
| `--anchor-splits` default | `cli.py` (`run`) | `5` |

Nota lateral sobre `BULK_TARGET_SUM`, porque toca al techo: `competition.py:43-44`
dice "#268: at the retired 1e6 the split-half ceiling is NEGATIVE on 6/6 real
lines, so an anchor built at that value is not merely different -- it is
unusable". El mismo hecho aparece en `config.py:316-321` y `moments.py:230`. Es
decir, el valor `bulk_target_sum = 50000.0` del preset `vcc2026`
(`configs/vcc2026.yaml`) **fue elegido en parte para que el techo de self-split
no saliera negativo**.

---

## Lo que no esta en el codigo

Para no inferirlo:

- **Por que `half_b` y no `half_a`** como lado predicho, y por que no las dos
  direcciones: no hay comentario en `ceiling.py` ni en `anchor.py`.
- **Por que `de_wilcoxon_sig_jaccard` queda fuera de `SB_METRICS`**: no esta
  nombrada en el comentario de exclusiones (`ceiling.py:95-113`). Solo se deduce
  de la regla general mas `metrics/de.py:195-196`.
- **Por que `de_wilcoxon_direction_fidelity_yield_raw` y
  `de_wilcoxon_direction_reach_raw` quedan fuera por su nombre exacto**: el
  comentario nombra otras tres metricas de la familia #187
  (`direction_precision`, `direction_sensitivity`,
  `direction_sensitivity_universe`), no estas dos.
- **El efecto del orden del `groupby`** sobre el reparto, dado que hay un solo
  `rng` para todos los grupos (`ceiling.py:151` + `155`).
- **Si `ceiling.py` sigue siendo codigo vivo o es un residuo** de antes de #276:
  el codigo no lo declara. Lo unico comprobable es que su salida no la consume
  nadie y que `anchor.py` describe su correccion como no aplicable a cinco de las
  seis metricas puntuadas.

---

## Procedencia

- Nada se ejecuto: ni `compute_ceiling`, ni `cell-eval2 run --ceiling`, ni nada
  que los llame. Solo lectura de fuente y `grep`.
- No se toco `bucle_techo.md`, `techo_arnes.py`, `bucle_techo.py` ni
  `techo_reporte.py`.
- No se leyo nada bajo `~/vcc/benchmark/` ni
  `/datos/benchmark_VEDADO_techo_metricas/`.
- Sin envios, sin `vcc sample`, sin commits.
- Archivos leidos: `ceiling.py`, `anchor.py`, `score.py`, `scoring.py`,
  `cli.py`, `competition.py`, `configs/vcc2026.yaml`, `io.py`, `de.py`,
  `metrics/de.py`, todos bajo `~/cell-eval2/src/cell_eval2/`.
