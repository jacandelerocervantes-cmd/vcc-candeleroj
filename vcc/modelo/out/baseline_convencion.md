# ¿Quien produce el `b`, y existe un camino nativo? -- lectura de codigo

Fuente: `~/cell-eval2/src/cell_eval2/`, version 0.20.0 (la misma que puntuo
`bucle_suelo.md` y `bucle_techo.md`). Nada se ejecuto salvo `--help` de los
subcomandos para transcribir sus banderas; nada se modifico; ningun arnes se
toco. Toda cita lleva archivo y linea.

Complementa `modelo/out/ceiling_convencion.md`, que resolvio el `r`.

---

## Respuesta corta

**A.** El `b` es una **emision medida que el marcador construye**: el
*generic-response baseline* de `baseline.py` -- una respuesta promedio de
perturbacion, emitida para TODAS las perturbaciones, y puntuada como una
submission cualquiera. No es constante y no la aporta el usuario como numero: el
usuario aporta el ARTEFACTO que el propio marcador genero con
`cell-eval2 baseline`.

**Y nuestro suelo NO es ese `b`.** Nuestro suelo emite el control **sin
escalar**; el `b` del marcador emite el control **escalado gen a gen** hasta la
respuesta promedio de perturbacion. Son dos comparadores distintos. Nuestro
suelo corresponde a otra cosa que el marcador tambien tiene nombrada: el extremo
0 de la escala constante `low-random_high-1_v10` (`scales.py:313`). Las citas
estan en A.4, y el cruce numerico en A.5 lo confirma.

**B.** Si existe el camino nativo, y son **cuatro invocaciones**. Deja sobrando
`suelo_arnes.py`, `suelo_floor.py`, `bucle_suelo.py` y `score_vs_self` por
completo. `techo_arnes.py` queda sobrando como productor de `r`, pero **no** como
lo que resulto ser: una medicion independiente de las seis metricas por
perturbacion, que el camino nativo no expone. Y si, el camino nativo resuelve
solo los dos huecos que habiamos identificado.

---

# A. ¿Quien produce el `b` de `s = (u - b)/(r - b)`?

## A.1 La cadena, de la formula al archivo

El `0` de "0 = baseline, 1 = replicate" (`score.py:395`) entra a la politica por
esta cadena, toda dentro de `score.py`:

| paso | linea | codigo |
|---|---|---|
| se lee la fila del estadistico comparado del frame de baseline | `score.py:731` | `b_row = results_base.filter(pl.col("statistic") == comparison_statistic).drop("statistic")` |
| se extraen los valores | `score.py:732` | `b_vals = b_row.row(0)` |
| se arma el mapa nombre -> `b` | `score.py:861` | `base_by_name = dict(zip(metric_names, b_vals))` |
| se pasa a la escala de replica | `score.py:969` | `entries = _replicate_entries(base_by_name, frame, exclude=ref_unscoreable)` |
| se canoniza | `score.py:414` | `base_by_canon = {_NAME_TO_CANONICAL.get(k, k): v for k, v in base_by_name.items()}` |
| se toma el `b` de la metrica | `score.py:443` | `b, problem, policy = base_by_canon.get(canon), None, None` |
| el `r` entra como `anchor` de la politica | `score.py:455` | `policy = replace(spec.scoring, anchor=float(rep), allow_negative_baseline=False)` |

Y el docstring que une las dos puntas, `score.py:395-399`:

> "The replicate scale as `(base, Scoring)` pairs: 0 = baseline, 1 = replicate.
>
> #276 comment 1 fixes the scale; part C makes it POLICY-APPLIED. `score_one`
> with `policy.anchor` set to the measured replicate has **`(u - b) / (r - b)`**
> as its unclamped linear core for both directions"

O sea: **`b` = `base`, que viene del frame `results_base`**; **`r` = `policy.anchor`,
que viene del artefacto de `anchor.py`**.

## A.2 ¿De donde sale `results_base`? Dos caminos, los dos con artefacto

**Camino suelto:** del CSV que el usuario pasa en `--baseline-agg`, descrito en
el `--help` de `score` como "wide agg CSV from `baseline`".

**Camino de competencia:** del bundle, `score.py:660`:

```
results_base = bundle.baseline_agg
```

y en el bundle ese archivo se llama `baseline_agg.csv`
(`real_bundle.py:41`: `BASELINE_AGG = "baseline_agg.csv"`), escrito en
`real_bundle.py:449`.

En los dos casos el contenido es el **agregado de una corrida real de metricas**
sobre una prediccion que el marcador construyo. **No es un numero declarado.**

## A.3 ¿Es constante? Si, pero en OTRA familia de escalas

Para la escala de replica (`(u-b)/(r-b)`) el `b` es medido, no constante. Pero el
paquete SI tiene `b` constantes, en `scales.py`, y conviene no confundirlos.
`scales.py:1-7`:

> "Named, frozen SCALES: fixed reference points that score against constant
> endpoints. ... It supplies, per metric, a constant zero point `base` plus a
> full `Scoring` policy. Nothing here is new arithmetic -- `scoring.score_one`
> already computes `(u - b) / (a - b)`; **a scale just hands it a constant `b`
> instead of a measured baseline**."

La escala que cubre las seis metricas de vcc2026 es `low-random_high-1_v10`
(`scales.py:313`), descrita en `scales.py:315-318` como "0 = the random minimum,
1 = real input (the real count matrix pasted as the prediction)". Sus seis bases
constantes:

| metrica vcc2026 | `base` constante | politica | linea |
|---|---|---|---|
| `expr_mse_unbiased_capped_norm` | **1.0** | `_lower(1.0, -6.0)` | `scales.py:369` |
| `de_wilcoxon_lfc_nmae` | **1.0** | `_lower(1.0, -1.0)` | `scales.py:418` |
| `pds_cosine` | **0.5** | `_higher(0.5, -1.0)` | `scales.py:428` |
| `de_wilcoxon_direction_fidelity_yield_raw` | **0.5** | `_higher(0.5, -1.0)` | `scales.py:430` |
| `de_wilcoxon_direction_reach_raw` | **0.0** | `_higher(0.0, 0.0)` | `scales.py:444` |
| `de_wilcoxon_sig_jaccard` | **0.0** | `_higher(0.0, 0.0)` | `scales.py:453` |

(`_lower` en `scales.py:303-310` fija `anchor=0.0`; `_higher` en
`scales.py:295-302` fija `anchor=1.0`.)

Asi que hay **tres** escalas distintas en el paquete, y solo la segunda es la de
la pregunta:

| escala | `b` | `a` (el 1) | como se pide |
|---|---|---|---|
| `from_baseline` | baseline medido | constante de la politica (`1.0` / `0.0`, `scoring.py:189`, `:161`) | `--baseline-agg` |
| **replica (#276)** | **baseline medido** | **replica medida (`anchor.py`)** | `--anchor` (diagnostico) o `--real-bundle` (competencia) |
| escala nombrada | **constante** (tabla de arriba) | constante | `--scale low-random_high-1_v10` |

## A.4 🔴 ¿Nuestro suelo corresponde a ese `b`? **No.**

El `b` del marcador se construye en tres pasos, todos en `baseline.py`:

**(1) El perfil.** `generic_response_profile` (`baseline.py:106`), definido
literalmente en `baseline.py:114-120`:

> "Mean per-perturbation pseudobulk of `real`, over non-control perturbations.
> With `exclude_target_gene=True` (the default) gene `g`'s mean omits the
> perturbation that targets `g`::
>
> `profile[g] = mean over { p : p != control, label(p) != g } of pseudobulk[p][g]`"

Con peso igual por perturbacion, no por celula (`baseline.py:133-134`). El
resultado es UN vector constante unico (`baseline.py:121-122`).

**(2) La emision.** `_emit_scaled_resample` (`baseline.py:312`), cuyo docstring
en `baseline.py:315` es "Resample control-pool row positions by group and scale
their genes". El remuestreo, `baseline.py:333`:

```
draws = [rng.integers(0, control_X.shape[0], size=size) for size in sizes]
```

y el escalado, `baseline.py:344`:

```
out.data = (out.data * scale[out.indices]).astype(np.float32)
```

donde `scale` es la razon `r = profile / ctrl` calculada en `_emission_scale`
(`baseline.py:272`). La propiedad que persigue, `baseline.py:372-374`:

> "In counts space, conditional on the fixed control pool, the construction has
> `E[group counts mean] = r * ctrl = profile` exactly on every reachable gene
> before float32 rounding."

**(3) La semilla.** `seed: int = 0` en `build_baseline_prediction`
(`baseline.py:362`), expuesta como `--seed` con default 0 en el `--help` de
`baseline` ("random seed for dispersed emission [default: 0]").

### El contraste con nuestro suelo

Nuestro `suelo_floor.predict_floor` (`modelo/src/suelo_floor.py:14-24`):

```
chosen = rng.integers(0, n_control, size=N_PER_LABEL)
out[label] = control_pool_X[chosen]
```

**Mismo remuestreo con reemplazo del pool de control. Sin el escalado.** Es el
paso (2) de `baseline.py` con `scale = 1` y sin el paso (1).

| | nuestro suelo | el `b` del marcador |
|---|---|---|
| que emite | el control, tal cual | el control **escalado** a la respuesta promedio |
| remuestreo | con reemplazo (`suelo_floor.py:22`) | con reemplazo (`baseline.py:333`) |
| escalado gen a gen | **ninguno** | `r = profile/ctrl` (`baseline.py:344`) |
| perfil de referencia | ninguno | media de pseudobulk por perturbacion (`baseline.py:114-120`) |
| que acota | emitir el control | la **TRIVIALIDAD** de la metrica |

Y el marcador es explicito en que su `b` **no** es un suelo alcanzable,
`baseline.py:6-10`:

> "It is an **ORACLE comparator, not a floor a submission could reach**: the
> profile is averaged from the evaluated real perturbations, so it is
> transductive and no real submission has access to that data. What it bounds is
> a metric's TRIVIALITY -- a score at or below it is not evidence of
> target-specific skill."

Nuestro suelo es lo contrario: una emision que **cualquiera** puede producir
(el pool de control se reparte). Son dos objetos con proposito distinto, no dos
versiones del mismo.

### Entonces, ¿que ES nuestro suelo en el vocabulario del marcador?

El extremo 0 de la escala constante, para las metricas donde esa escala lo dice.
`scales.py:321-322`, comentando la base `1.0` de
`expr_mse_unbiased_capped_norm`:

> "**base 1.0: a submission emitting the control unchanged reads 1.0**, and as of
> #257 that is a PROPERTY on any panel whatever the REFERENCE's depth rather than
> a declaration."

y `scales.py:324-326` describe exactamente nuestra construccion:

> "An arm emitting an **INDEPENDENT draw of the control's cells per perturbation**
> is at r ~ 0.97 and reads 1.0 to within that"

Para `de_wilcoxon_lfc_nmae`, base `1.0`, `scales.py:371-373`:

> "base 1.0: the `n` normalization makes a prediction whose log2FC is exactly zero
> read 1.0 -- at lfc_hat = 0 the numerator IS the denominator -- so this end of
> the scale needs no reference to the evaluation data."

**No esta en el codigo** una declaracion de que "el control remuestreado" sea el
0 de las otras cuatro metricas; para `pds_cosine` la base 0.5 se describe como
"the random minimum" (`scales.py:315`), que es una afirmacion mas general.

## A.5 El cruce numerico que lo confirma

Nuestro suelo medido (`bucle_suelo.md`) contra las bases constantes de
`low-random_high-1_v10`. Esto no es una cita de codigo, es la comprobacion de que
nuestro suelo cae donde el codigo dice que cae:

| metrica | base constante de la escala | nuestro suelo, 6 pliegues | ¿coincide? |
|---|---|---|---|
| `de_wilcoxon_lfc_nmae` | 1.0 (`scales.py:418`) | 0.996965, 1.000821, 1.028846, 1.036018, 1.044314, 1.056808 | **si**, +0.3% a +5.7% |
| `pds_cosine` | 0.5 (`scales.py:428`) | 0.474432, 0.490530, 0.495265, 0.512311, 0.512311, 0.517992 | **si**, -5.1% a +3.6% |
| `expr_mse_unbiased_capped_norm` | 1.0 (`scales.py:369`) | 1.332759, 1.774874, 1.868728, 1.931565, 2.562172, 2.913227 | no: 1.3x a 2.9x por encima |
| `de_wilcoxon_direction_reach_raw` | 0.0 (`scales.py:444`) | 0.013563, 0.142889, 0.144540, 0.192408, 0.200448, 0.216033 | cerca, no exacto |
| `de_wilcoxon_sig_jaccard` | 0.0 (`scales.py:453`) | 0.0, 0.0, 0.181818, 0.242424, 0.424242, 0.454545 | no (es la tasa de union vacia -- ver `bucle_techo.md`) |
| `de_wilcoxon_direction_fidelity_yield_raw` | 0.5 (`scales.py:430`) | 0.0, 0.0, 0.0, 0.0, 0.000026, 0.009259 | **no**: nuestro suelo da ~0, la base dice 0.5 |

Las dos metricas cuyo 0 el codigo deriva de "emitir el control" / "log2FC cero"
-- `lfc_nmae` y `pds_cosine` -- son justo las dos donde nuestro suelo aterriza
sobre la constante. Eso es evidencia de que nuestro suelo **es** ese punto de
referencia, y no el baseline oraculo.

Las desviaciones tienen explicacion en el propio codigo: para
`expr_mse_unbiased_capped_norm`, `scales.py:339-343` dice que la base es "a policy
constant ... it does not move with the panel" y que una submission que emite el
control "still scores slightly below 0 here rather than exactly 0"; nuestro
desvio es mucho mayor que el 2.49% que ahi se caracteriza, lo que **no esta en el
codigo** explicado para un panel como IFNB.

---

# B. El camino nativo

## B.1 Las cuatro invocaciones

Transcritas de los `--help` de los subcomandos y de las rutas de `cli.py`. Con
los valores que usamos (`pert_col=target_gene`, `control=NT`, IFNB como real).

```
# 1. el b -- la emision de referencia que el marcador construye
cell-eval2 baseline \
    -ar /datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad \
    --preset vcc2026 --profile vcc2026 \
    --pert-col target_gene --control NT \
    --emit dispersed --seed 0 \
    -o BASEDIR
# escribe: baseline_agg.csv, baseline_results.csv, baseline_meta.json  (cli.py:810-818)

# 2. el r, Y el empaquetado de los dos extremos con la identidad que los liga
cell-eval2 prep-real-bundle \
    --real /datos/ref/convertido/Seurat_object_IFNB_Perturb_seq.h5ad \
    --baseline BASEDIR \
    --preset vcc2026 --profile vcc2026 \
    --pert-col target_gene --control NT \
    --anchor-base-seed 0 --anchor-splits 5 \
    -o BUNDLE
# escribe: manifest.json, baseline_agg.csv, baseline_meta.json, config.yaml,
#          anchor_agg.parquet, anchor_splits.parquet, anchor_meta.json

# 3. el u -- la submission, igual que ya lo corremos
cell-eval2 run \
    -ap PRED.h5ad -ar REAL.h5ad \
    --preset vcc2026 --profile vcc2026 \
    --pert-col target_gene --control NT \
    -o RUNDIR

# 4. el s
cell-eval2 score \
    --user-agg RUNDIR/agg_results.csv \
    --user-meta RUNDIR/run_meta.json \
    --real-bundle BUNDLE \
    -o scored.csv
```

**El orden importa**: `prep-real-bundle` consume el directorio de `baseline`
(`--baseline BASELINE`, "the pre-built baseline arm, scored as the PREDICTION
(the 0 end)"), asi que (1) precede a (2). `score` consume el bundle y el run, asi
que (2) y (3) preceden a (4).

### `--anchor` NO es el camino de competencia

Del `--help` de `score`, para `--anchor`:

> "anchor directory written by `run --anchor`. ... Adds the from_replicate column
> **as a DIAGNOSTIC -- it never takes avg_score**. The competition score is
> `--real-bundle`, which supplies both ends and checks the submission against
> them."

Y para `--real-bundle`:

> "a real bundle from `prep-real-bundle`: BOTH ends of the scale plus the identity
> that binds them. Supplies the baseline and the anchor, so it is mutually
> exclusive with `--baseline-agg`/`--baseline-meta`/`--anchor`/`--anchor-cache`.
> ... With a competition bundle, **from_replicate takes avg_score** and
> from_baseline's goes null."

Asi que `run --anchor --anchor-splits 5 --anchor-base-seed 0` existe y produce el
`r`, pero por si solo da una columna diagnostica. Para el `s` enrolado hay que
pasar por `prep-real-bundle`.

## B.2 El flujo en `cli.py`

| subcomando | despacho | que llama | que escribe |
|---|---|---|---|
| `run` | `cli.py:515` | `compute_metrics`; con `--anchor`, `cli.py:741`; con `--ceiling`, `cli.py:747-753` | `agg_results.csv`, `results.csv`, `run_meta.json`, `run_params.yaml`; `anchor_*`; `ceiling_*` |
| `prep-real-bundle` | `cli.py:774` | `build_real_bundle` (`cli.py:775`, `780-783`) | el bundle; imprime si es de **competencia** o **DIAGNOSTICO** (`cli.py:788-793`) |
| `baseline` | `cli.py:794` | `build_generic_baseline` (`cli.py:795`, `799-807`) | `baseline_agg.csv` (`cli.py:809-810`), `baseline_results.csv` (`811`), `baseline_meta.json` (`815-816`) |
| `score` | `cli.py:823` | `score_metrics`; exige al menos uno de `--baseline-agg` / `--scale` / `--real-bundle` (`cli.py:825-829`) | el CSV de `-o` |

Dos guardas que conviene conocer antes de correr el camino:

- **`cli.py:788-793` imprime si el bundle quedo DIAGNOSTICO.** El comentario dice
  por que esta impreso: "A miscomputed rule check would otherwise un-enrol every
  submission scored against this bundle, silently and for weeks; here it is
  visible on the first build." Si sale "wrote DIAGNOSTIC bundle", el `avg_score`
  no se enrola.
- **`cache_strict` es obligatorio.** `configs/vcc2026.yaml:8-11` lo declara
  "LOAD-BEARING, not hygiene: `score.expect_from_run_meta` refuses a
  metadata-only fingerprint ... A competition run built without this cannot be
  scored against a real bundle at all". El preset ya lo pone en `true`, y
  **nuestras corridas lo cumplen**: `run_meta.json` de `run_A549` trae
  `source_fingerprint_strict = True`.

## B.3 Que de lo nuestro queda SOBRANDO

| nuestro | ¿sobra? | motivo |
|---|---|---|
| `modelo/src/suelo_floor.py` | **si, por completo** | `baseline.py` construye la emision de referencia, y la construye DISTINTA (escalada al perfil). Nuestra emision no es el `b` de nadie. |
| `modelo/src/suelo_arnes.py` | **si** como productor del `b`. `build_fold_adatas`, el submuestreo a 400 y el ensamblado de h5ad no tienen contraparte necesaria: `cell-eval2 baseline` toma el h5ad real directo. | El unico valor que queda es el lector de IFNB (`read_rows_csr` / `load_line_block`), util solo porque nosotros partimos por linea celular -- cosa que el camino nativo no hace. |
| `modelo/src/bucle_suelo.py` | **si** como productor del `b` | lo reemplazan las invocaciones 1 y 2 |
| `score_vs_self` (`suelo_arnes.py:229-246`) | **si, sin reservas** | puntuaba el agregado contra si mismo, lo que da `avg_score = 0.0` por construccion y un `score` que se niega cuando el denominador no es positivo (lo que nos paso en HAP1). El camino nativo puntua contra los DOS extremos reales. |
| `modelo/src/techo_arnes.py` | **si** como productor del `r` | `anchor.py` ya lo hace, con 5 particiones y la convencion oficial (ver `ceiling_convencion.md`) |
| `modelo/src/bucle_techo.py` | **parcialmente** | sobra como productor del `r`; **no** sobra como lo que resulto ser |
| `modelo/src/techo_reporte.py` + `contar_de_por_perturbacion` | **NO sobra** | las columnas `n_DE_mitadA`, `n_DE_mitadB`, `|union|` por perturbacion no las expone ningun artefacto nativo. `anchor_splits.parquet` guarda `value` y `n_perturbations` por (split, metrica) (`anchor.py:82-85`), no los conjuntos de DE. Fue lo que revelo que el `jac` del suelo es la tasa de union vacia. |

**Lo que el camino nativo NO da, y nosotros si medimos:** la particion por LINEA
CELULAR (6 pliegues leave-one-line-out). `baseline`, `prep-real-bundle` y
`anchor` toman UN h5ad real y lo tratan como una poblacion. Nuestra pregunta
-- ¿como varia el suelo y el techo entre lineas? -- no tiene bandera en el CLI.
**No esta en el codigo** ninguna estratificacion por lote o por linea: la unica
clave de agrupacion en `_disjoint_halves` es `pert_col` (`ceiling.py:155`), y el
perfil del baseline es un vector constante unico (`baseline.py:121-122`).

## B.4 ¿Resuelve solo los dos huecos que ya sabiamos?

### Hueco 1: `lfc_nmae` usa `FULL_GATE_RAW`, otro estimador

**Si, lo resuelve solo.** No hay que hacer nada: `compute_replicate_anchor` lo
sustituye internamente.

- La constante: `FULL_GATE_RAW = "full_gate_raw"` (`anchor.py:67`).
- El motivo, `anchor.py:57-66`: "`de_lfc_nmae` omits a perturbation whose
  real-side significance gate holds fewer than `min_gate_size` genes. A HALF calls
  far fewer genes significant than the full data (measured: the median gate is
  1.83-3.36x smaller), so the uniform core averages over 65-79% of the
  perturbations the MEMBER averages over. `compute_lfc_nmae_reference` takes its
  gate and denominator from the FULL real table and only the numerator's two
  vectors from the halves, so its cohort equals the member's EXACTLY (verified on
  six real cell lines)."
- La mecanica: una sola tabla DE del real COMPLETO, calculada una vez y reusada
  por los cinco splits (`anchor.py:222-232`), y la sustitucion por split en
  `anchor.py:244-250` -- hecha insertando en el dict, no con `with_columns`,
  justamente para que la metrica no desaparezca en silencio.
- El estampado: cada metrica lleva su productor en el artefacto,
  `anchor.py:296-299` (`estimator = FULL_GATE_RAW` para los nombres de lfc_nmae,
  `SPLIT_HALF_RAW` para el resto).

Y el riesgo que esto cierra, `anchor.py:64-66`: "Anchoring a metric on a
different population than it is measured on is the mismatch `score` cannot
detect -- it sees only the two scalars." Es decir: si nosotros replicaramos el
techo a mano para `lfc_nmae` con el nucleo uniforme, `score` **no lo detectaria**.

### Hueco 2: de las 6 metricas de vcc2026, `ceiling.py` solo cubre `pds_cosine`

**Si, lo resuelve -- porque el camino nativo no usa `ceiling.py`.** El `r` lo da
`anchor.py`, que cubre **todas** las metricas seleccionadas y **se niega** si
alguna falta, en vez de devolver `NaN`. `anchor.py:256-264`:

```
absent = [m for m in expected
          if vals.get(m) is None or not math.isfinite(float(vals[m]))]
if absent:
    raise ValueError(
        f"split {i} (seed {seed}) produced no usable value for {len(absent)} "
        f"selected metric(s): {absent}. The anchor must cover every selected "
        "metric: a member missing from the anchor cannot be scored on this scale, "
        "and silently dropping it would change the competition average's "
        "denominator with nothing in the artifact saying so."
    )
```

Comparar con `ceiling.py:216`, que para una metrica fuera de `SB_METRICS`
simplemente escribe `float("nan")`. Son politicas opuestas: `ceiling.py` calla,
`anchor.py` levanta. La `SB_METRICS` de `ceiling.py` **no participa** del camino
nativo.

---

## Lo que no esta en el codigo

- **Si nuestro suelo (control remuestreado, sin escalar) es el 0 de las cuatro
  metricas restantes** de `low-random_high-1_v10`. El codigo lo deriva solo para
  `expr_mse_unbiased_capped_norm` (`scales.py:321-322`) y `de_wilcoxon_lfc_nmae`
  (`scales.py:371-373`). Para `pds_cosine` dice "the random minimum"
  (`scales.py:315`), que no es lo mismo que "emitir el control".
- **Por que nuestro `fid` de suelo da ~0 contra una base constante de 0.5**
  (`scales.py:430`). No hay en el codigo una caracterizacion de esa base para una
  emision de control.
- **Por que nuestro `expr_mse_unbiased_capped_norm` de suelo se va a 1.3-2.9x la
  base** cuando `scales.py:339-343` caracteriza un desvio de ~2.5% para una
  submission que emite el control. El codigo no cubre un panel como IFNB.
- **Ninguna estratificacion por linea celular o por lote** en `baseline.py`,
  `anchor.py` ni `ceiling.py`. La unica clave de agrupacion es `pert_col`
  (`ceiling.py:155`), y el perfil del baseline es un vector constante unico
  (`baseline.py:121-122`).
- **Si `--emit tile`** reproduce algo que nos sirva: el codigo lo declara
  "known-biased" (`baseline.py:411`) y "tile exists only to reproduce pre-fix
  numbers" (`baseline.py:417`), sin decir contra que numeros previos.

---

## Procedencia

- No se ejecuto nada del paquete salvo `--help` de `run`, `baseline`, `score` y
  `prep-real-bundle`, para transcribir banderas y defaults literalmente. No se
  llamo a `compute_ceiling`, `compute_replicate_anchor`, `build_generic_baseline`
  ni `build_real_bundle`.
- No se toco ningun arnes: `suelo_arnes.py`, `suelo_floor.py`, `bucle_suelo.py`,
  `techo_arnes.py`, `bucle_techo.py`, `techo_reporte.py`, `bucle_techo.md` y
  `bucle_suelo.md` quedan como estaban.
- Se leyo un artefacto propio ya existente, `/datos/modelo/out/run_A549/run_meta.json`,
  solo para comprobar `source_fingerprint_strict`.
- No se leyo nada bajo `~/vcc/benchmark/` ni
  `/datos/benchmark_VEDADO_techo_metricas/`.
- Sin envios, sin `vcc sample`, sin commits.
- Archivos de fuente leidos, todos bajo `~/cell-eval2/src/cell_eval2/`:
  `score.py`, `scoring.py`, `scales.py`, `baseline.py`, `anchor.py`,
  `real_bundle.py`, `ceiling.py`, `cli.py`, `competition.py`,
  `configs/vcc2026.yaml`.
