# Censo de disco y de respaldo — qué hay, dónde vive, y qué está (o no) en S3

Esto es un censo, no un diagnóstico de qué arreglar. No se movió, copió, borró ni
renombró nada. No se tocó ninguna unidad de systemd. No se tocó `X` de ningún archivo.

## 1 · Árbol

### `~/vcc` completo hasta profundidad 3 — tamaño total **862 M**

| ruta de 1er nivel | tamaño | notas |
|---|---|---|
| `PLAN.md` | 9.4K | — |
| `mf.sh` | 128,128,881 (122 M) | archivo suelto en la raíz de `~/vcc`, no un script trivial |
| `corpus/` | 816 K | `inventario.parquet` (819,766) + `particion.md` (6,503) |
| `out/` | 739 M | 150+ archivos sueltos + `techo_metricas_work/` (la mayor parte del peso) |
| `registro/` | 24 K | `CONSTANTES_A_VERIFICAR.md`, `limitaciones.md` |
| `src/` | 312 K | scripts `.py`/`.sh`, más `__pycache__/` |
| `workflow/` | 24 K | `vigilantes/{spot_watch.sh, sync_out.sh, sync_out.pid, sync_out.log}` |
| `.claude/` | 4.0K (vacío) | ver directorios vacíos abajo |

**Directorios vacíos** (cualquier profundidad, dentro de `~/vcc`):
- `~/vcc/.claude`
- `~/vcc/out/techo_metricas_work/user_i1_cpm1e6_default_attempt`
- `~/vcc/out/techo_metricas_work/user_i2_invert_agg_default_attempt`
- `~/vcc/out/techo_metricas_work/user_i3_log1p_default_attempt`
- `~/vcc/out/techo_metricas_work/user_i4_log_default_attempt`

**Enlaces simbólicos dentro de `~/vcc`: ninguno** (búsqueda recursiva completa, no solo
hasta profundidad 3).

El árbol completo de `~/vcc` hasta profundidad 3 (tipo, ruta, tamaño en bytes, fecha de
modificación) tiene >160 líneas — no se repite aquí completo por espacio; el detalle
íntegro quedó capturado en esta sesión y puede regenerarse con:
```
find ~/vcc -maxdepth 3 -printf '%y %p %s %TY-%Tm-%Td %TH:%TM %l\n'
```

### `/datos` hasta profundidad 2 — tamaño total ~**242 G** usados

| ruta | tamaño | fecha de modificación (directorio) |
|---|---|---|
| `/datos/ref/` | 194 G | 2026-10-27 18:30 |
| `/datos/ref/convertido/` | (5 `.h5ad`) | 2026-10-27 18:32 |
| `/datos/ref/feng/` | (3 `.tsv.gz`) | 2026-10-27 18:32 |
| `/datos/ref/nadig/` | (2 `.h5ad`) | 2026-10-27 18:34 |
| `/datos/ref/orion/` | (2 `.parquet`) | 2026-10-27 18:35 |
| `/datos/ref/otros/` | — | 2026-10-27 18:40 |
| `/datos/ref/replogle/` | (6 `.h5ad`) | 2026-10-27 18:52 |
| `/datos/ref/tfact/` | (1 `.csv`) | 2026-10-27 18:48 |
| `/datos/ref/dl.log` | 1,388 | 2026-10-27 07:11 |
| `/datos/reto/` | 2.5 G | 2026-10-27 18:25 (`final/`, `val/`) |
| `/datos/vcc2025/` | 33 G | 2026-10-27 18:50 (`test/`, `train/`, `validation/`) + `gene_names.csv` |
| **`/datos/_scratch_techo_metricas/`** | **12 G** | 2026-10-28 06:02–06:04 — **`predA_base.h5ad` (5.97 G) + `realB.h5ad` (5.97 G) + `split_log.json`** |
| `/datos/lost+found/` | — | **no accesible, permiso denegado** (artefacto normal de ext4, no investigado más) |

**Directorios vacíos en `/datos` (hasta profundidad 4, excluyendo `lost+found` por
permiso)**: ninguno encontrado. **Enlaces simbólicos (profundidad 2)**: ninguno
encontrado.

`_scratch_techo_metricas/` es el único directorio de `/datos` **creado después** de la
carga inicial (27 oct) — es la partición train/test que usó `cell-eval2` para las
corridas de `techo_metricas_work`, y es la fuente de `de_real.parquet` citada en
`nreal_marcador.md`.

## 2 · Volúmenes

| dispositivo | modelo (vía `/sys/block/*/device/model`) | tamaño | punto de montaje | qué vive ahí |
|---|---|---|---|---|
| `/dev/nvme0n1p1` | **Amazon Elastic Block Store** | ~30 GiB | `/` (raíz) | **todo `~/vcc`** (PLAN.md, mf.sh, corpus/, out/, registro/, src/, workflow/) |
| `/dev/nvme1n1` | **Amazon Elastic Block Store** | ~300 GiB | `/datos` | `/datos/ref`, `/datos/reto`, `/datos/vcc2025`, `/datos/_scratch_techo_metricas` |
| `/dev/nvme2n1` | Amazon EC2 NVMe Instance Storage (**efímero real**) | ~280 GiB | **sin montar, sin usar** | nada |

**Corrección directa a la premisa de `PLAN.md` §6** ("`/datos` (NVMe local) = solo
datasets... nada irreemplazable aquí", "el volumen raíz sobrevive"): **medido ahora, es
al revés de lo que ese texto sugiere sobre cuál es efímero.** Según
`/sys/block/*/device/model` y el *block-device-mapping* de IMDSv2 (`root`→`sda1`,
`ebs1`→`sdf` ≡ `nvme1n1`), **tanto el volumen raíz (`nvme0n1`, donde vive `~/vcc`) como
`/datos` (`nvme1n1`) son EBS**, no almacenamiento de instancia — los dos sobreviven un
`stop`/`start` de la instancia (interrupción spot), no solo la raíz. El `r5ad.2xlarge`
sí trae un disco NVMe de instancia genuinamente efímero (`nvme2n1`, ~280 GiB), pero
**no está montado ni usado** — nada depende de él, así que tampoco hay nada en riesgo
por ese lado. No se pudo confirmar `DeleteOnTermination` por volumen (el rol IAM de
esta instancia no tiene permiso para `ec2:DescribeInstances`), pero esa bandera solo
importa en una terminación, no en el `stop` que `PLAN.md` dice temer.

## 3 · Respaldo — `syncout.service`, textual

```
$ systemctl cat syncout.service
# /etc/systemd/system/syncout.service
[Unit]
Description=Sync vcc out a S3

[Service]
User=ubuntu
ExecStart=/bin/bash -c "while true; do aws s3 sync /home/ubuntu/vcc/ s3://vc-challenge-data/vcc/ --exclude "corpus/*" --exclude "out/*.npz" --no-progress >/dev/null 2>&1; sleep 120; done"
Restart=always

[Install]
WantedBy=multi-user.target
```
`enabled` + `active (running)` desde las 17:31:36 de hoy (PID 696 actualmente, ciclo
cada 120 s). Fuente: `/home/ubuntu/vcc/` completo → destino: `s3://vc-challenge-data/vcc/`.

**Contraste contra el árbol del punto 1 — qué queda DENTRO y qué queda FUERA:**

| ruta | ¿dentro de `syncout.service`? |
|---|---|
| `~/vcc/PLAN.md`, `~/vcc/mf.sh` | **sí** |
| `~/vcc/out/` (todo excepto `.npz`) | **sí** — incluye `techo_metricas_work/` completo |
| `~/vcc/registro/`, `~/vcc/src/`, `~/vcc/workflow/` | **sí** |
| **`~/vcc/corpus/inventario.parquet`, `~/vcc/corpus/particion.md`** | **NO** — excluidos explícitamente (`--exclude "corpus/*"`) |
| `~/vcc/out/vcc2025_knockdown_raw_arrays.npz` | **NO** — excluido (`--exclude "out/*.npz"`) |
| `~/vcc/out/vcc2025_rejilla_raw_arrays_test_train.npz` | **NO** — excluido |
| `~/vcc/out/vcc2025_rejilla_raw_arrays_validation.npz` | **NO** — excluido |
| **cualquier cosa bajo `/datos`** | **NO** — `syncout.service` no referencia `/datos` en absoluto; no hay ninguna unidad ni cron que lo sincronice |

**Esta es la segunda lista** que pide el punto 3: `corpus/*` completo, los 3 `.npz` de
`out/`, y **todo `/datos`** (incluyendo los 12 G de `_scratch_techo_metricas`, creados
después de cualquier carga inicial a S3) quedan fuera de la sincronización automática.

**Hallazgo lateral, no buscado**: `~/vcc/workflow/vigilantes/sync_out.sh` y
`spot_watch.sh` son un **mecanismo de respaldo distinto y anterior**, con otro destino
(`s3://vc-challenge-data/out` — raíz del bucket, no `vcc/out/`, y **sin exclusión de
`.npz`**). Su `sync_out.pid` guarda PIDs (2982/2986) que **ya no corresponden a ningún
proceso vivo** (`ps -p` no encuentra ninguno), y su log (`sync_out.log`) tiene una sola
entrada, de ayer (2026-10-27 20:03). Pese a eso, **`s3://vc-challenge-data/out/`
contiene una copia de `t3_inventario.md` fechada hoy 2026-10-28T21:14:12Z** — 40
segundos antes de que `syncout.service` subiera la suya a `vcc/out/` (21:14:52Z). **No
se pudo identificar, desde este host, qué proceso escribe ahí ahora**: no aparece en
`ps aux`, no hay cron ni otra unidad de systemd con "sync" o "s3" en el nombre, no hay
gancho en `.bashrc`/`.profile`, y el bucket no tiene reglas de replicación ni
notificación configuradas. Se reporta la observación tal cual — resolverla excede
"reportar, no arreglar".

## 4 · S3 — por prefijo de primer nivel

| prefijo | objetos | tamaño total | fecha más reciente |
|---|---|---|---|
| `datos/` | 55 | 244.98 G | 2026-10-27 16:47:09 — **no se actualiza**, es la carga inicial |
| `out/` | 426 | 772.9 M | 2026-10-28 21:14:12 — mecanismo activo no identificado (ver §3) |
| `vcc/` | 470 | 885.7 M | 2026-10-28 21:16:29 — `syncout.service`, al día |
| `virtual cell/` | 4 | 18.9 K | 2026-10-27 06:07:32 — **legado abandonado**: `PLAN.md` + 2 scripts de `src/`, de la convención de nombres previa a que el proyecto se asentara en `vcc/` |

**¿Hay algo en S3 que ya no esté en local?** `datos/` en S3 (244.98 G, carga del
27-oct) corresponde a los archivos originales de `/datos/ref`, `/datos/reto` y
`/datos/vcc2025`, que **siguen en local** — no se detectó nada huérfano solo-en-S3 en
ese prefijo. `virtual cell/` sí es contenido que ya no se actualiza ahí, pero sus 4
archivos (`PLAN.md`, `CONSTANTES_A_VERIFICAR.md`, `cpm_stats.py`, `pseudobulk.py`)
siguen existiendo en local con los mismos nombres — no es una pérdida, es una ruta S3
duplicada y abandonada.

**¿Hay algo en local que no esté en S3?** Sí, lo que lista el contraste del punto 3:
`corpus/*`, los 3 `.npz`, y todo `/datos/_scratch_techo_metricas/` (12 G, los dos
`.h5ad` de 5.97 G cada uno) — ninguno de los tres aparece en ningún prefijo de S3.

## 5 · Duplicados y colisiones (mismo nombre base, más de un sitio en `~/vcc`)

**8 nombres de archivo duplicados — los 8 dentro de `out/techo_metricas_work/`**, y
los 8 son el mismo patrón: la grilla de experimentos (`baseline/`, `bundle/`,
`user_base/`, `user_i1_cpm1e6/`...`user_i5_1p0/`, `i6/user_i6_s*_f*/`) produce,
**por diseño del arnés**, un archivo con el mismo nombre (`agg_results.csv`,
`evaluation_results.json`, `metric_aggregation.csv`, `results.csv`, `run_meta.json`,
`run_params.yaml`) en cada subcarpeta — **no es una colisión accidental**: cada copia
tiene un tamaño ligeramente distinto (confirmado archivo por archivo), consistente con
ser la salida de una corrida distinta, no la misma corrida repetida.

**Única excepción que amerita nota aparte — `baseline_agg.csv`:**

| ruta | tamaño | fecha |
|---|---|---|
| `out/techo_metricas_work/baseline/baseline_agg.csv` | 1,229 | 2026-10-27 23:09 |
| `out/techo_metricas_work/bundle/baseline_agg.csv` | 1,229 | 2026-10-28 00:55 |

**Mismo tamaño, pero `diff` confirma que el CONTENIDO es distinto.** No se investigó
más allá (el encargo es reportar, no reconciliar), pero es la única de las 8 colisiones
donde "mismo nombre, tamaño igual" no implica "mismo contenido" — queda señalado.

**Ningún nombre de reporte `.md` está duplicado en ningún otro sitio de `~/vcc`.**
**Ningún nombre de archivo `.parquet` está duplicado** (solo existen dos:
`out/t13_replicas.parquet` y `corpus/inventario.parquet`, nombres distintos).

## 6 · Entregables de T3 — dónde quedaron exactamente

| archivo | ruta local | tamaño | fecha (local) | volumen | ¿dentro de la sincronización? |
|---|---|---|---|---|---|
| `inventario.parquet` | `/home/ubuntu/vcc/corpus/inventario.parquet` | 819,766 | 2026-10-28 21:07 | raíz (EBS, `nvme0n1p1`) | **NO** por `syncout.service` (excluido) — **SÍ** en `s3://vc-challenge-data/vcc/corpus/` (21:16:29Z) por una subida manual puntual hecha en esta sesión, fuera del mecanismo automático |
| `particion.md` | `/home/ubuntu/vcc/corpus/particion.md` | 6,503 | 2026-10-28 21:12 | raíz (EBS) | igual que arriba — **NO** automático, **SÍ** por la subida manual puntual |
| `t3_inventario.md` | `/home/ubuntu/vcc/out/t3_inventario.md` | 17,644 | 2026-10-28 21:14 | raíz (EBS) | **SÍ** — `syncout.service` lo subió a `vcc/out/t3_inventario.md` (21:14:52Z); **también** apareció en `out/t3_inventario.md` (bucket raíz, 21:14:12Z) por el mecanismo no identificado del §3 |

**Los dos primeros no tienen respaldo automático en curso.** La subida manual que los
puso en S3 hoy fue un evento único, no un mecanismo que se repita si estos archivos
cambian de nuevo.

## RIESGO — qué costaría trabajo rehacer y hoy está fuera del respaldo, en orden

1. **`/datos/_scratch_techo_metricas/predA_base.h5ad` + `realB.h5ad`** (5.97 G cada
   uno, 12 G total). No están bajo `~/vcc` (ningún mecanismo de sync los alcanzaría
   aunque se corrigieran los `--exclude`), y no están en ningún prefijo de S3. Rehacerlos
   significa repetir toda la partición train/test que usó `cell-eval2` para
   `techo_metricas_work` — el costo más alto de los tres, y la base de
   `nreal_marcador.md`/`nreal_transfiere.md`.
2. **`~/vcc/corpus/inventario.parquet` y `corpus/particion.md`** — excluidos a propósito
   de `syncout.service`. Hoy están en S3 solo por la subida manual de esta sesión; la
   próxima vez que cambien, no se subirán solos. Rehacerlos significa repetir el censo
   completo de T3 (~4 min de cómputo, más el diseño de esquema y las 7 consultas).
3. **Los 3 `.npz` de `~/vcc/out/`** (`vcc2025_knockdown_raw_arrays.npz`,
   `vcc2025_rejilla_raw_arrays_test_train.npz`,
   `vcc2025_rejilla_raw_arrays_validation.npz`, ~16.5 M total) — excluidos por
   `syncout.service`. Costo de rehacer: moderado, re-correr los scripts de
   `knockdown_grid*.py`/`b_techo_i6_ruido_gen.py` que los generaron.

No se lista `/datos/ref`, `/datos/reto` ni `/datos/vcc2025` (241 G) porque ya están
respaldados en `s3://vc-challenge-data/datos/` desde la carga inicial, y porque además
viven en un volumen EBS persistente (§2) — el riesgo ahí es bajo en ambos frentes, no
solo en uno.
