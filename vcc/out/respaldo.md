# Cierre de huecos de respaldo — qué se arregló, qué se determinó

Viene de `~/vcc/out/estructura.md`. Esta vez sí se modifica el sistema (una unidad de
systemd, el contenido de dos prefijos de S3) — todo lo demás (archivos locales, los
cuatro directorios vacíos de `techo_metricas_work`, el proceso de Claude Code y
cualquier otro cómputo en esta máquina) se dejó intacto, sin excepción.

## 1 · El segundo escritor a S3 — no identificado, nada detenido

Búsqueda exhaustiva en el host: `ps -ef` completo, `pgrep -af bash`, `tmux list-panes
-a`, `ls -l /proc/*/cwd`, `systemctl list-units --all | grep -i sync|s3|vcc`,
`crontab -l` (ubuntu y root), `/etc/cron.d`, `/etc/cron.daily`, `/etc/cron.hourly`.

**Resultado: un solo proceso de sincronización a S3 en todo el host — PID 696 (ahora
reemplazado, ver §2), el de `syncout.service`, destino `vcc/`.** Ningún otro proceso,
cron, unidad de systemd, ni gancho de shell hace `aws s3 sync`. `workflow/vigilantes/
sync_out.sh` y `spot_watch.sh` (destino `s3://vc-challenge-data/out`, el prefijo
antiguo) no tienen ningún proceso vivo correspondiente a los PIDs de su `.pid`
(2982/2986 — `ps -p` no encuentra ninguno).

**Por la regla de seguridad del encargo: no se identificó un bucle vivo en este host
que escriba al prefijo antiguo, así que no se detuvo nada ahí ni se tocó el `.pid`.**
Que `out/t3_inventario.md` apareciera a las 21:14:12Z — 40 s antes que la copia de
`syncout.service` en `vcc/out/` (21:14:52Z) — sigue sin explicación desde este host.
Puede ser un proceso en otra instancia con el mismo rol IAM, o un mecanismo que ya no
estaba corriendo para el momento en que se buscó. **Esto lo resuelve la consola de
AWS, no esta máquina.**

**Acción tomada de todos modos, pedida a mitad de esta tarea**: se subió un marcador
`LEGADO.txt` a los dos prefijos obsoletos, diciendo cuál es el vigente:
- `s3://vc-challenge-data/out/LEGADO.txt` (1,417 bytes, subido 21:35:31Z)
- `s3://vc-challenge-data/virtual cell/LEGADO.txt` (753 bytes, subido 21:35:32Z)

Ambos marcadores dicen explícitamente que el prefijo vigente es `vcc/`, y el de `out/`
dice también, sin adornos, que su escritor no fue identificado — el marcador no afirma
que el prefijo haya dejado de recibir escritura, solo que ya no es el destino que este
proyecto mantiene a propósito.

## 2 · `syncout.service` — las dos exclusiones, quitadas; el tercer directorio, añadido

### a) Unidad original, copiada a `~/vcc/registro/syncout.service.bak` antes de editar

```ini
[Unit]
Description=Sync vcc out a S3
[Service]
User=ubuntu
ExecStart=/bin/bash -c "while true; do aws s3 sync /home/ubuntu/vcc/ s3://vc-challenge-data/vcc/ --exclude "corpus/*" --exclude "out/*.npz" --no-progress >/dev/null 2>&1; sleep 120; done"
Restart=always
[Install]
WantedBy=multi-user.target
```

### b) + c) Unidad nueva (`/etc/systemd/system/syncout.service`)

```ini
[Unit]
Description=Sync vcc out a S3
[Service]
User=ubuntu
ExecStart=/bin/bash -c "while true; do aws s3 sync /home/ubuntu/vcc/ s3://vc-challenge-data/vcc/ --no-progress >/dev/null 2>&1; aws s3 sync /datos/_scratch_techo_metricas/ s3://vc-challenge-data/scratch_techo/ --no-progress >/dev/null 2>&1; sleep 120; done"
Restart=always
[Install]
WantedBy=multi-user.target
```

Sin las dos exclusiones; con una segunda línea de `aws s3 sync` dentro del mismo ciclo
de 120 s, para `/datos/_scratch_techo_metricas/` → `s3://vc-challenge-data/scratch_techo/`.

### d) `daemon-reload` + `restart`, confirmado activo

```
$ sudo systemctl daemon-reload && sudo systemctl restart syncout.service
$ systemctl status syncout.service
● syncout.service - Sync vcc out a S3
     Active: active (running) since Wed 2026-10-28 21:34:41 UTC
   Main PID: 12793 (bash)
     CGroup: ├─12793 /bin/bash -c "while true; do aws s3 sync ... ; aws s3 sync /datos/_scratch_techo_metricas/ ... ; sleep 120; done"
             └─12831 aws s3 sync /datos/_scratch_techo_metricas/ s3://vc-challenge-data/scratch_techo/ --no-progress
```

### e) Verificación real — no el acuse, el contenido en S3 después del primer ciclo

```
$ aws s3 ls s3://vc-challenge-data/scratch_techo/ --recursive
2026-10-28 21:34:45 5973295832 scratch_techo/predA_base.h5ad
2026-10-28 21:34:45 5971459718 scratch_techo/realB.h5ad
2026-10-28 21:34:45       8227 scratch_techo/split_log.json

$ ls -la /datos/_scratch_techo_metricas/
-rw-r--r-- 1 root   root   5973295832 predA_base.h5ad
-rw-r--r-- 1 root   root   5971459718 realB.h5ad
-rw-r--r-- 1 ubuntu ubuntu       8227 split_log.json

$ aws s3 ls s3://vc-challenge-data/vcc/corpus/
2026-10-28 21:16:29     819766 inventario.parquet
2026-10-28 21:16:29       6503 particion.md

$ aws s3 ls s3://vc-challenge-data/vcc/out/ | grep npz
2026-10-28 21:34:43    1884816 vcc2025_knockdown_raw_arrays.npz
2026-10-28 21:34:43   12113673 vcc2025_rejilla_raw_arrays_test_train.npz
2026-10-28 21:34:43    2503056 vcc2025_rejilla_raw_arrays_validation.npz
```

**Los tres tamaños de `scratch_techo/` coinciden byte a byte con los locales.** Los
tres `.npz` ya están en `vcc/out/`. `corpus/` sigue mostrando la marca de tiempo de la
subida manual de la sesión anterior (21:16:29Z) — correcto: `aws s3 sync` compara
contenido/tamaño y no vuelve a subir lo que no cambió; el archivo ya está dentro del
alcance del servicio desde ahora, así que el próximo cambio real sí se subirá solo. La
transferencia completa de los 12 G tardó menos de 90 s dentro del primer ciclo — no
hizo falta ningún procedimiento especial más allá de correr dentro de la sesión con
tmux disponible.

## 3 · ¿Es reproducible la mitad B? — **Sí, de forma determinista y verificada contra el código**

`split_log.json` completo:

```json
{
  "seed": 42,
  "n_obs_total": 170846,
  "n_A_total": 85448,
  "n_B_total": 85398,
  "odd_group_rule": "ceil to A: n_A = ceil(n/2)",
  "per_group": { "...": "100 perturbaciones + non-targeting, cada una con n_total/n_A/n_B" }
}
```
(las 101 entradas de `per_group` — 100 genes + `non-targeting`, 38,176 células — están
íntegras en `/datos/_scratch_techo_metricas/split_log.json`, no se recortan aquí por
espacio; el total cuadra: suma de `n_A` por grupo = 85,448, suma de `n_B` = 85,398).

Esto por sí solo ya es una semilla + una regla de desempate declaradas — pero lo que
hace la reproducción **verificable, no solo plausible**, es que el código que las
consumió existe y es legible: `~/vcc/src/techo_metricas_split.py`. Ahí:
- `SEED = 42`, un único `rng = np.random.default_rng(SEED)`.
- Orden de consumo **determinista**: `for g in sorted(tg.categories):` — alfabético,
  no el orden de aparición en el archivo ni ningún orden dependiente de hash.
- Por grupo: `perm = rng.permutation(idx)`; `n_a = (n + 1) // 2` — exactamente la regla
  que `split_log.json` declara (`"ceil to A"`).
- Fuente fija: `/datos/vcc2025/test/adata_Test.h5ad` — ya en S3 desde la carga inicial
  (`datos/` prefix, ver `estructura.md` §4).

**En una frase: la partición es determinista y reproducible desde cero — semilla fija,
orden de iteración fijo (alfabético sobre las categorías), regla de desempate
declarada y verificada línea por línea contra `techo_metricas_split.py`, y la fuente
(`adata_Test.h5ad`) ya está respaldada — así que `predA_base.h5ad`/`realB.h5ad` (y por
extensión `de_real.parquet` y toda la serie de `n_real`) se pueden regenerar
exactamente, y el respaldo que se subió en el punto 2c es conveniencia (ahorra volver a
correr la partición y sus ~80 s de I/O), no la única copia que exista de esos datos.**

(Advertencia menor, no bloqueante: la reproducción exacta asume la misma versión de
NumPy — `default_rng`/PCG64 es el generador documentado y estable de NumPy desde la
versión 1.17, no se esperan cambios de algoritmo entre versiones menores, pero no se
verificó explícitamente qué versión de NumPy corrió la partición original.)

## 4 · Los dos `baseline_agg.csv` — resuelto: es `baseline/`, con cita

### Diff completo (archivo de 1,229 bytes, cabe entero)

```diff
4c4
< mean,0.21014601054104973,0.08638446954761758,0.9466774556293243,0.014945130405859739,0.002685715804851687,0.0025865474680519546,0.0025934944352463407,0.9656622754206716,1.0,0.483232323232323
---
> mean,0.21014601054104973,0.08638446954761758,0.9466774556293238,0.014945130405859739,0.002685715804851687,0.0025865474680519546,0.0025934944352463407,0.9656622754206716,1.0,0.483232323232323
```
(`<` = `baseline/`, `>` = `bundle/`). **La única diferencia es la mean de
`de_wilcoxon_lfc_nmae` en el dígito 16: `...293243` vs `...293238` — ruido de punto
flotante (orden de suma no asociativo, probablemente por un paso de serialización
distinto), no dos resultados distintos.** Las otras 9 columnas y las otras 6 filas
(`count`, `null_count`, `std`, `min`, `max`, `median`) son **idénticas, carácter por
carácter**.

### Qué produjo cada carpeta, y en qué orden

| | `baseline/` | `bundle/` |
|---|---|---|
| `created_utc` (de su propio meta) | **2026-10-27T23:09:47Z** | 2026-10-27T23:10:27Z (40 s después) |
| contenido del directorio | `baseline_agg.csv`, `baseline_meta.json`, `run_meta.json`, `run_params.yaml` | `baseline_agg.csv`, `baseline_meta.json`, `anchor_agg.parquet`, `anchor_meta.json`, `anchor_splits.parquet`, `config.yaml`, `manifest.json` |
| `baseline_meta.json` → `source` | `"/datos/_scratch_techo_metricas/realB.h5ad"` (el dato real) | `"bundle"` (autorreferencial) |
| `config_digest` | `61d594ef...96067338` (nota: ver abajo) | **idéntico**: `61d594ef...96067338` |
| `source_fingerprint` | `249391849a...50efbc37979f` | **idéntico**: `249391849a...50efbc37979f` |
| qué es, según `techo_estado.md` (línea 15–17) | "Baseline (respuesta media) construido" | "Bundle real (baseline-leg + ancla de 5 splits, perfil vcc2026) ... (rule_digest de competencia confirmado en el manifest)" |

`bundle/` tiene el mismo `config_digest` y el mismo `source_fingerprint` que
`baseline/` — es la **misma computación de baseline**, empaquetada 40 segundos después
junto con el ancla de 5 splits (`derived_seeds` de `manifest.json`: 5 semillas
derivadas, `base_seed: 0`) para el chequeo de identidad de la competencia
(`--real-bundle`). No es un segundo baseline independiente; es el mismo, reexportado.

### Cuál suministró la `b` de las cifras reportadas

**`baseline/baseline_agg.csv`**, citado textualmente dos veces en los reportes ya
cerrados:

- `~/vcc/out/techo_metricas.md`, líneas 188–190: *"...se uso consistentemente el camino
  diagnostico `score --baseline-agg baseline/baseline_agg.csv --allow-config-mismatch`
  (sin `--anchor` ni `--real-bundle`), que reporta `from_baseline`..."* — usado para las
  9 de 10 condiciones diagnósticas (i1–i5) que no pasaron el gate de identidad del
  `--real-bundle`.
- `~/vcc/out/techo_estado.md`, línea 77: la misma cita, mismo comando.

**No hay ambigüedad real que resolver más allá de la ya resuelta**: el diff mostró que
ambos archivos son, para efectos prácticos, el mismo número (16 dígitos de acuerdo), y
la cita textual en los dos reportes ya cerrados confirma que el camino que efectivamente
se corrió usó `baseline/baseline_agg.csv`, no `bundle/`. Ninguna `s` reportada queda
pendiente de este punto.

## PARA

**Qué quedó cerrado:**

1. **No hay un segundo bucle de sincronización vivo en este host** (búsqueda
   exhaustiva, punto 1). El escritor del prefijo antiguo de S3 queda sin identificar —
   resultado reportado, no resuelto, tal como pedía el encargo. Ambos prefijos
   obsoletos (`out/`, `virtual cell/`) ahora tienen un `LEGADO.txt` señalando el
   prefijo vigente.
2. **`syncout.service` ya no excluye `corpus/*` ni `out/*.npz`, y ahora también
   sincroniza `/datos/_scratch_techo_metricas/`** — verificado con `aws s3 ls`, no solo
   con el estado del servicio. Unidad original respaldada en
   `~/vcc/registro/syncout.service.bak`.
3. **La mitad B es reproducible, no solo respaldada** — semilla + orden + regla
   verificados contra `techo_metricas_split.py`. El respaldo de 12 G ya subido es
   conveniencia.
4. **La `b` de las cifras de `techo_metricas.md`/`techo_estado.md` es
   `baseline/baseline_agg.csv`**, con cita textual — no queda ninguna `s` reportada sin
   baseline identificado.

**Qué sigue abierto, con su dependencia:**

- **La identidad del segundo escritor a S3** depende de la consola de AWS (listar
  instancias con el rol `ec2-vcc`, o revisar CloudTrail para la clave
  `out/t3_inventario.md` alrededor de las 21:14:12Z) — no de nada que se pueda hacer
  desde este host con el rol IAM actual.
- **La versión exacta de NumPy** que corrió `techo_metricas_split.py` originalmente no
  se verificó — si algún día `default_rng`/PCG64 cambiara de algoritmo entre versiones
  mayores de NumPy (no ha pasado históricamente), la reproducción exacta dependería de
  fijar esa versión también. Bajo riesgo, no verificado.
