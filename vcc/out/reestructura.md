# Reestructuración del bucket a dos prefijos — un prefijo por sistema de archivos

Regla fija: **un prefijo por sistema de archivos, y el nombre del prefijo es el punto
de montaje** — `datos/` = `/datos`, `vcc/` = `/home/ubuntu/vcc`. Nada en la raíz del
bucket que no sea un punto de montaje.

Bucket con versionado activo (todo delete es reversible), pero cada borrado de esta
tarea se hizo **solo después de verificar clave por clave**, nunca por conteo de
prefijo.

## 1 · La unidad, antes y después

**Antes** (`~/vcc/registro/syncout.service.bak2`, copiada antes de tocar nada):

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

**Después** (`systemctl cat syncout.service`, pegada ANTES de arrancar el servicio):

```ini
# Regla fija: UN PREFIJO POR SISTEMA DE ARCHIVOS, Y EL NOMBRE DEL PREFIJO ES EL PUNTO
# DE MONTAJE. datos/ = /datos ; vcc/ = /home/ubuntu/vcc. Nada en la raiz del bucket
# que no sea un punto de montaje. (fijado 2026-10-28, ver registro/limitaciones.md)
[Unit]
Description=Sync vcc out a S3
[Service]
User=ubuntu
ExecStart=/bin/bash -c "while true; do aws s3 sync /home/ubuntu/vcc/ s3://vc-challenge-data/vcc/ --no-progress >/dev/null 2>&1; aws s3 sync /datos/_scratch_techo_metricas/ s3://vc-challenge-data/datos/_scratch_techo_metricas/ --no-progress >/dev/null 2>&1; sleep 120; done"
Restart=always
[Install]
WantedBy=multi-user.target
```

Único cambio real: el destino del segundo `aws s3 sync` pasó de
`s3://vc-challenge-data/scratch_techo/` a
`s3://vc-challenge-data/datos/_scratch_techo_metricas/` — más la regla, como
comentario, antes de `[Unit]`. Secuencia real: `stop` → backup → editar → `cat` (pegado
arriba) → `daemon-reload` → `start` → confirmado `active (running)` (PID 14510).

## 2 · `scratch_techo/` movido del lado del servidor

```
$ aws s3 mv s3://vc-challenge-data/scratch_techo/ s3://vc-challenge-data/datos/_scratch_techo_metricas/ --recursive
move: s3://vc-challenge-data/scratch_techo/split_log.json to s3://vc-challenge-data/datos/_scratch_techo_metricas/split_log.json
move: s3://vc-challenge-data/scratch_techo/predA_base.h5ad to s3://vc-challenge-data/datos/_scratch_techo_metricas/predA_base.h5ad
move: s3://vc-challenge-data/scratch_techo/realB.h5ad to s3://vc-challenge-data/datos/_scratch_techo_metricas/realB.h5ad

$ aws s3 ls s3://vc-challenge-data/datos/_scratch_techo_metricas/
2026-10-28 21:57:04 5973295832 predA_base.h5ad
2026-10-28 21:57:08 5971459718 realB.h5ad
2026-10-28 21:57:04       8227 split_log.json

$ aws s3 ls s3://vc-challenge-data/scratch_techo/
(sin salida — vacío)
```

Las tres claves, tamaño exacto, coinciden con lo esperado (5,973,295,832 /
5,971,459,718 / 8,227) y con los archivos locales. `scratch_techo/` no lista nada.

## 3 · `out/` contra `vcc/out/` — diferencia por clave

`out/`: 429 objetos. `vcc/out/`: 431 objetos (ya incluía, antes de esta tarea, los 4
archivos de `~/vcc/out/mapa_*.{txt,md}` que no existían cuando `out/` se dejó de
escribir).

### Solo en `out/` (2 claves)

```
LEGADO.txt                  1417
techo_metricas_DRAFT.md      153
```

### Solo en `vcc/out/` (4 claves)

```
mapa.md                  79344
mapa_datos.txt            4974
mapa_s3.txt              79001
mapa_vcc.txt             46896
```

### En ambas pero con tamaño distinto

```
(ninguna — las 427 claves comunes tienen tamaño idéntico en las dos copias)
```

### Qué se hizo con las 2 exclusivas de `out/`

- **`techo_metricas_DRAFT.md`** (153 bytes): el borrador previo a `techo_metricas.md`
  (`# Techo de las metricas ... BORRADOR EN CONSTRUCCION`), sin equivalente local bajo
  ningún nombre. **Copiado del lado del servidor a `vcc/out/techo_metricas_DRAFT.md`
  antes de borrar nada** (verificado: 153 bytes en el destino).
- **`LEGADO.txt`** (1,417 bytes): el marcador de la tarea anterior. Por instrucción
  explícita de esta tarea, **no se copia** — se retira junto con el resto de `out/`.

Con las 427 claves comunes verificadas de tamaño idéntico y la única clave exclusiva
real ya copiada, se borró `out/` completo:

```
$ aws s3 rm s3://vc-challenge-data/out/ --recursive
delete: s3://vc-challenge-data/out/... (429 líneas de "delete:")

$ aws s3 ls s3://vc-challenge-data/out/
(sin salida — vacío)
```

## 4 · `virtual cell/` — verificación de las 5 claves, **2 no cuadran, NADA SE BORRA**

| clave en `virtual cell/` | tamaño en S3 | equivalente local | tamaño local | ¿cuadra? |
|---|---|---|---|---|
| `PLAN.md` | 9,637 | `/home/ubuntu/vcc/PLAN.md` | 9,637 | ✅ |
| `registro/CONSTANTES_A_VERIFICAR.md` | 5,333 | `/home/ubuntu/vcc/registro/CONSTANTES_A_VERIFICAR.md` | **7,699** | ❌ **tamaño distinto** |
| `src/cpm_stats.py` | 2,961 | `/home/ubuntu/vcc/src/cpm_stats.py` | 2,961 | ✅ |
| `src/pseudobulk.py` | 1,033 | `/home/ubuntu/vcc/src/pseudobulk.py` | 1,033 | ✅ |
| `LEGADO.txt` | 753 | — | **no existe ningún archivo local con ese nombre** | ❌ **sin equivalente local** |

**2 de 5 no cuadran.** Por la regla de seguridad ("si alguna verificación no cuadra, se
para ahí y se reporta; mejor cinco prefijos que una clave perdida"): **no se borró
ninguna de las 5 claves de `virtual cell/`.** El prefijo sigue existiendo, intacto,
exactamente como estaba.

(Nota aparte, no una excusa: `CONSTANTES_A_VERIFICAR.md` no cuadra porque el archivo
local creció después del 27 de octubre — ya se había señalado en `mapa.md` §6 como una
copia congelada que diverge. `LEGADO.txt` no cuadra porque nunca fue un archivo local —
se subió directo a S3 en la tarea anterior.)

## 5 · Verificación final

```
$ aws s3 ls s3://vc-challenge-data/
                           PRE datos/
                           PRE vcc/
                           PRE virtual cell/

$ aws s3 ls s3://vc-challenge-data/ --recursive | wc -l
545
```

**No son dos líneas, son tres.** `virtual cell/` quedó porque su verificación falló
(§4) — resultado correcto de la regla de seguridad, no un olvido.

**Conteo, reconciliado exactamente** (no solo "contrastado", cuadrado número a
número): snapshot previo (`mapa_s3.txt`) = 968 objetos.

| cambio | objetos |
|---|---|
| `out/` borrado | −429 |
| `scratch_techo/` movido (3 claves fuera) | −3 |
| `datos/_scratch_techo_metricas/` (las mismas 3, ahora ahí) | +3 |
| `vcc/out/`: +4 (`mapa.md`, `mapa_vcc.txt`, `mapa_datos.txt`, `mapa_s3.txt`, ya subidos por el ciclo de `syncout.service` antes de esta tarea) | +4 |
| `vcc/out/techo_metricas_DRAFT.md` (copiado en §3) | +1 |
| `vcc/registro/syncout.service.bak2` (creado en esta tarea, §1) | +1 |
| `virtual cell/` | 0 (sin tocar) |
| **total** | **968 − 429 + 6 = 545** |

Coincide exacto con el conteo medido. Desglose final por prefijo, confirmado por
`list-objects-v2`: `datos/` = 58, `vcc/` = 482, `virtual cell/` = 5 → 58+482+5=545.

## 6 · La regla, escrita donde no se pierda

- **En la unidad**: comentario de 3 líneas antes de `[Unit]` (pegado en §1, "Después").
- **En `~/vcc/registro/limitaciones.md`**: entrada nueva, "9. Regla fija de estructura
  de S3 — un prefijo por sistema de archivos", fecha 2026-10-28, con el estado abierto
  de `virtual cell/` anotado ahí también.

## Addendum — `out/` se repobló sola, ~3 minutos después de confirmarse vacía

Al cerrar el reporte (antes de limpiar los archivos temporales) se volvió a listar la
raíz del bucket por costumbre, no porque se esperara nada: **`out/` había reaparecido**,
con **432 objetos**, `LastModified` **2026-10-28T22:00:34Z** — la verificación de §5
(vacía) se hizo a las ~21:57Z. El mismo escritor no identificado de `estructura.md` y
`respaldo.md` sigue activo, y esta vez quedó agarrado en el acto de recrear contenido
recién borrado, no solo de escribir algo nuevo.

**Evidencia de que lee del disco local, no de otra copia en S3**: de los 432 objetos,
*solo uno* no tiene equivalente en `vcc/out/` — **`reestructura.md`**, este mismo
archivo, que se acababa de escribir en `~/vcc/out/` segundos antes de que apareciera
en `out/` (bucket raíz). Al mismo tiempo, **`techo_metricas_DRAFT.md`** (la clave que
se copió del lado del servidor, S3→S3, en §3 — nunca se escribió en disco local) **no
apareció** en la recreación de `out/`. Eso aísla el mecanismo: algo sincroniza
`~/vcc/out` (el directorio local) hacia `s3://vc-challenge-data/out`, igual que hacía
`workflow/vigilantes/sync_out.sh` en su diseño original — pero no desde ningún proceso
visible en `ps -ef` de este host (se volvió a comprobar: solo el PID 14510 de
`syncout.service`, que no referencia `out/` en absoluto).

**No se borró `out/` una segunda vez.** Borrar contra un escritor activo y no
identificado es perseguir un blanco en movimiento, no una limpieza — y el encargo no
pidió repetir el borrado, pidió reportar el estado. Esto queda como el hallazgo más
importante de la tarea: **la reestructura de `out/` no es durable mientras ese segundo
escritor siga activo**, venga de donde venga.

## PARA

**Qué quedó cerrado:**

1. La unidad ya no puede recrear `scratch_techo/` — el segundo `sync` apunta a
   `datos/_scratch_techo_metricas/`, verificado textual antes de arrancar.
2. `scratch_techo/` se movió del lado del servidor, verificado byte a byte; el prefijo
   quedó vacío.
3. `out/` se borró completo, pero solo después de verificar clave por clave (no por
   conteo) que las 427 claves comunes coinciden en tamaño y de copiar la única clave
   real exclusiva (`techo_metricas_DRAFT.md`) antes de borrar nada.
4. La regla queda escrita en dos sitios que no se pierden con el próximo reinicio ni
   con la próxima sesión: la unidad misma y `limitaciones.md`.

**Qué queda abierto, con su dependencia:**

- **`out/` no se quedó borrado — se repobló sola en minutos (ver addendum).** Esto
  depende enteramente de identificar al segundo escritor, que sigue sin aparecer en
  `ps -ef`, cron, ni systemd de este host. Hasta que se identifique y se detenga (o se
  confirme que es otro host que hay que apagar desde la consola de AWS), **cualquier
  borrado de `out/` es temporal**, no una limpieza definitiva.
- **La raíz del bucket tiene tres prefijos, no dos** (y, por el addendum, `out/`
  volvió a aparecer — son efectivamente cuatro ahora mismo, no tres). `virtual cell/` (5 claves,
  19.3 KiB) no se retiró porque 2 de sus 5 claves no superaron la verificación de
  seguridad: `registro/CONSTANTES_A_VERIFICAR.md` difiere en tamaño del archivo local
  actual (la copia en S3 quedó congelada en su versión del 27 de octubre), y
  `LEGADO.txt` nunca tuvo un archivo local que lo respalde. Resolverlo depende de una
  decisión que esta tarea no estaba autorizada a tomar sola: o se acepta borrar
  `registro/CONSTANTES_A_VERIFICAR.md` de `virtual cell/` sabiendo que la versión ahí
  es vieja y se pierde (el bucket tiene versionado, así que seguiría siendo
  recuperable), o se decide qué hacer con `LEGADO.txt` (borrarlo sin más, ya que no
  representa ningún archivo local). Ninguna de las dos decisiones es mecánica —
  requieren que alguien diga explícitamente que está bien perder esa cadena de
  verificación para esas 2 claves.
