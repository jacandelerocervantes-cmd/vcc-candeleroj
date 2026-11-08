#!/usr/bin/env bash
# Runner del SUELO DE REPRODUCIBILIDAD DE LA LLAMADA. Igual forma que la noche
# del 4-5 (modelo/src/noche/runner.sh):
#
#   * UN PASO POR LINEA, corridos EN ORDEN y cada uno AISLADO: el que falla
#     escribe su traza en este mismo log y el runner SIGUE con el siguiente.
#     Un fallo en HT29 no se lleva A549..HAP1 ni impide 2025.
#   * SALIDA INCREMENTAL: cada sonda escribe su json en sets/<tag>_j<j>.json AL
#     TERMINAR ELLA (no al final), y cada paso es reentrante: si el json ya
#     esta, se salta. Matar el runner y relanzarlo no repite trabajo hecho.
#   * TOPE DE DISCO: aborta si /datos baja de MIN_FREE_GIB (tambien lo mira
#     `comun.guarda_disco()` dentro de cada sonda).
#   * RELEVO: si sigue vivo el `sl_sondas.py v2025` adjunto a la sesion, se le
#     DEJA acabar antes de arrancar el paso de 2025 (tope ESPERA_MAX_S, luego
#     se le mata y la j en vuelo se rehace). Asi no hay dos procesos sobre el
#     mismo tmp.
#   * INTEGRIDAD: antes de 2025 y antes del analisis, todo sets/*.json tiene que
#     parsear; el que no parsea (corte a media escritura) se APARTA a
#     sets/corruptos_<ts>/ y se vuelve a calcular. No se borra nada.
#
# No se usa: ~/vcc/benchmark/, /datos/benchmark_VEDADO_*, vcc sample, envios,
# commits. Nada se ajusta mirando el `s`.

set -u
PY=/home/ubuntu/mf/envs/vcc/bin/python
DIR=/home/ubuntu/vcc/modelo/src/suelo_llamada
OUT=/home/ubuntu/vcc/modelo/out
WORK=/datos/modelo/out/suelo_llamada
SETS=$WORK/sets
MIN_FREE_GIB=10
ESPERA_MAX_S=3600

cd /home/ubuntu/vcc || exit 1
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$SETS" "$WORK/logs"

libre_gib() { df -BG --output=avail /datos | tail -1 | tr -dc '0-9'; }

# --- un solo runner a la vez
exec 9>"$WORK/.runner.lock"
if ! flock -n 9; then
  echo "YA HAY OTRO suelo_llamada.sh CORRIENDO (lock $WORK/.runner.lock). Salgo."
  exit 1
fi

# --- pasos: NOMBRE|argv de sl_sondas.py, o ANALISIS / REPORTE
PASOS=(
  "SONDAS A549|jiang A549"
  "SONDAS BXPC3|jiang BXPC3"
  "SONDAS HAP1|jiang HAP1"
  "SONDAS HT29|jiang HT29"
  "SONDAS K562|jiang K562"
  "SONDAS MCF7|jiang MCF7"
  "SONDAS v2025|v2025"
  "ANALISIS|ANALISIS"
  "REPORTE|REPORTE"
)

espera_relevo() {
  # Espera al sl_sondas.py v2025 que quedo adjunto a la sesion, si vive.
  local pids t0 n
  pids=$(pgrep -f "sl_sondas.py v2025" | grep -v "^$$\$" || true)
  if [ -z "$pids" ]; then echo "    relevo: no hay otro sl_sondas.py v2025 vivo"; return 0; fi
  echo "    relevo: vive(n) PID(s) [$(echo $pids | tr '\n' ' ')] sobre v2025; se les DEJA acabar (tope ${ESPERA_MAX_S}s)"
  t0=$(date +%s)
  while :; do
    pids=$(pgrep -f "sl_sondas.py v2025" || true)
    [ -z "$pids" ] && break
    n=$(( $(date +%s) - t0 ))
    if [ "$n" -ge "$ESPERA_MAX_S" ]; then
      echo "    relevo: tope de ${ESPERA_MAX_S}s alcanzado; se MATA [$(echo $pids | tr '\n' ' ')] y se rehace la j en vuelo"
      kill $pids 2>/dev/null; sleep 15; kill -9 $pids 2>/dev/null || true
      sleep 5; break
    fi
    if [ $((n % 300)) -lt 30 ]; then
      echo "    relevo: esperando ($((n / 60)) min) — sets: $(ls "$SETS" | grep -c '^v2025_j') de 8 de v2025"
    fi
    sleep 30
  done
  echo "    relevo: libre a las $(date -u '+%H:%M:%S') UTC; sigo yo"
}

integridad() {
  # Todo sets/*.json tiene que parsear. El que no, se APARTA (no se borra).
  "$PY" - "$SETS" <<'PYFIN'
import json, os, shutil, sys, time
sets = sys.argv[1]
malos = []
for f in sorted(os.listdir(sets)):
    if not f.endswith(".json"):
        continue
    p = os.path.join(sets, f)
    try:
        json.load(open(p))
    except Exception as e:
        malos.append((f, repr(e)[:120]))
if malos:
    d = os.path.join(sets, "corruptos_" + time.strftime("%Y%m%d_%H%M%S"))
    os.makedirs(d, exist_ok=True)
    for f, e in malos:
        shutil.move(os.path.join(sets, f), os.path.join(d, f))
        print(f"    INTEGRIDAD: {f} no parsea ({e}) -> apartado a {d}, se recalcula")
else:
    print(f"    integridad: los {len([f for f in os.listdir(sets) if f.endswith('.json')])} json de sets parsean")
PYFIN
}

echo "================================================================"
echo "SUELO DE LA LLAMADA  (sondas por mitades -> analisis -> reporte)"
echo "  arranque      : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
echo "  pasos         : ${#PASOS[@]} en orden (6 lineas Jiang, 2025, analisis, reporte)"
echo "  /datos libre  : $(libre_gib) GiB   (aborta por debajo de ${MIN_FREE_GIB})"
echo "  sondas        : $DIR/sl_sondas.py   -> $SETS/<tag>_j<j>.json"
echo "  analisis      : $DIR/sl_analisis.py -> $OUT/suelo_llamada.json"
echo "  reporte       : $DIR/sl_reporte.py  -> $OUT/suelo_llamada.md"
echo "  python        : $PY"
echo "  hecho ya      : $(ls "$SETS" | grep -c '_j[0-7].json$') sondas en disco de 56 (7 tags x 8)"
echo "================================================================"
echo

RES_OK=(); RES_FALLO=(); RES_SIN_EMPEZAR=()
T_TOTAL=$(date +%s)

for i in "${!PASOS[@]}"; do
  IFS='|' read -r NOMBRE CMD <<< "${PASOS[$i]}"

  # --- tope de disco: aborta, no sigue
  L=$(libre_gib)
  if [ "$L" -lt "$MIN_FREE_GIB" ]; then
    echo "################################################################"
    echo "ABORTO POR DISCO: /datos tiene ${L} GiB libres (< ${MIN_FREE_GIB})."
    echo "  no se arranca $NOMBRE ni los siguientes:"
    for j in $(seq "$i" $((${#PASOS[@]} - 1))); do
      IFS='|' read -r N2 _ <<< "${PASOS[$j]}"
      echo "    - $N2"; RES_SIN_EMPEZAR+=("$N2")
    done
    echo "################################################################"
    break
  fi

  echo "================================================================"
  echo ">>> $NOMBRE   INICIO : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
  echo "    /datos libre   : ${L} GiB"
  T0=$(date +%s)

  # --- aislamiento: el fallo de un paso no se lleva lo anterior
  set +e
  case "$CMD" in
    "v2025")
      espera_relevo
      integridad
      "$PY" "$DIR/sl_sondas.py" v2025
      EC=$?
      ;;
    "ANALISIS")
      integridad
      "$PY" "$DIR/sl_analisis.py"
      EC=$?
      ;;
    "REPORTE")
      "$PY" "$DIR/sl_reporte.py"
      EC=$?
      ;;
    *)
      "$PY" "$DIR/sl_sondas.py" $CMD
      EC=$?
      ;;
  esac
  set -e
  set +u; set -u

  T1=$(date +%s); D=$((T1 - T0))
  if [ "$EC" -eq 0 ]; then RES_OK+=("$NOMBRE"); ESTADO="OK"
  else RES_FALLO+=("$NOMBRE(exit=$EC)"); ESTADO="FALLO exit=$EC  <-- la traza esta arriba en este mismo log"; fi
  echo
  echo "<<< $NOMBRE   FIN    : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
  echo "    duracion       : ${D} s  ($((D / 60)) min)"
  echo "    estado         : $ESTADO"
  case "$CMD" in
    ANALISIS) [ -f "$OUT/suelo_llamada.json" ] \
        && echo "    salida         : $OUT/suelo_llamada.json ($(stat -c%s "$OUT/suelo_llamada.json") bytes, $(date -u -r "$OUT/suelo_llamada.json" '+%H:%M:%S'))" \
        || echo "    salida         : $OUT/suelo_llamada.json NO EXISTE" ;;
    REPORTE)  [ -f "$OUT/suelo_llamada.md" ] \
        && echo "    salida         : $OUT/suelo_llamada.md ($(wc -l < "$OUT/suelo_llamada.md") lineas, $(date -u -r "$OUT/suelo_llamada.md" '+%H:%M:%S'))" \
        || echo "    salida         : $OUT/suelo_llamada.md NO EXISTE" ;;
    *) TAG=$(echo "$CMD" | awk '{print $NF}')
       echo "    sondas en disco: $(ls "$SETS" | grep -c "^${TAG}_j[0-7].json$") de 8 de $TAG" ;;
  esac
  echo "    /datos libre   : $(libre_gib) GiB"
  echo
done

DT=$(( $(date +%s) - T_TOTAL ))
echo "================================================================"
echo "RESUMEN"
echo "  fin           : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
echo "  duracion total: ${DT} s  ($((DT / 60)) min)"
echo "  OK            : ${RES_OK[*]:-(ninguno)}"
echo "  FALLARON      : ${RES_FALLO[*]:-(ninguno)}"
echo "  SIN EMPEZAR   : ${RES_SIN_EMPEZAR[*]:-(ninguno)}"
echo "  /datos libre  : $(libre_gib) GiB"
echo "  sondas por tag:"
for T in A549 BXPC3 HAP1 HT29 K562 MCF7 v2025; do
  echo "    $T  $(ls "$SETS" | grep -c "^${T}_j[0-7].json$")/8"
done
for F in suelo_llamada.json suelo_llamada.md; do
  if [ -f "$OUT/$F" ]; then echo "  $OUT/$F  ($(stat -c%s "$OUT/$F") bytes, $(date -u -r "$OUT/$F" '+%H:%M:%S') UTC)"
  else echo "  $OUT/$F  (no existe)"; fi
done
echo "================================================================"
