#!/usr/bin/env bash
# Runner del MAPA DE SENSIBILIDAD DE LA ENTREGA. Misma forma que las dos
# noches anteriores (`noche/runner.sh`, `noche/suelo_llamada.sh`):
#
#   * PASOS EN ORDEN, cada uno AISLADO: el que falla deja su traza en este log
#     y el runner SIGUE con el siguiente.
#   * SALIDA INCREMENTAL: cada celda del banco escribe su json al terminar ELLA
#     (y brazo a brazo dentro de la celda); el .md se reescribe DESPUES DE CADA
#     PASO, asi que siempre hay informe de lo que ya hay.
#   * REENTRANTE: la celda cuyo json existe se salta. Relanzar no repite.
#   * TOPE DE DISCO: aborta si /datos baja de MIN_FREE_GIB; cada celda lo mira
#     otra vez antes de empezar y borra sus intermedios grandes al cerrar.
#   * 🔴 SIN CORTE POR RELOJ: corre hasta terminar los tres brazos, tarde lo que
#     tarde. Lo UNICO que lo para es el tope de disco.
#
# Argumentos: runner.sh [plan]   (plan: ligero | medio | completo)
# No se usa: ~/vcc/benchmark/, /datos/benchmark_VEDADO_*, vcc sample, envios,
# commits. Nada se ajusta mirando el `s`.

set -u
PY=/home/ubuntu/mf/envs/vcc/bin/python
DIR=/home/ubuntu/vcc/modelo/src/sensib
OUT=/home/ubuntu/vcc/modelo/out
WORK=/datos/modelo/out/sensib
MIN_FREE_GIB=10
PLAN=${1:-medio}

cd /home/ubuntu/vcc || exit 1
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$WORK/celdas" "$WORK/tmp"

libre_gib() { df -BG --output=avail /datos | tail -1 | tr -dc '0-9'; }

exec 9>"$WORK/.runner.lock"
if ! flock -n 9; then echo "YA HAY OTRO runner.sh CORRIENDO. Salgo."; exit 1; fi

# ── los pasos de cada plan. "nombre|comando"
case "$PLAN" in
  ligero)
    PASOS=(
      "A de|$PY $DIR/a_bio.py de"
      "A reparto|$PY $DIR/a_bio.py reparto"
      "A resumen|$PY $DIR/a_bio.py resumen"
      "B cruz D|$PY $DIR/b_sens.py D cruz"
      "C D|$PY $DIR/c_union.py D"
    ) ;;
  medio)
    PASOS=(
      "A de|$PY $DIR/a_bio.py de"
      "A reparto|$PY $DIR/a_bio.py reparto"
      "A resumen|$PY $DIR/a_bio.py resumen"
      "C D|$PY $DIR/c_union.py D"
      "B cruz D|$PY $DIR/b_sens.py D cruz"
      "C E|$PY $DIR/c_union.py E"
      "C F|$PY $DIR/c_union.py F"
      "B centro E|$PY $DIR/b_sens.py E 9"
      "B centro F|$PY $DIR/b_sens.py F 9"
      "B centro A (validacion)|$PY $DIR/b_sens.py A 9"
    ) ;;
  completo)
    PASOS=(
      "A de|$PY $DIR/a_bio.py de"
      "A reparto|$PY $DIR/a_bio.py reparto"
      "A resumen|$PY $DIR/a_bio.py resumen"
      "B rejilla D|$PY $DIR/b_sens.py D todo"
      "C D|$PY $DIR/c_union.py D"
      "C E|$PY $DIR/c_union.py E"
      "C F|$PY $DIR/c_union.py F"
      "B cruz E|$PY $DIR/b_sens.py E cruz"
      "B cruz F|$PY $DIR/b_sens.py F cruz"
      "B cruz A|$PY $DIR/b_sens.py A cruz"
      "B centro B|$PY $DIR/b_sens.py B 9"
      "B centro C|$PY $DIR/b_sens.py C 9"
    ) ;;
  *) echo "plan desconocido: $PLAN"; exit 1 ;;
esac

echo "================================================================"
echo "MAPA DE SENSIBILIDAD DE LA ENTREGA   (plan: $PLAN)"
echo "  arranque      : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
echo "  corte por reloj: NINGUNO (solo para el tope de disco)"
echo "  pasos         : ${#PASOS[@]}"
echo "  /datos libre  : $(libre_gib) GiB   (aborta por debajo de ${MIN_FREE_GIB})"
echo "  celdas ya en disco: $(ls "$WORK/celdas" 2>/dev/null | grep -c '.json$')"
echo "  salida        : $OUT/sensibilidad_entrega.md (+ .json)"
echo "================================================================"
echo

RES_OK=(); RES_FALLO=(); RES_SIN=()
T0T=$(date +%s)

for i in "${!PASOS[@]}"; do
  IFS='|' read -r NOMBRE CMD <<< "${PASOS[$i]}"
  L=$(libre_gib)
  if [ "$L" -lt "$MIN_FREE_GIB" ]; then
    echo "################################################################"
    echo "ABORTO POR DISCO: ${L} GiB libres (< ${MIN_FREE_GIB}). Sin empezar:"
    for j in $(seq "$i" $((${#PASOS[@]} - 1))); do
      IFS='|' read -r N2 _ <<< "${PASOS[$j]}"; echo "    - $N2"; RES_SIN+=("$N2"); done
    echo "################################################################"
    break
  fi
  echo "================================================================"
  echo ">>> $NOMBRE   INICIO : $(date -u '+%Y-%m-%d %H:%M:%S') UTC   (/datos ${L} GiB)"
  echo "    $CMD"
  T0=$(date +%s)
  set +e
  $CMD
  EC=$?
  set -e; set +u; set -u
  D=$(( $(date +%s) - T0 ))
  if [ "$EC" -eq 0 ]; then RES_OK+=("$NOMBRE"); EST="OK";
  else RES_FALLO+=("$NOMBRE(exit=$EC)"); EST="FALLO exit=$EC <-- traza arriba"; fi
  # ── el informe se reescribe tras CADA paso: siempre hay .md de lo que hay
  set +e
  "$PY" "$DIR/reporte.py"
  RC=$?
  set -e
  echo "<<< $NOMBRE   FIN    : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
  echo "    duracion       : ${D} s  ($((D / 60)) min)"
  echo "    estado         : $EST"
  echo "    celdas en disco: $(ls "$WORK/celdas" | grep -c '.json$')"
  echo "    reporte        : exit=$RC  $( [ -f "$OUT/sensibilidad_entrega.md" ] && \
       echo "$(wc -l < "$OUT/sensibilidad_entrega.md") lineas" || echo "NO EXISTE")"
  echo "    /datos libre   : $(libre_gib) GiB"
  echo
done

DT=$(( $(date +%s) - T0T ))
echo "================================================================"
echo "RESUMEN"
echo "  fin           : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
echo "  duracion total: ${DT} s  ($((DT / 60)) min)"
echo "  OK            : ${RES_OK[*]:-(ninguno)}"
echo "  FALLARON      : ${RES_FALLO[*]:-(ninguno)}"
echo "  SIN EMPEZAR   : ${RES_SIN[*]:-(ninguno)}"
echo "  celdas        : $(ls "$WORK/celdas" | grep -c '.json$')"
ls "$WORK/celdas" | sed 's/^/    /'
for F in sensibilidad_entrega.md sensibilidad_entrega.json; do
  [ -f "$OUT/$F" ] && echo "  $OUT/$F ($(stat -c%s "$OUT/$F") bytes)" || echo "  $OUT/$F (no existe)"
done
echo "  /datos libre  : $(libre_gib) GiB"
echo "================================================================"
