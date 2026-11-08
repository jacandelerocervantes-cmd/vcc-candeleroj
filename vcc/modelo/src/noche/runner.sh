#!/usr/bin/env bash
# Runner del trabajo nocturno. Corre los bloques EN ORDEN, cada uno aislado.
#
#  * un bloque que falla NO se lleva lo anterior: se escribe su traza al log y
#    el runner SIGUE con el siguiente
#  * SALIDA INCREMENTAL: cada bloque escribe su .json y su .md AL TERMINAR EL,
#    no al final de todo (B4c reescribe su .md tras CADA semilla)
#  * TOPE DE DISCO: el runner ABORTA si /datos baja de MIN_FREE_GIB
#  * COMPUERTA: pasada DEADLINE_UTC no se ARRANCA otro bloque; el que este
#    corriendo se DEJA terminar, y el log dice que quedo sin empezar
#
# No se usa: ~/vcc/benchmark/, /datos/benchmark_VEDADO_*, vcc sample, envios,
# commits. Nada se ajusta mirando el `s`.

set -u
PY=/home/ubuntu/mf/envs/vcc/bin/python
DIR=/home/ubuntu/vcc/modelo/src/noche
OUT=/home/ubuntu/vcc/modelo/out
MIN_FREE_GIB=10
DEADLINE_UTC=$(date -u -d "today 11:30" +%s)

cd /home/ubuntu/vcc || exit 1
export PYTHONDONTWRITEBYTECODE=1

# bloque:script:salidas
BLOQUES=(
  "B1:b1_venn_ed.py:venn_ed.md"
  "B2:b2_sobredispersion.py:sobredispersion.md"
  "B3:b3_posicion.py:posicion_emisiones.md"
  "B4a:b4a_inventario.py:vcc2025_inventario.md"
  "B4b:b4b_anclajes.py:vcc2025_anclajes.md"
  "B4c:b4c_suelo.py:vcc2025_suelo.md"
)

libre_gib() { df -BG --output=avail /datos | tail -1 | tr -dc '0-9'; }

echo "================================================================"
echo "TRABAJO NOCTURNO"
echo "  arranque      : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
echo "  compuerta     : $(date -u -d @"$DEADLINE_UTC" '+%Y-%m-%d %H:%M:%S') UTC"
echo "                  (pasada esa hora no se arranca otro bloque;"
echo "                   el que este corriendo se deja terminar)"
echo "  /datos libre  : $(libre_gib) GiB   (aborta por debajo de ${MIN_FREE_GIB})"
echo "  bloques       : ${#BLOQUES[@]} en orden: B1, B2, B3, B4a, B4b, B4c"
echo "  python        : $PY"
echo "================================================================"
echo

RES_OK=(); RES_FALLO=(); RES_SIN_EMPEZAR=()
T_TOTAL=$(date +%s)

for i in "${!BLOQUES[@]}"; do
  IFS=: read -r NOMBRE SCRIPT SALIDA <<< "${BLOQUES[$i]}"

  # --- tope de disco: aborta, no sigue
  L=$(libre_gib)
  if [ "$L" -lt "$MIN_FREE_GIB" ]; then
    echo "################################################################"
    echo "ABORTO POR DISCO: /datos tiene ${L} GiB libres (< ${MIN_FREE_GIB})."
    echo "  no se arranca $NOMBRE ni los siguientes."
    for j in $(seq "$i" $((${#BLOQUES[@]} - 1))); do
      IFS=: read -r N2 _ _ <<< "${BLOQUES[$j]}"
      RES_SIN_EMPEZAR+=("$N2")
    done
    echo "################################################################"
    break
  fi

  # --- compuerta horaria: no arranca otro bloque
  AHORA=$(date +%s)
  if [ "$AHORA" -ge "$DEADLINE_UTC" ]; then
    echo "################################################################"
    echo "COMPUERTA 11:30 UTC alcanzada a las $(date -u '+%H:%M:%S')."
    echo "  NO se arranca $NOMBRE. Bloques que quedan SIN EMPEZAR:"
    for j in $(seq "$i" $((${#BLOQUES[@]} - 1))); do
      IFS=: read -r N2 S2 O2 <<< "${BLOQUES[$j]}"
      echo "    - $N2 ($S2) -> habria escrito $OUT/$O2"
      RES_SIN_EMPEZAR+=("$N2")
    done
    echo "################################################################"
    break
  fi

  echo "================================================================"
  echo ">>> $NOMBRE   INICIO : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
  echo "    script         : $DIR/$SCRIPT"
  echo "    salida esperada: $OUT/$SALIDA"
  echo "    /datos libre   : ${L} GiB"
  echo "================================================================"
  T0=$(date +%s)

  # --- aislamiento: el fallo de un bloque no se lleva lo anterior.
  #     La traza de python va al log por stderr; el exitcode se registra.
  set +e
  "$PY" "$DIR/$SCRIPT"
  EC=$?
  set -e
  set +u; set -u

  T1=$(date +%s); D=$((T1 - T0))
  if [ "$EC" -eq 0 ]; then
    RES_OK+=("$NOMBRE")
    ESTADO="OK"
  else
    RES_FALLO+=("$NOMBRE(exit=$EC)")
    ESTADO="FALLO exit=$EC  <-- la traza esta arriba en este mismo log"
  fi
  echo
  echo "<<< $NOMBRE   FIN    : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
  echo "    duracion       : ${D} s  ($((D / 60)) min)"
  echo "    estado         : $ESTADO"
  if [ -f "$OUT/$SALIDA" ]; then
    echo "    salida         : $OUT/$SALIDA  ($(wc -l < "$OUT/$SALIDA") lineas, $(date -u -r "$OUT/$SALIDA" '+%H:%M:%S'))"
  else
    echo "    salida         : $OUT/$SALIDA  NO EXISTE"
  fi
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
echo "  salidas en disco:"
for b in "${BLOQUES[@]}"; do
  IFS=: read -r N2 _ O2 <<< "$b"
  if [ -f "$OUT/$O2" ]; then
    echo "    $N2  $OUT/$O2  ($(wc -l < "$OUT/$O2") lineas)"
  else
    echo "    $N2  $OUT/$O2  (no existe)"
  fi
done
echo "================================================================"
