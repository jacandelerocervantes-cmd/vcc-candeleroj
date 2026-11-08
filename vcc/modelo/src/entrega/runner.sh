#!/usr/bin/env bash
# Runner de LA ENTREGA. Misma forma que las noches anteriores:
#   * PASOS EN ORDEN, cada uno AISLADO: el que falla deja su traza en el log y
#     el runner SIGUE con el siguiente.
#   * SALIDA INCREMENTAL: cada paso escribe su json/.md al terminar EL; la
#     emision escribe su progreso cada 25 perturbaciones y es REANUDABLE.
#   * ABORTA si /datos baja de 10 GiB (y la emision lo vuelve a mirar sola).
#   * Marca de inicio, fin y duracion por paso.
#   * SIN corte por reloj.
# 🔴 NO ENVIA NADA: no hay `vcc submit` ni `vcc sample` en ningun paso.

set -u
PY=/home/ubuntu/mf/envs/vcc/bin/python
DIR=/home/ubuntu/vcc/modelo/src/entrega
OUT=/home/ubuntu/vcc/modelo/out
WORK=/datos/modelo/out/entrega
H5=$WORK/prediction.h5ad
VCC=$WORK/prediction.vcc
RETO=/datos/reto/final
MIN_FREE_GIB=10

cd /home/ubuntu/vcc || exit 1
export PYTHONDONTWRITEBYTECODE=1
mkdir -p "$WORK"
libre_gib() { df -BG --output=avail /datos | tail -1 | tr -dc '0-9'; }

exec 9>"$WORK/.runner.lock"
if ! flock -n 9; then echo "YA HAY OTRO runner de entrega CORRIENDO. Salgo."; exit 1; fi

PASOS=(
  "1 AJUSTE (la funcion en VCC 2025)|$PY $DIR/ajuste.py"
  "2 EMISION (D/E/F -> h5ad)|$PY $DIR/emite.py"
  "3 CONTRATO (comprobacion propia)|$PY $DIR/contrato.py"
  "4 PREP --dry-run|dryrun"
  "5 PREP (escribe el .vcc)|prep"
)

echo "================================================================"
echo "LA ENTREGA: entrenar -> emitir -> contrato -> empaquetar"
echo "  arranque      : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
echo "  pasos         : ${#PASOS[@]} (sin corte por reloj)"
echo "  /datos libre  : $(libre_gib) GiB   (aborta por debajo de ${MIN_FREE_GIB})"
echo "  h5ad          : $H5"
echo "  vcc           : $VCC"
echo "  reporte       : $OUT/entrega.md"
echo "  🔴 ningun paso envia nada"
echo "================================================================"
echo

OK=(); FALLO=(); SIN=()
T0T=$(date +%s)
for i in "${!PASOS[@]}"; do
  IFS='|' read -r NOMBRE CMD <<< "${PASOS[$i]}"
  L=$(libre_gib)
  if [ "$L" -lt "$MIN_FREE_GIB" ]; then
    echo "################################################################"
    echo "ABORTO POR DISCO: ${L} GiB libres (< ${MIN_FREE_GIB}). Sin empezar:"
    for k in $(seq "$i" $((${#PASOS[@]} - 1))); do
      IFS='|' read -r N2 _ <<< "${PASOS[$k]}"; echo "    - $N2"; SIN+=("$N2"); done
    echo "################################################################"
    break
  fi
  echo "================================================================"
  echo ">>> $NOMBRE   INICIO : $(date -u '+%Y-%m-%d %H:%M:%S') UTC  (/datos ${L} GiB)"
  T0=$(date +%s)
  set +e
  case "$CMD" in
    dryrun)
      ( cd "$RETO" && vcc prep "$H5" -g gene_names.csv --perts pert_counts.csv \
          -o "$VCC" --dry-run ) 2>&1 | tee "$OUT/entrega_prep_dryrun.log"
      EC=${PIPESTATUS[0]} ;;
    prep)
      ( cd "$RETO" && vcc prep "$H5" -g gene_names.csv --perts pert_counts.csv \
          -o "$VCC" -f ) 2>&1 | tee "$OUT/entrega_prep.log"
      EC=${PIPESTATUS[0]} ;;
    *) $CMD; EC=$? ;;
  esac
  set -e; set +u; set -u
  D=$(( $(date +%s) - T0 ))
  if [ "$EC" -eq 0 ]; then OK+=("$NOMBRE"); EST="OK"; else FALLO+=("$NOMBRE(exit=$EC)"); EST="FALLO exit=$EC <-- traza arriba"; fi
  set +e; "$PY" "$DIR/reporte2.py"; set -e
  echo "<<< $NOMBRE   FIN    : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
  echo "    duracion       : ${D} s  ($((D / 60)) min)"
  echo "    estado         : $EST"
  [ -f "$H5" ] && echo "    h5ad           : $(du -h "$H5" | cut -f1)"
  [ -f "$VCC" ] && echo "    vcc            : $(du -h "$VCC" | cut -f1)"
  echo "    /datos libre   : $(libre_gib) GiB"
  echo
done
DT=$(( $(date +%s) - T0T ))
echo "================================================================"
echo "RESUMEN"
echo "  fin           : $(date -u '+%Y-%m-%d %H:%M:%S') UTC"
echo "  duracion total: ${DT} s  ($((DT / 60)) min)"
echo "  OK            : ${OK[*]:-(ninguno)}"
echo "  FALLARON      : ${FALLO[*]:-(ninguno)}"
echo "  SIN EMPEZAR   : ${SIN[*]:-(ninguno)}"
[ -f "$H5" ] && echo "  h5ad: $H5  ($(stat -c%s "$H5") bytes)"
[ -f "$VCC" ] && echo "  🔴 ENTREGABLE LISTO (NO ENVIADO): $VCC  ($(stat -c%s "$VCC") bytes)"
echo "  /datos libre  : $(libre_gib) GiB"
echo "================================================================"
