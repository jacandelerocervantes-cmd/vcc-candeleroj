#!/usr/bin/env bash
# REEMISION A k=3.4. Desprendido. ARCHIVO NUEVO: prediction_k3.4.vcc
# 🔴 NO ENVIA NADA: no hay `vcc submit` ni `vcc sample` en ningun paso.
# 🔴 El respaldo /datos/modelo/out/entrega/prediction.vcc NO se toca: se le
#    verifica el sha256 al empezar y al terminar.
# Salida incremental por paso; un paso caido no mata a los demas; aborta bajo 10 GiB.
set -u
PY=/home/ubuntu/mf/envs/vcc/bin/python
SRC=/home/ubuntu/vcc/modelo/src/entrega
OUT=/home/ubuntu/vcc/modelo/out
RETO=/datos/reto/final
W=/datos/modelo/out/entrega
H5=$W/prediction_k3.4.h5ad
VCC=$W/prediction_k3.4.vcc
RESP=$W/prediction.vcc
SHA_RESP=66149bf6794c9f5d8cc062599eff9856451a62456e5a3f5ad710749656c5bff1
exec 9>"$W/.lock_k"; flock -n 9 || { echo "ya corre otra reemision"; exit 1; }
libre(){ df -P /datos | awk 'NR==2{printf "%.0f", $4/1048576}'; }
OK=(); FALLO=()
paso(){
  NOMBRE="$1"; shift
  echo ""; echo "================================================================"
  echo ">>> $NOMBRE   INICIO : $(date -u '+%F %T UTC')  (/datos $(libre) GiB)"
  T0=$(date +%s)
  if [ "$(libre)" -lt 10 ]; then echo "ABORTA: /datos bajo 10 GiB"; FALLO+=("$NOMBRE(disco)"); return; fi
  set +e; "$@"; EC=$?; set -e
  D=$(( $(date +%s) - T0 ))
  if [ "$EC" -eq 0 ]; then OK+=("$NOMBRE"); EST="OK"; else FALLO+=("$NOMBRE(exit=$EC)"); EST="FALLO exit=$EC"; fi
  echo "<<< $NOMBRE   FIN    : $(date -u '+%F %T UTC')"
  echo "    duracion    : ${D} s"
  echo "    estado      : $EST"
  [ -f "$H5" ] && echo "    h5ad        : $(du -h "$H5" | cut -f1)"
  [ -f "$VCC" ] && echo "    vcc         : $(du -h "$VCC" | cut -f1)"
  echo "    /datos libre: $(libre) GiB"
}
echo "ARRANQUE $(date -u '+%F %T UTC')   k=3.4   /datos $(libre) GiB"
echo "respaldo al empezar: $(sha256sum "$RESP" | cut -d' ' -f1)"
T00=$(date +%s)

paso "1 EMISION k=3.4 (D/E/F -> h5ad NUEVO)" "$PY" -u "$SRC/emite_k.py"
paso "2 CONTRATO (las 17, una por una)"      "$PY" -u "$SRC/contrato_k.py"
paso "3 PREP --dry-run"  bash -c "cd $RETO && vcc prep $H5 -g gene_names.csv --perts pert_counts.csv -o $VCC --dry-run 2>&1 | tee $OUT/entrega_k_prep_dryrun.log"
paso "4 PREP (escribe el .vcc NUEVO)" bash -c "cd $RETO && vcc prep $H5 -g gene_names.csv --perts pert_counts.csv -o $VCC -f 2>&1 | tee $OUT/entrega_k_prep.log"

echo ""; echo "================================================================"
echo "RESUMEN   fin $(date -u '+%F %T UTC')   total $(( $(date +%s) - T00 )) s"
echo "  OK      : ${OK[*]:-(ninguno)}"
echo "  FALLARON: ${FALLO[*]:-(ninguno)}"
if [ -f "$VCC" ]; then
  echo "  🔴 NUEVO LISTO (NO ENVIADO): $VCC"
  echo "     bytes : $(stat -c %s "$VCC")"
  echo "     sha256: $(sha256sum "$VCC" | cut -d' ' -f1)"
fi
S2=$(sha256sum "$RESP" | cut -d' ' -f1)
echo "  respaldo: $S2"
[ "$S2" = "$SHA_RESP" ] && echo "     ✓ intacto, coincide con 66149bf6...5bff1" || echo "     🔴 EL RESPALDO CAMBIO"
echo "  /datos libre: $(libre) GiB"
echo "================================================================"
