#!/usr/bin/env bash
# BARRIDO DE k sobre VCC 2025, blancos retenidos. Desprendido.
# 🔴 NO ENVIA NADA: no hay `vcc submit` ni `vcc sample` en ningun paso.
# 🔴 No toca /datos/modelo/out/entrega/prediction.vcc (el respaldo).
# Salida incremental por brazo; un brazo caido no mata a los demas (try/except
# dentro de bk.py, que guarda el json tras cada uno); aborta bajo 10 GiB.
set -u
PY=/home/ubuntu/mf/envs/vcc/bin/python
OUT=/home/ubuntu/vcc/modelo/out
WORK=/datos/modelo/out/barrido
mkdir -p "$WORK"
exec 9>"$WORK/.lock"; flock -n 9 || { echo "ya corre otro barrido"; exit 1; }
libre(){ df -P /datos | awk 'NR==2{printf "%.0f", $4/1048576}'; }
echo "================================================================"
echo ">>> BARRIDO k   INICIO : $(date -u '+%F %T UTC')  (/datos $(libre) GiB)"
echo "    respaldo intacto: $(stat -c %s /datos/modelo/out/entrega/prediction.vcc) bytes"
T0=$(date +%s)
if [ "$(libre)" -lt 10 ]; then echo "ABORTA: /datos bajo 10 GiB"; exit 2; fi
set +e
"$PY" -u /home/ubuntu/vcc/modelo/src/barrido/bk.py
EC=$?
set -e
echo "<<< BARRIDO k   FIN    : $(date -u '+%F %T UTC')   exit=$EC   $(( $(date +%s) - T0 ))s"
echo "    /datos libre: $(libre) GiB"
echo "    json        : $OUT/barrido_k.json"
"$PY" - <<'EOF'
import json
r = json.load(open("/home/ubuntu/vcc/modelo/out/barrido_k.json"))
b = r.get("brazos", {})
print("    brazos OK   :", [k for k, v in b.items() if v.get("estado") == "OK"] or "(ninguno)")
print("    brazos otros:", {k: v.get("estado") for k, v in b.items() if v.get("estado") != "OK"} or "(ninguno)")
EOF
echo "    respaldo intacto: $(stat -c %s /datos/modelo/out/entrega/prediction.vcc) bytes"
echo "================================================================"
