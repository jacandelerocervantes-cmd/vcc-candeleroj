#!/usr/bin/env bash
# Vigila el aviso de interrupción spot (IMDSv2). 404 = sin aviso.
# 200 = aviso: sync final de ~/vcc/out y salir.
# 401 = token vencido: se renueva. Otros códigos o fallo de red: se registra y se reintenta,
# para que un fallo transitorio del IMDS no dispare la salida.
IMDS=http://169.254.169.254
token() { curl -s -X PUT "$IMDS/latest/api/token" -H "X-aws-ec2-metadata-token-ttl-seconds: 21600"; }
TOK=$(token); T0=$(date +%s)
echo "[$(date -u +%FT%TZ)] inicio, token len=${#TOK}"
while true; do
  (( $(date +%s) - T0 > 18000 )) && { TOK=$(token); T0=$(date +%s); }
  CODE=$(curl -s -o /tmp/spot_action.json -w '%{http_code}' -H "X-aws-ec2-metadata-token: $TOK" \
         "$IMDS/latest/meta-data/spot/instance-action")
  case "$CODE" in
    404) ;;
    200)
      echo "[$(date -u +%FT%TZ)] AVISO DE INTERRUPCIÓN: $(cat /tmp/spot_action.json)"
      aws s3 sync "$HOME/vcc/out" s3://vc-challenge-data/out --only-show-errors
      echo "[$(date -u +%FT%TZ)] sync final exit=$?; saliendo"
      exit 0 ;;
    401) echo "[$(date -u +%FT%TZ)] 401, renovando token"; TOK=$(token); T0=$(date +%s) ;;
    *)   echo "[$(date -u +%FT%TZ)] código inesperado '$CODE', reintento" ;;
  esac
  sleep 5
done
