#!/usr/bin/env bash
# Sube ~/vcc/out a S3 cada 120 s.
while true; do
  echo "[$(date -u +%FT%TZ)] sync"
  aws s3 sync "$HOME/vcc/out" s3://vc-challenge-data/out --only-show-errors
  echo "[$(date -u +%FT%TZ)] exit=$?"
  sleep 120
done
