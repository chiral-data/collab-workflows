#!/bin/bash
set -e
# Stay in Processing for PARAM_SLEEP_SECONDS seconds (default 120).
SLEEP_SECONDS="${PARAM_SLEEP_SECONDS:-120}"
echo "sleep job started at $(date -u +%FT%TZ), sleeping ${SLEEP_SECONDS}s" | tee sleep.log
for i in $(seq "$SLEEP_SECONDS"); do
  echo "$i" | tee -a sleep.log
  sleep 1
done
echo "sleep job finished at $(date -u +%FT%TZ)" | tee -a sleep.log
