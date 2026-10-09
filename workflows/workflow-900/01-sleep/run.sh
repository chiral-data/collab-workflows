#!/bin/bash
set -e
# Test-only step: stay in Processing long enough for the cancellation E2E
# (chiral-test-e2e#74) to cancel it. A job that runs to completion means the
# cancel never reached the backend.
SLEEP_SECONDS="${PARAM_SLEEP_SECONDS:-120}"
echo "sleep job started at $(date -u +%FT%TZ), sleeping ${SLEEP_SECONDS}s" | tee sleep.log
for i in $(seq "$SLEEP_SECONDS"); do
  echo "$i" | tee -a sleep.log
  sleep 1
done
echo "sleep job finished at $(date -u +%FT%TZ)" | tee -a sleep.log
