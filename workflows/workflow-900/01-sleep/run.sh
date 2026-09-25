#!/bin/bash
set -e
# Test-only step: stay in Processing long enough for the cancellation E2E
# (chiral-test-e2e#74) to cancel it. A job that runs to completion means the
# cancel never reached the backend.
echo "sleep job started at $(date -u +%FT%TZ)" | tee sleep.log
for i in $(seq 1 60); do
  sleep 10
  echo "still sleeping: $((i * 10))s" | tee -a sleep.log
done
echo "sleep job finished at $(date -u +%FT%TZ)" | tee -a sleep.log
