#!/bin/bash
set -eo pipefail
# Print 1..PARAM_LINES on stdout, PARAM_LINES_PER_SECOND at a time with a
# one-second sleep after each burst. After each burst, "@ <n> <epoch ms>" says
# when line n was printed. Every PARAM_STDERR_EVERY lines, "stderr <k>" goes
# to stderr (0 = never).
LINES="${PARAM_LINES:-1000}"
PER_SECOND="${PARAM_LINES_PER_SECOND:-10}"
STDERR_EVERY="${PARAM_STDERR_EVERY:-0}"
count() {
  echo "count job started at $(date -u +%FT%T.%3NZ), ${LINES} lines at ${PER_SECOND}/s"
  for ((i = 1; i <= LINES; i++)); do
    echo "$i"
    if ((STDERR_EVERY > 0 && i % STDERR_EVERY == 0)); then
      echo "stderr $i" >&2
    fi
    if ((i % PER_SECOND == 0 || i == LINES)); then
      echo "@ $i $(date +%s%3N)"
      ((i == LINES)) || sleep 1
    fi
  done
  echo "count job finished at $(date -u +%FT%T.%3NZ)"
}
count | tee count.log
