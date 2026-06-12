#!/bin/bash
set -e
echo "Starting Node 04: BLAST Validation and Report"
# Exit code 2 means no primer sets passed all constraints — expected outcome, not an error
set +e
python blast_report.py
exit_code=$?
set -e
if [ $exit_code -eq 2 ]; then
  exit 0
fi
[ $exit_code -eq 0 ] || exit $exit_code
echo "Node 04 completed"
