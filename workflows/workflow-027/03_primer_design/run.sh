#!/bin/bash
set -e
echo "Starting Node 03: Primer Design"
# Exit code 2 means no primer sets could be designed from any ROI — expected outcome, not an error
set +e
python primer_design.py
exit_code=$?
set -e
if [ $exit_code -eq 2 ]; then
  exit 0
fi
[ $exit_code -eq 0 ] || exit $exit_code
echo "Node 03 completed"
