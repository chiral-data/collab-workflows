#!/bin/bash
# The entry point, kept at the name workflow-037's suite uses so the habit carries over.
#
# The work is in run_fixtures.py. Python rather than shell because this harness stages a
# tree per case, parses each fixture's harness tokens, applies named mutations to JSON
# between steps and records every exit code -- all of which 037 did in bash and which is
# where most of its harness defects lived, including a fixture that asserted a substring
# no correct implementation could ever produce.
#
# Every argument is passed through: --out DIR, --only ID,..., --section S, --no-docker.
set -e
exec python3 "$(dirname "$0")/run_fixtures.py" "$@"
