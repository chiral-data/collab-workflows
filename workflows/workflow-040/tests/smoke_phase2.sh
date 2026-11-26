#!/bin/bash
# A hand smoke test of the Phase 2 path, on the host interpreter, with no GPU.
#
# Not a replacement for the suite -- it decides nothing and asserts nothing. It exists to
# walk the complex path once end to end while it is being built, printing what each node
# wrote, so that the fixtures written against it are written against something that has
# actually run. The suite's own Phase 2 cases are the real check.
#
# Usage:  bash tests/smoke_phase2.sh [OUT_DIR]
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
WF="$(dirname "$HERE")"
WORKFLOWS="$(dirname "$WF")"
OUT="${1:-${TMPDIR:-/tmp}/wf040-smoke2/$(date +%Y%m%d-%H%M%S)}"

say() { printf '\n=== %s\n' "$*"; }
run_node() {        # run_node <step> <script> [extra env assignments...]
  local step="$1" script="$2"; shift 2
  mkdir -p "$OUT/$step/inputs" "$OUT/$step/outputs"
  ( cd "$OUT/$step" && env "$@" python3 "$script" ) 2>&1 | tee "$OUT/$step.log"
  return "${PIPESTATUS[0]}"
}

mkdir -p "$OUT/g"
say "generate: 2 designs with a target, via 037's generator"
python3 -I "$WORKFLOWS/workflow-037/tests/make_designs.py" "$OUT/g" 2 compact 30 \
  --manifest > "$OUT/g.log" 2>&1 || { cat "$OUT/g.log"; exit 1; }
ls "$OUT/g" | sed 's/^/  /'

say "037 node 04: the real design_metrics.json"
mkdir -p "$OUT/m/inputs"
cp "$OUT"/g/*_generated.cif "$OUT/g/generation_manifest.json" "$OUT/m/inputs/"
run_node m "$WORKFLOWS/workflow-037/04-analyze/script.py" >/dev/null || true

say "040 node 01: ingest"
mkdir -p "$OUT/n1/inputs"
cp "$OUT"/g/*_generated.cif "$OUT/m/outputs/design_metrics.json" "$OUT/n1/inputs/"
run_node n1 "$WF/01-ingest-designs/script.py" || exit 1

# A stand-in target alignment. The target sequence has to come out of the inventory, so
# the a3m is built from it here rather than hard-coded -- a hard-coded one would make the
# query mismatch and 01b would correctly refuse it, which would test the refusal instead
# of the path.
say "build a stand-in target alignment from the inventory"
mkdir -p "$OUT/n1b/inputs"
cp "$OUT/n1/outputs/design_inventory.json" "$OUT/n1b/inputs/"
python3 -I - "$OUT/n1b/inputs/design_inventory.json" "$OUT/n1b/inputs/target_msa.a3m" <<'PY'
import json
import sys

inventory = json.load(open(sys.argv[1]))
live = [r for r in inventory["designs"] if not r.get("skipped")]
target = "".join(c["sequence"] for c in live[0]["target_chains"])

# Four synthetic homologs: one point substitution each, so every row lines up with the
# query in exactly the way match_columns() requires. Synthetic on purpose -- nothing is
# sent anywhere to build a smoke-test fixture.
rows = [("query", target)]
for i in range(4):
    at = (i * 7) % len(target)
    swap = "A" if target[at] != "A" else "G"
    rows.append((f"hom{i + 1}", target[:at] + swap + target[at + 1:]))

with open(sys.argv[2], "w") as handle:
    for name, seq in rows:
        handle.write(f">{name}\n{seq}\n")
print(f"  target {len(target)} residues, {len(rows)} alignment row(s)")
PY

say "040 node 01b: align target (source=file, no network)"
run_node n1b "$WF/01b-align-target/script.py" \
  PARAM_SCORE_INTERFACE=true PARAM_MSA_SOURCE=file || exit 1
python3 -c "
import json,sys
d=json.load(open('$OUT/n1b/outputs/target_alignment.json'))
for k in ('aligned','source','n_sequences','n_homologs','target_n_residues','binder_aligned','network_used'):
    print(f'  {k}: {d.get(k)}')
"

say "040 node 02: build BOTH inputs"
mkdir -p "$OUT/n2/inputs"
cp "$OUT"/n1/outputs/* "$OUT"/n1b/outputs/* "$OUT/n2/inputs/"
run_node n2 "$WF/02-build-prediction-inputs/script.py" PARAM_SCORE_INTERFACE=true || exit 1
ls "$OUT/n2/outputs" | sed 's/^/  /'
say "the paired a3m, as written"
head -8 "$OUT/n2/outputs"/pairedmsa_design_0.a3m | cut -c1-100 | sed 's/^/  /'

say "040 node 03: DRY RUN, so the argv for both paths is visible without a GPU"
mkdir -p "$OUT/n3/inputs"
cp "$OUT"/n2/outputs/* "$OUT/n3/inputs/"
( cd "$OUT/n3" && env PARAM_SCORE_INTERFACE=true python3 "$WF/03-predict/script.py" \
  --dry-run ) 2>&1 | tee "$OUT/n3.log"
python3 -c "
import json
d=json.load(open('$OUT/n3/outputs/prediction_manifest.json'))
print(f\"  phase={d['phase']}  n_declared={d['n_declared']}\")
print(f\"  mono: model_type={d['model_type']} msa_mode={d['msa_mode']} recycles={d['num_recycle']}\")
c=d.get('cplx')
print(f\"  cplx: {c['model_type']} msa_mode={c['msa_mode']} recycles={c['num_recycle']} binder_aligned={c['binder_aligned']}\" if c else '  cplx: none')
for p in d['predictions']:
    print(f\"   {p['design']:10} {p['mode']:5} chains={p['chain_order']} lengths={p['chain_lengths']} msa={p.get('msa_source','-')[:40]}\")
"

say "done. tree under $OUT"
