#!/bin/bash
# Does --msa-mode single_sequence throw away a precomputed a3m?
#
# Read out of colabfold.batch.run()'s source (1.6.1), lines 137-142:
#
#     if msa_mode == "single_sequence":
#         num_seqs = 1
#         if is_complex and "multimer" not in model_type: num_seqs += max_num
#         if use_templates: num_seqs += 4
#         max_seq = min(num_seqs, max_seq)
#
# and lines 258-270, which show that when a3m_lines is not None the MSA is read by
# unserialize_msa() and msa_mode is never consulted for MSA GENERATION at all.
#
# Put together, those say something alarming: an a3m handed to colabfold_batch together
# with --msa-mode single_sequence is parsed, is re-saved as <job>.a3m (line 287), and gets
# a coverage.png -- every downstream signal says "an alignment was used" -- while max_seq=1
# means the model is fed ONE sequence. §4's whole decision would be silently undone, and
# §11.1's spike concluded "coverage.png written, so an MSA was really consumed" from
# exactly those signals.
#
# So it is measured here instead of argued. Same a3m, same seed, same everything, two
# msa_mode values, offline. If the clamp is real the two runs disagree.
#
# Usage:  bash tests/probe_msa_clamp.sh [OUT_DIR]
set -uo pipefail

HERE="$(cd "$(dirname "$0")" && pwd)"
A3M="$HERE/reference/complex_target_only_msa.a3m"
OUT="${1:-$HOME/wf040-msaclamp/$(date +%Y%m%d-%H%M%S)}"
IMAGE="colabfold:2026_10_08"

if [ -z "${DOCKER_HOST:-}" ] && [ -S "/run/user/$(id -u)/docker.sock" ]; then
  export DOCKER_HOST="unix:///run/user/$(id -u)/docker.sock"
fi

[ -f "$A3M" ] || { echo "no $A3M" >&2; exit 2; }
free_gb() { df --output=avail -BG / | tail -1 | tr -dc '0-9'; }
[ "$(free_gb)" -lt 8 ] && { echo "only $(free_gb)G free, refusing" >&2; exit 1; }

mkdir -p "$OUT"
echo "a3m: $(grep -c '^>' "$A3M") sequence(s), header $(head -1 "$A3M")"
echo "out: $OUT"

for mode in single_sequence mmseqs2_uniref_env; do
  echo
  echo "=== --msa-mode $mode"
  mkdir -p "$OUT/$mode/in" "$OUT/$mode/out"
  cp "$A3M" "$OUT/$mode/in/probe.a3m"
  # --network none throughout. With an a3m input colabfold never calls
  # get_msa_and_templates, so use_msa/use_env are only used to write the bibtex and the
  # config -- no egress either way. The flag proves it rather than asserting it.
  docker run --rm --gpus all --network none \
    -v "$OUT:$OUT" -w "$OUT/$mode" "$IMAGE" \
    colabfold_batch "$OUT/$mode/in/probe.a3m" "$OUT/$mode/out" \
      --model-type alphafold2_multimer_v3 \
      --msa-mode "$mode" \
      --num-models 1 \
      --num-recycle 3 \
    > "$OUT/$mode.log" 2>&1
  echo "  rc=$? $(grep -ciE 'error|traceback' "$OUT/$mode.log") error line(s)"
  grep -iE 'max_seq|num_seqs|sequences|running on|Query ' "$OUT/$mode.log" | head -8 \
    | sed 's/^/    /'
done

echo
echo "=== what the model actually saw, and what came out"
python3 -I "$HERE/report_msa_clamp.py" "$OUT" | tee "$OUT/VERDICT.txt"
