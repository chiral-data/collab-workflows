#!/usr/bin/env python3
"""
Node 01: Validate target and exclusion FASTA files.

Checks format, non-empty sequences, minimum length, and nucleotide content.
Passes validated FASTAs and global_params.json through to outputs/.
"""

import json
import logging
import os
import shutil
import sys

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")
log = logging.getLogger("validate")

VALID_BASES = set("ACGTNacgtn")
MIN_SEQUENCE_LENGTH = 100  # bp — sequences shorter than this can't yield a 500bp ROI window


def _validate_fasta(path, label, require_nonempty=True):
    """
    Validate a FASTA file. Returns (records_list, error_str).
    error_str is None on success.
    """
    if not os.path.exists(path):
        return [], f"{label}: file not found at {path}"

    try:
        from Bio import SeqIO
        records = list(SeqIO.parse(path, "fasta"))
    except Exception as exc:
        return [], f"{label}: failed to parse FASTA — {exc}"

    # Filter placeholder/empty sequences (e.g. the exclusion placeholder)
    real_records = [r for r in records if len(r.seq) >= MIN_SEQUENCE_LENGTH]

    if require_nonempty and not real_records:
        # Soft warning for exclusion; hard error for target
        return real_records, None  # caller decides

    issues = []
    for rec in real_records:
        seq = str(rec.seq).upper()
        bad = set(seq) - VALID_BASES
        if bad:
            issues.append(f"  {rec.id}: unexpected characters {bad}")

    if issues:
        return real_records, f"{label} contains invalid nucleotide characters:\n" + "\n".join(issues)

    return real_records, None


def main():
    os.makedirs("./outputs", exist_ok=True)

    # Forward global_params.json
    for fname in ("global_params.json",):
        src = f"./inputs/{fname}"
        if os.path.exists(src):
            shutil.copy(src, f"./outputs/{fname}")

    errors = []
    summary = {}

    # Validate target.fasta
    tgt_path = "./inputs/target.fasta"
    tgt_recs, tgt_err = _validate_fasta(tgt_path, "target.fasta", require_nonempty=True)
    if tgt_err:
        errors.append(tgt_err)
    elif not tgt_recs:
        errors.append("target.fasta: no sequences >= 100 bp found")
    else:
        log.info("target.fasta: %d sequence(s), lengths %s",
                 len(tgt_recs),
                 [len(r.seq) for r in tgt_recs])
        summary["target"] = {
            "num_sequences": len(tgt_recs),
            "sequence_ids":  [r.id for r in tgt_recs],
            "lengths_bp":    [len(r.seq) for r in tgt_recs],
        }
        shutil.copy(tgt_path, "./outputs/target.fasta")

    # Validate exclusion.fasta (allowed to be empty/minimal for offline testing)
    excl_path = "./inputs/exclusion.fasta"
    excl_recs, excl_err = _validate_fasta(excl_path, "exclusion.fasta", require_nonempty=False)
    if excl_err:
        errors.append(excl_err)
    else:
        if not excl_recs:
            log.warning("exclusion.fasta: no real sequences found — ROI uniqueness will not be screened (all windows score 1.0)")
        else:
            log.info("exclusion.fasta: %d sequence(s)", len(excl_recs))
        summary["exclusion"] = {
            "num_sequences": len(excl_recs),
            "sequence_ids":  [r.id for r in excl_recs],
            "lengths_bp":    [len(r.seq) for r in excl_recs],
        }
        shutil.copy(excl_path, "./outputs/exclusion.fasta")

    if errors:
        for e in errors:
            log.error(e)
        sys.exit(1)

    summary["status"] = "ok"
    with open("./outputs/validation_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    log.info("Validation complete.")


if __name__ == "__main__":
    main()
