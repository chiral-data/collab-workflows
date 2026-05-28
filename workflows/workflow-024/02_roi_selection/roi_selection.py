#!/usr/bin/env python3
"""
Node 02: ROI Selection via k-mer subtraction.

Slides a window across the target sequence(s) and scores each window by
the fraction of k-mers absent from the exclusion (close-relative) set.
Top-scoring non-overlapping windows are saved as rois.json.
"""

import json
import logging
import os
import shutil
import sys

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")
log = logging.getLogger("roi_selection")

KMER_K = 20


def _load_params():
    gp = {}
    gp_path = "./inputs/global_params.json"
    if os.path.exists(gp_path):
        with open(gp_path) as f:
            gp = json.load(f)
    return {
        "roi_window":     int(os.environ.get("PARAM_ROI_WINDOW",    str(gp.get("roi_window", 500)))),
        "roi_step":       int(os.environ.get("PARAM_ROI_STEP",      str(gp.get("roi_step", 100)))),
        "min_uniqueness": float(os.environ.get("PARAM_MIN_UNIQUENESS", str(gp.get("min_uniqueness", 0.80)))),
        "top_rois":       int(os.environ.get("PARAM_TOP_ROIS",      str(gp.get("top_rois", 3)))),
    }


def _build_kmer_set(records, k):
    kmer_set = set()
    for idx, rec in enumerate(records, 1):
        log.info("  Building k-mer set: seq %d/%d (%s, len=%d)", idx, len(records), rec.id, len(rec.seq))
        seq = str(rec.seq).upper().replace("N", "")
        rc  = str(rec.seq.reverse_complement()).upper().replace("N", "")
        for s in (seq, rc):
            for i in range(len(s) - k + 1):
                kmer_set.add(s[i:i + k])
    log.info("Exclusion k-mer set: %d k-mers (k=%d)", len(kmer_set), k)
    return kmer_set


def _score_window(seq, exclusion_kmers, k):
    total = unique = 0
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        total += 1
        if kmer not in exclusion_kmers:
            unique += 1
    return unique / total if total > 0 else 0.0


def _deduplicate(candidates, min_gap):
    kept = []
    for cand in candidates:
        overlap = False
        for kc in kept:
            if kc["record_id"] != cand["record_id"]:
                continue
            if min(cand["end"], kc["end"]) - max(cand["start"], kc["start"]) > min_gap:
                overlap = True
                break
        if not overlap:
            kept.append(cand)
    return kept


def _find_rois(target_records, exclusion_records, window_size, step, min_uniqueness, top_n):
    if not exclusion_records:
        log.warning("No exclusion sequences — all windows will score 1.0 (no exclusion screening)")
        exclusion_kmers = set()
    else:
        exclusion_kmers = _build_kmer_set(exclusion_records, KMER_K)

    candidates = []

    for rec in target_records:
        seq    = str(rec.seq).upper()
        seqlen = len(seq)

        # Auto-reduce window for short sequences (gene-targeted mode)
        win = window_size
        if seqlen < win:
            win = max(100, seqlen - 20)
            log.info("Short sequence (%dbp) — reducing window to %d", seqlen, win)

        log.info("Scanning %s (len=%d, window=%d, step=%d)", rec.id, seqlen, win, step)
        n_windows = (seqlen - win) // step + 1

        for i, start in enumerate(range(0, seqlen - win + 1, step)):
            if i > 0 and i % 10000 == 0:
                log.info("  %.0f%% (%d/%d windows)", i / n_windows * 100, i, n_windows)

            window = seq[start:start + win]
            if window.count("N") / len(window) > 0.05:
                continue

            score = _score_window(window, exclusion_kmers, KMER_K)
            if score >= min_uniqueness:
                candidates.append({
                    "record_id":        rec.id,
                    "start":            start,
                    "end":              start + win,
                    "length":           win,
                    "sequence":         window,
                    "uniqueness_score": round(score, 4),
                })

        log.info("  Scan complete for %s", rec.id)

    candidates.sort(key=lambda c: c["uniqueness_score"], reverse=True)

    # Cap before dedup to avoid O(n²) on large genomes
    pre_cap = top_n * 20
    if len(candidates) > pre_cap:
        candidates = candidates[:pre_cap]

    candidates = _deduplicate(candidates, min_gap=window_size // 2)
    return candidates[:top_n]


def main():
    params = _load_params()
    os.makedirs("./outputs", exist_ok=True)

    for fname in ("global_params.json", "target.fasta", "exclusion.fasta"):
        src = f"./inputs/{fname}"
        if os.path.exists(src):
            shutil.copy(src, f"./outputs/{fname}")

    from Bio import SeqIO

    target_records = list(SeqIO.parse("./inputs/target.fasta", "fasta"))
    if not target_records:
        log.error("target.fasta is empty")
        sys.exit(1)

    excl_path = "./inputs/exclusion.fasta"
    excl_records = [r for r in SeqIO.parse(excl_path, "fasta") if len(r.seq) >= 100] if os.path.exists(excl_path) else []

    log.info("Target: %d sequence(s), Exclusion: %d sequence(s)", len(target_records), len(excl_records))

    window  = params["roi_window"]
    step    = params["roi_step"]
    min_u   = params["min_uniqueness"]
    top_n   = params["top_rois"]

    rois = _find_rois(target_records, excl_records, window, step, min_u, top_n)

    if not rois:
        log.warning("No ROIs at uniqueness=%.2f — retrying at 0.60", min_u)
        rois = _find_rois(target_records, excl_records, window, step, 0.60, top_n)

    if not rois:
        log.error("No ROI candidates found. Try --min-uniqueness 0.5 or check input sequences.")
        sys.exit(1)

    log.info("Selected %d ROI candidate(s):", len(rois))
    for roi in rois:
        log.info("  [%s:%d-%d] uniqueness=%.3f", roi["record_id"], roi["start"], roi["end"], roi["uniqueness_score"])

    with open("./outputs/rois.json", "w") as f:
        json.dump({"n_rois": len(rois), "rois": rois}, f, indent=2)

    log.info("ROI selection complete.")


if __name__ == "__main__":
    main()
