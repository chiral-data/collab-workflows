#!/usr/bin/env python3
"""
Node 03: Primer/Probe Design via Primer3.

Reads rois.json, runs Primer3 on each ROI to design forward primer,
reverse primer, and hydrolysis probe. Validates thermodynamic constraints
(Tm, GC, probe delta-Tm, self-dimer, hairpin). Outputs primer_sets.json
ranked by constraint pass status and Primer3 penalty score.
"""

import json
import logging
import os
import shutil
import sys

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")
log = logging.getLogger("primer_design")

# Default constraints (overridden by params)
DEFAULTS = {
    "primer_min_size": 19, "primer_opt_size": 20, "primer_max_size": 26,
    "primer_min_tm": 59.0, "primer_opt_tm": 60.0, "primer_max_tm": 62.0,
    "primer_min_gc": 40.0, "primer_max_gc": 60.0,
    "probe_min_size": 22, "probe_opt_size": 25, "probe_max_size": 29,
    "probe_min_gc": 40.0, "probe_max_gc": 65.0,
    "probe_tm_delta_min": 5.0, "probe_tm_delta_max": 10.0,
    "amplicon_min": 70, "amplicon_max": 200,
}
DIMER_THRESHOLD = -6.0  # kcal/mol


def _load_params():
    gp = {}
    if os.path.exists("./inputs/global_params.json"):
        with open("./inputs/global_params.json") as f:
            gp = json.load(f)
    p = dict(DEFAULTS)
    for key in DEFAULTS:
        env_key = f"PARAM_{key.upper()}"
        if env_key in os.environ:
            val = os.environ[env_key]
            p[key] = float(val) if "." in val or key in ("primer_min_tm", "primer_opt_tm", "primer_max_tm") else int(val)
        elif key in gp:
            p[key] = gp[key]
    return p


def _build_primer3_globals(p):
    opt_tm = p["primer_opt_tm"]
    delta_min, delta_max = p["probe_tm_delta_min"], p["probe_tm_delta_max"]
    return {
        "PRIMER_MIN_SIZE":                      p["primer_min_size"],
        "PRIMER_OPT_SIZE":                      p["primer_opt_size"],
        "PRIMER_MAX_SIZE":                      p["primer_max_size"],
        "PRIMER_MIN_TM":                        p["primer_min_tm"],
        "PRIMER_OPT_TM":                        opt_tm,
        "PRIMER_MAX_TM":                        p["primer_max_tm"],
        "PRIMER_MIN_GC":                        p["primer_min_gc"],
        "PRIMER_OPT_GC":                        50.0,
        "PRIMER_MAX_GC":                        p["primer_max_gc"],
        "PRIMER_PRODUCT_SIZE_RANGE":            [[p["amplicon_min"], p["amplicon_max"]]],
        "PRIMER_INTERNAL_MIN_SIZE":             p["probe_min_size"],
        "PRIMER_INTERNAL_OPT_SIZE":             p["probe_opt_size"],
        "PRIMER_INTERNAL_MAX_SIZE":             p["probe_max_size"],
        "PRIMER_INTERNAL_MIN_GC":               p["probe_min_gc"],
        "PRIMER_INTERNAL_MAX_GC":               p["probe_max_gc"],
        "PRIMER_INTERNAL_MIN_TM":               opt_tm + delta_min,
        "PRIMER_INTERNAL_OPT_TM":               opt_tm + (delta_min + delta_max) / 2,
        "PRIMER_INTERNAL_MAX_TM":               p["primer_max_tm"] + delta_max,
        "PRIMER_PICK_INTERNAL_OLIGO":           1,
        "PRIMER_NUM_RETURN":                    5,
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "PRIMER_MAX_SELF_ANY":                  8,
        "PRIMER_MAX_SELF_END":                  3,
        "PRIMER_PAIR_MAX_COMPL_ANY":            8,
        "PRIMER_PAIR_MAX_COMPL_END":            3,
        "PRIMER_MAX_POLY_X":                    4,
        "PRIMER_INTERNAL_MAX_SELF_ANY":         8,
    }


def _self_dimer_dg(seq):
    try:
        import primer3
        return primer3.calc_homodimer(seq).dg / 1000.0
    except Exception:
        return 0.0


def _hairpin_dg(seq):
    try:
        import primer3
        return primer3.calc_hairpin(seq).dg / 1000.0
    except Exception:
        return 0.0


def _validate_set(ps, p):
    notes = []
    ok = True

    for name, tm in [("FWD", ps["fwd_tm"]), ("REV", ps["rev_tm"])]:
        if not (p["primer_min_tm"] <= tm <= p["primer_max_tm"]):
            notes.append(f"{name} Tm {tm:.1f}°C outside [{p['primer_min_tm']},{p['primer_max_tm']}]")
            ok = False

    mean_tm = (ps["fwd_tm"] + ps["rev_tm"]) / 2
    delta = ps["probe_tm"] - mean_tm
    if not (p["probe_tm_delta_min"] <= delta <= p["probe_tm_delta_max"]):
        notes.append(f"Probe Tm delta {delta:.1f}°C outside [{p['probe_tm_delta_min']},{p['probe_tm_delta_max']}]")
        ok = False

    if not (p["amplicon_min"] <= ps["amplicon_size"] <= p["amplicon_max"]):
        notes.append(f"Amplicon {ps['amplicon_size']}nt outside [{p['amplicon_min']},{p['amplicon_max']}]")
        ok = False

    for name, seq in [("FWD", ps["fwd_seq"]), ("REV", ps["rev_seq"]), ("PRB", ps["probe_seq"])]:
        if not seq:
            continue
        dg = _self_dimer_dg(seq)
        if dg < DIMER_THRESHOLD:
            notes.append(f"{name} self-dimer ΔG={dg:.1f} kcal/mol")
            ok = False
        hp = _hairpin_dg(seq)
        if hp < DIMER_THRESHOLD:
            notes.append(f"{name} hairpin ΔG={hp:.1f} kcal/mol")
            ok = False

    ps["passed_constraints"] = ok
    ps["constraint_notes"]   = notes
    return ps


def _design_for_roi(roi, global_args, p):
    import primer3
    seq_args = {
        "SEQUENCE_ID":       f"{roi['record_id']}_{roi['start']}_{roi['end']}",
        "SEQUENCE_TEMPLATE": roi["sequence"].upper(),
        "SEQUENCE_PRIMER_PAIR_OK_REGION_LIST": [],
    }
    try:
        result = primer3.bindings.design_primers(seq_args, global_args)
    except Exception as exc:
        log.error("Primer3 failed for ROI %s:%d-%d: %s", roi["record_id"], roi["start"], roi["end"], exc)
        return []

    n = result.get("PRIMER_PAIR_NUM_RETURNED", 0)
    log.info("  ROI %s:%d-%d → %d pair(s)", roi["record_id"], roi["start"], roi["end"], n)

    sets = []
    for i in range(n):
        try:
            ps = {
                "roi":           roi,
                "pair_index":    i,
                "fwd_seq":       result[f"PRIMER_LEFT_{i}_SEQUENCE"],
                "rev_seq":       result[f"PRIMER_RIGHT_{i}_SEQUENCE"],
                "probe_seq":     result.get(f"PRIMER_INTERNAL_{i}_SEQUENCE", ""),
                "fwd_tm":        result[f"PRIMER_LEFT_{i}_TM"],
                "rev_tm":        result[f"PRIMER_RIGHT_{i}_TM"],
                "probe_tm":      result.get(f"PRIMER_INTERNAL_{i}_TM", 0.0),
                "fwd_gc":        result[f"PRIMER_LEFT_{i}_GC_PERCENT"],
                "rev_gc":        result[f"PRIMER_RIGHT_{i}_GC_PERCENT"],
                "probe_gc":      result.get(f"PRIMER_INTERNAL_{i}_GC_PERCENT", 0.0),
                "fwd_start":     result[f"PRIMER_LEFT_{i}"][0],
                "rev_start":     result[f"PRIMER_RIGHT_{i}"][0],
                "probe_start":   result.get(f"PRIMER_INTERNAL_{i}", [0])[0],
                "amplicon_size": result[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                "penalty":       result[f"PRIMER_PAIR_{i}_PENALTY"],
                "blast_fwd":     None, "blast_rev": None, "blast_probe": None, "blast_pass": None,
            }
            ps = _validate_set(ps, p)
            sets.append(ps)
        except KeyError as exc:
            log.debug("Missing Primer3 key for pair %d: %s", i, exc)

    return sets


def main():
    p = _load_params()
    os.makedirs("./outputs", exist_ok=True)

    for fname in ("global_params.json", "target.fasta", "exclusion.fasta"):
        src = f"./inputs/{fname}"
        if os.path.exists(src):
            shutil.copy(src, f"./outputs/{fname}")

    with open("./inputs/rois.json") as f:
        roi_data = json.load(f)

    rois = roi_data.get("rois", [])
    if not rois:
        log.warning("rois.json contains no ROI candidates — no primers to design")
        with open("./outputs/primer_sets.json", "w") as f:
            json.dump({"n_sets": 0, "primer_sets": []}, f, indent=2)
        sys.exit(2)

    log.info("Designing primers for %d ROI(s)", len(rois))
    global_args = _build_primer3_globals(p)
    all_sets = []

    for idx, roi in enumerate(rois, 1):
        log.info("ROI %d/%d: %s:%d-%d (len=%d, uniqueness=%.3f)",
                 idx, len(rois), roi["record_id"], roi["start"], roi["end"],
                 roi["length"], roi["uniqueness_score"])
        sets = _design_for_roi(roi, global_args, p)
        all_sets.extend(sets)

    if not all_sets:
        log.warning("Primer3 returned no primer sets for any ROI — try relaxing constraints")
        with open("./outputs/primer_sets.json", "w") as f:
            json.dump({"n_sets": 0, "primer_sets": []}, f, indent=2)
        sys.exit(2)

    # Sort: passing sets first, then by penalty
    all_sets.sort(key=lambda ps: (not ps["passed_constraints"], ps["penalty"]))

    n_pass = sum(1 for ps in all_sets if ps["passed_constraints"])
    log.info("Total primer sets: %d  (%d pass all constraints)", len(all_sets), n_pass)

    with open("./outputs/primer_sets.json", "w") as f:
        json.dump({"n_sets": len(all_sets), "primer_sets": all_sets}, f, indent=2)

    log.info("Primer design complete.")


if __name__ == "__main__":
    main()
