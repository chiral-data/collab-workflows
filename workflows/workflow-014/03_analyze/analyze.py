#!/usr/bin/env python3
"""Compute MPO scores, apply configurable ADMET filters, rank top N candidates."""

import json
import os
import sys

import numpy as np
import pandas as pd

# =============================================================================
# CONFIGURATION — params come from PARAM_* environment variables (Silva)
# =============================================================================

INPUT_PATH = "./inputs/raw_predictions.csv"
OUTPUT_CSV = "./outputs/filtered_candidates.csv"
OUTPUT_JSON = "./outputs/analysis_summary.json"

TOP_N = int(os.environ.get("PARAM_TOP_N", "20"))

# Filter thresholds — parsed from PARAM_FILTER_* env vars
# Supported formats: "safe", "strict", ">0.5", "<0.3", ">=0.5", "<=0.3"
FILTER_PARAMS = {
    "hERG": os.environ.get("PARAM_FILTER_HERG", "safe"),
    "Caco2_Wang": os.environ.get("PARAM_FILTER_CACO2", ">-5.15"),
    "BBB_Martins": os.environ.get("PARAM_FILTER_BBB", ">0.5"),
    "HIA_Hou": os.environ.get("PARAM_FILTER_HIA", ">0.5"),
    "DILI": os.environ.get("PARAM_FILTER_DILI", "<0.5"),
    "AMES": os.environ.get("PARAM_FILTER_AMES", "<0.5"),
    "ClinTox": os.environ.get("PARAM_FILTER_CLINTOX", "<0.3"),
}

# Named presets for convenience
NAMED_PRESETS = {
    "hERG": {"safe": "<0.5", "strict": "<0.3"},
}

# MPO scoring weights — positive means higher is better, negative means lower is better
MPO_WEIGHTS = {
    "BBB_Martins": 1.0,
    "Bioavailability_Ma": 1.0,
    "HIA_Hou": 1.0,
    "hERG": -1.0,
    "DILI": -1.0,
    "AMES": -1.0,
    "ClinTox": -1.0,
    "Lipinski": 0.5,
    "QED": 1.0,
    "Caco2_Wang": 0.5,
}


def parse_filter(raw_value):
    """Parse a filter string into (operator, threshold).

    Supports: ">0.5", "<0.3", ">=0.5", "<=0.3", "safe", "strict", or plain number.
    Returns (op_func, threshold) or None if the filter should be skipped.
    """
    val = raw_value.strip()
    if not val or val.lower() == "none":
        return None

    if val.startswith(">="):
        return (lambda x, t: x >= t, float(val[2:]))
    elif val.startswith("<="):
        return (lambda x, t: x <= t, float(val[2:]))
    elif val.startswith(">"):
        return (lambda x, t: x > t, float(val[1:]))
    elif val.startswith("<"):
        return (lambda x, t: x < t, float(val[1:]))
    else:
        # Plain number — treat as max threshold
        return (lambda x, t: x <= t, float(val))


def resolve_filter(prop, raw_value):
    """Resolve named presets then parse the filter."""
    val = raw_value.strip()
    if prop in NAMED_PRESETS and val.lower() in NAMED_PRESETS[prop]:
        val = NAMED_PRESETS[prop][val.lower()]
    return parse_filter(val)


def compute_mpo_score(row):
    """Compute a multi-parameter optimization (MPO) score for a molecule."""
    scores = []
    weights = []
    for prop, weight in MPO_WEIGHTS.items():
        val = row.get(prop)
        if pd.notna(val):
            try:
                v = float(val)
                scores.append(v * weight)
                weights.append(abs(weight))
            except (ValueError, TypeError):
                pass
    if not weights:
        return np.nan
    return np.sum(scores) / np.sum(weights)


def apply_filters(df, filters):
    """Apply all configured filters and return the passing mask."""
    mask = pd.Series(True, index=df.index)
    applied = {}

    for prop, raw_value in filters.items():
        if prop not in df.columns:
            print(f"  Warning: filter property '{prop}' not found in data, skipping", flush=True)
            continue

        parsed = resolve_filter(prop, raw_value)
        if parsed is None:
            continue

        op_func, threshold = parsed
        col = pd.to_numeric(df[prop], errors="coerce")
        prop_mask = col.apply(lambda x: op_func(x, threshold) if pd.notna(x) else True)
        mask &= prop_mask
        applied[prop] = raw_value
        print(f"  Filter {prop} {raw_value}: {prop_mask.sum()}/{len(df)} pass", flush=True)

    return mask, applied


def main():
    if not os.path.exists(INPUT_PATH):
        print(f"ERROR: raw_predictions.csv not found in inputs/", flush=True)
        sys.exit(1)

    print("Loading raw predictions...", flush=True)
    df = pd.read_csv(INPUT_PATH)
    print(f"Total molecules: {len(df)}", flush=True)

    # Apply filters
    print("Applying filters...", flush=True)
    mask, applied_filters = apply_filters(df, FILTER_PARAMS)
    filtered = df[mask].copy()
    print(f"Molecules passing all filters: {len(filtered)}/{len(df)}", flush=True)

    # Compute MPO scores
    print("Computing MPO scores...", flush=True)
    filtered["mpo_score"] = filtered.apply(compute_mpo_score, axis=1)

    # Rank by MPO score descending
    filtered = filtered.sort_values("mpo_score", ascending=False).reset_index(drop=True)

    # Select top N
    top = filtered.head(TOP_N)
    print(f"Top {TOP_N} candidates selected (returning {len(top)}).", flush=True)

    # Write output CSV
    top.to_csv(OUTPUT_CSV, index=False)

    # Write analysis summary JSON
    summary = {
        "total_molecules": len(df),
        "passed_filters": int(mask.sum()),
        "top_n": TOP_N,
        "candidates_returned": len(top),
        "filters_applied": applied_filters,
        "mpo_weights": MPO_WEIGHTS,
        "top_candidates": [
            {
                "smiles": row.get("smiles", ""),
                "name": row.get("name", ""),
                "mpo_score": round(float(row["mpo_score"]), 4) if pd.notna(row["mpo_score"]) else None,
            }
            for _, row in top.iterrows()
        ],
    }
    with open(OUTPUT_JSON, "w") as f:
        json.dump(summary, f, indent=2)

    print("Analysis complete.", flush=True)


if __name__ == "__main__":
    main()
