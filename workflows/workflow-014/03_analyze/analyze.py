#!/usr/bin/env python3
"""Filter and rank molecules by ADMET properties."""

import csv
import json
import os
import sys

# ADMET property thresholds for favorable drug-likeness.
# Classification properties: threshold is the max acceptable probability of a negative outcome,
# or min acceptable probability of a positive outcome.
FILTER_CRITERIA = {
    # Positive outcomes (higher is better) — keep if predicted probability >= threshold
    "BBB_Martins": {"min": 0.5, "label": "BBB Permeability"},
    "Bioavailability_Ma": {"min": 0.5, "label": "Bioavailability"},
    "HIA_Hou": {"min": 0.5, "label": "Human Intestinal Absorption"},
    # Negative outcomes (lower is better) — keep if predicted probability <= threshold
    "hERG": {"max": 0.5, "label": "hERG Toxicity"},
    "DILI": {"max": 0.5, "label": "Drug-Induced Liver Injury"},
    "AMES": {"max": 0.5, "label": "AMES Mutagenicity"},
    "ClinTox": {"max": 0.3, "label": "Clinical Toxicity"},
}

# Properties used for ranking (higher is better for the composite score)
RANKING_WEIGHTS = {
    "BBB_Martins": 1.0,
    "Bioavailability_Ma": 1.0,
    "HIA_Hou": 1.0,
    "hERG": -1.0,        # lower is better
    "DILI": -1.0,        # lower is better
    "AMES": -1.0,        # lower is better
    "ClinTox": -1.0,     # lower is better
    "Lipinski": 0.5,     # higher is better
    "QED": 1.0,          # higher is better
}


def load_params():
    local_params = {}
    if os.path.exists("params.json"):
        with open("params.json") as f:
            local_params = json.load(f)
    return local_params


def compute_score(row):
    """Compute a composite ADMET favorability score."""
    score = 0.0
    count = 0
    for prop, weight in RANKING_WEIGHTS.items():
        val = row.get(prop)
        if val is not None and val != "":
            try:
                score += float(val) * weight
                count += 1
            except ValueError:
                pass
    return score / max(count, 1)


def passes_filters(row):
    """Check if a molecule passes all ADMET filter criteria."""
    for prop, criteria in FILTER_CRITERIA.items():
        val = row.get(prop)
        if val is None or val == "":
            continue
        try:
            val = float(val)
        except ValueError:
            continue
        if "min" in criteria and val < criteria["min"]:
            return False
        if "max" in criteria and val > criteria["max"]:
            return False
    return True


def main():
    params = load_params()
    top_n = int(params.get("top_n", 20))

    input_path = "./inputs/admet_predictions.csv"
    if not os.path.exists(input_path):
        print("ERROR: admet_predictions.csv not found in inputs/", flush=True)
        sys.exit(1)

    print(f"Loading ADMET predictions...", flush=True)

    with open(input_path, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        all_rows = list(reader)

    print(f"Total molecules: {len(all_rows)}", flush=True)

    # Apply filters
    filtered = [row for row in all_rows if passes_filters(row)]
    print(f"Molecules passing filters: {len(filtered)}", flush=True)

    # Score and rank
    for row in filtered:
        row["_admet_score"] = compute_score(row)

    filtered.sort(key=lambda r: r["_admet_score"], reverse=True)

    # Select top N
    top_candidates = filtered[:top_n]
    print(f"Top {top_n} candidates selected.", flush=True)

    # Write filtered candidates
    os.makedirs("outputs", exist_ok=True)
    output_fields = list(fieldnames) + ["_admet_score"]
    with open("outputs/filtered_candidates.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=output_fields)
        writer.writeheader()
        writer.writerows(top_candidates)

    # Write analysis summary
    summary = {
        "total_molecules": len(all_rows),
        "passed_filters": len(filtered),
        "top_n": top_n,
        "candidates_returned": len(top_candidates),
        "filter_criteria": {
            prop: {**criteria} for prop, criteria in FILTER_CRITERIA.items()
        },
        "top_candidates": [
            {
                "smiles": row.get("smiles", ""),
                "name": row.get("name", ""),
                "score": round(row["_admet_score"], 4),
            }
            for row in top_candidates
        ],
    }
    with open("outputs/analysis_summary.json", "w") as f:
        json.dump(summary, f, indent=2)

    print("Analysis complete.", flush=True)


if __name__ == "__main__":
    main()
