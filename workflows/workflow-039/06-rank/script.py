#!/usr/bin/env python3
"""Node 06 — Rank Panel.

Join every per-clone measurement into one table and sort it, so the output of this workflow is
a decision — which clones to express — rather than five CSVs the user has to reconcile.

Ranking policy, in order:

  1. **Red flags, ascending.** A clone with a value never observed among clinical-stage
     therapeutics sorts below one without.
  2. **Exposed high-severity liabilities, ascending.** An exposed N-glycosylation sequon or free
     cysteine counts; a buried one does not (node 03 supplies the exposure).
  3. **Amber flags, ascending.**
  4. **Peptide content, descending.** More human first.
  5. **Clone id**, so the order is deterministic.

Nothing is filtered. Every clone that entered the panel leaves it, including the ones that
failed annotation or modelling, with a status column saying so. The ranking is a sort over
columns that are all present in the output, so a user who disagrees with this weighting can
re-sort the CSV on their own terms — which is why there is no single opaque composite score.

`top_n` marks a shortlist; it does not truncate the table.
"""
import csv
import json
import os
import sys

FLAG_COLUMNS = [
    "flag_total_cdr_length",
    "flag_psh_cdr_vicinity",
    "flag_ppc_cdr_vicinity",
    "flag_pnc_cdr_vicinity",
    "flag_sfvcsp",
]


def fail(message):
    print(f"ERROR: {message}", flush=True)
    sys.exit(1)


def warn(message):
    print(f"WARNING: {message}", flush=True)


def read_csv(path, label, required=True):
    if not os.path.exists(path):
        if required:
            fail(f"{label} not found: {path}. The upstream node must run first.")
        warn(f"{label} not found ({path}); continuing without it.")
        return []
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


def read_json(path, label, required=False):
    if not os.path.exists(path):
        if required:
            fail(f"{label} not found: {path}.")
        return {}
    try:
        with open(path) as handle:
            return json.load(handle)
    except Exception as exc:
        warn(f"Could not read {label} ('{path}'): {exc}")
        return {}


def index_by_clone(rows):
    return {row["clone_id"]: row for row in rows if row.get("clone_id")}


def to_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def to_int(value, default=0):
    try:
        return int(float(value))
    except (TypeError, ValueError):
        return default


def main_sort_key(row):
    """The ranking policy, at module level so the unit tests can exercise it directly.

    Failed clones sort last but are never dropped. Within the scorable ones the order is red
    flags, then exposed high-severity liabilities, then amber flags, then human peptide content
    descending, then clone id for determinism.
    """
    scorable = row.get("status") != "failed"
    return (
        0 if scorable else 1,
        to_int(row.get("n_red_flags"), 99),
        to_int(row.get("n_exposed_high_liabilities"), 99),
        to_int(row.get("n_amber_flags"), 99),
        -(to_float(row.get("mean_peptide_content")) or 0.0),
        row.get("clone_id", ""),
    )


OUTPUT_FIELDS = [
    "rank",
    "shortlist",
    "clone_id",
    "receptor_type",
    "status",
    "n_red_flags",
    "n_amber_flags",
    "n_exposed_high_liabilities",
    "n_exposed_liabilities",
    "n_liabilities",
    "mean_peptide_content",
    "mean_v_identity",
    "total_cdr_length",
    "cdrh3",
    "cdrh3_length",
    "psh_cdr_vicinity",
    "ppc_cdr_vicinity",
    "pnc_cdr_vicinity",
    "sfvcsp",
    "isoelectric_point",
    "flag_total_cdr_length",
    "flag_psh_cdr_vicinity",
    "flag_ppc_cdr_vicinity",
    "flag_pnc_cdr_vicinity",
    "flag_sfvcsp",
    "thresholds_calibrated",
    "n_proposals",
    "model_file",
    "refinement",
    "mean_error_angstrom",
    "cdrh3_mean_error_angstrom",
    "vh_v_gene",
    "vl_v_gene",
    "germline_species",
    "reason",
]


def main():
    os.makedirs("./outputs", exist_ok=True)
    top_n = int(os.environ.get("PARAM_TOP_N", "10"))
    if top_n < 1:
        fail(f"top_n must be at least 1, got {top_n}")

    annotations = index_by_clone(read_csv("./inputs/annotations.csv", "annotations.csv"))
    developability = index_by_clone(
        read_csv("./inputs/developability.csv", "developability.csv", required=False)
    )
    flags = index_by_clone(read_csv("./inputs/flags.csv", "flags.csv", required=False))
    humanness = index_by_clone(read_csv("./inputs/humanness.csv", "humanness.csv", required=False))
    models = index_by_clone(
        read_csv("./inputs/model_summary.csv", "model_summary.csv", required=False)
    )
    liabilities = read_csv("./inputs/liabilities.csv", "liabilities.csv", required=False)
    proposals = read_csv("./inputs/proposals.csv", "proposals.csv", required=False)

    if not annotations:
        fail("annotations.csv contains no rows; node 01 must run first.")

    liability_counts = {}
    for hit in liabilities:
        clone_id = hit.get("clone_id")
        if clone_id:
            liability_counts[clone_id] = liability_counts.get(clone_id, 0) + 1

    proposal_counts = {}
    for row in proposals:
        if row.get("verdict") == "kept" and row.get("clone_id"):
            proposal_counts[row["clone_id"]] = proposal_counts.get(row["clone_id"], 0) + 1

    rows = []
    for clone_id, annotation in annotations.items():
        metrics = developability.get(clone_id, {})
        flag_row = flags.get(clone_id, {})
        human = humanness.get(clone_id, {})
        model = models.get(clone_id, {})

        statuses = {
            "annotate": annotation.get("status", ""),
            "model": model.get("status", ""),
            "developability": metrics.get("status", ""),
            "humanness": human.get("status", ""),
        }
        failed_stages = [stage for stage, state in statuses.items() if state == "failed"]
        reasons = [
            f"{stage}: {source.get('reason')}"
            for stage, source in (
                ("annotate", annotation),
                ("model", model),
                ("developability", metrics),
                ("humanness", human),
            )
            if source.get("status") == "failed" and source.get("reason")
        ]

        if annotation.get("status") == "failed":
            status = "failed"
        elif failed_stages:
            status = "partial"
        else:
            status = "ok"

        flag_values = [flag_row.get(column, "") for column in FLAG_COLUMNS]
        n_red = (
            to_int(flag_row.get("n_red"), 0)
            if flag_row.get("n_red") not in (None, "")
            else sum(1 for value in flag_values if value == "red")
        )
        n_amber = (
            to_int(flag_row.get("n_amber"), 0)
            if flag_row.get("n_amber") not in (None, "")
            else sum(1 for value in flag_values if value == "amber")
        )

        rows.append(
            {
                "clone_id": clone_id,
                "receptor_type": annotation.get("receptor_type") or metrics.get("receptor_type", ""),
                "status": status,
                "n_red_flags": n_red if status != "failed" else "",
                "n_amber_flags": n_amber if status != "failed" else "",
                "n_exposed_high_liabilities": metrics.get("n_exposed_high_liabilities", ""),
                "n_exposed_liabilities": metrics.get("n_exposed_liabilities", ""),
                "n_liabilities": liability_counts.get(clone_id, 0),
                "mean_peptide_content": human.get("mean_peptide_content", ""),
                "mean_v_identity": human.get("mean_v_identity", ""),
                "total_cdr_length": annotation.get("total_cdr_length", ""),
                "cdrh3": annotation.get("cdrh3", ""),
                "cdrh3_length": annotation.get("cdrh3_length", ""),
                "psh_cdr_vicinity": metrics.get("psh_cdr_vicinity", ""),
                "ppc_cdr_vicinity": metrics.get("ppc_cdr_vicinity", ""),
                "pnc_cdr_vicinity": metrics.get("pnc_cdr_vicinity", ""),
                "sfvcsp": metrics.get("sfvcsp", ""),
                "isoelectric_point": metrics.get("isoelectric_point", ""),
                "flag_total_cdr_length": flag_row.get("flag_total_cdr_length", ""),
                "flag_psh_cdr_vicinity": flag_row.get("flag_psh_cdr_vicinity", ""),
                "flag_ppc_cdr_vicinity": flag_row.get("flag_ppc_cdr_vicinity", ""),
                "flag_pnc_cdr_vicinity": flag_row.get("flag_pnc_cdr_vicinity", ""),
                "flag_sfvcsp": flag_row.get("flag_sfvcsp", ""),
                "thresholds_calibrated": metrics.get("thresholds_calibrated", ""),
                "n_proposals": proposal_counts.get(clone_id, 0),
                "model_file": model.get("model_file", ""),
                "refinement": model.get("refinement", ""),
                "mean_error_angstrom": model.get("mean_error_angstrom", ""),
                "cdrh3_mean_error_angstrom": model.get("cdrh3_mean_error_angstrom", ""),
                "vh_v_gene": annotation.get("vh_v_gene", ""),
                "vl_v_gene": annotation.get("vl_v_gene", ""),
                "germline_species": human.get("germline_species", ""),
                "reason": "; ".join(reasons),
                # Carried-through extra columns from the input panel.
                **{
                    key: value
                    for key, value in annotation.items()
                    if key not in set(OUTPUT_FIELDS) and key not in {"status", "reason"}
                },
            }
        )

    rows.sort(key=main_sort_key)
    for position, row in enumerate(rows, start=1):
        if row["status"] == "failed":
            row["rank"] = ""
            row["shortlist"] = ""
        else:
            row["rank"] = position
            row["shortlist"] = "yes" if position <= top_n else ""

    extra_fields = [
        key for row in rows for key in row if key not in set(OUTPUT_FIELDS)
    ]
    ordered_extras = sorted(set(extra_fields))
    fields = OUTPUT_FIELDS + ordered_extras

    with open("./outputs/ranked_panel.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    scorable = [row for row in rows if row["status"] != "failed"]
    shortlist = [row for row in rows if row["shortlist"] == "yes"]
    contents = [
        value
        for value in (to_float(row["mean_peptide_content"]) for row in scorable)
        if value is not None
    ]

    summary = {
        "top_n": top_n,
        "n_clones": len(rows),
        "n_scorable": len(scorable),
        "n_failed": sum(1 for row in rows if row["status"] == "failed"),
        "n_partial": sum(1 for row in rows if row["status"] == "partial"),
        "n_shortlisted": len(shortlist),
        "shortlist": [row["clone_id"] for row in shortlist],
        "clean_clones": [
            row["clone_id"]
            for row in scorable
            if to_int(row["n_red_flags"], 0) == 0
            and to_int(row["n_amber_flags"], 0) == 0
            and to_int(row["n_exposed_high_liabilities"], 0) == 0
        ],
        "flag_totals": {
            column: {
                state: sum(1 for row in scorable if row.get(column) == state)
                for state in ("green", "amber", "red")
            }
            for column in FLAG_COLUMNS
        },
        "n_red_total": sum(to_int(row["n_red_flags"], 0) for row in scorable),
        "n_amber_total": sum(to_int(row["n_amber_flags"], 0) for row in scorable),
        "peptide_content": (
            {
                "mean": round(sum(contents) / len(contents), 4),
                "min": round(min(contents), 4),
                "max": round(max(contents), 4),
            }
            if contents
            else None
        ),
        "n_proposals_total": sum(to_int(row["n_proposals"], 0) for row in rows),
        "ranking_policy": [
            "fewest red flags",
            "fewest exposed high-severity liabilities",
            "fewest amber flags",
            "highest human peptide content",
            "clone id, for determinism",
        ],
        "caveat": (
            "Nothing is filtered: every clone that entered the panel is in this table. The rank "
            "is a sort over columns that are all present here, so it can be recomputed or "
            "disagreed with. Flags are TAP-style risk indicators carrying a documented "
            "calibration offset, not predictions of clinical success."
        ),
    }
    with open("./outputs/panel_summary.json", "w") as handle:
        json.dump(summary, handle, indent=2)

    print(
        f"Ranked {len(scorable)} scorable clone(s) of {len(rows)}; "
        f"{summary['n_failed']} failed, {summary['n_partial']} partial",
        flush=True,
    )
    print(f"Shortlist (top {top_n}): {', '.join(summary['shortlist']) or 'none'}", flush=True)
    if summary["clean_clones"]:
        print(
            f"{len(summary['clean_clones'])} clone(s) carry no flag and no exposed "
            "high-severity liability",
            flush=True,
        )
    else:
        print(
            "No clone is completely clean. That is common and is not a reason to stop: read "
            "which flag each clone carries rather than treating the panel as failed.",
            flush=True,
        )


if __name__ == "__main__":
    main()
