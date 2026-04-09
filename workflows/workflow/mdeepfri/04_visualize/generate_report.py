#!/usr/bin/env python3
"""Generate a self-contained HTML report from mDeepFRI prediction results.

Reads results.tsv and alignment_summary.tsv, produces a single report.html
with per-protein GO term tables and score distribution charts (Bootstrap 5 + Plotly.js).
"""

import argparse
import csv
import json
import os
import sys

RESULTS_PATH = "./inputs/results.tsv"
ALIGNMENT_PATH = "./inputs/alignment_summary.tsv"
OUTPUT_HTML = "./outputs/report.html"

MODE_LABELS = {
    "mf": "Molecular Function",
    "bp": "Biological Process",
    "cc": "Cellular Component",
    "ec": "Enzyme Commission",
}


def load_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def safe_float(val, default=0.0):
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


def build_report(rows, align_rows, min_score, top_n):
    """Organise rows by protein and mode; apply score filter."""
    by_protein = {}
    for row in rows:
        # mDeepFRI results.tsv columns vary by version; handle common variants
        protein = (
            row.get("Protein_ID") or row.get("query") or row.get("protein") or "unknown"
        )
        mode = (row.get("Mode") or row.get("mode") or "").lower()
        score = safe_float(row.get("Score") or row.get("score") or row.get("confidence"))
        go_id = row.get("GO_ID") or row.get("go_term") or row.get("GO_term") or ""
        go_name = row.get("GO_Name") or row.get("go_name") or row.get("description") or ""
        network = (row.get("Network") or row.get("network") or "cnn").lower()

        if score < min_score:
            continue

        if protein not in by_protein:
            by_protein[protein] = {}
        if mode not in by_protein[protein]:
            by_protein[protein][mode] = []
        by_protein[protein][mode].append(
            {"go_id": go_id, "go_name": go_name, "score": score, "network": network}
        )

    # Sort by score descending and cap at top_n
    for protein in by_protein:
        for mode in by_protein[protein]:
            by_protein[protein][mode].sort(key=lambda x: x["score"], reverse=True)
            by_protein[protein][mode] = by_protein[protein][mode][:top_n]

    # Build alignment lookup
    align_lookup = {}
    for row in align_rows:
        qid = row.get("query") or row.get("Protein_ID") or ""
        if qid:
            align_lookup[qid] = row

    return by_protein, align_lookup


def score_color(score):
    if score >= 0.7:
        return "#198754"  # green
    if score >= 0.4:
        return "#fd7e14"  # orange
    return "#dc3545"  # red


def protein_table_html(protein, modes_data):
    rows_html = []
    for mode, preds in sorted(modes_data.items()):
        mode_label = MODE_LABELS.get(mode, mode.upper())
        for p in preds:
            color = score_color(p["score"])
            network_badge = (
                '<span class="badge bg-info text-dark">GCN</span>'
                if p["network"] == "gcn"
                else '<span class="badge bg-secondary">CNN</span>'
            )
            rows_html.append(
                f"<tr>"
                f"<td><small class='text-muted'>{mode_label}</small></td>"
                f"<td><code>{p['go_id']}</code></td>"
                f"<td>{p['go_name']}</td>"
                f"<td style='color:{color}; font-weight:600'>{p['score']:.3f}</td>"
                f"<td>{network_badge}</td>"
                f"</tr>"
            )
    if not rows_html:
        return "<p class='text-muted'>No predictions above threshold.</p>"
    table = (
        "<table class='table table-sm table-hover table-striped'>"
        "<thead><tr><th>Mode</th><th>GO ID</th><th>Name</th><th>Score</th><th>Model</th></tr></thead>"
        "<tbody>" + "".join(rows_html) + "</tbody></table>"
    )
    return table


def scores_chart_js(by_protein):
    """Return a Plotly.js <script> block with a score distribution histogram."""
    all_scores = []
    for modes in by_protein.values():
        for preds in modes.values():
            all_scores.extend(p["score"] for p in preds)

    if not all_scores:
        return ""

    scores_json = json.dumps(all_scores)
    return f"""
<div id="scoreChart" style="height:300px"></div>
<script>
  var scores = {scores_json};
  var trace = {{
    x: scores,
    type: 'histogram',
    nbinsx: 20,
    marker: {{ color: '#0d6efd', opacity: 0.75 }}
  }};
  var layout = {{
    title: 'Prediction Score Distribution',
    xaxis: {{ title: 'Confidence Score', range: [0, 1] }},
    yaxis: {{ title: 'Count' }},
    margin: {{ t: 40 }},
    plot_bgcolor: '#f8f9fa',
    paper_bgcolor: '#ffffff'
  }};
  Plotly.newPlot('scoreChart', [trace], layout, {{responsive: true}});
</script>
"""


def build_html(by_protein, align_lookup, min_score, top_n, results_path):
    total_proteins = len(by_protein)
    total_predictions = sum(
        len(preds) for modes in by_protein.values() for preds in modes.values()
    )

    protein_sections = []
    for protein, modes_data in sorted(by_protein.items()):
        align_info = align_lookup.get(protein, {})
        align_html = ""
        if align_info:
            fields = {k: v for k, v in align_info.items() if k not in ("query", "Protein_ID")}
            if fields:
                items = " &nbsp;|&nbsp; ".join(
                    f"<b>{k}</b>: {v}" for k, v in list(fields.items())[:6]
                )
                align_html = f"<p class='small text-muted mb-1'><i>Alignment:</i> {items}</p>"

        table = protein_table_html(protein, modes_data)
        n_preds = sum(len(p) for p in modes_data.values())
        protein_sections.append(f"""
<div class="card mb-3">
  <div class="card-header d-flex justify-content-between align-items-center">
    <strong>{protein}</strong>
    <span class="badge bg-primary">{n_preds} prediction(s)</span>
  </div>
  <div class="card-body">
    {align_html}
    {table}
  </div>
</div>""")

    protein_html = "\n".join(protein_sections) if protein_sections else (
        "<div class='alert alert-warning'>No predictions above the score threshold.</div>"
    )

    chart_js = scores_chart_js(by_protein)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>mDeepFRI Results</title>
<link rel="stylesheet"
  href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css">
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
<style>
  body {{ font-family: 'Segoe UI', sans-serif; background: #f8f9fa; }}
  .hero {{ background: linear-gradient(135deg, #0d6efd 0%, #6610f2 100%);
           color: white; padding: 2rem; border-radius: 0.5rem; }}
  code {{ background: #e9ecef; padding: 1px 4px; border-radius: 3px; }}
</style>
</head>
<body>
<div class="container py-4">
  <div class="hero mb-4">
    <h1 class="mb-1">mDeepFRI Results</h1>
    <p class="mb-0">Protein Function Prediction — GO Term Annotations</p>
  </div>

  <div class="row g-3 mb-4">
    <div class="col-md-4">
      <div class="card text-center">
        <div class="card-body">
          <h2 class="text-primary">{total_proteins}</h2>
          <p class="mb-0">Proteins annotated</p>
        </div>
      </div>
    </div>
    <div class="col-md-4">
      <div class="card text-center">
        <div class="card-body">
          <h2 class="text-success">{total_predictions}</h2>
          <p class="mb-0">GO term predictions</p>
        </div>
      </div>
    </div>
    <div class="col-md-4">
      <div class="card text-center">
        <div class="card-body">
          <h2 class="text-secondary">&ge;{min_score}</h2>
          <p class="mb-0">Score threshold</p>
        </div>
      </div>
    </div>
  </div>

  {chart_js}

  <h3 class="mt-4 mb-3">Per-Protein Predictions</h3>
  {protein_html}

  <footer class="text-center text-muted mt-4 small">
    Generated by mDeepFRI silva workflow &bull; Results from <code>{os.path.basename(results_path)}</code>
  </footer>
</div>
</body>
</html>"""


def main():
    parser = argparse.ArgumentParser(description="Generate mDeepFRI HTML report")
    parser.add_argument("--min-score", type=float, default=0.3)
    parser.add_argument("--top-n", type=int, default=10)
    args = parser.parse_args()

    os.makedirs("./outputs", exist_ok=True)

    print(f"Loading {RESULTS_PATH} ...", flush=True)
    rows = load_tsv(RESULTS_PATH)
    if not rows:
        print(f"WARNING: {RESULTS_PATH} is empty or missing", flush=True)

    align_rows = load_tsv(ALIGNMENT_PATH)
    print(f"  {len(rows)} prediction rows, {len(align_rows)} alignment rows", flush=True)

    by_protein, align_lookup = build_report(rows, align_rows, args.min_score, args.top_n)
    html = build_html(by_protein, align_lookup, args.min_score, args.top_n, RESULTS_PATH)

    with open(OUTPUT_HTML, "w") as f:
        f.write(html)

    print(f"Report written: {OUTPUT_HTML} ({os.path.getsize(OUTPUT_HTML):,} bytes)", flush=True)
    print(f"  {len(by_protein)} proteins | {sum(len(p) for m in by_protein.values() for p in m.values())} predictions shown", flush=True)
    print("Visualization complete.", flush=True)


if __name__ == "__main__":
    main()
