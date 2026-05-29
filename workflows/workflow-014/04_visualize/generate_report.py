#!/usr/bin/env python3
"""Generate a self-contained HTML dashboard with Bootstrap 5 + Plotly.js CDN."""

import csv
import json
import os
import sys

# =============================================================================
# CONFIGURATION
# =============================================================================

CANDIDATES_PATH = "./inputs/filtered_candidates.csv"
OUTPUT_HTML = "./report.html"

# ADMET properties for the radar chart (classification probabilities, 0-1 range)
RADAR_PROPS = [
    ("BBB_Martins", "BBB Permeability"),
    ("Bioavailability_Ma", "Bioavailability"),
    ("HIA_Hou", "Intestinal Absorption"),
    ("hERG", "hERG Risk"),
    ("DILI", "DILI Risk"),
    ("AMES", "AMES Mutagenicity"),
    ("ClinTox", "Clinical Toxicity"),
    ("Lipinski", "Lipinski"),
    ("QED", "QED"),
]

# All ADMET properties to show in the data table
TABLE_PROPS = [
    "AMES", "BBB_Martins", "Bioavailability_Ma", "ClinTox", "DILI",
    "HIA_Hou", "hERG", "Lipinski", "QED",
    "CYP1A2_Veith", "CYP2C19_Veith", "CYP2C9_Veith",
    "CYP2D6_Veith", "CYP3A4_Veith",
    "Caco2_Wang", "Clearance_Hepatocyte_AZ", "Half_Life_Obach",
    "Lipophilicity_AstraZeneca", "PPBR_AZ", "Solubility_AqSolDB",
]

# Properties where lower is better (for color coding)
LOWER_IS_BETTER = {"AMES", "hERG", "DILI", "ClinTox"}

# Key properties for the sortable summary table (column_key, display_label, tooltip)
SUMMARY_PROPS = [
    ("mpo_score",          "MPO Score",  "Multi-Parameter Optimization score — composite of all key ADMET properties. Higher is better."),
    ("hERG",               "hERG",       "hERG cardiac toxicity — probability of hERG channel blockade (IC50 ≤ 10 µM). Model classification boundary: 0.5. Green: <0.3 (low risk), Yellow: 0.3–0.5 (moderate), Red: ≥0.5 (high risk)."),
    ("DILI",               "DILI",       "Drug-Induced Liver Injury — probability based on DILIst dataset (Thakkar et al. 2020). Model classification boundary: 0.5. Green: <0.3, Yellow: 0.3–0.5, Red: ≥0.5."),
    ("AMES",               "AMES",       "Ames mutagenicity — probability of bacterial mutagenicity (ICH S2(R1)). Model classification boundary: 0.5. Green: <0.3, Yellow: 0.3–0.5, Red: ≥0.5."),
    ("ClinTox",            "ClinTox",    "Clinical Toxicity — probability of FDA clinical trial failure due to toxicity (Gayvert et al. 2016). Model classification boundary: 0.5. Green: <0.3, Yellow: 0.3–0.5, Red: ≥0.5."),
    ("BBB_Martins",        "BBB",        "Blood-Brain Barrier permeability — probability of CNS penetration (Martins et al. 2012, logBB ≥ −1 threshold). High needed for CNS drugs; low preferred for peripheral targets. Green: >0.7, Yellow: 0.5–0.7, Red: ≤0.5."),
    ("HIA_Hou",            "HIA",        "Human Intestinal Absorption — probability of oral absorption (Hou et al. 2007, ≥30% absorption threshold). Green: >0.7 (high confidence), Yellow: 0.5–0.7, Red: ≤0.5."),
    ("QED",                "QED",        "Quantitative Estimate of Druglikeness (Bickerton et al. 2012, Nature Chemistry). Green: ≥0.67 (drug-like), Yellow: 0.34–0.67 (borderline), Red: <0.34 (non-drug-like). Median of FDA-approved oral drugs: 0.49."),
    ("Solubility_AqSolDB", "Solubility", "Aqueous solubility (log mol/L) — higher (less negative) is better. Acceptable oral solubility: LogS > −4 (BCS guideline). No color coding applied as values are on a log scale."),
]


def load_candidates(path):
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        return list(reader)


def safe_float(val, default=None):
    if val is None or val == "":
        return default
    try:
        return float(val)
    except (ValueError, TypeError):
        return default


def build_radar_traces_js(candidates, max_traces=5):
    """Build Plotly.js radar trace objects for top N candidates."""
    traces = []
    for i, row in enumerate(candidates[:max_traces]):
        name = row.get("name", "") or row.get("smiles", f"Candidate {i+1}")
        if len(name) > 30:
            name = name[:27] + "..."

        r_vals = []
        theta_labels = []
        for prop_key, prop_label in RADAR_PROPS:
            val = safe_float(row.get(prop_key))
            if val is not None:
                # For risk properties, invert so that "good" = high on the radar
                if prop_key in LOWER_IS_BETTER:
                    val = 1.0 - val
                r_vals.append(round(val, 3))
            else:
                r_vals.append(0)
            theta_labels.append(prop_label)

        # Close the polygon
        r_vals.append(r_vals[0])
        theta_labels.append(theta_labels[0])

        traces.append({
            "type": "scatterpolar",
            "r": r_vals,
            "theta": theta_labels,
            "fill": "toself",
            "name": name,
            "opacity": 0.6,
        })
    return json.dumps(traces)


def build_table_html(candidates):
    """Build the HTML table rows for all candidates."""
    headers = "".join(f"<th>{p}</th>" for p in TABLE_PROPS)
    header_row = f"<tr><th>#</th><th>Name</th><th>SMILES</th><th>MPO Score</th>{headers}</tr>"

    rows = []
    for i, row in enumerate(candidates, 1):
        name = row.get("name", "N/A")
        smiles = row.get("smiles", "")
        score = safe_float(row.get("mpo_score"))
        score_str = f"{score:.4f}" if score is not None else "N/A"

        cells = ""
        for prop in TABLE_PROPS:
            val = safe_float(row.get(prop))
            if val is not None:
                if prop in LOWER_IS_BETTER:
                    color = "#198754" if val < 0.3 else ("#ffc107" if val < 0.5 else "#dc3545")
                elif prop in ("BBB_Martins", "Bioavailability_Ma", "HIA_Hou"):
                    color = "#198754" if val > 0.7 else ("#ffc107" if val > 0.5 else "#dc3545")
                else:
                    color = "#6c757d"
                cells += f'<td style="color:{color}">{val:.3f}</td>'
            else:
                cells += "<td class='text-muted'>-</td>"

        rows.append(
            f"<tr><td>{i}</td><td>{name}</td>"
            f"<td class='smiles-cell'>{smiles}</td>"
            f"<td><strong>{score_str}</strong></td>{cells}</tr>"
        )

    return header_row, "\n".join(rows)


def build_summary_table_html(candidates):
    """Build a compact sortable summary table with key ADMET properties and tooltips."""
    headers = "".join(
        f'<th data-tip="{tip}" style="cursor:pointer;white-space:nowrap">{label} ↕</th>'
        for _, label, tip in SUMMARY_PROPS
    )
    header_row = f"<tr><th>#</th><th>Name</th>{headers}</tr>"

    rows = []
    for i, row in enumerate(candidates, 1):
        name = row.get("name", "N/A")
        cells = ""
        for prop, _, _ in SUMMARY_PROPS:
            val = safe_float(row.get(prop))
            if val is not None:
                if prop in LOWER_IS_BETTER:
                    color = "#198754" if val < 0.3 else ("#ffc107" if val < 0.5 else "#dc3545")
                elif prop in ("BBB_Martins", "HIA_Hou"):
                    # Martins 2012 / Hou 2007: model boundary = 0.5; >0.7 = high confidence
                    color = "#198754" if val > 0.7 else ("#ffc107" if val > 0.5 else "#dc3545")
                elif prop == "QED":
                    # Bickerton 2012 (Nature Chemistry): >=0.67 drug-like, <0.34 non-drug-like
                    color = "#198754" if val >= 0.67 else ("#ffc107" if val >= 0.34 else "#dc3545")
                elif prop == "mpo_score":
                    color = "#198754" if val > 0.6 else ("#ffc107" if val > 0.4 else "#dc3545")
                else:
                    color = "#6c757d"
                cells += f'<td style="color:{color}">{val:.3f}</td>'
            else:
                cells += "<td class='text-muted'>-</td>"
        rows.append(f"<tr><td>{i}</td><td>{name}</td>{cells}</tr>")

    return header_row, "\n".join(rows)


def build_summary_cards(candidates):
    """Build summary statistics for the dashboard header."""
    total = len(candidates)
    scores = [safe_float(r.get("mpo_score")) for r in candidates]
    scores = [s for s in scores if s is not None]
    avg_score = sum(scores) / len(scores) if scores else 0
    max_score = max(scores) if scores else 0

    # Count how many pass common safety thresholds
    safe_count = 0
    for row in candidates:
        herg = safe_float(row.get("hERG"))
        ames = safe_float(row.get("AMES"))
        dili = safe_float(row.get("DILI"))
        if herg is not None and ames is not None and dili is not None:
            if herg < 0.5 and ames < 0.5 and dili < 0.5:
                safe_count += 1

    return total, avg_score, max_score, safe_count


def main():
    if not os.path.exists(CANDIDATES_PATH):
        print("ERROR: filtered_candidates.csv not found in inputs/", flush=True)
        sys.exit(1)

    candidates = load_candidates(CANDIDATES_PATH)
    print(f"Loaded {len(candidates)} candidates", flush=True)

    total, avg_score, max_score, safe_count = build_summary_cards(candidates)
    radar_traces_js = build_radar_traces_js(candidates)
    header_row, table_body = build_table_html(candidates)
    summary_header_row, summary_table_body = build_summary_table_html(candidates)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ADMET-AI Prediction Dashboard</title>
<link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  body {{ background: #f8f9fa; }}
  .smiles-cell {{ font-family: monospace; font-size: 0.8em; max-width: 220px;
                  overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
  .stat-value {{ font-size: 2rem; font-weight: 700; }}
  .table-sm td, .table-sm th {{ padding: 0.4rem 0.5rem; font-size: 0.85rem; }}
  .table thead th {{ position: sticky; top: 0; z-index: 1; }}
  #radar-chart {{ min-height: 500px; }}
</style>
</head>
<body>
<div class="container-fluid py-4">
  <div class="row mb-4">
    <div class="col">
      <h1 class="h3">ADMET-AI Prediction Dashboard</h1>
      <p class="text-muted mb-0">Filtered and ranked candidates with multi-parameter optimization scores</p>
    </div>
  </div>

  <!-- Summary Cards -->
  <div class="row g-3 mb-4">
    <div class="col-md-3">
      <div class="card text-center">
        <div class="card-body">
          <div class="stat-value text-primary">{total}</div>
          <div class="text-muted">Candidates</div>
        </div>
      </div>
    </div>
    <div class="col-md-3">
      <div class="card text-center">
        <div class="card-body">
          <div class="stat-value text-success">{avg_score:.3f}</div>
          <div class="text-muted">Avg MPO Score</div>
        </div>
      </div>
    </div>
    <div class="col-md-3">
      <div class="card text-center">
        <div class="card-body">
          <div class="stat-value text-info">{max_score:.3f}</div>
          <div class="text-muted">Best MPO Score</div>
        </div>
      </div>
    </div>
    <div class="col-md-3">
      <div class="card text-center">
        <div class="card-body">
          <div class="stat-value text-warning">{safe_count}</div>
          <div class="text-muted">Pass Safety Panel</div>
        </div>
      </div>
    </div>
  </div>

  <!-- Summary Table -->
  <div class="row mb-4">
    <div class="col-12">
      <div class="card">
        <div class="card-header d-flex justify-content-between align-items-center">
          <h5 class="card-title mb-0">Summary — Key ADMET Properties</h5>
          <small class="text-muted">Click column headers to sort &nbsp;·&nbsp; Hover headers for property descriptions</small>
        </div>
        <div class="card-body p-0">
          <div class="table-responsive">
            <table id="summary-table" class="table table-sm table-striped table-hover mb-0">
              <thead class="table-dark">
                {summary_header_row}
              </thead>
              <tbody>
                {summary_table_body}
              </tbody>
            </table>
          </div>
        </div>
      </div>
    </div>
  </div>

  <!-- Radar Chart -->
  <div class="row mb-4">
    <div class="col-12">
      <div class="card">
        <div class="card-header">
          <h5 class="card-title mb-0">Top Candidates — ADMET Radar</h5>
          <small class="text-muted">Risk properties (hERG, DILI, AMES, ClinTox) are inverted so higher = better for all axes</small>
        </div>
        <div class="card-body">
          <div id="radar-chart"></div>
        </div>
      </div>
    </div>
  </div>

  <!-- Data Table -->
  <div class="row">
    <div class="col-12">
      <div class="card">
        <div class="card-header">
          <h5 class="card-title mb-0">All Candidates</h5>
        </div>
        <div class="card-body p-0">
          <div class="table-responsive">
            <table class="table table-sm table-striped table-hover mb-0">
              <thead class="table-dark">
                {header_row}
              </thead>
              <tbody>
                {table_body}
              </tbody>
            </table>
          </div>
        </div>
      </div>
    </div>
  </div>
</div>

<script>
  var traces = {radar_traces_js};
  var layout = {{
    polar: {{
      radialaxis: {{
        visible: true,
        range: [0, 1]
      }}
    }},
    showlegend: true,
    legend: {{ orientation: "h", y: -0.1 }},
    margin: {{ t: 30, b: 60, l: 60, r: 60 }}
  }};
  Plotly.newPlot('radar-chart', traces, layout, {{ responsive: true }});

  // Sortable summary table — sync data-tip to title, then wire click handlers
  document.querySelectorAll('#summary-table thead th').forEach(function(th) {{
    if (th.dataset.tip) th.title = th.dataset.tip;
    th.addEventListener('click', function() {{
      var tbody = th.closest('table').querySelector('tbody');
      var col = th.cellIndex;
      var asc = th.dataset.dir !== 'asc';
      th.dataset.dir = asc ? 'asc' : 'desc';
      // Update sort indicator; restore title from data-tip after textContent change
      th.closest('tr').querySelectorAll('th').forEach(function(t) {{
        t.textContent = t.textContent.replace(/ [↕↑↓]$/, '') + ' ↕';
        if (t.dataset.tip) t.title = t.dataset.tip;
      }});
      th.textContent = th.textContent.replace(/ [↕↑↓]$/, '') + (asc ? ' ↑' : ' ↓');
      if (th.dataset.tip) th.title = th.dataset.tip;
      var rows = Array.from(tbody.rows);
      rows.sort(function(a, b) {{
        var av = parseFloat(a.cells[col].textContent);
        var bv = parseFloat(b.cells[col].textContent);
        if (isNaN(av)) av = a.cells[col].textContent;
        if (isNaN(bv)) bv = b.cells[col].textContent;
        return asc ? (av > bv ? 1 : -1) : (av < bv ? 1 : -1);
      }});
      rows.forEach(function(r) {{ tbody.appendChild(r); }});
    }});
  }});
</script>
</body>
</html>"""

    with open(OUTPUT_HTML, "w") as f:
        f.write(html)

    print(f"Dashboard generated: {OUTPUT_HTML} ({len(candidates)} candidates)", flush=True)


if __name__ == "__main__":
    main()
