#!/usr/bin/env python3
"""Generate an interactive HTML dashboard for ADMET predictions."""

import csv
import json
import os
import sys


def main():
    summary_path = "./inputs/analysis_summary.json"
    candidates_path = "./inputs/filtered_candidates.csv"

    if not os.path.exists(summary_path):
        print("ERROR: analysis_summary.json not found in inputs/", flush=True)
        sys.exit(1)
    if not os.path.exists(candidates_path):
        print("ERROR: filtered_candidates.csv not found in inputs/", flush=True)
        sys.exit(1)

    with open(summary_path) as f:
        summary = json.load(f)

    with open(candidates_path, newline="") as f:
        reader = csv.DictReader(f)
        candidates = list(reader)

    # ADMET properties to display in the dashboard
    admet_props = [
        "AMES", "BBB_Martins", "Bioavailability_Ma", "ClinTox", "DILI",
        "HIA_Hou", "hERG", "Lipinski", "QED",
        "CYP1A2_Veith", "CYP2C19_Veith", "CYP2C9_Veith",
        "CYP2D6_Veith", "CYP3A4_Veith",
        "Caco2_Wang", "Clearance_Hepatocyte_AZ", "Half_Life_Obach",
        "Lipophilicity_AstraZeneca", "PPBR_AZ", "Solubility_AqSolDB",
    ]

    # Build table rows
    table_rows = []
    for i, row in enumerate(candidates, 1):
        name = row.get("name", "N/A")
        smiles = row.get("smiles", "")
        score = row.get("_admet_score", "N/A")
        try:
            score = f"{float(score):.4f}"
        except (ValueError, TypeError):
            pass

        props_html = ""
        for prop in admet_props:
            val = row.get(prop, "")
            if val:
                try:
                    fval = float(val)
                    # Color-code: green for favorable, red for unfavorable
                    if prop in ("AMES", "hERG", "DILI", "ClinTox"):
                        color = "#28a745" if fval < 0.3 else ("#ffc107" if fval < 0.5 else "#dc3545")
                    elif prop in ("BBB_Martins", "Bioavailability_Ma", "HIA_Hou"):
                        color = "#28a745" if fval > 0.7 else ("#ffc107" if fval > 0.5 else "#dc3545")
                    else:
                        color = "#6c757d"
                    props_html += f'<td style="color:{color}">{fval:.3f}</td>'
                except ValueError:
                    props_html += f"<td>{val}</td>"
            else:
                props_html += "<td>-</td>"

        table_rows.append(
            f"<tr><td>{i}</td><td>{name}</td>"
            f"<td class='smiles'>{smiles}</td>"
            f"<td><b>{score}</b></td>{props_html}</tr>"
        )

    prop_headers = "".join(f"<th>{p}</th>" for p in admet_props)
    table_body = "\n".join(table_rows)

    # Filter criteria summary
    filter_html = ""
    for prop, criteria in summary.get("filter_criteria", {}).items():
        label = criteria.get("label", prop)
        if "min" in criteria:
            filter_html += f"<li><b>{label}</b> ({prop}): &ge; {criteria['min']}</li>"
        if "max" in criteria:
            filter_html += f"<li><b>{label}</b> ({prop}): &le; {criteria['max']}</li>"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>ADMET-AI Prediction Report</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
         margin: 0; padding: 20px; background: #f5f5f5; color: #333; }}
  .container {{ max-width: 1400px; margin: 0 auto; }}
  h1 {{ color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }}
  h2 {{ color: #34495e; margin-top: 30px; }}
  .stats {{ display: flex; gap: 20px; flex-wrap: wrap; margin: 20px 0; }}
  .stat-card {{ background: white; padding: 20px; border-radius: 8px;
               box-shadow: 0 2px 4px rgba(0,0,0,0.1); min-width: 180px; }}
  .stat-card .value {{ font-size: 2em; font-weight: bold; color: #3498db; }}
  .stat-card .label {{ color: #7f8c8d; font-size: 0.9em; }}
  table {{ width: 100%; border-collapse: collapse; background: white;
           box-shadow: 0 2px 4px rgba(0,0,0,0.1); border-radius: 8px;
           overflow: hidden; font-size: 0.85em; }}
  th {{ background: #2c3e50; color: white; padding: 10px 8px; text-align: left;
        position: sticky; top: 0; }}
  td {{ padding: 8px; border-bottom: 1px solid #ecf0f1; }}
  tr:hover {{ background: #f0f7ff; }}
  .smiles {{ font-family: monospace; font-size: 0.8em; max-width: 250px;
             overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }}
  ul {{ columns: 2; }}
  li {{ margin-bottom: 5px; }}
  .table-wrap {{ overflow-x: auto; }}
</style>
</head>
<body>
<div class="container">
  <h1>ADMET-AI Prediction Report</h1>

  <div class="stats">
    <div class="stat-card">
      <div class="value">{summary.get('total_molecules', 'N/A')}</div>
      <div class="label">Total Molecules</div>
    </div>
    <div class="stat-card">
      <div class="value">{summary.get('passed_filters', 'N/A')}</div>
      <div class="label">Passed Filters</div>
    </div>
    <div class="stat-card">
      <div class="value">{summary.get('candidates_returned', 'N/A')}</div>
      <div class="label">Top Candidates</div>
    </div>
  </div>

  <h2>Filter Criteria</h2>
  <ul>{filter_html}</ul>

  <h2>Top Candidates</h2>
  <div class="table-wrap">
  <table>
    <thead>
      <tr>
        <th>#</th><th>Name</th><th>SMILES</th><th>Score</th>
        {prop_headers}
      </tr>
    </thead>
    <tbody>
      {table_body}
    </tbody>
  </table>
  </div>
</div>
</body>
</html>"""

    os.makedirs("outputs", exist_ok=True)
    with open("outputs/report.html", "w") as f:
        f.write(html)

    # Write data.json for programmatic access
    data_json = {
        "summary": summary,
        "candidates": [
            {
                "rank": i + 1,
                "name": row.get("name", ""),
                "smiles": row.get("smiles", ""),
                "score": row.get("_admet_score", ""),
                **{prop: row.get(prop, "") for prop in admet_props},
            }
            for i, row in enumerate(candidates)
        ],
    }
    with open("outputs/data.json", "w") as f:
        json.dump(data_json, f, indent=2)

    print(f"Report generated: outputs/report.html ({len(candidates)} candidates)", flush=True)
    print(f"Data exported: outputs/data.json", flush=True)


if __name__ == "__main__":
    main()
