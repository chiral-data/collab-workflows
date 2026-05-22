#!/usr/bin/env python3
import argparse
import json
import os

import plotly.graph_objects as go
from plotly.subplots import make_subplots


def before_after_figure(data):
    names = ["Antibody (Receptor)", "Antigen (Ligand)"]
    residues_before = [
        data["antibody"]["totals"]["input_residues"],
        data["antigen"]["totals"]["input_residues"],
    ]
    residues_after = [
        data["antibody"]["totals"]["kept_residues"],
        data["antigen"]["totals"]["kept_residues"],
    ]
    atoms_before = [
        data["antibody"]["totals"]["input_atoms"],
        data["antigen"]["totals"]["input_atoms"],
    ]
    atoms_after = [
        data["antibody"]["totals"]["kept_atoms"],
        data["antigen"]["totals"]["kept_atoms"],
    ]

    fig = make_subplots(rows=1, cols=2, subplot_titles=("Residues: Before vs After", "Atoms: Before vs After"))
    fig.add_trace(go.Bar(x=names, y=residues_before, name="Residues Before", marker_color="#94a3b8"), row=1, col=1)
    fig.add_trace(go.Bar(x=names, y=residues_after, name="Residues After", marker_color="#16a34a"), row=1, col=1)
    fig.add_trace(go.Bar(x=names, y=atoms_before, name="Atoms Before", marker_color="#cbd5e1"), row=1, col=2)
    fig.add_trace(go.Bar(x=names, y=atoms_after, name="Atoms After", marker_color="#2563eb"), row=1, col=2)
    fig.update_layout(barmode="group", title="Structure Preparation Quality Control", margin=dict(l=20, r=20, t=60, b=40))
    return fig


def removal_figure(data):
    labels = ["Removed Water", "Removed Hetero", "Removed Nonstandard AA"]
    ab_values = [
        data["antibody"]["totals"]["removed_water_residues"],
        data["antibody"]["totals"]["removed_hetero_residues"],
        data["antibody"]["totals"]["removed_nonstandard_aa_residues"],
    ]
    ag_values = [
        data["antigen"]["totals"]["removed_water_residues"],
        data["antigen"]["totals"]["removed_hetero_residues"],
        data["antigen"]["totals"]["removed_nonstandard_aa_residues"],
    ]

    fig = make_subplots(
        rows=1,
        cols=2,
        specs=[[{"type": "domain"}, {"type": "domain"}]],
        subplot_titles=("Antibody Removal Composition", "Antigen Removal Composition"),
    )
    fig.add_trace(
        go.Pie(labels=labels, values=ab_values, textinfo="label+percent+value", marker=dict(colors=["#0ea5e9", "#f97316", "#ef4444"])),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Pie(labels=labels, values=ag_values, textinfo="label+percent+value", marker=dict(colors=["#22d3ee", "#fb923c", "#f87171"])),
        row=1,
        col=2,
    )
    fig.update_layout(title="Residue Filtering Breakdown", margin=dict(l=20, r=20, t=60, b=20))
    return fig


def chain_table_figure(data):
    rows = []
    for entity in ("antibody", "antigen"):
        role = data[entity]["role"]
        for row in data[entity]["chain_summaries"]:
            rows.append(
                {
                    "entity": entity,
                    "role": role,
                    **row,
                }
            )

    rows.sort(key=lambda x: (x["entity"], x["chain_id"]))
    fig = go.Figure(
        data=[
            go.Table(
                header=dict(
                    values=[
                        "Entity",
                        "Role",
                        "Chain",
                        "Input Residues",
                        "Kept Residues",
                        "Removed Water",
                        "Removed Hetero",
                        "Removed Nonstandard AA",
                    ],
                    fill_color="#1e293b",
                    font=dict(color="white", size=12),
                    align="left",
                ),
                cells=dict(
                    values=[
                        [r["entity"] for r in rows],
                        [r["role"] for r in rows],
                        [r["chain_id"] for r in rows],
                        [r["input_residues"] for r in rows],
                        [r["kept_residues"] for r in rows],
                        [r["removed_water_residues"] for r in rows],
                        [r["removed_hetero_residues"] for r in rows],
                        [r["removed_nonstandard_aa_residues"] for r in rows],
                    ],
                    align="left",
                    fill_color="#f8fafc",
                ),
            )
        ]
    )
    fig.update_layout(title="Per-Chain Preparation Summary", margin=dict(l=20, r=20, t=60, b=20), height=450)
    return fig


def main():
    parser = argparse.ArgumentParser(description="Generate interactive report for structure preparation (Node 2).")
    parser.add_argument("--data_json", default="2_structure_prep_outputs/data.json")
    parser.add_argument("--output_html", default="2_structure_prep_outputs/report.html")
    args = parser.parse_args()

    with open(args.data_json, encoding="utf-8") as f:
        data = json.load(f)

    fig_before_after = before_after_figure(data)
    fig_removal = removal_figure(data)
    fig_table = chain_table_figure(data)

    html = f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"UTF-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\" />
  <title>Node 2 - Structure Preparation Report</title>
  <style>
    body {{
      margin: 0;
      font-family: "Segoe UI", Tahoma, sans-serif;
      color: #0f172a;
      background: linear-gradient(145deg, #f8fafc 0%, #e2e8f0 100%);
    }}
    .container {{
      max-width: 1200px;
      margin: 0 auto;
      padding: 24px;
    }}
    .card {{
      background: #ffffff;
      border-radius: 14px;
      box-shadow: 0 12px 24px rgba(15, 23, 42, 0.08);
      margin-bottom: 18px;
      padding: 18px 20px;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
      margin-top: 10px;
    }}
    .metric {{
      border: 1px solid #e2e8f0;
      border-radius: 10px;
      padding: 12px;
      background: #f8fafc;
    }}
    .muted {{ color: #475569; }}
    a {{ color: #1d4ed8; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
  </style>
</head>
<body>
  <div class=\"container\">
    <div class=\"card\">
      <h1>Node 2 - DiffDock-PP Structure Preparation</h1>
      <p class=\"muted\">Standard amino acids retained, water/hetero removed, and residues renumbered from 1 per chain.</p>
      <div class=\"grid\">
        <div class=\"metric\"><strong>Antibody Role</strong><br><span class=\"muted\">{data['antibody']['role']} (receptor)</span></div>
        <div class=\"metric\"><strong>Antigen Role</strong><br><span class=\"muted\">{data['antigen']['role']} (ligand)</span></div>
        <div class=\"metric\"><strong>Geometric Changes</strong><br><span class=\"muted\">{data.get('geometric_modification', False)}</span></div>
        <div class=\"metric\"><strong>Status</strong><br><span class=\"muted\">{data.get('status', 'unknown')}</span></div>
      </div>
    </div>

    <div class=\"card\">{fig_before_after.to_html(full_html=False, include_plotlyjs='cdn')}</div>
    <div class=\"card\">{fig_removal.to_html(full_html=False, include_plotlyjs=False)}</div>
    <div class=\"card\">{fig_table.to_html(full_html=False, include_plotlyjs=False)}</div>

    <div class=\"card\">
      <h2>Processed Outputs</h2>
      <ul>
        <li><a href=\"processed_antibody.pdb\">processed_antibody.pdb</a></li>
        <li><a href=\"processed_antigen.pdb\">processed_antigen.pdb</a></li>
        <li><a href=\"data.json\">data.json</a></li>
      </ul>
    </div>
  </div>
</body>
</html>
"""

    os.makedirs(os.path.dirname(args.output_html) or ".", exist_ok=True)
    with open(args.output_html, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"HTML written to {args.output_html}")


if __name__ == "__main__":
    main()