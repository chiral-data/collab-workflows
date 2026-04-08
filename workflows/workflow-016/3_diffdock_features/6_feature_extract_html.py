#!/usr/bin/env python3
import argparse
import json
import os

import plotly.graph_objects as go
from plotly.subplots import make_subplots


def sequence_overview_figure(data):
        labels = ["Antibody", "Antigen"]
        lengths = [data["antibody"]["total_length"], data["antigen"]["total_length"]]
        unknowns = [data["antibody"]["unknown_residue_count"], data["antigen"]["unknown_residue_count"]]

        fig = make_subplots(rows=1, cols=2, subplot_titles=("Total Sequence Length", "Unknown Residues (X)"))
        fig.add_trace(go.Bar(x=labels, y=lengths, marker_color=["#1d4ed8", "#ea580c"], name="Length"), row=1, col=1)
        fig.add_trace(go.Bar(x=labels, y=unknowns, marker_color=["#60a5fa", "#fdba74"], name="Unknown"), row=1, col=2)
        fig.update_layout(title="Node 3 Sequence Quality Overview", margin=dict(l=20, r=20, t=60, b=40), showlegend=False)
        return fig


def chain_length_figure(data):
        ab_rows = data["antibody"]["chain_info"]
        ag_rows = data["antigen"]["chain_info"]
        labels = [f"Ab:{r['chain_id']}" for r in ab_rows] + [f"Ag:{r['chain_id']}" for r in ag_rows]
        values = [r["length"] for r in ab_rows] + [r["length"] for r in ag_rows]
        unknowns = [r["unknown_residue_count"] for r in ab_rows] + [r["unknown_residue_count"] for r in ag_rows]
        colors = ["#2563eb"] * len(ab_rows) + ["#f97316"] * len(ag_rows)

        fig = go.Figure()
        fig.add_trace(
                go.Bar(
                        x=labels,
                        y=values,
                        marker_color=colors,
                        name="Length",
                        hovertemplate="%{x}<br>Length: %{y}<extra></extra>",
                        visible=True,
                )
        )
        fig.add_trace(
                go.Bar(
                        x=labels,
                        y=unknowns,
                        marker_color=colors,
                        name="Unknown X",
                        hovertemplate="%{x}<br>Unknown residues: %{y}<extra></extra>",
                        visible=False,
                )
        )
        fig.update_layout(
                title="Per-Chain Sequence Metrics",
                xaxis_title="Chain",
                yaxis_title="Count",
                margin=dict(l=40, r=20, t=60, b=40),
                updatemenus=[
                        dict(
                                type="buttons",
                                direction="right",
                                x=0,
                                y=1.2,
                                buttons=[
                                        dict(label="Length", method="update", args=[{"visible": [True, False]}, {"yaxis": {"title": "Length"}}]),
                                        dict(label="Unknown X", method="update", args=[{"visible": [False, True]}, {"yaxis": {"title": "Unknown Residues"}}]),
                                ],
                        )
                ],
        )
        return fig


def chain_table_figure(data):
        rows = []
        for entity in ("antibody", "antigen"):
                for row in data[entity]["chain_info"]:
                        rows.append({"entity": entity, **row})

        rows.sort(key=lambda x: (x["entity"], x["chain_id"]))
        fig = go.Figure(
                data=[
                        go.Table(
                                header=dict(
                                        values=["Entity", "Chain", "Length", "Unknown Residues", "Known Fraction"],
                                        fill_color="#1f2937",
                                        font=dict(color="white"),
                                        align="left",
                                ),
                                cells=dict(
                                        values=[
                                                [r["entity"] for r in rows],
                                                [r["chain_id"] for r in rows],
                                                [r["length"] for r in rows],
                                                [r["unknown_residue_count"] for r in rows],
                                                [r["known_residue_fraction"] for r in rows],
                                        ],
                                        align="left",
                                        fill_color="#f8fafc",
                                ),
                        )
                ]
        )
        fig.update_layout(title="Chain-Level Sequence Summary", margin=dict(l=20, r=20, t=60, b=20), height=420)
        return fig


def main():
        parser = argparse.ArgumentParser(description="Generate interactive Node 3 report.")
        parser.add_argument("--data_json", default="3_feature_extract_outputs/data.json")
        parser.add_argument("--output_html", default="3_feature_extract_outputs/report.html")
        args = parser.parse_args()

        with open(args.data_json, encoding="utf-8") as f:
                data = json.load(f)

        ab_chains = ", ".join(data["antibody"]["chains"])
        ag_chains = ", ".join(data["antigen"]["chains"])

        fig_overview = sequence_overview_figure(data)
        fig_chain = chain_length_figure(data)
        fig_table = chain_table_figure(data)

        html = f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
    <meta charset=\"UTF-8\" />
    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\" />
    <title>Node 3 - Feature Extraction Report</title>
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
            box-shadow: 0 10px 24px rgba(15, 23, 42, 0.08);
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
        pre {{
            background: #0f172a;
            color: #e2e8f0;
            padding: 14px;
            border-radius: 10px;
            overflow-x: auto;
        }}
        a {{ color: #1d4ed8; text-decoration: none; }}
        a:hover {{ text-decoration: underline; }}
        .muted {{ color: #475569; }}
    </style>
</head>
<body>
    <div class=\"container\">
        <div class=\"card\">
            <h1>Node 3 - Feature Extraction and ESM-2 Preparation</h1>
            <p class=\"muted\">FASTA generation is complete; embedding files are placeholders until ESM-2 extraction is run.</p>
            <div class=\"grid\">
                <div class=\"metric\"><strong>Antibody Chains</strong><br><span class=\"muted\">{ab_chains}</span></div>
                <div class=\"metric\"><strong>Antigen Chains</strong><br><span class=\"muted\">{ag_chains}</span></div>
                <div class=\"metric\"><strong>Antibody Total Length</strong><br><span class=\"muted\">{data['antibody']['total_length']}</span></div>
                <div class=\"metric\"><strong>Antigen Total Length</strong><br><span class=\"muted\">{data['antigen']['total_length']}</span></div>
            </div>
        </div>

        <div class=\"card\">{fig_overview.to_html(full_html=False, include_plotlyjs='cdn')}</div>
        <div class=\"card\">{fig_chain.to_html(full_html=False, include_plotlyjs=False)}</div>
        <div class=\"card\">{fig_table.to_html(full_html=False, include_plotlyjs=False)}</div>

        <div class=\"card\">
            <h2>Generate Real ESM-2 Embeddings</h2>
            <p class=\"muted\">Run from the output directory where FASTA files are located:</p>
            <pre>python generate_embeddings.py --model esm2_t36_3B_UR50D --repr_layers 36</pre>
            <p class=\"muted\">Equivalent direct commands:</p>
            <pre>esm-extract esm2_t36_3B_UR50D antibody.fasta antibody_features.pt --repr_layers 36 --include mean per_tok
esm-extract esm2_t36_3B_UR50D antigen.fasta antigen_features.pt --repr_layers 36 --include mean per_tok</pre>
        </div>

        <div class=\"card\">
            <h2>Output Files</h2>
            <ul>
                <li><a href=\"antibody.fasta\">antibody.fasta</a></li>
                <li><a href=\"antigen.fasta\">antigen.fasta</a></li>
                <li><a href=\"antibody_features.pt\">antibody_features.pt</a></li>
                <li><a href=\"antigen_features.pt\">antigen_features.pt</a></li>
                <li><a href=\"generate_embeddings.py\">generate_embeddings.py</a></li>
                <li><a href=\"sequence_info.json\">sequence_info.json</a></li>
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