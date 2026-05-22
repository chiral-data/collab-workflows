#!/usr/bin/env python3
import argparse
import json
import os

import plotly.graph_objects as go
from plotly.subplots import make_subplots


def confidence_figure(data):
        poses = data.get("pose_details", [])
        ranks = [p["rank"] for p in poses]
        confidences = [p["confidence"] for p in poses]

        fig = make_subplots(
                rows=1,
                cols=2,
                subplot_titles=("Confidence by Final Rank", "Confidence Distribution"),
                specs=[[{"type": "xy"}, {"type": "xy"}]],
        )
        fig.add_trace(
                go.Bar(
                        x=ranks,
                        y=confidences,
                        marker_color="#2563eb",
                        name="Confidence",
                        hovertemplate="Rank %{x}<br>Confidence %{y:.4f}<extra></extra>",
                ),
                row=1,
                col=1,
        )
        fig.add_trace(
                go.Scatter(
                        x=ranks,
                        y=confidences,
                        mode="lines+markers",
                        marker=dict(color="#f97316", size=7),
                        line=dict(color="#f97316", width=2),
                        name="Trend",
                        hovertemplate="Rank %{x}<br>Confidence %{y:.4f}<extra></extra>",
                ),
                row=1,
                col=2,
        )
        fig.update_xaxes(title_text="Rank", row=1, col=1)
        fig.update_xaxes(title_text="Rank", row=1, col=2)
        fig.update_yaxes(title_text="Confidence", row=1, col=1)
        fig.update_yaxes(title_text="Confidence", row=1, col=2)
        fig.update_layout(title="DiffDock-PP Pose Confidence", margin=dict(l=20, r=20, t=60, b=30), showlegend=False)
        return fig


def pose_table_figure(data):
        rows = data.get("pose_details", [])
        rows = sorted(rows, key=lambda x: x["rank"])
        fig = go.Figure(
                data=[
                        go.Table(
                                header=dict(
                                        values=["Rank", "Confidence", "Pose File", "Source Name", "Mode"],
                                        fill_color="#1f2937",
                                        font=dict(color="white", size=12),
                                        align="left",
                                ),
                                cells=dict(
                                        values=[
                                                [r["rank"] for r in rows],
                                                [f"{r['confidence']:.4f}" for r in rows],
                                                [r["pose_filename"] for r in rows],
                                                [r.get("source_name", "") for r in rows],
                                                [r.get("mode", "") for r in rows],
                                        ],
                                        align="left",
                                        fill_color="#f8fafc",
                                ),
                        )
                ]
        )
        fig.update_layout(title="Ranked Pose Summary", margin=dict(l=20, r=20, t=60, b=20), height=420)
        return fig


def main():
        parser = argparse.ArgumentParser(description="Generate interactive report for Node 4 inference.")
        parser.add_argument("--data_json", default="4_diffdock_inference_outputs/data.json")
        parser.add_argument("--output_html", default="4_diffdock_inference_outputs/report.html")
        args = parser.parse_args()

        with open(args.data_json, encoding="utf-8") as f:
                data = json.load(f)

        conf_fig = confidence_figure(data)
        table_fig = pose_table_figure(data)
        pose_links = "\n".join(
                [
                        f'<li><a href="{row["pose_filename"]}">{row["pose_filename"]}</a> '
                        f'(confidence={row["confidence"]:.4f})</li>'
                        for row in sorted(data.get("pose_details", []), key=lambda x: x["rank"])
                ]
        )

        html = f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
    <meta charset=\"UTF-8\" />
    <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\" />
    <title>Node 4 - DiffDock-PP Inference Report</title>
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
            background: #f8fafc;
            padding: 12px;
        }}
        .muted {{ color: #475569; }}
        a {{ color: #1d4ed8; text-decoration: none; }}
        a:hover {{ text-decoration: underline; }}
    </style>
</head>
<body>
    <div class=\"container\">
        <div class=\"card\">
            <h1>Node 4 - DiffDock-PP Docking Inference</h1>
            <p class=\"muted\">Predicted antigen poses are ranked by confidence. Highest confidence is exported as rank1.pdb.</p>
            <div class=\"grid\">
                <div class=\"metric\"><strong>Status</strong><br><span class=\"muted\">{data.get('status', 'unknown')}</span></div>
                <div class=\"metric\"><strong>Run Mode</strong><br><span class=\"muted\">{data.get('mode', 'unknown')}</span></div>
                <div class=\"metric\"><strong>Num Poses</strong><br><span class=\"muted\">{data.get('num_poses', 0)}</span></div>
                <div class=\"metric\"><strong>Top Confidence</strong><br><span class=\"muted\">{(data.get('confidence_scores') or [0])[0]:.4f}</span></div>
            </div>
        </div>

        <div class=\"card\">{conf_fig.to_html(full_html=False, include_plotlyjs='cdn')}</div>
        <div class=\"card\">{table_fig.to_html(full_html=False, include_plotlyjs=False)}</div>

        <div class=\"card\">
            <h2>Pose Files</h2>
            <ul>{pose_links}</ul>
            <p><a href=\"confidence_scores.json\">confidence_scores.json</a></p>
            <p><a href=\"inference_log.txt\">inference_log.txt</a></p>
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