#!/usr/bin/env python3
"""
Example: Rebuild the workflow-025 Vina vs GNINA report using chiral_report.

Run:  python lib/chiral_report/example_025.py
Output: /tmp/chiral_report_demo.html
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from chiral_report import Report, RawHTML
from chiral_report import components as C
from chiral_report import theme

DEMO_DATA = Path(__file__).resolve().parent / "demo_data"

# ── Tool colors (consistent throughout the report) ──
VINA_COLOR = "#0072B2"   # Okabe-Ito blue
GNINA_COLOR = "#D55E00"  # Okabe-Ito vermillion

# ── Data (extracted from workflow-025 sample_output) ──

COMPOUNDS = [
    "4-chlorobenzenesulfonamide",
    "4-fluorobenzenesulfonamide",
    "4-methylbenzenesulfonamide",
    "acetazolamide",
    "benzene",
    "benzenesulfonamide",
    "ethoxzolamide",
    "ibuprofen",
    "indole",
    "naphthalene-2-sulfonamide",
]

VINA_DG = [-6.609, -6.699, -6.588, -6.512, -4.083, -6.465, -6.693, -6.439, -5.527, -7.829]
GNINA_CNN = [0.969, 0.970, 0.972, 0.956, 0.751, 0.965, 0.946, 0.945, 0.745, 0.976]
GNINA_AFF = [6.524, 6.478, 6.669, 6.900, 3.217, 6.343, 6.933, 5.757, 3.534, 6.517]
VINA_RANKS = [4, 2, 5, 6, 10, 7, 3, 8, 9, 1]
GNINA_RANKS = [4, 3, 2, 6, 9, 5, 7, 8, 10, 1]

RMSD_COMPOUNDS = [
    "4-chlorobenzenesulfonamide", "4-fluorobenzenesulfonamide",
    "4-methylbenzenesulfonamide", "indole", "benzenesulfonamide",
    "ethoxzolamide", "naphthalene-2-sulfonamide", "acetazolamide",
    "ibuprofen", "benzene",
]
RMSD_VALUES = [2.902, 3.038, 3.066, 3.126, 3.323, 4.014, 4.139, 4.718, 4.990, 9.392]

# ── Build report ──

report = Report(
    "Vina vs GNINA: Docking Comparison",
    subtitle="Target: Carbonic Anhydrase II (CA-II), PDB: 1OKL  |  10 Compounds  |  CPU-only",
)

# Stat cards
report.stat_cards([
    {"value": "10", "label": "Compounds Screened"},
    {"value": "0.806", "label": "Spearman ρ (ranks)"},
    {"value": "0.644", "label": "Kendall τ (ranks)"},
    {"value": "4 / 5", "label": "Top-5 Overlap"},
])

# Score correlation scatter
report.chart(
    "Score Correlation: Vina ΔG vs GNINA CNNscore",
    [{
        "div_id": "corr-chart",
        "traces": [{
            "type": "scatter", "mode": "markers",
            "x": VINA_DG, "y": GNINA_CNN, "text": COMPOUNDS,
            "marker": {
                "size": 12,
                "color": VINA_DG,
                "colorscale": "RdBu", "reversescale": True,
                "showscale": True,
                "colorbar": {"title": "Vina ΔG"},
            },
            "hovertemplate": "<b>%{text}</b><br>Vina ΔG: %{x:.2f} kcal/mol<br>GNINA CNN: %{y:.3f}<extra></extra>",
        }],
        "layout": {
            "xaxis": {"title": "Vina ΔG (kcal/mol)"},
            "yaxis": {"title": "GNINA CNNscore"},
            "height": 380,
        },
    }],
    note="Each point is one compound. Color encodes Vina ΔG (blue = stronger binding). "
         "Spearman ρ = 0.806 indicates strong rank agreement.",
)

# Rank agreement scatter
max_rank = 10
report.chart(
    "Ranking Agreement: Vina Rank vs GNINA Rank",
    [{
        "div_id": "rank-chart",
        "traces": [
            {
                "type": "scatter", "mode": "lines",
                "x": [1, max_rank], "y": [1, max_rank],
                "line": {"color": "#cbd5e1", "dash": "dot", "width": 1.5},
                "hoverinfo": "none", "showlegend": False,
            },
            {
                "type": "scatter", "mode": "markers",
                "x": VINA_RANKS, "y": GNINA_RANKS, "text": COMPOUNDS,
                "marker": {"size": 12, "color": VINA_COLOR},
                "hovertemplate": "<b>%{text}</b><br>Vina rank: %{x}<br>GNINA rank: %{y}<extra></extra>",
                "name": "Compounds",
            },
        ],
        "layout": {
            "xaxis": {"title": "Vina Rank (1 = best)", "dtick": 1},
            "yaxis": {"title": "GNINA Rank (1 = best)", "dtick": 1},
            "height": 380,
        },
    }],
    note="Points on the diagonal = perfect rank agreement. Kendall τ = 0.644.",
)

# Score distributions (two charts in one section via raw HTML)
report.chart(
    "Score Distributions",
    [
        {
            "div_id": "vina-hist",
            "traces": [{
                "type": "histogram", "x": VINA_DG,
                "marker": {"color": VINA_COLOR, "opacity": 0.8},
                "xbins": {"size": 0.5}, "name": "Vina ΔG",
            }],
            "layout": {
                "xaxis": {"title": "Vina ΔG (kcal/mol)"},
                "yaxis": {"title": "Count"},
                "height": 300, "showlegend": False,
            },
        },
        {
            "div_id": "gnina-hist",
            "traces": [{
                "type": "histogram", "x": GNINA_CNN,
                "marker": {"color": GNINA_COLOR, "opacity": 0.8},
                "xbins": {"size": 0.05}, "name": "GNINA CNNscore",
            }],
            "layout": {
                "xaxis": {"title": "GNINA CNNscore (0–1)", "range": [0, 1]},
                "yaxis": {"title": "Count"},
                "height": 300, "showlegend": False,
            },
        },
    ],
    height=300,
)

# Affinity correlation
report.chart(
    "Affinity Correlation: Vina ΔG vs GNINA CNN Affinity",
    [{
        "div_id": "aff-chart",
        "traces": [{
            "type": "scatter", "mode": "markers",
            "x": VINA_DG, "y": GNINA_AFF, "text": COMPOUNDS,
            "marker": {"size": 12, "color": "#009E73"},
            "hovertemplate": "<b>%{text}</b><br>Vina: %{x:.2f}<br>GNINA CNN aff: %{y:.2f}<extra></extra>",
        }],
        "layout": {
            "xaxis": {"title": "Vina ΔG (kcal/mol)"},
            "yaxis": {"title": "GNINA CNN Affinity (kcal/mol)"},
            "height": 340,
        },
    }],
    note="Spearman ρ = −0.661. Note: GNINA CNN affinity uses the opposite sign convention (positive = stronger).",
)

# Full comparison table
all_rows = []
for i in sorted(range(10), key=lambda i: VINA_RANKS[i]):
    rd = GNINA_RANKS[i] - VINA_RANKS[i]
    all_rows.append([
        COMPOUNDS[i], VINA_RANKS[i], f"{VINA_DG[i]:.3f}",
        GNINA_RANKS[i], f"{GNINA_CNN[i]:.3f}", f"{GNINA_AFF[i]:.2f}",
        f"+{rd}" if rd > 0 else str(rd),
    ])

report.table(
    "Full Comparison Table",
    ["Name", "Vina Rank", "Vina ΔG", "GNINA Rank", "CNN Score", "CNN Aff.", "ΔRank"],
    all_rows,
    sortable=True,
)

# 3D Pose Visualization (Mol*)
with open(DEMO_DATA / "manifest.json") as f:
    manifest = json.load(f)

structures = []
for entry in manifest:
    data_path = DEMO_DATA / entry["file"]
    structures.append({
        "label": entry["label"],
        "data": data_path.read_text(),
        "format": entry["format"],
        "color": entry["color"],
    })

report.section(
    "3D Pose Visualization (Mol*)",
    "Toggle structures to compare Vina (blue) and GNINA (vermillion) top-3 poses against the crystal ligand (amber).",
    molstar_viewer=structures,
    note="Receptor as cartoon. Ligand poses as ball-and-stick. "
         "Scroll to zoom. Drag to rotate. Right-click to pan.",
)

# Pose RMSD bar chart
rmsd_colors = [
    theme.SUCCESS if v < 2 else (theme.WARNING if v < 3 else theme.DANGER)
    for v in RMSD_VALUES
]

rmsd_stat_cards = C.stat_cards([
    {"value": "10", "label": "Compounds with both poses"},
    {"value": "0 / 10", "label": "RMSD < 2 Å (same binding mode)"},
    {"value": "0%", "label": "Agreement rate (2 Å cutoff)", "color": "warning"},
    {"value": "4.27", "label": "Mean RMSD (Å)"},
])

rmsd_chart_html = C.plotly_chart("rmsd-chart", height=320)
rmsd_js = C.plotly_script([{
    "div_id": "rmsd-chart",
    "traces": [
        {
            "type": "bar",
            "x": RMSD_COMPOUNDS, "y": RMSD_VALUES,
            "marker": {"color": rmsd_colors},
            "hovertemplate": "<b>%{x}</b><br>RMSD: %{y:.3f} Å<extra></extra>",
            "name": "Vina vs GNINA RMSD",
        },
        {
            "type": "scatter", "mode": "lines",
            "x": RMSD_COMPOUNDS, "y": [2.0] * 10,
            "line": {"color": "#64748b", "dash": "dash", "width": 1.5},
            "hoverinfo": "none", "showlegend": True, "name": "2 Å threshold",
        },
    ],
    "layout": {
        "xaxis": {"title": "Compound", "tickangle": -40, "automargin": True},
        "yaxis": {"title": "RMSD (Å)", "range": [0, 11]},
        "height": 320,
        "legend": {"orientation": "h", "y": 1.08},
    },
}])

report.raw(C.section(
    "Pose Agreement: Vina vs GNINA Top-Pose RMSD per Compound",
    f"{rmsd_stat_cards}\n{rmsd_chart_html}",
    note="RMSD between Vina and GNINA top-ranked pose for each compound (heavy atoms only). "
         "Green < 2 Å: tools agree on binding mode. Amber 2–3 Å: moderate deviation. "
         "Red > 3 Å: divergent poses.",
))
report._extra_js.append(rmsd_js)

# Runtime comparison
runtime_chart_html = C.plotly_chart("runtime-chart", height=220)
runtime_js = C.plotly_script([{
    "div_id": "runtime-chart",
    "traces": [{
        "type": "bar", "orientation": "h",
        "x": [12.0, 736.6], "y": ["Vina", "GNINA"],
        "text": ["12.0 s", "12.3 min (737 s)"],
        "textposition": "outside",
        "marker": {"color": [VINA_COLOR, GNINA_COLOR]},
        "hovertemplate": "%{y}: %{x:.1f} s<extra></extra>",
    }],
    "layout": {
        "xaxis": {"title": "Wall-clock time (seconds)", "zeroline": True},
        "yaxis": {"automargin": True},
        "height": 220,
        "showlegend": False,
    },
}])

runtime_table = C.data_table(
    ["", "Vina", "GNINA"],
    [
        ["Total wall-clock", "12.0 s", "12.3 min (737 s)"],
        ["Per-ligand avg", "1.2 s", "1.2 min (74 s)"],
        ["Compounds", "10", "10"],
    ],
)

runtime_body = C.side_by_side(runtime_chart_html, runtime_table)
report.raw(C.section(
    "Runtime Comparison: Vina vs GNINA",
    runtime_body,
    note="Timing covers the docking loop only (excludes library splitting and CSV/JSON writing). "
         "GNINA includes CNN rescoring overhead; both tools ran CPU-only.",
))
report._extra_js.append(runtime_js)

# Methods
report.methods([
    {
        "title": "AutoDock Vina",
        "badge": "Vina", "badge_color": VINA_COLOR,
        "rows": [
            ("Scoring function", "Empirical energy function (Vina weights)"),
            ("Exhaustiveness", "8"),
            ("Num. output modes", "9"),
            ("Hardware", "CPU-only"),
            ("Output format", "PDBQT"),
        ],
    },
    {
        "title": "GNINA",
        "badge": "GNINA", "badge_color": GNINA_COLOR,
        "rows": [
            ("CNN scoring mode", "<code>--cnn_scoring rescore</code>"),
            ("Exhaustiveness", "8"),
            ("Num. output modes", "9"),
            ("Hardware", "CPU-only (<code>--no_gpu</code>)"),
            ("Output format", "SDF"),
        ],
    },
    {
        "title": "Shared Experimental Settings",
        "rows": [
            ("Target protein", "Carbonic Anhydrase II (CA-II), PDB: 1OKL"),
            ("Compounds screened", "10"),
            ("Pose format (downstream)", "Top-5 poses per tool converted to MOL2 via OpenBabel"),
            ("Ranking metric — Vina", "ΔG (kcal/mol); lower (more negative) = stronger predicted binding"),
            ("Ranking metric — GNINA", "CNNscore (0–1); higher = CNN-predicted binder"),
        ],
    },
])

# Save
out = report.save("/tmp/chiral_report_demo.html")
print(f"Report saved to {out}")
print(f"Open in browser:  file://{out}")
