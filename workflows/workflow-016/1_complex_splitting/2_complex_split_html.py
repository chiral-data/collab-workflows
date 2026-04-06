#!/usr/bin/env python3
import argparse
import json
import os

import plotly.graph_objects as go
from plotly.subplots import make_subplots


def pie_composition_figure(data):
    ab_res = data["antibody"]["residue_count"]
    ag_res = data["antigen"]["residue_count"]
    ab_atoms = data["antibody"]["atom_count"]
    ag_atoms = data["antigen"]["atom_count"]
    ab_chain_n = data["antibody"]["chain_count"]
    ag_chain_n = data["antigen"]["chain_count"]

    fig = make_subplots(
        rows=1,
        cols=3,
        subplot_titles=("Residue Composition", "Atom Composition", "Chain Composition"),
        specs=[[{"type": "domain"}, {"type": "domain"}, {"type": "domain"}]],
    )
    fig.add_trace(
        go.Pie(
            labels=["Antibody", "Antigen"],
            values=[ab_res, ag_res],
            marker=dict(colors=["#1f77b4", "#ff7f0e"]),
            textinfo="label+percent+value",
            hovertemplate="%{label}<br>Residues: %{value}<extra></extra>",
            name="Residues",
        ),
        row=1,
        col=1,
    )
    fig.add_trace(
        go.Pie(
            labels=["Antibody", "Antigen"],
            values=[ab_atoms, ag_atoms],
            marker=dict(colors=["#2ca02c", "#d62728"]),
            textinfo="label+percent+value",
            hovertemplate="%{label}<br>Atoms: %{value}<extra></extra>",
            name="Atoms",
        ),
        row=1,
        col=2,
    )
    fig.add_trace(
        go.Pie(
            labels=["Antibody", "Antigen"],
            values=[ab_chain_n, ag_chain_n],
            marker=dict(colors=["#9467bd", "#8c564b"]),
            textinfo="label+percent+value",
            hovertemplate="%{label}<br>Chains: %{value}<extra></extra>",
            name="Chains",
        ),
        row=1,
        col=3,
    )
    fig.update_layout(
        title="Complex Composition Overview",
        margin=dict(l=20, r=20, t=60, b=20),
        legend=dict(orientation="h", yanchor="bottom", y=-0.2, x=0.5, xanchor="center"),
    )
    return fig


def per_chain_bar_figure(data):
    chain_records = data["antibody"]["chain_info"] + data["antigen"]["chain_info"]
    ids = [rec["id"] for rec in chain_records]
    residues = [rec["residues"] for rec in chain_records]
    atoms = [rec["atoms"] for rec in chain_records]
    classes = [rec["type"] for rec in chain_records]

    color_map = {"antibody": "#1f77b4", "antigen": "#ff7f0e"}
    colors = [color_map.get(c, "#7f7f7f") for c in classes]

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=ids,
            y=residues,
            marker_color=colors,
            name="Residues",
            hovertemplate="Chain %{x}<br>Residues: %{y}<extra></extra>",
            visible=True,
        )
    )
    fig.add_trace(
        go.Bar(
            x=ids,
            y=atoms,
            marker_color=colors,
            name="Atoms",
            hovertemplate="Chain %{x}<br>Atoms: %{y}<extra></extra>",
            visible=False,
        )
    )

    fig.update_layout(
        title="Per-Chain Statistics",
        xaxis_title="Chain ID",
        yaxis_title="Count",
        barmode="group",
        margin=dict(l=50, r=20, t=60, b=50),
        updatemenus=[
            dict(
                type="buttons",
                direction="right",
                x=0.0,
                y=1.2,
                buttons=[
                    dict(
                        label="Residues",
                        method="update",
                        args=[{"visible": [True, False]}, {"yaxis": {"title": "Residue Count"}}],
                    ),
                    dict(
                        label="Atoms",
                        method="update",
                        args=[{"visible": [False, True]}, {"yaxis": {"title": "Atom Count"}}],
                    ),
                ],
            )
        ],
    )
    return fig


def chain_table_figure(data):
    chain_records = data["antibody"]["chain_info"] + data["antigen"]["chain_info"]
    rows = sorted(chain_records, key=lambda x: x["id"])
    header = ["Chain", "Type", "Residues", "Atoms", "Hetero Residues"]
    fig = go.Figure(
        data=[
            go.Table(
                header=dict(values=header, fill_color="#1f2937", font=dict(color="white", size=12), align="left"),
                cells=dict(
                    values=[
                        [r["id"] for r in rows],
                        [r["type"] for r in rows],
                        [r["residues"] for r in rows],
                        [r["atoms"] for r in rows],
                        [r.get("hetero_residues", 0) for r in rows],
                    ],
                    fill_color="#f9fafb",
                    align="left",
                ),
            )
        ]
    )
    fig.update_layout(title="Chain-Level Scientific Summary", margin=dict(l=20, r=20, t=60, b=20), height=420)
    return fig


def main():
    parser = argparse.ArgumentParser(description="Generate interactive HTML report for Node 1 complex splitting.")
    parser.add_argument("--data_json", default="1_complex_split_outputs/data.json")
    parser.add_argument("--output_html", default="1_complex_split_outputs/report.html")
    args = parser.parse_args()

    with open(args.data_json, encoding="utf-8") as f:
        data = json.load(f)

    ab_chains = ", ".join(data["antibody"]["chains"])
    ag_chains = ", ".join(data["antigen"]["chains"])

    composition_fig = pie_composition_figure(data)
    chain_bar_fig = per_chain_bar_figure(data)
    chain_table_fig = chain_table_figure(data)

    html = f"""<!DOCTYPE html>
<html lang=\"en\">
<head>
  <meta charset=\"UTF-8\" />
  <meta name=\"viewport\" content=\"width=device-width, initial-scale=1.0\" />
  <title>Node 1 - Complex Splitting Report</title>
  <style>
    body {{
      margin: 0;
      padding: 0;
      font-family: "Segoe UI", Tahoma, sans-serif;
      color: #111827;
      background: linear-gradient(135deg, #f8fafc, #eef2ff);
    }}
    .container {{
      max-width: 1200px;
      margin: 0 auto;
      padding: 24px;
    }}
    .card {{
      background: #ffffff;
      border-radius: 14px;
      box-shadow: 0 10px 25px rgba(15, 23, 42, 0.08);
      padding: 18px 22px;
      margin-bottom: 18px;
    }}
    .grid {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
      gap: 12px;
      margin-top: 10px;
    }}
    .metric {{
      background: #f8fafc;
      border: 1px solid #e5e7eb;
      border-radius: 10px;
      padding: 12px;
    }}
    h1, h2 {{ margin: 0 0 8px 0; }}
    .muted {{ color: #4b5563; }}
    a {{ color: #1d4ed8; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    code {{ background: #f3f4f6; padding: 2px 6px; border-radius: 6px; }}
        .viewer-wrap {{
            border: 1px solid #e5e7eb;
            border-radius: 12px;
            overflow: hidden;
            background: #0b1220;
        }}
        .viewer-toolbar {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
            gap: 10px;
            padding: 12px;
            background: #f8fafc;
            border-bottom: 1px solid #e5e7eb;
        }}
        .viewer-toolbar label {{
            font-size: 12px;
            color: #374151;
            display: block;
            margin-bottom: 4px;
            font-weight: 600;
        }}
        .viewer-toolbar select,
        .viewer-toolbar input,
        .viewer-toolbar button {{
            width: 100%;
            padding: 8px 10px;
            border-radius: 8px;
            border: 1px solid #d1d5db;
            font-size: 13px;
            background: #fff;
        }}
        .viewer-toolbar button {{
            cursor: pointer;
            background: #1d4ed8;
            color: white;
            border: none;
            font-weight: 600;
        }}
        .viewer-toolbar button.secondary {{
            background: #4b5563;
        }}
        #proteinViewer {{
            width: 100%;
            height: 560px;
            min-height: 420px;
        }}
        .small {{ font-size: 12px; color: #6b7280; }}
  </style>
</head>
<body>
  <div class=\"container\">
    <div class=\"card\">
      <h1>Node 1 - Antibody/Antigen Complex Splitting</h1>
      <p class=\"muted\">Interactive and science-focused summary of chain assignment and structural composition.</p>
      <div class=\"grid\">
        <div class=\"metric\"><strong>Input PDB</strong><br><span class=\"muted\">{data['input_file']}</span></div>
        <div class=\"metric\"><strong>Detection Method</strong><br><span class=\"muted\">{data.get('detection_method', 'n/a')}</span></div>
        <div class=\"metric\"><strong>Antibody Chains</strong><br><span class=\"muted\">{ab_chains}</span></div>
        <div class=\"metric\"><strong>Antigen Chains</strong><br><span class=\"muted\">{ag_chains}</span></div>
      </div>
    </div>

    <div class=\"card\">{composition_fig.to_html(full_html=False, include_plotlyjs='cdn')}</div>
    <div class=\"card\">{chain_bar_fig.to_html(full_html=False, include_plotlyjs=False)}</div>
    <div class=\"card\">{chain_table_fig.to_html(full_html=False, include_plotlyjs=False)}</div>

        <div class=\"card\">
            <h2>Interactive 3D Protein Viewer (iCn3D)</h2>
            <p class=\"muted\">Reusable viewer API for generated HTML reports. Loads local PDB files from this output folder.</p>
            <div class=\"viewer-wrap\">
                <div class=\"viewer-toolbar\">
                    <div>
                        <label for=\"proteinFile\">Protein File</label>
                        <select id=\"proteinFile\">
                            <option value=\"antibody.pdb\">antibody.pdb</option>
                            <option value=\"antigen.pdb\">antigen.pdb</option>
                            <option value=\"original_complex.pdb\">original_complex.pdb</option>
                        </select>
                    </div>
                    <div>
                        <label for=\"renderStyle\">Style</label>
                        <select id=\"renderStyle\">
                            <option value=\"cartoon\" selected>cartoon</option>
                            <option value=\"surface\">surface</option>
                            <option value=\"stick\">stick</option>
                        </select>
                    </div>
                    <div>
                        <label for=\"colorScheme\">Color Scheme</label>
                        <select id=\"colorScheme\">
                            <option value=\"chain\" selected>chain</option>
                            <option value=\"spectrum\">spectrum</option>
                            <option value=\"secondary\">secondary</option>
                            <option value=\"atom\">atom</option>
                        </select>
                    </div>
                    <div>
                        <label for=\"highlightResidues\">Highlight Residues (comma-separated)</label>
                        <input id=\"highlightResidues\" placeholder=\"H:33,H:52,L:91\" />
                    </div>
                    <div>
                        <label>&nbsp;</label>
                        <button id=\"btnLoad\" type=\"button\">Load / Update</button>
                    </div>
                    <div>
                        <label>&nbsp;</label>
                        <button id=\"btnShot\" class=\"secondary\" type=\"button\">Export Screenshot</button>
                    </div>
                </div>
                <div id=\"proteinViewer\"></div>
            </div>
            <p class=\"small\">Uses iCn3D JS embedding API. If CDN is blocked, place local copies of <code>icn3d.min.js</code>/<code>icn3d.css</code> and dependencies in this folder.</p>
        </div>

    <div class=\"card\">
      <h2>Generated Files</h2>
      <ul>
        <li><a href=\"antibody.pdb\">antibody.pdb</a></li>
        <li><a href=\"antigen.pdb\">antigen.pdb</a></li>
        <li><a href=\"original_complex.pdb\">original_complex.pdb</a></li>
        <li><a href=\"chain_info.json\">chain_info.json</a></li>
        <li><a href=\"data.json\">data.json</a></li>
      </ul>
    </div>
  </div>

    <link rel=\"stylesheet\" href=\"https://www.ncbi.nlm.nih.gov/Structure/icn3d/lib/jquery-ui.min.css\">
    <link rel=\"stylesheet\" href=\"https://www.ncbi.nlm.nih.gov/Structure/icn3d/icn3d.css\">
    <script src=\"https://www.ncbi.nlm.nih.gov/Structure/icn3d/lib/jquery.min.js\"></script>
    <script src=\"https://www.ncbi.nlm.nih.gov/Structure/icn3d/lib/jquery-ui.min.js\"></script>
    <script src=\"https://www.ncbi.nlm.nih.gov/Structure/icn3d/lib/threeClass.min.js\"></script>
    <script src=\"https://www.ncbi.nlm.nih.gov/Structure/icn3d/icn3d.min.js\"></script>

    <script>
        // Reusable API for your HTML generator pipeline.
        // proteinFile: local pdb filename/path, width/height: CSS sizes, style: cartoon|surface|stick
        let icn3dUiInstance = null;

        async function loadProtein(proteinFile, width = '100%', height = '560px', style = 'cartoon', opts = {{}}) {{
            const container = document.getElementById('proteinViewer');
            container.style.width = width;
            container.style.height = height;

            const response = await fetch(proteinFile);
            if (!response.ok) {{
                throw new Error(`Failed to load ${{proteinFile}} (${{response.status}})`);
            }}
            const pdbStr = await response.text();

            const cfg = {{
                divid: 'proteinViewer',
                width,
                height,
                resize: true,
                mobilemenu: true,
                showcommand: false,
                showtitle: false,
            }};

            // Fresh init per load keeps it deterministic in generated pages.
            window.icn3dui = new icn3d.iCn3DUI(cfg);
            icn3dUiInstance = window.icn3dui;
            await icn3dUiInstance.show3DStructure(pdbStr);

            applyRenderStyle(style);
            if (opts.colorScheme) setColorScheme(opts.colorScheme);
            if (opts.highlightResidues) highlightResidues(opts.highlightResidues);

            return icn3dUiInstance;
        }}

        function applyRenderStyle(style) {{
            if (!icn3dUiInstance) return;
            const ic = icn3dUiInstance.icn3d;
            try {{
                const cmdMap = {{
                    cartoon: 'style proteins cartoon',
                    surface: 'surface proteins',
                    stick: 'style proteins stick',
                }};
                const cmd = cmdMap[style] || cmdMap.cartoon;
                ic.applyCommandCls.applyCommand(cmd);
            }} catch (e) {{
                // Fallback: keep default rendering if command syntax changes.
                console.warn('Style command failed, using default style', e);
            }}
        }}

        function setColorScheme(scheme) {{
            if (!icn3dUiInstance) return;
            try {{
                const ic = icn3dUiInstance.icn3d;
                const cmdMap = {{
                    chain: 'color chain',
                    spectrum: 'color spectrum',
                    secondary: 'color secondary structure',
                    atom: 'color atom',
                }};
                ic.applyCommandCls.applyCommand(cmdMap[scheme] || 'color chain');
            }} catch (e) {{
                console.warn('Color command failed', e);
            }}
        }}

        function highlightResidues(residueString) {{
            if (!icn3dUiInstance || !residueString) return;
            const list = residueString.split(',').map(s => s.trim()).filter(Boolean);
            if (!list.length) return;
            try {{
                const ic = icn3dUiInstance.icn3d;
                // Example selection syntax: H:33,H:52,L:91
                ic.applyCommandCls.applyCommand(`select ${{list.join(' or ')}}`);
                ic.applyCommandCls.applyCommand('color red');
                ic.applyCommandCls.applyCommand('style proteins stick');
            }} catch (e) {{
                console.warn('Residue highlight failed', e);
            }}
        }}

        function exportScreenshot(filename = 'protein_view.png') {{
            const canvas = document.querySelector('#proteinViewer canvas');
            if (!canvas) {{
                alert('No rendered canvas found yet. Load a protein first.');
                return;
            }}
            const a = document.createElement('a');
            a.href = canvas.toDataURL('image/png');
            a.download = filename;
            a.click();
        }}

        async function loadFromUi() {{
            const proteinFile = document.getElementById('proteinFile').value;
            const style = document.getElementById('renderStyle').value;
            const colorScheme = document.getElementById('colorScheme').value;
            const residues = document.getElementById('highlightResidues').value;
            try {{
                await loadProtein(proteinFile, '100%', '560px', style, {{
                    colorScheme,
                    highlightResidues: residues,
                }});
            }} catch (err) {{
                console.error(err);
                alert(`Viewer load failed: ${{err.message}}`);
            }}
        }}

        document.addEventListener('DOMContentLoaded', async () => {{
            document.getElementById('btnLoad').addEventListener('click', loadFromUi);
            document.getElementById('btnShot').addEventListener('click', () => exportScreenshot('node1_icn3d.png'));
            await loadFromUi();
        }});
    </script>
</body>
</html>
"""

    os.makedirs(os.path.dirname(args.output_html) or ".", exist_ok=True)
    with open(args.output_html, "w", encoding="utf-8") as f:
        f.write(html)

    print(f"HTML report written to {args.output_html}")


if __name__ == "__main__":
    main()