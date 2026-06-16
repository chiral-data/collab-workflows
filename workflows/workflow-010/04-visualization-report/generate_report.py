"""
Node 4 – PDB Export & Interactive HTML Visualization
=====================================================
Copies each PDB file from Node 3 into the output directory and
generates a single self-contained HTML report with an interactive
3Dmol viewer, pLDDT confidence coloring, per-residue pLDDT chart,
and structure summary statistics for every predicted structure.

The report works in any modern web browser with no Python runtime
required – all PDB data and pLDDT values are embedded directly in
the HTML.

Outputs
-------
outputs/report.html          – interactive viewer with pLDDT analysis
outputs/<pair_id>.pdb        – copied PDB files (for download links)
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path
from typing import Dict, List, Tuple

# ---------------------------------------------------------------------------
# PDB parsing
# ---------------------------------------------------------------------------

_AA3TO1: Dict[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def parse_pdb_chains(pdb_text: str) -> Dict[str, List[dict]]:
    """
    Parse CA atoms from PDB ATOM records.
    Returns dict: chain_id → list of {resnum, aa, plddt}.
    """
    chains: Dict[str, List[dict]] = {}
    seen: set = set()
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        atom_name = line[12:16].strip()
        if atom_name != "CA":
            continue
        chain = line[21]
        resnum_str = line[22:26].strip()
        resname = line[17:20].strip()
        key = (chain, resnum_str)
        if key in seen:
            continue
        seen.add(key)
        try:
            bfactor = float(line[60:66])
        except (ValueError, IndexError):
            bfactor = 0.0
        chains.setdefault(chain, []).append({
            "resnum": int(resnum_str),
            "aa": _AA3TO1.get(resname, "X"),
            "plddt": round(bfactor, 2),
        })
    return chains


def chain_stats(residues: List[dict]) -> Tuple[int, float, float]:
    """Return (length, avg_plddt, min_plddt)."""
    n = len(residues)
    if n == 0:
        return 0, 0.0, 0.0
    values = [r["plddt"] for r in residues]
    return n, round(sum(values) / n, 1), round(min(values), 1)


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------

_3DMOL_BUNDLED_PATH = Path("/workflow/3Dmol-min.js")
_3DMOL_CDN = "https://3dmol.org/build/3Dmol-min.js"


def _get_3dmol_script_tag() -> str:
    """Return an inline <script> tag if the bundled file exists, else a CDN tag."""
    if _3DMOL_BUNDLED_PATH.exists():
        js = _3DMOL_BUNDLED_PATH.read_text(encoding="utf-8")
        return f"<script>{js}</script>"
    return f'<script src="{_3DMOL_CDN}"></script>'


def escape_html(s: str) -> str:
    return (
        s.replace("&", "&amp;")
         .replace("<", "&lt;")
         .replace(">", "&gt;")
         .replace('"', "&quot;")
         .replace("'", "&#39;")
    )


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def build_summary_table(chains_data: Dict[str, List[dict]]) -> str:
    rows = ""
    for chain_id in sorted(chains_data):
        residues = chains_data[chain_id]
        n, avg, mn = chain_stats(residues)
        confidence = (
            "Very High" if avg >= 90 else
            "High" if avg >= 70 else
            "Medium" if avg >= 50 else "Low"
        )
        badge_class = (
            "badge-veryhigh" if avg >= 90 else
            "badge-high" if avg >= 70 else
            "badge-medium" if avg >= 50 else "badge-low"
        )
        rows += f"""
        <tr>
          <td><strong>Chain {escape_html(chain_id)}</strong></td>
          <td>{n}</td>
          <td>{avg}</td>
          <td>{mn}</td>
          <td><span class="badge {badge_class}">{confidence}</span></td>
        </tr>"""
    return f"""
    <table class="summary-table">
      <thead>
        <tr>
          <th>Chain</th><th>Residues</th>
          <th>Avg pLDDT</th><th>Min pLDDT</th><th>Confidence</th>
        </tr>
      </thead>
      <tbody>{rows}
      </tbody>
    </table>"""


def build_viewer_section(index: int, pdb_id: str, pdb_text: str) -> str:
    chains_data = parse_pdb_chains(pdb_text)
    summary_table = build_summary_table(chains_data)

    # Embed chain pLDDT data for the JS chart
    chart_data = {
        chain_id: [{"resnum": r["resnum"], "plddt": r["plddt"]}
                   for r in residues]
        for chain_id, residues in sorted(chains_data.items())
    }
    chart_data_json = json.dumps(chart_data)

    # Escape PDB for JS template literal
    pdb_js = pdb_text.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")
    preview = escape_html("\n".join(pdb_text.splitlines()[:40]))

    return f"""
    <section class="card">
      <div class="card-header">
        <h2 class="card-title">{escape_html(pdb_id)}</h2>
      </div>

      {summary_table}

      <div class="viewer-toolbar">
        <span class="toolbar-label">Color:</span>
        <div class="btn-group">
          <button class="btn active" data-viewer="{index}" data-mode="plddt">pLDDT</button>
          <button class="btn" data-viewer="{index}" data-mode="spectrum">Spectrum</button>
          <button class="btn" data-viewer="{index}" data-mode="surface">Surface</button>
        </div>
      </div>

      <div id="viewer_{index}" class="viewer"></div>

      <div class="plddt-legend">
        <span class="legend-item"><span class="legend-dot" style="background:#003F8A"></span>Very High (&ge;90)</span>
        <span class="legend-item"><span class="legend-dot" style="background:#4CB5F5"></span>High (70&ndash;90)</span>
        <span class="legend-item"><span class="legend-dot" style="background:#F5D76E"></span>Medium (50&ndash;70)</span>
        <span class="legend-item"><span class="legend-dot" style="background:#FF7D45"></span>Low (&lt;50)</span>
      </div>

      <div class="chart-section">
        <div class="section-title">Per-Residue pLDDT</div>
        <canvas id="chart_{index}" class="plddt-chart"></canvas>
      </div>

      <details class="pdb-preview">
        <summary>PDB preview (first 40 lines)</summary>
        <pre>{preview}</pre>
      </details>

      <script>
        (function () {{
          var pdb = `{pdb_js}`;
          var chainData = {chart_data_json};
          var idx = {index};

          function applyPlddt(viewer) {{
            viewer.setStyle({{}}, {{ cartoon: {{ colorfunc: function(atom) {{
              var b = atom.b;
              if (b >= 90) return '#003F8A';
              if (b >= 70) return '#4CB5F5';
              if (b >= 50) return '#F5D76E';
              return '#FF7D45';
            }} }} }});
            viewer.removeAllSurfaces();
            viewer.render();
          }}

          function drawChart() {{
            var canvas = document.getElementById("chart_" + idx);
            var dpr = window.devicePixelRatio || 1;
            var cssW = Math.round(canvas.getBoundingClientRect().width) || 800;
            var cssH = 200;
            canvas.width = cssW * dpr;
            canvas.height = cssH * dpr;
            canvas.style.width = cssW + 'px';
            canvas.style.height = cssH + 'px';
            var ctx = canvas.getContext("2d");
            ctx.scale(dpr, dpr);

            var W = cssW, H = cssH;

            var allResnums = [];
            Object.values(chainData).forEach(function(rs) {{
              rs.forEach(function(r) {{ allResnums.push(r.resnum); }});
            }});
            if (!allResnums.length) return;
            var minRes = Math.min.apply(null, allResnums);
            var maxRes = Math.max.apply(null, allResnums);

            var padL = 50, padR = 16, padT = 28, padB = 38;
            var plotW = W - padL - padR;
            var plotH = H - padT - padB;

            function xPos(r) {{ return padL + (r - minRes) / Math.max(maxRes - minRes, 1) * plotW; }}
            function yPos(v) {{ return padT + (1 - v / 100) * plotH; }}

            // Confidence bands
            var bands = [
              {{ min: 90, max: 100, color: 'rgba(0,63,138,0.08)' }},
              {{ min: 70, max: 90,  color: 'rgba(76,181,245,0.10)' }},
              {{ min: 50, max: 70,  color: 'rgba(245,215,110,0.12)' }},
              {{ min: 0,  max: 50,  color: 'rgba(255,125,69,0.08)' }},
            ];
            bands.forEach(function(b) {{
              ctx.fillStyle = b.color;
              ctx.fillRect(padL, yPos(b.max), plotW, yPos(b.min) - yPos(b.max));
            }});

            // Threshold lines
            ctx.setLineDash([4, 4]);
            ctx.lineWidth = 1;
            [90, 70, 50].forEach(function(v) {{
              ctx.strokeStyle = 'rgba(0,0,0,0.10)';
              ctx.beginPath();
              ctx.moveTo(padL, yPos(v));
              ctx.lineTo(padL + plotW, yPos(v));
              ctx.stroke();
            }});
            ctx.setLineDash([]);

            // Y-axis ticks and labels
            ctx.fillStyle = '#888';
            ctx.font = '11px sans-serif';
            ctx.textAlign = 'right';
            [0, 50, 70, 90, 100].forEach(function(v) {{
              ctx.fillText(v, padL - 6, yPos(v) + 4);
            }});

            // Y-axis title
            ctx.save();
            ctx.translate(14, padT + plotH / 2);
            ctx.rotate(-Math.PI / 2);
            ctx.textAlign = 'center';
            ctx.fillStyle = '#666';
            ctx.font = '11px sans-serif';
            ctx.fillText('pLDDT', 0, 0);
            ctx.restore();

            // X-axis ticks
            ctx.textAlign = 'center';
            ctx.fillStyle = '#888';
            ctx.font = '11px sans-serif';
            var step = Math.ceil((maxRes - minRes) / 8 / 10) * 10 || 10;
            for (var r = Math.ceil(minRes / step) * step; r <= maxRes; r += step) {{
              ctx.fillText(r, xPos(r), padT + plotH + 16);
            }}
            ctx.fillStyle = '#666';
            ctx.fillText('Residue position', padL + plotW / 2, H - 4);

            // Chain lines
            var chainIds = Object.keys(chainData).sort();
            var palette = ['#0f3460', '#c2410c'];
            chainIds.forEach(function(chainId, ci) {{
              var residues = chainData[chainId];
              if (!residues.length) return;
              var col = palette[ci % palette.length];
              ctx.strokeStyle = col;
              ctx.lineWidth = 1.5;
              ctx.beginPath();
              residues.forEach(function(r, i) {{
                var x = xPos(r.resnum), y = yPos(r.plddt);
                if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
              }});
              ctx.stroke();
            }});

            // Chain legend (top-left, outside plot data area)
            chainIds.forEach(function(chainId, ci) {{
              var col = palette[ci % palette.length];
              var lx = padL + 8;
              var ly = padT - 16 + ci * 0;
              if (chainIds.length > 1) {{
                lx = padL + 8 + ci * 80;
              }}
              ctx.strokeStyle = col; ctx.lineWidth = 2;
              ctx.beginPath(); ctx.moveTo(lx, ly); ctx.lineTo(lx + 18, ly); ctx.stroke();
              ctx.fillStyle = col; ctx.font = '11px sans-serif'; ctx.textAlign = 'left';
              ctx.fillText('Chain ' + chainId, lx + 22, ly + 4);
            }});

            // Plot border
            ctx.strokeStyle = '#e2e8f0'; ctx.lineWidth = 1;
            ctx.strokeRect(padL, padT, plotW, plotH);
          }}

          // Defer until DOM is laid out so viewer div has real dimensions.
          // $3Dmolpromise exists only when 3Dmol.js is loaded via CDN <script src>;
          // when inlined, $3Dmol is available immediately — use DOMContentLoaded instead.
          function initViewer() {{
            var el = document.getElementById("viewer_" + idx);
            var viewer = $3Dmol.createViewer(el, {{ backgroundColor: "#f8f9fa" }});
            viewer.resize();
            viewer.addModel(pdb, "pdb");
            applyPlddt(viewer);
            viewer.zoomTo();
            viewer.render();

            document.querySelectorAll('.btn[data-viewer="' + idx + '"]').forEach(function(btn) {{
              btn.addEventListener('click', function() {{
                document.querySelectorAll('.btn[data-viewer="' + idx + '"]').forEach(function(b) {{
                  b.classList.remove('active');
                }});
                this.classList.add('active');
                var mode = this.dataset.mode;
                if (mode === 'plddt') applyPlddt(viewer);
                else if (mode === 'spectrum') {{
                  viewer.setStyle({{}}, {{ cartoon: {{ colorscheme: 'spectrum' }} }});
                  viewer.removeAllSurfaces(); viewer.render();
                }} else if (mode === 'surface') {{
                  viewer.setStyle({{}}, {{ cartoon: {{ colorscheme: 'spectrum' }} }});
                  viewer.removeAllSurfaces();
                  viewer.addSurface($3Dmol.SurfaceType.VDW, {{ opacity: 0.7, colorscheme: 'spectrum' }});
                  viewer.render();
                }}
              }});
            }});

            drawChart();
          }}

          function deferInit() {{
            // requestAnimationFrame ensures the browser has painted at least once,
            // so the viewer div has its final CSS dimensions.
            requestAnimationFrame(function() {{ requestAnimationFrame(initViewer); }});
          }}

          if (typeof $3Dmolpromise !== 'undefined') {{
            $3Dmolpromise.then(deferInit);
          }} else if (document.readyState === 'loading') {{
            document.addEventListener('DOMContentLoaded', deferInit);
          }} else {{
            deferInit();
          }}
        }})();
      </script>
    </section>
"""


def generate_html_report(
    pdb_paths: List[Path],
    out_html: Path,
    title: str = "ABB3 Structure Predictions",
    use_plm: bool = True,
) -> None:
    sections = "".join(
        build_viewer_section(i, p.stem, read_text(p))
        for i, p in enumerate(pdb_paths)
    )

    if not sections:
        sections = "<p class='empty'>No PDB files were found.</p>"

    script_tag = _get_3dmol_script_tag()

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{escape_html(title)}</title>
  {script_tag}
  <style>
    :root {{
      --primary: #1a1a2e;
      --accent:  #0f3460;
      --bg:      #f0f2f5;
    }}

    *, *::before, *::after {{ box-sizing: border-box; }}

    body {{
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                   "Helvetica Neue", Arial, sans-serif;
      background: var(--bg);
      color: #1a1a2e;
    }}

    header {{
      background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
      color: #fff;
      padding: 28px 32px 22px;
    }}
    header h1 {{ margin: 0 0 8px; font-size: 1.8rem; font-weight: 700; letter-spacing: -0.02em; }}
    header p  {{ margin: 0; font-size: 0.93rem; opacity: 0.75; }}
    .mode-badge {{
      display: inline-block;
      margin-bottom: 10px;
      padding: 3px 10px;
      border-radius: 20px;
      font-size: 0.78rem;
      font-weight: 600;
      letter-spacing: 0.02em;
    }}
    .mode-lm   {{ background: rgba(76,181,245,0.25); color: #7dd3fc; border: 1px solid rgba(76,181,245,0.4); }}
    .mode-plain {{ background: rgba(255,255,255,0.12); color: #e2e8f0; border: 1px solid rgba(255,255,255,0.25); }}

    main {{
      max-width: 1100px;
      margin: 28px auto;
      padding: 0 20px 48px;
      display: grid;
      gap: 24px;
    }}

    .card {{
      background: #fff;
      border-radius: 14px;
      padding: 20px 24px;
      box-shadow: 0 2px 12px rgba(0,0,0,0.08);
    }}

    .card-header {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      margin-bottom: 16px;
    }}
    .card-title {{
      margin: 0;
      font-size: 1.15rem;
      font-weight: 600;
      color: var(--accent);
      word-break: break-all;
    }}

    /* Summary table */
    .summary-table {{
      width: 100%;
      border-collapse: collapse;
      font-size: 0.88rem;
      margin-bottom: 16px;
    }}
    .summary-table th {{
      background: #f7f9fc;
      color: #555;
      font-weight: 600;
      padding: 8px 12px;
      text-align: left;
      border-bottom: 2px solid #e2e8f0;
    }}
    .summary-table td {{
      padding: 7px 12px;
      border-bottom: 1px solid #f0f2f5;
    }}
    .badge {{
      display: inline-block;
      padding: 2px 8px;
      border-radius: 20px;
      font-size: 0.78rem;
      font-weight: 600;
    }}
    .badge-veryhigh {{ background: #dbeafe; color: #003F8A; }}
    .badge-high     {{ background: #e0f2fe; color: #0369a1; }}
    .badge-medium   {{ background: #fef9c3; color: #854d0e; }}
    .badge-low      {{ background: #ffedd5; color: #9a3412; }}

    /* Viewer toolbar */
    .viewer-toolbar {{
      display: flex;
      align-items: center;
      gap: 10px;
      margin-bottom: 10px;
    }}
    .toolbar-label {{ font-size: 0.85rem; color: #666; }}
    .btn-group {{ display: flex; gap: 4px; }}
    .btn {{
      padding: 5px 12px;
      font-size: 0.82rem;
      border: 1.5px solid #c6d4e8;
      background: #fff;
      color: #444;
      border-radius: 6px;
      cursor: pointer;
      transition: background 0.15s, border-color 0.15s;
    }}
    .btn:hover {{ background: #eef3fb; }}
    .btn.active {{ background: var(--accent); color: #fff; border-color: var(--accent); }}

    .viewer {{
      position: relative;
      width: 100%;
      height: 420px;
      border-radius: 10px;
      border: 1px solid #e2e8f0;
      overflow: hidden;
      margin-bottom: 10px;
    }}
    .viewer > div {{
      width: 100% !important;
      height: 100% !important;
    }}
    .viewer canvas {{
      width: 100% !important;
      height: 100% !important;
    }}

    /* pLDDT legend */
    .plddt-legend {{
      display: flex;
      flex-wrap: wrap;
      gap: 14px;
      margin-bottom: 16px;
      font-size: 0.82rem;
      color: #555;
    }}
    .legend-item {{ display: flex; align-items: center; gap: 5px; }}
    .legend-dot  {{ width: 10px; height: 10px; border-radius: 50%; flex-shrink: 0; }}

    /* pLDDT chart */
    .chart-section {{ margin-bottom: 16px; }}
    .section-title {{
      font-size: 0.9rem;
      font-weight: 600;
      color: var(--accent);
      border-left: 4px solid var(--accent);
      padding-left: 8px;
      margin-bottom: 8px;
    }}
    .plddt-chart {{
      width: 100%;
      height: 200px;
      display: block;
      border-radius: 8px;
      border: 1px solid #e2e8f0;
    }}

    /* PDB preview */
    .pdb-preview {{ margin-top: 8px; }}
    .pdb-preview summary {{ cursor: pointer; font-size: 0.85rem; color: #666; user-select: none; }}
    pre {{
      background: #f7f9fc;
      border: 1px solid #e2e8f0;
      border-radius: 8px;
      padding: 12px 14px;
      font-size: 0.78rem;
      overflow: auto;
      max-height: 240px;
      margin-top: 8px;
    }}

    .empty {{ text-align: center; color: #777; font-size: 1rem; padding: 48px 0; }}
  </style>
</head>
<body>
  <header>
    <div class="mode-badge {'mode-lm' if use_plm else 'mode-plain'}">
      {'ABB3-LM &mdash; language model enhanced' if use_plm else 'ABB3 &mdash; standard mode'}
    </div>
    <h1>{escape_html(title)}</h1>
    <p>{'Predicted with ABodyBuilder3-LM using ProtT5 embeddings. pLDDT confidence scores reflect the language-model-enhanced prediction.' if use_plm else 'Predicted with ABodyBuilder3 (plain mode, no language model). pLDDT confidence scores reflect standard structure prediction.'}</p>
  </header>

  <main>
    {sections}
  </main>
</body>
</html>
"""
    out_html.write_text(html, encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Node 4: Export PDBs and generate interactive HTML report."
    )
    parser.add_argument("--inputs", required=True, help="Directory of Node 3 .pdb files")
    parser.add_argument("--outputs", required=True, help="Directory to write report.html and PDBs")
    parser.add_argument("--title", default="ABB3 Structure Predictions", help="Report title")
    parser.add_argument("--use-plm", default="1", help="'1' for ABB3-LM mode, '0' for plain ABB3")
    args = parser.parse_args()

    inputs_dir = Path(args.inputs)
    outputs_dir = Path(args.outputs)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    pdb_in = sorted(inputs_dir.glob("*.pdb"))
    if not pdb_in:
        raise RuntimeError(f"No .pdb files found in {inputs_dir}")

    print(f"[Node 4] Found {len(pdb_in)} PDB file(s).")

    copied: List[Path] = []
    for src in pdb_in:
        dst = outputs_dir / src.name
        shutil.copy2(src, dst)
        copied.append(dst)
        print(f"[Node 4] Copied {src.name}")

    use_plm = args.use_plm.strip() == "1"
    out_html = outputs_dir / "report.html"
    generate_html_report(copied, out_html=out_html, title=args.title, use_plm=use_plm)

    print(f"[Node 4] HTML report written → {out_html}")
    print(f"[Node 4] To view: open report.html in any modern browser.")
    print(f"[Node 4] Done.")


if __name__ == "__main__":
    main()
