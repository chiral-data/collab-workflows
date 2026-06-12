"""
Node 4 – PDB Export & Interactive HTML Visualization
=====================================================
Copies each PDB file from Node 3 into the output directory and
generates a single self-contained HTML report with an interactive
py3Dmol viewer for every predicted structure.

The report works in any modern web browser with no Python runtime
required – all PDB data is embedded directly in the HTML.

Outputs
-------
results/node4/report.html          – interactive viewer
results/node4/<pair_id>.pdb        – copied PDB files (for download links)
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import List


# ---------------------------------------------------------------------------
# HTML helpers
# ---------------------------------------------------------------------------

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


def build_viewer_section(index: int, pdb_id: str, pdb_text: str) -> str:
    """Return the HTML block for one structure card."""
    # Escape backticks so the PDB can be embedded as a JS template literal
    pdb_js = pdb_text.replace("`", "\\`").replace("\\", "\\\\").replace("${", "\\${")
    preview = escape_html("\n".join(pdb_text.splitlines()[:40]))

    return f"""
    <section class="card">
      <div class="card-header">
        <h2 class="card-title">{escape_html(pdb_id)}</h2>
        <a class="btn-download" href="{escape_html(pdb_id)}.pdb" download>⬇ Download PDB</a>
      </div>
      <div id="viewer_{index}" class="viewer"></div>
      <details class="pdb-preview">
        <summary>PDB preview (first 40 lines)</summary>
        <pre>{preview}</pre>
      </details>
      <script>
        (function () {{
          const el = document.getElementById("viewer_{index}");
          const viewer = $3Dmol.createViewer(el, {{ backgroundColor: "#f8f9fa" }});
          const pdb = `{pdb_js}`;
          viewer.addModel(pdb, "pdb");
          viewer.setStyle({{}}, {{ cartoon: {{ colorscheme: "spectrum" }} }});
          viewer.zoomTo();
          viewer.render();
        }})();
      </script>
    </section>
"""


def generate_html_report(
    pdb_paths: List[Path],
    out_html: Path,
    title: str = "ABB3 Structure Predictions",
) -> None:
    sections = "".join(
        build_viewer_section(i, p.stem, read_text(p))
        for i, p in enumerate(pdb_paths)
    )

    if not sections:
        sections = "<p class='empty'>No PDB files were found.</p>"

    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>{escape_html(title)}</title>
  <script src="https://3dmol.org/build/3Dmol-min.js"></script>
  <style>
    *, *::before, *::after {{ box-sizing: border-box; }}

    body {{
      margin: 0;
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
                   "Helvetica Neue", Arial, sans-serif;
      background: #f0f2f5;
      color: #1a1a2e;
    }}

    header {{
      background: linear-gradient(135deg, #1a1a2e 0%, #16213e 60%, #0f3460 100%);
      color: #fff;
      padding: 28px 32px 22px;
    }}

    header h1 {{
      margin: 0 0 6px;
      font-size: 1.8rem;
      font-weight: 700;
      letter-spacing: -0.02em;
    }}

    header p {{
      margin: 0;
      font-size: 0.93rem;
      opacity: 0.75;
    }}

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
      box-shadow: 0 2px 12px rgba(0, 0, 0, 0.08);
    }}

    .card-header {{
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      margin-bottom: 14px;
    }}

    .card-title {{
      margin: 0;
      font-size: 1.15rem;
      font-weight: 600;
      color: #0f3460;
      word-break: break-all;
    }}

    .btn-download {{
      flex-shrink: 0;
      text-decoration: none;
      color: #0f3460;
      border: 1.5px solid #c6d4e8;
      border-radius: 8px;
      padding: 6px 14px;
      font-size: 0.88rem;
      font-weight: 500;
      transition: background 0.15s, border-color 0.15s;
    }}

    .btn-download:hover {{
      background: #eef3fb;
      border-color: #0f3460;
    }}

    .viewer {{
      width: 100%;
      height: 440px;
      border-radius: 10px;
      border: 1px solid #e2e8f0;
      overflow: hidden;
    }}

    .pdb-preview {{
      margin-top: 14px;
    }}

    .pdb-preview summary {{
      cursor: pointer;
      font-size: 0.88rem;
      color: #555;
      user-select: none;
    }}

    pre {{
      background: #f7f9fc;
      border: 1px solid #e2e8f0;
      border-radius: 8px;
      padding: 12px 14px;
      font-size: 0.8rem;
      overflow: auto;
      max-height: 260px;
      margin-top: 8px;
    }}

    .empty {{
      text-align: center;
      color: #777;
      font-size: 1rem;
      padding: 48px 0;
    }}
  </style>
</head>
<body>
  <header>
    <h1>{escape_html(title)}</h1>
    <p>Open this file in a browser to interactively explore predicted antibody structures.</p>
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
    args = parser.parse_args()

    inputs_dir = Path(args.inputs)
    outputs_dir = Path(args.outputs)
    outputs_dir.mkdir(parents=True, exist_ok=True)

    pdb_in = sorted(inputs_dir.glob("*.pdb"))
    if not pdb_in:
        raise RuntimeError(f"No .pdb files found in {inputs_dir}")

    print(f"[Node 4] Found {len(pdb_in)} PDB file(s).")

    # Copy PDBs to output dir so the HTML download links work alongside report.html
    copied: List[Path] = []
    for src in pdb_in:
        dst = outputs_dir / src.name
        shutil.copy2(src, dst)
        copied.append(dst)
        print(f"[Node 4] Copied {src.name}")

    out_html = outputs_dir / "report.html"
    generate_html_report(copied, out_html=out_html, title=args.title)

    print(f"[Node 4] HTML report written → {out_html}")
    print(f"[Node 4] To view: open report.html in any modern browser.")
    print(f"[Node 4] Done.")


if __name__ == "__main__":
    main()