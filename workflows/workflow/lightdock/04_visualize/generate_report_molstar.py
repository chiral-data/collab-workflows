#!/usr/bin/env python3
"""Generate a self-contained HTML docking report with 3D viewer and ranked pose table."""

import json
import os
import sys


def read_first_model_pdb(path):
    """Return ATOM/HETATM lines from the first MODEL and the ordered list of unique chain IDs."""
    lines = []
    chains = []
    seen_chains = set()
    in_model = False
    has_model = False

    with open(path) as f:
        for line in f:
            if line.startswith("MODEL"):
                has_model = True
                in_model = True
                continue
            if line.startswith("ENDMDL"):
                if in_model:
                    break
                continue
            if in_model or not has_model:
                if line.startswith(("ATOM", "HETATM")):
                    lines.append(line.rstrip())
                    chain = line[21] if len(line) > 21 else ""
                    if chain.strip() and chain not in seen_chains:
                        seen_chains.add(chain)
                        chains.append(chain)

    return "\n".join(lines), chains


def main():
    rank_path = os.path.join("inputs", "rank_results.json")
    pdb_path = os.path.join("inputs", "top_predictions.pdb")

    if not os.path.exists(rank_path):
        print(f"ERROR: {rank_path} not found", flush=True)
        sys.exit(1)
    if not os.path.exists(pdb_path):
        print(f"ERROR: {pdb_path} not found", flush=True)
        sys.exit(1)

    with open(rank_path) as f:
        rank_data = json.load(f)

    top_poses = rank_data.get("top_poses", [])
    top_pdb_content, chains = read_first_model_pdb(pdb_path)
    top_pdb_content = top_pdb_content.replace("`", "\\`").replace("${", "\\${")

    table_rows = ""
    for pose in top_poses:
        rank = pose.get("rank", "")
        score = pose.get("score", "")
        score_str = f"{score:.4f}" if isinstance(score, float) else str(score)
        table_rows += f"<tr><td>{rank}</td><td>{score_str}</td></tr>\n"

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>LightDock Docking Report</title>
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.css">
  <script src="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.js"></script>
  <style>
    body {{ font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; color: #333; }}
    h1 {{ color: #2c3e50; }}
    .container {{ max-width: 1100px; margin: 0 auto; }}
    .summary {{ background: #fff; border-radius: 8px; padding: 16px 24px; margin-bottom: 24px;
                box-shadow: 0 1px 4px rgba(0,0,0,0.1); }}
    .summary table {{ border-collapse: collapse; }}
    .summary td {{ padding: 4px 16px 4px 0; }}
    .summary td:first-child {{ font-weight: bold; color: #555; }}
    #viewer-container {{ background: #fff; border-radius: 8px; padding: 16px 24px;
                         margin-bottom: 24px; box-shadow: 0 1px 4px rgba(0,0,0,0.1); }}
    #viewer {{ width: 100%; height: 480px; position: relative; }}
    .rankings {{ background: #fff; border-radius: 8px; padding: 16px 24px;
                 box-shadow: 0 1px 4px rgba(0,0,0,0.1); }}
    .rankings table {{ width: 100%; border-collapse: collapse; }}
    .rankings th {{ background: #2c3e50; color: #fff; padding: 10px 14px; text-align: left; }}
    .rankings td {{ padding: 8px 14px; border-bottom: 1px solid #eee; }}
    .rankings tr:hover td {{ background: #f0f4f8; }}
    .best-score {{ color: #27ae60; font-weight: bold; }}
  </style>
</head>
<body>
<div class="container">
  <h1>LightDock Protein-Protein Docking Report</h1>

  <div class="summary">
    <h2>Run Summary</h2>
    <table>
      <tr><td>Scoring function</td><td>{rank_data.get('scoring_function', 'N/A')}</td></tr>
      <tr><td>Swarms</td><td>{rank_data.get('num_swarms', 'N/A')}</td></tr>
      <tr><td>Steps</td><td>{rank_data.get('steps', 'N/A')}</td></tr>
      <tr><td>Conformations per swarm</td><td>{rank_data.get('num_conformations', 'N/A')}</td></tr>
      <tr><td>Top poses reported</td><td>{len(top_poses)}</td></tr>
      {'<tr><td>Best score</td><td class="best-score">' + f"{top_poses[0]['score']:.4f}" + '</td></tr>' if top_poses else ''}
    </table>
  </div>

  <div id="viewer-container">
    <h2>Top-Ranked Docking Pose (3D View)</h2>
    <div id="viewer"></div>
  </div>

  <div class="rankings">
    <h2>Ranked Docking Poses</h2>
    <table>
      <thead>
        <tr>
          <th>Rank</th>
          <th>Score</th>
        </tr>
      </thead>
      <tbody>
        {table_rows}
      </tbody>
    </table>
  </div>
</div>

<script>
(function() {{
  var pdbData = `{top_pdb_content}`;
  if (!pdbData.trim()) {{
    document.getElementById('viewer-container').innerHTML +=
      '<p style="color:#999">No structure data available for visualization.</p>';
    return;
  }}
  molstar.Viewer.create('viewer', {{
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowLeftPanel: false,
    layoutShowSequence: false,
    layoutShowLog: false,
    layoutShowRemoteState: false,
    viewportShowAnimation: false,
    viewportShowExpand: true,
    viewportShowSelectionMode: false,
  }}).then(function(viewer) {{
    return viewer.loadStructureFromData(pdbData, 'pdb');
  }}).catch(function(err) {{
    console.error('Mol* error:', err);
    document.getElementById('viewer').innerHTML =
      '<p style="color:red;padding:16px">Viewer error: ' + err + '</p>';
  }});
}})();
</script>
</body>
</html>
"""

    with open("report.html", "w") as f:
        f.write(html)

    print(f"Report written to report.html ({len(top_poses)} poses).", flush=True)


if __name__ == "__main__":
    main()
