import json
import os
from datetime import datetime

import numpy as np

confidence = json.load(open("./inputs/confidence.json"))
cif_data   = open("./inputs/structure.cif").read()

plddt_raw  = confidence["plddt_mean"]
ptm        = confidence["ptm"]
iptm       = confidence["iptm"]
plddt_per  = confidence.get("plddt_per_residue", [])
pae        = confidence.get("pae")

# Normalise to 0-100 (model may return 0-1 or 0-100)
plddt_mean = plddt_raw * 100 if plddt_raw <= 1.0 else plddt_raw

# pLDDT colour scale (matches colab + AlphaFold convention)
def plddt_colour(v):
    if v >= 90: return "#0053D6"
    if v >= 70: return "#65CBF3"
    if v >= 50: return "#FFDB13"
    return "#FF7D45"

# Residue colours for 3Dmol pLDDT colouring (JS array)
plddt_arr = plddt_per if plddt_per else []
if plddt_arr and max(plddt_arr) <= 1.0:
    plddt_arr = [v * 100 for v in plddt_arr]

plddt_colours_js = json.dumps([plddt_colour(v) for v in plddt_arr])

# PAE heatmap (Plotly)
pae_plot_div = ""
if pae:
    pae_np = np.array(pae)
    z_json = json.dumps(pae_np.tolist())
    pae_plot_div = f"""
    <div id="pae-plot" style="width:100%;height:450px;"></div>
    <script>
    Plotly.newPlot('pae-plot', [{{
        z: {z_json},
        type: 'heatmap',
        colorscale: 'Greens_r',
        zmin: 0, zmax: 30,
        colorbar: {{title: 'PAE (Å)', len: 0.8}}
    }}], {{
        title: 'Predicted Aligned Error (PAE)',
        xaxis: {{title: 'Scored residue'}},
        yaxis: {{title: 'Aligned residue', autorange: 'reversed'}},
        margin: {{t: 50, b: 60, l: 70, r: 20}}
    }}, {{responsive: true}});
    </script>"""

# Quality badge (0-100 scale)
def quality_badge(v):
    if v >= 90: return ("Excellent", "#0f766e")
    if v >= 70: return ("Good",      "#0369a1")
    if v >= 50: return ("Moderate",  "#fb7c3c")
    return ("Low", "#dc2626")

q_label, q_colour = quality_badge(plddt_mean)

cif_escaped = cif_data.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>ESMFold2 Structure Report</title>
  <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
  <script src="https://3dmol.org/build/3Dmol-min.js"></script>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
  <style>
    body {{ background: linear-gradient(135deg, #075985, #0284c7); min-height: 100vh; font-family: system-ui, sans-serif; }}
    .wrap {{ max-width: 1200px; margin: 0 auto; padding: 24px; }}
    .card {{ background: rgba(255,255,255,0.96); border-radius: 14px; box-shadow: 0 8px 32px rgba(0,0,0,.12); margin-bottom: 24px; }}
    .card-header {{ background: linear-gradient(135deg,#075985,#0284c7); color:#fff; border-radius: 14px 14px 0 0; padding: 18px 24px; font-weight: 600; font-size: 1.05rem; }}
    .card-body {{ padding: 24px; }}
    .metric {{ text-align:center; }}
    .metric .val {{ font-size: 2rem; font-weight: 700; color: #075985; }}
    .metric .lbl {{ font-size: 0.8rem; text-transform: uppercase; letter-spacing: 1px; color: #64748b; }}
    #viewer {{ height: 500px; background: #1a1a2e; border-radius: 10px; }}
    .legend span {{ display:inline-block; width:14px; height:14px; border-radius:3px; margin-right:4px; vertical-align:middle; }}
  </style>
</head>
<body>
<div class="wrap">

  <div class="card">
    <div class="card-body">
      <h2 class="mb-1" style="color:#075985;">ESMFold2 Structure Prediction Report</h2>
      <p class="text-muted mb-0">Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    </div>
  </div>

  <div class="row g-3 mb-4">
    <div class="col-md-3">
      <div class="card h-100"><div class="card-body metric">
        <div class="val">{plddt_mean:.2f}</div>
        <div class="lbl">Mean pLDDT</div>
        <span class="badge mt-2" style="background:{q_colour}">{q_label}</span>
      </div></div>
    </div>
    <div class="col-md-3">
      <div class="card h-100"><div class="card-body metric">
        <div class="val">{ptm:.3f}</div>
        <div class="lbl">pTM</div>
      </div></div>
    </div>
    <div class="col-md-3">
      <div class="card h-100"><div class="card-body metric">
        <div class="val">{iptm:.3f}</div>
        <div class="lbl">ipTM</div>
      </div></div>
    </div>
    <div class="col-md-3">
      <div class="card h-100"><div class="card-body metric">
        <div class="val">{len(plddt_arr)}</div>
        <div class="lbl">Residues</div>
      </div></div>
    </div>
  </div>

  <div class="card">
    <div class="card-header">3D Structure Viewer — coloured by pLDDT</div>
    <div class="card-body">
      <div class="legend mb-3 small">
        <span style="background:#FF7D45"></span>&lt;50 &nbsp;
        <span style="background:#FFDB13"></span>50–70 &nbsp;
        <span style="background:#65CBF3"></span>70–90 &nbsp;
        <span style="background:#0053D6"></span>&gt;90
      </div>
      <div id="viewer"></div>
    </div>
  </div>

  {"<div class='card'><div class='card-header'>Predicted Aligned Error (PAE)</div><div class='card-body'>" + pae_plot_div + "</div></div>" if pae_plot_div else ""}

</div>
<script>
(function() {{
  var cif = `{cif_escaped}`;
  var colours = {plddt_colours_js};
  var viewer = $3Dmol.createViewer(document.getElementById('viewer'), {{backgroundColor:'#1a1a2e'}});
  viewer.addModel(cif, 'mmcif');
  if (colours.length > 0) {{
    viewer.setStyle({{}}, {{}});
    colours.forEach(function(c, i) {{
      viewer.setStyle({{resi: i+1}}, {{cartoon: {{color: c}}}});
    }});
  }} else {{
    viewer.setStyle({{}}, {{cartoon: {{colorscheme: 'chainHetatm'}}}});
  }}
  viewer.addStyle({{hetflag: true}}, {{stick: {{}}}});
  viewer.zoomTo();
  viewer.render();
}})();
</script>
</body>
</html>"""

os.makedirs("./outputs", exist_ok=True)
with open("./outputs/report.html", "w") as f:
    f.write(html)

print("Wrote ./outputs/report.html", flush=True)
