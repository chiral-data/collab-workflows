#!/usr/bin/env python3
"""
Node 05: Generate Pipeline Report.

Produces a self-contained report.html with five sections:

  1. Pipeline summary card  — sequence length, model count, pocket count,
                              QC-passed pockets, docked poses
  2. Structure confidence   — Boltz-2 pLDDT histogram + PAE heatmap (Plotly)
  3. Pocket discovery table — P2Rank pockets sortable by score, probability,
                              mean/min pLDDT, center coordinates
  4. 3D viewer              — PDBe-Molstar (CDN): protein cartoon coloured by
                              pLDDT B-factors, selected pocket sphere,
                              top docking pose as ball-and-stick
  5. Methods & caveats      — pipeline parameters + pLDDT caveat (Eguida &
                              Rognan 2023: ≥70 is necessary but not sufficient)

All Plotly traces are inlined as JSON; no server-side compute after generation.
The page is fully self-contained except for three CDN resources:
  - Bootstrap 5 CSS/JS
  - Plotly.js
  - PDBe-Molstar (viewer + CSS)
"""

import glob
import json
import os
import sys
from datetime import datetime
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Data loading helpers
# ---------------------------------------------------------------------------

def load_model_id() -> str:
    p = Path("selected_model_id.txt")
    if p.exists():
        return p.read_text().strip()
    # Fall back to first available plddt file
    candidates = sorted(glob.glob("plddt_*.npz"))
    if candidates:
        stem = Path(candidates[0]).stem  # e.g. "plddt_model_0"
        return stem[len("plddt_"):]       # e.g. "model_0"
    return "model_0"


def load_plddt(model_id: str) -> np.ndarray:
    """Return per-residue pLDDT on the 0-100 scale."""
    candidates = glob.glob(f"plddt_{model_id}.npz") + sorted(glob.glob("plddt_*.npz"))
    if not candidates:
        return np.array([])
    npz = np.load(candidates[0])
    key = "plddt" if "plddt" in npz.files else npz.files[0]
    arr = npz[key].astype(float).flatten()
    if arr.max() <= 1.01:
        arr = arr * 100.0
    return arr


def load_pae(model_id: str) -> np.ndarray | None:
    """Return the PAE matrix (N×N) or None if unavailable."""
    candidates = glob.glob(f"pae_{model_id}.npz") + sorted(glob.glob("pae_*.npz"))
    if not candidates:
        return None
    npz = np.load(candidates[0])
    key = "pae" if "pae" in npz.files else npz.files[0]
    arr = npz[key].astype(float)
    if arr.ndim == 1:
        # Some Boltz-2 versions store the upper triangle flattened
        n = int(round((-1 + (1 + 8 * len(arr)) ** 0.5) / 2))
        if n * (n + 1) // 2 == len(arr):
            mat = np.zeros((n, n))
            idx = np.triu_indices(n)
            mat[idx] = arr
            mat = mat + mat.T - np.diag(np.diag(mat))
            return mat
        # Square it if possible
        side = int(round(len(arr) ** 0.5))
        if side * side == len(arr):
            return arr.reshape(side, side)
        return None
    return arr


def load_affinity() -> dict:
    files = sorted(glob.glob("affinity_*.json"))
    if not files:
        return {}
    try:
        return json.loads(Path(files[0]).read_text())
    except Exception:
        return {}


def load_confidence() -> list[dict]:
    """Load all confidence_*.json files sorted by model index."""
    result = []
    for f in sorted(glob.glob("confidence_*.json")):
        try:
            d = json.loads(Path(f).read_text())
            d["_file"] = Path(f).name
            result.append(d)
        except Exception:
            pass
    return result


def load_input_summary() -> dict:
    p = Path("input_summary.json")
    if p.exists():
        try:
            return json.loads(p.read_text())
        except Exception:
            pass
    return {}


def load_pocket_qc() -> dict:
    p = Path("pocket_qc.json")
    if p.exists():
        try:
            return json.loads(p.read_text())
        except Exception:
            pass
    return {}


def load_docking_summary() -> dict:
    p = Path("docking_summary.json")
    if p.exists():
        try:
            return json.loads(p.read_text())
        except Exception:
            pass
    return {}


def load_top_pose_sdf() -> str | None:
    """Return the first (top-ranked) pose record from docked_poses.sdf, or None."""
    p = Path("docked_poses.sdf")
    if not p.exists() or p.stat().st_size == 0:
        return None
    text = p.read_text()
    end = text.find("$$$$")
    if end < 0:
        return text if text.strip() else None
    return text[:end + 4]


# ---------------------------------------------------------------------------
# pLDDT → B-factor PDB  (for Molstar colouring)
# ---------------------------------------------------------------------------

def receptor_pdb_with_bfactors(plddt: np.ndarray) -> str | None:
    """
    Patch the receptor.pdb B-factor column with pLDDT values so PDBe-Molstar
    can colour by pLDDT using its built-in 'b-factor' representation.
    Returns the patched PDB string, or None if receptor.pdb not found.
    """
    p = Path("receptor.pdb")
    if not p.exists() or len(plddt) == 0:
        return p.read_text() if p.exists() else None

    lines = p.read_text().splitlines()
    out = []
    res_seen: dict[tuple, int] = {}   # (chain, seqnum, icode) -> first atom index
    patched_res_idx = 0               # advances each time we see a new residue

    for line in lines:
        if not line.startswith(("ATOM", "HETATM")):
            out.append(line)
            continue
        chain  = line[21]
        seqnum = line[22:26].strip()
        icode  = line[26]
        key    = (chain, seqnum, icode)
        if key not in res_seen:
            res_seen[key] = patched_res_idx
            patched_res_idx += 1
        idx = res_seen[key]
        if idx < len(plddt):
            bfac = plddt[idx]
        else:
            bfac = 0.0
        # PDB B-factor column: cols 61-66 (1-based, 0-based 60:66), right-justified 6.2f
        patched = line[:60] + f"{bfac:6.2f}" + (line[66:] if len(line) > 66 else "")
        out.append(patched)

    return "\n".join(out)


def _escape_js(s: str) -> str:
    return s.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")


# ---------------------------------------------------------------------------
# HTML generation helpers
# ---------------------------------------------------------------------------

def _fmt(v, decimals=2):
    if v is None:
        return "—"
    try:
        return f"{float(v):.{decimals}f}"
    except (TypeError, ValueError):
        return str(v)


def _affinity_html(affinity: dict) -> str:
    if not affinity:
        return ""
    aff_val = affinity.get("affinity")
    if aff_val is None:
        return ""
    return (f'<p class="mt-2 mb-1" style="font-size:var(--font-sm);color:var(--text-secondary);">'
            f'<strong>Boltz-2 predicted affinity:</strong> '
            f'log<sub>10</sub>K<sub>D</sub> = {_fmt(aff_val, 2)}'
            f' <em style="color:var(--text-faint);">(orientation only — not validated for ranking)</em></p>')


def _badge(text: str, color: str) -> str:
    return f'<span class="pill" style="background:{color};">{text}</span>'


# ---------------------------------------------------------------------------
# Report sections
# ---------------------------------------------------------------------------

def _section_summary(input_summary, confidence_list, pocket_qc, docking_summary, plddt, affinity):
    n_models = len(confidence_list)
    entities = input_summary.get("entities", [])
    protein_entities = [e for e in entities if e.get("type") == "protein"]
    ligand_entities  = [e for e in entities if e.get("type") == "ligand"]

    seq_len = 0
    if protein_entities:
        seq_len = protein_entities[0].get("length", 0)

    n_pockets = len(pocket_qc.get("pockets", []))
    n_qc_pass = sum(1 for p in pocket_qc.get("pockets", []) if p.get("plddt_passes"))
    selected_rank = pocket_qc.get("selected_pocket_rank", "—")
    selected_pass = pocket_qc.get("selected_pocket_passes_qc")
    qc_badge = (_badge("QC PASS", "#16a34a") if selected_pass
                else _badge("QC WARN", "#f59e0b") if selected_pass is not None
                else _badge("N/A", "#94a3b8"))

    n_poses = docking_summary.get("num_poses_generated", "—")
    ligand_smiles = docking_summary.get("ligand_smiles", "")
    if not ligand_smiles and ligand_entities:
        ligand_smiles = ligand_entities[0].get("smiles", "")

    mean_plddt = f"{plddt.mean():.1f}" if len(plddt) > 0 else "—"

    return f"""
  <div class="section-card mb-4">
    <div class="section-header">Pipeline Summary</div>
    <div class="section-body">
      <div class="row g-3">
        <div class="col-6 col-md-2">
          <div class="stat-card">
            <div class="value">{seq_len if seq_len else "—"}</div>
            <div class="label">Residues</div>
          </div>
        </div>
        <div class="col-6 col-md-2">
          <div class="stat-card">
            <div class="value">{n_models}</div>
            <div class="label">Boltz-2 Models</div>
          </div>
        </div>
        <div class="col-6 col-md-2">
          <div class="stat-card">
            <div class="value">{mean_plddt}</div>
            <div class="label">Mean pLDDT</div>
          </div>
        </div>
        <div class="col-6 col-md-2">
          <div class="stat-card">
            <div class="value">{n_pockets}</div>
            <div class="label">Pockets Found</div>
          </div>
        </div>
        <div class="col-6 col-md-2">
          <div class="stat-card">
            <div class="value">{n_qc_pass} / {n_pockets}</div>
            <div class="label">Pockets QC ≥70</div>
          </div>
        </div>
        <div class="col-6 col-md-2">
          <div class="stat-card">
            <div class="value">{n_poses}</div>
            <div class="label">Docked Poses</div>
          </div>
        </div>
      </div>

      {"" if not ligand_smiles else f'<p class="mt-3 mb-1" style="font-size:var(--font-sm);color:var(--text-secondary);"><strong>Ligand SMILES:</strong> <code style="word-break:break-all;">{ligand_smiles}</code></p>'}
      {_affinity_html(affinity)}
      <p class="mt-2 mb-0" style="font-size:var(--font-sm);color:var(--text-secondary);">
        Selected pocket: rank {selected_rank}&nbsp;&nbsp;{qc_badge}
        &nbsp;&mdash;&nbsp;
        <em>QC badge reflects mean pocket pLDDT vs threshold {pocket_qc.get("threshold", 70.0)}</em>
      </p>
    </div>
  </div>"""


def _section_confidence(plddt, pae, confidence_list):
    """pLDDT histogram + PAE heatmap side by side."""
    if len(plddt) == 0 and pae is None:
        return ""

    plddt_vh_js = json.dumps(plddt[plddt >= 90].tolist())          if len(plddt) > 0 else "[]"
    plddt_h_js  = json.dumps(plddt[(plddt >= 70) & (plddt < 90)].tolist()) if len(plddt) > 0 else "[]"
    plddt_l_js  = json.dumps(plddt[(plddt >= 50) & (plddt < 70)].tolist()) if len(plddt) > 0 else "[]"
    plddt_vl_js = json.dumps(plddt[plddt < 50].tolist())           if len(plddt) > 0 else "[]"

    pae_section = ""
    pae_chart_js = ""
    if pae is not None:
        # Downsample large PAE matrices to ≤200×200 for the heatmap
        mat = pae
        if mat.shape[0] > 200:
            step = mat.shape[0] // 200 + 1
            mat = mat[::step, ::step]
        pae_z_js  = json.dumps(mat.tolist())
        pae_max   = round(float(pae.max()), 1)
        pae_section = '<div class="col-12 col-md-6"><div id="pae-chart" style="height:340px;"></div></div>'
        pae_chart_js = f"""Plotly.newPlot('pae-chart', [{{
  type: 'heatmap',
  z: {pae_z_js},
  colorscale: [[0,'#f0f7ff'],[0.25,'#a8cef0'],[0.5,'#4a9ade'],[0.75,'#1a6fba'],[1,'#075985']],
  zmin: 0, zmax: {pae_max},
  colorbar: {{ title: 'PAE (Å)', len: 0.8 }},
  hovertemplate: 'Res %{{x}},%{{y}}: %{{z:.1f}} Å<extra></extra>'
}}], {{
  xaxis: {{ title: 'Residue (aligned)' }},
  yaxis: {{ title: 'Residue (scored)', autorange: 'reversed' }},
  template: 'plotly_white', height: 340,
  margin: {{ t: 30, b: 50, l: 60, r: 20 }},
  title: {{ text: 'Predicted Aligned Error (PAE)', font: {{ size: 13 }} }}
}}, {{responsive: true}});"""

    plddt_chart_js = ""
    if len(plddt) > 0:
        bins_js = "{ start: 0, end: 100, size: 2 }"
        plddt_chart_js = f"""Plotly.newPlot('plddt-hist', [
  {{ type:'histogram', x:{plddt_vl_js}, xbins:{bins_js},
     marker:{{ color:'#dc2626' }}, name:'Very low <50',
     hovertemplate:'pLDDT %{{x:.0f}}: %{{y}} residues<extra></extra>' }},
  {{ type:'histogram', x:{plddt_l_js},  xbins:{bins_js},
     marker:{{ color:'#f59e0b' }}, name:'Low 50–70',
     hovertemplate:'pLDDT %{{x:.0f}}: %{{y}} residues<extra></extra>' }},
  {{ type:'histogram', x:{plddt_h_js},  xbins:{bins_js},
     marker:{{ color:'#16a34a' }}, name:'High 70–90',
     hovertemplate:'pLDDT %{{x:.0f}}: %{{y}} residues<extra></extra>' }},
  {{ type:'histogram', x:{plddt_vh_js}, xbins:{bins_js},
     marker:{{ color:'#1d4ed8' }}, name:'Very high ≥90',
     hovertemplate:'pLDDT %{{x:.0f}}: %{{y}} residues<extra></extra>' }}
], {{
  barmode: 'stack',
  xaxis: {{ title: 'pLDDT', range: [0, 100] }},
  yaxis: {{ title: 'Residue count' }},
  template: 'plotly_white', height: 340,
  margin: {{ t: 30, b: 50, l: 55, r: 20 }},
  title: {{ text: 'pLDDT Distribution (selected model)', font: {{ size: 13 }} }},
  legend: {{ orientation: 'h', y: -0.25 }},
  shapes: [
    {{ type:'line', x0:70, x1:70, y0:0, y1:1, yref:'paper',
       line:{{ color:'#94a3b8', dash:'dash', width:1.5 }} }},
    {{ type:'line', x0:90, x1:90, y0:0, y1:1, yref:'paper',
       line:{{ color:'#94a3b8', dash:'dot', width:1 }} }}
  ],
  annotations: [
    {{ x:70, y:1, yref:'paper', text:'70', showarrow:false,
       font:{{ size:10, color:'#94a3b8' }}, xanchor:'right', yanchor:'bottom' }},
    {{ x:90, y:1, yref:'paper', text:'90', showarrow:false,
       font:{{ size:10, color:'#94a3b8' }}, xanchor:'right', yanchor:'bottom' }}
  ]
}}, {{responsive: true}});"""

    plddt_col_width = "col-12 col-md-6" if pae is not None else "col-12"
    plddt_height    = "340px"

    return f"""
  <div class="section-card mb-4">
    <div class="section-header">Structure Confidence (Boltz-2)</div>
    <div class="section-body">
      <div class="row g-3">
        <div class="{plddt_col_width}">
          <div id="plddt-hist" style="height:{plddt_height};"></div>
        </div>
        {pae_section}
      </div>
      <p class="note mt-2">
        pLDDT colour bands: <span style="color:#1d4ed8;font-weight:600;">&#9632;</span> Very high ≥90 &nbsp;
        <span style="color:#16a34a;font-weight:600;">&#9632;</span> Confident 70–90 &nbsp;
        <span style="color:#f59e0b;font-weight:600;">&#9632;</span> Low 50–70 &nbsp;
        <span style="color:#dc2626;font-weight:600;">&#9632;</span> Very low &lt;50 &nbsp;
        (AlphaFold / ColabFold scale). Dashed line = 70 threshold used for pocket QC.
      </p>
    </div>
  </div>
__PLDDT_JS__{plddt_chart_js}
__PAE_JS__{pae_chart_js}"""


def _section_pocket_table(pocket_qc):
    pockets = pocket_qc.get("pockets", [])
    if not pockets:
        return ""

    selected_rank = pocket_qc.get("selected_pocket_rank")
    threshold = pocket_qc.get("threshold", 70.0)

    rows = ""
    for p in pockets:
        rank      = p.get("rank", "—")
        name      = p.get("name", "—")
        score     = _fmt(p.get("p2rank_score"), 3)
        prob      = _fmt(p.get("probability"), 3)
        mean_pl   = _fmt(p.get("plddt_mean"), 1)
        min_pl    = _fmt(p.get("plddt_min"), 1)
        n_res     = p.get("n_residues", "—")
        cx        = _fmt(p.get("center_x"), 2)
        cy        = _fmt(p.get("center_y"), 2)
        cz        = _fmt(p.get("center_z"), 2)
        passes    = p.get("plddt_passes")
        qc_badge  = (_badge("PASS", "#16a34a") if passes
                     else _badge("WARN", "#f59e0b") if passes is not None
                     else "—")
        selected  = rank == selected_rank
        row_style = ' style="background:#f0f7ff;"' if selected else ""
        sel_mark  = " ★" if selected else ""
        rows += (
            f"<tr{row_style}>"
            f"<td>{rank}{sel_mark}</td>"
            f"<td>{name}</td>"
            f"<td>{score}</td>"
            f"<td>{prob}</td>"
            f"<td>{mean_pl}</td>"
            f"<td>{min_pl}</td>"
            f"<td>{n_res}</td>"
            f"<td>{cx}</td>"
            f"<td>{cy}</td>"
            f"<td>{cz}</td>"
            f"<td>{qc_badge}</td>"
            f"</tr>\n"
        )

    return f"""
  <div class="section-card mb-4">
    <div class="section-header">Pocket Discovery (P2Rank) — sorted by score</div>
    <div class="section-body p-0">
      <div class="table-responsive">
        <table class="table table-sm table-hover mb-0" id="pocket-table">
          <thead>
            <tr>
              <th>Rank &#8597;</th>
              <th>Name</th>
              <th>Score &#8597;</th>
              <th>Probability &#8597;</th>
              <th>Mean pLDDT &#8597;</th>
              <th>Min pLDDT &#8597;</th>
              <th>Residues</th>
              <th>Center X</th>
              <th>Center Y</th>
              <th>Center Z</th>
              <th>QC (≥{threshold})</th>
            </tr>
          </thead>
          <tbody>{rows}</tbody>
        </table>
      </div>
      <p class="note px-4 pb-3 mt-2">
        &#9733; = selected pocket for docking.
        Click column headers to sort.
        pLDDT QC: mean per-residue pLDDT of pocket residues vs threshold {threshold}.
      </p>
    </div>
  </div>"""


def _section_viewer(plddt, pocket_qc, top_pose_sdf):
    """PDBe-Molstar viewer: protein coloured by pLDDT B-factors + top docking pose."""
    receptor_pdb_str = receptor_pdb_with_bfactors(plddt)
    if receptor_pdb_str is None and top_pose_sdf is None:
        return "", ""

    structures = []  # [(label, content, format, color)]

    if receptor_pdb_str is not None:
        structures.append(("Receptor (pLDDT B-factors)", _escape_js(receptor_pdb_str), "pdb", "#075985"))

    if top_pose_sdf is not None:
        structures.append(("Top docking pose (Uni-Mol)", _escape_js(top_pose_sdf), "sdf", "#0284c7"))

    if not structures:
        return "", ""

    toggle_buttons = "".join(
        f'<button id="toggle-struct-{i}" onclick="toggleStruct({i})" '
        f'class="pill" style="background:{color};border:none;padding:4px 14px;'
        f'cursor:pointer;margin-right:6px;margin-bottom:6px;">'
        f'{label}</button>'
        for i, (label, _, _, color) in enumerate(structures)
    )

    entries_js = ",\n".join(
        f'  {{label: {json.dumps(label)}, data: `{content}`, format: {json.dumps(fmt)}}}'
        for label, content, fmt, _ in structures
    )
    n = len(structures)

    # Pocket centre sphere annotation note
    pockets = pocket_qc.get("pockets", [])
    selected_rank = pocket_qc.get("selected_pocket_rank")
    selected_pocket = next((p for p in pockets if p.get("rank") == selected_rank), None)
    pocket_note = ""
    if selected_pocket:
        cx = _fmt(selected_pocket.get("center_x"), 2)
        cy = _fmt(selected_pocket.get("center_y"), 2)
        cz = _fmt(selected_pocket.get("center_z"), 2)
        pocket_note = (f" Pocket centre (rank {selected_rank}): "
                       f"({cx}, {cy}, {cz}) Å.")

    viewer_js = f"""(function() {{
  var structs = [
{entries_js}
  ];
  var visible = Array({n}).fill(true);
  var molPlugin = null;

  window.toggleStruct = function(idx) {{
    if (!molPlugin) return;
    visible[idx] = !visible[idx];
    var show = visible[idx];
    var hierarchy = molPlugin.managers.structure.hierarchy.current.structures;
    if (!hierarchy[idx]) return;
    var isHidden = !!hierarchy[idx].cell.state.isHidden;
    if ((show && isHidden) || (!show && !isHidden)) {{
      try {{
        molPlugin.managers.structure.hierarchy.toggleVisibility([hierarchy[idx]]);
      }} catch(e) {{
        try {{
          molstar.PluginCommands.State.ToggleVisibility(molPlugin, {{
            state: molPlugin.state.data,
            ref: hierarchy[idx].cell.transform.ref
          }});
        }} catch(e2) {{ console.warn('toggle failed:', e2); }}
      }}
    }}
    var btn = document.getElementById('toggle-struct-' + idx);
    if (btn) btn.style.opacity = show ? '1' : '0.35';
  }};

  molstar.Viewer.create('molstar-viewer', {{
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowLeftPanel: true,
    layoutShowSequence: false,
    layoutShowLog: false,
    layoutShowRemoteState: false,
    viewportShowAnimation: false,
    viewportShowExpand: true,
    viewportShowSelectionMode: false
  }}).then(function(viewer) {{
    molPlugin = viewer.plugin;
    var p = molPlugin;
    var chain = Promise.resolve();
    structs.forEach(function(s) {{
      chain = chain
        .then(function() {{ return p.builders.data.rawData({{data: s.data, label: s.label}}); }})
        .then(function(d) {{ return p.builders.structure.parseTrajectory(d, s.format); }})
        .then(function(t) {{ return p.builders.structure.hierarchy.applyPreset(t, 'default'); }});
    }});
    chain.then(function() {{
      var el = document.getElementById('molstar-loading');
      if (el) el.style.display = 'none';
    }}).catch(function(e) {{
      var el = document.getElementById('molstar-loading');
      if (el) el.style.display = 'none';
      document.getElementById('molstar-viewer').innerHTML =
        '<p style="color:#dc2626;padding:16px;">Mol* error: ' + e + '</p>';
    }});
  }}).catch(function(e) {{
    var el = document.getElementById('molstar-loading');
    if (el) el.style.display = 'none';
    document.getElementById('molstar-viewer').innerHTML =
      '<p style="color:#dc2626;padding:16px;">Mol* failed to initialize: ' + e + '</p>';
  }});
}})();"""

    html = f"""
  <div class="section-card mb-4">
    <div class="section-header">3D Structure Viewer (Mol*)</div>
    <div class="section-body">
      <p style="font-size:var(--font-sm);color:var(--text-secondary);margin-bottom:8px;">Toggle structures:</p>
      <div style="margin-bottom:14px;display:flex;flex-wrap:wrap;">{toggle_buttons}</div>
      <div style="position:relative;">
        <div id="molstar-loading" class="viewer-loading">
          <div class="spinner"></div> Loading 3D viewer&hellip;
        </div>
        <div id="molstar-viewer"
             style="width:100%;height:520px;position:relative;border-radius:var(--radius);
                    overflow:hidden;border:1px solid var(--border);"></div>
      </div>
      <p class="note mt-2">
        Receptor cartoon coloured by pLDDT B-factors (blue = high confidence).
        {pocket_note}
        Top docking pose (orange) from Uni-Mol Docking V2.
        Scroll to zoom &nbsp;·&nbsp; drag to rotate &nbsp;·&nbsp; right-click to pan.
      </p>
    </div>
  </div>"""

    return html, viewer_js


def _section_methods(input_summary, pocket_qc, docking_summary, model_id, timestamp):
    entities = input_summary.get("entities", [])
    protein  = next((e for e in entities if e.get("type") == "protein"), {})
    ligand   = next((e for e in entities if e.get("type") == "ligand"),  {})

    threshold      = pocket_qc.get("threshold", 70.0)
    sel_rank       = pocket_qc.get("selected_pocket_rank", "1")
    n_poses_req    = docking_summary.get("num_poses_requested", "—")
    grid           = docking_summary.get("grid", {})
    cx = _fmt(grid.get("center_x"), 2)
    cy = _fmt(grid.get("center_y"), 2)
    cz = _fmt(grid.get("center_z"), 2)
    box = _fmt(grid.get("size_x"), 1)

    return f"""
  <div class="section-card mb-4">
    <div class="section-header">Methods &amp; Caveats</div>
    <div class="section-body">

      <div class="row g-4">
        <div class="col-12 col-md-6">
          <h6 class="fw-semibold mb-2" style="color:var(--primary);">Node 01 — Boltz-2 Structure Prediction</h6>
          <table class="table table-sm table-borderless mb-0" style="font-size:var(--font-sm);">
            <tbody>
              <tr><td class="text-secondary pe-3">Mode</td>
                  <td>Holo (protein + ligand in YAML; <code>--output_format mmcif</code>)</td></tr>
              <tr><td class="text-secondary pe-3">Selected model</td>
                  <td><code>{model_id}</code> (highest overall confidence score)</td></tr>
              <tr><td class="text-secondary pe-3">Affinity</td>
                  <td>Boltz-2 affinity predicted for ligand chain {ligand.get("id", "B")}
                      (log<sub>10</sub>K<sub>D</sub>); for orientation — not validated for docking pose selection</td></tr>
            </tbody>
          </table>
        </div>

        <div class="col-12 col-md-6">
          <h6 class="fw-semibold mb-2" style="color:var(--primary);">Node 02 — P2Rank 2.4.2 Pocket Detection</h6>
          <table class="table table-sm table-borderless mb-0" style="font-size:var(--font-sm);">
            <tbody>
              <tr><td class="text-secondary pe-3">Profile</td>
                  <td><code>-c alphafold</code> (B-factor column = pLDDT; PrankWeb 4 protocol)</td></tr>
              <tr><td class="text-secondary pe-3">Input</td>
                  <td>Selected Boltz-2 mmCIF ({model_id})</td></tr>
            </tbody>
          </table>
        </div>

        <div class="col-12 col-md-6">
          <h6 class="fw-semibold mb-2" style="color:var(--primary);">Node 03 — Pocket QC + Grid</h6>
          <table class="table table-sm table-borderless mb-0" style="font-size:var(--font-sm);">
            <tbody>
              <tr><td class="text-secondary pe-3">pLDDT threshold</td>
                  <td>{threshold} (mean per-residue pLDDT of pocket residues)</td></tr>
              <tr><td class="text-secondary pe-3">Selected pocket</td>
                  <td>Rank {sel_rank}</td></tr>
              <tr><td class="text-secondary pe-3">Grid centre</td>
                  <td>({cx}, {cy}, {cz}) Å &nbsp;·&nbsp; box {box} Å</td></tr>
              <tr><td class="text-secondary pe-3">Receptor format</td>
                  <td>PDB (gemmi mmCIF→PDB; non-protein chains removed)</td></tr>
              <tr><td class="text-secondary pe-3">Ligand format</td>
                  <td>SDF (RDKit ETKDGv3 + MMFF from SMILES)</td></tr>
            </tbody>
          </table>
        </div>

        <div class="col-12 col-md-6">
          <h6 class="fw-semibold mb-2" style="color:var(--primary);">Node 04 — Uni-Mol Docking V2</h6>
          <table class="table table-sm table-borderless mb-0" style="font-size:var(--font-sm);">
            <tbody>
              <tr><td class="text-secondary pe-3">Engine</td>
                  <td>Uni-Mol Docking V2 (77.6% RMSD &lt; 2 Å on PoseBusters)</td></tr>
              <tr><td class="text-secondary pe-3">Poses requested</td>
                  <td>{n_poses_req}</td></tr>
              <tr><td class="text-secondary pe-3">Affinity</td>
                  <td>Not output — Uni-Mol outputs 3D poses only (no prmsd_score exposed)</td></tr>
            </tbody>
          </table>
        </div>

        <div class="col-12">
          <div class="alert alert-warning py-2 px-3 mb-2" style="font-size:var(--font-sm);">
            <strong>pLDDT caveat:</strong>
            pLDDT ≥ {threshold} is a <em>necessary but not sufficient</em> condition for reliable
            docking. High pLDDT filters disordered regions but does not guarantee correct
            side-chain geometry at the pocket.  Eguida &amp; Rognan (2023) found that 4 of the 5
            worst-performing docking targets had mean pLDDT ≥ 70 (JCIM, PMC9852548).
            Always inspect the selected pocket and consider experimental validation.
          </div>
          <p class="note mb-0">Report generated: {timestamp}</p>
        </div>
      </div>

    </div>
  </div>"""


# ---------------------------------------------------------------------------
# Full HTML assembly
# ---------------------------------------------------------------------------

def generate_html(
    model_id, input_summary, confidence_list, plddt, pae,
    pocket_qc, docking_summary, top_pose_sdf, timestamp, affinity
):
    entities      = input_summary.get("entities", [])
    protein_ent   = next((e for e in entities if e.get("type") == "protein"), {})
    target_name   = protein_ent.get("id", "A")

    # Build sections (some return (html, js) tuples for the viewer)
    sec_summary = _section_summary(input_summary, confidence_list, pocket_qc, docking_summary, plddt, affinity)

    # _section_confidence returns a single string with __PLDDT_JS__ / __PAE_JS__ markers
    conf_block  = _section_confidence(plddt, pae, confidence_list)
    # Extract the chart JS from the markers
    plddt_js, pae_js = "", ""
    if "__PLDDT_JS__" in conf_block:
        parts     = conf_block.split("__PLDDT_JS__")
        conf_html = parts[0]
        rest      = parts[1]
        if "__PAE_JS__" in rest:
            rest2    = rest.split("__PAE_JS__")
            plddt_js = rest2[0]
            pae_js   = rest2[1]
        else:
            plddt_js = rest
    else:
        conf_html = conf_block

    sec_pockets = _section_pocket_table(pocket_qc)
    viewer_html, viewer_js = _section_viewer(plddt, pocket_qc, top_pose_sdf)
    sec_methods = _section_methods(input_summary, pocket_qc, docking_summary, model_id, timestamp)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Structure → Pocket → Docking Pipeline Report</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/molstar@5.9.0/build/viewer/molstar.css">
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/molstar@5.9.0/build/viewer/molstar.js"></script>
  <style>
    :root {{
      --primary: #075985;
      --secondary: #0284c7;
      --text: #1e293b;
      --text-secondary: #475569;
      --text-muted: #64748b;
      --text-faint: #94a3b8;
      --bg-body: #f8fafc;
      --bg-card: #ffffff;
      --bg-section-header: #f8fafc;
      --border: #e2e8f0;
      --success: #16a34a;
      --warning: #d97706;
      --danger: #dc2626;
      --radius: 12px;
      --radius-pill: 20px;
      --shadow: 0 1px 6px rgba(0,0,0,0.06);
      --font-xs: 0.75rem;
      --font-sm: 0.84rem;
    }}
    body {{ background: var(--bg-body); color: var(--text); font-family: system-ui, sans-serif; }}
    .hero {{
      background: linear-gradient(135deg, var(--primary) 0%, var(--secondary) 100%);
      color: white; border-radius: var(--radius); padding: 32px 36px; margin-bottom: 28px;
    }}
    .hero h1 {{ font-size: 1.8rem; font-weight: 700; margin: 0 0 6px; }}
    .hero .sub {{ font-size: 0.9rem; opacity: 0.85; margin: 0; }}
    .stat-card {{
      background: var(--bg-card); border-radius: var(--radius); padding: 20px 16px;
      text-align: center; box-shadow: var(--shadow); height: 100%;
    }}
    .stat-card .value {{ font-size: 1.9rem; font-weight: 700; line-height: 1.1;
                         color: var(--primary); }}
    .stat-card .label {{ font-size: var(--font-xs); text-transform: uppercase;
                         letter-spacing: .7px; color: var(--text-muted); margin-top: 4px; }}
    .section-card {{
      background: var(--bg-card); border-radius: var(--radius);
      box-shadow: var(--shadow); overflow: hidden;
    }}
    .section-header {{
      padding: 14px 24px; font-weight: 600; font-size: 0.97rem; color: var(--text);
      border-bottom: 1px solid var(--border); background: var(--bg-section-header);
    }}
    .section-body {{ padding: 20px 24px; }}
    .table th {{ font-size: var(--font-sm); color: var(--text-secondary); font-weight: 600; cursor: pointer; }}
    .table td {{ font-size: var(--font-sm); vertical-align: middle; }}
    .note {{ font-size: var(--font-xs); color: var(--text-faint); }}
    .pill {{
      display: inline-block; padding: 2px 10px; border-radius: var(--radius-pill);
      font-size: var(--font-xs); font-weight: 600; color: white;
    }}
    #pocket-table tbody tr:nth-child(even) {{ background: var(--bg-body); }}
    #pocket-table tbody tr:hover {{ background: #f0f7ff; }}
    .viewer-loading {{
      position: absolute; inset: 0; display: flex; align-items: center;
      justify-content: center; background: var(--bg-body); z-index: 10;
      font-size: var(--font-sm); color: var(--text-muted);
    }}
    .viewer-loading .spinner {{
      width: 24px; height: 24px; border: 3px solid var(--border);
      border-top-color: var(--primary); border-radius: 50%;
      animation: spin 0.8s linear infinite; margin-right: 10px;
    }}
    @keyframes spin {{ to {{ transform: rotate(360deg); }} }}
    @media print {{
      body {{ background: white; }}
      .hero {{ background: var(--primary) !important; -webkit-print-color-adjust: exact; print-color-adjust: exact; }}
      .section-card {{ box-shadow: none; border: 1px solid var(--border); break-inside: avoid; }}
      #molstar-viewer, .viewer-loading {{ display: none !important; }}
      [onclick*="toggleStruct"] {{ display: none !important; }}
    }}
  </style>
</head>
<body>
<div class="container-fluid py-4" style="max-width:1280px; margin:0 auto;">

  <div class="hero">
    <h1>Structure &#8594; Pocket &#8594; Docking Pipeline Report</h1>
    <p class="sub">
      Boltz-2 &nbsp;&#8594;&nbsp; P2Rank 2.4.2 &nbsp;&#8594;&nbsp; Uni-Mol Docking V2
      &nbsp;&middot;&nbsp; Generated: {timestamp}
    </p>
  </div>

  {sec_summary}
  {conf_html}
  {sec_pockets}
  {viewer_html}
  {sec_methods}

</div>

<script>
// ── Plotly: pLDDT histogram ──
{plddt_js}

// ── Plotly: PAE heatmap ──
{pae_js}

// ── Mol* viewer ──
{viewer_js}

// ── Sortable pocket table ──
document.querySelectorAll('#pocket-table thead th').forEach(function(th) {{
  th.addEventListener('click', function() {{
    var tbody = th.closest('table').querySelector('tbody');
    var col = th.cellIndex;
    var asc = th.dataset.dir !== 'asc';
    th.dataset.dir = asc ? 'asc' : 'desc';
    var rows = Array.from(tbody.rows);
    rows.sort(function(a, b) {{
      var av = a.cells[col].textContent.replace(/[^0-9.\-]/g, '');
      var bv = b.cells[col].textContent.replace(/[^0-9.\-]/g, '');
      av = av !== '' ? parseFloat(av) : a.cells[col].textContent;
      bv = bv !== '' ? parseFloat(bv) : b.cells[col].textContent;
      if (typeof av === 'string') return asc ? av.localeCompare(bv) : bv.localeCompare(av);
      return asc ? av - bv : bv - av;
    }});
    rows.forEach(function(r) {{ tbody.appendChild(r); }});
  }});
}});
</script>
</body>
</html>"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    print("Starting Node 05: Generate Report", flush=True)

    model_id       = load_model_id()
    print(f"  Selected model: {model_id}", flush=True)

    plddt          = load_plddt(model_id)
    print(f"  pLDDT: {len(plddt)} residues, "
          f"range [{plddt.min():.1f}, {plddt.max():.1f}]" if len(plddt) > 0 else "  pLDDT: not found",
          flush=True)

    pae            = load_pae(model_id)
    print(f"  PAE: {pae.shape if pae is not None else 'not found'}", flush=True)

    confidence_list = load_confidence()
    print(f"  Confidence files: {len(confidence_list)}", flush=True)

    affinity       = load_affinity()
    print(f"  Affinity: {'found' if affinity else 'not found'}", flush=True)

    input_summary  = load_input_summary()
    pocket_qc      = load_pocket_qc()
    docking_summary = load_docking_summary()
    top_pose_sdf   = load_top_pose_sdf()
    print(f"  Top pose SDF: {'found' if top_pose_sdf else 'not found'}", flush=True)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M UTC")
    print("Generating HTML ...", flush=True)

    html = generate_html(
        model_id=model_id,
        input_summary=input_summary,
        confidence_list=confidence_list,
        plddt=plddt,
        pae=pae,
        pocket_qc=pocket_qc,
        docking_summary=docking_summary,
        top_pose_sdf=top_pose_sdf,
        timestamp=timestamp,
        affinity=affinity,
    )

    os.makedirs("./outputs", exist_ok=True)
    Path("./outputs/report.html").write_text(html)
    print(f"Wrote outputs/report.html  ({len(html):,} chars)", flush=True)
    print("Node 05 complete.", flush=True)


if __name__ == "__main__":
    main()
