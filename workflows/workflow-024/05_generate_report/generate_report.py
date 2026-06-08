#!/usr/bin/env python3
"""
Node 05: Generate Comparison Report

1. Load and merge Vina + GNINA score CSVs
2. Compute ranking statistics (Spearman ρ, Kendall τ, top-N overlap)
3. Convert top-5 poses to MOL2 via OpenBabel
4. Write comparative_screening_metrics.csv
5. Generate self-contained docking_performance_report.html (Bootstrap 5 + Plotly.js CDN)
"""

import csv
import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path


# ── Data loading and merging ──────────────────────────────────────────────────

def load_csv(path):
    """Load CSV as list of dicts; return [] if file missing."""
    p = Path(path)
    if not p.exists():
        print(f"  WARNING: {path} not found", flush=True)
        return []
    with open(p, newline="") as f:
        return list(csv.DictReader(f))


def safe_float(val):
    try:
        return float(val) if val not in (None, "") else None
    except (ValueError, TypeError):
        return None


def merge_scores(vina_rows, gnina_rows):
    """Merge on compound name; assign combined ranks."""
    gnina_by_name = {r["name"]: r for r in gnina_rows}
    vina_by_name  = {r["name"]: r for r in vina_rows}
    all_names = sorted(set(vina_by_name) | set(gnina_by_name))

    rows = []
    for name in all_names:
        v = vina_by_name.get(name, {})
        g = gnina_by_name.get(name, {})
        rows.append({
            "name":                         name,
            "vina_rank":                    safe_float(v.get("rank")),
            "vina_affinity_kcal_mol":       safe_float(v.get("best_affinity_kcal_mol")),
            "vina_status":                  v.get("status", "missing"),
            "gnina_rank":                   safe_float(g.get("rank")),
            "gnina_cnn_score":              safe_float(g.get("cnn_score")),
            "gnina_cnn_affinity_kcal_mol":  safe_float(g.get("cnn_affinity_kcal_mol")),
            "gnina_vina_affinity_kcal_mol": safe_float(g.get("vina_affinity_kcal_mol")),
            "gnina_status":                 g.get("status", "missing"),
        })

    # rank_difference: positive = Vina ranked it higher (lower number) than GNINA
    for r in rows:
        if r["vina_rank"] is not None and r["gnina_rank"] is not None:
            r["rank_difference"] = int(r["vina_rank"]) - int(r["gnina_rank"])
        else:
            r["rank_difference"] = None

    return rows


# ── Statistics ────────────────────────────────────────────────────────────────

def compute_statistics(merged):
    """Return dict of ranking agreement metrics."""
    both = [r for r in merged
            if r["vina_rank"] is not None and r["gnina_rank"] is not None]
    if len(both) < 2:
        return {"n_both": len(both), "error": "not enough matched compounds"}

    vina_ranks  = [r["vina_rank"]  for r in both]
    gnina_ranks = [r["gnina_rank"] for r in both]
    vina_aff    = [r["vina_affinity_kcal_mol"]      for r in both]
    gnina_cnn   = [r["gnina_cnn_score"]              for r in both]
    gnina_aff   = [r["gnina_cnn_affinity_kcal_mol"]  for r in both]

    def _rank_list(xs):
        sorted_idx = sorted(range(len(xs)), key=lambda i: xs[i])
        ranks = [0] * len(xs)
        for r, i in enumerate(sorted_idx, 1):
            ranks[i] = r
        return ranks

    def _spearmanr(x, y):
        rx, ry = _rank_list(x), _rank_list(y)
        n = len(x)
        d2 = sum((a - b) ** 2 for a, b in zip(rx, ry))
        return 1 - 6 * d2 / (n * (n ** 2 - 1))

    def _kendalltau(x, y):
        n = len(x)
        conc = disc = 0
        for i in range(n):
            for j in range(i + 1, n):
                s = (x[i] - x[j]) * (y[i] - y[j])
                if s > 0:
                    conc += 1
                elif s < 0:
                    disc += 1
        return (conc - disc) / (n * (n - 1) / 2)

    try:
        from scipy.stats import spearmanr as _sp, kendalltau as _kt
        def _spearmanr(x, y): return _sp(x, y)[0]
        def _kendalltau(x, y): return _kt(x, y)[0]
    except ImportError:
        print("  INFO: scipy not available; using pure-Python rank correlations", flush=True)

    rho_rank  = _spearmanr(vina_ranks, gnina_ranks)
    tau_rank  = _kendalltau(vina_ranks, gnina_ranks)
    aff_pairs = [(v, g) for v, g in zip(vina_aff, gnina_aff) if v is not None and g is not None]
    rho_aff   = _spearmanr([p[0] for p in aff_pairs], [p[1] for p in aff_pairs]) if aff_pairs else None
    cnn_pairs = [(v, g) for v, g in zip(vina_aff, gnina_cnn) if v is not None and g is not None]
    rho_score = _spearmanr([p[0] for p in cnn_pairs], [p[1] for p in cnn_pairs]) if cnn_pairs else None

    def topn_overlap(n):
        vt = {r["name"] for r in sorted(both, key=lambda r: r["vina_rank"])[:n]}
        gt = {r["name"] for r in sorted(both, key=lambda r: r["gnina_rank"])[:n]}
        return len(vt & gt)

    n = len(both)
    return {
        "n_compounds":          n,
        "spearman_rho_rank":    round(rho_rank,  4) if rho_rank  is not None else None,
        "spearman_p_rank":      None,
        "kendall_tau_rank":     round(tau_rank,  4) if tau_rank  is not None else None,
        "kendall_p_rank":       None,
        "spearman_rho_affinity": round(rho_aff,  4) if rho_aff   is not None else None,
        "spearman_rho_vina_cnn": round(rho_score,4) if rho_score is not None else None,
        "top3_overlap":  topn_overlap(min(3, n)),
        "top5_overlap":  topn_overlap(min(5, n)),
    }


# ── MOL2 conversion of top poses ─────────────────────────────────────────────

def convert_top_poses_mol2(vina_pdbqt, gnina_sdf, n=5):
    """Convert top-N poses from each tool to MOL2 via obabel."""
    mol2_dir = Path("top_poses_mol2")
    mol2_dir.mkdir(exist_ok=True)
    results = {"vina": [], "gnina": []}

    # Vina: split PDBQT into individual molecule files, convert each to MOL2
    if Path(vina_pdbqt).exists():
        content = Path(vina_pdbqt).read_text()
        # Split on REMARK Name markers (one block per compound's top pose)
        blocks = re.split(r'(?=REMARK Name = )', content)
        for i, block in enumerate(blocks[:n]):
            if not block.strip():
                continue
            m = re.search(r'REMARK Name = (.+)', block)
            name = re.sub(r'[^\w.-]', '_', m.group(1).strip()) if m else f"vina_{i+1}"
            tmp_pdbqt = mol2_dir / f"vina_{i+1}_{name}.pdbqt"
            out_mol2  = mol2_dir / f"vina_{i+1}_{name}.mol2"
            tmp_pdbqt.write_text(block)
            res = subprocess.run(
                ["obabel", str(tmp_pdbqt), "-O", str(out_mol2)],
                capture_output=True, text=True
            )
            if res.returncode == 0 and out_mol2.exists():
                results["vina"].append(str(out_mol2))
                print(f"  Vina top-{i+1}: {out_mol2.name}", flush=True)
            tmp_pdbqt.unlink(missing_ok=True)

    # GNINA: split SDF into individual molecules, convert each to MOL2
    if Path(gnina_sdf).exists():
        content = Path(gnina_sdf).read_text()
        records = re.split(r'(?<=\$\$\$\$)\n', content)
        for i, record in enumerate(records[:n]):
            if "$$$$" not in record:
                continue
            tmp_sdf  = mol2_dir / f"gnina_{i+1}.sdf"
            out_mol2 = mol2_dir / f"gnina_{i+1}.mol2"
            tmp_sdf.write_text(record)
            res = subprocess.run(
                ["obabel", str(tmp_sdf), "-O", str(out_mol2)],
                capture_output=True, text=True
            )
            if res.returncode == 0 and out_mol2.exists():
                results["gnina"].append(str(out_mol2))
                print(f"  GNINA top-{i+1}: {out_mol2.name}", flush=True)
            tmp_sdf.unlink(missing_ok=True)

    return results


# ── HTML report generation ────────────────────────────────────────────────────

def _fmt(v, decimals=3):
    if v is None:
        return "—"
    return f"{v:.{decimals}f}"


def generate_html(merged, stats, vina_meta, gnina_meta, timestamp):
    both = [r for r in merged
            if r["vina_rank"] is not None and r["gnina_rank"] is not None]
    both_sorted_vina  = sorted(both, key=lambda r: r["vina_rank"])
    both_sorted_gnina = sorted(both, key=lambda r: r["gnina_rank"])

    n_screened  = stats.get("n_compounds", len(merged))
    rho_rank    = stats.get("spearman_rho_rank")
    tau_rank    = stats.get("kendall_tau_rank")
    top5_ol     = stats.get("top5_overlap", "—")
    rho_aff     = stats.get("spearman_rho_affinity")

    # Plotly data — rank agreement scatter
    names_js   = json.dumps([r["name"] for r in both])
    vina_rk_js = json.dumps([r["vina_rank"]  for r in both])
    gnina_rk_js= json.dumps([r["gnina_rank"] for r in both])

    # Score correlation scatter (Vina ΔG vs GNINA CNNscore)
    corr_both = [r for r in both
                 if r["vina_affinity_kcal_mol"] is not None and r["gnina_cnn_score"] is not None]
    corr_x_js = json.dumps([r["vina_affinity_kcal_mol"] for r in corr_both])
    corr_y_js = json.dumps([r["gnina_cnn_score"]         for r in corr_both])
    corr_n_js = json.dumps([r["name"]                    for r in corr_both])

    # Affinity correlation scatter (Vina ΔG vs GNINA CNN affinity)
    aff_both = [r for r in both
                if r["vina_affinity_kcal_mol"] is not None and r["gnina_cnn_affinity_kcal_mol"] is not None]
    aff_x_js = json.dumps([r["vina_affinity_kcal_mol"]     for r in aff_both])
    aff_y_js = json.dumps([r["gnina_cnn_affinity_kcal_mol"] for r in aff_both])
    aff_n_js = json.dumps([r["name"]                        for r in aff_both])

    # Score distribution histograms
    vina_scores_js  = json.dumps([r["vina_affinity_kcal_mol"]
                                   for r in merged if r["vina_affinity_kcal_mol"] is not None])
    gnina_cnn_js    = json.dumps([r["gnina_cnn_score"]
                                   for r in merged if r["gnina_cnn_score"] is not None])
    gnina_aff_js    = json.dumps([r["gnina_cnn_affinity_kcal_mol"]
                                   for r in merged if r["gnina_cnn_affinity_kcal_mol"] is not None])

    # Top-5 hits tables
    def _top5_rows(rows, tool_ranks_by_vina=True):
        html = ""
        for r in rows[:5]:
            vina_rank  = int(r["vina_rank"])  if r.get("vina_rank")  else "—"
            gnina_rank = int(r["gnina_rank"]) if r.get("gnina_rank") else "—"
            html += (
                f'<tr><td>{r["name"]}</td>'
                f'<td>{vina_rank}</td>'
                f'<td>{_fmt(r["vina_affinity_kcal_mol"])} kcal/mol</td>'
                f'<td>{gnina_rank}</td>'
                f'<td>{_fmt(r["gnina_cnn_score"])}</td>'
                f'<td>{_fmt(r["gnina_cnn_affinity_kcal_mol"])} kcal/mol</td>'
                f'<td>{_fmt(r.get("rank_difference"), 0)}</td></tr>\n'
            )
        return html

    top5_vina_rows  = _top5_rows(both_sorted_vina)
    top5_gnina_rows = _top5_rows(both_sorted_gnina)

    # Full comparison table rows
    def _full_table_rows():
        html = ""
        for r in sorted(merged, key=lambda x: (x["vina_rank"] or 999)):
            rd = r.get("rank_difference")
            rd_str = f'+{rd}' if rd and rd > 0 else (str(rd) if rd is not None else "—")
            html += (
                f'<tr>'
                f'<td>{r["name"]}</td>'
                f'<td>{_fmt(r["vina_rank"], 0)}</td>'
                f'<td>{_fmt(r["vina_affinity_kcal_mol"])}</td>'
                f'<td>{_fmt(r["gnina_rank"], 0)}</td>'
                f'<td>{_fmt(r["gnina_cnn_score"])}</td>'
                f'<td>{_fmt(r["gnina_cnn_affinity_kcal_mol"])}</td>'
                f'<td>{rd_str}</td>'
                f'</tr>\n'
            )
        return html

    full_table_rows = _full_table_rows()

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>Vina vs GNINA Docking Comparison</title>
  <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css" rel="stylesheet">
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    body {{ background: #f1f5f9; color: #1e293b; }}
    .hero {{
      background: linear-gradient(135deg, #1d4ed8 0%, #7c3aed 100%);
      color: white; border-radius: 16px; padding: 32px 36px; margin-bottom: 28px;
    }}
    .hero h1 {{ font-size: 1.9rem; font-weight: 700; margin: 0 0 6px; }}
    .hero .sub {{ font-size: 0.92rem; opacity: 0.85; margin: 0; }}
    .stat-card {{
      background: white; border-radius: 14px; padding: 22px 20px; text-align: center;
      box-shadow: 0 2px 12px rgba(0,0,0,0.07); height: 100%;
    }}
    .stat-card .value {{ font-size: 2rem; font-weight: 700; line-height: 1.1; }}
    .stat-card .label {{ font-size: 0.78rem; text-transform: uppercase;
                         letter-spacing: .8px; color: #64748b; margin-top: 4px; }}
    .section-card {{
      background: white; border-radius: 14px;
      box-shadow: 0 2px 12px rgba(0,0,0,0.07); margin-bottom: 24px; overflow: hidden;
    }}
    .section-header {{
      padding: 16px 24px; font-weight: 600; font-size: 1rem; color: #1e293b;
      border-bottom: 1px solid #f1f5f9;
    }}
    .section-body {{ padding: 20px 24px; }}
    .vina-badge  {{ background:#3b82f6; color:white; padding:3px 10px;
                    border-radius:12px; font-size:11px; }}
    .gnina-badge {{ background:#8b5cf6; color:white; padding:3px 10px;
                    border-radius:12px; font-size:11px; }}
    .table th {{ font-size: 0.82rem; color: #475569; font-weight: 600; }}
    .table td {{ font-size: 0.85rem; vertical-align: middle; }}
    .note {{ font-size: 0.78rem; color: #94a3b8; margin-top: 8px; }}
  </style>
</head>
<body>
<div class="container-fluid py-4" style="max-width:1300px; margin:0 auto;">

  <!-- Hero header -->
  <div class="hero">
    <h1>&#9878; Vina vs GNINA: Docking Comparison</h1>
    <p class="sub">
      Target: Carbonic Anhydrase II (PDB: 1OKL) &nbsp;·&nbsp;
      {n_screened} compounds screened &nbsp;·&nbsp;
      Generated: {timestamp}
    </p>
  </div>

  <!-- Summary stats cards -->
  <div class="row g-3 mb-4">
    <div class="col-6 col-md-3">
      <div class="stat-card">
        <div class="value" style="color:#1d4ed8;">{n_screened}</div>
        <div class="label">Compounds Screened</div>
      </div>
    </div>
    <div class="col-6 col-md-3">
      <div class="stat-card">
        <div class="value" style="color:{'#16a34a' if rho_rank and rho_rank > 0.6 else '#f59e0b' if rho_rank and rho_rank > 0.3 else '#dc2626'}">
          {_fmt(rho_rank)}
        </div>
        <div class="label">Spearman &#961; (ranks)</div>
      </div>
    </div>
    <div class="col-6 col-md-3">
      <div class="stat-card">
        <div class="value" style="color:#7c3aed;">{_fmt(tau_rank)}</div>
        <div class="label">Kendall &#964; (ranks)</div>
      </div>
    </div>
    <div class="col-6 col-md-3">
      <div class="stat-card">
        <div class="value" style="color:#0891b2;">{top5_ol} / {min(5, n_screened)}</div>
        <div class="label">Top-5 Overlap</div>
      </div>
    </div>
  </div>

  <!-- Score correlation scatter -->
  <div class="section-card">
    <div class="section-header">&#128202; Score Correlation: Vina &#916;G vs GNINA CNNscore</div>
    <div class="section-body">
      <div id="corr-chart" style="height:380px;"></div>
      <p class="note">
        Vina &#916;G (kcal/mol, more negative = stronger predicted binding) vs
        GNINA CNNscore (0–1, higher = CNN-predicted binder). Compounds in the
        upper-left are ranked highly by GNINA but weakly by Vina — suggesting the
        CNN detects metal-coordination geometry that Vina&#8217;s additive terms miss.
      </p>
    </div>
  </div>

  <!-- Rank agreement scatter -->
  <div class="section-card">
    <div class="section-header">&#128257; Ranking Agreement: Vina Rank vs GNINA Rank</div>
    <div class="section-body">
      <div id="rank-chart" style="height:380px;"></div>
      <p class="note">
        Points on the diagonal (y&#8239;=&#8239;x) indicate both tools agree on the rank.
        Points above the diagonal are ranked higher by Vina; below by GNINA.
        Spearman &#961; = {_fmt(rho_rank)} &nbsp;·&nbsp; Kendall &#964; = {_fmt(tau_rank)}.
      </p>
    </div>
  </div>

  <!-- Score distributions -->
  <div class="section-card">
    <div class="section-header">&#128200; Score Distributions</div>
    <div class="section-body row g-3">
      <div class="col-12 col-md-6">
        <div id="vina-hist" style="height:300px;"></div>
        <p class="note">Vina &#916;G distribution (kcal/mol)</p>
      </div>
      <div class="col-12 col-md-6">
        <div id="gnina-hist" style="height:300px;"></div>
        <p class="note">GNINA CNNscore distribution (0–1)</p>
      </div>
    </div>
  </div>

  <!-- Affinity correlation -->
  <div class="section-card">
    <div class="section-header">&#9881; Affinity Correlation: Vina &#916;G vs GNINA CNN Affinity</div>
    <div class="section-body">
      <div id="aff-chart" style="height:340px;"></div>
      <p class="note">
        Both axes report kcal/mol; more negative = stronger predicted binding.
        Spearman &#961; = {_fmt(rho_aff)}.
        Divergence here reflects the scoring function, not pose generation
        (both tools share the same MCMC sampler).
      </p>
    </div>
  </div>

  <!-- Top-5 tables -->
  <div class="row g-3 mb-4">
    <div class="col-12 col-md-6">
      <div class="section-card h-100">
        <div class="section-header">
          <span class="vina-badge">Vina</span> Top-5 by &#916;G Affinity
        </div>
        <div class="section-body p-0">
          <table class="table table-hover mb-0">
            <thead><tr>
              <th>Name</th><th>Vina Rank</th><th>Vina &#916;G</th>
              <th>GNINA Rank</th><th>CNN Score</th><th>CNN Aff.</th><th>&#916;Rank</th>
            </tr></thead>
            <tbody>{top5_vina_rows}</tbody>
          </table>
        </div>
      </div>
    </div>
    <div class="col-12 col-md-6">
      <div class="section-card h-100">
        <div class="section-header">
          <span class="gnina-badge">GNINA</span> Top-5 by CNNscore
        </div>
        <div class="section-body p-0">
          <table class="table table-hover mb-0">
            <thead><tr>
              <th>Name</th><th>Vina Rank</th><th>Vina &#916;G</th>
              <th>GNINA Rank</th><th>CNN Score</th><th>CNN Aff.</th><th>&#916;Rank</th>
            </tr></thead>
            <tbody>{top5_gnina_rows}</tbody>
          </table>
        </div>
      </div>
    </div>
  </div>

  <!-- Full comparison table -->
  <div class="section-card">
    <div class="section-header">&#128203; Full Comparison Table</div>
    <div class="section-body p-0">
      <div class="table-responsive">
        <table class="table table-sm table-striped table-hover mb-0" id="full-table">
          <thead>
            <tr>
              <th>Name</th>
              <th>Vina Rank &#8597;</th><th>Vina &#916;G (kcal/mol)</th>
              <th>GNINA Rank &#8597;</th><th>CNN Score</th><th>CNN Aff. (kcal/mol)</th>
              <th title="Vina rank − GNINA rank (+ = Vina ranks it higher)">&#916;Rank</th>
            </tr>
          </thead>
          <tbody>{full_table_rows}</tbody>
        </table>
      </div>
    </div>
  </div>

  <!-- Methods -->
  <div class="section-card">
    <div class="section-header">&#128221; Methods</div>
    <div class="section-body">
      <div class="row g-4">

        <div class="col-12 col-md-6">
          <h6 class="fw-semibold mb-2" style="color:#1d4ed8;">
            <span class="vina-badge me-2">Vina</span> AutoDock Vina
          </h6>
          <table class="table table-sm table-borderless mb-0" style="font-size:0.84rem;">
            <tbody>
              <tr><td class="text-secondary pe-3" style="white-space:nowrap;">Scoring function</td>
                  <td>Empirical energy function (Vina weights)</td></tr>
              <tr><td class="text-secondary pe-3">Exhaustiveness</td>
                  <td>{vina_meta.get("exhaustiveness", 8)}</td></tr>
              <tr><td class="text-secondary pe-3">Num. output modes</td>
                  <td>{vina_meta.get("num_modes", 9)}</td></tr>
              <tr><td class="text-secondary pe-3">Hardware</td>
                  <td>CPU-only</td></tr>
              <tr><td class="text-secondary pe-3">Output format</td>
                  <td>PDBQT</td></tr>
            </tbody>
          </table>
        </div>

        <div class="col-12 col-md-6">
          <h6 class="fw-semibold mb-2" style="color:#7c3aed;">
            <span class="gnina-badge me-2">GNINA</span> GNINA
          </h6>
          <table class="table table-sm table-borderless mb-0" style="font-size:0.84rem;">
            <tbody>
              <tr><td class="text-secondary pe-3" style="white-space:nowrap;">CNN scoring mode</td>
                  <td><code>--cnn_scoring rescore</code> — Vina poses rescored by CNN; CNN score used for ranking</td></tr>
              <tr><td class="text-secondary pe-3">Exhaustiveness</td>
                  <td>{gnina_meta.get("exhaustiveness", 8)}</td></tr>
              <tr><td class="text-secondary pe-3">Num. output modes</td>
                  <td>{gnina_meta.get("num_modes", 9)}</td></tr>
              <tr><td class="text-secondary pe-3">Hardware</td>
                  <td>CPU-only (<code>--no_gpu</code>)</td></tr>
              <tr><td class="text-secondary pe-3">Output format</td>
                  <td>SDF</td></tr>
            </tbody>
          </table>
        </div>

        <div class="col-12">
          <h6 class="fw-semibold mb-2" style="color:#0f172a;">Shared Experimental Settings</h6>
          <table class="table table-sm table-borderless mb-0" style="font-size:0.84rem;">
            <tbody>
              <tr><td class="text-secondary pe-3" style="width:200px;">Target protein</td>
                  <td>Carbonic Anhydrase II (CA-II), PDB: 1OKL</td></tr>
              <tr><td class="text-secondary pe-3">Compounds screened</td>
                  <td>{n_screened}</td></tr>
              <tr><td class="text-secondary pe-3">Pose format (downstream)</td>
                  <td>Top-5 poses per tool converted to MOL2 via OpenBabel</td></tr>
              <tr><td class="text-secondary pe-3">Ranking metric — Vina</td>
                  <td>&#916;G (kcal/mol); lower (more negative) = stronger predicted binding</td></tr>
              <tr><td class="text-secondary pe-3">Ranking metric — GNINA</td>
                  <td>CNNscore (0&#8211;1); higher = CNN-predicted binder; CNN affinity (kcal/mol) also reported</td></tr>
              <tr><td class="text-secondary pe-3">Report generated</td>
                  <td>{timestamp}</td></tr>
            </tbody>
          </table>
        </div>

      </div>
    </div>
  </div>

</div>

<script>
// ── Score correlation scatter ──
Plotly.newPlot('corr-chart', [{{
  type: 'scatter', mode: 'markers',
  x: {corr_x_js}, y: {corr_y_js}, text: {corr_n_js},
  textposition: 'top center', textfont: {{ size: 10 }},
  marker: {{ size: 12, color: {corr_x_js},
             colorscale: 'RdBu', reversescale: true,
             showscale: true, colorbar: {{ title: 'Vina ΔG' }} }},
  hovertemplate: '<b>%{{text}}</b><br>Vina ΔG: %{{x:.2f}} kcal/mol<br>GNINA CNN: %{{y:.3f}}<extra></extra>'
}}], {{
  xaxis: {{ title: 'Vina ΔG (kcal/mol)' }},
  yaxis: {{ title: 'GNINA CNNscore' }},
  template: 'plotly_white', height: 380,
  margin: {{ t: 20, b: 50, l: 60, r: 20 }}
}}, {{responsive: true}});

// ── Rank agreement scatter ──
var maxRank = Math.max(...{vina_rk_js}.concat({gnina_rk_js}));
Plotly.newPlot('rank-chart', [
  {{
    type: 'scatter', mode: 'lines',
    x: [1, maxRank], y: [1, maxRank],
    line: {{ color: '#cbd5e1', dash: 'dot', width: 1.5 }},
    hoverinfo: 'none', showlegend: false
  }},
  {{
    type: 'scatter', mode: 'markers',
    x: {vina_rk_js}, y: {gnina_rk_js}, text: {names_js},
    textposition: 'top center', textfont: {{ size: 10 }},
    marker: {{ size: 12, color: '#6366f1' }},
    hovertemplate: '<b>%{{text}}</b><br>Vina rank: %{{x}}<br>GNINA rank: %{{y}}<extra></extra>',
    name: 'Compounds'
  }}
], {{
  xaxis: {{ title: 'Vina Rank (1 = best)', dtick: 1 }},
  yaxis: {{ title: 'GNINA Rank (1 = best)', dtick: 1 }},
  template: 'plotly_white', height: 380,
  margin: {{ t: 20, b: 50, l: 60, r: 20 }}
}}, {{responsive: true}});

// ── Vina ΔG histogram ──
Plotly.newPlot('vina-hist', [{{
  type: 'histogram', x: {vina_scores_js},
  marker: {{ color: '#3b82f6', opacity: 0.8 }},
  xbins: {{ size: 0.5 }}, name: 'Vina ΔG'
}}], {{
  xaxis: {{ title: 'Vina ΔG (kcal/mol)' }},
  yaxis: {{ title: 'Count' }},
  template: 'plotly_white', height: 300,
  margin: {{ t: 20, b: 50, l: 50, r: 20 }},
  showlegend: false
}}, {{responsive: true}});

// ── GNINA CNNscore histogram ──
Plotly.newPlot('gnina-hist', [{{
  type: 'histogram', x: {gnina_cnn_js},
  marker: {{ color: '#8b5cf6', opacity: 0.8 }},
  xbins: {{ size: 0.05 }}, name: 'GNINA CNNscore'
}}], {{
  xaxis: {{ title: 'GNINA CNNscore (0–1)', range: [0, 1] }},
  yaxis: {{ title: 'Count' }},
  template: 'plotly_white', height: 300,
  margin: {{ t: 20, b: 50, l: 50, r: 20 }},
  showlegend: false
}}, {{responsive: true}});

// ── Affinity correlation scatter ──
Plotly.newPlot('aff-chart', [{{
  type: 'scatter', mode: 'markers',
  x: {aff_x_js}, y: {aff_y_js}, text: {aff_n_js},
  textposition: 'top center', textfont: {{ size: 10 }},
  marker: {{ size: 12, color: '#0891b2' }},
  hovertemplate: '<b>%{{text}}</b><br>Vina: %{{x:.2f}}<br>GNINA CNN aff: %{{y:.2f}}<extra></extra>'
}}], {{
  xaxis: {{ title: 'Vina ΔG (kcal/mol)' }},
  yaxis: {{ title: 'GNINA CNN Affinity (kcal/mol)' }},
  template: 'plotly_white', height: 340,
  margin: {{ t: 20, b: 50, l: 70, r: 20 }}
}}, {{responsive: true}});

// ── Sortable full table ──
document.querySelectorAll('#full-table thead th').forEach(function(th) {{
  th.style.cursor = 'pointer';
  th.addEventListener('click', function() {{
    var tbody = th.closest('table').querySelector('tbody');
    var col = th.cellIndex;
    var asc = th.dataset.dir !== 'asc';
    th.dataset.dir = asc ? 'asc' : 'desc';
    var rows = Array.from(tbody.rows);
    rows.sort(function(a, b) {{
      var av = parseFloat(a.cells[col].textContent.replace(/[^0-9.-]/g,''));
      var bv = parseFloat(b.cells[col].textContent.replace(/[^0-9.-]/g,''));
      if (isNaN(av)) av = a.cells[col].textContent;
      if (isNaN(bv)) bv = b.cells[col].textContent;
      return asc ? (av > bv ? 1 : -1) : (av < bv ? 1 : -1);
    }});
    rows.forEach(function(r) {{ tbody.appendChild(r); }});
  }});
}});
</script>
</body>
</html>"""


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    print("Loading scores ...", flush=True)
    vina_rows  = load_csv("vina_screening_scores.csv")
    gnina_rows = load_csv("gnina_screening_scores.csv")

    if not vina_rows and not gnina_rows:
        print("ERROR: no score files found", flush=True)
        sys.exit(1)

    print(f"  Vina:  {len(vina_rows)} compounds", flush=True)
    print(f"  GNINA: {len(gnina_rows)} compounds", flush=True)

    merged = merge_scores(vina_rows, gnina_rows)
    print(f"  Merged: {len(merged)} unique compounds", flush=True)

    print("\nComputing statistics ...", flush=True)
    stats = compute_statistics(merged)
    print(f"  Spearman ρ (rank): {stats.get('spearman_rho_rank')}", flush=True)
    print(f"  Kendall τ (rank):  {stats.get('kendall_tau_rank')}", flush=True)
    print(f"  Top-5 overlap:     {stats.get('top5_overlap')}", flush=True)

    print("\nConverting top-5 poses to MOL2 ...", flush=True)
    convert_top_poses_mol2("vina_screening_poses.pdbqt", "gnina_screening_poses.sdf", n=5)

    # Write comparative CSV
    fieldnames = ["name", "vina_rank", "vina_affinity_kcal_mol", "vina_status",
                  "gnina_rank", "gnina_cnn_score", "gnina_cnn_affinity_kcal_mol",
                  "gnina_vina_affinity_kcal_mol", "gnina_status", "rank_difference"]
    with open("comparative_screening_metrics.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for r in sorted(merged, key=lambda x: (x["vina_rank"] or 999)):
            writer.writerow({k: r.get(k, "") for k in fieldnames})
    print("\nWrote comparative_screening_metrics.csv", flush=True)

    # Load optional metadata
    vina_meta  = json.loads(Path("vina_docking_report.json").read_text()) \
                 if Path("vina_docking_report.json").exists() else {}
    gnina_meta = json.loads(Path("gnina_docking_report.json").read_text()) \
                 if Path("gnina_docking_report.json").exists() else {}

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print("Generating HTML report ...", flush=True)
    html = generate_html(merged, stats, vina_meta, gnina_meta, timestamp)
    Path("docking_performance_report.html").write_text(html)
    print(f"Wrote docking_performance_report.html ({len(html):,} chars)", flush=True)

    print("\nNode 05 completed", flush=True)


if __name__ == "__main__":
    main()
