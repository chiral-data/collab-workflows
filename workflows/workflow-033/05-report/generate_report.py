"""Node 05 — Report.

Reads the campaign manifest from node 04, re-applies filter thresholds,
ranks survivors by ipTM, generates scrambled negative controls (sequences
only — not co-folded; marked as pending for a follow-up run if desired),
computes the success rate, and produces three outputs:

  outputs/candidates.csv  — full ranked table of all shortlisted candidates
  outputs/summary.json    — success rate, filter params, top-20 designs
  outputs/report.html     — self-contained HTML report

Controls are composition-preserving scrambles of the top passed sequences.
They are recorded in the manifest with is_control=True but have no Boltz2
scores — a subsequent co-fold pass is needed to benchmark the distribution.
"""
from __future__ import annotations

import csv
import io
import json
import os
import pathlib
import random

# ── inputs ────────────────────────────────────────────────────────────────────

def _load_params() -> dict:
    """Load global params from PARAM_* env vars (silva) or inputs/global_params.json (local)."""
    try:
        base = json.loads(pathlib.Path("inputs/global_params.json").read_text())
    except FileNotFoundError:
        base = {}
    for key, val in os.environ.items():
        if not key.startswith("PARAM_"):
            continue
        param = key[6:].lower()
        if param in base:
            t = type(base[param])
            try:
                base[param] = t(val)
            except (ValueError, TypeError):
                base[param] = val
        else:
            try:
                base[param] = int(val)
            except ValueError:
                try:
                    base[param] = float(val)
                except ValueError:
                    base[param] = val
    return base

PARAMS   = _load_params()
SCORES   = json.loads(pathlib.Path("inputs/scores.json").read_text())
MANIFEST = json.loads(pathlib.Path("inputs/manifest.json").read_text())

FILTER_IPTM  = PARAMS["filter_iptm_min"]
FILTER_PLDDT = PARAMS["filter_binder_plddt_min"]
FILTER_RMSD  = PARAMS["filter_self_consistency_rmsd_max"]

N_CONTROLS = 5   # scrambled negatives to generate
N_TOP_HTML = 20  # candidates shown in the top-designs table

OUT = pathlib.Path("outputs")

# ── filter + rank ─────────────────────────────────────────────────────────────

def _passes(c: dict) -> bool:
    s = c.get("scores", {})
    checks = []
    if s.get("iptm") is not None:
        checks.append(s["iptm"] >= FILTER_IPTM)
    if s.get("binder_plddt") is not None:
        checks.append(s["binder_plddt"] >= FILTER_PLDDT)
    if s.get("self_consistency_rmsd") is not None:
        checks.append(s["self_consistency_rmsd"] <= FILTER_RMSD)
    return bool(checks) and all(checks)


def apply_filters(candidates: list[dict]) -> list[dict]:
    for c in candidates:
        c["passed_filter"] = _passes(c)
    return candidates


def rank_candidates(candidates: list[dict], by: str = "iptm",
                    descending: bool = True) -> list[dict]:
    real = [c for c in candidates if not c.get("is_control")]
    real = [c for c in real if c.get("scores", {}).get(by) is not None]
    return sorted(real, key=lambda c: c["scores"][by], reverse=descending)

# ── scrambled controls ────────────────────────────────────────────────────────

def _scramble(seq: str, seed: int) -> str:
    rng = random.Random(seed)
    chars = list(seq)
    rng.shuffle(chars)
    return "".join(chars)


def make_controls(passed: list[dict], n: int = N_CONTROLS, seed: int = 42) -> list[dict]:
    """Return N scrambled-sequence control records (not co-folded)."""
    if not passed:
        return []
    rng = random.Random(seed)
    controls = []
    for i in range(n):
        base = passed[i % len(passed)]
        ctrl_seq = _scramble(base["sequence"], seed=rng.randint(0, 2**31 - 1))
        controls.append({
            "id": f"ctrl_scrambled_{i:02d}",
            "backbone_id": base.get("backbone_id"),
            "sequence": ctrl_seq,
            "scores": {},
            "artifacts": {},
            "passed_filter": None,
            "is_control": True,
            "control_type": "scrambled",
            "created": base.get("created"),
            "note": "composition-preserving scramble; not co-folded",
        })
    return controls

# ── CSV ───────────────────────────────────────────────────────────────────────

SCORE_KEYS = ["iptm", "binder_plddt", "self_consistency_rmsd", "proteinmpnn_score"]

def write_csv(ranked: list[dict], controls: list[dict], path: pathlib.Path) -> None:
    cols = ["rank", "id", "backbone_id", "passed_filter",
            "is_control", "control_type"] + SCORE_KEYS + ["sequence"]
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        all_rows = ranked + controls
        for i, c in enumerate(all_rows, start=1):
            s = c.get("scores", {})
            w.writerow([
                i if not c.get("is_control") else "",
                c.get("id", ""),
                c.get("backbone_id", ""),
                c.get("passed_filter", ""),
                c.get("is_control", False),
                c.get("control_type", ""),
                *[_fmt(s.get(k)) for k in SCORE_KEYS],
                c.get("sequence", ""),
            ])


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.4f}"
    return str(v)

# ── summary JSON ──────────────────────────────────────────────────────────────

def write_summary(ranked: list[dict], all_real: list[dict], controls: list[dict],
                  path: pathlib.Path) -> dict:
    n_total   = len(all_real)
    n_scored  = sum(1 for c in all_real if c.get("scores", {}).get("iptm") is not None)
    n_passed  = sum(1 for c in all_real if c.get("passed_filter"))
    rate      = n_passed / max(n_scored, 1)

    top20 = []
    for c in ranked[:N_TOP_HTML]:
        s = c.get("scores", {})
        top20.append({
            "id": c["id"],
            "backbone_id": c.get("backbone_id"),
            "iptm": s.get("iptm"),
            "binder_plddt": s.get("binder_plddt"),
            "self_consistency_rmsd": s.get("self_consistency_rmsd"),
            "proteinmpnn_score": s.get("proteinmpnn_score"),
            "sequence": c.get("sequence", ""),
        })

    summary = {
        "target": MANIFEST.get("target", {}),
        "params": MANIFEST.get("params", {}),
        "filters": {
            "iptm_min": FILTER_IPTM,
            "binder_plddt_min": FILTER_PLDDT,
            "self_consistency_rmsd_max": FILTER_RMSD,
        },
        "results": {
            "n_shortlisted": n_total,
            "n_scored": n_scored,
            "n_passed_filter": n_passed,
            "success_rate": round(rate, 4),
            "n_controls_generated": len(controls),
        },
        "top_designs": top20,
    }
    path.write_text(json.dumps(summary, indent=2))
    return summary

# ── HTML report ───────────────────────────────────────────────────────────────

_CSS = """
:root {
  --bg: #0f1117; --surface: #1a1d27; --border: #2a2d3a;
  --text: #e2e5ef; --muted: #7b8099; --accent: #4f8ef7;
  --pass: #34c97a; --fail: #e05d5d; --warn: #f0a03c;
  --font: "Consolas","Menlo","Liberation Mono",monospace;
}
*, *::before, *::after { box-sizing: border-box; margin: 0; padding: 0; }
body {
  background: var(--bg); color: var(--text);
  font-family: var(--font); font-size: 13px; line-height: 1.6;
  padding: 24px 32px;
}
h1 { font-size: 20px; color: var(--accent); margin-bottom: 4px; }
h2 { font-size: 14px; color: var(--accent); margin: 28px 0 10px;
     border-bottom: 1px solid var(--border); padding-bottom: 4px; }
.meta { color: var(--muted); font-size: 12px; margin-bottom: 24px; }
.kv-grid {
  display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr));
  gap: 12px; margin-bottom: 20px;
}
.kv { background: var(--surface); border: 1px solid var(--border);
      border-radius: 6px; padding: 12px 16px; }
.kv .label { color: var(--muted); font-size: 11px; text-transform: uppercase;
             letter-spacing: .06em; margin-bottom: 2px; }
.kv .value { font-size: 22px; font-weight: 700; }
.kv .value.pass { color: var(--pass); }
.kv .value.fail { color: var(--fail); }
.kv .value.neutral { color: var(--accent); }
.table-wrap { overflow-x: auto; margin-bottom: 24px; }
table { border-collapse: collapse; width: 100%; min-width: 700px; }
th { background: var(--surface); color: var(--muted); font-size: 11px;
     text-transform: uppercase; letter-spacing: .05em;
     padding: 8px 12px; border-bottom: 1px solid var(--border);
     white-space: nowrap; text-align: left; }
td { padding: 7px 12px; border-bottom: 1px solid var(--border);
     white-space: nowrap; vertical-align: middle; }
tr:hover td { background: var(--surface); }
.pass-tag { color: var(--pass); font-weight: 700; }
.fail-tag { color: var(--fail); }
.bar-cell { min-width: 100px; }
.bar-wrap { background: var(--border); border-radius: 3px; height: 6px;
            overflow: hidden; width: 80px; display: inline-block;
            vertical-align: middle; margin-right: 6px; }
.bar-fill { height: 100%; border-radius: 3px; background: var(--accent); }
.bar-fill.pass { background: var(--pass); }
.seq { font-size: 11px; color: var(--muted); max-width: 220px;
       overflow: hidden; text-overflow: ellipsis; }
.note { background: var(--surface); border: 1px solid var(--border);
        border-left: 3px solid var(--warn); border-radius: 4px;
        padding: 10px 14px; color: var(--muted); font-size: 12px;
        margin-bottom: 16px; }
.ctrl-tag { color: var(--warn); font-size: 11px; }
"""


def _bar(value: float | None, lo: float, hi: float, pass_color: bool = True) -> str:
    if value is None:
        return "<span style='color:var(--muted)'>—</span>"
    pct = max(0.0, min(1.0, (value - lo) / max(hi - lo, 1e-9))) * 100
    cls = "pass" if pass_color else ""
    return (
        f'<span class="bar-wrap"><span class="bar-fill {cls}" '
        f'style="width:{pct:.0f}%"></span></span>'
        f'<span>{value:.3f}</span>'
    )


def _rmsd_bar(value: float | None) -> str:
    if value is None:
        return "<span style='color:var(--muted)'>—</span>"
    # For RMSD: low is good. Invert the bar so full bar = 0 Å, empty = 4+ Å
    cap = 4.0
    pct = max(0.0, min(1.0, 1.0 - value / cap)) * 100
    cls = "pass" if value <= FILTER_RMSD else ""
    return (
        f'<span class="bar-wrap"><span class="bar-fill {cls}" '
        f'style="width:{pct:.0f}%"></span></span>'
        f'<span>{value:.2f} Å</span>'
    )


def _score_color(value: float | None, threshold: float, higher_better: bool = True) -> str:
    if value is None:
        return f"<span style='color:var(--muted)'>—</span>"
    ok = (value >= threshold) if higher_better else (value <= threshold)
    col = "var(--pass)" if ok else "var(--fail)"
    return f"<span style='color:{col}'>{value:.4f}</span>"


def _top_table(ranked: list[dict]) -> str:
    rows = []
    for i, c in enumerate(ranked[:N_TOP_HTML], start=1):
        s = c.get("scores", {})
        iptm   = s.get("iptm")
        plddt  = s.get("binder_plddt")
        rmsd   = s.get("self_consistency_rmsd")
        nll    = s.get("proteinmpnn_score")
        passed = c.get("passed_filter")
        tag    = '<span class="pass-tag">PASS</span>' if passed else '<span class="fail-tag">FAIL</span>'
        seq    = c.get("sequence", "")
        seq_td = f'<span class="seq" title="{seq}">{seq[:30]}{"…" if len(seq) > 30 else ""}</span>'

        rows.append(
            f"<tr>"
            f"<td>{i}</td>"
            f"<td>{c['id']}</td>"
            f"<td>{c.get('backbone_id','')}</td>"
            f"<td>{tag}</td>"
            f'<td class="bar-cell">{_bar(iptm, 0.0, 1.0)}</td>'
            f'<td class="bar-cell">{_bar(plddt, 0.0, 100.0)}</td>'
            f'<td class="bar-cell">{_rmsd_bar(rmsd)}</td>'
            f"<td>{_score_color(nll, 1.5, higher_better=False)}</td>"
            f"<td>{seq_td}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


def _all_table(all_cands: list[dict]) -> str:
    rows = []
    for c in all_cands:
        s      = c.get("scores", {})
        iptm   = s.get("iptm")
        plddt  = s.get("binder_plddt")
        rmsd   = s.get("self_consistency_rmsd")
        nll    = s.get("proteinmpnn_score")
        passed = c.get("passed_filter")
        if passed is True:
            tag = '<span class="pass-tag">PASS</span>'
        elif passed is False:
            tag = '<span class="fail-tag">FAIL</span>'
        else:
            tag = '<span style="color:var(--muted)">—</span>'
        ctrl_tag = (
            f' <span class="ctrl-tag">[{c["control_type"]}]</span>'
            if c.get("is_control") else ""
        )
        rows.append(
            f"<tr>"
            f"<td>{c['id']}{ctrl_tag}</td>"
            f"<td>{tag}</td>"
            f"<td>{_fmt(iptm)}</td>"
            f"<td>{_fmt(plddt)}</td>"
            f"<td>{_fmt(rmsd)}</td>"
            f"<td>{_fmt(nll)}</td>"
            f"</tr>"
        )
    return "\n".join(rows)


def write_html(summary: dict, ranked: list[dict],
               all_real: list[dict], controls: list[dict],
               path: pathlib.Path) -> None:
    r = summary["results"]
    tgt = summary["target"]
    flt = summary["filters"]
    success_cls = "pass" if r["success_rate"] >= 0.1 else "fail"

    params_html = "".join(
        f'<div class="kv"><div class="label">{k}</div>'
        f'<div class="value neutral">{v}</div></div>'
        for k, v in summary["params"].items()
    )
    filter_html = "".join([
        f'<div class="kv"><div class="label">ipTM ≥</div>'
        f'<div class="value neutral">{flt["iptm_min"]}</div></div>',
        f'<div class="kv"><div class="label">pLDDT ≥</div>'
        f'<div class="value neutral">{flt["binder_plddt_min"]}</div></div>',
        f'<div class="kv"><div class="label">sc-RMSD ≤</div>'
        f'<div class="value neutral">{flt["self_consistency_rmsd_max"]} Å</div></div>',
    ])

    ctrl_note = ""
    if controls:
        ctrl_note = (
            '<div class="note">'
            f'<strong>{len(controls)} scrambled control sequences</strong> were generated '
            '(composition-preserving shuffles of top passed designs). '
            'They are listed in candidates.csv but have not been co-folded — '
            'run them through node 04 to benchmark the ipTM/sc-RMSD distribution.'
            '</div>'
        )

    all_sorted = sorted(
        all_real + controls,
        key=lambda c: (c.get("is_control", False),
                       -(c.get("scores", {}).get("iptm") or 0.0))
    )

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Binder Design Report — {tgt.get('pdb_id','?')} chain {tgt.get('chain','?')}</title>
<meta name="description" content="RFdiffusion → ProteinMPNN → Boltz2 binder design campaign results">
<style>{_CSS}</style>
</head>
<body>
<h1>Binder Design Report</h1>
<p class="meta">Target: {tgt.get('pdb_id','?')} chain {tgt.get('chain','?')} &nbsp;·&nbsp;
Hotspots: {tgt.get('hotspots','?')} &nbsp;·&nbsp;
Pipeline: RFdiffusion → ProteinMPNN → Boltz2</p>

<h2>Results</h2>
<div class="kv-grid">
  <div class="kv"><div class="label">Shortlisted</div>
    <div class="value neutral">{r['n_shortlisted']}</div></div>
  <div class="kv"><div class="label">Scored</div>
    <div class="value neutral">{r['n_scored']}</div></div>
  <div class="kv"><div class="label">Passed filter</div>
    <div class="value {success_cls}">{r['n_passed_filter']}</div></div>
  <div class="kv"><div class="label">Success rate</div>
    <div class="value {success_cls}">{r['success_rate']*100:.1f}%</div></div>
</div>

<h2>Filter thresholds</h2>
<div class="kv-grid">{filter_html}</div>

<h2>Campaign params</h2>
<div class="kv-grid">{params_html}</div>

<h2>Top designs (ranked by ipTM, passed only)</h2>
{ctrl_note}
<div class="table-wrap">
<table>
<thead><tr>
  <th>#</th><th>Candidate</th><th>Backbone</th><th>Filter</th>
  <th>ipTM ↑</th><th>pLDDT ↑</th><th>sc-RMSD ↓</th>
  <th>NLL ↓</th><th>Sequence (30 aa)</th>
</tr></thead>
<tbody>{_top_table(ranked)}</tbody>
</table>
</div>

<h2>All shortlisted candidates</h2>
<div class="table-wrap">
<table>
<thead><tr>
  <th>Candidate</th><th>Filter</th>
  <th>ipTM</th><th>pLDDT</th><th>sc-RMSD (Å)</th><th>NLL</th>
</tr></thead>
<tbody>{_all_table(all_sorted)}</tbody>
</table>
</div>

<p class="meta" style="margin-top:32px">
All scores are in-silico. Report a <strong>success rate</strong>, not just top
scores; always corroborate with self-consistency RMSD and NLL before
prioritising designs for experimental follow-up.
</p>
</body>
</html>"""

    path.write_text(html, encoding="utf-8")

# ── main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    candidates = MANIFEST.get("candidates", [])
    real_cands = [c for c in candidates if not c.get("is_control")]

    # Re-apply filters (ensures params changes in global_params.json take effect)
    apply_filters(real_cands)

    ranked  = rank_candidates(real_cands)
    passed  = [c for c in ranked if c.get("passed_filter")]
    controls = make_controls(passed)

    n_total  = len(real_cands)
    n_scored = sum(1 for c in real_cands if c.get("scores", {}).get("iptm") is not None)
    n_passed = len(passed)
    rate     = n_passed / max(n_scored, 1)

    print(
        f"[05] {n_total} candidates, {n_scored} scored, "
        f"{n_passed} passed ({rate*100:.1f}% success rate)",
        flush=True,
    )
    print(f"[05] generated {len(controls)} scrambled control sequences", flush=True)

    # CSV: ranked real candidates first, then controls
    write_csv(ranked, controls, OUT / "candidates.csv")
    print("[05] candidates.csv written", flush=True)

    summary = write_summary(ranked, real_cands, controls, OUT / "summary.json")
    print("[05] summary.json written", flush=True)

    write_html(summary, ranked, real_cands, controls, OUT / "report.html")
    print("[05] report.html written", flush=True)

    if ranked:
        top = ranked[0]["scores"]
        print(
            f"[05] top design: {ranked[0]['id']}  "
            f"ipTM={top.get('iptm','n/a')}  "
            f"pLDDT={top.get('binder_plddt','n/a')}  "
            f"sc-RMSD={top.get('self_consistency_rmsd','n/a')}",
            flush=True,
        )

    print("[05] report complete", flush=True)


if __name__ == "__main__":
    main()
