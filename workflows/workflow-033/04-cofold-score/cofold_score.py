"""Node 04 — Co-fold and Score.

Shortlists top candidates by ProteinMPNN score (lower NLL = better), then
co-folds each binder–target complex with Boltz2 NIM. Extracts ipTM, binder
pLDDT, and self-consistency RMSD from the mmCIF output. Writes the run
manifest used by node 05 for filtering and ranking.

Chain assignment in every Boltz2 call:
  id "A" — binder (designed sequence; single-sequence, no fabricated MSA)
  id "B" — target (target sequence + target A3M from node 01)

Inputs (mounted under inputs/ by the platform):
  sequence_manifest.json  — flat list of candidates from node 03
  chain_seq.txt           — target chain one-letter sequence
  target_a3m.txt          — target MSA (A3M format) from node 01
  global_params.json      — workflow params

Outputs written to outputs/:
  complexes/{candidate_id}.cif  — Boltz2 co-folded complex (mmCIF)
  scores.json                   — {candidate_id: {iptm, binder_plddt, sc_rmsd, ...}}
  manifest.json                 — full campaign manifest (consumed by node 05)
  cofold_report.json            — run summary
"""
from __future__ import annotations

import json
import os
import pathlib
import time
import urllib.error
import urllib.request
from datetime import datetime, timezone

import metrics  # local: metrics.py alongside this script

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

PARAMS       = _load_params()
SEQ_MANIFEST = json.loads(pathlib.Path("inputs/sequence_manifest.json").read_text())
TARGET_SEQ   = pathlib.Path("inputs/chain_seq.txt").read_text().strip()
TARGET_A3M   = pathlib.Path("inputs/target_a3m.txt").read_text()

SHORTLIST_N      = PARAMS["cofold_shortlist_n"]
RECYCLING_STEPS  = PARAMS.get("boltz2_recycling_steps", 3)
SAMPLING_STEPS   = PARAMS.get("boltz2_sampling_steps", 50)
NIM_MODE         = PARAMS.get("nim_mode", "hosted")
THROTTLE_S       = 5.0          # seconds between hosted calls to avoid 429s

FILTER_IPTM      = PARAMS["filter_iptm_min"]
FILTER_PLDDT     = PARAMS["filter_binder_plddt_min"]
FILTER_RMSD      = PARAMS["filter_self_consistency_rmsd_max"]

OUT         = pathlib.Path("outputs")
COMPLEX_DIR = OUT / "complexes"

# ── API key ───────────────────────────────────────────────────────────────────

def _api_key() -> str | None:
    return os.environ.get("NGC_API_KEY") or os.environ.get("NVIDIA_API_KEY")

# ── shortlist ─────────────────────────────────────────────────────────────────

def build_shortlist(candidates: list[dict], n: int) -> list[dict]:
    """Return top-N candidates by ProteinMPNN score (lower = better NLL).

    Candidates without a score (score is None or missing) are placed at the
    end and still included if slots remain.
    """
    scored   = [c for c in candidates if c.get("proteinmpnn_score") is not None]
    unscored = [c for c in candidates if c.get("proteinmpnn_score") is None]
    scored.sort(key=lambda c: c["proteinmpnn_score"])
    combined = scored + unscored
    return combined[:n]

# ── Boltz2 NIM call ───────────────────────────────────────────────────────────

def _nim_post(url: str, payload: dict, headers: dict,
              retries: int = 5, base_delay: float = 10.0) -> dict:
    body = json.dumps(payload).encode()
    for attempt in range(1, retries + 1):
        req = urllib.request.Request(url, data=body, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=1200) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if attempt == retries or e.code not in (429, 500, 502, 503, 504):
                raise
            ra = e.headers.get("Retry-After") if e.headers else None
            wait = float(ra) if (ra and str(ra).isdigit()) else base_delay * (2 ** attempt)
            wait = min(wait, 120.0)
            print(f"[04] HTTP {e.code}; retry {attempt}/{retries} in {wait:.0f}s", flush=True)
            time.sleep(wait)
        except (urllib.error.URLError, TimeoutError) as e:
            if attempt == retries:
                raise
            wait = min(base_delay * (2 ** attempt), 120.0)
            print(f"[04] {type(e).__name__}; retry {attempt}/{retries} in {wait:.0f}s", flush=True)
            time.sleep(wait)
    raise RuntimeError("_nim_post exhausted retries")


def _boltz2_url_headers() -> tuple[str, dict]:
    if NIM_MODE == "hosted":
        key = _api_key()
        if not key:
            raise EnvironmentError("NGC_API_KEY / NVIDIA_API_KEY not set; required in hosted mode")
        url  = "https://health.api.nvidia.com/v1/biology/mit/boltz2/predict"
        hdrs = {"Content-Type": "application/json", "Authorization": f"Bearer {key}"}
    elif NIM_MODE == "local":
        url  = "http://localhost:8000/biology/mit/boltz2/predict"
        hdrs = {"Content-Type": "application/json"}
    else:
        raise ValueError(f"unknown nim_mode {NIM_MODE!r}")
    return url, hdrs


def cofold(binder_seq: str) -> dict:
    """Co-fold binder (chain A) + target (chain B + MSA) with Boltz2."""
    url, hdrs = _boltz2_url_headers()
    payload = {
        "polymers": [
            {
                "id": "A",
                "molecule_type": "protein",
                "sequence": binder_seq,
                # de novo binder: single-sequence only, no fabricated MSA
            },
            {
                "id": "B",
                "molecule_type": "protein",
                "sequence": TARGET_SEQ,
                "msa": {
                    "msa_search": {
                        "a3m": {
                            "alignment": TARGET_A3M,
                            "format": "a3m",
                            "rank": 0,
                        }
                    }
                },
            },
        ],
        "recycling_steps":  RECYCLING_STEPS,
        "sampling_steps":   SAMPLING_STEPS,
        "diffusion_samples": 1,
        "step_scale": 1.638,
        "output_format": "mmcif",
    }
    return _nim_post(url, payload, hdrs)

# ── manifest helpers ──────────────────────────────────────────────────────────

def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _make_manifest(shortlist: list[dict]) -> dict:
    return {
        "schema_version": "1.0",
        "campaign": "protein-binder-design",
        "created": _now(),
        "target": {
            "pdb_id": PARAMS["target_pdb_id"],
            "chain": PARAMS["target_chain"],
            "hotspots": PARAMS["hotspot_residues"],
        },
        "mode": NIM_MODE,
        "params": {
            "n_backbones": PARAMS["n_backbones"],
            "seqs_per_backbone": PARAMS["seqs_per_backbone"],
            "binder_len": f"{PARAMS['binder_length_min']}-{PARAMS['binder_length_max']}",
            "cofold_shortlist_n": SHORTLIST_N,
            "boltz2_recycling_steps": RECYCLING_STEPS,
            "boltz2_sampling_steps": SAMPLING_STEPS,
        },
        "filters": {
            "iptm_min": FILTER_IPTM,
            "binder_plddt_min": FILTER_PLDDT,
            "self_consistency_rmsd_max": FILTER_RMSD,
        },
        "candidates": [],
    }


def _apply_filters(manifest: dict) -> None:
    f = manifest["filters"]
    for c in manifest["candidates"]:
        s = c.get("scores", {})
        checks = []
        if f.get("iptm_min") is not None and s.get("iptm") is not None:
            checks.append(s["iptm"] >= f["iptm_min"])
        if f.get("binder_plddt_min") is not None and s.get("binder_plddt") is not None:
            checks.append(s["binder_plddt"] >= f["binder_plddt_min"])
        if f.get("self_consistency_rmsd_max") is not None and s.get("self_consistency_rmsd") is not None:
            checks.append(s["self_consistency_rmsd"] <= f["self_consistency_rmsd_max"])
        c["passed_filter"] = bool(checks) and all(checks)

# ── main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    COMPLEX_DIR.mkdir(parents=True, exist_ok=True)

    shortlist = build_shortlist(SEQ_MANIFEST, SHORTLIST_N)
    print(
        f"[04] shortlisted {len(shortlist)}/{len(SEQ_MANIFEST)} candidates "
        f"for Boltz2 co-folding",
        flush=True,
    )

    manifest    = _make_manifest(shortlist)
    scores_out  = {}
    succeeded   = 0
    failed_ids: list[str] = []

    for idx, cand in enumerate(shortlist):
        cid          = cand["candidate_id"]
        bb_id        = cand["backbone_id"]
        binder_chain = cand["binder_chain"]
        binder_seq   = cand["sequence"]
        # cand["backbone_path"] is root-relative (e.g. "backbones/bb_0000.pdb"),
        # forwarded unchanged from node 02 via node 03.
        backbone_path = pathlib.Path("inputs") / cand["backbone_path"]
        cif_path     = COMPLEX_DIR / f"{cid}.cif"

        # Seed manifest entry for resumability
        mc = {
            "id": cid,
            "backbone_id": bb_id,
            "sequence": binder_seq,
            "scores": {
                "proteinmpnn_score": cand.get("proteinmpnn_score"),
            },
            "artifacts": {
                "backbone_pdb": str(backbone_path),
            },
            "passed_filter": None,
            "is_control": False,
            "control_type": None,
            "created": _now(),
        }

        # Resume: load existing CIF if already computed
        nim_result: dict | None = None
        if cif_path.exists():
            print(f"[04] {cid} skipped (cached)", flush=True)
            cif_content = cif_path.read_text()
        else:
            if not backbone_path.exists():
                print(f"[04] {cid} SKIP: backbone not found {backbone_path}", flush=True)
                failed_ids.append(cid)
                manifest["candidates"].append(mc)
                continue

            try:
                print(
                    f"[04] {cid} ({idx+1}/{len(shortlist)}) "
                    f"binder={len(binder_seq)}aa target={len(TARGET_SEQ)}aa …",
                    flush=True,
                )
                nim_result  = cofold(binder_seq)
                structures  = nim_result.get("structures", [])
                if not structures:
                    raise ValueError("Boltz2 returned no structures")
                cif_content = structures[0]["structure"]
                cif_path.write_text(cif_content)

                elapsed = nim_result.get("metrics", {}).get("elapsed_ms")
                if elapsed:
                    print(f"[04]   elapsed {elapsed} ms", flush=True)

                # Throttle between hosted calls to avoid 429s
                if NIM_MODE == "hosted" and idx < len(shortlist) - 1:
                    time.sleep(THROTTLE_S)

            except Exception as e:
                print(f"[04] {cid} FAILED: {e}", flush=True)
                failed_ids.append(cid)
                manifest["candidates"].append(mc)
                continue

        # Extract scores
        backbone_pdb = backbone_path.read_text() if backbone_path.exists() else ""

        # ipTM is only in the NIM response, not persisted in the CIF. On resume
        # paths (nim_result is None) ipTM will be None; node 05 handles missing scores.
        conf_scores = nim_result.get("confidence_scores", []) if nim_result else []
        iptm        = conf_scores[0] if conf_scores else None

        try:
            binder_plddt = metrics.mean_plddt_from_cif(cif_content, chain_id="A")
        except Exception as e:
            print(f"[04]   pLDDT extraction failed: {e}", flush=True)
            binder_plddt = None

        try:
            sc_rmsd = metrics.self_consistency_rmsd(
                cif_content, backbone_pdb,
                binder_cif_chain="A", binder_pdb_chain=binder_chain,
            )
        except Exception as e:
            print(f"[04]   sc-RMSD failed: {e}", flush=True)
            sc_rmsd = None

        sc = {"iptm": iptm, "binder_plddt": binder_plddt, "self_consistency_rmsd": sc_rmsd}

        print(
            f"[04] {cid} OK  iptm={sc['iptm']:.3f if sc['iptm'] else 'n/a'}  "
            f"plddt={sc['binder_plddt']:.1f if sc['binder_plddt'] else 'n/a'}  "
            f"sc_rmsd={sc['self_consistency_rmsd']:.2f if sc['self_consistency_rmsd'] else 'n/a'}",
            flush=True,
        )

        mc["scores"].update(sc)
        mc["artifacts"]["complex_cif"] = str(cif_path)
        scores_out[cid] = mc["scores"]
        manifest["candidates"].append(mc)
        succeeded += 1

    # Apply filters to all candidates
    _apply_filters(manifest)

    n_passed = sum(1 for c in manifest["candidates"] if c.get("passed_filter"))

    # Write outputs
    (OUT / "manifest.json").write_text(json.dumps(manifest, indent=2))
    (OUT / "scores.json").write_text(json.dumps(scores_out, indent=2))

    report = {
        "n_shortlisted": len(shortlist),
        "n_cofolded": succeeded,
        "n_failed": len(failed_ids),
        "n_passed_filter": n_passed,
        "failed_ids": failed_ids,
        "filters": manifest["filters"],
        "nim_mode": NIM_MODE,
    }
    (OUT / "cofold_report.json").write_text(json.dumps(report, indent=2))

    print(
        f"[04] done: {succeeded}/{len(shortlist)} co-folded, "
        f"{n_passed} passed filters, {len(failed_ids)} failed",
        flush=True,
    )
    if failed_ids:
        print(f"[04] failed: {failed_ids}", flush=True)


if __name__ == "__main__":
    main()
