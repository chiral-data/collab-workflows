"""Node 03 — Design Sequences.

Calls ProteinMPNN NIM once per backbone to design binder sequences.
Redesigns the binder chain only; the target chain stays fixed.
Drops the native/WT row from the mfasta output before pairing with scores.

Inputs (mounted under inputs/ by the platform):
  backbone_list.json  — [{id, path, binder_chain, binder_length, seed, contig}]
  global_params.json  — workflow params

Outputs written to outputs/:
  sequences/{bb_id}/seqs.fa       — full mfasta per backbone (retained for audit)
  sequence_manifest.json          — flat list of all candidates for node 04
  seq_report.json                 — summary (backbones processed, total candidates)
"""
from __future__ import annotations

import json
import os
import pathlib
import time
import urllib.error
import urllib.request

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

PARAMS  = _load_params()
BB_LIST = json.loads(pathlib.Path("inputs/backbone_list.json").read_text())

SEQS_PER_BACKBONE = PARAMS["seqs_per_backbone"]
SAMPLING_TEMP     = PARAMS.get("proteinmpnn_sampling_temp", 0.1)
NIM_MODE          = PARAMS.get("nim_mode", "hosted")

OUT     = pathlib.Path("outputs")
SEQ_DIR = OUT / "sequences"

# ── API key ───────────────────────────────────────────────────────────────────

def _api_key() -> str | None:
    return os.environ.get("NGC_API_KEY") or os.environ.get("NVIDIA_API_KEY")

# ── mfasta helpers ────────────────────────────────────────────────────────────

def _parse_mfasta(mfasta: str) -> list[tuple[str, str]]:
    """Return list of (header, sequence) from a multi-FASTA string."""
    entries: list[tuple[str, str]] = []
    header = seq_lines = None
    for line in mfasta.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                entries.append((header, "".join(seq_lines)))
            header = line[1:]  # strip leading >
            seq_lines = []
        elif header is not None:
            seq_lines.append(line)
    if header is not None and seq_lines is not None:
        entries.append((header, "".join(seq_lines)))
    return entries


def _is_native(header: str) -> bool:
    """Return True if this mfasta entry is the native/WT scaffold sequence."""
    h = header.lower()
    return "native" in h or ", wt," in h or h.startswith("wt,") or h.endswith(", wt")


def parse_designed(mfasta: str, scores: list[float]) -> list[dict]:
    """Extract designed (non-native) sequences and pair with NIM scores.

    The mfasta may prepend one native/WT row. Scores from the response
    correspond only to designed rows — we skip native rows before pairing.
    """
    all_entries = _parse_mfasta(mfasta)
    designed = [(h, s) for h, s in all_entries if not _is_native(h)]

    if len(designed) != len(scores):
        # Tolerate a count mismatch with a warning rather than hard-failing;
        # take the shorter side so we never index out of range.
        print(
            f"[03] warning: {len(designed)} designed entries vs "
            f"{len(scores)} scores — taking min({len(designed)}, {len(scores)})",
            flush=True,
        )
        n = min(len(designed), len(scores))
        designed = designed[:n]
        scores   = scores[:n]

    out = []
    for (header, sequence), score in zip(designed, scores):
        out.append({"header": header, "sequence": sequence, "proteinmpnn_score": score})
    return out


def _extract_temp(header: str) -> float | None:
    """Parse sampling temperature from mfasta header like 'T=0.1, score=...'."""
    for part in header.split(","):
        part = part.strip()
        if part.startswith("T="):
            try:
                return float(part[2:])
            except ValueError:
                pass
    return None

# ── ProteinMPNN NIM call ──────────────────────────────────────────────────────

def _nim_post(url: str, payload: dict, headers: dict, retries: int = 3) -> dict:
    body = json.dumps(payload).encode()
    for attempt in range(1, retries + 1):
        req = urllib.request.Request(url, data=body, headers=headers, method="POST")
        try:
            with urllib.request.urlopen(req, timeout=300) as r:
                return json.loads(r.read().decode())
        except urllib.error.HTTPError as e:
            if attempt == retries or e.code not in (429, 500, 502, 503, 504):
                raise
            wait = 2 ** attempt
            print(f"[03] HTTP {e.code}; retry {attempt}/{retries} in {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError("unreachable")


def _proteinmpnn_url() -> tuple[str, dict]:
    if NIM_MODE == "hosted":
        key = _api_key()
        if not key:
            raise EnvironmentError("NGC_API_KEY / NVIDIA_API_KEY not set; required in hosted mode")
        url = "https://health.api.nvidia.com/v1/biology/ipd/proteinmpnn/predict"
        hdrs = {"Content-Type": "application/json", "Authorization": f"Bearer {key}"}
    elif NIM_MODE == "local":
        url = "http://localhost:8000/biology/ipd/proteinmpnn/predict"
        hdrs = {"Content-Type": "application/json"}
    else:
        raise ValueError(f"unknown nim_mode {NIM_MODE!r}; expected 'hosted' or 'local'")
    return url, hdrs


def design_sequences(backbone_pdb: str, binder_chain: str) -> tuple[str, list[float]]:
    """Call ProteinMPNN NIM and return (mfasta, scores)."""
    url, hdrs = _proteinmpnn_url()
    payload = {
        "input_pdb": backbone_pdb,
        "input_pdb_chains": [binder_chain],
        "num_seq_per_target": SEQS_PER_BACKBONE,
        "sampling_temp": [SAMPLING_TEMP],
        "use_soluble_model": True,
    }
    result = _nim_post(url, payload, hdrs)
    mfasta = result.get("mfasta", "")
    scores = result.get("scores", [])
    if not mfasta:
        raise ValueError("ProteinMPNN returned empty mfasta")
    return mfasta, scores

# ── main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    SEQ_DIR.mkdir(parents=True, exist_ok=True)

    print(
        f"[03] designing sequences for {len(BB_LIST)} backbones "
        f"({SEQS_PER_BACKBONE} seqs each, T={SAMPLING_TEMP})",
        flush=True,
    )

    sequence_manifest: list[dict] = []
    succeeded = 0
    failed_ids: list[str] = []

    for bb in BB_LIST:
        bb_id        = bb["id"]
        # bb["path"] is root-relative (e.g. "backbones/bb_0000.pdb"); silva
        # copies node 02's outputs/ tree into this node's inputs/.
        bb_path      = pathlib.Path("inputs") / bb["path"]
        binder_chain = bb["binder_chain"]
        binder_len   = bb["binder_length"]

        bb_seq_dir = SEQ_DIR / bb_id
        seqs_fa    = bb_seq_dir / "seqs.fa"

        # Resume support: load cached mfasta if present
        if seqs_fa.exists():
            print(f"[03] {bb_id} skipped (cached)", flush=True)
            mfasta = seqs_fa.read_text()
            # Scores can't be recovered from mfasta alone; mark as None so
            # node 04 falls back to sequence-order shortlisting.
            scores = []
        else:
            if not bb_path.exists():
                print(f"[03] {bb_id} SKIP backbone not found: {bb_path}", flush=True)
                failed_ids.append(bb_id)
                continue
            try:
                backbone_pdb = bb_path.read_text()
                print(
                    f"[03] {bb_id} binder_chain={binder_chain} len={binder_len} …",
                    flush=True,
                )
                mfasta, scores = design_sequences(backbone_pdb, binder_chain)
                bb_seq_dir.mkdir(parents=True, exist_ok=True)
                seqs_fa.write_text(mfasta)
                print(
                    f"[03] {bb_id} OK  {len(scores)} scores, "
                    f"{mfasta.count('>')} FASTA entries",
                    flush=True,
                )
            except Exception as e:
                print(f"[03] {bb_id} FAILED: {e}", flush=True)
                failed_ids.append(bb_id)
                continue

        # Parse designed entries and build manifest rows
        designed = parse_designed(mfasta, scores)
        if not designed:
            print(f"[03] {bb_id} warning: no designed sequences parsed", flush=True)
            failed_ids.append(bb_id)
            continue

        for j, entry in enumerate(designed):
            candidate_id = f"{bb_id}_seq{j:03d}"
            sequence_manifest.append({
                "candidate_id": candidate_id,
                "backbone_id": bb_id,
                # Root-relative (e.g. "backbones/bb_0000.pdb"), matching bb["path"]
                # above — node 04 resolves it under its own inputs/.
                "backbone_path": bb["path"],
                "binder_chain": binder_chain,
                "binder_length": binder_len,
                "sequence": entry["sequence"],
                "proteinmpnn_score": entry["proteinmpnn_score"],
                "mfasta_header": entry["header"],
                "sampling_temp": _extract_temp(entry["header"]) or SAMPLING_TEMP,
            })

        succeeded += 1

    if not sequence_manifest:
        raise RuntimeError("no sequences produced — cannot continue to co-folding")

    (OUT / "sequence_manifest.json").write_text(
        json.dumps(sequence_manifest, indent=2)
    )

    report = {
        "n_backbones_attempted": len(BB_LIST),
        "n_backbones_succeeded": succeeded,
        "n_backbones_failed": len(failed_ids),
        "failed_backbone_ids": failed_ids,
        "n_total_candidates": len(sequence_manifest),
        "seqs_per_backbone": SEQS_PER_BACKBONE,
        "sampling_temp": SAMPLING_TEMP,
        "nim_mode": NIM_MODE,
    }
    (OUT / "seq_report.json").write_text(json.dumps(report, indent=2))

    print(
        f"[03] done: {succeeded}/{len(BB_LIST)} backbones, "
        f"{len(sequence_manifest)} total candidates",
        flush=True,
    )
    if failed_ids:
        print(f"[03] failed backbones: {failed_ids}", flush=True)


if __name__ == "__main__":
    main()
