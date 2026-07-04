"""Node 02 — Generate Backbones.

Calls RFdiffusion NIM N times with distinct seeds to generate de novo
binder backbone PDBs conditioned on the target hotspots.

Inputs (mounted under inputs/ by the platform):
  target.pdb         — extracted target chain (from node 01)
  chain_seq.txt      — target chain one-letter sequence (for chain length)
  hotspots.json      — {"author": ["C56", ...], "seq_indices": [...]}
  global_params.json — workflow params

Outputs written to outputs/:
  backbones/bb_{i:04d}.pdb   — one backbone PDB per successful call
  backbone_list.json          — [{id, path, binder_chain, binder_length, seed}]
  gen_report.json             — summary (attempted, succeeded, failed)
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

PARAMS   = _load_params()
HOTSPOTS = json.loads(pathlib.Path("inputs/hotspots.json").read_text())
TARGET   = pathlib.Path("inputs/target.pdb").read_text()
TARGET_SEQ = pathlib.Path("inputs/chain_seq.txt").read_text().strip()

CHAIN          = PARAMS["target_chain"]
N_BACKBONES    = PARAMS["n_backbones"]
BINDER_LEN_MIN = PARAMS["binder_length_min"]
BINDER_LEN_MAX = PARAMS["binder_length_max"]
DIFFUSION_STEPS = PARAMS["diffusion_steps"]
NIM_MODE       = PARAMS.get("nim_mode", "hosted")

OUT      = pathlib.Path("outputs")
BB_DIR   = OUT / "backbones"

# ── API key ───────────────────────────────────────────────────────────────────

def _api_key() -> str | None:
    return os.environ.get("NGC_API_KEY") or os.environ.get("NVIDIA_API_KEY")

# ── contig + hotspot helpers ───────────────────────────────────────────────────

def chain_residue_range(pdb_text: str, chain: str) -> tuple[int, int]:
    """First and last author (PDB) residue numbers for a chain's CA atoms.

    RFdiffusion contigs reference actual PDB residue numbering, which rarely
    starts at 1 in real structures.
    """
    resnums = [
        int(line[22:26])
        for line in pdb_text.splitlines()
        if line.startswith("ATOM") and len(line) > 21
        and line[21] == chain and line[12:16].strip() == "CA"
    ]
    if not resnums:
        raise ValueError(f"no CA atoms found for chain {chain!r}")
    return resnums[0], resnums[-1]


def build_contig(target_chain: str, first_res: int, last_res: int,
                 binder_min: int, binder_max: int) -> str:
    """Build RFdiffusion contig string for binder design.

    Format: '{chain}{first}-{last}/0 {min}-{max}'
    Keeps the full target chain and appends a chain break followed by the
    generated binder segment.
    """
    return f"{target_chain}{first_res}-{last_res}/0 {binder_min}-{binder_max}"


def hotspot_res_list(hotspot_author: list[str]) -> list[str]:
    """Return hotspot_res list for RFdiffusion.

    Author specs like 'C56' are already in 'ChainResidue' format — pass through.
    Pure numeric specs (e.g. '56') get the target chain prepended.
    """
    out = []
    for spec in hotspot_author:
        if spec and spec[0].isalpha():
            out.append(spec)  # already has chain prefix
        else:
            out.append(f"{CHAIN}{spec}")
    return out

# ── PDB chain parsing ─────────────────────────────────────────────────────────

def _chain_ids(pdb_text: str) -> list[str]:
    """Distinct chain IDs from ATOM records, in order of first appearance."""
    seen: list[str] = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and len(line) > 21:
            cid = line[21]
            if cid not in seen:
                seen.append(cid)
    return seen


def _chain_ca_count(pdb_text: str, chain: str) -> int:
    """Number of CA atoms in the given chain (= residue count)."""
    n = 0
    for line in pdb_text.splitlines():
        if (line.startswith("ATOM") and len(line) > 21
                and line[21] == chain
                and line[12:16].strip() == "CA"):
            n += 1
    return n


def identify_binder_chain(backbone_pdb: str, target_chain: str) -> tuple[str, int]:
    """Return (binder_chain_id, binder_length) from an RFdiffusion output PDB.

    RFdiffusion places the target chain first and appends the generated binder
    as a new chain. We take the first chain that is not the target chain.
    """
    chains = _chain_ids(backbone_pdb)
    binder_chains = [c for c in chains if c != target_chain]
    if not binder_chains:
        raise ValueError(
            f"no binder chain found in backbone; chains present: {chains}"
        )
    binder_chain = binder_chains[0]
    binder_len = _chain_ca_count(backbone_pdb, binder_chain)
    return binder_chain, binder_len

# ── RFdiffusion NIM call ──────────────────────────────────────────────────────

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
            print(f"[02] HTTP {e.code}; retry {attempt}/{retries} in {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError("unreachable")


def _rfdiffusion_url() -> tuple[str, dict]:
    if NIM_MODE == "hosted":
        key = _api_key()
        if not key:
            raise EnvironmentError("NGC_API_KEY / NVIDIA_API_KEY not set; required in hosted mode")
        url = "https://health.api.nvidia.com/v1/biology/ipd/rfdiffusion/generate"
        headers = {"Content-Type": "application/json", "Authorization": f"Bearer {key}"}
    elif NIM_MODE == "local":
        url = "http://localhost:8000/biology/ipd/rfdiffusion/generate"
        headers = {"Content-Type": "application/json"}
    else:
        raise ValueError(f"unknown nim_mode {NIM_MODE!r}; expected 'hosted' or 'local'")
    return url, headers


def generate_backbone(contig: str, hotspot_res: list[str], seed: int) -> str:
    """Call RFdiffusion NIM once and return the output PDB string."""
    url, headers = _rfdiffusion_url()
    payload: dict = {
        "input_pdb": TARGET,
        "contigs": contig,
        "hotspot_res": hotspot_res,
        "diffusion_steps": DIFFUSION_STEPS,
        "random_seed": seed,
    }
    result = _nim_post(url, payload, headers)
    pdb = result.get("output_pdb")
    if not pdb:
        raise ValueError(f"RFdiffusion returned no output_pdb: {list(result.keys())}")
    elapsed = result.get("elapsed_ms")
    if elapsed is not None:
        print(f"[02]   server elapsed {elapsed} ms", flush=True)
    return pdb

# ── main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    BB_DIR.mkdir(parents=True, exist_ok=True)

    target_len  = len(TARGET_SEQ)
    first_res, last_res = chain_residue_range(TARGET, CHAIN)
    contig      = build_contig(CHAIN, first_res, last_res, BINDER_LEN_MIN, BINDER_LEN_MAX)
    hotspot_res = hotspot_res_list(HOTSPOTS["author"])

    print(f"[02] target chain={CHAIN} len={target_len}", flush=True)
    print(f"[02] contig: {contig}", flush=True)
    print(f"[02] hotspot_res: {hotspot_res}", flush=True)
    print(f"[02] generating {N_BACKBONES} backbones (diffusion_steps={DIFFUSION_STEPS})", flush=True)

    backbone_list: list[dict] = []
    succeeded = 0
    failed_ids: list[str] = []

    for i in range(N_BACKBONES):
        bb_id   = f"bb_{i:04d}"
        bb_path = BB_DIR / f"{bb_id}.pdb"
        seed    = i  # deterministic, distinct per backbone

        # Resume support: skip if already generated
        if bb_path.exists():
            pdb_text = bb_path.read_text()
            try:
                binder_chain, binder_len = identify_binder_chain(pdb_text, CHAIN)
                backbone_list.append({
                    "id": bb_id,
                    "path": str(bb_path),
                    "binder_chain": binder_chain,
                    "binder_length": binder_len,
                    "seed": seed,
                    "contig": contig,
                })
                succeeded += 1
                print(f"[02] {bb_id} skipped (cached) binder_chain={binder_chain} len={binder_len}", flush=True)
            except ValueError as e:
                print(f"[02] {bb_id} cached but unparseable: {e}", flush=True)
                failed_ids.append(bb_id)
            continue

        try:
            print(f"[02] {bb_id} seed={seed} …", flush=True)
            pdb_text = generate_backbone(contig, hotspot_res, seed)
            bb_path.write_text(pdb_text)

            binder_chain, binder_len = identify_binder_chain(pdb_text, CHAIN)
            backbone_list.append({
                "id": bb_id,
                "path": str(bb_path),
                "binder_chain": binder_chain,
                "binder_length": binder_len,
                "seed": seed,
                "contig": contig,
            })
            succeeded += 1
            print(f"[02] {bb_id} OK  binder_chain={binder_chain} len={binder_len}", flush=True)

        except Exception as e:
            print(f"[02] {bb_id} FAILED: {e}", flush=True)
            failed_ids.append(bb_id)

    if not backbone_list:
        raise RuntimeError(
            f"all {N_BACKBONES} RFdiffusion calls failed — cannot continue"
        )

    (OUT / "backbone_list.json").write_text(
        json.dumps(backbone_list, indent=2)
    )

    report = {
        "n_requested": N_BACKBONES,
        "n_succeeded": succeeded,
        "n_failed": len(failed_ids),
        "failed_ids": failed_ids,
        "contig": contig,
        "hotspot_res": hotspot_res,
        "diffusion_steps": DIFFUSION_STEPS,
        "nim_mode": NIM_MODE,
    }
    (OUT / "gen_report.json").write_text(json.dumps(report, indent=2))

    print(
        f"[02] done: {succeeded}/{N_BACKBONES} succeeded, "
        f"{len(failed_ids)} failed",
        flush=True,
    )
    if failed_ids:
        print(f"[02] failed: {failed_ids}", flush=True)


if __name__ == "__main__":
    main()
