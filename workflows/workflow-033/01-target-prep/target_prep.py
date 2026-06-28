"""Node 01 — Target Prep.

Fetches target PDB from RCSB, extracts the target chain, remaps hotspot
author residue numbers to 1-based sequence indices, calls MSA-search NIM
(hosted or local) for the target A3M.

Outputs written to outputs/:
  target.pdb        — extracted chain (ATOM records only)
  chain_seq.txt     — one-letter amino acid sequence
  hotspots.json     — {"author": [...], "seq_indices": [...]}
  target_a3m.txt    — MSA in A3M format (Uniref30 + env databases merged)
  prep_report.json  — summary metadata
"""
from __future__ import annotations

import io
import json
import os
import pathlib
import sys
import tarfile
import time
import urllib.error
import urllib.parse
import urllib.request

import pdb_utils

# ── inputs ────────────────────────────────────────────────────────────────────

PARAMS   = json.loads(pathlib.Path("inputs/global_params.json").read_text())
PDB_ID   = PARAMS["target_pdb_id"].upper()
CHAIN    = PARAMS["target_chain"]
HOTSPOTS = [h.strip() for h in PARAMS["hotspot_residues"].split(",")]
NIM_MODE = PARAMS.get("nim_mode", "hosted")
OUT      = pathlib.Path("outputs")

# ── API key ───────────────────────────────────────────────────────────────────

def _api_key() -> str | None:
    return os.environ.get("NGC_API_KEY") or os.environ.get("NVIDIA_API_KEY")

# ── hotspot helpers ────────────────────────────────────────────────────────────

def _strip_chain_prefix(spec: str) -> str:
    """Strip leading chain letter(s) from a hotspot spec.

    'C56' → '56', 'C56A' → '56A', '56' → '56'.
    """
    i = 0
    while i < len(spec) and spec[i].isalpha():
        i += 1
    return spec[i:]


def remap_hotspots(chain_pdb: str, specs: list[str]) -> list[int]:
    mapping = pdb_utils.residue_index_map(chain_pdb, CHAIN)
    out, missing = [], []
    for spec in specs:
        key = _strip_chain_prefix(spec)
        if key in mapping:
            out.append(mapping[key])
        else:
            missing.append(spec)
    if missing:
        raise ValueError(f"hotspot residues not found in chain {CHAIN!r}: {missing}")
    return out

# ── RCSB fetch ────────────────────────────────────────────────────────────────

def fetch_pdb(pdb_id: str) -> str:
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"[01] fetching {url}", flush=True)
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read().decode("utf-8", errors="replace")

# ── MSA-search NIM ────────────────────────────────────────────────────────────

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
            print(f"[01] HTTP {e.code}; retry {attempt}/{retries} in {wait}s", flush=True)
            time.sleep(wait)
    raise RuntimeError("unreachable")


def _merge_a3m(alignments: dict) -> str:
    """Concatenate a3m alignments from multiple databases into one block.

    Keeps one query header (from the first database) and appends sequences
    from subsequent databases, skipping their redundant query headers.
    """
    parts: list[str] = []
    for _db, formats in alignments.items():
        a3m = formats.get("a3m", {}).get("alignment", "")
        if a3m:
            parts.append(a3m.rstrip("\n"))
    if not parts:
        return ""
    merged = parts[0]
    for extra in parts[1:]:
        lines = extra.splitlines()
        # Skip the query header + sequence line from subsequent databases
        skip = 2 if (lines and lines[0].startswith(">")) else 0
        tail = "\n".join(lines[skip:]).strip()
        if tail:
            merged += "\n" + tail
    return merged + "\n"


def _sanitize_a3m(a3m: str) -> str:
    """Replace characters Boltz2 rejects with X/x; preserve gaps and case."""
    UPPER_OK = set("ACDEFGHIKLMNPQRSTVWYX")
    LOWER_OK = set("acdefghiklmnpqrstvwyx")
    out = []
    for ln in a3m.splitlines():
        if ln.startswith(">") or not ln:
            out.append(ln)
            continue
        fixed = []
        for c in ln:
            if c in UPPER_OK or c in LOWER_OK or c in "-.":
                fixed.append(c)
            elif c.islower():
                fixed.append("x")
            else:
                fixed.append("X")
        out.append("".join(fixed))
    return "\n".join(out) + "\n"


_MSA_PAYLOAD = {
    "databases": ["Uniref30_2302", "colabfold_envdb_202108"],
    "e_value": 0.0001,
    "output_alignment_formats": ["a3m"],
}


def call_msa_nim_hosted(sequence: str) -> str:
    key = _api_key()
    if not key:
        raise EnvironmentError("NGC_API_KEY / NVIDIA_API_KEY not set; required in hosted mode")
    url = "https://health.api.nvidia.com/v1/biology/colabfold/msa-search/predict"
    headers = {"Content-Type": "application/json", "Authorization": f"Bearer {key}"}
    payload = {"sequence": sequence, **_MSA_PAYLOAD}
    print("[01] calling MSA-search NIM (hosted) …", flush=True)
    result = _nim_post(url, payload, headers)
    return _sanitize_a3m(_merge_a3m(result.get("alignments", {})))


def call_msa_nim_local(sequence: str) -> str:
    url = "http://localhost:8000/biology/colabfold/msa-search/predict"
    headers = {"Content-Type": "application/json"}
    payload = {"sequence": sequence, **_MSA_PAYLOAD}
    print("[01] calling MSA-search NIM (local) …", flush=True)
    result = _nim_post(url, payload, headers)
    return _sanitize_a3m(_merge_a3m(result.get("alignments", {})))


def call_msa_colabfold_public(sequence: str) -> str:
    """Fallback: public api.colabfold.com (rate-limited, for dev/testing only)."""
    HOST = "https://api.colabfold.com"
    UA   = "chiral-binder-design/1.0"

    def _post(path, data):
        req = urllib.request.Request(
            f"{HOST}/{path}",
            data=urllib.parse.urlencode(data).encode(),
            headers={"User-Agent": UA}, method="POST",
        )
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.loads(r.read().decode())

    def _get_json(path):
        req = urllib.request.Request(f"{HOST}/{path}", headers={"User-Agent": UA})
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.loads(r.read().decode())

    def _download(ticket):
        req = urllib.request.Request(
            f"{HOST}/result/download/{ticket}", headers={"User-Agent": UA}
        )
        with urllib.request.urlopen(req, timeout=300) as r:
            return r.read()

    print("[01] calling ColabFold public API (fallback) …", flush=True)
    sub = _post("ticket/msa", {"q": f">target\n{sequence}\n", "mode": "env"})
    tid = sub.get("id")
    status = sub.get("status")
    if not tid:
        raise RuntimeError(f"ColabFold submission failed: {sub}")
    print(f"[01] ColabFold ticket={tid} status={status}", flush=True)
    waited = 0
    while status in ("PENDING", "RUNNING", "UNKNOWN", None):
        time.sleep(10)
        waited += 10
        if waited > 900:
            raise RuntimeError("ColabFold timed out after 900s")
        status = _get_json(f"ticket/{tid}").get("status")
        print(f"[01]   … {waited}s status={status}", flush=True)
    if status != "COMPLETE":
        raise RuntimeError(f"ColabFold ended with status={status}")
    tar_bytes = _download(tid)
    parts = []
    with tarfile.open(fileobj=io.BytesIO(tar_bytes)) as tf:
        for m in sorted(tf.getmembers(), key=lambda m: m.name):
            if m.name.endswith(".a3m"):
                parts.append(tf.extractfile(m).read().decode(errors="replace"))
    if not parts:
        raise RuntimeError("no .a3m in ColabFold result tar")
    merged = "\n".join(p.rstrip("\n") for p in parts) + "\n"
    return _sanitize_a3m(merged)


def call_msa(sequence: str) -> str:
    if NIM_MODE == "hosted":
        return call_msa_nim_hosted(sequence)
    elif NIM_MODE == "local":
        return call_msa_nim_local(sequence)
    else:
        raise ValueError(f"unknown nim_mode {NIM_MODE!r}; expected 'hosted' or 'local'")

# ── main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)

    # 1. Fetch full PDB
    raw_pdb = fetch_pdb(PDB_ID)

    # 2. Extract target chain
    chain_pdb = pdb_utils.extract_chain(raw_pdb, CHAIN)
    if not chain_pdb:
        raise ValueError(f"chain {CHAIN!r} not found in {PDB_ID}")
    (OUT / "target.pdb").write_text(chain_pdb)
    print(f"[01] target.pdb: {len(chain_pdb.splitlines())} lines", flush=True)

    # 3. One-letter sequence
    seq = pdb_utils.sequence(chain_pdb, CHAIN)
    if not seq:
        raise ValueError(f"no CA atoms for chain {CHAIN!r} in {PDB_ID}")
    (OUT / "chain_seq.txt").write_text(seq)
    print(f"[01] chain {CHAIN}: {len(seq)} aa", flush=True)

    # 4. Remap hotspots → 1-based sequence indices
    seq_indices = remap_hotspots(chain_pdb, HOTSPOTS)
    hotspot_data = {"author": HOTSPOTS, "seq_indices": seq_indices}
    (OUT / "hotspots.json").write_text(json.dumps(hotspot_data, indent=2))
    print(f"[01] hotspots: {list(zip(HOTSPOTS, seq_indices))}", flush=True)

    # 5. MSA for target chain
    a3m = call_msa(seq)
    n_seqs = a3m.count(">")
    (OUT / "target_a3m.txt").write_text(a3m)
    print(f"[01] target_a3m.txt: ~{n_seqs} sequences, {len(a3m)} chars", flush=True)

    # 6. Prep report
    report = {
        "pdb_id": PDB_ID,
        "chain": CHAIN,
        "seq_length": len(seq),
        "hotspots_author": HOTSPOTS,
        "hotspots_seq_indices": seq_indices,
        "msa_sequences": n_seqs,
        "nim_mode": NIM_MODE,
    }
    (OUT / "prep_report.json").write_text(json.dumps(report, indent=2))
    print("[01] target-prep complete", flush=True)


if __name__ == "__main__":
    main()
