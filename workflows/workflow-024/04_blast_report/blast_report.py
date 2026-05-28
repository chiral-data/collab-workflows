#!/usr/bin/env python3
"""
Node 04: BLAST Validation and Report Generation.

1. Builds a local BLAST+ database from target.fasta + exclusion.fasta.
2. BLASTs the top N primer sets (forward, reverse, probe) for target specificity.
3. Writes results.json, results.tsv, and report.txt to ./outputs/.

Set PARAM_SKIP_BLAST=true to skip BLAST and report Primer3-only results.
"""

import json
import logging
import os
import shutil
import subprocess
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from datetime import datetime

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")
log = logging.getLogger("blast_report")


def _load_params():
    gp = {}
    if os.path.exists("./inputs/global_params.json"):
        with open("./inputs/global_params.json") as f:
            gp = json.load(f)
    organism  = os.environ.get("PARAM_ORGANISM",    gp.get("organism", "unknown"))
    blast_n   = int(os.environ.get("PARAM_BLAST_SETS", str(gp.get("blast_sets", 3))))
    skip_str  = os.environ.get("PARAM_SKIP_BLAST", str(gp.get("skip_blast", False))).lower()
    skip_blast = skip_str in ("true", "1", "yes")
    return organism, blast_n, skip_blast


# ── BLAST db builder ──────────────────────────────────────────────────────────

def _build_blast_db(target_fasta, exclusion_fasta):
    db_dir = "/tmp/blast_db"
    os.makedirs(db_dir, exist_ok=True)
    combined_fasta = os.path.join(db_dir, "combined.fasta")
    db_prefix      = os.path.join(db_dir, "qpcr_db")

    if os.path.exists(db_prefix + ".nsi") or os.path.exists(db_prefix + ".nin"):
        log.info("Local BLAST db already exists: %s", db_prefix)
        return db_prefix

    with open(combined_fasta, "w") as out:
        for fa in [target_fasta, exclusion_fasta]:
            if os.path.exists(fa):
                with open(fa) as f:
                    out.write(f.read())

    cmd = ["makeblastdb", "-in", combined_fasta, "-dbtype", "nucl", "-out", db_prefix, "-title", "qpcr_db"]
    log.info("Building local BLAST db from target + exclusion sequences...")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if result.returncode != 0:
            log.error("makeblastdb failed: %s", result.stderr)
            return None
        log.info("BLAST db built: %s", db_prefix)
        return db_prefix
    except Exception as exc:
        log.error("makeblastdb error: %s", exc)
        return None


# ── BLAST runner ──────────────────────────────────────────────────────────────

def _run_blast(sequences, db_path):
    """Run blastn-short for a dict of {label: sequence}. Returns {label: [hits]}."""
    results = {label: [] for label in sequences}

    with tempfile.NamedTemporaryFile(mode="w", suffix=".fasta", delete=False) as qf:
        for label, seq in sequences.items():
            qf.write(f">{label}\n{seq}\n")
        query_path = qf.name

    out_path = query_path + ".xml"
    cmd = [
        "blastn", "-task", "blastn-short",
        "-db", db_path, "-query", query_path,
        "-out", out_path, "-outfmt", "5",
        "-word_size", "7", "-evalue", "1000",
        "-dust", "no", "-num_alignments", "10", "-num_threads", "2",
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if proc.returncode != 0:
            log.error("blastn error: %s", proc.stderr)
            return results
    except subprocess.TimeoutExpired:
        log.error("blastn timed out")
        return results
    except FileNotFoundError:
        log.error("blastn not found — is BLAST+ installed?")
        return results
    finally:
        try: os.unlink(query_path)
        except Exception: pass

    try:
        root = ET.parse(out_path).getroot()
    except Exception as exc:
        log.error("XML parse error: %s", exc)
        return results
    finally:
        try: os.unlink(out_path)
        except Exception: pass

    labels = list(sequences.keys())
    for iteration in root.findall(".//Iteration"):
        query_def = iteration.findtext("Iteration_query-def", "").strip()
        query_len = int(iteration.findtext("Iteration_query-len", "1"))
        label = next((l for l in labels if l in query_def or query_def in l), query_def)

        hits = []
        for rank, hit_el in enumerate(iteration.findall(".//Hit"), start=1):
            title     = hit_el.findtext("Hit_def", "")
            accession = hit_el.findtext("Hit_accession", "")
            hsp       = hit_el.find(".//Hsp")
            if hsp is None:
                continue
            align_len  = int(hsp.findtext("Hsp_align-len", "1"))
            identity   = int(hsp.findtext("Hsp_identity",  "0"))
            evalue     = float(hsp.findtext("Hsp_evalue",  "999"))
            bitscore   = float(hsp.findtext("Hsp_bit-score","0"))
            hits.append({
                "rank": rank, "accession": accession, "title": title,
                "identity":  round(identity / align_len * 100, 1) if align_len else 0,
                "coverage":  round(align_len / query_len * 100, 1) if query_len else 0,
                "evalue": evalue, "bitscore": bitscore,
            })
            if rank >= 10:
                break
        results[label] = hits
    return results


def _evaluate_specificity(hits, target_organism):
    """Returns True if all significant hits map back to the target organism."""
    parts = target_organism.lower().split()
    genus, species = (parts[0] if parts else ""), (parts[1] if len(parts) > 1 else "")

    for hit in hits:
        if hit["evalue"] > 0.01:
            hit["is_target"] = True
            continue
        title = hit["title"].lower()
        species_in = species and species in title
        genus_in   = genus in title
        if not genus_in and not species_in:
            hit["is_target"] = False
        elif species_in:
            hit["is_target"] = True
        elif genus_in:
            idx = title.find(genus)
            next_w = title[idx + len(genus):].split()
            nw = next_w[0].strip(".,;") if next_w else ""
            hit["is_target"] = not (nw and nw != species and nw not in ("sp", "sp.", "spp", "spp."))
        else:
            hit["is_target"] = True

    return all(h["is_target"] for h in hits), hits


def _validate_primer_sets(primer_sets, organism, db_path, max_sets):
    cache = {}
    to_blast_sets = primer_sets[:max_sets]

    for ps in to_blast_sets:
        seqs = {
            "fwd": ps["fwd_seq"],
            "rev": ps["rev_seq"],
            "prb": ps["probe_seq"],
        }
        uncached = {k: v for k, v in seqs.items() if v and v not in cache}
        if uncached:
            hit_map = _run_blast(uncached, db_path)
            for label, seq in uncached.items():
                hits = hit_map.get(label, [])
                specific, hits = _evaluate_specificity(hits, organism)
                cache[seq] = {"specific": specific, "hits": hits}

        blast_results = {}
        for key, seq in seqs.items():
            if seq and seq in cache:
                blast_results[key] = cache[seq]

        ps["blast_fwd"]   = blast_results.get("fwd")
        ps["blast_rev"]   = blast_results.get("rev")
        ps["blast_probe"] = blast_results.get("prb")
        checks = [v for v in blast_results.values() if v is not None]
        ps["blast_pass"]  = all(v["specific"] for v in checks) if checks else None
        log.info("  Pair #%d: blast_pass=%s", ps["pair_index"], ps["blast_pass"])

    return primer_sets


# ── Report writers ────────────────────────────────────────────────────────────

def _set_to_dict(ps, organism):
    def _blast_summary(br):
        if br is None: return {"status": "not_run"}
        return {
            "status":   "ok",
            "specific": br["specific"],
            "top_hits": [
                {k: h[k] for k in ("rank", "accession", "title", "identity", "coverage", "evalue", "is_target")}
                for h in br["hits"][:5]
            ],
        }
    return {
        "pair_index": ps["pair_index"],
        "roi": {
            "source_id":        ps["roi"]["record_id"],
            "start":            ps["roi"]["start"],
            "end":              ps["roi"]["end"],
            "uniqueness_score": ps["roi"]["uniqueness_score"],
        },
        "primers": {
            "forward": {"sequence": ps["fwd_seq"], "tm": round(ps["fwd_tm"], 2), "gc": round(ps["fwd_gc"], 1), "length": len(ps["fwd_seq"])},
            "reverse": {"sequence": ps["rev_seq"], "tm": round(ps["rev_tm"], 2), "gc": round(ps["rev_gc"], 1), "length": len(ps["rev_seq"])},
            "probe":   {"sequence": ps["probe_seq"], "tm": round(ps["probe_tm"], 2), "gc": round(ps["probe_gc"], 1), "length": len(ps.get("probe_seq") or "")},
        },
        "amplicon_size":      ps["amplicon_size"],
        "primer3_penalty":    round(ps["penalty"], 4),
        "passed_constraints": ps["passed_constraints"],
        "constraint_notes":   ps["constraint_notes"],
        "blast": {
            "forward_primer": _blast_summary(ps.get("blast_fwd")),
            "reverse_primer": _blast_summary(ps.get("blast_rev")),
            "probe":          _blast_summary(ps.get("blast_probe")),
            "overall_pass":   ps.get("blast_pass"),
        },
    }


def _write_json(primer_sets, organism, outdir):
    path = os.path.join(outdir, "results.json")
    payload = {
        "pipeline":        "qpcr-primer-design",
        "version":         "1.0.0",
        "timestamp":       datetime.utcnow().isoformat() + "Z",
        "target_organism": organism,
        "n_sets":          len(primer_sets),
        "primer_sets":     [_set_to_dict(ps, organism) for ps in primer_sets],
    }
    with open(path, "w") as f:
        json.dump(payload, f, indent=2)
    log.info("JSON report → %s", path)
    return path


def _write_tsv(primer_sets, outdir):
    path = os.path.join(outdir, "results.tsv")
    headers = [
        "pair_index", "roi_id", "roi_start", "roi_end", "roi_uniqueness",
        "fwd_seq", "fwd_len", "fwd_tm", "fwd_gc",
        "rev_seq", "rev_len", "rev_tm", "rev_gc",
        "probe_seq", "probe_len", "probe_tm", "probe_gc",
        "amplicon_size", "primer3_penalty",
        "passed_constraints", "constraint_notes",
        "blast_fwd_specific", "blast_rev_specific", "blast_probe_specific", "blast_overall_pass",
    ]

    def _spec(br):
        if br is None: return "N/A"
        return str(br.get("specific", "N/A"))

    rows = []
    for ps in primer_sets:
        rows.append([
            ps["pair_index"],
            ps["roi"]["record_id"], ps["roi"]["start"], ps["roi"]["end"],
            ps["roi"]["uniqueness_score"],
            ps["fwd_seq"], len(ps["fwd_seq"]), round(ps["fwd_tm"], 2), round(ps["fwd_gc"], 1),
            ps["rev_seq"], len(ps["rev_seq"]), round(ps["rev_tm"], 2), round(ps["rev_gc"], 1),
            ps.get("probe_seq", ""), len(ps.get("probe_seq") or ""), round(ps.get("probe_tm", 0), 2), round(ps.get("probe_gc", 0), 1),
            ps["amplicon_size"], round(ps["penalty"], 4),
            ps["passed_constraints"], "; ".join(ps["constraint_notes"]),
            _spec(ps.get("blast_fwd")), _spec(ps.get("blast_rev")), _spec(ps.get("blast_probe")),
            str(ps.get("blast_pass")),
        ])

    with open(path, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in rows:
            f.write("\t".join(str(v) for v in row) + "\n")

    log.info("TSV report → %s", path)
    return path


def _write_txt(primer_sets, organism, skip_blast, outdir):
    path = os.path.join(outdir, "report.txt")
    lines = [
        "=" * 70,
        "  qPCR Primer/Probe Pipeline Report",
        f"  Target: {organism}",
        f"  Generated: {datetime.utcnow().strftime('%Y-%m-%d %H:%M UTC')}",
        "=" * 70, "",
    ]

    passed = [ps for ps in primer_sets
              if ps["passed_constraints"] and (skip_blast or ps.get("blast_pass"))]
    lines.append(f"Total primer sets designed : {len(primer_sets)}")
    lines.append(f"Passed ALL constraints     : {len(passed)}")
    if not skip_blast:
        lines.append(f"BLAST-validated            : {sum(1 for ps in primer_sets if ps.get('blast_pass'))}")
    lines.append("")

    for rank, ps in enumerate(primer_sets, start=1):
        blast_ok = ps.get("blast_pass") if ps.get("blast_pass") is not None else True
        status = "PASS" if (ps["passed_constraints"] and (skip_blast or blast_ok)) else "FAIL"
        lines.append("─" * 70)
        lines.append(f"  Rank #{rank}  (P3 pair #{ps['pair_index']})  [{status}]  "
                     f"ROI {ps['roi']['record_id']}:{ps['roi']['start']}-{ps['roi']['end']}")
        lines.append(f"  ROI uniqueness: {ps['roi']['uniqueness_score']:.3f}")
        lines.append("")
        lines.append(f"  FWD  5'-{ps['fwd_seq']}-3'")
        lines.append(f"       Tm={ps['fwd_tm']:.1f}°C  GC={ps['fwd_gc']:.1f}%  len={len(ps['fwd_seq'])}nt")
        lines.append(f"  REV  5'-{ps['rev_seq']}-3'")
        lines.append(f"       Tm={ps['rev_tm']:.1f}°C  GC={ps['rev_gc']:.1f}%  len={len(ps['rev_seq'])}nt")
        if ps.get("probe_seq"):
            lines.append(f"  PRB  5'-{ps['probe_seq']}-3'")
            lines.append(f"       Tm={ps['probe_tm']:.1f}°C  GC={ps['probe_gc']:.1f}%  len={len(ps['probe_seq'])}nt")
            delta = ps["probe_tm"] - (ps["fwd_tm"] + ps["rev_tm"]) / 2
            lines.append(f"       ΔTm vs primers = {delta:+.1f}°C")
        lines.append(f"  Amplicon: {ps['amplicon_size']} nt  |  P3 penalty: {ps['penalty']:.4f}")
        if ps["constraint_notes"]:
            lines.append("  Constraint issues:")
            for note in ps["constraint_notes"]:
                lines.append(f"    • {note}")
        if not skip_blast and ps.get("blast_pass") is not None:
            label = "SPECIFIC" if ps["blast_pass"] else "NOT SPECIFIC"
            lines.append(f"  BLAST: {label}")
            for bl_label, br in [("FWD", ps.get("blast_fwd")), ("REV", ps.get("blast_rev")), ("PRB", ps.get("blast_probe"))]:
                if br and br.get("hits"):
                    for h in br["hits"][:3]:
                        flag = "+" if h.get("is_target") else "x"
                        lines.append(f"    {bl_label}  [{flag}] {h['title'][:55]:<55} "
                                     f"id={h['identity']:.1f}% e={h['evalue']:.1e}")
        lines.append("")

    with open(path, "w") as f:
        f.write("\n".join(lines))

    log.info("TXT report → %s", path)
    return path


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    organism, blast_n, skip_blast = _load_params()
    os.makedirs("./outputs", exist_ok=True)

    with open("./inputs/primer_sets.json") as f:
        data = json.load(f)
    primer_sets = data.get("primer_sets", [])

    if not primer_sets:
        log.error("primer_sets.json contains no primer sets")
        sys.exit(1)

    log.info("Loaded %d primer set(s) for organism: %s", len(primer_sets), organism)

    if not skip_blast:
        db_path = _build_blast_db("./inputs/target.fasta", "./inputs/exclusion.fasta")
        if db_path:
            log.info("Running BLAST validation on top %d set(s)...", blast_n)
            primer_sets = _validate_primer_sets(primer_sets, organism, db_path, blast_n)
        else:
            log.warning("BLAST db build failed — skipping BLAST validation")
            skip_blast = True
    else:
        log.info("BLAST validation skipped (skip_blast=true)")

    _write_json(primer_sets, organism, "./outputs")
    _write_tsv(primer_sets, "./outputs")
    _write_txt(primer_sets, organism, skip_blast, "./outputs")

    fully_passed = [ps for ps in primer_sets
                    if ps["passed_constraints"] and (skip_blast or ps.get("blast_pass"))]
    log.info("Pipeline complete. %d/%d sets fully passed.", len(fully_passed), len(primer_sets))

    if not fully_passed:
        log.warning("No sets passed all filters. Review report.txt and relax constraints if needed.")
        sys.exit(2)


if __name__ == "__main__":
    main()
