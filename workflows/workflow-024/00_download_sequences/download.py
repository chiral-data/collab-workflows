#!/usr/bin/env python3
"""
Node 00: Download or accept FASTA sequences for qPCR primer design.

Checks for user-provided FASTAs in ./inputs/; if both target.fasta and
exclusion.fasta are present, copies them through (offline mode).
Otherwise downloads from NCBI Entrez based on the organism parameter.
"""

import json
import logging
import os
import shutil
import sys
import time

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)-8s %(message)s")
log = logging.getLogger("download")


def _load_params():
    gp = {}
    # Probe multiple locations: inputs/ (standard Chiral forwarded file), then
    # the workflow root mounted at the container CWD, then the local dev path.
    for gp_path in ("./inputs/global_params.json", "./global_params.json", "../global_params.json"):
        if os.path.exists(gp_path):
            with open(gp_path) as f:
                gp = json.load(f)
            break

    return {
        "organism":          os.environ.get("PARAM_ORGANISM") or      gp.get("organism", ""),
        "email":             os.environ.get("PARAM_EMAIL",             gp.get("email", "user@example.com")),
        "ncbi_api_key":      os.environ.get("PARAM_NCBI_API_KEY",      gp.get("ncbi_api_key", "")),
        "target_fasta":      os.environ.get("PARAM_TARGET_FASTA",      gp.get("target_fasta", "")),
        "target_gene":       os.environ.get("PARAM_TARGET_GENE",       gp.get("target_gene", "")),
        "max_target_seqs":   int(os.environ.get("PARAM_MAX_TARGET_SEQS",  str(gp.get("max_target_seqs", 3)))),
        "max_relative_seqs": int(os.environ.get("PARAM_MAX_RELATIVE_SEQS", str(gp.get("max_relative_seqs", 5)))),
    }


def _configure_entrez(email, api_key=""):
    from Bio import Entrez
    Entrez.email  = email
    Entrez.tool   = "qpcr_pipeline_silva"
    if api_key:
        Entrez.api_key = api_key


def _search_accessions(organism, max_seqs, email, api_key):
    from Bio import Entrez
    _configure_entrez(email, api_key)
    queries = [
        f'"{organism}"[Organism] AND "complete genome"[Title] AND refseq[filter]',
        f'"{organism}"[Organism] AND "complete sequence"[Title]',
        f'"{organism}"[Organism]',
    ]
    for q in queries:
        try:
            h = Entrez.esearch(db="nucleotide", term=q, retmax=max_seqs, usehistory="y")
            rec = Entrez.read(h); h.close()
            ids = rec.get("IdList", [])
            if ids:
                log.info("Found %d accession(s) with query: %s", len(ids), q)
                return ids
        except Exception as exc:
            log.warning("Search failed: %s", exc)
            time.sleep(1)
    return []


def _fetch_sequences(ids, email, api_key):
    from Bio import Entrez, SeqIO
    _configure_entrez(email, api_key)
    records = []
    for i in range(0, len(ids), 5):
        batch = ids[i:i + 5]
        try:
            h = Entrez.efetch(db="nucleotide", id=",".join(batch), rettype="fasta", retmode="text")
            for rec in SeqIO.parse(h, "fasta-blast"):
                if len(rec.seq) > 0:
                    records.append(rec)
            h.close()
            time.sleep(0.4)
        except Exception as exc:
            log.error("Fetch failed for %s: %s", batch, exc)
            time.sleep(2)
    return records


def _search_gene_sequences(organism, gene, max_seqs, email, api_key):
    from Bio import Entrez
    _configure_entrez(email, api_key)
    queries = [
        f'"{organism}"[Organism] AND "{gene}"[Gene] AND refseq[filter]',
        f'"{organism}"[Organism] AND "{gene}"[Title]',
        f'"{organism}"[Organism] AND "{gene}"[All Fields]',
    ]
    for q in queries:
        try:
            h = Entrez.esearch(db="nucleotide", term=q, retmax=max_seqs)
            rec = Entrez.read(h); h.close()
            ids = rec.get("IdList", [])
            if ids:
                recs = _fetch_sequences(ids, email, api_key)
                if recs:
                    return recs
        except Exception as exc:
            log.warning("Gene search failed: %s", exc)
            time.sleep(1)
    return []


def _get_close_relative_gene_sequences(organism, gene, n, email, api_key):
    from Bio import Entrez
    parts = organism.strip().split()
    if len(parts) < 2:
        return []
    genus, species_q = parts[0], " ".join(parts[:2])
    queries = [
        f'"{genus}"[Organism] NOT "{species_q}"[Organism] AND "{gene}"[Gene] AND refseq[filter]',
        f'"{genus}"[Organism] NOT "{species_q}"[Organism] AND "{gene}"[Title]',
    ]
    _configure_entrez(email, api_key)
    for q in queries:
        try:
            h = Entrez.esearch(db="nucleotide", term=q, retmax=n)
            rec = Entrez.read(h); h.close()
            ids = rec.get("IdList", [])
            if ids:
                recs = _fetch_sequences(ids, email, api_key)
                if recs:
                    return recs
        except Exception as exc:
            log.warning("Relative gene search failed: %s", exc)
            time.sleep(1)
    return []


def _get_close_relatives(organism, n, email, api_key):
    from Bio import Entrez
    parts = organism.strip().split()
    if len(parts) < 2:
        return []
    genus, species_q = parts[0], " ".join(parts[:2])
    query = f'"{genus}"[Organism] NOT "{species_q}"[Organism] AND "complete genome"[Title]'
    _configure_entrez(email, api_key)
    try:
        h = Entrez.esearch(db="nucleotide", term=query, retmax=n)
        rec = Entrez.read(h); h.close()
        ids = rec.get("IdList", [])
    except Exception as exc:
        log.warning("Relative search failed: %s", exc)
        return []
    return _fetch_sequences(ids, email, api_key) if ids else []


def main():
    params       = _load_params()
    organism     = params["organism"]
    email        = params["email"]
    api_key      = params["ncbi_api_key"]
    target_fasta = params["target_fasta"]
    gene         = params["target_gene"]
    max_tgt      = params["max_target_seqs"]
    max_rel      = params["max_relative_seqs"]

    os.makedirs("./outputs", exist_ok=True)

    # Always forward global_params.json so downstream nodes can read it
    gp_src = "./inputs/global_params.json"
    if os.path.exists(gp_src):
        shutil.copy(gp_src, "./outputs/global_params.json")
    else:
        # Write defaults so downstream nodes have something to read
        with open("./outputs/global_params.json", "w") as f:
            json.dump({"organism": organism, "email": email}, f, indent=2)

    excl_in = "./inputs/exclusion.fasta"

    # target_fasta mode: use a named FASTA from inputs/ as the target sequence
    if target_fasta:
        tgt_path = f"./inputs/{target_fasta}"
        if not os.path.exists(tgt_path):
            log.error("target_fasta '%s' not found in inputs/ — check the filename in global_params.json", target_fasta)
            sys.exit(1)
        log.info("Using user-provided target FASTA: %s", target_fasta)
        shutil.copy(tgt_path, "./outputs/target.fasta")

        # Exclusion: use inputs/exclusion.fasta if present, otherwise fetch from NCBI
        if os.path.exists(excl_in):
            log.info("Using user-provided exclusion.fasta from inputs/")
            shutil.copy(excl_in, "./outputs/exclusion.fasta")
        elif organism:
            log.info("Fetching exclusion sequences from NCBI for organism: '%s'", organism)
            from Bio import SeqIO
            if gene:
                excl_recs = _get_close_relative_gene_sequences(organism, gene, max_rel, email, api_key)
            else:
                excl_recs = _get_close_relatives(organism, max_rel, email, api_key)
            total_bp = sum(len(r.seq) for r in excl_recs)
            if total_bp > 10_000_000 and len(excl_recs) > 5:
                log.warning("Large exclusion set (%.1f MB) — capping to 5 sequences", total_bp / 1e6)
                excl_recs = excl_recs[:5]
            log.info("Fetched %d exclusion sequence(s)", len(excl_recs))
            with open("./outputs/exclusion.fasta", "w") as fh:
                SeqIO.write(excl_recs, fh, "fasta")
        else:
            log.warning("No exclusion source — writing empty exclusion.fasta (all ROIs will score 1.0)")
            open("./outputs/exclusion.fasta", "w").close()

        log.info("Download complete.")
        return

    # Offline / no-organism mode: skip NCBI entirely when organism is not set
    tgt_in = "./inputs/target.fasta"
    if not organism:
        if not os.path.exists(tgt_in):
            log.error(
                "No organism specified and no target.fasta found in inputs/ — "
                "provide 'organism' in global_params.json for NCBI download mode, "
                "or place target.fasta (and optionally exclusion.fasta) in input_files/ "
                "for offline mode."
            )
            sys.exit(1)
        log.info("No organism specified — using local FASTAs from inputs/ (offline mode)")
        shutil.copy(tgt_in, "./outputs/target.fasta")
        if os.path.exists(excl_in):
            shutil.copy(excl_in, "./outputs/exclusion.fasta")
            log.info("Copied target.fasta and exclusion.fasta to outputs/")
        else:
            open("./outputs/exclusion.fasta", "w").close()
            log.warning("No exclusion.fasta found in inputs/ — writing empty exclusion.fasta (all ROIs will score 1.0)")
        return

    log.info("Downloading sequences for organism: '%s'", organism)
    from Bio import SeqIO

    # Target sequences
    if gene:
        log.info("Gene-targeted mode: fetching '%s' gene sequences", gene)
        target_recs = _search_gene_sequences(organism, gene, max_tgt, email, api_key)
    else:
        ids = _search_accessions(organism, max_tgt, email, api_key)
        if not ids:
            log.error("No sequences found for '%s' — check organism name and NCBI connectivity", organism)
            sys.exit(1)
        target_recs = _fetch_sequences(ids, email, api_key)

    if not target_recs:
        log.error("Failed to fetch any target sequences")
        sys.exit(1)

    log.info("Fetched %d target sequence(s)", len(target_recs))
    with open("./outputs/target.fasta", "w") as fh:
        SeqIO.write(target_recs, fh, "fasta")

    # Exclusion (close-relative) sequences
    if gene:
        excl_recs = _get_close_relative_gene_sequences(organism, gene, max_rel, email, api_key)
    else:
        excl_recs = _get_close_relatives(organism, max_rel, email, api_key)

    # Cap large exclusion sets to keep runtime reasonable
    total_bp = sum(len(r.seq) for r in excl_recs)
    if total_bp > 10_000_000 and len(excl_recs) > 5:
        log.warning("Large exclusion set (%.1f MB) — capping to 5 sequences", total_bp / 1e6)
        excl_recs = excl_recs[:5]

    log.info("Fetched %d exclusion sequence(s)", len(excl_recs))
    with open("./outputs/exclusion.fasta", "w") as fh:
        SeqIO.write(excl_recs, fh, "fasta")

    log.info("Download complete.")


if __name__ == "__main__":
    main()
