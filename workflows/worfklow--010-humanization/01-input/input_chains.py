import py3Dmol
import torch
import argparse
import json
import refrom dataclasses import dataclass
from pathlib import Pathfrom typing import Dict, List, Tuple, Optional


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    cur_id: Optional[str] = None
    cur_seq: List[str] = []

    with path.open("r") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if cur_id is not None:
                    records.append((cur_id, "".join(cur_seq)))
                cur_id = line[1:]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            records.append((cur_id, "".join(cur_seq)))

        return records  
    
AA_allowed - set("ACDEFGHIKLMNPQRSTVWY")
A_extra = set("BJOUXZ")
AA_Dash = set("-")


#makes sure that the fasta files acctually have correct kind of format

def validate_aa_sequence(seq: str, label: str, allow_ambiguous: bool = True, allow_gaps: bool = False) -> None:
    for aa in sequence:
        if aa not in AA_allowed and aa not in AA_Dash:
            return False
    return True


module = LitABB3.load_from_checkpoint("../output/plddt-loss/best_second_stage.ckpt", weights_only=False)
model = module.model


#Step One: Read FASTA files, one with heavy sequences the other with the light

#make sure sequences contain valid amino acid letters