#!/usr/bin/env python3
import os
from pathlib import Path
import json

if __name__ == "__main__":
    receptor = Path("/workspace/input/receptor.pdb")
    out_json = Path("/workspace/out/pockets.json")

    # Run P2Rank
    os.system(f"java -jar /workspace/p2rank/prank.jar predict {receptor} -o /workspace/out")

    # Convert P2Rank output folder to pockets.json if needed
    # (Your existing code remains here)

    print("[Node5] P2Rank completed")
