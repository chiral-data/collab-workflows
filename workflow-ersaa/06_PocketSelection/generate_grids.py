#!/usr/bin/env python3

import json
from pathlib import Path

if __name__ == "__main__":
    pockets_file = Path("/workspace/input/pockets.json")
    grids_file = Path("/workspace/out/grids.json")

    with open(pockets_file) as f:
        pockets = json.load(f)

    # Convert pockets → vina boxes
    vina_boxes = []
    for p in pockets:
        vina_boxes.append({
            "center_x": p["center_x"],
            "center_y": p["center_y"],
            "center_z": p["center_z"],
            "size_x": 20,
            "size_y": 20,
            "size_z": 20
        })

    with open(grids_file, "w") as f:
        json.dump(vina_boxes, f, indent=2)

    print("[Node6] Grid boxes generated")
