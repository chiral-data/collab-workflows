# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/ doc: https://pocketeer.readthedocs.io/en/latest/

import json
import os
import urllib.request

import imageio.v2 as imageio
import numpy as np
import pocketeer as pt
from pocketeer.core.types import AlphaSphere


def load_pockets_json(json_path: str) -> list[pt.Pocket]:
    """Load pockets from a JSON file created by pt.write_pockets_json()."""
    with open(json_path, "r") as f:
        data = json.load(f)

    pockets = []
    for p in data:
        # Reconstruct AlphaSphere objects
        spheres = [
            AlphaSphere(
                sphere_id=s["sphere_id"],
                center=np.array(s["center"]),
                radius=s["radius"],
                mean_sasa=s["mean_sasa"],
                atom_indices=s["atom_indices"],
            )
            for s in p["spheres"]
        ]
        # Create Pocket object
        pocket = pt.Pocket(
            pocket_id=p["pocket_id"],
            spheres=spheres,
            centroid=np.array(p["centroid"]),
            volume=p["volume"],
            score=p["score"],
        )
        pockets.append(pocket)

    return pockets


print("Initializing PyMOL...", flush=True)
from pymol import cmd

print("PyMOL initialized.", flush=True)

# Visualization representation: "cartoon" or "surface"
representation = "surface"

pdb_code = "4tos"

# Download the pdb file for demonstration
print(f"\nDownloading PDB file: {pdb_code}...", flush=True)
pdb_filename = f"{pdb_code.upper()}.pdb"
url = f"https://files.rcsb.org/download/{pdb_code.upper()}.pdb"
urllib.request.urlretrieve(url, pdb_filename)
print(f"Downloaded {pdb_filename}", flush=True)

# Load structure
print("\nLoading structure...", flush=True)
atomarray = pt.load_structure(pdb_filename)

# Detect pockets
print("Detecting pockets (this may take a moment)...", flush=True)
pockets = pt.find_pockets(atomarray)

# Display results
print(f"\nFound {len(pockets)} pockets:")
for pocket in pockets[:5]:  # Show top 5
    print(
        f"  Pocket {pocket.pocket_id}: score={pocket.score:.2f}, "
        f"volume={pocket.volume:.1f} Å³, "
        f"spheres={pocket.n_spheres}"
    )

# Create visualization
viewer = pt.view_pockets(atomarray, pockets)

# Save visualization to HTML file (works in non-notebook environments)
html_content = viewer._make_html()
with open("pocket_visualization.html", "w") as f:
    f.write(html_content)
print("\nVisualization saved to pocket_visualization.html")

# Generate rotating GIF using PyMOL
print("\n" + "=" * 50, flush=True)
print("Generating rotating GIF...", flush=True)
print("=" * 50, flush=True)

# Initialize PyMOL in headless mode
print("Setting up PyMOL scene...", flush=True)
cmd.reinitialize()

# Enable multi-threading for faster ray tracing
import multiprocessing

num_cpus = multiprocessing.cpu_count()
cmd.set("max_threads", num_cpus)
print(f"Using {num_cpus} CPU threads for ray tracing", flush=True)

cmd.load(pdb_filename, "structure")

# Style the protein based on representation choice
print(f"Using '{representation}' representation", flush=True)
cmd.hide("everything", "structure")
cmd.show(representation, "structure")
cmd.color("cyan", "structure")

if representation == "cartoon":
    cmd.set("cartoon_fancy_helices", 1)
elif representation == "surface":
    cmd.set("surface_quality", 1)
    cmd.set("transparency", 0.3, "structure")

cmd.set("ray_opaque_background", 1)
cmd.bg_color("white")

# Highlight top pockets by rendering all alpha spheres that define each pocket
colors = ["green", "yellow", "cyan", "magenta", "orange"]
print(f"Adding {len(pockets[:5])} pockets to visualization...", flush=True)

for i, pocket in enumerate(pockets[:5]):
    pocket_name = f"pocket_{i}"
    # Create a pseudoatom for each alpha sphere in the pocket
    for j, sphere in enumerate(pocket.spheres):
        atom_name = f"{pocket_name}_s{j}"
        cmd.pseudoatom(atom_name, pos=list(sphere.center), vdw=sphere.radius)

    # Group all spheres of this pocket and show as surface
    cmd.group(pocket_name, f"{pocket_name}_s*")
    cmd.show("spheres", pocket_name)
    cmd.color(colors[i], pocket_name)
    # Make pocket spheres slightly transparent
    cmd.set("sphere_transparency", 0.3, pocket_name)
    print(f"  Pocket {i+1}: {len(pocket.spheres)} spheres added", flush=True)

cmd.zoom("all", buffer=5)
print("Scene setup complete.", flush=True)

# Generate frames for rotation
step = 10  # degrees per frame
total_frames = 360 // step
images = []
frame_dir = "/tmp/frames"
os.makedirs(frame_dir, exist_ok=True)

print(f"\nRendering {total_frames} frames (using OpenGL draw)...", flush=True)
for i, angle in enumerate(range(0, 360, step)):
    cmd.rotate("y", float(step))
    # Use draw() instead of ray() for faster OpenGL rendering
    cmd.draw(512, 512)
    frame_path = f"{frame_dir}/frame_{angle:03d}.png"
    cmd.png(frame_path, ray=0)  # ray=0 to use the drawn image
    images.append(imageio.imread(frame_path))
    print(f"  Frame {i + 1}/{total_frames} ({angle}°) complete", flush=True)

# Create GIF
print("\nAssembling GIF...", flush=True)
gif_path = f"protein_pockets_rotation_{representation}.gif"
imageio.mimsave(gif_path, images, duration=0.1, loop=0)
print(f"\nRotating GIF saved to {gif_path}", flush=True)

# Cleanup temporary frames
for f in os.listdir(frame_dir):
    os.remove(os.path.join(frame_dir, f))
os.rmdir(frame_dir)
