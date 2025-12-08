# Part 3: Generate visualization output (HTML/GIF)
# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/

import json
import os

import imageio.v2 as imageio
import numpy as np
import pocketeer as pt
from pocketeer.core.types import AlphaSphere

# =============================================================================
# CONFIGURATION FROM ENVIRONMENT VARIABLES
# =============================================================================

# PDB ID (from global workflow parameter)
pdb_id = os.environ.get("PARAM_PDB_ID", "4TOS")
pdb_filename = f"{pdb_id.upper()}.pdb"
pockets_json = "pockets.json"

# Visualization parameters (from job parameters)
# Pocket visualization style: "filled_surfaces" or "single_sphere"
pocket_style = os.environ.get("PARAM_POCKET_STYLE", "filled_surfaces")

# Rendering method: "ray" (higher quality, slower) or "draw" (faster)
render_method = os.environ.get("PARAM_RENDER_METHOD", "draw")

# Protein visualization representation: "cartoon" or "surface"
representation = os.environ.get("PARAM_REPRESENTATION", "surface")

# Output format: "html", "gif", or "both"
output_format = os.environ.get("PARAM_OUTPUT_FORMAT", "html")

# =============================================================================


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


# Print configuration
print(f"Visualization configuration:", flush=True)
print(f"  pocket_style={pocket_style}", flush=True)
print(f"  render_method={render_method}", flush=True)
print(f"  representation={representation}", flush=True)
print(f"  output_format={output_format}", flush=True)

# Load structure and pockets
print("Loading structure...", flush=True)
atomarray = pt.load_structure(pdb_filename)

print("Loading pockets from JSON...", flush=True)
pockets = load_pockets_json(pockets_json)
print(f"Loaded {len(pockets)} pockets", flush=True)

# Generate HTML visualization if requested
if output_format in ("html", "both"):
    # Create visualization
    viewer = pt.view_pockets(atomarray, pockets)

    # Save visualization to HTML file (works in non-notebook environments)
    html_content = viewer._make_html()
    with open("pocket_visualization.html", "w") as f:
        f.write(html_content)
    print("\nVisualization saved to pocket_visualization.html")

# Generate rotating GIF using PyMOL if requested
if output_format in ("gif", "both"):
    print("Initializing PyMOL...", flush=True)
    from pymol import cmd

    print("PyMOL initialized.", flush=True)

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

    # Highlight top pockets based on pocket_style
    colors = ["red", "orange", "yellow", "green", "magenta"]
    print(
        f"Adding {len(pockets[:5])} pockets using '{pocket_style}' style...", flush=True
    )

    if pocket_style == "filled_surfaces":
        # Show all alpha spheres that define each pocket
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
            print(f"  Pocket {i + 1}: {len(pocket.spheres)} spheres added", flush=True)

    elif pocket_style == "single_sphere":
        # Show a single sphere at each pocket's centroid
        for i, pocket in enumerate(pockets[:5]):
            centroid = pocket.centroid
            pocket_name = f"pocket_{i}"
            # Create a pseudoatom at pocket centroid
            cmd.pseudoatom(pocket_name, pos=list(centroid))
            cmd.show("spheres", pocket_name)
            cmd.color(colors[i], pocket_name)
            cmd.set("sphere_scale", 2.0, pocket_name)
            print(f"  Pocket {i + 1}: centroid sphere added", flush=True)

    cmd.zoom("all", buffer=5)
    print("Scene setup complete.", flush=True)

    # Generate frames for rotation
    step = 10  # degrees per frame
    total_frames = 360 // step
    images = []
    frame_dir = "/tmp/frames"
    os.makedirs(frame_dir, exist_ok=True)

    if render_method == "ray":
        print(
            f"\nRendering {total_frames} frames (ray tracing - this takes time)...",
            flush=True,
        )
    else:
        print(f"\nRendering {total_frames} frames (using OpenGL draw)...", flush=True)

    for i, angle in enumerate(range(0, 360, step)):
        cmd.rotate("y", float(step))

        if render_method == "ray":
            cmd.ray(512, 512)
            frame_path = f"{frame_dir}/frame_{angle:03d}.png"
            cmd.png(frame_path)
        else:  # draw
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
