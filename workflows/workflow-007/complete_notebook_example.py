# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/ doc: https://pocketeer.readthedocs.io/en/latest/

import os
import urllib.request

import imageio.v2 as imageio
import pocketeer as pt
from pymol import cmd

pdb_code = "4tos"

# Download the pdb file for demonstration
pdb_filename = f"{pdb_code.upper()}.pdb"
url = f"https://files.rcsb.org/download/{pdb_code.upper()}.pdb"
urllib.request.urlretrieve(url, pdb_filename)

# Load structure
atomarray = pt.load_structure(pdb_filename)

# Detect pockets
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
print("\nGenerating rotating GIF...")

# Initialize PyMOL in headless mode
cmd.reinitialize()
cmd.load(pdb_filename, "protein")

# Style the protein
cmd.show("cartoon", "protein")
cmd.color("cyan", "protein")
cmd.set("cartoon_fancy_helices", 1)
cmd.set("ray_opaque_background", 1)
cmd.bg_color("white")

# Highlight top pockets with spheres
for i, pocket in enumerate(pockets[:5]):
    # Get pocket centroid coordinates
    centroid = pocket.centroid
    pocket_name = f"pocket_{i}"
    # Create a pseudoatom at pocket centroid
    cmd.pseudoatom(pocket_name, pos=list(centroid))
    cmd.show("spheres", pocket_name)
    # Color pockets differently
    colors = ["red", "orange", "yellow", "green", "magenta"]
    cmd.color(colors[i], pocket_name)
    cmd.set("sphere_scale", 2.0, pocket_name)

cmd.zoom("all", buffer=5)

# Generate frames for rotation
step = 10  # degrees per frame
images = []
frame_dir = "/tmp/frames"
os.makedirs(frame_dir, exist_ok=True)

for angle in range(0, 360, step):
    cmd.rotate("y", float(step))
    cmd.ray(512, 512)
    frame_path = f"{frame_dir}/frame_{angle:03d}.png"
    cmd.png(frame_path)
    images.append(imageio.imread(frame_path))

# Create GIF
gif_path = "protein_pockets_rotation.gif"
imageio.mimsave(gif_path, images, duration=0.1, loop=0)
print(f"Rotating GIF saved to {gif_path}")

# Cleanup temporary frames
for f in os.listdir(frame_dir):
    os.remove(os.path.join(frame_dir, f))
os.rmdir(frame_dir)
