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

# Input/Output directories (silva 0.4.0+)
input_dir = "inputs"
output_dir = "outputs"

# Input files
pdb_filename = f"{pdb_id.upper()}.pdb"
pdb_path = os.path.join(input_dir, pdb_filename)
pockets_json = os.path.join(input_dir, "pockets.json")

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

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Load structure and pockets
print(f"Loading structure from {pdb_path}...", flush=True)
atomarray = pt.load_structure(pdb_path)

print("Loading pockets from JSON...", flush=True)
pockets = load_pockets_json(pockets_json)
print(f"Loaded {len(pockets)} pockets", flush=True)

# Generate HTML visualization if requested
if output_format in ("html", "both"):
    # Configure visualization based on representation option
    receptor_cartoon = representation == "cartoon"
    receptor_surface = representation == "surface"

    # Configure sphere scale based on pocket_style
    # For single_sphere style, use larger spheres; for filled_surfaces, use normal scale
    sphere_scale = 2.0 if pocket_style == "single_sphere" else 1.0

    print(f"\nGenerating HTML with:", flush=True)
    print(f"  receptor_cartoon={receptor_cartoon}, receptor_surface={receptor_surface}", flush=True)
    print(f"  sphere_scale={sphere_scale} (pocket_style={pocket_style})", flush=True)

    # Create visualization with options
    viewer = pt.view_pockets(
        atomarray,
        pockets,
        receptor_cartoon=receptor_cartoon,
        receptor_surface=receptor_surface,
        sphere_scale=sphere_scale,
    )

    # Save visualization to HTML file (works in non-notebook environments)
    viewer_html = viewer._make_html()

    # Wrap in a complete HTML document with full-page layout and controls
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pocket Visualization - {pdb_id}</title>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        html, body {{
            height: 100%;
            width: 100%;
            overflow: hidden;
        }}
        body {{
            display: flex;
            flex-direction: column;
            background-color: #f5f5f5;
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
        }}
        header {{
            padding: 12px 20px;
            background: white;
            border-bottom: 1px solid #e0e0e0;
            box-shadow: 0 1px 3px rgba(0, 0, 0, 0.1);
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 12px;
        }}
        h1 {{
            color: #333;
            font-size: 1.25rem;
            font-weight: 500;
        }}
        .controls {{
            display: flex;
            gap: 8px;
            align-items: center;
        }}
        .controls label {{
            font-size: 0.875rem;
            color: #666;
            margin-right: 4px;
        }}
        .btn-group {{
            display: flex;
            gap: 4px;
        }}
        .btn {{
            padding: 6px 12px;
            font-size: 0.8rem;
            border: 1px solid #ddd;
            background: #fff;
            color: #333;
            border-radius: 4px;
            cursor: pointer;
            transition: all 0.15s ease;
        }}
        .btn:hover {{
            background: #f0f0f0;
            border-color: #ccc;
        }}
        .btn.active {{
            background: #0066cc;
            color: white;
            border-color: #0066cc;
        }}
        .viewer-container {{
            flex: 1;
            display: flex;
            justify-content: center;
            align-items: center;
            padding: 16px;
            min-height: 0;
        }}
        .viewer-wrapper {{
            width: 100%;
            height: 100%;
            max-width: 1400px;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
            overflow: hidden;
        }}
        /* Override 3Dmol viewer size to fill wrapper */
        .viewer-wrapper > div {{
            width: 100% !important;
            height: 100% !important;
        }}
    </style>
</head>
<body>
    <header>
        <h1>Protein Pocket Visualization: {pdb_id}</h1>
        <div class="controls">
            <label>Style:</label>
            <div class="btn-group">
                <button class="btn active" data-style="cartoon">Cartoon</button>
                <button class="btn" data-style="stick">Stick</button>
                <button class="btn" data-style="sphere">Sphere</button>
                <button class="btn" data-style="line">Line</button>
                <button class="btn" data-style="surface">Surface</button>
            </div>
        </div>
    </header>
    <div class="viewer-container">
        <div class="viewer-wrapper">
            {viewer_html}
        </div>
    </div>
    <script>
        var globalViewer = null;

        // Initialize viewer and set up controls
        if (typeof $3Dmolpromise !== 'undefined') {{
            $3Dmolpromise.then(function() {{
                setTimeout(function() {{
                    // Find the viewer
                    var viewers = document.querySelectorAll('[id^="3dmolviewer_"]');
                    viewers.forEach(function(el) {{
                        var viewerId = el.id.replace('3dmolviewer_', '');
                        globalViewer = window['viewer_' + viewerId];
                        if (globalViewer) {{
                            globalViewer.zoom(0.7);
                            // Set default representation to Cartoon
                            setStyle('cartoon');
                        }}
                    }});
                }}, 100);
            }});
        }}

        // Style change handler
        function setStyle(style) {{
            if (!globalViewer) return;

            // Clear all styles and surfaces first
            globalViewer.removeAllSurfaces();
            globalViewer.setStyle({{}}, {{}});

            // Apply new style to protein
            var styleSpec = {{}};
            switch(style) {{
                case 'cartoon':
                    styleSpec = {{cartoon: {{color: 'spectrum'}}}};
                    break;
                case 'stick':
                    styleSpec = {{stick: {{colorscheme: 'Jmol'}}}};
                    break;
                case 'sphere':
                    styleSpec = {{sphere: {{colorscheme: 'Jmol', scale: 0.3}}}};
                    break;
                case 'line':
                    styleSpec = {{line: {{colorscheme: 'Jmol'}}}};
                    break;
                case 'surface':
                    styleSpec = {{cartoon: {{color: 'spectrum'}}}};
                    globalViewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.7, color: 'white'}});
                    break;
            }}
            globalViewer.setStyle({{}}, styleSpec);
            globalViewer.render();
        }}

        // Set up button click handlers
        document.querySelectorAll('.btn[data-style]').forEach(function(btn) {{
            btn.addEventListener('click', function() {{
                // Update active state
                document.querySelectorAll('.btn[data-style]').forEach(function(b) {{
                    b.classList.remove('active');
                }});
                this.classList.add('active');

                // Apply style
                setStyle(this.dataset.style);
            }});
        }});
    </script>
</body>
</html>"""

    html_output = os.path.join(output_dir, "pocket_visualization.html")
    with open(html_output, "w") as f:
        f.write(html_content)
    print(f"Visualization saved to {html_output}")

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

    cmd.load(pdb_path, "structure")

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
    gif_path = os.path.join(output_dir, f"protein_pockets_rotation_{representation}.gif")
    imageio.mimsave(gif_path, images, duration=0.1, loop=0)
    print(f"\nRotating GIF saved to {gif_path}", flush=True)

    # Cleanup temporary frames
    for f in os.listdir(frame_dir):
        os.remove(os.path.join(frame_dir, f))
    os.rmdir(frame_dir)
