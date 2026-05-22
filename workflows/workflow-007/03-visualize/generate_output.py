# Part 3: Generate multi-protein visualization (HTML)
# github repo: https://github.com/cch1999/pocketeer
# doc: https://pocketeer.readthedocs.io/en/latest/

import json
import os

import numpy as np
import pocketeer as pt
from pocketeer.core.types import AlphaSphere

# =============================================================================
# CONFIGURATION
# =============================================================================

input_dir = "inputs"
output_dir = "outputs"

# Read PDB IDs from config.json
config_path = os.path.join(input_dir, "config.json")
with open(config_path, "r") as f:
    config = json.load(f)
pdb_ids = config["pdb_ids"]

# Visualization parameters
pocket_style = os.environ.get("PARAM_POCKET_STYLE", "filled_surfaces")
render_method = os.environ.get("PARAM_RENDER_METHOD", "draw")
representation = os.environ.get("PARAM_REPRESENTATION", "surface")
output_format = os.environ.get("PARAM_OUTPUT_FORMAT", "html")

# =============================================================================


def load_pockets_json(json_path: str, atomarray) -> list[pt.Pocket]:
    """Load pockets from a JSON file created by pt.write_pockets_json()."""
    with open(json_path, "r") as f:
        data = json.load(f)

    pockets = []
    for p in data:
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
        residues = [(r["chain_id"], r["res_id"], r["res_name"]) for r in p["residues"]]
        residue_set = set(residues)
        mask = np.array(
            [(c, r, n) in residue_set
             for c, r, n in zip(atomarray.chain_id, atomarray.res_id, atomarray.res_name)],
            dtype=bool,
        )
        pocket = pt.Pocket(
            pocket_id=p["pocket_id"],
            spheres=spheres,
            centroid=np.array(p["centroid"]),
            volume=p["volume"],
            score=p["score"],
            residues=residues,
            mask=mask,
        )
        pockets.append(pocket)
    return pockets


# =============================================================================

os.makedirs(output_dir, exist_ok=True)

print(f"Generating visualization for {len(pdb_ids)} protein(s)", flush=True)
print(f"  pocket_style={pocket_style}, representation={representation}", flush=True)
print(f"  render_method={render_method}, output_format={output_format}", flush=True)

# Load all protein data and generate viewer HTML fragments
protein_data = []  # list of { id, pockets_count, top_score, viewer_html | None }

receptor_cartoon = representation == "cartoon"
receptor_surface = representation == "surface"
sphere_scale = 2.0 if pocket_style == "single_sphere" else 1.0

for pdb_id in pdb_ids:
    pdb_path = os.path.join(input_dir, f"{pdb_id}.pdb")
    pockets_path = os.path.join(input_dir, f"{pdb_id}_pockets.json")

    print(f"  Loading {pdb_id}...", flush=True)
    atomarray = pt.load_structure(pdb_path)

    pockets = load_pockets_json(pockets_path, atomarray)
    pocket_count = len(pockets)
    top_score = f"{pockets[0].score:.2f}" if pockets else "—"

    viewer_html = None
    if pockets:
        viewer = pt.view_pockets(
            atomarray,
            pockets,
            receptor_cartoon=receptor_cartoon,
            receptor_surface=receptor_surface,
            sphere_scale=sphere_scale,
        )
        viewer_html = viewer._make_html()

    protein_data.append({
        "id": pdb_id,
        "pocket_count": pocket_count,
        "top_score": top_score,
        "viewer_html": viewer_html,
    })
    print(f"    {pocket_count} pockets (top score: {top_score})", flush=True)

# Build protein cards HTML
cards_html = ""
for p in protein_data:
    if p["viewer_html"]:
        viewer_block = f'<div class="viewer-content">{p["viewer_html"]}</div>'
    else:
        viewer_block = '<div class="viewer-content no-pockets"><p>No pockets detected</p><p class="hint">Try adjusting parameters</p></div>'

    cards_html += f"""
    <div class="protein-card" data-id="{p['id']}">
        <div class="card-header">
            <div class="card-info">
                <input type="checkbox" class="compare-check" data-id="{p['id']}" title="Select for compare">
                <h2>{p['id']}</h2>
                <span class="badge">{p['pocket_count']} pockets</span>
                <span class="score">top: {p['top_score']}</span>
            </div>
            <div class="card-actions">
                <button class="btn btn-focus" onclick="toggleFocus('{p['id']}')" title="Focus view">&#x2922;</button>
            </div>
        </div>
        {viewer_block}
    </div>"""

title = f"Pocket Analysis: {', '.join(pdb_ids)}"

html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <style>
        * {{ margin: 0; padding: 0; box-sizing: border-box; }}
        html, body {{ height: 100%; overflow: auto; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: #f0f2f5;
        }}
        header {{
            padding: 12px 20px;
            background: white;
            border-bottom: 1px solid #e0e0e0;
            box-shadow: 0 1px 3px rgba(0,0,0,0.1);
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 10px;
            position: sticky;
            top: 0;
            z-index: 100;
        }}
        h1 {{ color: #333; font-size: 1.15rem; font-weight: 500; }}
        .toolbar {{
            display: flex;
            gap: 8px;
            align-items: center;
        }}
        .btn {{
            padding: 6px 14px;
            font-size: 0.8rem;
            border: 1px solid #ddd;
            background: #fff;
            color: #333;
            border-radius: 4px;
            cursor: pointer;
            transition: all 0.15s;
        }}
        .btn:hover {{ background: #f0f0f0; border-color: #ccc; }}
        .btn.active {{ background: #0066cc; color: white; border-color: #0066cc; }}
        .btn-focus {{ font-size: 1rem; padding: 4px 8px; line-height: 1; }}

        /* Grid layout */
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(480px, 1fr));
            gap: 16px;
            padding: 16px;
        }}
        .grid.compare-mode {{
            grid-template-columns: repeat(auto-fill, minmax(400px, 1fr));
        }}

        /* Cards */
        .protein-card {{
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.08);
            overflow: hidden;
            display: flex;
            flex-direction: column;
        }}
        .protein-card.hidden {{ display: none; }}
        .protein-card.focused {{
            position: fixed;
            top: 0; left: 0; right: 0; bottom: 0;
            z-index: 200;
            border-radius: 0;
            margin: 0;
        }}
        .card-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding: 10px 14px;
            border-bottom: 1px solid #eee;
            background: #fafafa;
        }}
        .card-info {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        .card-info h2 {{ font-size: 1rem; font-weight: 600; color: #333; }}
        .badge {{
            font-size: 0.75rem;
            padding: 2px 8px;
            border-radius: 10px;
            background: #e8f4fd;
            color: #0066cc;
        }}
        .score {{ font-size: 0.75rem; color: #888; }}
        .compare-check {{ width: 16px; height: 16px; cursor: pointer; }}

        /* Viewer */
        .viewer-content {{
            flex: 1;
            min-height: 400px;
            position: relative;
        }}
        .viewer-content > div {{
            width: 100% !important;
            height: 100% !important;
            min-height: 400px;
        }}
        .focused .viewer-content {{
            min-height: calc(100vh - 50px);
        }}
        .no-pockets {{
            display: flex;
            flex-direction: column;
            justify-content: center;
            align-items: center;
            color: #999;
            font-size: 1.1rem;
            min-height: 400px;
            background: #fafafa;
        }}
        .no-pockets .hint {{ font-size: 0.85rem; margin-top: 8px; color: #bbb; }}

        /* Compare mode indicator */
        .compare-bar {{
            display: none;
            padding: 8px 20px;
            background: #0066cc;
            color: white;
            font-size: 0.85rem;
            align-items: center;
            gap: 12px;
        }}
        .compare-bar.visible {{ display: flex; }}
    </style>
</head>
<body>
    <header>
        <h1>{title}</h1>
        <div class="toolbar">
            <label style="font-size:0.85rem;color:#666;">Style:</label>
            <div class="btn-group" style="display:flex;gap:4px;">
                <button class="btn style-btn active" data-style="cartoon">Cartoon</button>
                <button class="btn style-btn" data-style="stick">Stick</button>
                <button class="btn style-btn" data-style="sphere">Sphere</button>
                <button class="btn style-btn" data-style="line">Line</button>
                <button class="btn style-btn" data-style="surface">Surface</button>
            </div>
            <span style="color:#ddd;">|</span>
            <button class="btn" onclick="showAll()">Grid View</button>
            <button class="btn" id="compareBtn" onclick="enterCompare()">Compare Selected</button>
        </div>
    </header>
    <div class="compare-bar" id="compareBar">
        <span id="compareCount">0 selected</span>
        <button class="btn" style="background:rgba(255,255,255,0.2);color:white;border-color:rgba(255,255,255,0.3);" onclick="exitCompare()">Exit Compare</button>
    </div>
    <div class="grid" id="grid">
        {cards_html}
    </div>
    <script>
        var focusedId = null;
        var compareMode = false;
        var allViewers = [];

        // Collect all 3Dmol viewers once loaded
        function initViewers() {{
            allViewers = [];
            document.querySelectorAll('[id^="3dmolviewer_"]').forEach(function(el) {{
                var vid = el.id.replace('3dmolviewer_', '');
                var v = window['viewer_' + vid];
                if (v) {{
                    allViewers.push(v);
                    v.zoom(0.7);
                }}
            }});
            // Set default style to cartoon
            setStyle('cartoon');
        }}

        // Apply style to all viewers
        function setStyle(style) {{
            allViewers.forEach(function(viewer) {{
                viewer.removeAllSurfaces();
                viewer.setStyle({{}}, {{}});
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
                        viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.7, color: 'white'}});
                        break;
                }}
                viewer.setStyle({{}}, styleSpec);
                viewer.render();
            }});
        }}

        // Style button handlers
        document.querySelectorAll('.style-btn').forEach(function(btn) {{
            btn.addEventListener('click', function() {{
                document.querySelectorAll('.style-btn').forEach(function(b) {{ b.classList.remove('active'); }});
                this.classList.add('active');
                setStyle(this.dataset.style);
            }});
        }});

        // Initialize viewers when 3Dmol is ready
        if (typeof $3Dmolpromise !== 'undefined') {{
            $3Dmolpromise.then(function() {{ setTimeout(initViewers, 200); }});
        }}

        function toggleFocus(id) {{
            var card = document.querySelector('.protein-card[data-id="' + id + '"]');
            if (focusedId === id) {{
                card.classList.remove('focused');
                focusedId = null;
                document.body.style.overflow = 'auto';
            }} else {{
                if (focusedId) {{
                    document.querySelector('.protein-card[data-id="' + focusedId + '"]').classList.remove('focused');
                }}
                card.classList.remove('hidden');
                card.classList.add('focused');
                focusedId = id;
                document.body.style.overflow = 'hidden';
            }}
            resizeViewers();
        }}

        function showAll() {{
            exitCompare();
            if (focusedId) {{
                document.querySelector('.protein-card[data-id="' + focusedId + '"]').classList.remove('focused');
                focusedId = null;
                document.body.style.overflow = 'auto';
            }}
            document.querySelectorAll('.protein-card').forEach(function(c) {{ c.classList.remove('hidden'); }});
            document.getElementById('grid').classList.remove('compare-mode');
            resizeViewers();
        }}

        function enterCompare() {{
            var checked = document.querySelectorAll('.compare-check:checked');
            if (checked.length < 2) {{ alert('Select at least 2 proteins to compare'); return; }}

            compareMode = true;
            var selectedIds = Array.from(checked).map(function(cb) {{ return cb.dataset.id; }});

            document.querySelectorAll('.protein-card').forEach(function(card) {{
                if (selectedIds.indexOf(card.dataset.id) >= 0) {{
                    card.classList.remove('hidden');
                }} else {{
                    card.classList.add('hidden');
                }}
            }});
            document.getElementById('grid').classList.add('compare-mode');
            document.getElementById('compareBar').classList.add('visible');
            document.getElementById('compareCount').textContent = checked.length + ' selected';
            resizeViewers();
        }}

        function exitCompare() {{
            compareMode = false;
            document.querySelectorAll('.protein-card').forEach(function(c) {{ c.classList.remove('hidden'); }});
            document.querySelectorAll('.compare-check').forEach(function(cb) {{ cb.checked = false; }});
            document.getElementById('grid').classList.remove('compare-mode');
            document.getElementById('compareBar').classList.remove('visible');
            resizeViewers();
        }}

        function resizeViewers() {{
            setTimeout(function() {{
                if (typeof $3Dmol !== 'undefined') {{
                    document.querySelectorAll('[id^="3dmolviewer_"]').forEach(function(el) {{
                        var vid = el.id.replace('3dmolviewer_', '');
                        var v = window['viewer_' + vid];
                        if (v) v.resize();
                    }});
                }}
            }}, 100);
        }}

        // ESC to exit focus
        document.addEventListener('keydown', function(e) {{
            if (e.key === 'Escape') {{
                if (focusedId) toggleFocus(focusedId);
                else if (compareMode) exitCompare();
            }}
        }});

        // Update compare count on checkbox change
        document.querySelectorAll('.compare-check').forEach(function(cb) {{
            cb.addEventListener('change', function() {{
                var count = document.querySelectorAll('.compare-check:checked').length;
                document.getElementById('compareBtn').textContent = count > 0 ? 'Compare (' + count + ')' : 'Compare Selected';
            }});
        }});
    </script>
</body>
</html>"""

if output_format in ("html", "both"):
    html_output = os.path.join(output_dir, "pocket_visualization.html")
    with open(html_output, "w") as f:
        f.write(html_content)
    print(f"Visualization saved to {html_output}", flush=True)

# Generate rotating GIF using PyMOL if requested
if output_format in ("gif", "both"):
    import multiprocessing
    import shutil

    import imageio.v2 as imageio
    from pymol import cmd

    print("\nInitializing PyMOL...", flush=True)
    num_cpus = multiprocessing.cpu_count()
    cmd.set("max_threads", num_cpus)
    print(f"Using {num_cpus} CPU threads for ray tracing", flush=True)

    for p in protein_data:
        pdb_id = p["id"]
        gif_atomarray = pt.load_structure(os.path.join(input_dir, f"{pdb_id}.pdb"))
        pockets = load_pockets_json(os.path.join(input_dir, f"{pdb_id}_pockets.json"), gif_atomarray)

        if not pockets:
            print(f"\nSkipping GIF for {pdb_id} (no pockets)", flush=True)
            continue

        print(f"\n{'='*50}", flush=True)
        print(f"Generating GIF for {pdb_id}", flush=True)
        print(f"{'='*50}", flush=True)

        cmd.reinitialize()
        pdb_path = os.path.join(input_dir, f"{pdb_id}.pdb")
        cmd.load(pdb_path, "structure")

        # Style the protein
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

        # Highlight top pockets
        colors = ["red", "orange", "yellow", "green", "magenta"]
        if pocket_style == "filled_surfaces":
            for i, pocket in enumerate(pockets[:5]):
                pocket_name = f"pocket_{i}"
                for j, sphere in enumerate(pocket.spheres):
                    atom_name = f"{pocket_name}_s{j}"
                    cmd.pseudoatom(atom_name, pos=list(sphere.center), vdw=sphere.radius)
                cmd.group(pocket_name, f"{pocket_name}_s*")
                cmd.show("spheres", pocket_name)
                cmd.color(colors[i], pocket_name)
                cmd.set("sphere_transparency", 0.3, pocket_name)
        elif pocket_style == "single_sphere":
            for i, pocket in enumerate(pockets[:5]):
                pocket_name = f"pocket_{i}"
                cmd.pseudoatom(pocket_name, pos=list(pocket.centroid))
                cmd.show("spheres", pocket_name)
                cmd.color(colors[i], pocket_name)
                cmd.set("sphere_scale", 2.0, pocket_name)

        cmd.zoom("all", buffer=5)

        # Generate frames
        step = 10
        total_frames = 360 // step
        images = []
        frame_dir = "/tmp/frames"
        os.makedirs(frame_dir, exist_ok=True)

        print(f"Rendering {total_frames} frames ({render_method})...", flush=True)
        for i, angle in enumerate(range(0, 360, step)):
            cmd.rotate("y", float(step))
            frame_path = f"{frame_dir}/frame_{angle:03d}.png"
            if render_method == "ray":
                cmd.ray(512, 512)
                cmd.png(frame_path)
            else:
                cmd.draw(512, 512)
                cmd.png(frame_path, ray=0)
            images.append(imageio.imread(frame_path))

        # Save GIF
        gif_path = os.path.join(output_dir, f"{pdb_id}_pockets_{representation}.gif")
        imageio.mimsave(gif_path, images, duration=0.1, loop=0)
        print(f"GIF saved to {gif_path}", flush=True)

        # Cleanup frames
        shutil.rmtree(frame_dir)
