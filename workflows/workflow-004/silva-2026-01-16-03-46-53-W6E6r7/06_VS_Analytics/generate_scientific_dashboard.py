#!/usr/bin/env python3
"""
Node 6: Virtual Screening Dashboard
Creates comprehensive interactive dashboard for docking results
"""

import json
from pathlib import Path

def generate_report():
    print("[VS Analytics] Generating Unified Scientific Dashboard...")
    
    # Load Data
    json_path = Path("outputs/virtual_screening_metadata.json")
    if not json_path.exists():
        json_path = Path("virtual_screening_metadata.json")
    if not json_path.exists():
        print("[VS Analytics] ERROR: virtual_screening_metadata.json not found")
        return
    
    with open(json_path) as f:
        data = json.load(f)
    
    # Load Receptor
    receptor_file = data.get("receptor_file", "receptor.pdbqt")
    receptor_path = Path("outputs") / receptor_file
    if not receptor_path.exists():
        receptor_path = Path(receptor_file)
    
    receptor_content = ""
    if receptor_path.exists():
        with open(receptor_path) as f:
            receptor_content = json.dumps(f.read())
    else:
        receptor_content = '""'

    # Load Ligands and Poses
    ligands_db = {}  # Pose content
    ligands_meta = {}  # Ligand metadata
    cards_html = ""
    
    # Support both "docking_results" and "ligand_results" keys
    ligand_results = data.get("ligand_results", data.get("docking_results", []))
    
    for res in ligand_results:
        if res.get("status") == "failed":
            continue
            
        lig_name = res.get("ligand")
        best_aff = res.get("best_affinity", 0)
        rank = res.get("rank", 0)
        poses = res.get("poses", [])
        
        ligands_meta[lig_name] = {
            "rank": rank,
            "best_affinity": best_aff,
            "poses": poses
        }
        
        for p in poses:
            p_num = p["pose"]
            p_key = f"{lig_name}_pose{p_num}"
            p_file = Path("outputs") / p["file"]
            if not p_file.exists():
                p_file = Path(p["file"])
            if p_file.exists():
                with open(p_file) as f:
                    ligands_db[p_key] = f.read()
        
        # Build Sidebar Cards
        affinity_color = "#22c55e" if best_aff < -8 else "#fbbf24" if best_aff < -6 else "#ef4444"
        cards_html += f"""
        <div class="ligand-card" onclick="app.loadLigand('{lig_name}')" id="card-{lig_name}">
            <div class="card-header">
                <div class="rank-badge">#{rank}</div>
                <strong class="lig-name">{lig_name}</strong>
            </div>
            <div class="affinity-display" style="background: {affinity_color};">
                {best_aff:.2f} kcal/mol
            </div>
            <div class="card-meta">{len(poses)} poses</div>
        </div>
        """

    # Summary Stats
    summary = data.get("summary", {})
    stats_html = f"""
    <div class="stat-grid">
        <div class="stat-item">
            <div class="stat-val">{summary.get('total_compounds', 0)}</div>
            <div class="stat-label">Total Ligands</div>
        </div>
        <div class="stat-item">
            <div class="stat-val">{summary.get('successful_dockings', 0)}</div>
            <div class="stat-label">Successes</div>
        </div>
        <div class="stat-item">
            <div class="stat-val">{summary.get('best_affinity', 0):.2f}</div>
            <div class="stat-label">Best Energy</div>
        </div>
    </div>
    """

    # Parameters
    params = data.get("parameters", {})
    grid = data.get("grid_config", {})
    protocol_html = f"""
    <div class="prop-row"><span class="label">Software</span><span class="value">{data.get('software', 'AutoDock Vina')}</span></div>
    <div class="prop-row"><span class="label">Exhaustiveness</span><span class="value">{params.get('exhaustiveness', 'N/A')}</span></div>
    <div class="prop-row"><span class="label">Energy Range</span><span class="value">{params.get('energy_range', 'N/A')} kcal/mol</span></div>
    <div class="prop-row"><span class="label">Grid Size</span><span class="value">{grid.get('size_x', 0)} × {grid.get('size_y', 0)} × {grid.get('size_z', 0)} Å</span></div>
    """

    # JSON Packing
    js_ligands_db = json.dumps(ligands_db)
    js_ligands_meta = json.dumps(ligands_meta)
    js_grid = json.dumps(grid)

    # HTML Template
    html_template = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Virtual Screening Dashboard</title>
    <script src="https://cdn.jsdelivr.net/npm/ngl@2.0.0-dev.37/dist/ngl.js"></script>
    <style>
        :root {
            --bg: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            --panel: rgba(255,255,255,0.95);
            --text: #212529;
            --text-dim: #6c757d;
            --accent: #667eea;
            --border: #e9ecef;
            --success: #22c55e;
            --warning: #fbbf24;
            --error: #ef4444;
        }
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { 
            margin: 0;
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            display: grid;
            grid-template-columns: 320px 1fr 360px;
            height: 100vh;
            background: var(--bg);
            color: var(--text);
        }
        
        #sidebar { 
            background: var(--panel);
            border-right: 3px solid var(--accent);
            display: flex;
            flex-direction: column;
            overflow: hidden;
            box-shadow: 5px 0 15px rgba(0,0,0,0.1);
        }
        #viewer-container { 
            position: relative;
            background: #1a1a2e;
            overflow: hidden;
        }
        #inspector { 
            background: var(--panel);
            border-left: 3px solid var(--accent);
            display: flex;
            flex-direction: column;
            overflow: hidden;
            box-shadow: -5px 0 15px rgba(0,0,0,0.1);
        }
        
        .hdr { 
            padding: 25px;
            border-bottom: 2px solid var(--border);
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
        }
        .hdr h1 { 
            font-size: 20px;
            margin: 0 0 15px 0;
        }
        .stat-grid { 
            display: grid;
            grid-template-columns: 1fr 1fr 1fr;
            gap: 10px;
        }
        .stat-item { 
            background: rgba(255,255,255,0.2);
            padding: 12px;
            border-radius: 8px;
            text-align: center;
        }
        .stat-val { 
            font-size: 20px;
            font-weight: bold;
            color: white;
        }
        .stat-label { 
            font-size: 10px;
            color: rgba(255,255,255,0.9);
            text-transform: uppercase;
            margin-top: 5px;
            letter-spacing: 0.5px;
        }
        
        #hit-list { 
            overflow-y: auto;
            padding: 15px;
            height: 40%;
            border-bottom: 2px solid var(--border);
        }
        .ligand-card { 
            margin-bottom: 12px;
            border: 2px solid var(--border);
            border-radius: 12px;
            cursor: pointer;
            transition: all 0.3s;
            overflow: hidden;
            background: white;
        }
        .ligand-card:hover { 
            border-color: var(--accent);
            transform: translateX(5px);
            box-shadow: 0 5px 15px rgba(102, 126, 234, 0.3);
        }
        .ligand-card.active { 
            border-color: var(--accent);
            background: linear-gradient(to right, rgba(102, 126, 234, 0.1), white);
            box-shadow: 0 5px 15px rgba(102, 126, 234, 0.4);
        }
        .card-header { 
            display: flex;
            align-items: center;
            gap: 12px;
            padding: 15px;
        }
        .rank-badge { 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            width: 35px;
            height: 35px;
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            font-weight: bold;
            font-size: 14px;
        }
        .lig-name {
            font-size: 15px;
            color: var(--text);
        }
        .affinity-display { 
            padding: 10px 15px;
            color: white;
            font-weight: bold;
            font-size: 16px;
            text-align: center;
        }
        .card-meta { 
            font-size: 11px;
            color: var(--text-dim);
            padding: 8px 15px;
            background: #f8f9fa;
            border-top: 1px solid var(--border);
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        
        #viewport { 
            width: 100%;
            height: 100%;
        }
        #hud { 
            position: absolute;
            top: 20px;
            left: 20px;
            background: rgba(30,41,59,0.9);
            padding: 20px;
            border-radius: 12px;
            border: 2px solid var(--accent);
            backdrop-filter: blur(10px);
            pointer-events: none;
            min-width: 250px;
        }
        #hud h2 { 
            font-size: 16px;
            margin: 0 0 12px 0;
            color: white;
        }
        .hud-row { 
            font-size: 13px;
            margin-bottom: 8px;
            color: rgba(255,255,255,0.9);
        }
        .hud-row strong {
            color: #38bdf8;
        }
        
        .controls { 
            position: absolute;
            bottom: 20px;
            left: 50%;
            transform: translateX(-50%);
            background: rgba(255,255,255,0.95);
            padding: 15px 20px;
            border-radius: 50px;
            display: flex;
            gap: 12px;
            border: 2px solid var(--accent);
            box-shadow: 0 5px 20px rgba(0,0,0,0.3);
        }
        .btn { 
            background: white;
            border: 2px solid var(--border);
            color: var(--text);
            padding: 8px 18px;
            border-radius: 25px;
            cursor: pointer;
            font-size: 13px;
            font-weight: 600;
            transition: all 0.3s;
        }
        .btn:hover { 
            border-color: var(--accent);
            color: var(--accent);
            transform: translateY(-2px);
        }
        .btn.active { 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-color: var(--accent);
        }
        
        .panel { 
            padding: 20px;
            border-bottom: 2px solid var(--border);
        }
        .p-title { 
            font-size: 12px;
            font-weight: bold;
            color: var(--text-dim);
            text-transform: uppercase;
            letter-spacing: 1px;
            margin-bottom: 15px;
        }
        .prop-row { 
            display: flex;
            justify-content: space-between;
            font-size: 13px;
            margin-bottom: 10px;
            padding-bottom: 10px;
            border-bottom: 1px solid var(--border);
        }
        .prop-row:last-child {
            border-bottom: none;
        }
        .label { 
            color: var(--text-dim);
        }
            font-weight: 600;
            font-family: 'Courier New', monospace;
        }
        
        .pose-row { 
            display: flex;
            gap: 8px;
            flex-wrap: wrap;
        }
        .pose-chip { 
            padding: 8px 14px;
            border: 2px solid var(--border);
            border-radius: 8px;
            background: white;
            font-size: 12px;
            cursor: pointer;
            font-weight: 600;
            transition: all 0.3s;
        }
        .pose-chip:hover { 
            border-color: var(--accent);
            background: rgba(102, 126, 234, 0.1);
        }
        .pose-chip.active { 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-color: var(--accent);
        }
        
        .int-table { 
            width: 100%;
            border-collapse: collapse;
            font-size: 12px;
        }
        .int-table th { 
            text-align: left;
            padding-bottom: 10px;
            border-bottom: 2px solid var(--border);
            color: var(--text-dim);
            font-weight: 700;
            text-transform: uppercase;
            font-size: 11px;
        }
        .int-table td { 
            padding: 10px 5px;
            border-bottom: 1px solid var(--border);
        }
        .res-tag { 
            font-family: 'Courier New', monospace;
            background: rgba(102, 126, 234, 0.1);
            padding: 3px 8px;
            border-radius: 4px;
            font-weight: 600;
        }
        .dist-tag { 
            font-weight: bold;
            color: var(--accent);
        }
    </style>
</head>
<body>

    <div id="sidebar">
        <div class="hdr">
            <h1>🎯 Virtual Screening</h1>
            __SUMMARY_STATS_HTML__
        </div>
        <div id="hit-list">
            __CARDS_HTML__
        </div>
        <div class="panel" style="flex:1; display:flex; flex-direction:column; overflow:hidden; background: var(--panel);">
            <div class="p-title" style="padding: 15px 15px 5px 15px;">Interaction Analysis (<span id="int-count">0</span>)</div>
            <div style="overflow-y:auto; flex:1; padding: 0 15px 15px 15px;">
                <table class="int-table">
                    <thead><tr><th>Residue</th><th>Dist (Å)</th></tr></thead>
                    <tbody id="int-body"></tbody>
                </table>
            </div>
        </div>
    </div>

    <div id="viewer-container">
        <div id="viewport"></div>
        <div id="hud">
            <h2>📊 Active Result</h2>
            <div class="hud-row"><strong>Ligand:</strong> <span id="hud-name">-</span></div>
            <div class="hud-row"><strong>Affinity:</strong> <span id="hud-aff">-</span></div>
            <div class="hud-row"><strong>Pose:</strong> <span id="hud-pose">-</span></div>
        </div>
        <div class="controls">
            <button class="btn active" id="btn-surface" onclick="app.toggleSurface()">Surface</button>
            <button class="btn" id="btn-backbone" onclick="app.toggleBackbone()">Backbone</button>
            <button class="btn active" id="btn-grid" onclick="app.toggleBox()">Grid Box</button>
            <div style="width: 1px; height: 20px; background: #dee2e6; margin: 0 5px;"></div>
            <button class="btn" onclick="app.toggleStyle()">Style: Licorice</button>
            <button class="btn" onclick="app.resetView()">Reset View</button>
        </div>
    </div>

    <div id="inspector">
        <div class="panel" style="background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white;">
            <div class="p-title" style="color: rgba(255,255,255,0.9);">Ligand Context</div>
            <div class="prop-row" style="border-color: rgba(255,255,255,0.2);">
                <span class="label" style="color: rgba(255,255,255,0.9);">Rank</span>
                <span id="p-rank" class="value" style="color: white;">-</span>
            </div>
            <div class="prop-row" style="border-color: rgba(255,255,255,0.2);">
                <span class="label" style="color: rgba(255,255,255,0.9);">Best Energy</span>
                <span id="p-aff" class="value" style="color: white;">-</span>
            </div>
        </div>
        
        <div class="panel">
            <div class="p-title">Switch Poses</div>
            <div id="pose-container" class="pose-row"></div>
        </div>
        
        <div class="panel">
            <div class="p-title">Protocol Parameters</div>
            __PROTOCOL_HTML__
        </div>
    </div>

    <script>
        const RECEPTOR_RAW = __RECEPTOR_JSON__;
        const LIGANDS_DB = __LIGANDS_DB__;
        const LIGANDS_META = __LIGANDS_META__;
        const GRID = __GRID_JSON__;

        class UnifiedApp {
            constructor() {
                this.stage = new NGL.Stage("viewport", { backgroundColor: "#1a1a2e" });
                this.comps = { rec: null, lig: null, box: null, back: null };
                this.current = { name: null, pose: 1 };
                this.showSurface = true;
                this.showBackbone = false;
                this.ligandStyle = 'licorice';
                window.addEventListener("resize", () => this.stage.handleResize());
            }

            async init() {
                if (!RECEPTOR_RAW) return;
                const blob = new Blob([RECEPTOR_RAW], { type: 'text/plain' });
                this.comps.rec = await this.stage.loadFile(blob, { ext: "pdbqt" });
                this.comps.rec.addRepresentation("cartoon", { color: "#64748b" });
                this.surfaceRep = this.comps.rec.addRepresentation("surface", {
                    color: "white", opacity: 0.15, side: "front", visible: this.showSurface
                });
                this.backboneRep = this.comps.rec.addRepresentation("backbone", {
                    color: "lightgray", radiusScale: 2.0, visible: this.showBackbone
                });
                if (GRID && GRID.center_x !== undefined) this.initBox();
                this.resetView();
            }

            initBox() {
                const shape = new NGL.Shape("box");
                const cx = GRID.center_x, cy = GRID.center_y, cz = GRID.center_z;
                const sx = GRID.size_x, sy = GRID.size_y, sz = GRID.size_z;
                const c = [1, 0.75, 0];
                const x = sx/2, y = sy/2, z = sz/2;
                const p = [[cx-x, cy-y, cz-z], [cx+x, cy-y, cz-z], [cx+x, cy+y, cz-z], [cx-x, cy+y, cz-z],
                           [cx-x, cy-y, cz+z], [cx+x, cy-y, cz+z], [cx+x, cy+y, cz+z], [cx-x, cy+y, cz+z]];
                const edges = [[0,1],[1,2],[2,3],[3,0], [4,5],[5,6],[6,7],[7,4], [0,4],[1,5],[2,6],[3,7]];
                edges.forEach(e => shape.addLine(p[e[0]], p[e[1]], c));
                this.comps.box = this.stage.addComponentFromObject(shape);
                this.comps.box.addRepresentation("buffer", {linewidth: 4});
                this.comps.box.setVisibility(true);
            }

            async loadLigand(name) {
                this.current.name = name;
                document.querySelectorAll('.ligand-card').forEach(c => c.classList.remove('active'));
                document.getElementById('card-' + name).classList.add('active');
                const meta = LIGANDS_META[name];
                document.getElementById('hud-name').innerText = name;
                document.getElementById('hud-aff').innerText = meta.best_affinity + " kcal/mol";
                document.getElementById('p-rank').innerText = "#" + meta.rank;
                document.getElementById('p-aff').innerText = meta.best_affinity + " kcal/mol";
                const cont = document.getElementById('pose-container');
                cont.innerHTML = "";
                meta.poses.forEach(p => {
                    const chip = document.createElement('div');
                    chip.className = 'pose-chip';
                    chip.innerText = "Pose " + p.pose + " (" + p.affinity + ")";
                    chip.id = `pose-${p.pose}`;
                    chip.onclick = () => this.loadPose(p.pose);
                    cont.appendChild(chip);
                });
                await this.loadPose(1);
            }

            async loadPose(num) {
                this.current.pose = num;
                document.querySelectorAll('.pose-chip').forEach(c => c.classList.remove('active'));
                document.getElementById(`pose-${num}`).classList.add('active');
                document.getElementById('hud-pose').innerText = num;
                if (this.comps.lig) this.stage.removeComponent(this.comps.lig);
                const key = `${this.current.name}_pose${num}`;
                const blob = new Blob([LIGANDS_DB[key]], { type: 'text/plain' });
                const params = this.ligandStyle === 'licorice' ? 
                    { colorScheme: "element", multipleBond: true } : 
                    { colorScheme: "element", multipleBond: true }; // ball+stick handled differently by type
                
                this.comps.lig = await this.stage.loadFile(blob, { ext: "pdbqt" });
                this.comps.lig.addRepresentation(this.ligandStyle, params);
                this.analyze();
            }

            analyze() {
                const interactions = [];
                const rec = this.comps.rec.structure;
                const lig = this.comps.lig.structure;
                lig.eachAtom(la => {
                    if (la.element === "H") return;
                    rec.eachResidue(rp => {
                        let minD = 999;
                        rp.eachAtom(ra => {
                            if (ra.element === "H") return;
                            const d = Math.sqrt((ra.x-la.x)**2 + (ra.y-la.y)**2 + (ra.z-la.z)**2);
                            if (d < minD) minD = d;
                        });
                        if (minD < 4.0) {
                            if (!interactions.find(i => i.resNo === rp.resno && i.chain === rp.chainname)) {
                                interactions.push({ 
                                    res: rp.resname, 
                                    resNo: rp.resno,
                                    chain: rp.chainname,
                                    dist: minD 
                                });
                            }
                        }
                    });
                });
                interactions.sort((a,b) => a.dist - b.dist);
                const body = document.getElementById('int-body');
                body.innerHTML = "";
                document.getElementById('int-count').innerText = interactions.length;
                interactions.slice(0, 20).forEach(i => {
                    const tr = document.createElement('tr');
                    tr.innerHTML = `<td><span class="res-tag">${i.res}${i.resNo}:${i.chain}</span></td><td class="dist-tag">${i.dist.toFixed(2)}</td>`;
                    body.appendChild(tr);
                });
                this.comps.lig.autoView();
            }

            toggleSurface() {
                this.showSurface = !this.showSurface;
                this.surfaceRep.setVisibility(this.showSurface);
                document.getElementById('btn-surface').classList.toggle('active', this.showSurface);
            }
            
            toggleBackbone() {
                this.showBackbone = !this.showBackbone;
                this.backboneRep.setVisibility(this.showBackbone);
                document.getElementById('btn-backbone').classList.toggle('active', this.showBackbone);
            }
            
            toggleStyle() {
                this.ligandStyle = this.ligandStyle === 'licorice' ? 'ball+stick' : 'licorice';
                event.target.innerText = "Style: " + (this.ligandStyle === 'licorice' ? 'Licorice' : 'Ball+Stick');
                if (this.current.name) this.loadPose(this.current.pose);
            }

            toggleBox() {
                if(this.comps.box) {
                    this.comps.box.setVisibility(!this.comps.box.visible);
                    document.getElementById('btn-grid').classList.toggle('active', this.comps.box.visible);
                }
            }

            resetView() {
                if (GRID && GRID.center_x !== undefined) {
                    this.stage.animationControls.zoomMove([GRID.center_x, GRID.center_y, GRID.center_z], -30, 1000);
                } else { 
                    this.stage.autoView();
                }
            }
        }

        const app = new UnifiedApp();
        document.addEventListener("DOMContentLoaded", () => {
            app.init();
            const keys = Object.keys(LIGANDS_META);
            if(keys.length > 0) app.loadLigand(keys[0]);
        });
    </script>
</body>
</html>
"""

    final_html = html_template.replace("__SUMMARY_STATS_HTML__", stats_html)
    final_html = final_html.replace("__CARDS_HTML__", cards_html)
    final_html = final_html.replace("__PROTOCOL_HTML__", protocol_html)
    final_html = final_html.replace("__RECEPTOR_JSON__", receptor_content)
    final_html = final_html.replace("__LIGANDS_DB__", js_ligands_db)
    final_html = final_html.replace("__LIGANDS_META__", js_ligands_meta)
    final_html = final_html.replace("__GRID_JSON__", js_grid)
    
    # Save to outputs
    output_path = Path("outputs/report.html")
    output_path.parent.mkdir(exist_ok=True)
    with open(output_path, "w") as f:
        f.write(final_html)
        
    print(f"[VS Analytics] ✅ Dashboard generated: {output_path}")

if __name__ == "__main__":
    generate_report()