#!/usr/bin/env python3
"""
Node 5: Pocket Discovery Report
Creates interactive visualization of predicted binding pockets
"""

import json
from pathlib import Path

def generate_report():
    print("[Pocket Discovery] Generating HTML report from metadata...")
    
    # Read data from JSON
    json_file = Path("outputs/pocket_discovery_metadata.json")
    if not json_file.exists():
        json_file = Path("pocket_discovery_metadata.json")
    if not json_file.exists():
        print("[Pocket Discovery] ERROR: pocket_discovery_metadata.json not found")
        return
    
    with open(json_file) as f:
        data = json.load(f)
    
    # Extract Grid Selection
    grid_sel = data.get("grid_selection", {})
    grid_params = grid_sel.get("grid_params", {})
    grid_js = json.dumps(grid_params)
    
    # Read protein PDB for visualization
    protein_pdb_file = data.get("protein_pdb", "protein.pdb")
    protein_pdb_path = Path("outputs") / protein_pdb_file
    if not protein_pdb_path.exists():
        protein_pdb_path = Path(protein_pdb_file)
    if protein_pdb_path.exists():
        with open(protein_pdb_path) as f:
            protein_pdb_content = f.read()
        protein_pdb_js = json.dumps(protein_pdb_content)
    else:
        print(f"[Pocket Discovery] WARNING: {protein_pdb_file} not found")
        protein_pdb_js = '""'
    
    # Build pockets table and coordinates for JS
    pockets_rows = ""
    pockets_js_list = []
    
    for pocket in data.get("pockets", []):
        rank = pocket.get("rank", 0)
        rank_class = "top-pocket" if rank == 1 else ""
        pockets_rows += f"""
        <tr class="{rank_class}" onclick="focusPocket({rank-1})">
            <td><strong>{pocket['pocket_name']}</strong></td>
            <td>#{rank}</td>
            <td>{pocket['score']:.3f}</td>
            <td>{pocket.get('probability', 0):.1%}</td>
            <td>{pocket.get('residue_count', 0)}</td>
            <td>{pocket.get('surface_atoms', 0)}</td>
        </tr>"""
        
        # Add to JS list for sphere drawing
        if pocket.get("center_x") is not None:
            pockets_js_list.append({
                "name": pocket["pocket_name"],
                "rank": rank,
                "x": pocket["center_x"],
                "y": pocket["center_y"],
                "z": pocket["center_z"],
                "score": pocket["score"]
            })
    
    pockets_json_str = json.dumps(pockets_js_list)
    top_pocket = data.get("top_pocket")
    if not top_pocket:
        top_pocket = {"pocket_name": "None", "score": 0, "probability": 0, "residue_count": 0, "surface_atoms": 0, "rank": 0}
    
    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Pocket Discovery Report</title>
    <script src="https://cdn.jsdelivr.net/npm/ngl@2.0.0-dev.37/dist/ngl.js"></script>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body { 
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }
        .container { max-width: 1400px; margin: 0 auto; }
        .header { 
            background: rgba(255,255,255,0.95);
            border-radius: 12px;
            padding: 30px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            margin-bottom: 20px;
        }
        .header-content { 
            display: flex; 
            align-items: center; 
            justify-content: space-between; 
        }
        .title-section { 
            display: flex; 
            align-items: center; 
            gap: 20px; 
        }
        .page-title { 
            font-size: 32px;
            font-weight: 700;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        .node-label { 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 8px 16px;
            border-radius: 20px;
            font-size: 14px;
            font-weight: 600;
        }
        .legend { 
            font-size: 12px; 
            display: flex; 
            gap: 20px; 
        }
        .legend-item { 
            display: flex; 
            align-items: center; 
            gap: 8px;
            background: white;
            padding: 8px 12px;
            border-radius: 8px;
        }
        .dot { 
            width: 12px;
            height: 12px;
            border-radius: 50%;
        }
        .main-content { 
            display: grid;
            grid-template-columns: 450px 1fr;
            gap: 20px;
        }
        .sidebar { 
            background: rgba(255,255,255,0.95);
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            overflow: hidden;
            height: 800px;
            display: flex;
            flex-direction: column;
        }
        .section { 
            padding: 25px;
            border-bottom: 1px solid #e9ecef;
        }
        .section-title { 
            font-size: 12px;
            font-weight: 700;
            text-transform: uppercase;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin-bottom: 20px;
            letter-spacing: 1px;
        }
        .top-pocket-card { 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            border-radius: 10px;
            box-shadow: 0 5px 15px rgba(102, 126, 234, 0.3);
        }
        .pocket-name { 
            font-size: 20px;
            font-weight: 700;
            color: white;
            margin-bottom: 12px;
        }
        .pocket-stats { 
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 12px;
        }
        .pocket-stat { 
            font-size: 11px;
            color: rgba(255,255,255,0.9);
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }
        .pocket-stat strong { 
            display: block;
            font-size: 18px;
            color: white;
            margin-top: 5px;
        }
        .table-container { 
            flex: 1;
            overflow-y: auto;
            padding: 0 25px 25px 25px;
        }
        .pockets-table { 
            width: 100%;
            border-collapse: collapse;
            font-size: 13px;
        }
        .pockets-table th { 
            position: sticky;
            top: 0;
            background: #f8f9fa;
            z-index: 10;
            padding: 12px 8px;
            text-align: left;
            border-bottom: 2px solid #dee2e6;
            color: #667eea;
            font-weight: 700;
            text-transform: uppercase;
            font-size: 11px;
            letter-spacing: 0.5px;
        }
        .pockets-table td { 
            padding: 12px 8px;
            border-bottom: 1px solid #eee;
            cursor: pointer;
        }
        .pockets-table tr:hover { 
            background: #f8f9fa;
        }
        .pockets-table tr.top-pocket { 
            background: linear-gradient(to right, rgba(102, 126, 234, 0.1), transparent);
            font-weight: 600;
        }
        .viewer-panel { 
            background: rgba(255,255,255,0.95);
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            display: flex;
            flex-direction: column;
        }
        .viewer-header { 
            padding: 20px 25px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        .viewer-title { 
            font-size: 18px;
            font-weight: 600;
        }
        #viewport { 
            flex: 1;
            background: #1a1a2e;
            min-height: 600px;
        }
        .controls { 
            padding: 20px 25px;
            background: #f8f9fa;
            display: flex;
            gap: 10px;
            align-items: center;
            flex-wrap: wrap;
        }
        .btn { 
            padding: 10px 20px;
            border: 2px solid #dee2e6;
            background: white;
            border-radius: 8px;
            font-size: 13px;
            cursor: pointer;
            transition: all 0.3s;
            font-weight: 600;
        }
        .btn:hover { 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-color: #667eea;
            transform: translateY(-2px);
        }
        .btn.active { 
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-color: #667eea;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="header-content">
                <div class="title-section">
                    <div class="page-title">🎯 Pocket Discovery</div>
                    <div class="node-label">Node 5: Binding Site Analysis</div>
                </div>
                <div class="legend">
                    <div class="legend-item">
                        <span class="dot" style="background: #22c55e;"></span>
                        <span>Top Pocket</span>
                    </div>
                    <div class="legend-item">
                        <span class="dot" style="background: #fbbf24;"></span>
                        <span>Other Pockets</span>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="main-content">
            <div class="sidebar">
                <div class="section">
                    <div class="section-title">🏆 Top Binding Site</div>
                    <div class="top-pocket-card">
                        <div class="pocket-name">""" + top_pocket['pocket_name'] + """</div>
                        <div class="pocket-stats">
                            <div class="pocket-stat">Rank <strong>#""" + str(top_pocket.get('rank', 1)) + """</strong></div>
                            <div class="pocket-stat">Score <strong>""" + f"{top_pocket['score']:.2f}" + """</strong></div>
                            <div class="pocket-stat">Atoms <strong>""" + str(top_pocket.get('surface_atoms', 0)) + """</strong></div>
                            <div class="pocket-stat">Probability <strong>""" + f"{top_pocket.get('probability', 0):.0%}" + """</strong></div>
                        </div>
                    </div>
                </div>
                
                <div class="section">
                    <div class="section-title">📦 Docking Grid Box</div>
                    <div style="background: #fff3cd; padding: 15px; border-radius: 8px; border: 2px solid #fbbf24;">
                        <div class="pocket-stats">
                            <div class="pocket-stat" style="color: #856404;">Center X <strong style="color: #212529;">""" + f"{grid_params.get('center_x', 0):.2f}" + """</strong></div>
                            <div class="pocket-stat" style="color: #856404;">Size X <strong style="color: #212529;">""" + str(grid_params.get("size_x", "-")) + """</strong></div>
                            <div class="pocket-stat" style="color: #856404;">Center Y <strong style="color: #212529;">""" + f"{grid_params.get('center_y', 0):.2f}" + """</strong></div>
                            <div class="pocket-stat" style="color: #856404;">Size Y <strong style="color: #212529;">""" + str(grid_params.get("size_y", "-")) + """</strong></div>
                            <div class="pocket-stat" style="color: #856404;">Center Z <strong style="color: #212529;">""" + f"{grid_params.get('center_z', 0):.2f}" + """</strong></div>
                            <div class="pocket-stat" style="color: #856404;">Size Z <strong style="color: #212529;">""" + str(grid_params.get("size_z", "-")) + """</strong></div>
                        </div>
                    </div>
                </div>

                <div class="section-title" style="padding: 25px 25px 10px 25px;">📊 All Pockets (""" + str(data.get('total_pockets', 0)) + """)</div>
                <div class="table-container">
                    <table class="pockets-table">
                        <thead>
                            <tr>
                                <th>Name</th>
                                <th>Rank</th>
                                <th>Score</th>
                                <th>Prob</th>
                                <th>Res</th>
                                <th>Atoms</th>
                            </tr>
                        </thead>
                        <tbody>""" + pockets_rows + """
                        </tbody>
                    </table>
                </div>
            </div>
            
            <div class="viewer-panel">
                <div class="viewer-header">
                    <div class="viewer-title">🔬 3D Visualization</div>
                </div>
                <div id="viewport"></div>
                <div class="controls">
                    <button id="btn-cartoon" class="btn active" onclick="setVisual('cartoon')">Cartoon</button>
                    <button id="btn-backbone" class="btn" onclick="setVisual('backbone')">Backbone</button>
                    <button id="btn-surface" class="btn" onclick="setVisual('surface')">Surface</button>
                    <button id="btn-stick" class="btn" onclick="setVisual('ball+stick')">Ball+Stick</button>
                    <div style="width: 1px; height: 20px; background: #dee2e6; margin: 0 5px;"></div>
                    <button id="btn-top" class="btn" onclick="toggleTopPocket()">Top Pocket Only</button>
                    <button id="btn-grid" class="btn active" onclick="toggleGrid()">Toggle Grid Box</button>
                    <div style="flex: 1"></div>
                    <button class="btn" onclick="stage.autoView();">Reset View</button>
                </div>
            </div>
        </div>
    </div>
    
    <script>
        var stage = new NGL.Stage('viewport', {backgroundColor: '#1a1a2e'});
        var component;
        var shapeComp;
        var gridComp;
        var proteinData = """ + protein_pdb_js + """;
        var pockets = """ + pockets_json_str + """;
        var grid = """ + grid_js + """;
        
        function init() {
            if (!proteinData) return;
            var blob = new Blob([proteinData], {type: 'text/plain'});
            stage.loadFile(blob, {ext: 'pdb'}).then(function(comp) {
                component = comp;
                renderPockets();
                renderGrid();
                setVisual('cartoon');
                component.autoView();
            });
        }
        
        function renderPockets() {
            // Draw all pockets
            var shape = new NGL.Shape('pockets');
            pockets.forEach(function(p, i) {
                var color = (p.rank === 1) ? [0.13, 0.77, 0.34] : [0.98, 0.75, 0.14];
                var radius = 2.5 + (p.score / 15.0);
                if (radius > 6.0) radius = 6.0;
                shape.addSphere([p.x, p.y, p.z], color, radius, p.name);
            });
            shapeComp = stage.addComponentFromObject(shape);
            shapeComp.addRepresentation('buffer', {opacity: 0.7});
            
            // Draw ONLY top pocket (hidden by default, toggled later)
            var topShape = new NGL.Shape('top-pocket');
            var topP = pockets.find(p => p.rank === 1);
            if (topP) {
                var radius = 2.5 + (topP.score / 15.0);
                if (radius > 6.0) radius = 6.0;
                topShape.addSphere([topP.x, topP.y, topP.z], [0.13, 0.77, 0.34], radius, topP.name);
            }
            topPocketComp = stage.addComponentFromObject(topShape);
            topPocketComp.addRepresentation('buffer', {opacity: 0.7});
            topPocketComp.setVisibility(false);
        }

        function renderGrid() {
            if (!grid.center_x) return;
            var shape = new NGL.Shape("gridbox");
            
            var cx = grid.center_x, cy = grid.center_y, cz = grid.center_z;
            var sx = grid.size_x, sy = grid.size_y, sz = grid.size_z;
            
            var x = sx/2, y = sy/2, z = sz/2;
            var p = [[cx-x, cy-y, cz-z], [cx+x, cy-y, cz-z], [cx+x, cy+y, cz-z], [cx-x, cy+y, cz-z],
                     [cx-x, cy-y, cz+z], [cx+x, cy-y, cz+z], [cx+x, cy+y, cz+z], [cx-x, cy+y, cz+z]];
            var edges = [[0,1],[1,2],[2,3],[3,0], [4,5],[5,6],[6,7],[7,4], [0,4],[1,5],[2,6],[3,7]];
            var c = [1, 0.75, 0]; // Orange
            
            edges.forEach(e => shape.addLine(p[e[0]], p[e[1]], c));
            
            gridComp = stage.addComponentFromObject(shape);
            gridComp.addRepresentation('buffer', {linewidth: 4});
            gridComp.setVisibility(true);
        }

        function toggleGrid() {
            if (gridComp) {
                gridComp.setVisibility(!gridComp.visible);
                document.getElementById('btn-grid').classList.toggle('active', gridComp.visible);
            }
        }
        
        function toggleTopPocket() {
            var btn = document.getElementById('btn-top');
            var showOnlyTop = !btn.classList.contains('active');
            
            if (showOnlyTop) {
                btn.classList.add('active');
                shapeComp.setVisibility(false); // Hide all
                topPocketComp.setVisibility(true); // Show top only
            } else {
                btn.classList.remove('active');
                shapeComp.setVisibility(true); // Show all
                topPocketComp.setVisibility(false); // Hide top duplicate
            }
        }

        function setVisual(mode) {
            component.removeAllRepresentations();
            
            // Reset buttons (except grid and top pocket)
            var btns = document.querySelectorAll('.controls .btn');
            for (var i = 0; i < btns.length; i++) {
                if (btns[i].id !== 'btn-grid' && btns[i].id !== 'btn-top') {
                    btns[i].classList.remove('active');
                }
            }
            
            if (mode === 'cartoon') {
                component.addRepresentation('cartoon', {color: '#64748b', opacity: 1.0});
                document.getElementById('btn-cartoon').classList.add('active');
            } else if (mode === 'backbone') {
                component.addRepresentation('backbone', {color: 'lightgray', radiusScale: 2.0});
                document.getElementById('btn-backbone').classList.add('active');
            } else if (mode === 'surface') {
                component.addRepresentation('cartoon', {color: '#64748b', opacity: 1.0});
                component.addRepresentation('surface', {color: 'white', opacity: 0.15, probeRadius: 1.4, side: 'front'});
                document.getElementById('btn-surface').classList.add('active');
            } else if (mode === 'ball+stick') {
                component.addRepresentation('ball+stick', {colorScheme: 'element', multipleBond: true});
                document.getElementById('btn-stick').classList.add('active');
            }
        }
        
        function focusPocket(index) {
            var p = pockets[index];
            if (p) {
                stage.animationControls.zoomMove([p.x, p.y, p.z], -20, 1000);
            }
        }
        
        init();
        window.addEventListener('resize', function(){ stage.handleResize(); });
    </script>
</body>
</html>"""
    
    output_file = Path("outputs/report.html")
    output_file.parent.mkdir(exist_ok=True)
    with open(output_file, "w") as f:
        f.write(html_content)
    
    print(f"[Pocket Discovery] Report saved to {output_file}")

if __name__ == "__main__":
    generate_report()