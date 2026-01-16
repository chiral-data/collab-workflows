#!/usr/bin/env python3
"""
Node 1: Receptor Report Generation
Creates interactive HTML visualization of downloaded receptor
"""

import json
from pathlib import Path

def generate_report():
    print("[Receptor Acquisition] Generating HTML report from metadata...")
    
    # Read data from JSON
    json_file = Path("outputs/receptor_metadata.json")
    if not json_file.exists():
        json_file = Path("receptor_metadata.json")
    if not json_file.exists():
        print("[Receptor Acquisition] ERROR: receptor_metadata.json not found")
        return
    
    with open(json_file) as f:
        data = json.load(f)
    
    pdb_file = data["pdb_file"]
    
    # Read PDB content for NGL viewer
    pdb_path = Path(pdb_file)
    if not pdb_path.exists() and Path("outputs") / pdb_file:
        pdb_path = Path("outputs") / pdb_file
    if not pdb_path.exists():
        print(f"[Receptor Acquisition] ERROR: {pdb_file} not found")
        return
    
    with open(pdb_path) as f:
        pdb_content = f.read()
    
    # Escape for JavaScript
    pdb_content_js = json.dumps(pdb_content)
    
    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Protein Structure - """ + data['pdb_id'] + """</title>
    <script src="https://cdn.jsdelivr.net/npm/ngl@2.0.0-dev.37/dist/ngl.js"></script>
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, Oxygen, Ubuntu, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: #212529;
            min-height: 100vh;
            padding: 20px;
        }
        .header {
            background: rgba(255,255,255,0.95);
            border-radius: 12px;
            padding: 30px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            margin-bottom: 20px;
        }
        .container {
            max-width: 1400px;
            margin: 0 auto;
        }
        .header-content {
            display: flex;
            align-items: center;
            gap: 20px;
        }
        .pdb-id {
            font-size: 42px;
            font-weight: 700;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            font-family: 'Courier New', monospace;
        }
        .node-label {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 8px 16px;
            border-radius: 20px;
            font-size: 14px;
            font-weight: 600;
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);
        }
        .main-content {
            display: grid;
            grid-template-columns: 350px 1fr;
            gap: 20px;
        }
        .sidebar {
            background: rgba(255,255,255,0.95);
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            height: fit-content;
        }
        .section {
            padding: 25px;
            border-bottom: 1px solid #e9ecef;
        }
        .section:last-child {
            border-bottom: none;
        }
        .section-title {
            font-size: 12px;
            font-weight: 700;
            text-transform: uppercase;
            color: #667eea;
            margin-bottom: 20px;
            letter-spacing: 1px;
        }
        .info-row {
            display: flex;
            justify-content: space-between;
            padding: 12px 0;
            border-bottom: 1px solid #f1f3f5;
        }
        .info-row:last-child {
            border-bottom: none;
        }
        .info-label {
            font-size: 13px;
            color: #6c757d;
        }
        .info-value {
            font-size: 14px;
            font-weight: 600;
            color: #212529;
        }
        .viewer-panel {
            background: rgba(255,255,255,0.95);
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            overflow: hidden;
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
            width: 100%;
            height: 600px;
            background: #000;
        }
        .controls {
            padding: 20px 25px;
            background: #f8f9fa;
            display: flex;
            gap: 10px;
            flex-wrap: wrap;
        }
        .btn {
            padding: 10px 20px;
            border: 2px solid #dee2e6;
            background: #fff;
            border-radius: 8px;
            font-size: 14px;
            cursor: pointer;
            transition: all 0.3s;
            font-weight: 600;
        }
        .btn:hover {
            background: #667eea;
            color: white;
            border-color: #667eea;
            transform: translateY(-2px);
            box-shadow: 0 5px 15px rgba(102, 126, 234, 0.3);
        }
        .btn.active {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: #fff;
            border-color: #667eea;
        }
        .stats-grid {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 15px;
        }
        .stat-box {
            text-align: center;
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            border-radius: 10px;
            box-shadow: 0 5px 15px rgba(102, 126, 234, 0.3);
        }
        .stat-number {
            font-size: 32px;
            font-weight: 700;
            color: white;
        }
        .stat-label {
            font-size: 11px;
            color: rgba(255,255,255,0.9);
            margin-top: 5px;
            text-transform: uppercase;
            letter-spacing: 1px;
        }
        .control-group {
            display: flex;
            align-items: center;
            gap: 8px;
            padding-right: 20px;
            border-right: 1px solid #dee2e6;
        }
        .control-group:last-of-type {
            border-right: none;
        }
        .control-label {
            font-size: 11px;
            font-weight: 700;
            text-transform: uppercase;
            color: #6c757d;
            margin-right: 5px;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="header-content">
                <div class="pdb-id">""" + data['pdb_id'] + """</div>
                <div class="node-label">Node 1: Receptor Acquisition</div>
            </div>
        </div>
        
        <div class="main-content">
            <div class="sidebar">
                <div class="section">
                    <div class="section-title">Structure Summary</div>
                    <div class="stats-grid">
                        <div class="stat-box">
                            <div class="stat-number">""" + f"{data['num_atoms']:,}" + """</div>
                            <div class="stat-label">Atoms</div>
                        </div>
                        <div class="stat-box">
                            <div class="stat-number">""" + f"{data['num_residues']:,}" + """</div>
                            <div class="stat-label">Residues</div>
                        </div>
                        <div class="stat-box">
                            <div class="stat-number">""" + str(data['num_chains']) + """</div>
                            <div class="stat-label">Chains</div>
                        </div>
                        <div class="stat-box">
                            <div class="stat-number">""" + f"{data['file_size_bytes'] / 1024:.1f}" + """</div>
                            <div class="stat-label">KB</div>
                        </div>
                    </div>
                </div>
                
                <div class="section">
                    <div class="section-title">File Information</div>
                    <div class="info-row">
                        <span class="info-label">PDB ID</span>
                        <span class="info-value">""" + data['pdb_id'] + """</span>
                    </div>
                    <div class="info-row">
                        <span class="info-label">Filename</span>
                        <span class="info-value">""" + data['pdb_file'] + """</span>
                    </div>
                    <div class="info-row">
                        <span class="info-label">Downloaded</span>
                        <span class="info-value">""" + data['download_timestamp'].split('T')[0] + """</span>
                    </div>
                </div>
            </div>
            
            <div class="viewer-panel">
                <div class="viewer-header">
                    <div class="viewer-title">🔬 3D Structure Viewer</div>
                </div>
                <div id="viewport"></div>
                <div class="controls">
                    <div class="control-group">
                        <span class="control-label">Representation</span>
                        <button class="btn active" onclick="setRep('cartoon')">Cartoon</button>
                        <button class="btn" onclick="setRep('backbone')">Backbone</button>
                        <button class="btn" onclick="setRep('surface')">Surface</button>
                        <button class="btn" onclick="setRep('ball+stick')">Ball+Stick</button>
                    </div>
                    <div class="control-group">
                        <span class="control-label">Color Scheme</span>
                        <button class="btn active" id="btn-col-chain" onclick="setColor('chainname')">Chain</button>
                        <button class="btn" id="btn-col-elem" onclick="setColor('element')">Element</button>
                        <button class="btn" id="btn-col-bfactor" onclick="setColor('bfactor')">B-Factor</button>
                    </div>
                    <button class="btn" style="margin-left: auto;" onclick="stage.autoView()">Reset View</button>
                </div>
            </div>
        </div>
    </div>
    
    <script>
        var stage = new NGL.Stage('viewport', {backgroundColor: '#1a1a2e'});
        var component;
        var currentRep = 'cartoon';
        var currentColor = 'chainname';
        
        var pdbData = """ + pdb_content_js + """;
        var blob = new Blob([pdbData], {type: 'text/plain'});
        
        stage.loadFile(blob, {ext: 'pdb'}).then(function(comp) {
            component = comp;
            updateVisual();
            component.autoView();
        });
        
        function setRep(type) {
            currentRep = type;
            updateVisual();
            updateButtons();
        }
        
        function setColor(scheme) {
            currentColor = scheme;
            updateVisual();
            updateButtons();
        }
        
        function updateVisual() {
            if (!component) return;
            component.removeAllRepresentations();
            
            var params = { color: currentColor };
            if (currentRep === 'surface') {
                params = { color: currentColor, opacity: 0.7, probeRadius: 1.4 };
            } else if (currentRep === 'ball+stick') {
                params = { colorScheme: currentColor }; 
            }
            
            component.addRepresentation(currentRep, params);
        }
        
        function updateButtons() {
            // Update Representation buttons
            document.querySelectorAll('.control-group:first-child .btn').forEach(btn => {
                if (btn.innerText.toLowerCase().replace('&','').replace(' ','').includes(currentRep.replace('+',''))) 
                    btn.classList.add('active');
                else 
                    btn.classList.remove('active');
            });
            
            // Update Color buttons
            document.getElementById('btn-col-chain').classList.toggle('active', currentColor === 'chainname');
            document.getElementById('btn-col-elem').classList.toggle('active', currentColor === 'element');
            document.getElementById('btn-col-bfactor').classList.toggle('active', currentColor === 'bfactor');
        }
        
        window.addEventListener('resize', function() {
            stage.handleResize();
        });
    </script>
</body>
</html>"""
    
    output_file = Path("outputs/report.html")
    output_file.parent.mkdir(exist_ok=True)
    with open(output_file, "w") as f:
        f.write(html_content)
    
    print(f"[Receptor Acquisition] Report saved to {output_file}")

if __name__ == "__main__":
    generate_report()