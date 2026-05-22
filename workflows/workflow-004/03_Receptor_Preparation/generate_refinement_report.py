#!/usr/bin/env python3
"""
Node 3: Receptor Preparation Report
Creates interactive visualization of prepared receptor
"""

import json
from pathlib import Path

def generate_report():
    print("[Receptor Preparation] Generating HTML report from metadata...")
    
    # Read data from JSON
    json_file = Path("outputs/refined_receptor_metadata.json")
    if not json_file.exists():
        json_file = Path("refined_receptor_metadata.json")
    if not json_file.exists():
        print("[Receptor Preparation] ERROR: refined_receptor_metadata.json not found")
        return
    
    with open(json_file) as f:
        data = json.load(f)
    
    # Read PDBQT content for visualization
    pdbqt_path = Path("outputs") / data["output_pdbqt"]
    if not pdbqt_path.exists():
        pdbqt_path = Path(data["output_pdbqt"])
    if pdbqt_path.exists():
        with open(pdbqt_path) as f:
            pdbqt_content = f.read()
        pdbqt_content_js = json.dumps(pdbqt_content)
    else:
        pdbqt_content_js = '""'
    
    # Build preparation steps list
    steps_html = ""
    for i, step in enumerate(data["preparation_steps"], 1):
        steps_html += f"<li><span class='step-num'>{i}</span>{step}</li>"
    
    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Receptor Preparation Report</title>
    <script src="https://cdn.jsdelivr.net/npm/ngl@2.0.0-dev.37/dist/ngl.js"></script>
    <style>
        * { margin: 0; padding: 0; box-sizing: border-box; }
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
            padding: 20px;
        }
        .container { 
            max-width: 1400px; 
            margin: 0 auto; 
        }
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
            grid-template-columns: 400px 1fr;
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
        .section-title {
            font-size: 12px;
            font-weight: 700;
            text-transform: uppercase;
            color: #667eea;
            margin-bottom: 20px;
            letter-spacing: 1px;
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
        }
        .btn.active {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: #fff;
            border-color: #667eea;
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
        /* Steps List Styling */
        .steps-list {
            list-style: none;
            padding: 0;
            margin: 0;
        }
        .steps-list li {
            padding: 12px;
            margin-bottom: 8px;
            background: #f8f9fa;
            border-left: 3px solid #667eea;
            border-radius: 4px;
            font-size: 13px;
            display: flex;
            gap: 10px;
            align-items: center;
        }
        .step-num {
            background: #667eea;
            color: white;
            width: 20px;
            height: 20px;
            border-radius: 50%;
            display: flex;
            align-items: center;
            justify-content: center;
            font-size: 10px;
            font-weight: bold;
            flex-shrink: 0;
        }
        .info-row {
            display: flex;
            justify-content: space-between;
            padding: 10px 0;
            border-bottom: 1px solid #f1f3f5;
            font-size: 13px;
        }
        .info-label { color: #6c757d; }
        .info-value { font-weight: 600; color: #212529; }
        .status-badge {
            display: inline-block;
            padding: 4px 10px;
            border-radius: 12px;
            font-size: 11px;
            font-weight: 600;
            background: #22c55e;
            color: white;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="header-content">
                <div>
                    <h1>⚙️ Receptor Preparation</h1>
                    <div style="color: #6c757d; font-size: 14px; margin-top: 5px;">Refining Structure for Docking</div>
                </div>
                <div class="node-label">Node 3: Preparation</div>
            </div>
        </div>
        
        <div class="main-content">
            <div class="sidebar">
                <div class="section">
                    <div class="section-title">Status Report</div>
                    <div class="info-row">
                        <span class="info-label">Input File</span>
                        <span class="info-value">""" + data['input_pdb'] + """</span>
                    </div>
                    <div class="info-row">
                        <span class="info-label">Output File</span>
                        <span class="info-value">""" + data['output_pdbqt'] + """</span>
                    </div>
                    <div class="info-row">
                        <span class="info-label">Status</span>
                        <span class="status-badge">✓ """ + data['conversion_status'].upper() + """</span>
                    </div>
                </div>

                <div class="section">
                    <div class="section-title">File Size Statistics</div>
                    <div class="info-row">
                        <span class="info-label">Input PDB</span>
                        <span class="info-value">""" + f"{data['file_sizes']['input_pdb_bytes'] / 1024:.1f}" + """ KB</span>
                    </div>
                    <div class="info-row">
                        <span class="info-label">Cleaned PDB</span>
                        <span class="info-value">""" + f"{data['file_sizes']['fixed_pdb_bytes'] / 1024:.1f}" + """ KB</span>
                    </div>
                    <div class="info-row">
                        <span class="info-label">Final PDBQT</span>
                        <span class="info-value">""" + f"{data['file_sizes']['pdbqt_bytes'] / 1024:.1f}" + """ KB</span>
                    </div>
                </div>

                <div class="section">
                    <div class="section-title">Preparation Pipeline</div>
                    <ul class="steps-list">""" + steps_html + """</ul>
                </div>
            </div>
            
            <div class="viewer-panel">
                <div class="viewer-header">
                    <div class="viewer-title">🔬 Prepared Receptor Viewer</div>
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
        
        var pdbqtData = """ + pdbqt_content_js + """;
        var blob = new Blob([pdbqtData], {type: 'text/plain'});
        
        stage.loadFile(blob, {ext: 'pdbqt'}).then(function(comp) {
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
    
    print(f"[Receptor Preparation] Report saved to {output_file}")

if __name__ == "__main__":
    generate_report()