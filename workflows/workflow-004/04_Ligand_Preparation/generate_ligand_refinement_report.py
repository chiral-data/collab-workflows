#!/usr/bin/env python3
"""
Node 4: Ligand Preparation Report
Creates interactive visualization of prepared ligands
"""

import json
from pathlib import Path

def generate_report():
    print("[Ligand Preparation] Generating HTML report from metadata...")
    
    # Read data from JSON
    json_file = Path("outputs/refined_ligands_metadata.json")
    if not json_file.exists():
        json_file = Path("refined_ligands_metadata.json")
    if not json_file.exists():
        print("[Ligand Preparation] ERROR: refined_ligands_metadata.json not found")
        return
    
    with open(json_file) as f:
        data = json.load(f)
    
    # Build ligand cards with NGL viewers
    cards_html = ""
    viewer_scripts = ""
    
    for i, ligand in enumerate(data["ligands"]):
        if ligand["conversion_status"] == "success":
            pdbqt_file = ligand["output_pdbqt"]
            name = ligand["ligand_name"]
            
            # Read PDBQT content and embed it
            pdbqt_path = Path("outputs") / pdbqt_file
            if not pdbqt_path.exists():
                pdbqt_path = Path(pdbqt_file)
            if pdbqt_path.exists():
                with open(pdbqt_path) as f:
                    pdbqt_content = f.read()
                pdbqt_content_js = json.dumps(pdbqt_content)
            else:
                pdbqt_content_js = '""'
            
            cards_html += f"""
            <div class="ligand-card">
                <div class="card-header">
                    <h3>{name}</h3>
                    <span class="badge success">✓ Success</span>
                </div>
                <div class="card-note">
                    📊 PDBQT format with {ligand['num_rotatable_bonds']} rotatable bonds
                </div>
                <div class="viewer-container">
                    <div id='viewport-{i}' class='ngl-viewport'></div>
                </div>
                <div class="card-stats">
                    <div class="stat-item">
                        <span class="stat-label">Atoms</span>
                        <span class="stat-value">{ligand['num_atoms']}</span>
                    </div>
                    <div class="stat-item">
                        <span class="stat-label">Rotatable Bonds</span>
                        <span class="stat-value">{ligand['num_rotatable_bonds']}</span>
                    </div>

                </div>
            </div>"""
            
            # NGL Viewer script for this ligand
            viewer_scripts += f"""
            (function() {{
                var stage{i} = new NGL.Stage('viewport-{i}', {{backgroundColor: '#1a1a2e'}});
                var pdbqtData{i} = {pdbqt_content_js};
                var blob{i} = new Blob([pdbqtData{i}], {{type: 'text/plain'}});
                stage{i}.loadFile(blob{i}, {{ext: 'pdbqt'}}).then(function(component) {{
                    component.addRepresentation('licorice', {{colorScheme: 'element'}});
                    component.autoView();
                }});
            }})();
            """
        else:
            cards_html += f"""
            <div class="ligand-card error-card">
                <div class="card-header">
                    <h3>{ligand['ligand_name']}</h3>
                    <span class="badge error">✗ Failed</span>
                </div>
                <p class="error-message">Error: {ligand.get('error', 'Unknown error')}</p>
            </div>"""
    
    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Ligand Preparation Report</title>
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
        .header h1 {
            font-size: 32px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            margin: 0 0 10px 0;
        }
        .node-label {
            display: inline-block;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 8px 16px;
            border-radius: 20px;
            font-size: 14px;
            font-weight: 600;
            margin-top: 10px;
        }
        .info-note {
            background: rgba(255,255,255,0.95);
            border-left: 4px solid #667eea;
            padding: 20px;
            margin-bottom: 20px;
            border-radius: 8px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }
        .summary {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 20px;
            margin-bottom: 20px;
        }
        .summary-card {
            background: rgba(255,255,255,0.95);
            padding: 25px;
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            text-align: center;
        }
        .summary-card h3 {
            margin: 0 0 15px 0;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
            font-size: 14px;
            text-transform: uppercase;
            letter-spacing: 1px;
        }
        .summary-card .value {
            font-size: 42px;
            font-weight: bold;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        .ligands-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(380px, 1fr));
            gap: 20px;
        }
        .ligand-card {
            background: rgba(255,255,255,0.95);
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            overflow: hidden;
            transition: transform 0.3s;
        }
        .ligand-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 15px 50px rgba(0,0,0,0.3);
        }
        .card-header {
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        .card-header h3 {
            margin: 0;
            color: white;
            font-size: 18px;
        }
        .card-note {
            padding: 12px 20px;
            background: #e7f1ff;
            border-bottom: 1px solid #b3d9ff;
            color: #004085;
            font-size: 13px;
        }
        .viewer-container {
            padding: 15px;
        }
        .ngl-viewport {
            width: 100%;
            height: 300px;
            border-radius: 8px;
            overflow: hidden;
        }
        .card-stats {
            display: grid;
            grid-template-columns: repeat(2, 1fr);
            gap: 10px;
            padding: 20px;
            background: #f8f9fa;
        }
        .stat-item {
            text-align: center;
            padding: 15px;
            background: white;
            border-radius: 8px;
        }
        .stat-label {
            display: block;
            font-size: 11px;
            color: #6c757d;
            text-transform: uppercase;
            margin-bottom: 8px;
            letter-spacing: 0.5px;
        }
        .stat-value {
            display: block;
            font-size: 20px;
            font-weight: bold;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        .badge {
            padding: 6px 14px;
            border-radius: 20px;
            font-size: 12px;
            font-weight: 600;
        }
        .badge.success {
            background: rgba(255,255,255,0.3);
            color: white;
        }
        .badge.error {
            background: #dc3545;
            color: white;
        }
        .error-card {
            border: 2px solid #dc3545;
        }
        .error-message {
            padding: 30px;
            color: #dc3545;
            text-align: center;
            font-weight: 600;
        }
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>⚗️ Ligand Preparation</h1>
            <p>Molecular Structure Processing for Docking</p>
            <div class="node-label">Node 4: Ligand Preparation</div>
        </div>
        
        <div class="info-note">
            <strong>ℹ️ PDBQT Format:</strong> Ligands converted to PDBQT with preserved hydrogens, 
            rotatable bonds identified, and charges assigned for accurate docking simulations.
        </div>
        
        <div class="summary">
            <div class="summary-card">
                <h3>Total Ligands</h3>
                <div class="value">""" + str(data['total_count']) + """</div>
            </div>
            <div class="summary-card">
                <h3>Successful</h3>
                <div class="value">""" + str(data['successful_conversions']) + """</div>
            </div>
            <div class="summary-card">
                <h3>Failed</h3>
                <div class="value">""" + str(data['failed_conversions']) + """</div>
            </div>
        </div>
        
        <div class="ligands-grid">""" + cards_html + """
        </div>
    </div>
    
    <script>
        document.addEventListener('DOMContentLoaded', function() {
            """ + viewer_scripts + """
        });
    </script>
</body>
</html>"""
    
    output_file = Path("outputs/report.html")
    output_file.parent.mkdir(exist_ok=True)
    with open(output_file, "w") as f:
        f.write(html_content)
    
    print(f"[Ligand Preparation] Report saved to {output_file}")

if __name__ == "__main__":
    generate_report()