#!/usr/bin/env python3
"""
Node 2: Ligand Report Generation
Creates interactive HTML visualization of downloaded ligands
"""

import json
from pathlib import Path

def generate_report():
    print("[Ligand Collection] Generating HTML report from metadata...")
    
    # Read data from JSON
    json_file = Path("outputs/ligands_metadata.json")
    if not json_file.exists():
        json_file = Path("ligands_metadata.json")
    if not json_file.exists():
        print("[Ligand Collection] ERROR: ligands_metadata.json not found")
        return
    
    with open(json_file) as f:
        data = json.load(f)
    
    # Build ligand cards with NGL viewers
    ligand_cards = ""
    viewer_scripts = ""
    
    for i, ligand in enumerate(data["ligands"]):
        if ligand["download_status"] == "success":
            sdf_file = ligand["sdf_file"]
            cid = ligand["cid"]
            file_size = ligand["file_size_bytes"]
            
            # Read SDF content and embed it
            sdf_path = Path("outputs") / sdf_file
            if not sdf_path.exists():
                sdf_path = Path(sdf_file)
            if sdf_path.exists():
                with open(sdf_path) as f:
                    sdf_content = f.read()
                sdf_content_js = json.dumps(sdf_content)
            else:
                sdf_content_js = '""'
            
            ligand_cards += f"""
            <div class="ligand-card">
                <div class="card-header">
                    <div class="cid-label">CID {cid}</div>
                    <span class="status-badge success">✓ Downloaded</span>
                </div>
                <div id="viewport-{i}" class="viewport"></div>
                <div class="card-footer">
                    <div class="info-item">
                        <span class="label">File</span>
                        <span class="value">{sdf_file}</span>
                    </div>

                </div>
            </div>"""
            
            # NGL Viewer script
            viewer_scripts += f"""
            (function() {{
                var stage{i} = new NGL.Stage('viewport-{i}', {{backgroundColor: '#1a1a2e'}});
                var sdfData{i} = {sdf_content_js};
                var blob{i} = new Blob([sdfData{i}], {{type: 'text/plain'}});
                stage{i}.loadFile(blob{i}, {{ext: 'sdf'}}).then(function(component) {{
                    component.addRepresentation('ball+stick', {{colorScheme: 'element'}});
                    component.autoView();
                }});
            }})();
            """
        else:
            ligand_cards += f"""
            <div class="ligand-card error">
                <div class="card-header">
                    <div class="cid-label">CID {ligand['cid']}</div>
                    <span class="status-badge error">✗ Failed</span>
                </div>
                <div class="error-message">
                    Download failed: {ligand.get('error', 'Unknown error')}
                </div>
            </div>"""
    
    html_content = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Ligand Collection Report</title>
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
            box-shadow: 0 4px 15px rgba(102, 126, 234, 0.4);
        }
        .stats-summary {
            display: flex;
            gap: 30px;
        }
        .stat-item {
            text-align: center;
        }
        .stat-number {
            font-size: 36px;
            font-weight: 700;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            -webkit-background-clip: text;
            -webkit-text-fill-color: transparent;
        }
        .stat-label {
            font-size: 11px;
            color: #6c757d;
            text-transform: uppercase;
            letter-spacing: 1px;
            margin-top: 5px;
        }
        .main-content {
            padding: 20px 0;
        }
        .ligands-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(350px, 1fr));
            gap: 30px;
            justify-content: center;
        }
        .ligand-card {
            background: rgba(255,255,255,0.95);
            border-radius: 12px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            overflow: hidden;
            transition: transform 0.3s, box-shadow 0.3s;
        }
        .ligand-card:hover {
            transform: translateY(-5px);
            box-shadow: 0 15px 50px rgba(0,0,0,0.3);
        }
        .ligand-card.error {
            border: 2px solid #dc3545;
        }
        .card-header {
            padding: 20px;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            display: flex;
            justify-content: space-between;
            align-items: center;
        }
        .cid-label {
            font-size: 20px;
            font-weight: 700;
            color: white;
            font-family: 'Courier New', monospace;
        }
        .status-badge {
            padding: 6px 14px;
            border-radius: 20px;
            font-size: 12px;
            font-weight: 600;
        }
        .status-badge.success {
            background: rgba(255,255,255,0.3);
            color: white;
        }
        .status-badge.error {
            background: #dc3545;
            color: white;
        }
        .viewport {
            width: 100%;
            height: 320px;
            background: #1a1a2e;
        }
        .card-footer {
            padding: 20px;
            background: #f8f9fa;
        }
        .info-item {
            display: flex;
            justify-content: space-between;
            padding: 8px 0;
            border-bottom: 1px solid #e9ecef;
        }
        .info-item:last-child {
            border-bottom: none;
        }
        .info-item .label {
            font-size: 13px;
            color: #6c757d;
        }
        .info-item .value {
            font-size: 14px;
            font-weight: 600;
            color: #212529;
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
            <div class="header-content">
                <div class="title-section">
                    <div class="page-title">💊 Ligand Structures</div>
                    <div class="node-label">Node 2: Ligand Collection</div>
                </div>
                <div class="stats-summary">
                    <div class="stat-item">
                        <div class="stat-number">""" + str(data['total_count']) + """</div>
                        <div class="stat-label">Total</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-number">""" + str(data['successful_downloads']) + """</div>
                        <div class="stat-label">Success</div>
                    </div>
                    <div class="stat-item">
                        <div class="stat-number">""" + str(data['failed_downloads']) + """</div>
                        <div class="stat-label">Failed</div>
                    </div>
                </div>
            </div>
        </div>
        
        <div class="main-content">
            <div class="ligands-grid">""" + ligand_cards + """
            </div>
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
    
    print(f"[Ligand Collection] Report saved to {output_file}")

if __name__ == "__main__":
    generate_report()