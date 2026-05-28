import json
import os

def generate_report():
    print("Generating Report for Node 4...")
    if not os.path.exists("outputs/data.json"):
        return

    with open("outputs/data.json", "r") as f:
        data = json.load(f)

    js_data = json.dumps(data)

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <title>Predicted Binding Affinities - Interactive Dashboard</title>
    <meta charset="UTF-8">
    <script src="https://unpkg.com/@rdkit/rdkit/Code/MinimalLib/dist/RDKit_minimal.js"></script>
    <style>
        * {{ box-sizing: border-box; }}
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            padding: 0;
            margin: 0;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }}
        .header {{
            background: rgba(255,255,255,0.95);
            padding: 25px 30px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            position: sticky;
            top: 0;
            z-index: 100;
        }}
        .header h1 {{
            margin: 0 0 15px 0;
            color: #2c3e50;
            font-size: 2em;
        }}
        .controls {{
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
            align-items: center;
        }}
        .control-group {{
            display: flex;
            align-items: center;
            gap: 10px;
        }}
        .control-group label {{
            font-weight: 600;
            color: #555;
        }}
        select, button {{
            padding: 10px 15px;
            border: 2px solid #ddd;
            border-radius: 6px;
            background: white;
            cursor: pointer;
            font-size: 14px;
            transition: all 0.2s;
        }}
        select:hover, button:hover {{ border-color: #667eea; }}
        button {{
            background: #667eea;
            color: white;
            border-color: #667eea;
            font-weight: 600;
        }}
        button:hover {{ background: #5568d3; }}
        .stats {{
            display: flex;
            gap: 20px;
            margin-top: 15px;
            padding-top: 15px;
            border-top: 2px solid #eee;
        }}
        .stat {{ display: flex; align-items: center; gap: 8px; }}
        .stat-label {{ color: #666; font-size: 14px; }}
        .stat-value {{ font-weight: bold; font-size: 18px; color: #667eea; }}
        .container {{
            padding: 30px;
            max-width: 1600px;
            margin: 0 auto;
        }}
        .grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
            gap: 25px;
        }}
        .card {{
            background: white;
            border-radius: 15px;
            overflow: hidden;
            box-shadow: 0 10px 30px rgba(0,0,0,0.2);
            transition: transform 0.3s ease, box-shadow 0.3s ease;
        }}
        .card:hover {{
            transform: translateY(-5px);
            box-shadow: 0 15px 40px rgba(0,0,0,0.3);
        }}
        .mol-canvas {{
            width: 100%;
            height: 220px;
            background: #f8f9fa;
            display: flex;
            align-items: center;
            justify-content: center;
        }}
        .mol-canvas svg {{
            width: 100%;
            height: 100%;
        }}
        .card-content {{ padding: 20px; }}
        h3 {{
            font-size: 1.3em;
            margin: 0 0 15px 0;
            color: #2c3e50;
            word-wrap: break-word;
        }}
        .smiles-box {{
            background: #f8f9fa;
            border: 2px dashed #dee2e6;
            border-radius: 8px;
            padding: 12px;
            margin: 15px 0;
            position: relative;
            cursor: pointer;
            transition: all 0.2s ease;
            overflow-x: auto;
            max-height: 100px;
            overflow-y: auto;
        }}
        .smiles-box:hover {{
            background: #e9ecef;
            border-color: #667eea;
        }}
        .smiles-box code {{
            font-family: 'Courier New', monospace;
            font-size: 0.85em;
            color: #495057;
            word-break: break-all;
            user-select: all;
        }}
        .copy-indicator {{
            position: absolute;
            top: 8px;
            right: 8px;
            font-size: 1.2em;
            opacity: 0.5;
        }}
        .smiles-box:hover .copy-indicator {{ opacity: 1; }}
        .info-row {{
            margin: 10px 0;
            padding: 8px 0;
            border-bottom: 1px solid #e9ecef;
        }}
        .info-row:last-child {{ border-bottom: none; }}
        .prediction {{ font-size: 1.6em; color: #667eea; font-weight: bold; }}
        .ad-badge {{
            display: inline-block;
            padding: 6px 14px;
            border-radius: 12px;
            font-weight: bold;
            font-size: 0.9em;
        }}
        .ad-in {{ background: #d4edda; color: #155724; }}
        .ad-out {{ background: #f8d7da; color: #721c24; }}
        .toast {{
            position: fixed;
            bottom: 30px;
            right: 30px;
            background: #28a745;
            color: white;
            padding: 15px 25px;
            border-radius: 8px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.3);
            opacity: 0;
            transition: opacity 0.3s ease;
            z-index: 1000;
            font-weight: bold;
        }}
        .toast.show {{ opacity: 1; }}
        .hidden {{ display: none !important; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Predicted Binding Affinities - Interactive Dashboard</h1>
        <div class="controls">
            <div class="control-group">
                <label for="filterAD">Filter by AD:</label>
                <select id="filterAD" onchange="applyFilters()">
                    <option value="all">All Compounds</option>
                    <option value="IN">IN (Within Domain)</option>
                    <option value="OUT">OUT (Outside Domain)</option>
                </select>
            </div>
            <div class="control-group">
                <label for="sortBy">Sort by Score:</label>
                <select id="sortBy" onchange="applyFilters()">
                    <option value="default">Default Order</option>
                    <option value="high">Highest First</option>
                    <option value="low">Lowest First</option>
                </select>
            </div>
            <button onclick="resetFilters()">Reset All</button>
        </div>
        <div class="stats">
            <div class="stat">
                <span class="stat-label">Total:</span>
                <span class="stat-value" id="totalCount">0</span>
            </div>
            <div class="stat">
                <span class="stat-label">Showing:</span>
                <span class="stat-value" id="visibleCount">0</span>
            </div>
            <div class="stat">
                <span class="stat-label">IN Domain:</span>
                <span class="stat-value" id="inCount">0</span>
            </div>
            <div class="stat">
                <span class="stat-label">OUT Domain:</span>
                <span class="stat-value" id="outCount">0</span>
            </div>
        </div>
    </div>
    <div class="container">
        <div class="grid" id="cardsGrid"></div>
    </div>
    <div id="toast" class="toast">SMILES copied to clipboard!</div>

    <script>
        const allData = {js_data};
        let filteredData = [...allData];
        let RDKitInstance = null;

        let observer = null;

        window.initRDKitModule().then(RDKit => {{
            RDKitInstance = RDKit;
            renderCards();
        }});

        function renderMol(placeholder) {{
            if (!RDKitInstance) return;
            const smiles = placeholder.dataset.smiles;
            try {{
                const mol = RDKitInstance.get_mol(smiles);
                if (mol) {{
                    placeholder.innerHTML = mol.get_svg();
                    mol.delete();
                }} else {{
                    placeholder.textContent = 'No structure';
                    placeholder.style.color = '#aaa';
                }}
            }} catch (e) {{
                placeholder.textContent = 'No structure';
                placeholder.style.color = '#aaa';
            }}
        }}

        function createCard(item, index) {{
            const drugName = item.drug_name || `Compound ${{index + 1}}`;
            const adClass = item.ad_status === 'IN' ? 'ad-in' : 'ad-out';
            const escapedSmiles = item.smiles.replace(/\\/g, '\\\\').replace(/'/g, "\\'");

            return `
                <div class="card" data-ad="${{item.ad_status}}" data-score="${{item.prediction}}">
                    <div class="mol-canvas mol-placeholder" data-smiles="${{item.smiles.replace(/"/g, '&quot;')}}"></div>
                    <div class="card-content">
                        <h3>${{drugName}}</h3>
                        <div class="smiles-box" onclick="copyToClipboard('${{escapedSmiles}}')" title="Click to copy SMILES">
                            <code>${{item.smiles}}</code>
                            <span class="copy-indicator">&#128203;</span>
                        </div>
                        <div class="info-row">
                            <strong>Prediction:</strong>
                            <span class="prediction">${{item.prediction.toFixed(4)}}</span>
                        </div>
                        <div class="info-row">
                            <strong>Applicability Domain:</strong>
                            <span class="ad-badge ${{adClass}}">${{item.ad_status}}</span>
                        </div>
                    </div>
                </div>
            `;
        }}

        function renderCards() {{
            const grid = document.getElementById('cardsGrid');
            grid.innerHTML = filteredData.map((item, idx) => createCard(item, idx)).join('');
            updateStats();
            observePlaceholders();
        }}

        function observePlaceholders() {{
            if (observer) observer.disconnect();
            observer = new IntersectionObserver((entries) => {{
                entries.forEach(entry => {{
                    if (entry.isIntersecting) {{
                        renderMol(entry.target);
                        observer.unobserve(entry.target);
                    }}
                }});
            }}, {{ rootMargin: '200px' }});
            document.querySelectorAll('.mol-placeholder').forEach(el => observer.observe(el));
        }}

        function applyFilters() {{
            const adFilter = document.getElementById('filterAD').value;
            const sortBy = document.getElementById('sortBy').value;

            filteredData = allData.filter(item => adFilter === 'all' || item.ad_status === adFilter);

            if (sortBy === 'high') {{
                filteredData.sort((a, b) => b.prediction - a.prediction);
            }} else if (sortBy === 'low') {{
                filteredData.sort((a, b) => a.prediction - b.prediction);
            }}

            renderCards();
        }}

        function resetFilters() {{
            document.getElementById('filterAD').value = 'all';
            document.getElementById('sortBy').value = 'default';
            filteredData = [...allData];
            renderCards();
        }}

        function updateStats() {{
            const inCount = allData.filter(item => item.ad_status === 'IN').length;
            const outCount = allData.filter(item => item.ad_status === 'OUT').length;
            document.getElementById('totalCount').textContent = allData.length;
            document.getElementById('visibleCount').textContent = filteredData.length;
            document.getElementById('inCount').textContent = inCount;
            document.getElementById('outCount').textContent = outCount;
        }}

        function copyToClipboard(text) {{
            navigator.clipboard.writeText(text).then(() => {{
                const toast = document.getElementById('toast');
                toast.classList.add('show');
                setTimeout(() => toast.classList.remove('show'), 2000);
            }}).catch(err => console.error('Failed to copy:', err));
        }}
    </script>
</body>
</html>
"""

    with open("outputs/report.html", "w", encoding="utf-8") as f:
        f.write(html_content)

    print("Interactive dashboard generated at outputs/report.html")

if __name__ == "__main__":
    generate_report()
