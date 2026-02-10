import json
import os
import sys

def generate_report():
    print("Generating Interactive Dashboard for Node 4...")
    json_path = "outputs/data.json"
    if not os.path.exists(json_path):
        print("Error: outputs/data.json not found. Run 7_predict_from_csv.py first.")
        sys.exit(1)

    with open(json_path, "r") as f:
        data = json.load(f)

    # Unpack data
    summary = data.get("summary", {})
    predictions = data.get("predictions", [])
    
    # Sort predictions by score (best/lowest first) by default
    predictions.sort(key=lambda x: x['score'])

    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Node 4: Prediction Results</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
        <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&display=swap" rel="stylesheet">
        <style>
            :root {{ 
                --primary: #4338ca;       /* Indigo 700 */
                --primary-light: #6366f1; /* Indigo 500 */
                --bg-body: #f8fafc;       /* Slate 50 */
                --bg-card: #ffffff;
                --text-main: #1e293b;     /* Slate 800 */
                --text-muted: #64748b;    /* Slate 500 */
                --border: #e2e8f0;
                --shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);
                --success: #10b981;
                --danger: #ef4444;
            }}

            body {{ 
                font-family: 'Outfit', sans-serif; 
                background: var(--bg-body); 
                color: var(--text-main); 
                margin: 0; 
                padding: 0; 
                height: 100vh; 
                display: flex; 
                flex-direction: column; 
                overflow: hidden;
            }}
            
            /* TOP NAVIGATION BAR */
            .navbar {{ 
                background: var(--bg-card); 
                border-bottom: 1px solid var(--border); 
                padding: 0 2rem; 
                height: 64px; 
                display: flex; 
                align-items: center; 
                justify-content: space-between;
                box-shadow: 0 1px 2px 0 rgba(0, 0, 0, 0.05);
                z-index: 50;
            }}

            .brand {{ display: flex; align-items: center; gap: 1rem; }}
            .brand h1 {{ margin: 0; font-size: 1.25rem; font-weight: 700; color: var(--text-main); letter-spacing: -0.025em; }}
            .brand .badge {{ 
                padding: 0.25rem 0.75rem; 
                border-radius: 9999px; 
                font-size: 0.75rem; 
                font-weight: 600; 
                background: var(--primary); 
                color: white; 
            }}
            
            /* MAIN LAYOUT */
            .main {{ 
                flex: 1; 
                padding: 2rem; 
                overflow-y: auto; 
                max-width: 1600px; 
                width: 100%; 
                margin: 0 auto; 
                box-sizing: border-box;
            }}
            
            /* CONTROLS BAR */
            .controls-bar {{
                display: flex;
                gap: 1rem;
                margin-bottom: 1rem;
                background: white;
                padding: 1rem;
                border-radius: 0.75rem;
                border: 1px solid var(--border);
                box-shadow: var(--shadow);
                align-items: center;
                flex-wrap: wrap;
            }}
            
            .control-group {{ display: flex; align-items: center; gap: 0.5rem; }}
            .control-label {{ font-weight: 600; font-size: 0.9rem; color: var(--text-muted); }}
            
            select, input {{
                padding: 0.5rem 0.75rem;
                border: 1px solid var(--border);
                border-radius: 0.375rem;
                font-family: inherit;
                font-size: 0.9rem;
                color: var(--text-main);
                background-color: #f8fafc;
                transition: all 0.2s;
            }}
            select:focus, input:focus {{ outline: none; border-color: var(--primary); background: white; }}
            
            .search-bar {{ width: 250px; }}

            /* DATA TABLE */
            .table-container {{ overflow-x: auto; border-radius: 0.75rem; border: 1px solid var(--border); box-shadow: var(--shadow); background: white; }}
            table {{ width: 100%; border-collapse: collapse; background: white; font-size: 0.9rem; }}
            th {{ background: #f8fafc; padding: 1rem; text-align: left; font-weight: 600; color: var(--text-muted); border-bottom: 1px solid var(--border); position: sticky; top: 0; z-index: 10; cursor: pointer; user-select: none; }}
            th:hover {{ background: #f1f5f9; color: var(--primary); }}
            td {{ padding: 1rem; border-bottom: 1px solid #f1f5f9; vertical-align: middle; color: var(--text-main); }}
            tr:hover {{ background: #f8fafc; }}
            
            .mol-img {{ width: 120px; height: 80px; object-fit: contain; background: white; padding: 0.25rem; border-radius: 0.5rem; border: 1px solid #f1f5f9; }}
            .score-val {{ font-weight: 700; color: var(--primary); font-size: 1.1em; }}
            .tag {{ padding: 0.25rem 0.75rem; border-radius: 9999px; font-size: 0.75rem; font-weight: 600; display: inline-block; }}
            .tag.in {{ background: #dcfce7; color: #166534; }}
            .tag.out {{ background: #fee2e2; color: #991b1b; }}
            
            /* COPY BUTTON */
            .copy-btn {{
                background: transparent;
                border: 1px solid var(--border);
                color: var(--text-muted);
                cursor: pointer;
                padding: 0.4rem 0.6rem;
                border-radius: 0.375rem;
                font-size: 0.8rem;
                margin-left: 0.5rem;
                transition: all 0.2s;
            }}
            .copy-btn:hover {{
                background: #f1f5f9;
                color: var(--primary);
                border-color: var(--primary-light);
            }}
            .copy-btn:active {{ transform: scale(0.95); }}
            
            .toast {{
                position: fixed;
                bottom: 2rem;
                right: 2rem;
                background: #1e293b;
                color: white;
                padding: 0.75rem 1.5rem;
                border-radius: 0.5rem;
                box-shadow: 0 10px 15px -3px rgba(0, 0, 0, 0.1);
                opacity: 0;
                transform: translateY(1rem);
                transition: all 0.3s;
                z-index: 100;
            }}
            .toast.show {{ opacity: 1; transform: translateY(0); }}
        </style>
    </head>
    <body>

        <!-- TOP NAVIGATION -->
        <div class="navbar">
            <div class="brand">
                <i class="fas fa-atom" style="color: var(--primary); font-size: 1.5rem;"></i>
                <div>
                    <h1>Predictions Data Grid</h1>
                    <div style="font-size: 0.75rem; color: var(--text-muted);">Virtual Screening Results</div>
                </div>
                <span class="badge" title="Total Molecules">{summary.get('total', 0)}</span>
            </div>
        </div>

        <div class="main">
            
            <!-- CONTROLS -->
            <div class="controls-bar">
                <div class="control-group">
                    <span class="control-label"><i class="fas fa-search"></i> Search</span>
                    <input type="text" id="searchInput" class="search-bar" placeholder="Name or SMILES..." onkeyup="applyFilters()">
                </div>
                
                <div class="control-group" style="border-left: 1px solid var(--border); padding-left: 1rem;">
                    <span class="control-label"><i class="fas fa-filter"></i> AD Status</span>
                    <select id="filterAD" onchange="applyFilters()">
                        <option value="ALL">All Compounds</option>
                        <option value="IN">In Domain</option>
                        <option value="OUT">Out of Domain</option>
                    </select>
                </div>
                
                <div class="control-group">
                    <span class="control-label"><i class="fas fa-sort"></i> Sort By</span>
                    <select id="sortBy" onchange="applyFilters()">
                        <option value="SCORE_ASC">Score (Best First)</option>
                        <option value="SCORE_DESC">Score (Worst First)</option>
                        <option value="MW_ASC">Mol. Weight (Low-High)</option>
                        <option value="MW_DESC">Mol. Weight (High-Low)</option>
                        <option value="LOGP_ASC">LogP (Low-High)</option>
                        <option value="LOGP_DESC">LogP (High-Low)</option>
                    </select>
                </div>
                
                <div style="margin-left: auto; font-size: 0.9rem; color: var(--text-muted);">
                    Showing <span id="visibleCount" style="font-weight: 700; color: var(--primary);">0</span> records
                </div>
            </div>
            
            <!-- TABLE -->
            <div class="table-container">
                <table id="dataTable">
                    <thead>
                        <tr>
                            <th style="width: 50px;">#</th>
                            <th style="width: 140px;">Structure</th>
                            <th>Name / ID</th>
                            <th>SMILES</th>
                            <th>Score</th>
                            <th>AD Status</th>
                            <th>MW</th>
                            <th>LogP</th>
                        </tr>
                    </thead>
                    <tbody></tbody>
                </table>
            </div>

        </div>

        <div id="toast" class="toast"><i class="fas fa-check-circle"></i> SMILES copied to clipboard!</div>

        <script>
            // DATA INJECTION
            const allData = {json.dumps(predictions[:1000])}; // Load up to 1000 records
            
            // RENDER FUNCTION
            function renderTable(data) {{
                const tbody = document.querySelector('#dataTable tbody');
                tbody.innerHTML = '';
                
                data.forEach((p, i) => {{
                    const tr = document.createElement('tr');
                    tr.innerHTML = `
                        <td>${{i+1}}</td>
                        <td><img src="${{p.image}}" class="mol-img" loading="lazy"></td>
                        <td style="font-weight: 500;">${{p.name}}</td>
                        <td style="font-family: monospace; font-size: 0.85rem; color: #475569;">
                            <div style="display: flex; align-items: center; max-width: 300px;">
                                <div style="overflow: hidden; text-overflow: ellipsis; white-space: nowrap; margin-right: 5px;">${{p.smiles}}</div>
                                <button class="copy-btn" onclick="copyToClipboard('${{p.smiles.replace(/'/g, "\\\\'") }}')" title="Copy SMILES">
                                    <i class="far fa-copy"></i>
                                </button>
                            </div>
                        </td>
                        <td class="score-val">${{p.score.toFixed(3)}}</td>
                        <td><span class="tag ${{p.ad === 'IN' ? 'in' : 'out'}}">${{p.ad}}</span></td>
                        <td>${{p.mw.toFixed(1)}}</td>
                        <td>${{p.logp.toFixed(2)}}</td>
                    `;
                    tbody.appendChild(tr);
                }});
                
                document.getElementById('visibleCount').textContent = data.length;
            }}

            // FILTER & SORT ENGINE
            function applyFilters() {{
                let result = [...allData];
                
                // 1. SEARCH
                const term = document.getElementById('searchInput').value.toUpperCase();
                if (term) {{
                    result = result.filter(p => 
                        (p.name && p.name.toUpperCase().includes(term)) || 
                        (p.smiles && p.smiles.toUpperCase().includes(term))
                    );
                }}
                
                // 2. FILTER AD
                const adFilter = document.getElementById('filterAD').value;
                if (adFilter !== 'ALL') {{
                    result = result.filter(p => p.ad === adFilter);
                }}
                
                // 3. SORT
                const sortType = document.getElementById('sortBy').value;
                result.sort((a, b) => {{
                    switch(sortType) {{
                        case 'SCORE_ASC': return a.score - b.score;
                        case 'SCORE_DESC': return b.score - a.score;
                        case 'MW_ASC':     return a.mw - b.mw;
                        case 'MW_DESC':    return b.mw - a.mw;
                        case 'LOGP_ASC':   return a.logp - b.logp;
                        case 'LOGP_DESC':  return b.logp - a.logp;
                        default: return a.score - b.score;
                    }}
                }});
                
                renderTable(result);
            }}
            
            // COPY FUNCTION
            function copyToClipboard(text) {{
                navigator.clipboard.writeText(text).then(() => {{
                    const toast = document.getElementById('toast');
                    toast.classList.add('show');
                    setTimeout(() => {{
                        toast.classList.remove('show');
                    }}, 2500);
                }});
            }}
            
            // INITIAL RENDER
            applyFilters();

        </script>
    </body>
    </html>
    """

    with open("outputs/report.html", "w") as f:
        f.write(html_content)
    
    print(f"Interactive Dashboard generated: outputs/report.html")

if __name__ == "__main__":
    generate_report()