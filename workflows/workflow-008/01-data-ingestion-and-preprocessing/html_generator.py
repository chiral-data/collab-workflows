"""HTML Generator for Node 1: Preprocessing Report
Visualizes Raw vs Standardized Amino Acid distributions + Data Tables.
"""
import json
import os

def generate_preprocessing_report(json_path):
    """Generate HTML report comparing Raw vs Processed data"""
    
    # Load Data
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
    except Exception as e:
        return f"<h1>Error loading JSON: {e}</h1>"

    # Embed JSON directly into HTML for the frontend to consume
    json_str = json.dumps(data)

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Data Preprocessing Report (Node 1)</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {{
            --primary: #2c3e50;
            --accent: #e74c3c;
            --processed: #3498db;
            --bg: #f8f9fa;
        }}
        
        body {{ font-family: 'Segoe UI', sans-serif; margin: 0; padding: 20px; background: var(--bg); }}
        
        .container {{ max-width: 1400px; margin: 0 auto; background: white; padding: 30px; border-radius: 12px; box-shadow: 0 4px 12px rgba(0,0,0,0.05); }}
        
        header {{ border-bottom: 2px solid #eee; padding-bottom: 20px; margin-bottom: 30px; display: flex; justify-content: space-between; align-items: center; }}
        h1 {{ margin: 0; color: var(--primary); }}
        .stats {{ color: #7f8c8d; font-size: 0.9em; }}
        
        .section-title {{ font-size: 1.2rem; font-weight: bold; margin: 30px 0 15px 0; color: var(--primary); border-left: 5px solid var(--accent); padding-left: 10px; }}
        
        /* PLOTS */
        .plot-container {{ display: flex; gap: 20px; height: 600px; margin-bottom: 30px; }}
        .plot-box {{ flex: 1; border: 1px solid #eee; border-radius: 8px; padding: 10px; display: flex; flex-direction: column; }}
        .plot-header {{ text-align: center; font-weight: bold; font-size: 1.1em; margin-bottom: 10px; padding: 10px; background: #f8f9fa; border-radius: 4px; }}
        .raw-header {{ color: var(--accent); border-left: 4px solid var(--accent); }}
        .proc-header {{ color: var(--processed); border-left: 4px solid var(--processed); }}
        .chart {{ flex: 1; min-height: 0; }}
        
        /* TABS */
        .tabs {{ display: flex; gap: 5px; margin-bottom: 10px; }}
        .tab-btn {{ padding: 10px 20px; border: none; background: #eee; cursor: pointer; border-radius: 5px 5px 0 0; font-weight: 600; color: #555; }}
        .tab-btn.active {{ background: var(--primary); color: white; }}
        
        /* TABLES */
        .table-container {{ overflow-x: auto; height: 500px; border: 1px solid #ddd; border-radius: 8px; font-size: 0.85em; }}
        table {{ width: 100%; border-collapse: collapse; white-space: nowrap; }}
        th, td {{ padding: 8px 12px; border-bottom: 1px solid #eee; text-align: left; }}
        th {{ background: #f1f2f6; position: sticky; top: 0; z-index: 10; font-weight: 600; color: #444; cursor: pointer; user-select: none; }}
        th:hover {{ background: #e1e2e6; color: var(--primary); }}
        th::after {{ content: ' ↕'; opacity: 0.3; font-size: 0.8em; }}
        tr:hover {{ background: #f8f9fa; }}
        .numeric {{ text-align: right; font-family: monospace; }}
        
        .controls-bar {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px; }}
        .search-box {{ padding: 8px 12px; border: 1px solid #ccc; border-radius: 4px; width: 250px; }}
        
        .explanation {{ background: #e8f4fd; padding: 15px; border-radius: 6px; border-left: 4px solid #3498db; margin-top: 20px; }}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <div>
                <h1>📊 Preprocessing Report: Raw vs Standardized</h1>
                <p>Node 1 Data Ingestion Pipeline</p>
            </div>
            <div class="stats">
                <strong>Total Records:</strong> {data['metadata']['total_records']}<br>
                <strong>Features:</strong> {len(data['metadata']['amino_acids'])} Amino Acids
            </div>
        </header>

        <div class="explanation">
            <strong>Why Standardization Matters:</strong> 
            Raw biological data often varies by orders of magnitude (e.g., GLN vs ORN). 
            Standardization (Z-score normalization) centers all features at 0 with a variance of 1, 
            ensuring that high-concentration metabolites do not dominate downstream analyses.
        </div>
        
        <!-- SECTION 1: VISUALIZATION -->
        <div class="section-title">1. Distribution Overview</div>
        <div class="plot-container">
            <div class="plot-box">
                <div class="plot-header raw-header">Before: Raw Concentrations (µM)</div>
                <div id="raw_plot" class="chart"></div>
            </div>
            <div class="plot-box">
                <div class="plot-header proc-header">After: Standardized (Z-Scores)</div>
                <div id="proc_plot" class="chart"></div>
            </div>
        </div>

        <!-- SECTION 2: DATA INSPECTION -->
        <div class="section-title">2. Data Inspection (Table View)</div>
        <div class="controls-bar">
            <div class="tabs">
                <button class="tab-btn active" onclick="showTable('raw')">Raw Data</button>
                <button class="tab-btn" onclick="showTable('processed')">Processed Data</button>
            </div>
            <input type="text" id="searchInput" class="search-box" placeholder="Search data..." onkeyup="filterTable()">
        </div>
        
        <div id="table_raw" class="table-container active-table"></div>
        <div id="table_processed" class="table-container" style="display:none;"></div>

    </div>

    <script>
        const DATA = {json_str};
        let currentTableId = 'table_raw'; // Track active table for filtering
        
        function init() {{
            renderPlots();
            renderTable('raw', DATA.tables.raw);
            renderTable('processed', DATA.tables.processed);
        }}
        
        function showTable(type) {{
            // Toggle Buttons
            document.querySelectorAll('.tab-btn').forEach(btn => btn.classList.remove('active'));
            event.target.classList.add('active');
            
            // Toggle Tables
            document.getElementById('table_raw').style.display = 'none';
            document.getElementById('table_processed').style.display = 'none';
            document.getElementById('table_raw').classList.remove('active-table');
            document.getElementById('table_processed').classList.remove('active-table');
            
            // Activate selected
            currentTableId = 'table_' + type;
            document.getElementById(currentTableId).style.display = 'block';
            document.getElementById(currentTableId).classList.add('active-table');
            
            // Re-apply filter
            filterTable();
        }}
        
        function renderTable(id, records) {{
            if(!records || records.length === 0) return;
            
            const cols = Object.keys(records[0]);
            let html = `<table id="tbl_${{id}}"><thead><tr>`;
            
            // Headers with Sort OnClick
            cols.forEach((c, idx) => {{
                html += `<th onclick="sortTable('${{id}}', ${{idx}})">${{c}}</th>`;
            }});
            html += '</tr></thead><tbody>';
            
            // Rows
            records.forEach(row => {{
                html += '<tr>';
                cols.forEach(c => {{
                    let val = row[c];
                    // Format numeric values
                    if (typeof val === 'number') {{
                        val = val.toFixed(4);
                        html += `<td class="numeric">${{val}}</td>`;
                    }} else {{
                        html += `<td>${{val}}</td>`;
                    }}
                }});
                html += '</tr>';
            }});
            
            html += '</tbody></table>';
            document.getElementById('table_' + id).innerHTML = html;
        }}

        // --- FILTERING ---
        function filterTable() {{
            const input = document.getElementById('searchInput');
            const filter = input.value.toUpperCase();
            if(!currentTableId) return;
            
            const table = document.getElementById(currentTableId).querySelector('table');
            if(!table) return;
            
            const tr = table.getElementsByTagName('tr');

            for (let i = 1; i < tr.length; i++) {{ // Skip header
                let visible = false;
                const tds = tr[i].getElementsByTagName('td');
                for (let j = 0; j < tds.length; j++) {{
                    if (tds[j]) {{
                        const txtValue = tds[j].textContent || tds[j].innerText;
                        if (txtValue.toUpperCase().indexOf(filter) > -1) {{
                            visible = true;
                            break;
                        }}
                    }}
                }}
                tr[i].style.display = visible ? "" : "none";
            }}
        }}

        // --- SORTING ---
        function sortTable(type, colIndex) {{
            const table = document.getElementById('table_' + type).querySelector('table');
            let switching = true;
            let shouldSwitch, i;
            let dir = "asc"; 
            let switchcount = 0;
            
            while (switching) {{
                switching = false;
                const rows = table.rows;
                
                // Loop through all table rows (except header)
                for (i = 1; i < (rows.length - 1); i++) {{
                    shouldSwitch = false;
                    const x = rows[i].getElementsByTagName("TD")[colIndex];
                    const y = rows[i + 1].getElementsByTagName("TD")[colIndex];
                    
                    // Check if numeric or string
                    let xContent = x.innerHTML.toLowerCase();
                    let yContent = y.innerHTML.toLowerCase();
                    const isNum = !isNaN(parseFloat(xContent)) && !isNaN(parseFloat(yContent));
                    
                    if(isNum) {{
                        xContent = parseFloat(xContent);
                        yContent = parseFloat(yContent);
                    }}

                    if (dir == "asc") {{
                        if (xContent > yContent) {{ shouldSwitch = true; break; }}
                    }} else if (dir == "desc") {{
                        if (xContent < yContent) {{ shouldSwitch = true; break; }}
                    }}
                }}
                
                if (shouldSwitch) {{
                    rows[i].parentNode.insertBefore(rows[i + 1], rows[i]);
                    switching = true;
                    switchcount++;
                }} else {{
                    if (switchcount == 0 && dir == "asc") {{
                        dir = "desc";
                        switching = true;
                    }}
                }}
            }}
        }}

        function renderPlots() {{
            const amino_acids = DATA.metadata.amino_acids;
            
            // 1. Raw Traces
            const rawTraces = [];
            amino_acids.forEach(aa => {{
                if (DATA.raw[aa]) {{
                    rawTraces.push({{
                        y: DATA.raw[aa],
                        type: 'box',
                        name: aa,
                        boxpoints: 'outliers',
                        marker: {{ color: '#e74c3c', size: 3 }}
                    }});
                }}
            }});

            Plotly.newPlot('raw_plot', rawTraces, {{
                margin: {{ t: 20, r: 20, b: 100, l: 60 }},
                yaxis: {{ title: 'Concentration (Log)', gridcolor: '#eee', type: 'log' }},
                xaxis: {{ tickangle: -90 }},
                showlegend: false
            }});

            // 2. Processed Traces
            const procTraces = [];
            amino_acids.forEach(aa => {{
                if (DATA.processed[aa]) {{
                    procTraces.push({{
                        y: DATA.processed[aa],
                        type: 'box',
                        name: aa,
                        boxpoints: 'outliers',
                        marker: {{ color: '#3498db', size: 3 }}
                    }});
                }}
            }});
            
            Plotly.newPlot('proc_plot', procTraces, {{
                margin: {{ t: 20, r: 20, b: 100, l: 60 }},
                yaxis: {{ title: 'Z-Score', gridcolor: '#eee' }},
                xaxis: {{ tickangle: -90 }},
                showlegend: false,
                shapes: [{{ type: 'line', x0: 0, x1: 1, xref: 'paper', y0: 0, y1: 0, line: {{ color: '#333', dash: 'dot', width: 1 }} }}]
            }});
        }}

        init();
    </script>
</body>
</html>"""

    return html
