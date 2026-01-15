"""HTML Generator for Node 02: Cohort Demographics (Python 3)
Generates interactive HTML Dashboard matching Node 10 style.
"""

import json
import os

def generate_demographics_html(json_filename: str = 'demographics_data.json') -> str:
    """Generate HTML for Demographics Dashboard"""
    
    # Locate JSON
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(SCRIPT_DIR, 'outputs', json_filename)
    
    if not os.path.exists(json_path):
        if os.path.exists(os.path.join(os.getcwd(), 'outputs', json_filename)):
            json_path = os.path.join(os.getcwd(), 'outputs', json_filename)
        elif os.path.exists(json_filename):
             json_path = json_filename

    json_content = '{}'

    try:
        print(f"DEBUG: Reading JSON from {json_path}")
        with open(json_path, 'r', encoding='utf-8') as f:
            json_content = f.read()
            json.loads(json_content) # Verify
    except Exception as e:
        print(f"Warning: Could not read JSON: {e}")
        json_content = '{}'

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Cohort Demographics</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {{
            --primary: #2c3e50;
            --accent: #3498db;
            --bg: #f8f9fa;
            --card-bg: #ffffff;
            --text: #333;
        }}
        
        body {{ font-family: 'Segoe UI', Roboto, sans-serif; margin: 0; padding: 0; background: var(--bg); color: var(--text); height: 100vh; display: flex; flex-direction: column; }}
        
        header {{ background: var(--card-bg); padding: 15px 30px; border-bottom: 1px solid #ddd; display: flex; justify-content: space-between; align-items: center; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }}
        h1 {{ margin: 0; font-size: 1.2rem; color: var(--primary); }}
        .timestamp {{ font-size: 0.8rem; color: #777; }}

        .controls {{ background: #fff; padding: 10px 30px; border-bottom: 1px solid #eee; display: flex; gap: 20px; align-items: center; flex-wrap: wrap; }}
        
        .tab-group {{ display: flex; gap: 5px; background: #eee; padding: 4px; border-radius: 6px; }}
        .tab-btn {{ border: none; background: none; padding: 8px 16px; cursor: pointer; border-radius: 4px; font-weight: 500; color: #666; transition: 0.2s; }}
        .tab-btn.active {{ background: white; color: var(--accent); box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}

        .toggle-group {{ display: flex; gap: 5px; background: #eee; padding: 4px; border-radius: 6px; margin-left: 0; }}
        .toggle-btn {{ border: none; background: none; padding: 8px 16px; cursor: pointer; border-radius: 4px; font-weight: 500; color: #666; transition: 0.2s; }}
        .toggle-btn.active {{ background: white; color: var(--accent); box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}

        #app {{ flex: 1; overflow-y: auto; padding: 20px 30px; }}
        
        /* Layouts */
        .dashboard-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .full-width {{ grid-column: 1 / -1; }}
        
        .chart-card {{
            background: var(--card-bg);
            border-radius: 8px;
            padding: 15px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            display: flex;
            flex-direction: column;
            height: 400px;
        }}
        
        .table-card {{
            background: var(--card-bg);
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            overflow-x: auto;
            height: auto;
        }}
        
        .card-title {{ font-weight: bold; color: var(--primary); margin-bottom: 15px; border-bottom: 1px solid #eee; padding-bottom: 10px; }}
        
        /* Table Styles */
        table {{ width: 100%; border-collapse: collapse; font-size: 0.9em; }}
        th {{ background: #f8f9fa; padding: 12px; text-align: left; font-weight: 600; color: #555; border-bottom: 2px solid #ddd; }}
        td {{ padding: 12px; border-bottom: 1px solid #eee; color: #333; }}
        tr:hover {{ background: #fcfcfc; }}
        
        .loading {{ text-align: center; padding: 50px; font-size: 1.2rem; color: #666; }}
    </style>
    <script>
        const RAW_DATA = __JSON_DATA__;
    </script>
</head>
<body>
    <header>
        <h1>Cohort Demographics</h1>
        <span class="timestamp" id="timestamp">Generated: Just now</span>
    </header>

    <div class="controls">
        <label style="font-weight:600; color:#555;">Dashboard:</label>
        <div class="tab-group">
            <button class="tab-btn active" onclick="switchTab('tab1')">MS vs Controls</button>
            <button class="tab-btn" onclick="switchTab('tab2')">MS vs MG</button>
        </div>

        <label style="font-weight:600; color:#555; margin-left:15px;">Chart:</label>
        <div class="toggle-group">
             <button class="toggle-btn active" onclick="setChartType('box')">Box</button>
             <button class="toggle-btn" onclick="setChartType('violin')">Violin</button>
        </div>
    </div>

    <div id="app">
        <div class="loading">Initializing Dashboard...</div>
    </div>

    <script>
        let currentTab = 'tab1';
        let currentType = 'box';
        
        const COLORS = {{
            'MS': '#4682b4',     // SteelBlue
            'Control': '#fa8072', // Salmon
            'MG': '#ffa500'      // Orange
        }};

        function init() {{
            document.getElementById('timestamp').textContent = new Date().toLocaleString();
            if (!RAW_DATA || Object.keys(RAW_DATA).length === 0) {{
                document.getElementById('app').innerHTML = '<div style="padding:20px">No Data Available</div>';
                return;
            }}
            switchTab('tab1');
        }}
        
        function switchTab(tab) {{
            currentTab = tab;
            document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
            if(tab === 'tab1') document.querySelector('.tab-btn:first-child').classList.add('active');
            else document.querySelector('.tab-btn:last-child').classList.add('active');
            renderDashboard();
        }}
        
        function setChartType(type) {{
            currentType = type;
            document.querySelectorAll('.toggle-btn').forEach(b => b.classList.remove('active'));
            if(type === 'box') document.querySelector('.toggle-btn:first-child').classList.add('active');
            else document.querySelector('.toggle-btn:last-child').classList.add('active');
            renderDashboard();
        }}

        function renderDashboard() {{
            const app = document.getElementById('app');
            const plots = RAW_DATA.plots;
            
            let html = '<div class="dashboard-grid">';
            
            if (currentTab === 'tab1') {{
                // --- Tab 1 Visuals ---
                // 1. Age Distribution (Box)
                html += `
                    <div class="chart-card">
                        <div class="card-title">Age Distribution</div>
                        <div id="plot_age" style="flex:1"></div>
                    </div>
                `;
                // 2. Sex Ratio (Bar)
                html += `
                    <div class="chart-card">
                        <div class="card-title">Sex Distribution</div>
                        <div id="plot_sex" style="flex:1"></div>
                    </div>
                `;
                html += '</div>'; // End Grid
                
                // 3. Table 1 (Full Width)
                html += `
                    <div class="table-card">
                        <div class="card-title">Detailed Demographics (Table 1)</div>
                        ${{generateTableHTML(RAW_DATA.table1)}}
                    </div>
                `;
                
                app.innerHTML = html;
                
                // Render Plots
                plotDist('plot_age', plots.age, ['MS', 'Control']);
                plotSex('plot_sex', plots.sex, ['MS', 'Control']);

            }} else {{
                // --- Tab 2 Visuals ---
                // 1. Age Distribution
                html += `
                    <div class="chart-card">
                        <div class="card-title">Age Comparison</div>
                        <div id="plot_age_2" style="flex:1"></div>
                    </div>
                `;
                // 2. Duration Distribution
                html += `
                    <div class="chart-card">
                        <div class="card-title">Disease Duration</div>
                        <div id="plot_dur_2" style="flex:1"></div>
                    </div>
                `;
                // 3. Subtype Counts
                html += `
                    <div class="chart-card full-width" style="height:350px">
                        <div class="card-title">Subtype Composition</div>
                        <div id="plot_sub_2" style="flex:1"></div>
                    </div>
                `;
                html += '</div>'; // End Grid
                
                // 4. Table 2
                html += `
                    <div class="table-card">
                        <div class="card-title">MS vs MG Statistics (Table 2)</div>
                        ${{generateTableHTML(RAW_DATA.table2)}}
                    </div>
                `;
                
                app.innerHTML = html;
                
                plotDist('plot_age_2', plots.age, ['MS', 'MG']);
                plotDist('plot_dur_2', plots.duration, ['MS', 'MG']);
                plotSubtypes('plot_sub_2', plots.subtypes);
            }}
        }}

        // --- Plotting Helpers ---
        function plotDist(divId, dataObj, groups) {{
            const traces = [];
            groups.forEach(g => {{
                if (dataObj[g]) {{
                    traces.push({{
                        y: dataObj[g],
                        type: currentType,
                        name: g,
                        marker: {{ color: COLORS[g] }},
                        boxpoints: 'all',
                        jitter: 0.3,
                        pointpos: -1.8
                    }});
                }}
            }});
            const layout = {{ margin: {{t:10, b:30, l:30, r:10}}, showlegend: true }};
            Plotly.newPlot(divId, traces, layout, {{displayModeBar: false}});
        }}

        function plotSex(divId, sexData, groups) {{
            const males = [];
            const females = [];
            groups.forEach(g => {{
                if(sexData[g]) {{
                    males.push(sexData[g].Male);
                    females.push(sexData[g].Female);
                }}
            }});
            
            const traces = [
                {{ x: groups, y: males, name: 'Male', type: 'bar', marker: {{color: '#34495e'}} }},
                {{ x: groups, y: females, name: 'Female', type: 'bar', marker: {{color: '#e74c3c'}} }}
            ];
            const layout = {{ barmode: 'group', margin: {{t:10, b:30, l:30, r:10}} }};
            Plotly.newPlot(divId, traces, layout, {{displayModeBar: false}});
        }}
        
        function plotSubtypes(divId, subData) {{
            // Nested Pies
             const labelsMS = Object.keys(subData.MS);
             const valuesMS = Object.values(subData.MS);
             const labelsMG = Object.keys(subData.MG);
             const valuesMG = Object.values(subData.MG);
             
             const trace1 = {{
                 values: valuesMS, labels: labelsMS,
                 type: 'pie', name: 'MS Subtypes', title: 'MS',
                 domain: {{ column: 0 }},
                 hole: 0.4
             }};
             const trace2 = {{
                 values: valuesMG, labels: labelsMG,
                 type: 'pie', name: 'MG Subtypes', title: 'MG',
                 domain: {{ column: 1 }},
                 hole: 0.4
             }};
             
             const layout = {{ 
                 grid: {{rows: 1, columns: 2}},
                 margin: {{t:0, b:10, l:0, r:0}},
                 showlegend: true
             }};
             Plotly.newPlot(divId, [trace1, trace2], layout, {{displayModeBar: false}});
        }}

        function generateTableHTML(tableData) {{
            if(!tableData) return '';
            let h = '<table><thead><tr>';
            tableData.columns.forEach(c => h += `<th>${{c}}</th>`);
            h += '</tr></thead><tbody>';
            tableData.rows.forEach(r => {{
                h += '<tr>';
                r.forEach(c => h += `<td>${{c}}</td>`);
                h += '</tr>';
            }});
            h += '</tbody></table>';
            return h;
        }}

        setTimeout(init, 100);
    </script>
</body>
</html>"""
    
    # Inject JSON data
    html = html.replace('__JSON_DATA__', json_content)
    
    return html
