"""HTML Generator for Node 08: Global Metabolic Load Analysis (Python 3)
Generates interactive HTML visualization matching Node 4's Card Grid style.
"""

import json
import os

def generate_metabolic_load_html(json_filename: str = 'metabolic_load_data.json') -> str:
    """Generate HTML for Global Metabolic Load visualization"""
    
    # Locate JSON
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(SCRIPT_DIR, 'outputs', json_filename)
    
    if not os.path.exists(json_path):
        # Try looking in current working directory as fallback
        if os.path.exists(os.path.join(os.getcwd(), 'outputs', json_filename)):
            json_path = os.path.join(os.getcwd(), 'outputs', json_filename)
        elif os.path.exists(json_filename):
             json_path = json_filename

    json_content = '{}'
    error_msg = ''
    
    try:
        print(f"DEBUG: Reading JSON from {json_path}")
        with open(json_path, 'r', encoding='utf-8') as f:
            json_content = f.read()
            # Verify valid JSON
            json.loads(json_content)
    except Exception as e:
        print(f"Warning: Could not read JSON: {e}")
        error_msg = str(e)
        json_content = '{}'

    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Global Metabolic Load Analysis</title>
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
        
        /* Header */
        header {{ background: var(--card-bg); padding: 15px 30px; border-bottom: 1px solid #ddd; display: flex; justify-content: space-between; align-items: center; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }}
        h1 {{ margin: 0; font-size: 1.2rem; color: var(--primary); }}
        .timestamp {{ font-size: 0.8rem; color: #777; }}

        /* Controls Bar */
        .controls {{ background: #fff; padding: 10px 30px; border-bottom: 1px solid #eee; display: flex; gap: 20px; align-items: center; flex-wrap: wrap; }}
        
        .toggle-group {{ display: flex; gap: 5px; background: #eee; padding: 4px; border-radius: 6px; }}
        .toggle-btn {{ border: none; background: none; padding: 8px 16px; cursor: pointer; border-radius: 4px; font-weight: 500; color: #666; transition: 0.2s; }}
        .toggle-btn.active {{ background: white; color: var(--accent); box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}

        .toggle-container {{ display: flex; align-items: center; gap: 8px; font-size: 0.9rem; margin-left: auto; }}
        input[type="checkbox"] {{ accent-color: var(--accent); cursor: pointer; transform: scale(1.2); }}

        .legend-bar {{ display: flex; gap: 15px; font-size: 0.85rem; margin-left: 20px; }}
        .legend-item {{ display: flex; align-items: center; gap: 5px; }}
        .dot {{ width: 10px; height: 10px; border-radius: 50%; display: inline-block; }}

        /* Main Content */
        #app {{ flex: 1; overflow-y: auto; padding: 20px 30px; display: flex; justify-content: center; }}
        
        .hero-container {{
             width: 100%;
             max-width: 1000px;
             margin-top: 20px;
        }}

        .chart-card {{
            background: var(--card-bg);
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.08);
            border: 1px solid #eee;
            height: 600px;
            display: flex;
            flex-direction: column;
        }}
        
        .card-header {{ display: flex; justify-content: space-between; padding-bottom: 10px; margin-bottom: 10px; border-bottom: 1px solid #f0f0f0; font-size: 1.1rem; font-weight: bold; color: var(--primary); }}
        
        .stats-badge {{ font-size: 0.9rem; background: #e8f6fd; color: var(--accent); padding: 4px 10px; border-radius: 4px; font-weight: 600; }}
        .sig-badge {{ background: #fee; color: #c0392b; }} /* Red for significant */

        .plot-div {{ flex: 1; width: 100%; }}

        .loading {{ text-align: center; padding: 50px; font-size: 1.2rem; color: #666; }}
        
        .error-msg {{
            background: #ffebee; color: #c62828; padding: 20px; border-radius: 8px; border: 1px solid #ffcdd2;
            margin: 20px; text-align: center;
        }}
        
        .info-msg {{
            background: #e3f2fd; color: #0d47a1; padding: 20px; border-radius: 8px; border: 1px solid #bbdefb;
            margin: 20px; text-align: center;
        }}
    </style>
    <script>
        const RAW_DATA = __JSON_DATA__;
    </script>
</head>
<body>
    <header>
        <h1>Global Metabolic Load Analysis</h1>
        <span class="timestamp" id="timestamp">Generated: Just now</span>
    </header>

    <div class="controls">
        <label style="font-weight:600; color:#555; margin-right:10px;">Chart Type:</label>
        <div class="toggle-group">
            <button class="toggle-btn active" onclick="setChartType('box')">Box Plot</button>
            <button class="toggle-btn" onclick="setChartType('violin')">Violin Plot</button>
        </div>

        <div class="legend-bar">
            <div class="legend-item"><span class="dot" style="background:#1f77b4"></span> PPMS</div>
            <div class="legend-item"><span class="dot" style="background:#ff7f0e"></span> SPMS</div>
            <div class="legend-item"><span class="dot" style="background:#2ca02c"></span> RRMS</div>
            <div class="legend-item"><span class="dot" style="background:#d62728"></span> GMG</div>
        </div>
        
        <div class="toggle-container">
            <input type="checkbox" id="sigFilter" onchange="renderCurrentView()">
            <label for="sigFilter" title="Show only results with p < 0.05">Significant Only (p<0.05)</label>
        </div>

        <div class="toggle-container">
            <input type="checkbox" id="showPoints" checked onchange="renderCurrentView()">
            <label for="showPoints">Show Data Points</label>
        </div>
    </div>

    <div id="app">
        <div class="hero-container" id="hero">
            <div class="loading">Initializing Dashboard...</div>
        </div>
    </div>

    <script>
        let currentType = 'box';
        
        const COLORS = {{
            'PPMS': '#1f77b4', 
            'SPMS': '#ff7f0e', 
            'RRMS': '#2ca02c', 
            'GMG': '#d62728'
        }};
        
        // Order for x-axis
        const GROUPS = ['PPMS', 'SPMS', 'RRMS', 'GMG'];

        function init() {{
            document.getElementById('timestamp').textContent = new Date().toLocaleString();
            
            // Check for empty data
            if (!RAW_DATA || Object.keys(RAW_DATA).length === 0) {{
                document.getElementById('app').innerHTML = `
                    <div class="error-msg">
                        <h3>Dashboard Error</h3>
                        <p>No data available to display.</p>
                        <p>Debug info: JSON could not be loaded. {error_msg}</p>
                    </div>
                `;
                return;
            }}
            
            renderCurrentView();
        }}
        
        function setChartType(type) {{
            currentType = type;
            document.querySelectorAll('.toggle-btn').forEach(b => b.classList.remove('active'));
            if(type === 'box') document.querySelector('.toggle-btn:first-child').classList.add('active');
            else document.querySelector('.toggle-btn:last-child').classList.add('active');
            renderCurrentView();
        }}

        function renderCurrentView() {{
            const container = document.getElementById('hero');
            
            if (!RAW_DATA['fig8'] || !RAW_DATA['fig8'].traces) {{
                  container.innerHTML = '<div class="error-msg">Missing Figure 8 Data</div>';
                  return;
            }}
            
            const data = RAW_DATA['fig8'].traces[0];
            const pVal = data.stats.p_value;
            const testName = data.stats.test || 'Test';
            const showPoints = document.getElementById('showPoints').checked;
            const isSigOnly = document.getElementById('sigFilter').checked;

            container.innerHTML = '';
            
            // Filter logic
            if (isSigOnly && pVal >= 0.05) {{
                container.innerHTML = `
                    <div class="info-msg">
                        <h3>Result Hidden</h3>
                        <p>The result is not statistically significant (p = ${{pVal.toExponential(2)}}).</p>
                        <p>Uncheck "Significant Only" to view.</p>
                    </div>
                `;
                return;
            }}
            
            const card = document.createElement('div');
            card.className = 'chart-card';
            
            // Add sig-badge class if significant
            const badgeClass = pVal < 0.05 ? 'stats-badge sig-badge' : 'stats-badge';
            
            card.innerHTML = `
                <div class="card-header">
                    <span>Total Amino Acid Concentration by Disease Type</span>
                    <span class="${{badgeClass}}">${{testName}} P = ${{pVal.toExponential(2)}}</span>
                </div>
                <div class="plot-div" id="plot_hero"></div>
            `;
            container.appendChild(card);
            
            const plotData = [];
            
            GROUPS.forEach(g => {{
                if(data[g] && data[g].y) {{
                    plotData.push({{
                        y: data[g].y,
                        type: currentType,
                        name: g,
                        marker: {{ color: COLORS[g] }},
                        boxpoints: showPoints ? 'all' : false,
                        points: showPoints ? 'all' : false, // for violin
                        jitter: 0.3,
                        pointpos: 0,
                        box: {{ visible: true }}, // inside violin
                        meanline: {{ visible: true }}
                    }});
                }}
            }});
            
            const layout = {{
                margin: {{ t: 20, r: 20, b: 40, l: 60 }},
                yaxis: {{ title: 'Concentration [nmol/ml]', showgrid: true, gridcolor: '#eee', zeroline: false }},
                xaxis: {{ showgrid: false }},
                paper_bgcolor: 'rgba(0,0,0,0)',
                plot_bgcolor: 'rgba(0,0,0,0)',
                showlegend: false,
                font: {{ family: 'Segoe UI, sans-serif' }}
            }};
            
            Plotly.newPlot('plot_hero', plotData, layout, {{responsive: true, displayModeBar: false}});
        }}

        // Run
        setTimeout(init, 100);
    </script>
</body>
</html>"""
    
    # Inject JSON data
    html = html.replace('__JSON_DATA__', json_content)
    
    return html
