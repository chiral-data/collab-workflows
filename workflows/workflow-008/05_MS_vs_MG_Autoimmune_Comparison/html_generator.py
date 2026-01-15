"""HTML Generator for Node 05: Autoimmune Comparison Dashboard (Python 3)
Generates interactive HTML visualization matching Node 10 style (Grid Only).
"""

import json
import os

def generate_autoimmune_html(json_filename: str = 'autoimmune_data.json') -> str:
    """Generate HTML for Autoimmune Comparison visualization"""
    
    # Locate JSON
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(SCRIPT_DIR, 'outputs', json_filename)
    
    if not os.path.exists(json_path):
        # Checks for local execution fallbacks
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
    <title>Autoimmune Comparison: MS+MG vs Controls</title>
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
        
        .toggle-group {{ display: flex; gap: 5px; background: #eee; padding: 4px; border-radius: 6px; margin-left: 0; }}
        .toggle-btn {{ border: none; background: none; padding: 8px 16px; cursor: pointer; border-radius: 4px; font-weight: 500; color: #666; transition: 0.2s; }}
        .toggle-btn.active {{ background: white; color: var(--accent); box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}

        .toggle-container {{ display: flex; align-items: center; gap: 8px; font-size: 0.9rem; margin-left: auto; }}
        input[type="checkbox"] {{ accent-color: var(--accent); cursor: pointer; transform: scale(1.2); }}

        .legend-bar {{ display: flex; gap: 15px; font-size: 0.85rem; margin-left: 20px; }}
        .legend-item {{ display: flex; align-items: center; gap: 5px; }}
        .dot {{ width: 10px; height: 10px; border-radius: 50%; display: inline-block; }}

        /* Main Content Grid */
        #app {{ flex: 1; overflow-y: auto; padding: 20px 30px; }}
        
        .grid-container {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
            gap: 20px;
            padding-bottom: 50px;
        }}

        .chart-card {{
            background: var(--card-bg);
            border-radius: 8px;
            padding: 10px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            transition: transform 0.2s;
            border: 1px solid transparent;
            height: 350px;
            display: flex;
            flex-direction: column;
        }}
        
        .chart-card:hover {{ transform: translateY(-2px); box-shadow: 0 5px 15px rgba(0,0,0,0.1); }}
        
        .card-header {{ display: flex; justify-content: space-between; padding: 0 10px; margin-bottom: 5px; font-size: 0.95rem; font-weight: bold; color: var(--primary); }}
        .stats-badge {{ font-size: 0.75rem; background: #eee; padding: 2px 6px; border-radius: 4px; font-weight: normal; color: #666; }}
        .sig-badge {{ background: #fee; color: #c0392b; font-weight: bold; }}
        
        .plot-div {{ flex: 1; width: 100%; }}

        .loading {{ text-align: center; padding: 50px; font-size: 1.2rem; color: #666; }}
        
        .error-msg {{
            background: #ffebee; color: #c62828; padding: 20px; border-radius: 8px; border: 1px solid #ffcdd2;
            margin: 20px; text-align: center;
        }}
    </style>
    <script>
        const RAW_DATA = __JSON_DATA__;
    </script>
</head>
<body>
    <header>
        <h1>Autoimmune Comparison: MS+MG vs Controls</h1>
        <span class="timestamp" id="timestamp">Generated: Just now</span>
    </header>

    <div class="controls">
        <label style="font-weight:600; color:#555; margin-right:5px;">Chart:</label>
        <div class="toggle-group">
             <button class="toggle-btn active" onclick="setChartType('box')">Box</button>
             <button class="toggle-btn" onclick="setChartType('violin')">Violin</button>
        </div>

        <div class="legend-bar">
            <div class="legend-item"><span class="dot" style="background:#87ceeb"></span> MS+MG</div>
            <div class="legend-item"><span class="dot" style="background:#808080"></span> Controls</div>
        </div>
        
        <div class="toggle-container">
            <input type="checkbox" id="sigFilter" onchange="renderCurrentView()">
            <label for="sigFilter" title="Show only results with p < 0.05">Significant Only (p<0.05)</label>
        </div>
        
        <div class="toggle-container">
            <input type="checkbox" id="showPoints" checked onchange="renderCurrentView()">
            <label for="showPoints">Show Data</label>
        </div>
    </div>

    <div id="app">
        <div class="loading">Initializing Dashboard...</div>
    </div>

    <script>
        let currentType = 'box';
        
        const COLORS = {{
            'MS+MG': '#87ceeb', // SkyBlue
            'Controls': '#808080' // Grey
        }};

        function init() {{
            document.getElementById('timestamp').textContent = new Date().toLocaleString();
            
             // Check for empty data
            if (!RAW_DATA || Object.keys(RAW_DATA).length === 0) {{
                document.getElementById('app').innerHTML = `
                    <div class="error-msg">
                        <h3>Dashboard Error</h3>
                        <p>No data available to display.</p>
                        <p>Debug info: JSON could not be loaded.</p>
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
            const container = document.getElementById('app');
            
            // Default to 'grid' view, ignoring 'hero'
            if (!RAW_DATA['grid'] || !RAW_DATA['grid'].traces) {{
                  container.innerHTML = '<div class="error-msg">Missing Grid Data</div>';
                  return;
            }}

            const traces = RAW_DATA['grid'].traces;
            const isSigOnly = document.getElementById('sigFilter').checked;
            const showPoints = document.getElementById('showPoints').checked;
            
            // Clean container
            container.innerHTML = '<div class="grid-container" id="grid"></div>';
            const grid = document.getElementById('grid');
            
            const groups = ['MS+MG', 'Controls'];

            traces.forEach((item, index) => {{
                const pVal = item.stats ? item.stats.p_value : 1.0;
                
                if (isSigOnly && pVal >= 0.05) return;
                
                const card = document.createElement('div');
                card.className = 'chart-card';
                const badgeClass = pVal < 0.05 ? 'stats-badge sig-badge' : 'stats-badge';
                
                card.innerHTML = `
                    <div class="card-header">
                        <span>${{item.aa}}</span>
                        <span class="${{badgeClass}}">P = ${{pVal.toExponential(2)}}</span>
                    </div>
                    <div class="plot-div" id="plot_${{index}}"></div>
                `;
                grid.appendChild(card);

                const plotData = [];
                
                groups.forEach(g => {{
                    if(item[g] && item[g].y) {{
                        plotData.push({{
                            y: item[g].y,
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
                
                // Adjust layout for Grid
                const layout = {{
                    margin: {{ t: 20, r: 20, b: 30, l: 40 }},
                    yaxis: {{ title: 'Concentration [nmol/ml]', showgrid: true, gridcolor: '#eee' }},
                    paper_bgcolor: 'rgba(0,0,0,0)',
                    plot_bgcolor: 'rgba(0,0,0,0)',
                    showlegend: false,
                    font: {{ family: 'Segoe UI, sans-serif' }}
                }};

                Plotly.newPlot(`plot_${{index}}`, plotData, layout, {{responsive: true, displayModeBar: false}});
            }});
        }}

        // Run
        setTimeout(init, 100);
    </script>
</body>
</html>"""
    
    # Inject JSON data
    html = html.replace('__JSON_DATA__', json_content)
    
    return html
