"""HTML Generator for Node 10: Clinical Mimicry Dashboard
Generates interactive HTML visualization matching Node 4's Card Grid style.
"""

import json
import os
import io

def generate_mimicry_html(json_filename='mimicry_data.json'):
    """Generate HTML for Clinical Mimicry visualization"""
    
    # Locate JSON
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(SCRIPT_DIR, 'outputs', json_filename)
    
    if not os.path.exists(json_path):
        json_path = json_filename # Fallback

    try:
        with io.open(json_path, 'r', encoding='utf-8') as f:
            json_content = f.read()
    except Exception as e:
        print("Warning: Could not read JSON: " + str(e))
        json_content = '{}'

    html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Clinical Mimicry Test: RRMS vs MG</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {
            --primary: #2c3e50;
            --accent: #3498db;
            --bg: #f8f9fa;
            --card-bg: #ffffff;
            --text: #333;
        }
        
        body { font-family: 'Segoe UI', Roboto, sans-serif; margin: 0; padding: 0; background: var(--bg); color: var(--text); height: 100vh; display: flex; flex-direction: column; }
        
        /* Header */
        header { background: var(--card-bg); padding: 15px 30px; border-bottom: 1px solid #ddd; display: flex; justify-content: space-between; align-items: center; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }
        h1 { margin: 0; font-size: 1.2rem; color: var(--primary); }
        .timestamp { font-size: 0.8rem; color: #777; }

        /* Controls Bar */
        .controls { background: #fff; padding: 10px 30px; border-bottom: 1px solid #eee; display: flex; gap: 20px; align-items: center; flex-wrap: wrap; }
        
        .tab-group { display: flex; gap: 5px; background: #eee; padding: 4px; border-radius: 6px; }
        .tab-btn { border: none; background: none; padding: 8px 16px; cursor: pointer; border-radius: 4px; font-weight: 500; color: #666; transition: 0.2s; }
        .tab-btn.active { background: white; color: var(--accent); box-shadow: 0 1px 3px rgba(0,0,0,0.1); }

        .toggle-container { display: flex; align-items: center; gap: 8px; font-size: 0.9rem; margin-left: auto; }
        input[type="checkbox"] { accent-color: var(--accent); cursor: pointer; transform: scale(1.2); }

        .legend-bar { display: flex; gap: 15px; font-size: 0.85rem; margin-left: 20px; }
        .legend-item { display: flex; align-items: center; gap: 5px; }
        .dot { width: 10px; height: 10px; border-radius: 50%; display: inline-block; }

        /* Main Content Grid */
        #app { flex: 1; overflow-y: auto; padding: 20px 30px; }
        
        .grid-container {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(320px, 1fr));
            gap: 20px;
            padding-bottom: 50px;
        }
        
        .hero-container {
             max-width: 800px;
             margin: 0 auto;
             height: 500px;
        }

        .chart-card {
            background: var(--card-bg);
            border-radius: 8px;
            padding: 10px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.05);
            transition: transform 0.2s;
            border: 1px solid transparent;
            height: 350px;
            display: flex;
            flex-direction: column;
        }
        
        .chart-card:hover { transform: translateY(-2px); box-shadow: 0 5px 15px rgba(0,0,0,0.1); }
        
        .card-header { display: flex; justify-content: space-between; padding: 0 10px; margin-bottom: 5px; font-size: 0.9rem; font-weight: bold; color: var(--primary); }
        .stats-badge { font-size: 0.75rem; background: #eee; padding: 2px 6px; border-radius: 4px; font-weight: normal; color: #666; }
        .sig-badge { background: #fee; color: #c0392b; font-weight: bold; }
        
        .plot-div { flex: 1; width: 100%; }

        .loading { text-align: center; padding: 50px; font-size: 1.2rem; color: #666; }
    </style>
    <script>
        const RAW_DATA = __JSON_DATA__;
    </script>
</head>
<body>
    <header>
        <h1>Clinical Mimicry Test: RRMS vs MG</h1>
        <span class="timestamp" id="timestamp">Generated: Just now</span>
    </header>

    <div class="controls">
        <div class="tab-group">
            <button class="tab-btn active" onclick="switchTab('overview')">Total Overview</button>
            <button class="tab-btn" onclick="switchTab('detailed')">Specific Amino Acids</button>
        </div>

        <div class="legend-bar">
            <div class="legend-item"><span class="dot" style="background:grey"></span> RRMS</div>
            <div class="legend-item"><span class="dot" style="background:#87CEEB"></span> MG</div>
        </div>

        <div class="toggle-container">
            <input type="checkbox" id="sigFilter" onchange="renderCurrentView()">
            <label for="sigFilter" title="Show only correlations with p < 0.05">Significant Only (p<0.05)</label>
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
        let currentTab = 'overview';
        
        // Colors from Python script
        const COLORS = {
            'RRMS': 'grey', 
            'MG': '#87CEEB' // Skyblue
        };

        function init() {
            document.getElementById('timestamp').textContent = new Date().toLocaleString();
            renderCurrentView();
        }
        
        function switchTab(tab) {
            currentTab = tab;
            document.querySelectorAll('.tab-btn').forEach(b => b.classList.remove('active'));
            if(tab === 'overview') document.querySelector('.tab-btn:first-child').classList.add('active');
            else document.querySelector('.tab-btn:last-child').classList.add('active');
            renderCurrentView();
        }

        function renderCurrentView() {
            const container = document.getElementById('app');
            container.innerHTML = ''; // basic clean
            
            const isSigOnly = document.getElementById('sigFilter').checked;
            const showPoints = document.getElementById('showPoints').checked;

            if (currentTab === 'overview') {
                 // Hero Plot (Total AA)
                 const data = RAW_DATA['fig10'].traces[0];
                 const pVal = data.stats.p_value;
                 
                 // Create hero container
                 const heroDiv = document.createElement('div');
                 heroDiv.className = 'chart-card hero-container';
                 heroDiv.style.height = '600px';
                 
                 heroDiv.innerHTML = `
                    <div class="card-header" style="font-size:1.2rem;">
                        <span>Total Amino Acid Load</span>
                        <span class="stats-badge ${pVal < 0.05 ? 'sig-badge' : ''}" style="font-size:1rem;">P = ${pVal.toExponential(2)}</span>
                    </div>
                    <div class="plot-div" id="plot_total"></div>
                 `;
                 container.appendChild(heroDiv);
                 
                 drawBoxPlot('plot_total', data, showPoints, true);

            } else {
                 // Grid View
                 const traces = RAW_DATA['fig11'].traces;
                 container.innerHTML = '<div class="grid-container" id="grid"></div>';
                 const grid = document.getElementById('grid');

                 traces.forEach((item, index) => {
                     const pVal = item.stats ? item.stats.p_value : 1.0;
                     
                     if (isSigOnly && pVal >= 0.05) return;
                     
                     const card = document.createElement('div');
                     card.className = 'chart-card';
                     card.innerHTML = `
                        <div class="card-header">
                            <span>${item.aa}</span>
                            <span class="stats-badge ${pVal < 0.05 ? 'sig-badge' : ''}">P = ${pVal.toExponential(2)}</span>
                        </div>
                        <div class="plot-div" id="plot_${index}"></div>
                     `;
                     grid.appendChild(card);
                     
                     drawBoxPlot(`plot_${index}`, item, showPoints, false);
                 });
            }
        }
        
        function drawBoxPlot(divId, data, showPoints, showLegend) {
             const plotData = [];
             
             ['RRMS', 'MG'].forEach(group => {
                 if(data[group] && data[group].y) {
                     plotData.push({
                         y: data[group].y,
                         type: 'box',
                         name: group,
                         marker: { color: COLORS[group] },
                         boxpoints: showPoints ? 'all' : false,
                         jitter: 0.3,
                         pointpos: 0, 
                         showlegend: showLegend
                     });
                 }
             });
             
             const layout = {
                 margin: { t: 10, r: 10, b: 30, l: 40 },
                 yaxis: { title: 'Concentration', showgrid: true, gridcolor: '#eee' },
                 paper_bgcolor: 'rgba(0,0,0,0)',
                 plot_bgcolor: 'rgba(0,0,0,0)',
                 showlegend: showLegend
             };
             
             Plotly.newPlot(divId, plotData, layout, {responsive: true, displayModeBar: false});
        }

        // Run
        setTimeout(init, 100);
    </script>
</body>
</html>"""
    
    # Inject JSON data
    html = html.replace('__JSON_DATA__', json_content)
    
    return html
