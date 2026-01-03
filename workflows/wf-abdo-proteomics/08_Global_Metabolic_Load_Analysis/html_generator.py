"""HTML Generator for Node 08: Global Metabolic Load Analysis"""
import os
import json
import io

def generate_metabolic_load_html(json_filename='metabolic_load_data.json'):
    
    # Locate JSON (Logic from your original file)
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(SCRIPT_DIR, 'outputs', json_filename)
    
    if not os.path.exists(json_path):
        json_path = json_filename # Fallback

    try:
        # Use io.open for Python 2 compatibility
        with io.open(json_path, 'r', encoding='utf-8') as f:
            json_content = f.read()
    except Exception as e:
        print("Warning: Could not read JSON: " + str(e))
        json_content = '{}'

    html = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Global Metabolic Load Analysis</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.27.0/plotly.min.js"></script>
    <style>
        :root {
            --sidebar-width: 280px;
            --primary: #2c3e50;
            --accent: #3498db;
            --bg: #f8f9fa;
        }
        body { margin: 0; padding: 0; font-family: 'Segoe UI', sans-serif; background: var(--bg); display: flex; height: 100vh; overflow: hidden; }
        
        /* Sidebar Styles */
        .sidebar {
            width: var(--sidebar-width);
            background: white;
            border-right: 1px solid #ddd;
            padding: 20px;
            display: flex;
            flex-direction: column;
            box-shadow: 2px 0 5px rgba(0,0,0,0.05);
            z-index: 10;
        }
        .sidebar h2 { font-size: 1.2rem; color: var(--primary); margin-bottom: 20px; border-bottom: 2px solid var(--accent); padding-bottom: 10px; }
        .control-group { margin-bottom: 25px; }
        .control-group label { display: block; font-weight: 600; margin-bottom: 10px; color: #555; font-size: 0.9rem; }
        
        select, input[type="checkbox"] { margin-bottom: 10px; }
        select {
            width: 100%; padding: 10px; border: 1px solid #ddd; border-radius: 6px; font-size: 1rem;
            background-color: white; cursor: pointer; transition: border 0.2s;
        }
        select:hover { border-color: var(--accent); }
        
        .stats-panel {
            margin-top: auto; padding: 15px; background: #eef7fb; border-radius: 8px; font-size: 0.9rem; color: #444;
            border-left: 4px solid var(--accent);
        }
        .stats-panel h3 { margin: 0 0 10px 0; font-size: 1rem; color: var(--primary); }
        .stat-item { margin-bottom: 5px; display: flex; justify-content: space-between; }
        .stat-value { font-weight: bold; color: var(--accent); }

        /* Main Content */
        .main { flex: 1; padding: 20px; display: flex; flex-direction: column; position: relative; }
        #plot-container { flex: 1; width: 100%; background: white; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,0.05); overflow: hidden; }
        
        .loading { position: absolute; top: 50%; left: 55%; transform: translate(-50%, -50%); font-size: 1.2rem; color: #777; }
    </style>
</head>
<body>

    <div class="sidebar">
        <h2>Metabolic Load</h2>
        
        <div class="control-group">
            <label for="chart-type">Visualization Type:</label>
            <select id="chart-type" onchange="updatePlot()">
                <option value="box">Box Plot</option>
                <option value="violin">Violin Plot</option>
            </select>
        </div>

        <div class="control-group">
            <label for="points-toggle">Data Points:</label>
            <select id="points-toggle" onchange="updatePlot()">
                <option value="outliers">Outliers Only (Clean)</option>
                <option value="all">Show All Points (Jitter)</option>
                <option value="false" selected>Hidden</option>
            </select>
        </div>

        <div class="control-group" style="border-top: 1px solid #eee; padding-top: 15px;">
            <label class="checkbox-label" style="display: flex; align-items: center; cursor: pointer;">
                <input type="checkbox" id="sig-filter" onchange="updatePlot()" style="margin: 0 10px 0 0;">
                <span>Filter Non-Significant</span>
            </label>
        </div>

        <div class="stats-panel" id="stats-panel">
            <h3>Statistical Significance</h3>
            <div class="stat-item">
                <span>Test:</span> <span id="stat-test" style="font-weight: 600;">Loading...</span>
            </div>
            <div class="stat-item">
                <span>P-Value:</span> <span id="stat-p" class="stat-value" style="font-size: 1.2rem;">...</span>
            </div>
        </div>
    </div>

    <div class="main">
        <div id="plot-container"></div>
        <div id="loading" class="loading">Loading Data...</div>
        <div id="filter-msg" style="display:none; text-align: center; margin-top: 50px; color: #7f8c8d;">
            <h3>Result Not Significant</h3>
            <p>Chart hidden by filter (P > 0.05)</p>
        </div>
    </div>

    <script>
        const RAW_DATA = __JSON_DATA__;
        
        // Helper to parse p-value
        function isSignificant(pValStr) {
            if (typeof pValStr !== 'string') return false;
            if (pValStr.includes('<')) return true; 
            return parseFloat(pValStr) < 0.05;
        }

        function init() {
            document.getElementById('loading').style.display = 'none';
            if(!RAW_DATA.fig8 || !RAW_DATA.fig8.traces) {
                alert("Error: Data structure invalid.");
                return;
            }
            
            // Populate Stats
            if (RAW_DATA.stats) {
                const pVal = RAW_DATA.stats.p_value_fmt;
                const pElem = document.getElementById('stat-p');
                document.getElementById('stat-test').textContent = RAW_DATA.stats.test;
                pElem.textContent = pVal;
                
                // Color logic
                if (isSignificant(pVal)) {
                    pElem.style.color = '#27ae60'; // Green
                    pElem.textContent += ' (Sig)';
                } else {
                    pElem.style.color = '#7f8c8d'; // Grey
                    pElem.textContent += ' (ns)';
                }
            }

            updatePlot();
        }

        function updatePlot() {
            var chartType = document.getElementById('chart-type').value;
            var pointsMode = document.getElementById('points-toggle').value;
            var filterOn = document.getElementById('sig-filter').checked;
            
            // P-Value Filter Logic
            var pVal = RAW_DATA.stats ? RAW_DATA.stats.p_value_fmt : "1.0";
            var isSig = isSignificant(pVal);
            
            if (filterOn && !isSig) {
                document.getElementById('plot-container').style.display = 'none';
                document.getElementById('filter-msg').style.display = 'block';
                return; // Stop rendering
            } else {
                document.getElementById('plot-container').style.display = 'block';
                document.getElementById('filter-msg').style.display = 'none';
            }

            var dataConfig = RAW_DATA.fig8;
            
            var traces = dataConfig.traces.map(function(t) {
                var trace = {
                    y: t.y,
                    name: t.name,
                    type: chartType,
                    marker: { color: t.color },
                    line: { color: 'black', width: 1.5 }
                };

                if (chartType === 'box') {
                    trace.boxpoints = pointsMode;
                    trace.jitter = 0.3;
                    trace.pointpos = 0; // Centered
                    trace.fillcolor = t.color; 
                } else if (chartType === 'violin') {
                    trace.points = pointsMode;
                    trace.side = 'positive';
                    trace.meanline = { visible: true };
                    trace.fillcolor = t.color; 
                    trace.line = { color: 'black', width: 1 };
                }
                
                return trace;
            });

            var layout = {
                title: { text: dataConfig.title, font: {size: 20} },
                yaxis: { title: dataConfig.yaxis, zeroline: false },
                margin: { l: 60, r: 30, t: 80, b: 60 },
                showlegend: true,
                legend: { orientation: 'h', x: 0.5, xanchor: 'center', y: -0.1 }
            };

            Plotly.newPlot('plot-container', traces, layout, {responsive: true});
        }

        init();
    </script>
</body>
</html>'''

    # Inject JSON data
    html = html.replace('__JSON_DATA__', json_content)
    
    return html
