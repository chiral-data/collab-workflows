"""HTML Generator for Node 05: Autoimmune Dashboard
Generates interactive HTML visualization with Box/Violin toggles and Significance filters.
"""

import json
import os

def generate_autoimmune_html(json_filename='autoimmune_data.json'):
    """Generate HTML for MS vs MG Autoimmune visualization"""
    
    json_path = os.path.join('outputs', json_filename)
    if not os.path.exists(json_path):
        json_path = json_filename
        
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            json_content = f.read()
    except Exception as e:
        print(f"Warning: Could not read JSON file {json_path}: {e}")
        json_content = '{}'

    html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Autoimmune Analysis Dashboard</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {
            --primary: #2c3e50;
            --accent: #3498db;
            --bg: #f8f9fa;
            --card-bg: #ffffff;
            --text: #333;
        }
        
        body { font-family: 'Segoe UI', Roboto, Helvetica, Arial, sans-serif; margin: 0; padding: 0; background: var(--bg); color: var(--text); height: 100vh; display: flex; flex-direction: column; }
        
        /* Header */
        header { background: var(--card-bg); padding: 15px 30px; border-bottom: 1px solid #ddd; display: flex; justify-content: space-between; align-items: center; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }
        h1 { margin: 0; font-size: 1.4rem; color: var(--primary); }
        .subtitle { font-size: 0.9rem; color: #666; margin-left: 10px; }

        /* Controls Bar */
        .controls { background: #fff; padding: 12px 30px; border-bottom: 1px solid #eee; display: flex; gap: 25px; align-items: center; flex-wrap: wrap; }
        
        .control-group { display: flex; align-items: center; gap: 8px; font-size: 0.9rem; }
        .control-label { font-weight: 600; color: #555; margin-right: 5px; }
        
        /* Styled Radio Buttons as Tabs */
        .radio-group { display: flex; background: #eee; border-radius: 6px; padding: 3px; }
        .radio-group label {
            padding: 6px 14px; cursor: pointer; border-radius: 4px; font-weight: 500; color: #666; transition: 0.2s; user-select: none;
        }
        .radio-group input { display: none; }
        .radio-group input:checked + label { background: white; color: var(--accent); box-shadow: 0 1px 3px rgba(0,0,0,0.1); }

        /* Checkbox */
        input[type="checkbox"] { accent-color: var(--accent); transform: scale(1.2); cursor: pointer; }

        /* Main Content Grid */
        #app { flex: 1; overflow-y: auto; padding: 20px 30px; }
        
        .grid-container {
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
            gap: 20px;
            padding-bottom: 50px;
        }

        /* Chart Card */
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
        .chart-card.highlight { border-color: var(--accent); }
        
        .card-header { display: flex; justify-content: space-between; align-items: center; padding: 0 10px; margin-bottom: 5px; }
        .card-title { font-weight: bold; color: var(--primary); }
        
        /* P-value Badge */
        .stats-badge { font-size: 0.75rem; background: #eee; padding: 2px 8px; border-radius: 10px; color: #666; }
        .sig-badge { background: #fee; color: #c0392b; font-weight: bold; border: 1px solid #f5b7b1; }
        
        .plot-div { flex: 1; width: 100%; }
        .loading { text-align: center; padding: 50px; font-size: 1.2rem; color: #666; }
    </style>
    <script>
        const RAW_DATA = __JSON_DATA__;
    </script>
</head>
<body>
    <header>
        <div>
            <h1>Autoimmune Comparison <span class="subtitle">MS+MG vs Controls</span></h1>
        </div>
        <div style="font-size: 0.8rem; color: #888;">Generated: <span id="timestamp"></span></div>
    </header>

    <div class="controls">
        <div class="control-group">
            <span class="control-label">Chart Type:</span>
            <div class="radio-group">
                <input type="radio" id="viewBox" name="viewType" value="box" checked onchange="renderDashboard()">
                <label for="viewBox">Box Plot</label>
                
                <input type="radio" id="viewViolin" name="viewType" value="violin" onchange="renderDashboard()">
                <label for="viewViolin">Violin Plot</label>
            </div>
        </div>

        <div class="control-group" style="margin-left: auto;">
            <input type="checkbox" id="sigFilter" onchange="renderDashboard()">
            <label for="sigFilter" title="Show only results with p < 0.05">Significant Only (p < 0.05)</label>
        </div>
    </div>

    <div id="app">
        <div class="loading">Loading analysis data...</div>
    </div>

    <script>
        function init() {
            document.getElementById('timestamp').textContent = new Date().toLocaleString();
            renderDashboard();
        }

        function renderDashboard() {
            const container = document.getElementById('app');
            const data = RAW_DATA.fig3.subplots;
            const isSigOnly = document.getElementById('sigFilter').checked;
            const chartType = document.querySelector('input[name="viewType"]:checked').value;
            
            // Clear previous content
            container.innerHTML = '<div class="grid-container" id="grid"></div>';
            const grid = document.getElementById('grid');

            // Render each chart
            data.forEach((item, index) => {
                // Filter Logic
                const pVal = item.p_value ?? 1.0;
                const isSig = pVal < 0.05;
                
                if (isSigOnly && !isSig) return;

                // Create Card Element
                const card = document.createElement('div');
                card.className = `chart-card ${isSig ? 'highlight' : ''}`;
                
                const pText = pVal < 0.001 ? 'p < 0.001' : `p = ${pVal.toFixed(3)}`;
                
                card.innerHTML = `
                    <div class="card-header">
                        <span class="card-title">${item.title}</span>
                        <span class="stats-badge ${isSig ? 'sig-badge' : ''}">${pText}</span>
                    </div>
                    <div class="plot-div" id="plot_${index}"></div>
                `;
                grid.appendChild(card);

                // Prepare Plotly Traces
                const traces = item.traces.map(t => ({
                    y: t.y,
                    name: t.name,
                    type: chartType, // 'box' or 'violin'
                    boxpoints: 'outliers', // Show only outliers (Standard Box Plot)
                    marker: { color: t.color, size: 3, opacity: 0.6 },
                    line: { width: 1.5 },
                    // Violin specific settings
                    side: 'positive',
                    meanline: { visible: true }
                }));

                const layout = {
                    margin: { t: 30, r: 10, b: 30, l: 50 },
                    yaxis: { title: 'Log C', zeroline: false },  // Original Axis Label
                    showlegend: (index === 0), // Only show legend on first chart
                    legend: { x: 1, y: 1, xanchor: 'right' },
                    paper_bgcolor: 'rgba(0,0,0,0)',
                    plot_bgcolor: 'rgba(0,0,0,0)'
                };
                
                // Remove X-axis labels to clean up the grid
                layout.xaxis = { showticklabels: true };

                Plotly.newPlot(`plot_${index}`, traces, layout, {responsive: true, displayModeBar: false});
            });
            
            // Empty State Message
            if (grid.children.length === 0) {
                container.innerHTML = '<div class="loading">No significant results found matching criteria.</div>';
            }
        }

        // Run initialization
        setTimeout(init, 100);
    </script>
</body>
</html>"""
    
    # Inject JSON data
    html = html.replace('__JSON_DATA__', json_content)
    
    return html
