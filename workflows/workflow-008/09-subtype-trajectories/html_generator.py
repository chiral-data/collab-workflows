"""HTML Generator for Node 09: Subtype Trajectories Dashboard
Generates interactive HTML visualization matching Node 4's Card Grid style.
"""

import json
import os
import io

def generate_trajectories_html(json_filename='trajectories_data.json'):
    """Generate HTML for Subtype Trajectories visualization"""
    
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
    <title>Subtype Trajectories Dashboard</title>
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
        
        .toggle-container { display: flex; align-items: center; gap: 8px; font-size: 0.9rem; margin-left: auto; }
        input[type="checkbox"] { accent-color: var(--accent); cursor: pointer; transform: scale(1.2); }

        .legend-bar { display: flex; gap: 15px; font-size: 0.85rem; }
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
        <h1>Subtype Trajectories Analysis</h1>
        <span class="timestamp" id="timestamp">Generated: Just now</span>
    </header>

    <div class="controls">
        <div class="legend-bar">
            <!-- Legend matching python colors -->
            <div class="legend-item"><span class="dot" style="background:#E74C3C"></span> PPMS</div>
            <div class="legend-item"><span class="dot" style="background:#3498DB"></span> SPMS</div>
            <div class="legend-item"><span class="dot" style="background:#2ECC71"></span> RRMS</div>
            <div class="legend-item"><span class="dot" style="background:#9B59B6"></span> GMG</div>
            <div class="legend-item"><span class="dot" style="background:#F39C12"></span> OMG</div>
        </div>

        <div class="toggle-container">
            <input type="checkbox" id="sigFilter" onchange="renderCurrentView()">
            <label for="sigFilter" title="Show only correlations with p < 0.05">Significant Only (p<0.05)</label>
        </div>
        
        <div class="toggle-container">
            <input type="checkbox" id="trendlineToggle" checked onchange="renderCurrentView()">
            <label for="trendlineToggle">Show Trendlines</label>
        </div>
    </div>

    <div id="app">
        <div class="loading">Initializing Dashboard...</div>
    </div>

    <script>
        // Colors from Python script
        const COLORS = {
            'PPMS': '#E74C3C', 
            'SPMS': '#3498DB', 
            'RRMS': '#2ECC71', 
            'GMG': '#9B59B6', 
            'OMG': '#F39C12'
        };

        function init() {
            document.getElementById('timestamp').textContent = new Date().toLocaleString();
            renderCurrentView();
        }

        function renderCurrentView() {
            const container = document.getElementById('app');
            const dataSection = RAW_DATA['fig9']; // Only one fig currently
            const traces = dataSection.traces;
            const isSigOnly = document.getElementById('sigFilter').checked;
            const showLines = document.getElementById('trendlineToggle').checked;
            
            // Clean container
            container.innerHTML = '<div class="grid-container" id="grid"></div>';
            const grid = document.getElementById('grid');

            // Process each AA
            traces.forEach((item, index) => {
                const aaName = item.aa;
                const subtypes = ['PPMS', 'SPMS', 'RRMS', 'GMG', 'OMG'];
                
                // Determine significance (if ANY subtype has p < 0.05)
                let isSig = false;
                let minP = 1.0;
                
                subtypes.forEach(st => {
                    if (item[st] && item[st].stats && item[st].stats.p_value < 0.05) {
                        isSig = true;
                        if(item[st].stats.p_value < minP) minP = item[st].stats.p_value;
                    }
                });

                if (isSigOnly && !isSig) return; // Skip if filter is on

                // Create Card
                const card = document.createElement('div');
                card.className = 'chart-card';
                
                // Header info
                let statsText = isSig ? `Min P: ${minP.toExponential(2)}` : 'ns';
                
                card.innerHTML = `
                    <div class="card-header">
                        <span>${aaName}</span>
                        <span class="stats-badge ${isSig ? 'sig-badge' : ''}">${statsText}</span>
                    </div>
                    <div class="plot-div" id="plot_${index}"></div>
                `;
                grid.appendChild(card);

                // Prepare Plotly Data
                const plotData = [];
                
                subtypes.forEach(st => {
                     const groupData = item[st];
                     if(!groupData) return;

                     // Scatter points
                     plotData.push({
                         x: groupData.x,
                         y: groupData.y,
                         mode: 'markers',
                         type: 'scatter',
                         name: st,
                         marker: { color: COLORS[st], size: 5, opacity: 0.6 },
                         showlegend: false
                     });

                     // Trendline
                     if(showLines && groupData.stats) {
                         plotData.push({
                             x: groupData.stats.line_x,
                             y: groupData.stats.line_y,
                             mode: 'lines',
                             type: 'scatter',
                             name: `${st} Fit`,
                             line: { color: COLORS[st], width: 2 },
                             showlegend: false,
                             hoverinfo: 'skip'
                         });
                     }
                });

                // Layout
                const layout = {
                    margin: { t: 10, r: 10, b: 30, l: 40 },
                    xaxis: { title: dataSection.variable, showgrid: true, gridcolor: '#eee' },
                    yaxis: { title: 'Conc', showgrid: true, gridcolor: '#eee' },
                    paper_bgcolor: 'rgba(0,0,0,0)',
                    plot_bgcolor: 'rgba(0,0,0,0)',
                    showlegend: false, // Legend is in header
                    hovermode: 'closest'
                };

                Plotly.newPlot(`plot_${index}`, plotData, layout, {responsive: true, displayModeBar: false});
            });
        }

        // Run
        setTimeout(init, 100);
    </script>
</body>
</html>"""
    
    # Inject JSON data
    html = html.replace('__JSON_DATA__', json_content)
    
    return html
