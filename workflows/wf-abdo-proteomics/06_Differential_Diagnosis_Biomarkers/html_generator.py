"""HTML Generator for Node 06: Biomarkers Dashboard
Fixes the 'locking' issue by correctly handling Fig4/Fig5 data structures.
"""

import json
import os

def generate_biomarkers_html(json_filename='biomarkers_data.json'):
    """Generate HTML for Differential Diagnosis Biomarkers visualization"""
    
    # Locate the JSON file
    script_dir = os.path.dirname(os.path.abspath(__file__))
    json_path = os.path.join(script_dir, 'outputs', json_filename)
    if not os.path.exists(json_path):
        json_path = json_filename
        
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            json_content = f.read()
    except Exception as e:
        print(f"Warning: Could not read JSON file {json_path}: {e}")
        json_content = '{}'

    # The HTML Template
    html = """<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Biomarker Discovery Dashboard</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <style>
        :root {
            --primary: #2c3e50;
            --accent: #008080; /* Teal for MG */
            --bg: #f8f9fa;
            --card-bg: #ffffff;
            --text: #333;
        }
        
        body { font-family: 'Segoe UI', Roboto, sans-serif; margin: 0; padding: 0; background: var(--bg); color: var(--text); height: 100vh; display: flex; flex-direction: column; }
        
        /* Header */
        header { background: var(--card-bg); padding: 15px 30px; border-bottom: 1px solid #ddd; display: flex; justify-content: space-between; align-items: center; box-shadow: 0 2px 5px rgba(0,0,0,0.05); }
        h1 { margin: 0; font-size: 1.4rem; color: var(--primary); }
        .tag { font-size: 0.8rem; background: var(--accent); color: white; padding: 3px 8px; border-radius: 4px; vertical-align: middle; margin-left: 10px; }

        /* Navigation Tabs */
        .nav-tabs { display: flex; background: #fff; padding: 0 30px; border-bottom: 1px solid #ddd; gap: 20px; }
        .nav-item { 
            padding: 15px 5px; cursor: pointer; color: #666; font-weight: 500; border-bottom: 3px solid transparent; transition: 0.2s; 
        }
        .nav-item:hover { color: var(--accent); }
        .nav-item.active { color: var(--accent); border-bottom-color: var(--accent); }

        /* Controls */
        .controls { padding: 15px 30px; display: flex; justify-content: flex-end; align-items: center; }
        .radio-group { display: flex; background: #e9ecef; border-radius: 6px; padding: 3px; }
        .radio-group label { padding: 6px 14px; cursor: pointer; border-radius: 4px; font-weight: 500; color: #666; font-size: 0.9rem; transition: 0.2s; }
        .radio-group input { display: none; }
        .radio-group input:checked + label { background: white; color: var(--accent); box-shadow: 0 1px 3px rgba(0,0,0,0.1); }

        /* Content Grid */
        #app { flex: 1; overflow-y: auto; padding: 0 30px 30px 30px; }
        
        .grid-container {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(320px, 1fr));
            gap: 25px;
            max-width: 1400px;
            margin: 0 auto;
        }

        .chart-card {
            background: var(--card-bg);
            border-radius: 10px;
            padding: 15px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.05);
            height: 400px;
            display: flex;
            flex-direction: column;
            border: 1px solid #eee;
            min-width: 0; /* Prevents grid blowout */
        }
        
        .card-header { display: flex; justify-content: space-between; align-items: center; margin-bottom: 10px; padding: 0 5px; }
        .card-title { font-weight: 700; font-size: 1.1rem; color: var(--primary); }
        .p-badge { font-size: 0.85rem; background: #f1f3f5; padding: 4px 8px; border-radius: 6px; color: #555; }
        .p-sig { background: #ffe3e3; color: #c0392b; font-weight: bold; }

        .plot-div { flex: 1; width: 100%; min-width: 0; }
        
        .loading { text-align: center; padding: 50px; font-size: 1.2rem; color: #666; }
    </style>
    <script>
        const RAW_DATA = __JSON_DATA__;
    </script>
</head>
<body>
    <header>
        <div>
            <h1>Differential Diagnosis <span class="tag">MS vs MG</span></h1>
        </div>
        <div style="font-size: 0.8rem; color: #888;">Generated: <span id="timestamp"></span></div>
    </header>

    <div class="nav-tabs">
        <div class="nav-item active" onclick="switchTab('fig4', this)">General Cohort</div>
        <div class="nav-item" onclick="switchTab('fig5', this)">Female Subgroup</div>
    </div>

    <div class="controls">
        <div class="radio-group">
            <input type="radio" id="viewBox" name="viewType" value="box" checked onchange="renderDashboard()">
            <label for="viewBox">Box Plot</label>
            
            <input type="radio" id="viewViolin" name="viewType" value="violin" onchange="renderDashboard()">
            <label for="viewViolin">Violin Plot</label>
        </div>
    </div>

    <div id="app">
        <div class="loading">Initializing Dashboard...</div>
    </div>

    <script>
        // Start with Fig 4 (General Cohort)
        let currentSection = 'fig4';

        function init() {
            document.getElementById('timestamp').textContent = new Date().toLocaleString();
            
            // Check if data loaded correctly
            if (!RAW_DATA || !RAW_DATA.fig4) {
                document.getElementById('app').innerHTML = '<div class="loading">Error: Data not loaded correctly.</div>';
                return;
            }
            
            renderDashboard();
        }

        function switchTab(sectionKey, element) {
            currentSection = sectionKey;
            
            // Update Tab UI
            document.querySelectorAll('.nav-item').forEach(el => el.classList.remove('active'));
            element.classList.add('active');
            
            renderDashboard();
        }

        function renderDashboard() {
            const container = document.getElementById('app');
            
            // Get data for current section (fig4 or fig5)
            const sectionData = RAW_DATA[currentSection];
            
            if (!sectionData || !sectionData.subplots) {
                 container.innerHTML = '<div class="loading">No data found for this section.</div>';
                 return;
            }

            const data = sectionData.subplots;
            const chartType = document.querySelector('input[name="viewType"]:checked').value;
            
            // Clear container
            container.innerHTML = '<div class="grid-container" id="grid"></div>';
            const grid = document.getElementById('grid');

            // Render Charts
            data.forEach((item, index) => {
                const pVal = item.p_value ?? 1;
                const isSig = pVal < 0.05;
                const pText = pVal < 0.001 ? 'p < 0.001' : `p = ${pVal.toFixed(4)}`;

                // Create Card HTML
                const card = document.createElement('div');
                card.className = 'chart-card';
                card.innerHTML = `
                    <div class="card-header">
                        <span class="card-title">${item.title}</span>
                        <span class="p-badge ${isSig ? 'p-sig' : ''}">${pText}</span>
                    </div>
                    <div class="plot-div" id="plot_${index}"></div>
                `;
                grid.appendChild(card);

                // Create Traces for Plotly
                const traces = item.traces.map(t => ({
                    y: t.y,
                    name: t.name,
                    type: chartType,
                    boxpoints: 'outliers',
                    marker: { color: t.color, size: 4, opacity: 0.7 },
                    line: { width: 1.5 },
                    meanline: { visible: true },
                    side: 'positive'
                }));

                // Layout settings
                const layout = {
                    margin: { t: 10, r: 10, b: 30, l: 40 },
                    yaxis: { title: 'Concentration (Std)', zeroline: false },
                    showlegend: (index === 0), // Show legend only on the first chart
                    legend: { x: 1, y: 1, xanchor: 'right' },
                    paper_bgcolor: 'rgba(0,0,0,0)',
                    plot_bgcolor: 'rgba(0,0,0,0)'
                };

                Plotly.newPlot(`plot_${index}`, traces, layout, {responsive: true, displayModeBar: false});
            });
        }

        // Run initialization
        setTimeout(init, 100);
    </script>
</body>
</html>"""
    
    # Inject JSON data
    html = html.replace('__JSON_DATA__', json_content)
    
    return html
