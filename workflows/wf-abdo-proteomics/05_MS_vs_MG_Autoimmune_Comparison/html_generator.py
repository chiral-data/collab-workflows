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
                    boxpoints: 'all', // Show all points
                    jitter: 0.3,
                    pointpos: -1.8,
                    marker: { color: t.color, size: 3, opacity: 0.6 },
                    line: { width: 1.5 },
                    // Violin specific settings
                    side: 'positive',
                    meanline: { visible: true }
                }));

                const layout = {
                    margin: { t: 10, r: 10, b: 30, l: 40 },
                    yaxis: { title: 'Standardized Level', zeroline: false },
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
Generates interactive HTML visualization
Based on JavaScript + Plotly.js architecture
"""

import json
import os

def generate_autoimmune_html(json_filename='autoimmune_data.json'):
    """Generate HTML for MS vs MG Autoimmune Comparison visualization"""
    
    # Read JSON data to embed directly
    json_path = os.path.join('outputs', json_filename)
    # Fallback if running from inside output dir or elsewhere
    if not os.path.exists(json_path):
        json_path = json_filename
        
    try:
        with open(json_path, 'r', encoding='utf-8') as f:
            json_content = f.read()
    except Exception as e:
        print(f"Warning: Could not read JSON file {json_path} for embedding: {e}")
        json_content = '{}'

    html = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>MS vs MG Autoimmune Comparison</title>
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js"></script>
    <script>
        // Embedded JSON data
        const RAW_DATA = __JSON_DATA__;
    </script>
    <style>
        :root {
            --primary: #2c3e50;
            --secondary: #667eea;
            --gradient-start: #667eea;
            --gradient-end: #764ba2;
            --bg-light: #f5f5f5;
            --text-muted: #7f8c8d;
        }
        
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
            background: linear-gradient(135deg, var(--gradient-start) 0%, var(--gradient-end) 100%);
            padding: 20px;
            line-height: 1.6;
            min-height: 100vh;
        }
        
        .container {
            max-width: 1800px;
            margin: 0 auto;
            background: white;
            padding: 40px;
            border-radius: 12px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.15);
        }
        
        header {
            text-align: center;
            margin-bottom: 50px;
            padding-bottom: 30px;
            border-bottom: 3px solid var(--secondary);
        }
        
        h1 {
            color: var(--primary);
            font-size: 2.5em;
            margin-bottom: 10px;
            font-weight: 700;
        }
        
        .subtitle {
            color: var(--text-muted);
            font-size: 1.1em;
            margin-top: 10px;
        }
        
        .figure-section {
            margin: 60px 0;
            animation: fadeIn 0.6s ease-in;
        }
        
        @keyframes fadeIn {
            from {
                opacity: 0;
                transform: translateY(20px);
            }
            to {
                opacity: 1;
                transform: translateY(0);
            }
        }
        
        .figure-title {
            font-size: 1.6em;
            color: var(--primary);
            font-weight: 600;
            margin-bottom: 25px;
            padding-left: 20px;
            border-left: 5px solid var(--secondary);
        }
        
        .plot-container {
            background: var(--bg-light);
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.08);
            min-height: 300px;
        }
        
        .loading {
            text-align: center;
            padding: 100px 20px;
        }
        
        .spinner {
            border: 4px solid #f3f3f3;
            border-top: 4px solid var(--secondary);
            border-radius: 50%;
            width: 60px;
            height: 60px;
            animation: spin 1s linear infinite;
            margin: 0 auto 20px;
        }
        
        @keyframes spin {
            100% { transform: rotate(360deg); }
        }
        
        .error {
            background: linear-gradient(135deg, #ff6b6b 0%, #ee5a6f 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin: 20px 0;
            box-shadow: 0 4px 15px rgba(255, 107, 107, 0.3);
        }
        
        footer {
            margin-top: 80px;
            text-align: center;
            color: var(--text-muted);
            padding-top: 30px;
            border-top: 2px solid #ecf0f1;
            font-size: 0.95em;
        }
        
        @media (max-width: 1400px) {
            .plot-container[style*="grid-template-columns"] {
                grid-template-columns: repeat(3, 1fr) !important;
            }
        }
        
        @media (max-width: 900px) {
            .container {
                padding: 25px;
            }
            h1 {
                font-size: 2em;
            }
            .plot-container[style*="grid-template-columns"] {
                grid-template-columns: repeat(2, 1fr) !important;
            }
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>🔬 MS vs MG Autoimmune Comparison</h1>
            <p class="subtitle">MS & MG Amino Acid Analysis</p>
        </header>
        
        <div id="app">
            <div class="loading">
                <div class="spinner"></div>
                <p style="color: var(--text-muted); font-size: 1.1em;">Loading data and generating visualizations...</p>
            </div>
        </div>
        
        <footer>
            <p><strong>Dataset:</strong> Multiple Sclerosis & Myasthenia Gravis Proteomics</p>
            <p><strong>Generated:</strong> <span id="timestamp"></span></p>
            <p><strong>Node:</strong> 05 - MS vs MG Autoimmune Comparison</p>
        </footer>
    </div>
    
    <script>
        'use strict';
        
        const CONFIG = {
            colors: {
                case: '#4682b4',
                control: '#fa8072',
                ms: '#4682b4',
                mg: '#20b2aa',
                ppms: '#800080',
                spms: '#ffa500',
                rrms: '#90ee90',
                gmg: '#f17cb0',
                omg: '#b2912f'
            },
            plotly: {
                responsive: true,
                displayModeBar: true,
                modeBarButtons: [
                    ['zoom2d', 'pan2d', 'autoScale2d', 'resetScale2d'],
                    ['toImage']
                ],
                displaylogo: false
            }
        };
        
        function renderGridPlot(containerId, data) {
            if (!data || !data.subplots) {
                console.error('No subplot data for ' + containerId);
                return;
            }
            
            const subplots = data.subplots;
            const container = document.getElementById(containerId);
            
            container.style.display = 'grid';
            container.style.gridTemplateColumns = 'repeat(5, 1fr)';
            container.style.gap = '12px';
            container.style.padding = '10px';
            container.style.background = '#f5f5f5';
            
            subplots.forEach((subplot, idx) => {
                const div = document.createElement('div');
                div.id = containerId + '_' + idx;
                div.style.minHeight = '250px';
                div.style.background = 'white';
                div.style.borderRadius = '8px';
                div.style.border = '1px solid #d0d0d0';
                div.style.boxShadow = '0 2px 4px rgba(0,0,0,0.1)';
                container.appendChild(div);
                
                const traces = subplot.traces.map((t, tIdx) => ({
                    y: t.y,
                    name: t.name,
                    type: 'box',
                    marker: { 
                        color: t.color || CONFIG.colors[t.name.toLowerCase()] || '#3498db',
                        line: { color: 'black', width: 1 }
                    },
                    boxmean: 'sd',
                    showlegend: idx === 0
                }));
                
                const layout = {
                    title: { text: subplot.title, font: { size: 12 }, y: 0.95 },
                    showlegend: idx === 0,
                    legend: { x: 0, y: -0.25, orientation: 'h', font: { size: 9 } },
                    height: 250,
                    margin: { t: 45, r: 15, b: 50, l: 50 },
                    plot_bgcolor: '#ebebeb',
                    paper_bgcolor: 'white',
                    xaxis: { showticklabels: true, tickfont: { size: 9 }, gridcolor: '#d0d0d0' },
                    yaxis: { 
                        title: 'Log C',
                        titlefont: { size: 10 },
                        tickfont: { size: 9 },
                        gridcolor: '#d0d0d0',
                        zeroline: true
                    }
                };
                
                Plotly.newPlot(div.id, traces, layout, CONFIG.plotly);
            });
        }
        
        function renderBoxPlot(containerId, data) {
            if (!data || !data.traces) {
                console.error('No data for ' + containerId);
                return;
            }
            
            const traces = data.traces.map(trace => ({
                y: trace.y,
                name: trace.name,
                type: 'box',
                boxmean: 'sd',
                marker: { 
                    color: trace.color || CONFIG.colors[trace.name.toLowerCase()] || '#3498db',
                    line: { color: 'black', width: 1 }
                }
            }));
            
            const layout = {
                title: data.title || '',
                xaxis: { title: data.xaxis || '' },
                yaxis: { title: data.yaxis || 'Value' },
                showlegend: true,
                height: 500,
                plot_bgcolor: '#ebebeb',
                paper_bgcolor: 'white',
                margin: { t: 50, r: 20, b: 40, l: 50 }
            };
            
            Plotly.newPlot(containerId, traces, layout, CONFIG.plotly);
        }
        
        function init() {
            try {
                const data = RAW_DATA;
                console.log('Data loaded successfully:', Object.keys(data));
                
                // Build HTML based on data structure
                let html = '';
                for (const key in data) {
                    if (key !== 'metadata') {
                        html += `
                            <div class="figure-section">
                                <h2 class="figure-title">${key.toUpperCase().replace('_', ' ')}</h2>
                                <div class="plot-container" id="${key}"></div>
                            </div>
                        `;
                    }
                }
                
                document.getElementById('app').innerHTML = html;
                document.getElementById('timestamp').textContent = new Date().toLocaleString();
                
                setTimeout(() => {
                    console.log('Rendering visualizations...');
                    for (const key in data) {
                        if (key !== 'metadata') {
                            if (data[key].subplots) {
                                renderGridPlot(key, data[key]);
                            } else if (data[key].traces) {
                                renderBoxPlot(key, data[key]);
                            }
                        }
                    }
                    console.log('All visualizations rendered successfully!');
                }, 150);
                
            } catch (error) {
                console.error('Error during initialization:', error);
                document.getElementById('app').innerHTML = `
                    <div class="error">
                        <h2>⚠️ Error Loading Data</h2>
                        <p><strong>Message:</strong> ${error.message}</p>
                        <p><strong>File:</strong> autoimmune_data.json</p>
                        <p>Please ensure the JSON file exists in the same directory.</p>
                    </div>
                `;
            }
        }
        
        init();
    </script>
</body>
</html>'''
    
    # Inject JSON data
    html = html.replace('__JSON_DATA__', json_content)
    
    return html

