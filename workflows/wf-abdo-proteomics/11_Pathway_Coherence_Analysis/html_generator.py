"""HTML Generator for Node 11: Pathway Coherence Analysis
Generates interactive HTML visualization
Based on JavaScript + Plotly.js architecture
"""

import json
import os

def generate_coherence_html(json_filename='coherence_data.json'):
    """Generate HTML for Pathway Coherence Analysis visualization"""
    
    # Read JSON data to embed directly
    json_path = os.path.join('output', json_filename)
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
    <title>Pathway Coherence Analysis</title>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/plotly.js/2.27.0/plotly.min.js" crossorigin="anonymous"></script>
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
            <h1>🔬 Pathway Coherence Analysis</h1>
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
            <p><strong>Node:</strong> 11 - Pathway Coherence Analysis</p>
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
        
        function renderHeatmap(containerId, data) {
            if (!data || !data.z) {
                console.error('No z data for ' + containerId);
                return;
            }
            
            const trace = {
                z: data.z,
                x: data.x,
                y: data.y,
                type: 'heatmap',
                colorscale: data.colorscale || 'RdBu',
                zmin: -0.2,
                zmax: 1,
                xgap: 1,
                ygap: 1,
                text: data.z,
                texttemplate: '%{text:.2f}',
                textfont: { size: 10 }
            };
            
            const layout = {
                title: data.title || '',
                height: 900,
                width: 900,
                xaxis: { tickangle: -45 },
                yaxis: { automargin: true, autorange: 'reversed' },
                margin: { t: 50, r: 50, b: 100, l: 100 }
            };
            
            Plotly.newPlot(containerId, [trace], layout, CONFIG.plotly);
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
                            } else if (data[key].z) {
                                renderHeatmap(key, data[key]);
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
                        <p><strong>File:</strong> coherence_data.json</p>
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

