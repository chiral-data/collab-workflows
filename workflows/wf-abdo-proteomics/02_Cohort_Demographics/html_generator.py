"""HTML Generator for Node 02: Cohort Demographics
Generates interactive HTML visualization for demographic tables
Based on JavaScript + Plotly.js architecture
"""

import json
import os

def generate_demographics_html(json_filename='demographics_data.json'):
    """Generate HTML for demographics tables visualization"""
    
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
    <title>Cohort Demographics - Tables 1 & 2</title>
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
            max-width: 1400px;
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
        
        .data-table {
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
            font-size: 0.95em;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }
        
        .data-table th,
        .data-table td {
            padding: 14px 16px;
            text-align: left;
            border-bottom: 1px solid #d0d0d0;
        }
        
        .data-table th {
            background: linear-gradient(135deg, var(--gradient-start) 0%, var(--gradient-end) 100%);
            color: white;
            font-weight: 600;
            text-transform: uppercase;
            letter-spacing: 0.5px;
            font-size: 0.9em;
        }
        
        .data-table tbody tr:hover {
            background-color: var(--bg-light);
            transition: background-color 0.2s;
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
        
        .error h2 {
            margin-bottom: 15px;
            font-size: 1.5em;
        }
        
        footer {
            margin-top: 80px;
            text-align: center;
            color: var(--text-muted);
            padding-top: 30px;
            border-top: 2px solid #ecf0f1;
            font-size: 0.95em;
        }
        
        footer p {
            margin: 8px 0;
        }
        
        @media (max-width: 768px) {
            .container {
                padding: 25px;
            }
            h1 {
                font-size: 2em;
            }
        }
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1>📊 Cohort Demographics</h1>
            <p class="subtitle">MS & MG Amino Acid Analysis - Table 1 & Table 2</p>
        </header>
        
        <div id="app">
            <div class="loading">
                <div class="spinner"></div>
                <p style="color: var(--text-muted); font-size: 1.1em;">Loading demographic data...</p>
            </div>
        </div>
        
        <footer>
            <p><strong>Dataset:</strong> Multiple Sclerosis & Myasthenia Gravis Proteomics</p>
            <p><strong>Generated:</strong> <span id="timestamp"></span></p>
            <p><strong>Node:</strong> 02 - Cohort Demographics</p>
        </footer>
    </div>
    
    <script>
        'use strict';
        
        function renderTable(containerId, data) {
            if (!data || !data.columns || !data.rows) {
                console.error('No table data for ' + containerId);
                return;
            }
            
            const container = document.getElementById(containerId);
            let html = '<table class="data-table"><thead><tr>';
            
            data.columns.forEach(col => {
                html += `<th>${col}</th>`;
            });
            html += '</tr></thead><tbody>';
            
            data.rows.forEach(row => {
                html += '<tr>';
                row.forEach(cell => {
                    html += `<td>${cell !== null && cell !== undefined ? cell : '-'}</td>`;
                });
                html += '</tr>';
            });
            html += '</tbody></table>';
            
            container.innerHTML = html;
        }
        
        function init() {
            try {
                const data = RAW_DATA;
                console.log('Data loaded successfully:', Object.keys(data));
                
                document.getElementById('app').innerHTML = `
                    <div class="figure-section">
                        <h2 class="figure-title">Table 1: MS vs Controls Demographics</h2>
                        <div id="table1"></div>
                    </div>
                    
                    <div class="figure-section">
                        <h2 class="figure-title">Table 2: MS and MG Demographic Comparison</h2>
                        <div id="table2"></div>
                    </div>
                `;
                
                document.getElementById('timestamp').textContent = new Date().toLocaleString();
                
                setTimeout(() => {
                    console.log('Rendering tables...');
                    renderTable('table1', data.table1);
                    renderTable('table2', data.table2);
                    console.log('Tables rendered successfully!');
                }, 150);
                
            } catch (error) {
                console.error('Error during initialization:', error);
                document.getElementById('app').innerHTML = `
                    <div class="error">
                        <h2>⚠️ Error Loading Data</h2>
                        <p><strong>Message:</strong> ${error.message}</p>
                        <p><strong>File:</strong> demographics_data.json</p>
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
