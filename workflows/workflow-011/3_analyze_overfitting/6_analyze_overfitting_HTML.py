import json
import os
import sys

def generate_report():
    print("Generating Interactive Dashboard for Node 3...")
    json_path = "outputs/data.json"
    if not os.path.exists(json_path):
        print("Error: outputs/data.json not found. Run 5_analyze_overfitting.py first.")
        sys.exit(1)

    with open(json_path, "r") as f:
        data = json.load(f)

    # Unpack data
    metrics = data.get("metrics", {})
    history = data.get("history", {})
    res_train = data.get("residuals", {}).get("train", {})
    res_test = data.get("residuals", {}).get("test", {})
    verdict = metrics.get("verdict", "UNKNOWN")

    # Define color based on verdict
    verdict_color = "#10b981" # Green
    if "MODERATE" in verdict: verdict_color = "#f59e0b" # Yellow
    if "SEVERE" in verdict: verdict_color = "#ef4444" # Red

    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <meta name="viewport" content="width=device-width, initial-scale=1.0">
        <title>Node 3: Overfitting Analysis</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.0.0/css/all.min.css" rel="stylesheet">
        <link href="https://fonts.googleapis.com/css2?family=Outfit:wght@300;400;500;600;700&display=swap" rel="stylesheet">
        <style>
            :root {{ 
                --primary: #4338ca;       /* Indigo 700 */
                --primary-light: #6366f1; /* Indigo 500 */
                --bg-body: #f8fafc;       /* Slate 50 */
                --bg-card: #ffffff;
                --text-main: #1e293b;     /* Slate 800 */
                --text-muted: #64748b;    /* Slate 500 */
                --border: #e2e8f0;
                --shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1), 0 2px 4px -1px rgba(0, 0, 0, 0.06);
            }}

            body {{ 
                font-family: 'Outfit', sans-serif; 
                background: var(--bg-body); 
                color: var(--text-main); 
                margin: 0; 
                padding: 0; 
                height: 100vh; 
                display: flex; 
                flex-direction: column; 
                overflow: hidden;
            }}
            
            /* TOP NAVIGATION BAR */
            .navbar {{ 
                background: var(--bg-card); 
                border-bottom: 1px solid var(--border); 
                padding: 0 2rem; 
                height: 64px; 
                display: flex; 
                align-items: center; 
                justify-content: space-between;
                box-shadow: 0 1px 2px 0 rgba(0, 0, 0, 0.05);
                z-index: 50;
            }}

            .brand {{ display: flex; align-items: center; gap: 1rem; }}
            .brand h1 {{ margin: 0; font-size: 1.25rem; font-weight: 700; color: var(--text-main); letter-spacing: -0.025em; }}
            .brand .badge {{ 
                padding: 0.25rem 0.75rem; 
                border-radius: 9999px; 
                font-size: 0.75rem; 
                font-weight: 600; 
                background: {verdict_color}15; 
                color: {verdict_color}; 
                border: 1px solid {verdict_color}40; 
            }}

            .nav-links {{ display: flex; height: 100%; gap: 0.5rem; }}
            .nav-item {{ 
                display: flex; 
                align-items: center; 
                gap: 0.5rem; 
                padding: 0 1.25rem; 
                height: 100%; 
                cursor: pointer; 
                color: var(--text-muted); 
                font-weight: 500; 
                font-size: 0.9rem;
                transition: all 0.2s; 
                border-bottom: 3px solid transparent;
            }}
            .nav-item:hover {{ color: var(--primary); background: #f1f5f9; }}
            .nav-item.active {{ color: var(--primary); border-bottom-color: var(--primary); font-weight: 600; }}
            
            /* MAIN LAYOUT */
            .main {{ 
                flex: 1; 
                padding: 2rem; 
                overflow-y: auto; 
                max-width: 1600px; 
                width: 100%; 
                margin: 0 auto; 
                box-sizing: border-box;
            }}

            /* FILTERS BAR */
            .filters-bar {{
                display: flex;
                align-items: center;
                gap: 1.5rem;
                margin-bottom: 1.5rem;
                padding: 1rem;
                background: var(--bg-card);
                border: 1px solid var(--border);
                border-radius: 0.75rem;
                box-shadow: var(--shadow);
            }}

            .filter-group {{ display: flex; align-items: center; gap: 0.75rem; }}
            .filter-label {{ font-size: 0.85rem; font-weight: 600; color: var(--text-muted); }}
            
            input[type=range] {{
                -webkit-appearance: none; width: 150px; height: 6px; background: #e2e8f0; border-radius: 5px; outline: none;
            }}
            input[type=range]::-webkit-slider-thumb {{
                -webkit-appearance: none; appearance: none; width: 18px; height: 18px; border-radius: 50%; background: var(--primary); cursor: pointer; border: 2px solid white; box-shadow: 0 1px 3px rgba(0,0,0,0.3);
            }}

            /* CARDS & GRID */
            .grid {{ display: grid; gap: 1.5rem; }}
            .grid-2 {{ grid-template-columns: repeat(2, 1fr); }}
            .grid-4 {{ grid-template-columns: repeat(4, 1fr); }}

            .card {{ 
                background: var(--bg-card); 
                border-radius: 1rem; 
                border: 1px solid var(--border); 
                padding: 1.5rem; 
                box-shadow: var(--shadow); 
                transition: transform 0.2s;
            }}
            .card:hover {{ transform: translateY(-2px); }}

            .stat-label {{ font-size: 0.85rem; font-weight: 600; text-transform: uppercase; letter-spacing: 0.05em; color: var(--text-muted); margin-bottom: 0.5rem; }}
            .stat-value {{ font-size: 1.75rem; font-weight: 700; color: var(--text-main); }}

            .plot-header {{ display: flex; flex-direction: column; justify-content: center; align-items: center; margin-bottom: 1rem; border-bottom: 1px solid #f1f5f9; padding-bottom: 0.75rem; text-align: center; gap: 0.25rem; }}
            .plot-title {{ font-size: 1.1rem; font-weight: 600; color: var(--text-main); margin: 0; }}

            .plot-container {{ height: 400px; width: 100%; border-radius: 0.5rem; overflow: hidden; }}

            /* TABS */
            .view {{ display: none; animation: fadeIn 0.4s cubic-bezier(0.16, 1, 0.3, 1); }}
            .view.active {{ display: block; }}
            @keyframes fadeIn {{ from {{ opacity: 0; transform: translateY(10px) scale(0.98); }} to {{ opacity: 1; transform: translateY(0) scale(1); }} }}
        </style>
    </head>
    <body>

        <!-- TOP NAVIGATION -->
        <div class="navbar">
            <div class="brand">
                <i class="fas fa-project-diagram" style="color: var(--primary); font-size: 1.5rem;"></i>
                <div>
                    <h1>Node 3 Analysis</h1>
                    <div style="font-size: 0.75rem; color: var(--text-muted);">Ultimate Hybrid Model Checkpoint</div>
                </div>
                <span class="badge">{verdict}</span>
            </div>
            
            <div class="nav-links">
                <div class="nav-item active" onclick="switchTab('overview', this)"><i class="fas fa-chart-pie"></i> Overview</div>
                <div class="nav-item" onclick="switchTab('learning', this)"><i class="fas fa-chart-line"></i> Learning Curves</div>
                <div class="nav-item" onclick="switchTab('residuals', this)"><i class="fas fa-poll"></i> Diagnostics</div>
                <div class="nav-item" onclick="switchTab('distribution', this)"><i class="fas fa-chart-bar"></i> Error Dist</div>
            </div>
        </div>

        <div class="main">
            
            <!-- VIEW 1: OVERVIEW -->
            <div id="overview" class="view active">
                <div class="grid grid-4" style="margin-bottom: 2rem;">
                    <div class="card">
                        <div class="stat-label">Train R² Score</div>
                        <div class="stat-value">{metrics.get('train_r2', 0):.4f}</div>
                    </div>
                    <div class="card">
                        <div class="stat-label">Test R² Score</div>
                        <div class="stat-value">{metrics.get('test_r2', 0):.4f}</div>
                    </div>
                    <div class="card">
                        <div class="stat-label">Generalization Gap</div>
                        <div class="stat-value" style="color: {verdict_color}">{metrics.get('r2_gap', 0):.2f}%</div>
                    </div>
                    <div class="card">
                        <div class="stat-label">Data Efficiency</div>
                        <div class="stat-value">{metrics.get('data_ratio', 0):.1f}<span style="font-size: 1rem; color: var(--text-muted);"> ratio</span></div>
                    </div>
                </div>

                <div class="grid grid-2">
                    <div class="card">
                        <div class="plot-header"><h3 class="plot-title">Performance Radar</h3></div>
                        <div id="radarPlot" class="plot-container" style="height: 350px;"></div>
                    </div>
                    <div class="card" style="background: linear-gradient(135deg, var(--bg-card) 0%, #f8fafc 100%);">
                        <div class="plot-header"><h3 class="plot-title">Verdict Analysis</h3></div>
                        <div style="height: 350px; display: flex; flex-direction: column; justify-content: center; align-items: center; text-align: center;">
                            <div style="background: {verdict_color}10; padding: 2rem; border-radius: 50%; margin-bottom: 1.5rem;">
                                <i class="fas fa-shield-alt" style="font-size: 4rem; color: {verdict_color};"></i>
                            </div>
                            <h2 style="margin: 0; color: {verdict_color}; font-size: 2rem;">{verdict}</h2>
                            <p style="color: var(--text-muted); max-width: 320px; margin-top: 1rem; line-height: 1.6;">
                                The model demonstrates a generalization gap of <strong>{metrics.get('r2_gap', 0):.2f}%</strong>.
                                <br>
                                { "Excellent stability across splits." if metrics.get('r2_gap', 0) < 5 else "Consider regularization or more data." }
                            </p>
                        </div>
                    </div>
                </div>
            </div>

            <!-- VIEW 2: LEARNING CURVES -->
            <div id="learning" class="view">
                <div class="card" style="margin-bottom: 1.5rem;">
                    <div class="plot-header">
                        <h3 class="plot-title">Loss Evolution (Log Scale)</h3>
                        <div style="font-size: 0.8rem; color: var(--text-muted);">Tracks convergence over {len(history.get('loss', []))} epochs</div>
                    </div>
                    <div id="lossPlot" class="plot-container"></div>
                </div>
                <div class="grid grid-2">
                    <div class="card">
                        <div class="plot-header"><h3 class="plot-title">R² Score Evolution</h3></div>
                        <div id="r2Plot" class="plot-container"></div>
                    </div>
                    <div class="card">
                        <div class="plot-header"><h3 class="plot-title">RMSE Evolution</h3></div>
                        <div id="rmsePlot" class="plot-container"></div>
                    </div>
                </div>
            </div>

            <!-- VIEW 3: DIAGNOSTICS -->
            <div id="residuals" class="view">
                
                <!-- INTERACTIVE FILTERS -->
                <div class="filters-bar">
                    <div class="filter-group">
                        <i class="fas fa-filter" style="color: var(--primary);"></i>
                        <span class="filter-label">Highlight High Error:</span>
                    </div>
                    <div class="filter-group">
                        <span style="font-size: 0.8rem;">Min Residual:</span>
                        <input type="range" id="residSlider" min="0" max="2" step="0.05" value="0" oninput="updateDiagnostics(this.value)">
                        <span id="residVal" style="font-weight: 600; color: var(--primary); width: 40px;">0.00</span>
                    </div>
                    <div style="flex: 1;"></div>
                    <div style="font-size: 0.8rem; color: var(--text-muted);">Adjust slider to isolate outliers</div>
                </div>

                <div class="grid grid-2">
                    <div class="card">
                        <div class="plot-header"><h3 class="plot-title">Predicted vs Actual</h3></div>
                        <div id="parityPlot" class="plot-container"></div>
                    </div>
                    <div class="card">
                        <div class="plot-header"><h3 class="plot-title">Residuals vs Predicted</h3></div>
                        <div id="residPlot" class="plot-container"></div>
                    </div>
                </div>
            </div>

            <!-- VIEW 4: DISTRIBUTION -->
            <div id="distribution" class="view">
                <div class="card">
                    <div class="plot-header"><h3 class="plot-title">Error Distribution (Residuals)</h3></div>
                    <div id="histPlot" class="plot-container"></div>
                </div>
            </div>

        </div>

        <script>
            // DATA INJECTION
            const history = {json.dumps(history)};
            const resTrain = {json.dumps(res_train)};
            const resTest = {json.dumps(res_test)};
            const metrics = {json.dumps(metrics)};
            const verdictColor = '{verdict_color}';

            // TABS LOGIC
            function switchTab(viewId, navElement) {{
                document.querySelectorAll('.view').forEach(el => el.classList.remove('active'));
                document.getElementById(viewId).classList.add('active');
                document.querySelectorAll('.nav-item').forEach(el => el.classList.remove('active'));
                navElement.classList.add('active');
                
                // Force Plotly resize after layout update
                setTimeout(() => {{
                    window.dispatchEvent(new Event('resize'));
                    const activeView = document.getElementById(viewId);
                    const plots = activeView.querySelectorAll('.plot-container');
                    plots.forEach(div => {{
                        if (div.data) Plotly.Plots.resize(div);
                    }});
                }}, 10); 
            }}

            // ------------------------
            // PLOT GENERATION
            // ------------------------
            
            // Common Layout Config
            const layoutConfig = {{
                font: {{ family: 'Outfit, sans-serif' }},
                plot_bgcolor: 'rgba(0,0,0,0)',
                paper_bgcolor: 'rgba(0,0,0,0)',
                margin: {{ t: 20, b: 40, l: 60, r: 20 }},
                xaxis: {{ showgrid: true, gridcolor: '#f1f5f9' }},
                yaxis: {{ showgrid: true, gridcolor: '#f1f5f9' }}
            }};

            // 1. RADAR PLOT
            Plotly.newPlot('radarPlot', [{{
                type: 'scatterpolar',
                r: [metrics.train_r2, 1 - (metrics.r2_gap/100), metrics.test_r2, metrics.test_rmse < 1 ? 1 : 0.5],
                theta: ['Train R²', 'Generalization', 'Test R²', 'RMSE Stability'],
                fill: 'toself',
                line: {{color: verdictColor}},
                fillcolor: verdictColor + '33'
            }}], {{
                ...layoutConfig,
                polar: {{ radialaxis: {{ visible: true, range: [0, 1] }} }},
                showlegend: false
            }});

            // 2. LOSS PLOT
            const epochs = history.loss.map((_, i) => i + 1);
            Plotly.newPlot('lossPlot', [
                {{ x: epochs, y: history.loss, name: 'Train Loss', line: {{width: 3, color: '#00008b'}} }},
                {{ x: epochs, y: history.val_loss, name: 'Val Loss', line: {{width: 3, dash: 'dash', color: '#f59e0b'}} }}
            ], {{
                ...layoutConfig,
                yaxis: {{ type: 'log', title: 'Loss (Log)', gridcolor: '#f1f5f9' }},
                xaxis: {{ title: 'Epoch' }},
                legend: {{ orientation: 'h', x: 0.5, xanchor: 'center', y: 1.1 }}
            }});

            // 3. R2 & RMSE
            Plotly.newPlot('r2Plot', [
                {{ x: epochs, y: history.train_r2_sklearn, name: 'Train R²', line: {{color: '#00008b'}} }},
                {{ x: epochs, y: history.val_r2_sklearn, name: 'Val R²', line: {{color: '#f59e0b'}} }}
            ], {{ ...layoutConfig, legend: {{ orientation: 'h', x: 0.5, xanchor: 'center', y: 1.1 }} }});

            Plotly.newPlot('rmsePlot', [
                {{ x: epochs, y: history.train_rmse_sklearn, name: 'Train RMSE', line: {{color: '#00008b'}} }},
                {{ x: epochs, y: history.val_rmse_sklearn, name: 'Val RMSE', line: {{color: '#f59e0b'}} }}
            ], {{ ...layoutConfig, legend: {{ orientation: 'h', x: 0.5, xanchor: 'center', y: 1.1 }} }});

            // 4. PREPARE DIAGNOSTICS DATA
            function drawDiagnostics(threshold) {{
                // Filter points based on absolute residual > threshold
                const trainMask = resTrain.residual.map(r => Math.abs(r) >= threshold);
                const testMask = resTest.residual.map(r => Math.abs(r) >= threshold);
                
                // Helper to filter arrays
                const filter = (arr, mask) => arr.filter((_, i) => mask[i]);

                // Parity Plot
                Plotly.react('parityPlot', [
                    {{ 
                        x: filter(resTrain.actual, trainMask), 
                        y: filter(resTrain.predicted, trainMask), 
                        mode: 'markers', name: 'Train', 
                        marker: {{opacity: 0.5, size: 5, color: '#00008b'}} 
                    }},
                    {{ 
                        x: filter(resTest.actual, testMask), 
                        y: filter(resTest.predicted, testMask), 
                        mode: 'markers', name: 'Test', 
                        marker: {{opacity: 0.7, size: 6, color: '#ef4444'}} 
                    }},
                    {{ 
                        x: [Math.min(...resTrain.actual), Math.max(...resTrain.actual)], 
                        y: [Math.min(...resTrain.actual), Math.max(...resTrain.actual)], 
                        mode: 'lines', name: 'Ideal', 
                        line: {{color: '#ef4444', dash: 'dash'}} 
                    }}
                ], {{ ...layoutConfig, xaxis: {{ title: 'Actual' }}, yaxis: {{ title: 'Predicted' }}, legend: {{ orientation: 'h', x: 0.5, xanchor: 'center', y: 1.1 }} }});

                // Residual Plot
                Plotly.react('residPlot', [
                    {{ 
                        x: filter(resTrain.predicted, trainMask), 
                        y: filter(resTrain.residual, trainMask), 
                        mode: 'markers', name: 'Train', 
                        marker: {{opacity: 0.5, size: 5, color: '#00008b'}} 
                    }},
                    {{ 
                        x: filter(resTest.predicted, testMask), 
                        y: filter(resTest.residual, testMask), 
                        mode: 'markers', name: 'Test', 
                        marker: {{opacity: 0.7, size: 6, color: '#ef4444', symbol: 'diamond'}} 
                    }}
                ], {{ 
                    ...layoutConfig, 
                    xaxis: {{ title: 'Predicted' }}, yaxis: {{ title: 'Residuals' }},
                    shapes: [{{ type: 'line', x0: Math.min(...resTrain.predicted), x1: Math.max(...resTrain.predicted), y0: 0, y1: 0, line: {{color: '#ef4444', width: 2}} }}],
                    legend: {{ orientation: 'h', x: 0.5, xanchor: 'center', y: 1.1 }}
                }});
            }}

            // Initial Draw
            drawDiagnostics(0);

            // Update function for Slider
            window.updateDiagnostics = function(val) {{
                document.getElementById('residVal').innerText = parseFloat(val).toFixed(2);
                drawDiagnostics(parseFloat(val));
            }}

            // 5. HISTOGRAM
            Plotly.newPlot('histPlot', [
                {{ x: resTrain.residual, type: 'histogram', name: 'Train Error', opacity: 0.6, marker: {{color: '#00008b'}} }},
                {{ x: resTest.residual, type: 'histogram', name: 'Test Error', opacity: 0.6, marker: {{color: '#ef4444'}} }}
            ], {{ ...layoutConfig, barmode: 'overlay', title: 'Residual Distribution', legend: {{ orientation: 'h', x: 0.5, xanchor: 'center', y: 1.1 }} }});

        </script>
    </body>
    </html>
    """
    
    with open("outputs/report.html", "w") as f:
        f.write(html_content)
    
    print(f"Interactive Dashboard generated: outputs/report.html")

if __name__ == "__main__":
    generate_report()