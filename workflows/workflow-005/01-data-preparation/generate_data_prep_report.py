import json
import os

def generate_report():
    print("Generating Report for Node 1...")
    if not os.path.exists("outputs/data.json"):
        print("Error: outputs/data.json not found.")
        return

    with open("outputs/data.json", "r") as f:
        data = json.load(f)

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Data Prep Report</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {{ font-family: sans-serif; padding: 20px; background: #f4f6f9; }}
            .card {{ background: white; padding: 20px; margin-bottom: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }}
            .grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(150px, 1fr)); gap: 10px; }}
            .mol-item {{ text-align: center; border: 1px solid #eee; padding: 5px; border-radius: 4px; }}
            img {{ max-width: 100%; height: auto; }}
        </style>
    </head>
    <body>
        <h1>Data Preparation Overview</h1>
        
        <div class="card">
            <h2>Summary</h2>
            <p><strong>Total Samples:</strong> {data['total_samples']}</p>
            <p><strong>Features Computed:</strong> {data['feature_count']}</p>
        </div>

        <div class="card">
            <h2>Feature Correlation (First 20 Features)</h2>
            <div id="heatmap"></div>
        </div>

        <div class="card">
            <h2>Sample Molecules</h2>
            <div class="grid">
                {''.join([f'<div class="mol-item"><img src="{item["img"]}"><p>{item["smiles"][:15]}...</p></div>' for item in data['sample_images']])}
            </div>
        </div>

        <script>
            var data = {json.dumps(data.get('correlation', {}))};
            if(data.z) {{
                var plotData = [{{
                    z: data.z,
                    x: data.x,
                    y: data.y,
                    type: 'heatmap',
                    colorscale: 'Viridis'
                }}];
                var layout = {{
                    width: 700,
                    height: 700,
                    margin: {{ l: 150, r: 50, t: 50, b: 150 }},
                    yaxis: {{ scaleanchor: 'x', scaleratio: 1 }}
                }};
                Plotly.newPlot('heatmap', plotData, layout);
            }}
        </script>
    </body>
    </html>
    """

    with open("outputs/report.html", "w") as f:
        f.write(html_content)
    print("Report generated at outputs/report.html")

if __name__ == "__main__":
    generate_report()