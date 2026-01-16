import json
import os

def generate_report():
    print("Generating Report for Node 2...")
    if not os.path.exists("outputs/data.json"):
        return

    with open("outputs/data.json", "r") as f:
        data = json.load(f)

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Feature Engineering Report</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>
            body {{ font-family: sans-serif; padding: 20px; }}
            .card {{ background: white; padding: 20px; box-shadow: 0 0 10px rgba(0,0,0,0.1); margin-bottom: 20px; }}
        </style>
    </head>
    <body>
        <h1>Feature Engineering Analysis</h1>
        <div class="card">
            <h3>Outlier Removal Stats</h3>
            <p>Final Training Samples: {data['train_samples']}</p>
            <p>Final Test Samples: {data['test_samples']}</p>
            <p>Outliers Removed from Train: {data['removed_train']}</p>
        </div>
        <div class="card">
            <h3>PCA Space Distribution (Training Set)</h3>
            <div id="pca_plot"></div>
        </div>

        <script>
            var pca = {json.dumps(data['pca_data'])};
            
            var traceIn = {{
                x: pca.train_x.filter((_, i) => pca.train_label[i] === 'Inlier'),
                y: pca.train_y.filter((_, i) => pca.train_label[i] === 'Inlier'),
                mode: 'markers',
                type: 'scatter',
                name: 'Inliers',
                marker: {{color: 'blue', size: 6}}
            }};
            
            var traceOut = {{
                x: pca.train_x.filter((_, i) => pca.train_label[i] === 'Outlier'),
                y: pca.train_y.filter((_, i) => pca.train_label[i] === 'Outlier'),
                mode: 'markers',
                type: 'scatter',
                name: 'Outliers',
                marker: {{color: 'red', symbol: 'x', size: 8}}
            }};

            Plotly.newPlot('pca_plot', [traceIn, traceOut], {{title: 'PCA Outlier Detection'}});
        </script>
    </body>
    </html>
    """
    with open("outputs/report.html", "w") as f:
        f.write(html_content)

if __name__ == "__main__":
    generate_report()