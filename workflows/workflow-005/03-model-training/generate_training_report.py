import json
import os

def generate_report():
    print("Generating Report for Node 3...")
    if not os.path.exists("outputs/data.json"):
        return

    with open("outputs/data.json", "r") as f:
        data = json.load(f)

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Training Report</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
        <style>body {{ font-family: sans-serif; padding: 20px; }} .plot {{ height: 400px; width: 100%; }}</style>
    </head>
    <body>
        <h1>Model Training Performance</h1>
        <h3>Test R² Score: {data['metrics']['r2_score']:.4f}</h3>
        
        <div id="loss_plot" class="plot"></div>
        <div id="parity_plot" class="plot"></div>

        <script>
            var hist = {json.dumps(data['history'])};
            var preds = {json.dumps(data['predictions'])};

            // Loss Plot
            Plotly.newPlot('loss_plot', [
                {{ y: hist.loss, name: 'Train Loss', type: 'scatter' }},
                {{ y: hist.val_loss, name: 'Val Loss', type: 'scatter' }}
            ], {{title: 'Loss over Epochs'}});

            // Parity Plot
            var trace = {{
                x: preds.y_test,
                y: preds.y_pred,
                mode: 'markers',
                type: 'scatter',
                marker: {{opacity: 0.6}}
            }};
            var line = {{
                x: [Math.min(...preds.y_test), Math.max(...preds.y_test)],
                y: [Math.min(...preds.y_test), Math.max(...preds.y_test)],
                mode: 'lines', name: 'Ideal', line: {{dash: 'dash', color: 'gray'}}
            }};
            Plotly.newPlot('parity_plot', [trace, line], {{
                title: 'Predicted vs Actual (Test Set)',
                xaxis: {{title: 'Actual'}}, yaxis: {{title: 'Predicted'}}
            }});
        </script>
    </body>
    </html>
    """
    with open("outputs/report.html", "w") as f:
        f.write(html_content)

if __name__ == "__main__":
    generate_report()