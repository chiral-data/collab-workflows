import json
import os

def generate_report():
    print("Generating Report for Node 4...")
    if not os.path.exists("outputs/data.json"):
        return

    with open("outputs/data.json", "r") as f:
        data = json.load(f)

    cards_html = ""
    for item in data:
        color = "#d4edda" if item['ad_status'] == "IN" else "#f8d7da"
        cards_html += f"""
        <div class="card">
            <img src="{item['image']}">
            <div class="container">
                <h4><b>{item['smiles'][:20]}...</b></h4> 
                <p>Prediction: <b>{item['prediction']:.4f}</b></p>
                <p style="background-color: {color}; padding: 5px; border-radius: 4px;">
                    Domain: {item['ad_status']}
                </p>
            </div>
        </div>
        """

    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Prediction Results</title>
        <style>
            body {{ font-family: sans-serif; padding: 20px; background: #f0f2f5; }}
            .grid {{ display: flex; flex-wrap: wrap; gap: 20px; justify-content: center; }}
            .card {{ 
                box-shadow: 0 4px 8px 0 rgba(0,0,0,0.2); 
                transition: 0.3s; 
                width: 250px; 
                background: white; 
                border-radius: 10px; 
                overflow: hidden;
                text-align: center;
            }}
            .card:hover {{ box-shadow: 0 8px 16px 0 rgba(0,0,0,0.2); }}
            .container {{ padding: 10px; }}
            img {{ width: 100%; height: 200px; object-fit: contain; background: white; }}
        </style>
    </head>
    <body>
        <h1 style="text-align:center;">Predicted Binding Affinities</h1>
        <div class="grid">
            {cards_html}
        </div>
    </body>
    </html>
    """
    with open("outputs/report.html", "w") as f:
        f.write(html_content)

if __name__ == "__main__":
    generate_report()