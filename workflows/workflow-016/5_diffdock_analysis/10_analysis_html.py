#!/usr/bin/env python3
"""
Node 5: HTML Report Generation for Interface Analysis

Generates an interactive HTML dashboard with Plotly visualizations
showing interface contacts, interaction types, and residue details.
"""

import json
import argparse
import os
from pathlib import Path
from collections import defaultdict


def load_analysis_data(data_json_path):
    """Load the analysis data from JSON."""
    with open(data_json_path) as f:
        return json.load(f)


def load_contact_residues(contact_json_path):
    """Load contact residues JSON."""
    if not os.path.exists(contact_json_path):
        return {'receptor_interface_residues': [], 'ligand_interface_residues': []}
    
    with open(contact_json_path) as f:
        return json.load(f)


def load_interface_analysis_txt(txt_path):
    """Parse the interface analysis text file."""
    if not os.path.exists(txt_path):
        return ""
    
    with open(txt_path) as f:
        return f.read()


def create_html_report(data, contact_residues, interface_txt, output_html_path):
    """
    Generate an interactive HTML report with Plotly visualizations.
    """
    
    # Default empty data if not available
    if data.get('status') != 'success':
        interaction_types = {}
        total_contacts = 0
        hbond_count = 0
    else:
        interaction_types = data.get('interaction_types', {})
        total_contacts = data.get('total_contacts', 0)
        hbond_count = data.get('hydrogen_bonds', 0)
    
    receptor_residues = contact_residues.get('receptor_interface_residues', [])
    ligand_residues = contact_residues.get('ligand_interface_residues', [])
    
    # Prepare data for Plotly
    interaction_names = list(interaction_types.keys())
    interaction_values = list(interaction_types.values())
    
    # Color mapping for interaction types
    color_map = {
        'hydrophobic': '#FF6B6B',
        'electrostatic': '#4ECDC4',
        'hydrogen_bond': '#45B7D1',
        'polar': '#FFA07A',
        'other': '#95E1D3'
    }
    colors = [color_map.get(itype, '#CCCCCC') for itype in interaction_names]
    
    # Create Plotly script for pie chart
    pie_chart_data = []
    for name, value, color in zip(interaction_names, interaction_values, colors):
        pie_chart_data.append(f"{{x: '{name}', y: {value}, marker: {{color: '{color}'}}}}")
    
    pie_chart_json = "[" + ",".join(pie_chart_data) + "]"
    
    # Create residue tables HTML
    receptor_table_html = _generate_residue_table(receptor_residues, "Antibody (Receptor)")
    ligand_table_html = _generate_residue_table(ligand_residues, "Antigen (Ligand)")
    
    # Format summary statistics
    unique_receptor = data.get('unique_receptor_residues', len(receptor_residues))
    unique_ligand = data.get('unique_ligand_residues', len(ligand_residues))
    
    # Escape text for HTML
    interface_txt_escaped = interface_txt.replace('<', '&lt;').replace('>', '&gt;').replace('\n', '<br>')
    
    # Build HTML
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Interface Analysis Report</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: #333;
            line-height: 1.6;
            padding: 20px;
        }}
        
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 10px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }}
        
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        
        .header p {{
            font-size: 1.1em;
            opacity: 0.9;
        }}
        
        .content {{
            padding: 40px;
        }}
        
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .summary-card {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 20px;
            border-radius: 8px;
            border-left: 5px solid #667eea;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        .summary-card h3 {{
            color: #667eea;
            margin-bottom: 10px;
            font-size: 0.9em;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        
        .summary-card .value {{
            font-size: 2.2em;
            font-weight: bold;
            color: #333;
        }}
        
        .section {{
            margin-bottom: 40px;
            border-bottom: 2px solid #eee;
            padding-bottom: 30px;
        }}
        
        .section:last-child {{
            border-bottom: none;
        }}
        
        .section h2 {{
            color: #667eea;
            margin-bottom: 20px;
            font-size: 1.8em;
            display: flex;
            align-items: center;
        }}
        
        .section h2 span {{
            display: inline-block;
            width: 40px;
            height: 40px;
            background: #667eea;
            color: white;
            border-radius: 50%;
            text-align: center;
            line-height: 40px;
            margin-right: 15px;
            font-size: 0.8em;
        }}
        
        .chart-container {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }}
        
        .chart {{
            background: white;
            border: 1px solid #eee;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        .pie-chart {{
            width: 100%;
            height: 400px;
        }}
        
        .table-container {{
            overflow-x: auto;
            margin-bottom: 20px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            background: white;
            box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        }}
        
        table thead {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
        }}
        
        table th {{
            padding: 15px;
            text-align: left;
            font-weight: 600;
            letter-spacing: 0.5px;
        }}
        
        table td {{
            padding: 12px 15px;
            border-bottom: 1px solid #eee;
        }}
        
        table tbody tr:hover {{
            background: #f9f9f9;
            transition: background 0.2s ease;
        }}
        
        table tbody tr:nth-child(even) {{
            background: #f5f5f5;
        }}
        
        .residue-list {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(150px, 1fr));
            gap: 10px;
            margin-bottom: 20px;
        }}
        
        .residue-badge {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 10px 15px;
            border-radius: 5px;
            text-align: center;
            font-weight: 600;
            box-shadow: 0 2px 5px rgba(0,0,0,0.1);
        }}
        
        .interaction-badge {{
            display: inline-block;
            padding: 5px 10px;
            border-radius: 3px;
            font-size: 0.85em;
            font-weight: 600;
            color: white;
            margin: 2px;
        }}
        
        .hydrophobic {{
            background: #FF6B6B;
        }}
        
        .electrostatic {{
            background: #4ECDC4;
        }}
        
        .hydrogen_bond {{
            background: #45B7D1;
        }}
        
        .polar {{
            background: #FFA07A;
        }}
        
        .other {{
            background: #95E1D3;
        }}
        
        .footer {{
            background: #f5f5f5;
            padding: 20px 40px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
            border-top: 1px solid #eee;
        }}
        
        pre {{
            background: #f5f5f5;
            border: 1px solid #ddd;
            border-radius: 5px;
            padding: 15px;
            overflow-x: auto;
            font-size: 0.85em;
            line-height: 1.4;
        }}
        
        .stats-bar {{
            display: flex;
            gap: 10px;
            flex-wrap: wrap;
            margin-bottom: 20px;
        }}
        
        .stat-item {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 10px 15px;
            border-radius: 5px;
            font-weight: 600;
        }}
        
        @media (max-width: 768px) {{
            .header {{
                padding: 20px;
            }}
            
            .header h1 {{
                font-size: 1.8em;
            }}
            
            .content {{
                padding: 20px;
            }}
            
            .chart-container {{
                grid-template-columns: 1fr;
            }}
            
            .summary-grid {{
                grid-template-columns: 1fr;
            }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🔬 Interface Analysis Report</h1>
            <p>Antibody-Antigen Binding Interface Characterization</p>
        </div>
        
        <div class="content">
            <!-- Summary Statistics -->
            <div class="summary-grid">
                <div class="summary-card">
                    <h3>Total Contacts</h3>
                    <div class="value">{total_contacts}</div>
                </div>
                <div class="summary-card">
                    <h3>Hydrogen Bonds</h3>
                    <div class="value">{hbond_count}</div>
                </div>
                <div class="summary-card">
                    <h3>Antibody Interface Residues</h3>
                    <div class="value">{unique_receptor}</div>
                </div>
                <div class="summary-card">
                    <h3>Antigen Interface Residues</h3>
                    <div class="value">{unique_ligand}</div>
                </div>
            </div>
            
            <!-- Interaction Type Distribution -->
            <div class="section">
                <h2><span>1</span>Interaction Type Distribution</h2>
                <div class="chart-container">
                    <div class="chart">
                        <div id="pie-chart" class="pie-chart"></div>
                    </div>
                </div>
                <div class="stats-bar">
                    {_generate_interaction_badges(interaction_types)}
                </div>
            </div>
            
            <!-- Interface Residues -->
            <div class="section">
                <h2><span>2</span>Interface Residues</h2>
                
                <h3 style="margin-top: 25px; margin-bottom: 15px; color: #555;">Antibody (Receptor) Interface Residues</h3>
                {receptor_table_html}
                
                <h3 style="margin-top: 25px; margin-bottom: 15px; color: #555;">Antigen (Ligand) Interface Residues</h3>
                {ligand_table_html}
            </div>
            
            <!-- Detailed Analysis -->
            <div class="section">
                <h2><span>3</span>Detailed Analysis</h2>
                <pre style="max-height: 500px; overflow-y: auto;">{interface_txt_escaped}</pre>
            </div>
        </div>
        
        <div class="footer">
            <p>Generated by DiffDock_AbAg_Pipeline Node 5 - Interface Analysis</p>
            <p>For details, see interface_analysis.txt and contact_residues.json</p>
        </div>
    </div>
    
    <script>
        // Pie chart data
        var data = {pie_chart_json};
        var layout = {{
            title: 'Interaction Types',
            font: {{
                family: "'Segoe UI', Tahoma, Geneva, Verdana, sans-serif",
                size: 12
            }},
            showlegend: true,
            margin: {{
                l: 50,
                r: 50,
                t: 50,
                b: 50
            }},
            paper_bgcolor: 'white',
            plot_bgcolor: 'white'
        }};
        
        Plotly.newPlot('pie-chart', data, layout, {{responsive: true}});
    </script>
</body>
</html>"""
    
    with open(output_html_path, 'w') as f:
        f.write(html_content)


def _generate_residue_table(residues, title):
    """Generate an HTML table for residues."""
    if not residues:
        return f"<p><em>No {title} interface residues found.</em></p>"
    
    # Group by chain
    by_chain = defaultdict(list)
    for res in residues:
        by_chain[res['chain']].append(res)
    
    html = f'<div class="table-container"><table><thead><tr><th>Chain</th><th>Residue #</th><th>Residue Name</th></tr></thead><tbody>'
    
    for chain_id in sorted(by_chain.keys()):
        residues_in_chain = by_chain[chain_id]
        for res in sorted(residues_in_chain, key=lambda r: r['residue_number']):
            html += f"<tr><td><strong>{res['chain']}</strong></td><td>{res['residue_number']}</td><td>{res['residue_name']}</td></tr>"
    
    html += '</tbody></table></div>'
    return html


def _generate_interaction_badges(interaction_types):
    """Generate HTML badges for interaction types."""
    badges = []
    for itype, count in sorted(interaction_types.items(), key=lambda x: x[1], reverse=True):
        class_name = itype.replace('_', ' ').lower()
        badge = f'<div class="stat-item"><span class="interaction-badge {class_name}">{itype.upper()}</span> {count}</div>'
        badges.append(badge)
    return ''.join(badges)


def main():
    parser = argparse.ArgumentParser(
        description='Generate interactive HTML report for interface analysis'
    )
    parser.add_argument(
        '--data_json',
        default='5_diffdock_analysis/outputs/data.json',
        help='Path to data.json from analysis script'
    )
    parser.add_argument(
        '--contact_residues_json',
        default=None,
        help='Path to contact_residues.json'
    )
    parser.add_argument(
        '--interface_analysis_txt',
        default=None,
        help='Path to interface_analysis.txt'
    )
    parser.add_argument(
        '--output_html',
        default='5_diffdock_analysis/outputs/report.html',
        help='Output HTML file path'
    )
    
    args = parser.parse_args()
    
    # Determine default paths if not provided
    output_dir = os.path.dirname(args.output_html) or '5_diffdock_analysis/outputs'
    
    if not args.contact_residues_json:
        args.contact_residues_json = os.path.join(output_dir, 'contact_residues.json')
    
    if not args.interface_analysis_txt:
        args.interface_analysis_txt = os.path.join(output_dir, 'interface_analysis.txt')
    
    # Load data
    data = load_analysis_data(args.data_json)
    contact_residues = load_contact_residues(args.contact_residues_json)
    interface_txt = load_interface_analysis_txt(args.interface_analysis_txt)
    
    # Generate HTML
    create_html_report(data, contact_residues, interface_txt, args.output_html)
    print(f"Generated report: {args.output_html}")


if __name__ == '__main__':
    main()
