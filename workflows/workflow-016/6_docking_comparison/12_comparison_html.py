#!/usr/bin/env python3
"""
Node 6 HTML Report Generator: Comparison with Wet‑Lab Structure

Produces an interactive HTML dashboard summarizing docking accuracy assessment,
quality metrics, and visualization links for structure comparison.
"""

import json, argparse, os
from datetime import datetime

def get_quality_color(quality):
    """Return HTML color code for quality category."""
    colors = {
        'excellent': '#28a745',  # Green
        'good': '#17a2b8',        # Teal
        'moderate': '#ffc107',    # Orange/Yellow
        'poor': '#dc3545',        # Red
        'unknown': '#6c757d'      # Gray
    }
    return colors.get(quality, '#6c757d')

def get_quality_icon(quality):
    """Return emoji indicator for quality category."""
    icons = {
        'excellent': '✓✓',
        'good': '✓',
        'moderate': '⚠',
        'poor': '✗',
        'unknown': '?'
    }
    return icons.get(quality, '?')

def format_rmsd_metrics(data):
    """Format RMSD metrics as HTML table rows."""
    html = ""
    if 'full_rmsd' in data and data['full_rmsd'] is not None:
        html += f"""
    <tr>
        <td><strong>Full Antigen RMSD</strong></td>
        <td><code>{data['full_rmsd']:.2f}</code> Å</td>
        <td>All Cα atoms of antigen after antibody alignment</td>
    </tr>
"""
    
    if 'interface_rmsd' in data and data.get('interface_rmsd'):
        html += f"""
    <tr>
        <td><strong>Interface RMSD</strong></td>
        <td><code>{data['interface_rmsd']:.2f}</code> Å</td>
        <td>Only contact residues at antibody‑antigen interface</td>
    </tr>
"""
    
    if 'num_atoms_compared' in data:
        html += f"""
    <tr>
        <td><strong>Atoms Compared</strong></td>
        <td><code>{data['num_atoms_compared']}</code></td>
        <td>Number of Cα atoms in antigen</td>
    </tr>
"""
    
    return html

def create_html_report(data):
    """Generate comprehensive HTML report."""
    
    # Determine report status
    if not data.get('comparison_performed', False):
        return create_error_report(data)
    
    quality = data.get('quality', 'unknown')
    quality_color = get_quality_color(quality)
    quality_icon = get_quality_icon(quality)
    timestamp = data.get('timestamp', 'Unknown')
    full_rmsd = data.get('full_rmsd', 'N/A')
    
    # Format RMSD metrics table
    metrics_html = format_rmsd_metrics(data)
    
    # Prepare PyMOL script link
    pymol_link = ""
    if data.get('pymol_script'):
        pymol_link = f"""
    <div class="card mt-4">
        <div class="card-header bg-info text-white">
            <h5 class="mb-0">🔬 PyMOL Alignment Script</h5>
        </div>
        <div class="card-body">
            <p>A PyMOL script has been generated for interactive visualization and alignment:</p>
            <code>{data['pymol_script']}</code>
            <p class="mt-3"><strong>Usage:</strong></p>
            <pre><code>pymol {data['pymol_script']}</code></pre>
            <p class="text-muted small">The script loads both structures, aligns the antibodies, and applies surface representation to visualize the interface.</p>
        </div>
    </div>
"""
    
    # Build interpretation details
    interpretation = data.get('interpretation', 'Assessment not available')
    quality_description = {
        'excellent': 'The predicted pose shows excellent spatial agreement with the experimental structure. The docking model has successfully captured the correct binding mode.',
        'good': 'The predicted pose shows good agreement with the experimental structure. Some conformational variations are visible, but key interactions are well‑predicted.',
        'moderate': 'The predicted pose shows moderate agreement with the experimental structure. Significant conformational differences exist; some aspects of binding may need refinement.',
        'poor': 'The predicted pose shows poor agreement with the experimental structure. Consider re‑running with different parameters or verifying input data.'
    }
    quality_details = quality_description.get(quality, 'Assessment details not available.')
    
    # Chain information
    chain_info = ""
    if data.get('antibody_chains') and data.get('antigen_chains'):
        chain_info = f"""
    <tr>
        <td><strong>Antibody chains</strong></td>
        <td><code>{', '.join(data['antibody_chains'])}</code></td>
    </tr>
    <tr>
        <td><strong>Antigen chains</strong></td>
        <td><code>{', '.join(data['antigen_chains'])}</code></td>
    </tr>
"""
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Antibody‑Antigen Complex Comparison</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        :root {{
            --primary-color: #2c3e50;
            --secondary-color: #3498db;
            --success-color: #27ae60;
            --warning-color: #f39c12;
            --danger-color: #e74c3c;
        }}
        
        body {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #2c3e50;
        }}
        
        .header {{
            background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
            color: white;
            padding: 2rem 0;
            margin-bottom: 2rem;
            box-shadow: 0 4px 6px rgba(0, 0, 0, 0.1);
        }}
        
        .header h1 {{
            margin-bottom: 0.5rem;
            font-weight: 700;
            letter-spacing: 0.5px;
        }}
        
        .header-subtitle {{
            font-size: 0.95rem;
            opacity: 0.95;
        }}
        
        .container {{
            max-width: 1200px;
        }}
        
        .card {{
            border: none;
            box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1);
            border-radius: 0.5rem;
            transition: transform 0.2s, box-shadow 0.2s;
        }}
        
        .card:hover {{
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15);
        }}
        
        .quality-badge {{
            display: inline-block;
            padding: 0.5rem 1rem;
            border-radius: 0.25rem;
            font-weight: 600;
            font-size: 1.1rem;
            background-color: {quality_color};
            color: white;
        }}
        
        .rmsd-value {{
            font-size: 2.5rem;
            font-weight: 700;
            color: {quality_color};
            font-family: 'Courier New', monospace;
        }}
        
        .rmsd-unit {{
            font-size: 1.2rem;
            color: #666;
            margin-left: 0.3rem;
        }}
        
        .metric-box {{
            background: #f8f9fa;
            padding: 1.5rem;
            border-radius: 0.5rem;
            border-left: 4px solid {quality_color};
        }}
        
        .table-striped tbody tr:hover {{
            background-color: #f8f9fa;
        }}
        
        .table td {{
            vertical-align: middle;
            padding: 1rem 0.75rem;
        }}
        
        .code-block {{
            background-color: #2c3e50;
            color: #ecf0f1;
            padding: 1rem;
            border-radius: 0.5rem;
            font-family: 'Courier New', monospace;
            font-size: 0.9rem;
            overflow-x: auto;
        }}
        
        .interpretation-text {{
            font-size: 1.05rem;
            line-height: 1.6;
            padding: 1.5rem;
            background-color: #f0f8ff;
            border-left: 4px solid {quality_color};
            border-radius: 0.25rem;
            color: #2c3e50;
        }}
        
        .section-header {{
            border-bottom: 3px solid {quality_color};
            padding-bottom: 0.75rem;
            margin-top: 2rem;
            margin-bottom: 1.5rem;
        }}
        
        .icon-column {{
            font-size: 1.5rem;
            font-weight: 700;
            color: {quality_color};
        }}
        
        .info-badge {{
            background-color: #d4edff;
            border-left: 4px solid #17a2b8;
            padding: 0.75rem 1rem;
            margin-bottom: 1rem;
            border-radius: 0.25rem;
        }}
        
        .timestamp {{
            color: #7f8c8d;
            font-size: 0.9rem;
            font-style: italic;
        }}
        
        .footer {{
            margin-top: 3rem;
            padding-top: 2rem;
            border-top: 1px solid #bdc3c7;
            text-align: center;
            color: #7f8c8d;
            font-size: 0.9rem;
        }}
        
        @media (max-width: 768px) {{
            .rmsd-value {{
                font-size: 2rem;
            }}
            .container {{
                padding: 1rem;
            }}
        }}
    </style>
</head>
<body>
    <div class="header">
        <div class="container">
            <h1>🧬 Antibody‑Antigen Complex Comparison</h1>
            <p class="header-subtitle">Structural validation of DiffDock‑PP predictions against experimental data</p>
        </div>
    </div>
    
    <div class="container mb-5">
        <!-- Quality Assessment Card -->
        <div class="card mb-4">
            <div class="card-header" style="background-color: {quality_color}; color: white; border-bottom: none;">
                <h5 class="mb-0">🎯 Quality Assessment</h5>
            </div>
            <div class="card-body">
                <div class="row align-items-center">
                    <div class="col-md-6">
                        <span class="quality-badge">{quality_icon} {quality.upper()}</span>
                        <p class="mt-3 interpretation-text">
                            {quality_details}
                        </p>
                    </div>
                    <div class="col-md-6 text-center metric-box">
                        <p class="text-muted mb-2 small">Full Antigen RMSD</p>
                        <div class="rmsd-value">{full_rmsd:.2f}<span class="rmsd-unit">Å</span></div>
                        <p class="text-muted small mt-2">{interpretation}</p>
                    </div>
                </div>
            </div>
        </div>
        
        <!-- RMSD Interpretation Guide -->
        <div class="card mb-4">
            <div class="card-header bg-secondary text-white">
                <h5 class="mb-0">📊 RMSD Interpretation Scale</h5>
            </div>
            <div class="card-body">
                <table class="table table-striped mb-0">
                    <tbody>
                        <tr style="border-left: 4px solid #28a745;">
                            <td class="icon-column" style="color: #28a745;">✓✓</td>
                            <td><strong>Excellent</strong></td>
                            <td>RMSD &lt; 2.0 Å</td>
                            <td>Exceptional agreement with experimental structure</td>
                        </tr>
                        <tr style="border-left: 4px solid #17a2b8;">
                            <td class="icon-column" style="color: #17a2b8;">✓</td>
                            <td><strong>Good</strong></td>
                            <td>RMSD 2.0–5.0 Å</td>
                            <td>Reliable prediction with acceptable differences</td>
                        </tr>
                        <tr style="border-left: 4px solid #ffc107;">
                            <td class="icon-column" style="color: #ffc107;">⚠</td>
                            <td><strong>Moderate</strong></td>
                            <td>RMSD 5.0–10.0 Å</td>
                            <td>Significant structural differences; refinement needed</td>
                        </tr>
                        <tr style="border-left: 4px solid #dc3545;">
                            <td class="icon-column" style="color: #dc3545;">✗</td>
                            <td><strong>Poor</strong></td>
                            <td>RMSD &gt; 10.0 Å</td>
                            <td>Poor agreement; model may have failed</td>
                        </tr>
                    </tbody>
                </table>
            </div>
        </div>
        
        <!-- Detailed Metrics -->
        <div class="card mb-4">
            <div class="card-header bg-info text-white">
                <h5 class="mb-0">📈 Detailed RMSD Metrics</h5>
            </div>
            <div class="card-body">
                <table class="table table-hover mb-0">
                    <tbody>
                        {metrics_html}
                        {chain_info}
                        <tr>
                            <td><strong>Reference complex</strong></td>
                            <td><code>{data.get('reference_complex', 'N/A')}</code></td>
                        </tr>
                        <tr>
                            <td><strong>Predicted pose</strong></td>
                            <td><code>{data.get('predicted_pose', 'N/A')}</code></td>
                        </tr>
                    </tbody>
                </table>
            </div>
        </div>
        
        <!-- PyMOL Script Section -->
        {pymol_link}
        
        <!-- Residue-level Analysis Info -->
        <div class="card mb-4">
            <div class="card-header bg-warning text-dark">
                <h5 class="mb-0">🔍 Interface Analysis</h5>
            </div>
            <div class="card-body">
                <div class="info-badge">
                    <strong>ℹ️ Interface Residues:</strong> See Node 5 analysis output for detailed contact information and interface residues.
                </div>
                <p>Residue-level interface analysis (hydrogen bonds, hydrophobic contacts, electrostatic interactions) is computed in Node 5.</p>
                <p class="text-muted small">The current comparison focuses on quantitative RMSD-based accuracy metrics using the backbone atoms (Cα) of the antigen.</p>
            </div>
        </div>
        
        <!-- Technical Details -->
        <div class="card mb-4">
            <div class="card-header bg-dark text-white">
                <h5 class="mb-0">⚙️ Technical Details</h5>
            </div>
            <div class="card-body">
                <h6>Alignment Method</h6>
                <p>Predicted antibody Cα atoms were rigidly superposed onto reference antibody Cα atoms using Biopython's SVD-based Superimposer. The same rotation/translation was applied to the antigen for RMSD computation.</p>
                
                <h6 class="mt-3">RMSD Definition</h6>
                <p class="code-block">RMSD = √( Σ(dist²) / N )</p>
                <p>where <code>dist</code> is the distance between paired atoms and <code>N</code> is the number of atoms.</p>
                
                <h6 class="mt-3">Interpretation</h6>
                <p>RMSD thresholds are based on common standards in structural biology for docking validation:</p>
                <ul>
                    <li>&lt; 2.0 Å: Excellent match (typically &lt; 0.5 Å error per residue)</li>
                    <li>2.0–5.0 Å: Good match (recovers main features)</li>
                    <li>5.0–10.0 Å: Moderate match (requires refinement)</li>
                    <li>&gt; 10.0 Å: Poor match (docking failure)</li>
                </ul>
            </div>
        </div>
        
        <div class="footer">
            <p>Report generated on <strong>{timestamp}</strong></p>
            <p>DiffDock‑PP Antibody‑Antigen Docking Pipeline | Node 6: Structure Comparison</p>
        </div>
    </div>
    
    <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/js/bootstrap.bundle.min.js"></script>
</body>
</html>
"""
    return html

def create_error_report(data):
    """Create an error report when comparison cannot be performed."""
    error = data.get('error', 'Unknown error occurred')
    timestamp = data.get('timestamp', 'Unknown')
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Comparison Failed</title>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.1.3/dist/css/bootstrap.min.css" rel="stylesheet">
    <style>
        body {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: #2c3e50;
        }}
        .header {{
            background: linear-gradient(135deg, #2c3e50 0%, #e74c3c 100%);
            color: white;
            padding: 2rem 0;
            margin-bottom: 2rem;
        }}
        .error-box {{
            background-color: #fff3cd;
            border-left: 4px solid #dc3545;
            padding: 1.5rem;
            border-radius: 0.25rem;
            color: #856404;
        }}
    </style>
</head>
<body>
    <div class="header">
        <div class="container">
            <h1>⚠️ Comparison Failed</h1>
            <p class="header-subtitle">Unable to complete antibody‑antigen structure comparison</p>
        </div>
    </div>
    
    <div class="container mt-4">
        <div class="card mb-4">
            <div class="card-header bg-danger text-white">
                <h5 class="mb-0">Error Report</h5>
            </div>
            <div class="card-body">
                <div class="error-box">
                    <h5>Error Details</h5>
                    <p><strong>Status:</strong> Comparison could not be performed</p>
                    <p><strong>Reason:</strong> {error}</p>
                    <p><strong>Timestamp:</strong> {timestamp}</p>
                </div>
                
                <h6 class="mt-4">Troubleshooting Steps</h6>
                <ol>
                    <li>Verify that both PDB files exist and are properly formatted</li>
                    <li>Check that chain identifiers are correctly specified</li>
                    <li>Ensure Cα atoms are present in both structures</li>
                    <li>Review logs for detailed error messages</li>
                    <li>Re-run Node 6 with valid inputs</li>
                </ol>
            </div>
        </div>
    </div>
</body>
</html>
"""
    return html

def main():
    parser = argparse.ArgumentParser(
        description='Generate HTML report for wet‑lab structure comparison'
    )
    parser.add_argument('--data_json', default='outputs/data.json',
                        help='Path to data.json from comparison science script')
    parser.add_argument('--output_html', default='outputs/report.html',
                        help='Output path for HTML report')
    args = parser.parse_args()
    
    try:
        with open(args.data_json) as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"Error: {args.data_json} not found", file=__import__('sys').stderr)
        return 1
    except json.JSONDecodeError:
        print(f"Error: {args.data_json} is not valid JSON", file=__import__('sys').stderr)
        return 1
    
    html_content = create_html_report(data)
    
    os.makedirs(os.path.dirname(args.output_html) or '.', exist_ok=True)
    with open(args.output_html, 'w') as f:
        f.write(html_content)
    
    print(f"HTML report written to {args.output_html}")
    return 0

if __name__ == '__main__':
    import sys
    sys.exit(main())