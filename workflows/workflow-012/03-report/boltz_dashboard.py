#!/usr/bin/env python3
"""
Boltz-2 Structure Prediction Results Dashboard Generator

This script generates an interactive HTML dashboard for analyzing Boltz-2 structure prediction results.
It processes confidence scores, pLDDT values, and other quality metrics to provide comprehensive
analysis and recommendations for structure prediction workflows.

Usage:
    python boltz_dashboard.py <results_directory>

The script expects the following files in the results directory:
    - confidence_*_model_*.json (confidence scores and metrics)
    - *_model_*.pdb (structure files)
    - pae_*_model_*.npz (Predicted Aligned Error matrices)
    - pde_*_model_*.npz (Predicted Distance Error matrices)
    - plddt_*_model_*.npz (pLDDT confidence scores)

Output:
    - boltz_dashboard_<protein_name>_<timestamp>.html
"""

import json
import sys
import re
from pathlib import Path
from datetime import datetime
from dataclasses import dataclass
from typing import List, Dict, Any, Optional, Tuple
import numpy as np

# Try importing required libraries
try:
    import pandas as pd
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    import plotly.offline as pyo
except ImportError as e:
    print(f"Error importing required libraries: {e}")
    print("Please ensure pandas, plotly, and numpy are installed:")
    print("pip install pandas plotly numpy")
    sys.exit(1)

@dataclass
class BoltzModel:
    """Data class for storing Boltz-2 model results."""
    model_id: int
    confidence_score: float
    ptm: float
    iptm: float
    plddt_mean: float
    pde_mean: float
    structure_file: str
    confidence_file: str
    pae_file: str
    pde_file: str
    plddt_file: str
    rank: int = 0

class BoltzDashboard:
    """Generate interactive HTML dashboard for Boltz-2 structure prediction results."""
    
    def __init__(self, results_dir: str):
        self.results_dir = Path(results_dir)
        self.models: List[BoltzModel] = []
        self.protein_name = ""
        self.job_id = ""
    
    def _parse_results(self):
        """Parse Boltz-2 results from the directory."""
        confidence_files = sorted(self.results_dir.glob('confidence_*_model_*.json'))
        
        if not confidence_files:
            raise ValueError("No confidence files found. Expected format: confidence_*_model_*.json")
            
        for conf_file in confidence_files:
            try:
                # Extract protein name and model ID from filename
                match = re.search(r'confidence_(.+)_model_(\d+)\.json', conf_file.name)
                if not match:
                    continue
                    
                protein_name = match.group(1)
                model_id = int(match.group(2))
                
                if not self.protein_name:
                    self.protein_name = protein_name
                elif self.protein_name != protein_name:
                    continue  # Skip if different protein
                    
                # Load confidence data
                with open(conf_file, 'r') as f:
                    conf_data = json.load(f)
                    
                # Find corresponding files
                structure_file = f"{protein_name}_model_{model_id}.pdb"
                pae_file = f"pae_{protein_name}_model_{model_id}.npz"
                pde_file = f"pde_{protein_name}_model_{model_id}.npz"
                plddt_file = f"plddt_{protein_name}_model_{model_id}.npz"

                # Calculate pLDDT mean from complex_plddt or load from file
                plddt_mean = conf_data.get('complex_plddt', 0.0)
                if plddt_mean == 0.0 and (self.results_dir / plddt_file).exists():
                    plddt_data = np.load(self.results_dir / plddt_file)
                    if 'plddt' in plddt_data:
                        plddt_mean = float(np.mean(plddt_data['plddt']) / 100.0)
                            
                # Calculate PDE mean
                pde_mean = conf_data.get('complex_pde', 0.0)
                
                model = BoltzModel(
                    model_id=model_id,
                    confidence_score=conf_data.get('confidence_score', 0.0),
                    ptm=conf_data.get('ptm', 0.0),
                    iptm=conf_data.get('iptm', 0.0),
                    plddt_mean=plddt_mean,
                    pde_mean=pde_mean,
                    structure_file=structure_file,
                    confidence_file=conf_file.name,
                    pae_file=pae_file,
                    pde_file=pde_file,
                    plddt_file=plddt_file
                )
                
                self.models.append(model)
                
            except (json.JSONDecodeError, KeyError, ValueError) as e:
                print(f"Warning: Could not parse {conf_file.name}: {e}")
                continue
                
        if not self.protein_name:
            self.protein_name = "unknown"
    
    def generate_dashboard(self) -> str:
        """Generate the complete dashboard and return the output filename."""
        print("Parsing Boltz-2 results...")
        self._parse_results()
        
        if not self.models:
            raise ValueError("No valid Boltz-2 results found in the specified directory")
            
        print(f"Found {len(self.models)} models for protein: {self.protein_name}")
        
        # Sort models by confidence score (descending)
        self.models.sort(key=lambda x: x.confidence_score, reverse=True)
        
        # Assign ranks
        for i, model in enumerate(self.models):
            model.rank = i + 1
            
        # Generate HTML
        html_content = self._generate_html_dashboard()
        
        # Save to file
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"boltz_dashboard_{self.protein_name}_{timestamp}.html"
        output_path = self.results_dir / output_file
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
            
        print(f"Dashboard generated: {output_file}")
        return output_file

    # GROMACS-exact color scheme matching  
    COLORS = {
        'page_bg': 'linear-gradient(135deg, #075985 0%, #0284c7 100%)',  # GROMACS blue gradient
        'background': '#ffffff',           # White content areas like GROMACS
        'card_bg': 'linear-gradient(135deg, #f8f9fa 0%, #e9ecef 100%)',  # GROMACS card gradient
        'excellent': '#0f766e',           # GROMACS teal for excellent (energy/success)
        'good': '#0369a1',               # GROMACS dark blue for good (structural/strong)
        'moderate': '#fb7c3c',           # Complementary orange for moderate
        'poor': '#dc2626',               # Red for poor
        'text_primary': '#2c3e50',       # GROMACS dark gray for body text
        'text_secondary': '#6c757d',     # GROMACS medium gray for labels
        'header_color': '#0369a1',       # GROMACS header blue
        'table_header': '#0284c7',       # GROMACS medium blue for table headers
        'border': '#e9ecef'              # Light gray borders
    }
    
    # Boltz-2 quality thresholds based on research
    QUALITY_THRESHOLDS = {
        'confidence': {
            'excellent': 0.9,
            'good': 0.7, 
            'moderate': 0.5
        },
        'plddt': {
            'excellent': 0.9,
            'good': 0.7,
            'moderate': 0.5
        }
    }

    def _get_confidence_quality(self, confidence: float) -> str:
        """Determine quality category based on confidence score."""
        if confidence >= self.QUALITY_THRESHOLDS['confidence']['excellent']:
            return 'excellent'
        elif confidence >= self.QUALITY_THRESHOLDS['confidence']['good']:
            return 'good'
        elif confidence >= self.QUALITY_THRESHOLDS['confidence']['moderate']:
            return 'moderate'
        else:
            return 'poor'
            
    def _get_plddt_quality(self, plddt: float) -> str:
        """Determine quality category based on pLDDT score (0-1 scale)."""
        if plddt >= self.QUALITY_THRESHOLDS['plddt']['excellent']:
            return 'excellent'
        elif plddt >= self.QUALITY_THRESHOLDS['plddt']['good']:
            return 'good'
        elif plddt >= self.QUALITY_THRESHOLDS['plddt']['moderate']:
            return 'moderate'
        else:
            return 'poor'
            
    def _assess_structure_quality(self, confidence: float, plddt: float) -> str:
        """Assess overall structure quality using both confidence and pLDDT."""
        # Both metrics should be good for overall excellent quality
        if (confidence >= self.QUALITY_THRESHOLDS['confidence']['excellent'] and 
            plddt >= self.QUALITY_THRESHOLDS['plddt']['excellent']):
            return 'excellent'
        elif (confidence >= self.QUALITY_THRESHOLDS['confidence']['good'] and 
              plddt >= self.QUALITY_THRESHOLDS['plddt']['good']):
            return 'good'
        elif (confidence >= self.QUALITY_THRESHOLDS['confidence']['moderate'] and 
              plddt >= self.QUALITY_THRESHOLDS['plddt']['moderate']):
            return 'moderate'
        else:
            return 'poor'
        
    def _calculate_summary_metrics(self) -> Dict[str, Any]:
        """Calculate summary metrics for the dashboard."""
        if not self.models:
            return {}
            
        best_model = self.models[0]  # Already sorted by confidence
        confidence_scores = [m.confidence_score for m in self.models]
        plddt_scores = [m.plddt_mean for m in self.models]
        
        # Calculate consistency (coefficient of variation)
        confidence_cv = np.std(confidence_scores) / np.mean(confidence_scores) if np.mean(confidence_scores) > 0 else 0
        
        # Overall structure quality assessment (conservative approach)
        structure_quality = self._assess_structure_quality(best_model.confidence_score, best_model.plddt_mean)
        
        return {
            'total_models': len(self.models),
            'best_model_id': best_model.model_id,
            'best_confidence': best_model.confidence_score,
            'best_plddt': max(plddt_scores),
            'mean_confidence': np.mean(confidence_scores),
            'mean_plddt': np.mean(plddt_scores),
            'confidence_cv': confidence_cv,
            'model_consistency': confidence_cv,
            'structure_quality': structure_quality,
            'protein_name': self.protein_name
        }
        
    def _generate_recommendation(self, metrics: Dict[str, Any]) -> str:
        """Generate evidence-based recommendations based on Boltz-2 quality standards."""
        quality = metrics['structure_quality']
        confidence = metrics['best_confidence']
        plddt = metrics['best_plddt']
        consistency = metrics['confidence_cv']
        
        if quality == 'excellent':
            if consistency < 0.1:
                return f"Near-Atomic Quality Structure - Confidence: {confidence:.3f}, pLDDT: {plddt:.3f}. Excellent for molecular dynamics, drug design, functional analysis, and publication. Both backbone and side-chains highly accurate."
            else:
                return f"High-Quality Structure with Model Diversity - Best model shows near-atomic quality (Conf: {confidence:.3f}, pLDDT: {plddt:.3f}). Consider ensemble approaches for critical applications."
        elif quality == 'good':
            if consistency < 0.15:
                return f"High-Confidence Structure - Confidence: {confidence:.3f}, pLDDT: {plddt:.3f}. Suitable for most applications including binding site analysis, comparative modeling, and initial drug screening."
            else:
                return f"Good Quality with Variability - Top model reliable (Conf: {confidence:.3f}, pLDDT: {plddt:.3f}). Validate critical regions and consider experimental confirmation for detailed applications."
        elif quality == 'moderate':
            return f"Moderate Confidence Prediction - Confidence: {confidence:.3f}, pLDDT: {plddt:.3f}. Use with caution. Consider experimental validation, focus on high-confidence regions, or use for preliminary analysis only."
        else:
            return f"Low Confidence Structure - Confidence: {confidence:.3f}, pLDDT: {plddt:.3f}. Not recommended for critical applications. Consider alternative prediction methods, experimental approaches, or sequence-based analysis."
    
    def _create_confidence_plot(self) -> str:
        """Create horizontal bar chart of confidence scores."""
        color_map = {
            'excellent': self.COLORS['excellent'],
            'good': self.COLORS['good'],
            'moderate': self.COLORS['moderate'],
            'poor': self.COLORS['poor']
        }

        # Build plain lists per model to avoid binary data encoding issues
        model_labels = []
        model_scores = []
        model_colors = []
        model_texts = []
        for model in reversed(self.models):  # reversed so rank 1 is on top
            quality = self._get_confidence_quality(model.confidence_score)
            model_labels.append(f"Model {model.model_id} (Rank {model.rank})")
            model_scores.append(round(model.confidence_score, 6))
            model_colors.append(color_map[quality])
            model_texts.append(f"{model.confidence_score:.3f}")

        fig = go.Figure()
        fig.add_trace(go.Bar(
            y=model_labels,
            x=model_scores,
            orientation='h',
            marker_color=model_colors,
            text=model_texts,
            textposition='inside',
            textfont=dict(color='white', size=12),
            width=0.6,
        ))

        # Zoom x-axis to data range for better visual differentiation
        x_min = max(0, min(model_scores) - 0.05)
        x_max = min(1.0, max(model_scores) + 0.02)

        fig.update_layout(
            title=dict(
                text="Model Confidence Scores",
                x=0.5,
                font=dict(size=18, family="Arial, sans-serif")
            ),
            xaxis=dict(
                title="Confidence Score",
                range=[x_min, x_max],
                tickformat='.3f',
                gridcolor='lightgray'
            ),
            yaxis=dict(
                title="Models",
            ),
            height=max(300, 120 + len(self.models) * 50),
            template="plotly_white",
            showlegend=False,
            margin=dict(l=120, r=20, t=80, b=60)
        )

        return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    
    def _create_quality_analysis_plots(self) -> str:
        """Create improved quality analysis plots with better visualization."""
        # Create subplot with 2 rows, 2 columns for better analysis
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('Confidence Score by Model', 'pLDDT Score by Model',
                           'Model Quality Overview', 'Confidence vs pLDDT Correlation'),
            specs=[[{"secondary_y": False}, {"secondary_y": False}],
                   [{"type": "domain"}, {"secondary_y": False}]]
        )
        
        # Data preparation - use plain lists to avoid binary serialization
        model_labels = [f"Model {m.model_id}" for m in self.models]
        confidence_scores = [round(m.confidence_score, 6) for m in self.models]
        plddt_scores = [round(m.plddt_mean, 6) for m in self.models]

        # 1. Confidence Score per Model (top-left)
        colors_conf = [self.COLORS[self._get_confidence_quality(c)] for c in confidence_scores]
        fig.add_trace(
            go.Bar(
                x=model_labels,
                y=confidence_scores,
                name="Confidence",
                marker_color=colors_conf,
                text=[f"{c:.3f}" for c in confidence_scores],
                textposition='outside',
                showlegend=False
            ),
            row=1, col=1
        )

        # 2. pLDDT Score per Model (top-right)
        colors_plddt = [self.COLORS[self._get_plddt_quality(p)] for p in plddt_scores]
        fig.add_trace(
            go.Bar(
                x=model_labels,
                y=plddt_scores,
                name="pLDDT",
                marker_color=colors_plddt,
                text=[f"{p:.3f}" for p in plddt_scores],
                textposition='outside',
                showlegend=False
            ),
            row=1, col=2
        )
        
        # 3. Model Quality Overview (bottom-left) - Quality distribution pie chart
        quality_counts = {}
        for model in self.models:
            quality = self._get_confidence_quality(model.confidence_score)
            quality_counts[quality] = quality_counts.get(quality, 0) + 1
            
        colors = [self.COLORS[quality] for quality in quality_counts.keys()]
        
        fig.add_trace(
            go.Pie(
                labels=[q.title() for q in quality_counts.keys()],
                values=list(quality_counts.values()),
                marker=dict(colors=colors),
                showlegend=True
            ),
            row=2, col=1
        )
        
        # 4. Confidence vs pLDDT Correlation (bottom-right)
        # Color by quality for better insight
        colors_scatter = [self.COLORS[self._get_confidence_quality(conf)] for conf in confidence_scores]
        
        fig.add_trace(
            go.Scatter(
                x=confidence_scores,
                y=plddt_scores,
                mode='markers',
                name="Models",
                marker=dict(
                    size=12,
                    color=colors_scatter,
                    line=dict(width=2, color='white')
                ),
                text=[f"Model {model.model_id}" for model in self.models],
                hovertemplate="<b>%{text}</b><br>" +
                            "Confidence: %{x:.3f}<br>" +
                            "pLDDT: %{y:.3f}<br>" +
                            "<extra></extra>",
                showlegend=False
            ),
            row=2, col=2
        )
        
        # Update subplot layouts
        conf_min = max(0, min(confidence_scores) - 0.05)
        conf_max = min(1.0, max(confidence_scores) + 0.02)
        plddt_min = max(0, min(plddt_scores) - 0.05)
        plddt_max = min(1.0, max(plddt_scores) + 0.02)

        fig.update_xaxes(title_text="Model", row=1, col=1)
        fig.update_yaxes(title_text="Confidence Score", range=[conf_min, conf_max], tickformat='.3f', row=1, col=1)
        fig.update_xaxes(title_text="Model", row=1, col=2)
        fig.update_yaxes(title_text="pLDDT Score", range=[plddt_min, plddt_max], tickformat='.3f', row=1, col=2)
        fig.update_xaxes(title_text="Confidence Score", range=[conf_min, conf_max], tickformat='.3f', row=2, col=2)
        fig.update_yaxes(title_text="pLDDT Score", range=[plddt_min, plddt_max], tickformat='.3f', row=2, col=2)
        
        fig.update_layout(
            height=600,
            showlegend=True,
            template="plotly_white",
            title=dict(
                text="Quality Analysis Dashboard",
                x=0.5,
                font=dict(size=20, family="Arial, sans-serif")
            )
        )
        
        return pyo.plot(fig, output_type='div', include_plotlyjs=False)
    
    def _generate_html_dashboard(self) -> str:
        """Generate the complete HTML dashboard using GROMACS styling."""
        metrics = self._calculate_summary_metrics()
        recommendation = self._generate_recommendation(metrics)

        # Generate plots
        confidence_plot = self._create_confidence_plot()
        quality_plots = self._create_quality_analysis_plots()

        # Generate detailed results table
        table_rows = []
        for model in self.models:
            confidence_quality = self._get_confidence_quality(model.confidence_score)
            plddt_quality = self._get_plddt_quality(model.plddt_mean)
            
            # Quality badges with GROMACS colors
            conf_badge_color = self.COLORS[confidence_quality]
            plddt_badge_color = self.COLORS[plddt_quality]
            
            table_rows.append(f"""
                <tr>
                    <td><strong>Model {model.model_id}</strong></td>
                    <td>{model.rank}</td>
                    <td>{model.confidence_score:.4f}</td>
                    <td><span class="badge" style="background-color: {conf_badge_color}; color: white; padding: 3px 8px; border-radius: 12px; font-size: 11px;">{confidence_quality.title()}</span></td>
                    <td>{model.plddt_mean:.3f}</td>
                    <td><span class="badge" style="background-color: {plddt_badge_color}; color: white; padding: 3px 8px; border-radius: 12px; font-size: 11px;">{plddt_quality.title()}</span></td>
                    <td>{model.ptm:.4f}</td>
                    <td>{model.iptm:.4f}</td>
                    <td>{model.pde_mean:.4f}</td>
                    <td><a href="{model.structure_file}" class="btn btn-sm btn-outline-primary">PDB</a></td>
                </tr>
            """)
        
        table_html = "".join(table_rows)
        
        # Quality assessment panel color using GROMACS scheme
        quality_colors = {
            'excellent': '#d4edda', 'good': '#d1ecf1', 'moderate': '#fff3cd', 'poor': '#f8d7da'
        }
        assessment_bg_color = quality_colors.get(metrics['structure_quality'], '#f8f9fa')
        
        # HTML template with GROMACS styling (converted from AutoDock dashboard)
        html_template = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Boltz-2 Structure Prediction Results - {metrics['protein_name']}</title>
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    <link href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.0/dist/css/bootstrap.min.css" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css" rel="stylesheet">
    <style>
        :root {{
            --primary-color: #075985;
            --secondary-color: #0284c7;
            --excellent-color: #0f766e;
            --good-color: #0369a1;
            --moderate-color: #fb7c3c;
            --poor-color: #dc2626;
        }}
        
        body {{
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
            background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
            min-height: 100vh;
            margin: 0;
            color: #1f2937;
        }}
        
        .main-container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        
        .header {{
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(10px);
            border-radius: 15px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
            border: 1px solid rgba(255, 255, 255, 0.2);
        }}
        
        .header h1 {{
            color: var(--primary-color);
            font-size: 2.5rem;
            font-weight: 700;
            margin: 0 0 10px 0;
            display: flex;
            align-items: center;
            gap: 15px;
        }}
        
        .header .subtitle {{
            color: #64748b;
            font-size: 1.1rem;
            margin: 0;
        }}
        
        .card {{
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(10px);
            border: 1px solid rgba(255, 255, 255, 0.2);
            border-radius: 15px;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
            margin-bottom: 30px;
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }}
        
        .card:hover {{
            transform: translateY(-2px);
            box-shadow: 0 12px 40px rgba(0, 0, 0, 0.15);
        }}
        
        .card-header {{
            background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
            color: white;
            border-radius: 15px 15px 0 0 !important;
            border: none;
            padding: 20px 25px;
            font-weight: 600;
            font-size: 1.1rem;
        }}
        
        .summary-cards {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .summary-card {{
            background: rgba(255, 255, 255, 0.95);
            backdrop-filter: blur(10px);
            border-radius: 15px;
            padding: 25px;
            text-align: center;
            box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
            border: 1px solid rgba(255, 255, 255, 0.2);
            transition: transform 0.2s ease;
        }}
        
        .summary-card:hover {{
            transform: translateY(-3px);
        }}
        
        .summary-card .icon {{
            font-size: 2.5rem;
            margin-bottom: 15px;
            opacity: 0.8;
        }}
        
        .summary-card .value {{
            font-size: 2rem;
            font-weight: 700;
            margin: 10px 0;
            color: var(--primary-color);
        }}
        
        .summary-card .label {{
            color: #64748b;
            font-size: 0.9rem;
            text-transform: uppercase;
            letter-spacing: 1px;
        }}
        
        .badge-excellent {{ background-color: var(--excellent-color) !important; }}
        .badge-good {{ background-color: var(--good-color) !important; }}
        .badge-moderate {{ background-color: var(--moderate-color) !important; }}
        .badge-poor {{ background-color: var(--poor-color) !important; }}
        
        .table th {{
            background-color: #f8fafc;
            border-top: none;
            font-weight: 600;
            color: var(--primary-color);
        }}
        
        .table tbody tr:hover {{
            background-color: rgba(7, 89, 133, 0.05);
        }}
        
        .assessment-panel {{
            background: {assessment_bg_color};
            border-left: 5px solid var(--{metrics['structure_quality']}-color);
            border-radius: 10px;
            padding: 20px;
            margin: 20px 0;
        }}
        
        .plot-container {{
            background: white;
            border-radius: 10px;
            padding: 20px;
            margin: 15px 0;
            box-shadow: 0 2px 10px rgba(0, 0, 0, 0.1);
        }}
        
        @media (max-width: 768px) {{
            .summary-cards {{
                grid-template-columns: 1fr;
            }}
            
            .header h1 {{
                font-size: 2rem;
            }}
            
            .main-container {{
                padding: 15px;
            }}
        }}
    </style>
</head>
<body>
    <div class="main-container">
        <div class="header">
            <h1><i class="fas fa-cubes"></i>Boltz-2 Structure Prediction Results</h1>
            <p class="subtitle">
                <strong>Protein:</strong> {metrics['protein_name']} | 
                <strong>Models:</strong> {len(self.models)} | 
                <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
            </p>
        </div>
        
        <div class="summary-cards">
            <div class="summary-card">
                <div class="icon"><i class="fas fa-trophy" style="color: var(--excellent-color);"></i></div>
                <div class="value">{metrics['best_confidence']:.3f}</div>
                <div class="label">Best Confidence</div>
            </div>
            
            <div class="summary-card">
                <div class="icon"><i class="fas fa-chart-line" style="color: var(--good-color);"></i></div>
                <div class="value">{metrics['best_plddt']:.1f}</div>
                <div class="label">Best pLDDT</div>
            </div>
            
            <div class="summary-card">
                <div class="icon"><i class="fas fa-medal" style="color: var(--moderate-color);"></i></div>
                <div class="value">{metrics['structure_quality'].title()}</div>
                <div class="label">Structure Quality</div>
            </div>
            
            <div class="summary-card">
                <div class="icon"><i class="fas fa-balance-scale" style="color: var(--primary-color);"></i></div>
                <div class="value">{metrics['model_consistency']:.3f}</div>
                <div class="label">Model Consistency</div>
            </div>
        </div>
        
        <div class="card">
            <div class="card-header">
                <i class="fas fa-lightbulb"></i> Quality Assessment & Recommendations
            </div>
            <div class="card-body">
                <div class="assessment-panel">
                    <h5><i class="fas fa-microscope"></i> Structure Quality: {metrics['structure_quality'].title()}</h5>
                    <p class="mb-3">{recommendation}</p>
                    
                    <div class="row">
                        <div class="col-md-6">
                            <strong>Quality Metrics:</strong>
                            <ul class="mt-2">
                                <li>Best Confidence: {metrics['best_confidence']:.3f}</li>
                                <li>Best pLDDT: {metrics['best_plddt']:.1f}</li>
                                <li>Model Consistency: {metrics['model_consistency']:.3f}</li>
                            </ul>
                        </div>
                        <div class="col-md-6">
                            <strong>Recommended Actions:</strong>
                            <ul class="mt-2">
                                <li>Validate with experimental data if available</li>
                                <li>Consider molecular dynamics simulation</li>
                                <li>Assess functional relevance of predictions</li>
                            </ul>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <div class="card">
            <div class="card-header">
                <i class="fas fa-chart-bar"></i> Confidence Scores
            </div>
            <div class="card-body">
                <div class="plot-container">{confidence_plot}</div>
            </div>
        </div>

        <div class="card">
            <div class="card-header">
                <i class="fas fa-chart-area"></i> Quality Analysis
            </div>
            <div class="card-body">
                <div class="plot-container">{quality_plots}</div>
            </div>
        </div>

        <div class="card">
            <div class="card-header">
                <i class="fas fa-table"></i> Detailed Results
            </div>
            <div class="card-body">
                <div class="table-responsive">
                    <table class="table table-striped table-hover">
                        <thead>
                            <tr>
                                <th>Model</th>
                                <th>Rank</th>
                                <th>Confidence</th>
                                <th>Conf. Quality</th>
                                <th>pLDDT</th>
                                <th>pLDDT Quality</th>
                                <th>PTM</th>
                                <th>iPTM</th>
                                <th>PDE</th>
                                <th>Structure</th>
                            </tr>
                        </thead>
                        <tbody>
                            {table_html}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>
    </div>
</body>
</html>"""
        
        return html_template


def main():
    """Main function to generate Boltz-2 dashboard"""
    import sys
    
    if len(sys.argv) > 1:
        results_dir = sys.argv[1]
    else:
        results_dir = "."
    
    print(f"Processing Boltz-2 results from: {results_dir}")
    dashboard = BoltzDashboard(results_dir)
    output_file = dashboard.generate_dashboard()
    print(f"Dashboard saved as: {output_file}")


if __name__ == "__main__":
    main()
