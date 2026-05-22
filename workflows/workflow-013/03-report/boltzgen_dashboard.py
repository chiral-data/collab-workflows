#!/usr/bin/env python3
"""
BoltzGen Binder Design Results Dashboard Generator

Generates an interactive HTML dashboard for analyzing BoltzGen binder design results.
Processes designed binder CIF files, aggregate metrics CSV, and quality scores.

Modeled after boltz_dashboard.py (Boltz-2 report) with adaptations for binder design metrics.

Usage:
    python boltzgen_dashboard.py <results_directory>

Expected files:
    - *.cif (designed binder structures)
    - aggregate_metrics_analyze.csv (design metrics)
    - results_overview.pdf (optional)

Output:
    - boltzgen_dashboard_<target_name>_<timestamp>.html
"""

import csv
import json
import os
import re
import sys
from datetime import datetime
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np

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
class DesignResult:
    """Data class for a single BoltzGen design result."""
    design_id: str
    cif_file: str
    confidence_score: float = 0.0
    plddt_mean: float = 0.0
    iptm: float = 0.0
    ptm: float = 0.0
    binding_energy: float = 0.0
    sequence_recovery: float = 0.0
    rmsd: float = 0.0
    rank: int = 0
    extra_metrics: Dict[str, Any] = field(default_factory=dict)


class BoltzGenDashboard:
    """Generate interactive HTML dashboard for BoltzGen binder design results."""

    COLORS = {
        'excellent': '#0f766e',
        'good': '#0369a1',
        'moderate': '#fb7c3c',
        'poor': '#dc2626',
    }

    QUALITY_THRESHOLDS = {
        'confidence': {'excellent': 0.9, 'good': 0.7, 'moderate': 0.5},
        'plddt': {'excellent': 0.9, 'good': 0.7, 'moderate': 0.5},
        'iptm': {'excellent': 0.8, 'good': 0.6, 'moderate': 0.4},
    }

    def __init__(self, results_dir: str):
        self.results_dir = Path(results_dir)
        self.designs: List[DesignResult] = []
        self.target_name = ""
        self.metrics_df: Optional[pd.DataFrame] = None

    def _parse_results(self):
        """Parse BoltzGen results from directory."""
        # Try loading aggregate metrics CSV first
        csv_files = list(self.results_dir.glob('aggregate_metrics*.csv'))
        if csv_files:
            self._parse_metrics_csv(csv_files[0])
        else:
            # Fall back to parsing CIF files directly
            self._parse_cif_files()

        if not self.target_name:
            # Infer target name from directory or file names
            cif_files = list(self.results_dir.glob('*.cif'))
            if cif_files:
                name = cif_files[0].stem
                name = re.sub(r'_design_?\d+.*', '', name)
                name = re.sub(r'_\d+$', '', name)
                self.target_name = name or "binder_design"
            else:
                self.target_name = "binder_design"

    def _parse_metrics_csv(self, csv_path: Path):
        """Parse aggregate metrics CSV from BoltzGen."""
        print(f"Loading metrics from: {csv_path}")
        self.metrics_df = pd.read_csv(csv_path)

        for _, row in self.metrics_df.iterrows():
            design_id = str(row.get('id', row.get('design_id', row.get('name', row.name))))

            # Find corresponding CIF file
            cif_candidates = list(self.results_dir.glob(f'*{design_id}*.cif'))
            cif_file = cif_candidates[0].name if cif_candidates else ''

            # Map BoltzGen CSV columns to DesignResult fields
            # BoltzGen outputs: id, iptm, ptm, complex_plddt, bb_rmsd, delta_sasa_refolded
            confidence = float(row.get('confidence_score',
                               row.get('design_to_target_iptm',
                               row.get('iptm', 0.0))))
            plddt = float(row.get('complex_plddt',
                          row.get('plddt',
                          row.get('plddt_mean', 0.0))))
            iptm = float(row.get('iptm',
                         row.get('design_iptm', 0.0)))
            ptm = float(row.get('ptm', 0.0))
            rmsd = float(row.get('bb_rmsd',
                         row.get('rmsd', 0.0)))
            binding_energy = float(row.get('binding_energy',
                                   row.get('delta_sasa_refolded', 0.0)))
            sequence_recovery = float(row.get('sequence_recovery', 0.0))

            extra = {}
            standard_cols = {'id', 'design_id', 'name', 'file_name',
                             'confidence_score', 'plddt', 'plddt_mean',
                             'complex_plddt', 'iptm', 'design_iptm',
                             'design_to_target_iptm', 'ptm',
                             'binding_energy', 'delta_sasa_refolded',
                             'sequence_recovery', 'rmsd', 'bb_rmsd'}
            for col in row.index:
                if col not in standard_cols:
                    val = row[col]
                    if pd.notna(val):
                        extra[col] = val

            design = DesignResult(
                design_id=design_id,
                cif_file=cif_file,
                confidence_score=confidence,
                plddt_mean=plddt,
                iptm=iptm,
                ptm=ptm,
                binding_energy=binding_energy,
                sequence_recovery=sequence_recovery,
                rmsd=rmsd,
                extra_metrics=extra,
            )
            self.designs.append(design)

    def _parse_cif_files(self):
        """Parse CIF files directly when no CSV is available."""
        cif_files = sorted(self.results_dir.glob('*.cif'))
        if not cif_files:
            raise ValueError("No CIF files or metrics CSV found in results directory")

        for i, cif_path in enumerate(cif_files):
            design = DesignResult(
                design_id=cif_path.stem,
                cif_file=cif_path.name,
            )
            self.designs.append(design)

    def generate_dashboard(self) -> str:
        """Generate the complete dashboard and return the output filename."""
        print("Parsing BoltzGen results...")
        self._parse_results()

        if not self.designs:
            raise ValueError("No valid BoltzGen results found")

        print(f"Found {len(self.designs)} designs for target: {self.target_name}")

        # Sort by confidence score (or iptm if confidence not available)
        self.designs.sort(
            key=lambda x: x.confidence_score if x.confidence_score > 0 else x.iptm,
            reverse=True,
        )

        for i, design in enumerate(self.designs):
            design.rank = i + 1

        html_content = self._generate_html_dashboard()

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"boltzgen_dashboard_{self.target_name}_{timestamp}.html"
        output_path = self.results_dir / output_file

        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)

        print(f"Dashboard generated: {output_file}")
        return output_file

    def _get_quality(self, value: float, metric: str = 'confidence') -> str:
        """Determine quality category based on metric value."""
        thresholds = self.QUALITY_THRESHOLDS.get(metric, self.QUALITY_THRESHOLDS['confidence'])
        if value >= thresholds['excellent']:
            return 'excellent'
        elif value >= thresholds['good']:
            return 'good'
        elif value >= thresholds['moderate']:
            return 'moderate'
        return 'poor'

    def _calculate_summary_metrics(self) -> Dict[str, Any]:
        """Calculate summary metrics for the dashboard."""
        if not self.designs:
            return {}

        confidences = [d.confidence_score for d in self.designs if d.confidence_score > 0]
        plddts = [d.plddt_mean for d in self.designs if d.plddt_mean > 0]
        iptms = [d.iptm for d in self.designs if d.iptm > 0]

        best = self.designs[0]

        # Primary scoring metric is iPTM (design_to_target_iptm)
        if confidences:
            best_score = max(confidences)
            score_label = "iPTM"
        elif iptms:
            best_score = max(iptms)
            score_label = "iPTM"
        else:
            best_score = 0.0
            score_label = "iPTM"

        quality = self._get_quality(best_score)

        return {
            'total_designs': len(self.designs),
            'best_design_id': best.design_id,
            'best_score': best_score,
            'score_label': score_label,
            'best_plddt': (max(plddts) if plddts else 0.0),
            'best_rmsd': min([d.rmsd for d in self.designs if d.rmsd > 0]) if any(d.rmsd > 0 for d in self.designs) else 0.0,
            'mean_score': np.mean(confidences) if confidences else (np.mean(iptms) if iptms else 0.0),
            'design_quality': quality,
            'target_name': self.target_name,
        }

    def _generate_recommendation(self, metrics: Dict[str, Any]) -> str:
        """Generate evidence-based recommendation."""
        quality = metrics['design_quality']
        score = metrics['best_score']
        label = metrics['score_label']

        if quality == 'excellent':
            return (f"Excellent binder candidates identified - {label}: {score:.3f}. "
                    f"Top designs show strong predicted binding. Recommended for experimental "
                    f"validation via SPR, ITC, or co-crystallization.")
        elif quality == 'good':
            return (f"Good binder candidates - {label}: {score:.3f}. "
                    f"Consider affinity maturation or sequence optimization for top hits. "
                    f"Validate binding experimentally before downstream use.")
        elif quality == 'moderate':
            return (f"Moderate quality designs - {label}: {score:.3f}. "
                    f"Results may benefit from additional design rounds with refined parameters, "
                    f"alternative scaffolds, or expanded sampling budget.")
        else:
            return (f"Low confidence designs - {label}: {score:.3f}. "
                    f"Consider alternative target binding sites, different design protocols, "
                    f"or increasing the sampling budget.")

    def _create_ranking_plot(self) -> str:
        """Create bar chart of design scores."""
        labels = [f"Design {d.design_id}" if len(d.design_id) < 10 else f"#{d.rank}" for d in self.designs[:20]]
        scores = [d.confidence_score if d.confidence_score > 0 else d.iptm for d in self.designs[:20]]
        qualities = [self._get_quality(s) for s in scores]
        colors = [self.COLORS[q] for q in qualities]

        fig = go.Figure()
        fig.add_trace(go.Bar(
            y=labels,
            x=scores,
            orientation='h',
            marker_color=colors,
            text=[f"{s:.3f}" for s in scores],
            textposition='inside',
            textfont=dict(color='white', size=12),
        ))

        fig.update_layout(
            title=dict(text=f"Design Ranking (Top {min(20, len(labels))})", x=0.5, font=dict(size=18)),
            xaxis=dict(title="iPTM Score", range=[0, 1.0], tickformat='.2f'),
            yaxis=dict(title="Designs", categoryorder="total ascending"),
            height=max(400, len(labels) * 25),
            template="plotly_white",
            margin=dict(l=120, r=20, t=80, b=60),
        )

        return pyo.plot(fig, output_type='div', include_plotlyjs=False)

    def _create_analysis_plots(self) -> str:
        """Create quality analysis subplot grid."""
        scores = [round(d.confidence_score if d.confidence_score > 0 else d.iptm, 6) for d in self.designs]
        rmsds = [round(d.rmsd, 6) for d in self.designs]
        design_labels = [d.design_id for d in self.designs]

        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=('iPTM Score by Design', 'RMSD by Design',
                            'Quality Overview', 'iPTM vs RMSD'),
            specs=[[{}, {}], [{"type": "domain"}, {}]],
        )

        # iPTM per design bar chart
        score_colors = [self.COLORS[self._get_quality(s)] for s in scores]
        fig.add_trace(go.Bar(
            x=design_labels, y=scores, name="iPTM",
            marker_color=score_colors,
            text=[f"{s:.3f}" for s in scores],
            textposition='outside', showlegend=False,
        ), row=1, col=1)

        # RMSD per design bar chart (lower is better)
        rmsd_colors = []
        for r in rmsds:
            if r < 2.0:
                rmsd_colors.append(self.COLORS['excellent'])
            elif r < 4.0:
                rmsd_colors.append(self.COLORS['good'])
            elif r < 6.0:
                rmsd_colors.append(self.COLORS['moderate'])
            else:
                rmsd_colors.append(self.COLORS['poor'])
        fig.add_trace(go.Bar(
            x=design_labels, y=rmsds, name="RMSD",
            marker_color=rmsd_colors,
            text=[f"{r:.2f}" for r in rmsds],
            textposition='outside', showlegend=False,
        ), row=1, col=2)

        # Quality pie
        quality_counts = {}
        for s in scores:
            q = self._get_quality(s)
            quality_counts[q] = quality_counts.get(q, 0) + 1

        fig.add_trace(go.Pie(
            labels=[q.title() for q in quality_counts],
            values=list(quality_counts.values()),
            marker=dict(colors=[self.COLORS[q] for q in quality_counts]),
        ), row=2, col=1)

        # iPTM vs RMSD scatter (high iPTM + low RMSD = best)
        fig.add_trace(go.Scatter(
            x=scores, y=rmsds, mode='markers+text', showlegend=False,
            marker=dict(size=10, color=score_colors, line=dict(width=1, color='white')),
            text=design_labels,
            textposition='top center',
            hovertemplate="<b>%{text}</b><br>iPTM: %{x:.3f}<br>RMSD: %{y:.2f}<extra></extra>",
        ), row=2, col=2)

        fig.update_xaxes(title_text="Design", row=1, col=1)
        fig.update_yaxes(title_text="iPTM", range=[0, max(scores) * 1.15], tickformat='.2f', row=1, col=1)
        fig.update_xaxes(title_text="Design", row=1, col=2)
        rmsd_max = max(rmsds) if rmsds else 10
        fig.update_yaxes(title_text="RMSD (\u00c5)", range=[0, rmsd_max * 1.2], row=1, col=2)
        fig.update_xaxes(title_text="iPTM", range=[0, max(scores) * 1.15], tickformat='.2f', row=2, col=2)
        fig.update_yaxes(title_text="RMSD (\u00c5)", range=[0, rmsd_max * 1.2], row=2, col=2)

        fig.update_layout(
            height=700, template="plotly_white",
            title=dict(text="Design Quality Analysis", x=0.5, font=dict(size=20)),
        )

        return pyo.plot(fig, output_type='div', include_plotlyjs=False)

    def _read_structure_files(self, top_n: int = 5) -> List[Dict[str, Any]]:
        """Read top N CIF structure files for 3D viewer embedding."""
        structures = []
        for design in self.designs[:top_n]:
            if not design.cif_file:
                continue
            cif_path = self.results_dir / design.cif_file
            if not cif_path.exists():
                print(f"Warning: Structure file not found: {cif_path}")
                continue
            try:
                raw = cif_path.read_text(encoding='utf-8')
                # Escape for safe JS template literal embedding
                escaped = raw.replace('\\', '\\\\').replace('`', '\\`').replace('${', '\\${')
                score = design.confidence_score if design.confidence_score > 0 else design.iptm
                structures.append({
                    'design_id': design.design_id,
                    'rank': design.rank,
                    'score': score,
                    'label': f"Rank {design.rank} - {design.design_id} (iPTM: {score:.3f})",
                    'data': escaped,
                    'format': 'cif',
                })
            except Exception as e:
                print(f"Warning: Could not read {cif_path}: {e}")
        return structures

    def _generate_viewer_section(self) -> str:
        """Generate 3Dmol.js viewer card HTML for embedding in dashboard."""
        structures = self._read_structure_files()
        if not structures:
            return """
        <div class="card">
            <div class="card-header">
                <i class="fas fa-cube"></i> 3D Structure Viewer
            </div>
            <div class="card-body">
                <p class="text-muted text-center py-4">No structure files (CIF) found for 3D viewing.</p>
            </div>
        </div>
"""

        # Build dropdown options
        options_html = ""
        for s in structures:
            selected = ' selected' if s['rank'] == 1 else ''
            options_html += f'<option value="{s["rank"]}"{selected}>{s["label"]}</option>\n'

        # Build JS structure data object
        struct_entries = []
        for s in structures:
            struct_entries.append(
                f'            {s["rank"]}: {{data: `{s["data"]}`, format: "{s["format"]}"}}'
            )
        struct_js = ',\n'.join(struct_entries)

        return f"""
        <div class="card">
            <div class="card-header">
                <i class="fas fa-cube"></i> 3D Structure Viewer
            </div>
            <div class="card-body">
                <div class="d-flex flex-wrap align-items-center gap-3 mb-3">
                    <div>
                        <label for="structureSelect" class="form-label mb-1" style="font-size: 0.85rem; color: #64748b;">Select Design</label>
                        <select id="structureSelect" class="form-select form-select-sm" style="min-width: 300px;" onchange="loadStructure(this.value)">
                            {options_html}
                        </select>
                    </div>
                    <div>
                        <label class="form-label mb-1" style="font-size: 0.85rem; color: #64748b;">Style</label>
                        <div class="btn-group btn-group-sm" role="group" id="styleButtons">
                            <button type="button" class="btn btn-outline-primary active" onclick="setViewerStyle('cartoon', this)">Cartoon</button>
                            <button type="button" class="btn btn-outline-primary" onclick="setViewerStyle('stick', this)">Stick</button>
                            <button type="button" class="btn btn-outline-primary" onclick="setViewerStyle('sphere', this)">Sphere</button>
                            <button type="button" class="btn btn-outline-primary" onclick="setViewerStyle('line', this)">Line</button>
                            <button type="button" class="btn btn-outline-primary" onclick="setViewerStyle('surface', this)">Surface</button>
                        </div>
                    </div>
                </div>
                <div id="viewer3d" style="height: 500px; width: 100%; position: relative; background: #1a1a2e; border-radius: 10px;"></div>
            </div>
        </div>

        <script>
        (function() {{
            var structureData = {{
{struct_js}
            }};

            var viewer = null;
            var currentStyle = 'cartoon';

            function initViewer() {{
                var container = document.getElementById('viewer3d');
                viewer = $3Dmol.createViewer(container, {{backgroundColor: '#1a1a2e'}});
                // Load the first structure (rank 1)
                var firstRank = Object.keys(structureData)[0];
                loadStructure(firstRank);
            }}

            window.loadStructure = function(rank) {{
                if (!viewer || !structureData[rank]) return;
                viewer.removeAllModels();
                viewer.removeAllSurfaces();
                viewer.addModel(structureData[rank].data, structureData[rank].format);
                applyStyle(currentStyle);
                viewer.zoomTo();
                viewer.render();
            }};

            function applyStyle(style) {{
                if (!viewer) return;
                viewer.removeAllSurfaces();
                viewer.setStyle({{}}, {{}});
                var styleSpec = {{}};
                switch(style) {{
                    case 'cartoon':
                        styleSpec = {{cartoon: {{colorscheme: 'chainHetatm'}}}};
                        break;
                    case 'stick':
                        styleSpec = {{stick: {{colorscheme: 'chainHetatm', radius: 0.12}}}};
                        break;
                    case 'sphere':
                        styleSpec = {{sphere: {{colorscheme: 'chainHetatm', scale: 0.3}}}};
                        break;
                    case 'line':
                        styleSpec = {{line: {{colorscheme: 'chainHetatm'}}}};
                        break;
                    case 'surface':
                        styleSpec = {{cartoon: {{colorscheme: 'chainHetatm'}}}};
                        viewer.addSurface($3Dmol.SurfaceType.VDW, {{opacity: 0.7, color: 'white'}});
                        break;
                }}
                viewer.setStyle({{}}, styleSpec);
                viewer.render();
            }}

            window.setViewerStyle = function(style, btn) {{
                currentStyle = style;
                document.querySelectorAll('#styleButtons .btn').forEach(function(b) {{
                    b.classList.remove('active');
                }});
                if (btn) btn.classList.add('active');
                applyStyle(style);
            }};

            // Initialize when DOM and 3Dmol.js are ready
            if (document.readyState === 'loading') {{
                document.addEventListener('DOMContentLoaded', initViewer);
            }} else {{
                initViewer();
            }}
        }})();
        </script>
"""

    def _generate_html_dashboard(self) -> str:
        """Generate the complete HTML dashboard."""
        metrics = self._calculate_summary_metrics()
        recommendation = self._generate_recommendation(metrics)

        # Build table rows
        table_rows = []
        for d in self.designs:
            score = d.confidence_score if d.confidence_score > 0 else d.iptm
            quality = self._get_quality(score)
            badge_color = self.COLORS[quality]

            cif_link = (f'<a href="{d.cif_file}" class="btn btn-sm btn-outline-primary">CIF</a>'
                        if d.cif_file else '-')

            sequence = d.extra_metrics.get('designed_chain_sequence', '')

            table_rows.append(f"""
                <tr>
                    <td><strong>{d.design_id}</strong></td>
                    <td>{d.rank}</td>
                    <td>{score:.4f}</td>
                    <td><span class="badge" style="background-color: {badge_color}; color: white; padding: 3px 8px; border-radius: 12px; font-size: 11px;">{quality.title()}</span></td>
                    <td>{d.plddt_mean:.3f}</td>
                    <td>{d.ptm:.4f}</td>
                    <td>{d.rmsd:.3f}</td>
                    <td>{d.binding_energy:.1f}</td>
                    <td><code style="font-size: 11px;">{sequence}</code></td>
                    <td>{cif_link}</td>
                </tr>
            """)

        table_html = "".join(table_rows)

        ranking_plot = self._create_ranking_plot()
        analysis_plots = self._create_analysis_plots()
        viewer_section = self._generate_viewer_section()

        quality_bg = {
            'excellent': '#d4edda', 'good': '#d1ecf1',
            'moderate': '#fff3cd', 'poor': '#f8d7da',
        }
        assessment_bg = quality_bg.get(metrics['design_quality'], '#f8f9fa')

        html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BoltzGen Binder Design Results - {metrics['target_name']}</title>
    <script src="https://cdn.plot.ly/plotly-2.24.1.min.js"></script>
    <script src="https://3dmol.org/build/3Dmol-min.js"></script>
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
            min-height: 100vh; margin: 0; color: #1f2937;
        }}
        .main-container {{ max-width: 1400px; margin: 0 auto; padding: 20px; }}
        .header {{
            background: rgba(255,255,255,0.95); backdrop-filter: blur(10px);
            border-radius: 15px; padding: 30px; margin-bottom: 30px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.1); border: 1px solid rgba(255,255,255,0.2);
        }}
        .header h1 {{ color: var(--primary-color); font-size: 2.5rem; font-weight: 700; margin: 0 0 10px 0; display: flex; align-items: center; gap: 15px; }}
        .header .subtitle {{ color: #64748b; font-size: 1.1rem; margin: 0; }}
        .card {{
            background: rgba(255,255,255,0.95); backdrop-filter: blur(10px);
            border: 1px solid rgba(255,255,255,0.2); border-radius: 15px;
            box-shadow: 0 8px 32px rgba(0,0,0,0.1); margin-bottom: 30px;
            transition: transform 0.2s ease, box-shadow 0.2s ease;
        }}
        .card:hover {{ transform: translateY(-2px); box-shadow: 0 12px 40px rgba(0,0,0,0.15); }}
        .card-header {{
            background: linear-gradient(135deg, var(--primary-color) 0%, var(--secondary-color) 100%);
            color: white; border-radius: 15px 15px 0 0 !important; border: none;
            padding: 20px 25px; font-weight: 600; font-size: 1.1rem;
        }}
        .summary-cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 20px; margin-bottom: 30px; }}
        .summary-card {{
            background: rgba(255,255,255,0.95); backdrop-filter: blur(10px);
            border-radius: 15px; padding: 25px; text-align: center;
            box-shadow: 0 8px 32px rgba(0,0,0,0.1); border: 1px solid rgba(255,255,255,0.2);
            transition: transform 0.2s ease;
        }}
        .summary-card:hover {{ transform: translateY(-3px); }}
        .summary-card .icon {{ font-size: 2.5rem; margin-bottom: 15px; opacity: 0.8; }}
        .summary-card .value {{ font-size: 2rem; font-weight: 700; margin: 10px 0; color: var(--primary-color); }}
        .summary-card .label {{ color: #64748b; font-size: 0.9rem; text-transform: uppercase; letter-spacing: 1px; }}
        .table th {{ background-color: #f8fafc; border-top: none; font-weight: 600; color: var(--primary-color); }}
        .table tbody tr:hover {{ background-color: rgba(7,89,133,0.05); }}
        .assessment-panel {{
            background: {assessment_bg};
            border-left: 5px solid var(--{metrics['design_quality']}-color);
            border-radius: 10px; padding: 20px; margin: 20px 0;
        }}
        .plot-container {{ background: white; border-radius: 10px; padding: 20px; margin: 15px 0; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }}
        @media (max-width: 768px) {{
            .summary-cards {{ grid-template-columns: 1fr; }}
            .header h1 {{ font-size: 2rem; }}
            .main-container {{ padding: 15px; }}
        }}
    </style>
</head>
<body>
    <div class="main-container">
        <div class="header">
            <h1><i class="fas fa-dna"></i> BoltzGen Binder Design Results</h1>
            <p class="subtitle">
                <strong>Target:</strong> {metrics['target_name']} |
                <strong>Designs:</strong> {metrics['total_designs']} |
                <strong>Generated:</strong> {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
            </p>
        </div>

        <div class="summary-cards">
            <div class="summary-card">
                <div class="icon"><i class="fas fa-trophy" style="color: var(--excellent-color);"></i></div>
                <div class="value">{metrics['best_score']:.3f}</div>
                <div class="label">Best {metrics['score_label']}</div>
            </div>
            <div class="summary-card">
                <div class="icon"><i class="fas fa-chart-line" style="color: var(--good-color);"></i></div>
                <div class="value">{metrics['best_plddt']:.3f}</div>
                <div class="label">Best pLDDT</div>
            </div>
            <div class="summary-card">
                <div class="icon"><i class="fas fa-ruler" style="color: var(--moderate-color);"></i></div>
                <div class="value">{metrics['best_rmsd']:.2f} &Aring;</div>
                <div class="label">Lowest RMSD</div>
            </div>
            <div class="summary-card">
                <div class="icon"><i class="fas fa-medal" style="color: var(--primary-color);"></i></div>
                <div class="value">{metrics['design_quality'].title()}</div>
                <div class="label">Design Quality</div>
            </div>
        </div>

        <div class="card">
            <div class="card-header">
                <i class="fas fa-lightbulb"></i> Quality Assessment & Recommendations
            </div>
            <div class="card-body">
                <div class="assessment-panel">
                    <h5><i class="fas fa-microscope"></i> Design Quality: {metrics['design_quality'].title()}</h5>
                    <p class="mb-3">{recommendation}</p>
                    <div class="row">
                        <div class="col-md-6">
                            <strong>Quality Metrics:</strong>
                            <ul class="mt-2">
                                <li>Best {metrics['score_label']}: {metrics['best_score']:.3f}</li>
                                <li>Best pLDDT: {metrics['best_plddt']:.3f}</li>
                                <li>Lowest RMSD: {metrics['best_rmsd']:.2f} &Aring;</li>
                            </ul>
                        </div>
                        <div class="col-md-6">
                            <strong>Recommended Next Steps:</strong>
                            <ul class="mt-2">
                                <li>Validate top designs experimentally</li>
                                <li>Assess binding affinity via SPR or ITC</li>
                                <li>Consider structural refinement with MD</li>
                            </ul>
                        </div>
                    </div>
                </div>
            </div>
        </div>

        <div class="card">
            <div class="card-header">
                <i class="fas fa-chart-bar"></i> Design Ranking
            </div>
            <div class="card-body">
                <div class="plot-container">{ranking_plot}</div>
            </div>
        </div>

        <div class="card">
            <div class="card-header">
                <i class="fas fa-chart-area"></i> Quality Analysis
            </div>
            <div class="card-body">
                <div class="plot-container">{analysis_plots}</div>
            </div>
        </div>

        {viewer_section}

        <div class="card">
            <div class="card-header">
                <i class="fas fa-table"></i> Detailed Results
            </div>
            <div class="card-body">
                <div class="table-responsive">
                    <table class="table table-striped table-hover">
                        <thead>
                            <tr>
                                <th>Design</th>
                                <th>Rank</th>
                                <th>iPTM</th>
                                <th>Quality</th>
                                <th>pLDDT</th>
                                <th>PTM</th>
                                <th>RMSD (&Aring;)</th>
                                <th>dSASA (&Aring;&sup2;)</th>
                                <th>Sequence</th>
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

    <script>
        // Plots are embedded as Plotly divs above
    </script>
</body>
</html>"""

        return html


def main():
    """Main function to generate BoltzGen dashboard."""
    if len(sys.argv) > 1:
        results_dir = sys.argv[1]
    else:
        results_dir = "."

    print(f"Processing BoltzGen results from: {results_dir}")
    dashboard = BoltzGenDashboard(results_dir)
    output_file = dashboard.generate_dashboard()
    print(f"Dashboard saved as: {output_file}")


if __name__ == "__main__":
    main()
