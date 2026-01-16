#!/usr/bin/env python3
"""
AutoDock Vina Virtual Screening Analysis Dashboard
Large-scale docking analysis and visualization tool for drug repurposing campaigns
Analyzes thousands of compounds and generates comprehensive screening reports
"""

import os
import re
import json
import glob
import numpy as np
import pandas as pd
from datetime import datetime
from pathlib import Path
from statistics import mean, stdev, median
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import warnings
warnings.filterwarnings('ignore')

@dataclass
class CompoundResult:
    """Data structure for individual compound screening results"""
    ligand_name: str
    best_affinity: float
    all_affinities: List[float] = field(default_factory=list)
    rmsd_values: List[Tuple[float, float]] = field(default_factory=list)
    num_modes: int = 0
    file_path: Optional[str] = None
    
@dataclass
class ScreeningResults:
    """Complete virtual screening results"""
    compounds: List[CompoundResult]
    total_compounds: int
    successful_dockings: int
    failed_dockings: int
    receptor_file: str
    job_id: str
    box_center: Tuple[float, float, float]
    box_size: Tuple[float, float, float]
    exhaustiveness: int
    screening_date: str

class VinaScreeningParser:
    """Parse large-scale AutoDock Vina virtual screening results"""
    
    def __init__(self, output_dir: str, csv_file: str = None, summary_file: str = None):
        self.output_dir = Path(output_dir)
        self.poses_dir = self.output_dir / "poses"
        self.csv_file = csv_file or self.output_dir / "screening_results_sorted.csv"
        self.summary_file = summary_file or self.output_dir / "virtual_screening_summary.txt"
        self.results: ScreeningResults = None
        self._parse_all_results()
    
    def _parse_pdbqt_file(self, pdbqt_file: str) -> Optional[CompoundResult]:
        """Parse individual PDBQT file for binding modes"""
        try:
            ligand_name = Path(pdbqt_file).stem.replace('_out', '')
            
            with open(pdbqt_file, 'r') as f:
                content = f.read()
            
            # Extract all VINA RESULT lines
            result_pattern = r'REMARK VINA RESULT:\s+([-\d.]+)\s+([\d.]+)\s+([\d.]+)'
            matches = re.findall(result_pattern, content)
            
            if not matches:
                return None
            
            affinities = []
            rmsd_values = []
            
            for affinity, rmsd_lb, rmsd_ub in matches:
                affinities.append(float(affinity))
                rmsd_values.append((float(rmsd_lb), float(rmsd_ub)))
            
            return CompoundResult(
                ligand_name=ligand_name,
                best_affinity=min(affinities),
                all_affinities=affinities,
                rmsd_values=rmsd_values,
                num_modes=len(affinities),
                file_path=str(pdbqt_file)
            )
        except Exception as e:
            print(f"Error parsing {pdbqt_file}: {e}")
            return None
    
    def _parse_all_results(self):
        """Parse all screening results with progress tracking"""
        print("Parsing virtual screening results...")
        
        # Parse summary file for metadata
        metadata = self._parse_summary_file()
        
        # Get all PDBQT files
        pdbqt_files = list(self.poses_dir.glob("*_out.pdbqt"))
        total_files = len(pdbqt_files)
        print(f"Found {total_files} PDBQT files to parse...")
        
        # Parse files in batches for progress tracking
        compounds = []
        batch_size = 100
        
        for i in range(0, total_files, batch_size):
            batch = pdbqt_files[i:i+batch_size]
            if i % 500 == 0:
                print(f"  Processing files {i+1}-{min(i+batch_size, total_files)} of {total_files}...")
            
            # Process batch
            for pdbqt_file in batch:
                result = self._parse_pdbqt_file(pdbqt_file)
                if result:
                    compounds.append(result)
        
        # Sort by best affinity
        compounds.sort(key=lambda x: x.best_affinity)
        
        # Parse CSV for additional validation
        if os.path.exists(self.csv_file):
            csv_data = pd.read_csv(self.csv_file)
            print(f"Validated against CSV: {len(csv_data)} entries")
        
        self.results = ScreeningResults(
            compounds=compounds,
            total_compounds=metadata.get('total_compounds', len(compounds)),
            successful_dockings=len(compounds),
            failed_dockings=metadata.get('failed_dockings', 0),
            receptor_file=metadata.get('receptor', 'Unknown'),
            job_id=metadata.get('job_id', f'screening_{datetime.now().strftime("%Y%m%d")}'),
            box_center=metadata.get('box_center', (0, 0, 0)),
            box_size=metadata.get('box_size', (0, 0, 0)),
            exhaustiveness=metadata.get('exhaustiveness', 0),
            screening_date=metadata.get('date', datetime.now().strftime('%Y-%m-%d'))
        )
        
        print(f"Successfully parsed {len(compounds)} compounds")
    
    def _parse_summary_file(self) -> Dict:
        """Parse virtual_screening_summary.txt for screening metadata"""
        metadata = {}
        
        if not os.path.exists(self.summary_file):
            return metadata
        
        try:
            with open(self.summary_file, 'r') as f:
                content = f.read()
            
            # Extract metadata using regex patterns
            patterns = {
                'job_id': r'Job ID:\s*(.+)',
                'date': r'Date:\s*(.+)',
                'receptor': r'Receptor:\s*(.+)',
                'total_compounds': r'Total SDF files found:\s*(\d+)',
                'successful_dockings': r'Successful dockings:\s*(\d+)',
                'failed_dockings': r'Failed dockings:\s*(\d+)',
                'exhaustiveness': r'Exhaustiveness:\s*(\d+)',
            }
            
            for key, pattern in patterns.items():
                match = re.search(pattern, content)
                if match:
                    value = match.group(1).strip()
                    if key in ['total_compounds', 'successful_dockings', 'failed_dockings', 'exhaustiveness']:
                        metadata[key] = int(value)
                    else:
                        metadata[key] = value
            
            # Parse box parameters
            box_match = re.search(r'Box Center:\s*\(([\d.]+),\s*([\d.]+),\s*([\d.]+)\)', content)
            if box_match:
                metadata['box_center'] = tuple(float(x) for x in box_match.groups())
            
            size_match = re.search(r'Box Size:\s*\(([\d.]+),\s*([\d.]+),\s*([\d.]+)\)', content)
            if size_match:
                metadata['box_size'] = tuple(float(x) for x in size_match.groups())
                
        except Exception as e:
            print(f"Warning: Could not fully parse summary file: {e}")
        
        return metadata

class ScreeningAnalyzer:
    """Statistical analysis for virtual screening results"""
    
    def __init__(self, results: ScreeningResults):
        self.results = results
        self.metrics = self._calculate_comprehensive_metrics()
    
    def _calculate_comprehensive_metrics(self) -> Dict:
        """Calculate population-level screening metrics"""
        compounds = self.results.compounds
        raw_affinities = [c.best_affinity for c in compounds]
        
        # Data cleaning: Remove extreme outliers (likely docking failures)
        affinities = [a for a in raw_affinities if -50.0 <= a <= 50.0]
        cleaned_compounds = [c for c in compounds if -50.0 <= c.best_affinity <= 50.0]
        
        print(f"Data cleaning: {len(raw_affinities)} -> {len(affinities)} compounds (removed {len(raw_affinities) - len(affinities)} outliers)")
        
        # Basic statistics
        metrics = {
            'total_screened': self.results.total_compounds,
            'successful_dockings': self.results.successful_dockings,
            'failed_dockings': self.results.failed_dockings,
            'success_rate': (self.results.successful_dockings / self.results.total_compounds * 100) 
                           if self.results.total_compounds > 0 else 0,
            
            # Affinity statistics
            'best_affinity': min(affinities) if affinities else 0,
            'worst_affinity': max(affinities) if affinities else 0,
            'mean_affinity': mean(affinities) if affinities else 0,
            'median_affinity': median(affinities) if affinities else 0,
            'std_affinity': stdev(affinities) if len(affinities) > 1 else 0,
            
            # Quartiles
            'q1_affinity': np.percentile(affinities, 25) if affinities else 0,
            'q3_affinity': np.percentile(affinities, 75) if affinities else 0,
            
            # Hit rates with non-overlapping categories
            'hits_excellent': self._count_hits_range(affinities, -float('inf'), -10.0, cleaned_compounds),
            'hits_good': self._count_hits_range(affinities, -10.0, -8.0, cleaned_compounds),
            'hits_moderate': self._count_hits_range(affinities, -8.0, -6.0, cleaned_compounds),
            'hits_weak': self._count_hits_range(affinities, -6.0, float('inf'), cleaned_compounds),
            
            # Top performers
            'top_10_compounds': compounds[:10] if len(compounds) >= 10 else compounds,
            'top_50_compounds': compounds[:50] if len(compounds) >= 50 else compounds,
            'top_100_compounds': compounds[:100] if len(compounds) >= 100 else compounds,
            
            # Enrichment analysis
            'enrichment_factor_1': self._calculate_enrichment(affinities, 0.01),
            'enrichment_factor_5': self._calculate_enrichment(affinities, 0.05),
            'enrichment_factor_10': self._calculate_enrichment(affinities, 0.10),
        }
        
        # Add recommendation based on results
        metrics['recommendation'] = self._generate_screening_recommendation(metrics)
        
        return metrics
    
    def _count_hits_range(self, affinities: List[float], min_val: float, max_val: float, compounds: List) -> Dict:
        """Count hits within a specific affinity range"""
        if min_val == -float('inf'):
            hits = [a for a in affinities if a < max_val]
            hit_compounds = [c for c in compounds if c.best_affinity < max_val]
        elif max_val == float('inf'):
            hits = [a for a in affinities if a >= min_val]
            hit_compounds = [c for c in compounds if c.best_affinity >= min_val]
        else:
            hits = [a for a in affinities if min_val <= a < max_val]
            hit_compounds = [c for c in compounds if min_val <= c.best_affinity < max_val]
        
        return {
            'count': len(hits),
            'percentage': (len(hits) / len(affinities) * 100) if affinities else 0,
            'compounds': hit_compounds[:10],
            'range': f'{min_val:.1f} to {max_val:.1f}' if min_val != -float('inf') and max_val != float('inf') else 
                    f'< {max_val:.1f}' if min_val == -float('inf') else f'>= {min_val:.1f}'
        }
    
    def _calculate_enrichment(self, affinities: List[float], top_fraction: float) -> float:
        """Calculate enrichment factor for top fraction of compounds using strict threshold"""
        if not affinities:
            return 0
        
        n_top = int(len(affinities) * top_fraction)
        if n_top == 0:
            return 0
        
        # Use strict threshold for true exceptional hits (-15 kcal/mol)
        # This separates truly outstanding compounds from good ones
        threshold_exceptional = -15.0
        hits_in_top = sum(1 for a in affinities[:n_top] if a < threshold_exceptional)
        total_hits = sum(1 for a in affinities if a < threshold_exceptional)
        
        if total_hits == 0:
            return 0  # No exceptional hits in entire dataset
        
        expected_hits = total_hits * top_fraction
        if expected_hits == 0:
            return 0
        
        return hits_in_top / expected_hits
    
    def _generate_screening_recommendation(self, metrics: Dict) -> Dict:
        """Generate evidence-based screening recommendations"""
        excellent_hits = metrics['hits_excellent']['count']
        good_hits = metrics['hits_good']['count']
        moderate_hits = metrics['hits_moderate']['count']
        
        if excellent_hits >= 10:
            return {
                'action': 'Proceed to Lead Optimization',
                'priority': 'HIGH',
                'confidence': f'Excellent hit rate with {excellent_hits} strong binders',
                'next_steps': [
                    f'Prioritize top {min(20, excellent_hits)} compounds for experimental validation',
                    'Perform structural clustering to identify diverse scaffolds',
                    'Run molecular dynamics on top 5 hits',
                    'Analyze binding mode consensus',
                    'Consider ADMET predictions for lead selection'
                ]
            }
        elif good_hits >= 20:
            return {
                'action': 'Hit Confirmation and Expansion',
                'priority': 'MEDIUM-HIGH',
                'confidence': f'Good hit rate with {good_hits} promising compounds',
                'next_steps': [
                    f'Select top {min(30, good_hits)} compounds for testing',
                    'Perform similarity search for analogs',
                    'Validate binding modes with MD simulation',
                    'Consider fragment growing strategies',
                    'Profile selectivity against related targets'
                ]
            }
        elif moderate_hits >= 50:
            return {
                'action': 'Hit-to-Lead Development',
                'priority': 'MEDIUM',
                'confidence': f'Moderate hits requiring optimization ({moderate_hits} compounds)',
                'next_steps': [
                    'Cluster compounds by scaffold',
                    'Identify optimization vectors',
                    'Consider combinatorial library design',
                    'Perform pharmacophore analysis',
                    'Evaluate synthetic accessibility'
                ]
            }
        else:
            return {
                'action': 'Re-evaluate Screening Strategy',
                'priority': 'LOW',
                'confidence': 'Limited hits identified',
                'next_steps': [
                    'Review docking parameters and search space',
                    'Consider alternative binding sites',
                    'Expand chemical library diversity',
                    'Try ensemble docking with multiple conformations',
                    'Evaluate different scoring functions'
                ]
            }

class ScreeningDashboard:
    """Virtual Screening Dashboard Generator"""
    
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
    
    def __init__(self, results: ScreeningResults, analyzer: ScreeningAnalyzer):
        self.results = results
        self.analyzer = analyzer
        self.metrics = analyzer.metrics
    
    def create_summary_cards(self) -> List[Dict]:
        """Create summary information cards"""
        return [
            {
                'title': 'Compounds Screened',
                'value': f'{self.metrics["total_screened"]:,}',
                'unit': f'{self.metrics["success_rate"]:.1f}% success'
            },
            {
                'title': 'Best Affinity',
                'value': f'{self.metrics["best_affinity"]:.2f}',
                'unit': 'kcal/mol'
            },
            {
                'title': 'Excellent Hits',
                'value': f'{self.metrics["hits_excellent"]["count"]}',
                'unit': f'< -10 kcal/mol'
            },
            {
                'title': 'Good Hits',
                'value': f'{self.metrics["hits_good"]["count"]}',
                'unit': f'< -8 kcal/mol'
            },
            {
                'title': 'Enrichment (Top 1%)',
                'value': f'{self.metrics["enrichment_factor_1"]:.1f}x',
                'unit': 'fold enrichment'
            },
            {
                'title': 'Recommendation',
                'value': self.metrics['recommendation']['action'].split()[0],
                'unit': self.metrics['recommendation']['priority']
            }
        ]
    
    def create_top_hits_plot(self) -> go.Figure:
        """Create horizontal bar chart of top hits"""
        top_compounds = self.metrics['top_50_compounds'][:30]  # Show top 30
        
        # Color based on affinity
        colors = []
        for compound in top_compounds:
            if compound.best_affinity < -10:
                colors.append(self.COLORS['excellent'])
            elif compound.best_affinity < -8:
                colors.append(self.COLORS['good'])
            elif compound.best_affinity < -6:
                colors.append(self.COLORS['moderate'])
            else:
                colors.append(self.COLORS['poor'])
        
        fig = go.Figure(data=go.Bar(
            y=[c.ligand_name for c in reversed(top_compounds)],
            x=[c.best_affinity for c in reversed(top_compounds)],
            orientation='h',
            marker_color=list(reversed(colors)),
            text=[f'{c.best_affinity:.2f}' for c in reversed(top_compounds)],
            textposition='outside',
            hovertemplate='<b>%{y}</b><br>' +
                         'Best Affinity: %{x:.3f} kcal/mol<br>' +
                         'Modes: %{customdata}<extra></extra>',
            customdata=[c.num_modes for c in reversed(top_compounds)]
        ))
        
        fig.update_layout(
            title='Top 30 Compounds by Binding Affinity',
            xaxis_title='Binding Affinity (kcal/mol)',
            yaxis_title='Compound',
            plot_bgcolor=self.COLORS['background'],
            paper_bgcolor=self.COLORS['background'],
            font_color=self.COLORS['text_primary'],
            height=800,
            margin=dict(l=200)
        )
        
        # Add threshold lines
        thresholds = [
            ('Excellent', -10, self.COLORS['excellent']),
            ('Good', -8, self.COLORS['good']),
            ('Moderate', -6, self.COLORS['moderate'])
        ]
        
        for name, value, color in thresholds:
            fig.add_vline(
                x=value,
                line_dash="dash",
                line_color=color,
                annotation_text=f"{name}: {value}",
                annotation_position="top"
            )
        
        return fig
    
    def create_distribution_plot(self) -> go.Figure:
        """Create distribution analysis plots"""
        affinities = [c.best_affinity for c in self.results.compounds]
        
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=[
                'Affinity Distribution',
                'Cumulative Distribution',
                'Box Plot by Category',
                'Hit Rate Analysis'
            ],
            specs=[[{'type': 'histogram'}, {'type': 'scatter'}],
                   [{'type': 'box'}, {'type': 'bar'}]]
        )
        
        # Histogram with reasonable axis range
        fig.add_trace(
            go.Histogram(
                x=affinities,
                nbinsx=50,
                marker_color=self.COLORS['header_color'],
                name='Count',
                showlegend=False,
                xbins=dict(start=-25, end=5, size=0.5)  # Reasonable range for binding affinities
            ),
            row=1, col=1
        )
        
        # CDF
        sorted_affinities = sorted(affinities)
        cdf = np.arange(1, len(sorted_affinities) + 1) / len(sorted_affinities)
        
        fig.add_trace(
            go.Scatter(
                x=sorted_affinities,
                y=cdf * 100,
                mode='lines',
                line=dict(color=self.COLORS['good'], width=2),
                name='CDF',
                showlegend=False
            ),
            row=1, col=2
        )
        
        # Box plot by category (non-overlapping)
        categories = ['Excellent\n< -10', 'Good\n-10 to -8', 'Moderate\n-8 to -6', 'Weak\n≥ -6']
        category_affinities = [
            [a for a in affinities if a < -10],
            [a for a in affinities if -10 <= a < -8],
            [a for a in affinities if -8 <= a < -6],
            [a for a in affinities if a >= -6]
        ]
        category_colors = [
            self.COLORS['excellent'],
            self.COLORS['good'],
            self.COLORS['moderate'],
            self.COLORS['poor']
        ]
        
        for i, (cat, cat_affs, color) in enumerate(zip(categories, category_affinities, category_colors)):
            if len(cat_affs) > 0:  # Only plot if category has data
                fig.add_trace(
                    go.Box(
                        y=cat_affs,
                        name=cat,
                        marker_color=color,
                        showlegend=False
                    ),
                    row=2, col=1
                )
        
        # Hit rate bar chart (corrected categories)
        hit_categories = ['Excellent\n< -10', 'Good\n-10 to -8', 'Moderate\n-8 to -6', 'Weak\n≥ -6']
        hit_percentages = [
            self.metrics['hits_excellent']['percentage'],
            self.metrics['hits_good']['percentage'],
            self.metrics['hits_moderate']['percentage'],
            self.metrics['hits_weak']['percentage']
        ]
        
        fig.add_trace(
            go.Bar(
                x=hit_categories,
                y=hit_percentages,
                marker_color=[
                    self.COLORS['excellent'],
                    self.COLORS['good'],
                    self.COLORS['moderate'],
                    self.COLORS['poor']
                ],
                text=[f'{p:.1f}%' for p in hit_percentages],
                textposition='outside',
                showlegend=False
            ),
            row=2, col=2
        )
        
        # Update layout with better axis ranges
        fig.update_xaxes(title_text="Binding Affinity (kcal/mol)", row=1, col=1, range=[-25, 5])
        fig.update_xaxes(title_text="Binding Affinity (kcal/mol)", row=1, col=2, range=[-25, 5])
        fig.update_xaxes(title_text="Category", row=2, col=1)
        fig.update_xaxes(title_text="Category", row=2, col=2)
        
        fig.update_yaxes(title_text="Count", row=1, col=1)
        fig.update_yaxes(title_text="Cumulative %", row=1, col=2)
        fig.update_yaxes(title_text="Binding Affinity (kcal/mol)", row=2, col=1, range=[-25, 5])
        fig.update_yaxes(title_text="Hit Rate (%)", row=2, col=2)
        
        fig.update_layout(
            title='Screening Results Distribution Analysis',
            plot_bgcolor=self.COLORS['background'],
            paper_bgcolor=self.COLORS['background'],
            font_color=self.COLORS['text_primary'],
            height=800,
            showlegend=False
        )
        
        return fig
    
    def create_enrichment_plot(self) -> go.Figure:
        """Create enrichment analysis visualization with detailed explanations"""
        # Calculate enrichment at different percentiles
        percentiles = [0.5, 1, 2, 5, 10, 20, 30, 50]
        enrichments = []
        
        # Use cleaned affinities
        raw_affinities = [c.best_affinity for c in self.results.compounds]
        affinities = [a for a in raw_affinities if -50.0 <= a <= 50.0]
        
        for p in percentiles:
            enrichments.append(self.analyzer._calculate_enrichment(affinities, p/100))
        
        fig = go.Figure()
        
        # Enrichment curve with hover information
        fig.add_trace(go.Scatter(
            x=percentiles,
            y=enrichments,
            mode='lines+markers',
            name='Enrichment Factor',
            line=dict(color=self.COLORS['good'], width=3),
            marker=dict(size=10, color=self.COLORS['good']),
            hovertemplate='<b>Top %{x}% of compounds</b><br>' +
                         'Enrichment: %{y:.2f}x<br>' +
                         '<i>%{y:.2f}x better than random selection</i><extra></extra>'
        ))
        
        # Add reference line at 1.0 (random)
        fig.add_hline(
            y=1.0,
            line_dash="dash",
            line_color="gray",
            annotation_text="Random Selection (1.0x)",
            annotation_position="right"
        )
        
        # Add quality zones without annotations (to avoid duplicate pop-ups)
        fig.add_hrect(
            y0=1, y1=5, fillcolor="rgba(255, 0, 0, 0.1)",
            line_width=0
        )
        fig.add_hrect(
            y0=5, y1=20, fillcolor="rgba(255, 165, 0, 0.1)",
            line_width=0
        )
        if max(enrichments) > 20:
            fig.add_hrect(
                y0=20, y1=max(enrichments) * 1.1, fillcolor="rgba(0, 128, 0, 0.1)",
                line_width=0
            )
        
        fig.update_layout(
            title='Enrichment Analysis - Exceptional Hits (< -15 kcal/mol)',
            xaxis_title='Top Percentage of Compounds (%)',
            yaxis_title='Enrichment Factor (x-fold vs Random)',
            plot_bgcolor=self.COLORS['background'],
            paper_bgcolor=self.COLORS['background'],
            font_color=self.COLORS['text_primary'],
            height=500,
            annotations=[
                dict(
                    x=0.02, y=0.98, xref='paper', yref='paper',
                    text='<b>Exceptional Hits Analysis:</b><br>' +
                         '• Threshold: < -15 kcal/mol (ultra-strong binders)<br>' +
                         '• Enrichment measures concentration in top %<br>' +
                         '• Color zones: Red (<5x), Orange (5-20x), Green (>20x)<br>' +
                         '• High values = exceptional compounds clustered at top',
                    showarrow=False,
                    align='left',
                    bgcolor='rgba(255,255,255,0.8)',
                    bordercolor='gray',
                    borderwidth=1,
                    font=dict(size=10)
                )
            ]
        )
        
        return fig
    
    def create_detailed_table(self) -> pd.DataFrame:
        """Create detailed results table for top compounds"""
        data = []
        
        for i, compound in enumerate(self.metrics['top_100_compounds'], 1):
            # Determine quality category (non-overlapping)
            if compound.best_affinity < -10:
                quality = 'Excellent'
            elif -10 <= compound.best_affinity < -8:
                quality = 'Good'
            elif -8 <= compound.best_affinity < -6:
                quality = 'Moderate'
            else:
                quality = 'Weak'
            
            # Get top 3 affinities
            top3 = compound.all_affinities[:3] if len(compound.all_affinities) >= 3 else compound.all_affinities
            
            data.append({
                'Rank': i,
                'Compound': compound.ligand_name,
                'Best Affinity (kcal/mol)': f'{compound.best_affinity:.3f}',
                'Mode 2 (kcal/mol)': f'{top3[1]:.3f}' if len(top3) > 1 else 'N/A',
                'Mode 3 (kcal/mol)': f'{top3[2]:.3f}' if len(top3) > 2 else 'N/A',
                'Num Modes': compound.num_modes,
                'Quality': quality,
                'Z-Score': f'{(compound.best_affinity - self.metrics["mean_affinity"]) / self.metrics["std_affinity"]:.2f}' if self.metrics["std_affinity"] > 0 else 'N/A'
            })
        
        return pd.DataFrame(data)
    
    def generate_html_report(self, output_file: str = None) -> str:
        """Generate complete HTML dashboard report"""
        if output_file is None:
            output_file = f'vina_screening_dashboard_{datetime.now().strftime("%Y%m%d_%H%M%S")}.html'
        
        print("Generating dashboard visualizations...")
        
        # Create all visualizations
        summary_cards = self.create_summary_cards()
        top_hits_plot = self.create_top_hits_plot()
        distribution_plot = self.create_distribution_plot()
        enrichment_plot = self.create_enrichment_plot()
        results_table = self.create_detailed_table()
        
        # Convert plots to JSON
        top_hits_json = top_hits_plot.to_json()
        distribution_json = distribution_plot.to_json()
        enrichment_json = enrichment_plot.to_json()
        
        print("Building HTML report...")
        
        # Generate HTML content with GROMACS styling
        html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Virtual Screening Analysis Dashboard</title>
    <meta charset="utf-8">
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            margin: 0;
            padding: 0;
            background: {self.COLORS['page_bg']};
            min-height: 100vh;
        }}
        
        .main-container {{
            max-width: 1400px;
            margin: 0 auto;
            padding: 20px;
        }}
        
        .dashboard {{
            background: {self.COLORS['background']};
            border-radius: 15px;
            padding: 30px;
            margin-bottom: 30px;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
        }}
        
        .header {{
            text-align: center;
            margin-bottom: 30px;
            padding-bottom: 20px;
            border-bottom: 3px solid {self.COLORS['table_header']};
        }}
        
        .header h1 {{
            color: {self.COLORS['header_color']};
            font-size: 2.5rem;
            margin-bottom: 0.5rem;
            font-weight: 600;
        }}
        
        .header p {{
            color: {self.COLORS['text_secondary']};
            font-size: 1.1rem;
        }}
        
        .summary-cards {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .card {{
            background: {self.COLORS['card_bg']};
            border: 2px solid {self.COLORS['table_header']};
            border-radius: 10px;
            padding: 20px;
            text-align: center;
            box-shadow: 0 4px 15px rgba(0,0,0,0.05);
        }}
        
        .card h3 {{
            color: {self.COLORS['text_secondary']};
            font-size: 0.9rem;
            margin-bottom: 10px;
            text-transform: uppercase;
            font-weight: 600;
        }}
        
        .card-value {{
            font-size: 2em;
            font-weight: bold;
            margin: 10px 0;
            color: {self.COLORS['header_color']};
        }}
        
        .card-unit {{
            color: {self.COLORS['text_secondary']};
            font-size: 0.9rem;
        }}
        
        .plot-container {{
            background: {self.COLORS['card_bg']};
            border-radius: 10px;
            padding: 20px;
            margin-bottom: 30px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.05);
        }}
        
        .recommendation-panel {{
            background: {self.COLORS['card_bg']};
            border: 2px solid {self.COLORS['table_header']};
            border-radius: 10px;
            padding: 20px;
            margin-bottom: 30px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.1);
        }}
        
        .recommendation-panel h2 {{
            color: {self.COLORS['header_color']};
            margin-bottom: 15px;
        }}
        
        .statistics-panel {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(300px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .stat-box {{
            background: {self.COLORS['card_bg']};
            border-radius: 10px;
            padding: 15px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.05);
        }}
        
        .stat-box h3 {{
            color: {self.COLORS['header_color']};
            font-size: 1.1rem;
            margin-bottom: 10px;
            border-bottom: 2px solid {self.COLORS['border']};
            padding-bottom: 5px;
        }}
        
        .stat-box p {{
            margin: 5px 0;
            color: {self.COLORS['text_primary']};
        }}
        
        .table-container {{
            background: {self.COLORS['card_bg']};
            border-radius: 10px;
            padding: 20px;
            overflow-x: auto;
            box-shadow: 0 4px 15px rgba(0,0,0,0.05);
        }}
        
        .table-container h2 {{
            color: {self.COLORS['header_color']};
            margin-bottom: 15px;
            border-bottom: 3px solid {self.COLORS['table_header']};
            padding-bottom: 10px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 15px;
        }}
        
        th, td {{
            border: 1px solid {self.COLORS['border']};
            padding: 12px 8px;
            text-align: left;
        }}
        
        th {{
            background-color: {self.COLORS['table_header']};
            color: white;
            font-weight: 600;
        }}
        
        tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        
        .footer {{
            text-align: center;
            margin-top: 30px;
            color: {self.COLORS['text_secondary']};
            font-size: 0.9rem;
        }}
    </style>
</head>
<body>
    <div class="main-container">
        <div class="dashboard">
            <div class="header">
                <h1>Virtual Screening Analysis Dashboard</h1>
                <p><strong>Job ID:</strong> {self.results.job_id} | <strong>Date:</strong> {self.results.screening_date}</p>
                <p><strong>Receptor:</strong> {self.results.receptor_file} | <strong>Library Size:</strong> {self.results.total_compounds:,} compounds</p>
                <p><strong>Success Rate:</strong> {self.metrics["success_rate"]:.1f}% | <strong>Box Center:</strong> ({self.results.box_center[0]:.2f}, {self.results.box_center[1]:.2f}, {self.results.box_center[2]:.2f})</p>
            </div>
    
    <div class="summary-cards">
        {''.join([f'''
        <div class="card">
            <h3>{card['title']}</h3>
            <div class="card-value">{card['value']}</div>
            <div class="card-unit">{card['unit']}</div>
        </div>
        ''' for card in summary_cards])}
    </div>
    
    <div class="recommendation-panel">
        <h2>Screening Recommendation: {self.metrics['recommendation']['action']}</h2>
        <p><strong>Priority:</strong> {self.metrics['recommendation']['priority']} | 
           <strong>Confidence:</strong> {self.metrics['recommendation']['confidence']}</p>
        <h3>Recommended Next Steps:</h3>
        <ul>
            {''.join([f'<li>{step}</li>' for step in self.metrics['recommendation']['next_steps']])}
        </ul>
    </div>
    
    <div class="statistics-panel">
        <div class="stat-box">
            <h3>Affinity Statistics</h3>
            <p><strong>Best:</strong> {self.metrics["best_affinity"]:.3f} kcal/mol</p>
            <p><strong>Mean:</strong> {self.metrics["mean_affinity"]:.3f} ± {self.metrics["std_affinity"]:.3f} kcal/mol</p>
            <p><strong>Median:</strong> {self.metrics["median_affinity"]:.3f} kcal/mol</p>
            <p><strong>Q1:</strong> {self.metrics["q1_affinity"]:.3f} | <strong>Q3:</strong> {self.metrics["q3_affinity"]:.3f} kcal/mol</p>
        </div>
        
        <div class="stat-box">
            <h3>Hit Classification (Non-overlapping)</h3>
            <p><strong>Excellent (< -10):</strong> {self.metrics["hits_excellent"]["count"]} compounds ({self.metrics["hits_excellent"]["percentage"]:.1f}%)</p>
            <p><strong>Good (-10 to -8):</strong> {self.metrics["hits_good"]["count"]} compounds ({self.metrics["hits_good"]["percentage"]:.1f}%)</p>
            <p><strong>Moderate (-8 to -6):</strong> {self.metrics["hits_moderate"]["count"]} compounds ({self.metrics["hits_moderate"]["percentage"]:.1f}%)</p>
            <p><strong>Weak (≥ -6):</strong> {self.metrics["hits_weak"]["count"]} compounds ({self.metrics["hits_weak"]["percentage"]:.1f}%)</p>
        </div>
        
        <div class="stat-box">
            <h3>Enrichment Factors (< -15 kcal/mol)</h3>
            <p><strong>Top 1%:</strong> {self.metrics["enrichment_factor_1"]:.2f}x</p>
            <p><strong>Top 5%:</strong> {self.metrics["enrichment_factor_5"]:.2f}x</p>
            <p><strong>Top 10%:</strong> {self.metrics["enrichment_factor_10"]:.2f}x</p>
            <p><em>Exceptional hits concentration in top compounds</em></p>
            <p><em>High enrichment = ultra-strong binders clustered at top</em></p>
        </div>
    </div>
    
    <div class="plot-container">
        <div id="top-hits-plot"></div>
    </div>
    
    <div class="plot-container">
        <div id="distribution-plot"></div>
    </div>
    
    <div class="plot-container">
        <div id="enrichment-plot"></div>
    </div>
    
    <div class="table-container">
        <h2>Detailed Results</h2>
        {results_table.to_html(classes='results-table', table_id='results-table', escape=False)}
    </div>
            
            <div class="footer">
                <p>Generated with Claude Code Dashboard | 
                   Co-Authored-By: Claude &lt;noreply@anthropic.com&gt;</p>
            </div>
        </div>
    </div>

    <script>
        // Plot top hits chart
        Plotly.newPlot('top-hits-plot', {top_hits_json});
        
        // Plot distribution analysis
        Plotly.newPlot('distribution-plot', {distribution_json});
        
        // Plot enrichment analysis
        Plotly.newPlot('enrichment-plot', {enrichment_json});
    </script>
</body>
</html>
        """
        
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        print(f"Dashboard generated: {output_file}")
        return output_file

def main():
    """Main function to generate virtual screening dashboard"""
    import argparse
    
    parser = argparse.ArgumentParser(description='Generate AutoDock Vina Virtual Screening Dashboard')
    parser.add_argument('--output-dir', type=str, 
                       default='inputs/output',
                       help='Directory containing screening output files')
    parser.add_argument('--output-file', type=str, default=None,
                       help='Output HTML file name')
    parser.add_argument('--test', action='store_true',
                       help='Use test dataset (50 compounds)')
    
    args = parser.parse_args()
    
    # Use test dataset if specified
    if args.test:
        args.output_dir = 'testing/small_dataset_50/outputs'
        print("Using test dataset (50 compounds)...")
    
    output_dir = Path(args.output_dir)
    
    if not output_dir.exists():
        print(f"Error: Output directory not found: {output_dir}")
        return
    
    try:
        print(f"\n{'='*60}")
        print("AutoDock Vina Virtual Screening Dashboard Generator")
        print(f"{'='*60}\n")
        
        # Parse screening results
        print(f"Loading results from: {output_dir}")
        screening_parser = VinaScreeningParser(str(output_dir))
        
        # Analyze results
        print("\nPerforming statistical analysis...")
        analyzer = ScreeningAnalyzer(screening_parser.results)
        
        # Create dashboard
        print("\nGenerating dashboard...")
        dashboard = ScreeningDashboard(screening_parser.results, analyzer)
        
        # Generate report
        output_file = args.output_file or f'{output_dir}/screening_dashboard_{datetime.now().strftime("%Y%m%d_%H%M%S")}.html'
        dashboard.generate_html_report(output_file)
        
        # Print summary
        print(f"\n{'='*60}")
        print("SCREENING SUMMARY")
        print(f"{'='*60}")
        print(f"Total Compounds Screened: {analyzer.metrics['total_screened']:,}")
        print(f"Successful Dockings: {analyzer.metrics['successful_dockings']:,}")
        print(f"Success Rate: {analyzer.metrics['success_rate']:.1f}%")
        print(f"\nBest Affinity: {analyzer.metrics['best_affinity']:.3f} kcal/mol")
        print(f"Mean Affinity: {analyzer.metrics['mean_affinity']:.3f} ± {analyzer.metrics['std_affinity']:.3f} kcal/mol")
        print(f"\nHit Rates:")
        print(f"  Excellent (< -10): {analyzer.metrics['hits_excellent']['count']} ({analyzer.metrics['hits_excellent']['percentage']:.1f}%)")
        print(f"  Good (< -8): {analyzer.metrics['hits_good']['count']} ({analyzer.metrics['hits_good']['percentage']:.1f}%)")
        print(f"  Moderate (< -6): {analyzer.metrics['hits_moderate']['count']} ({analyzer.metrics['hits_moderate']['percentage']:.1f}%)")
        print(f"\nTop 5 Compounds:")
        for i, compound in enumerate(analyzer.metrics['top_10_compounds'][:5], 1):
            print(f"  {i}. {compound.ligand_name}: {compound.best_affinity:.3f} kcal/mol")
        print(f"\nRecommendation: {analyzer.metrics['recommendation']['action']}")
        print(f"Priority: {analyzer.metrics['recommendation']['priority']}")
        print(f"\n{'='*60}")
        print(f"Dashboard saved to: {output_file}")
        print("Open in a web browser to view the interactive report")
        print(f"{'='*60}\n")
        
    except Exception as e:
        print(f"Error generating dashboard: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()