"""
Chiral Report — CSS theme and color constants.

Design language: "Laboratory Precision"
Primary: Nature/Springer science blue (#025e8d)
Brand accent: Chiral orange (#ff6b35) — logo area only
"""

# ---------------------------------------------------------------------------
# Color tokens
# ---------------------------------------------------------------------------
PRIMARY = "#025e8d"
PRIMARY_DARK = "#01324b"
PRIMARY_LIGHT = "#e8f4f8"
BRAND = "#ff6b35"
BRAND_DARK = "#e55a2b"

SLATE_900 = "#0f172a"
SLATE_700 = "#334155"
SLATE_600 = "#475569"
SLATE_400 = "#94a3b8"
SLATE_200 = "#e2e8f0"
SLATE_100 = "#f1f5f9"

BG = "#f8f9fa"
WHITE = "#ffffff"

SUCCESS = "#16a34a"
WARNING = "#d97706"
DANGER = "#dc2626"
INFO = "#0891b2"

# Okabe-Ito colorblind-safe palette (Nature-recommended)
CHART_COLORS = [
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # bluish green
    "#E69F00",  # orange
    "#56B4E9",  # sky blue
    "#CC79A7",  # reddish purple
    "#F0E442",  # yellow
]

# ---------------------------------------------------------------------------
# CDN dependencies
# ---------------------------------------------------------------------------
BOOTSTRAP_CSS = "https://cdn.jsdelivr.net/npm/bootstrap@5.3.3/dist/css/bootstrap.min.css"
PLOTLY_JS = "https://cdn.plot.ly/plotly-2.35.2.min.js"
MOLSTAR_CSS = "https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.css"
MOLSTAR_JS = "https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.js"

# ---------------------------------------------------------------------------
# Master stylesheet
# ---------------------------------------------------------------------------
CSS = """
/* ===================================================================
   Chiral Report Theme — Laboratory Precision
   Science-blue foundation · brand-orange whisper · austere clarity
   =================================================================== */

:root {
  --cr-primary:       #025e8d;
  --cr-primary-dark:  #01324b;
  --cr-primary-mid:   #0284c7;
  --cr-primary-light: #e8f4f8;
  --cr-brand:         #ff6b35;
  --cr-brand-dark:    #e55a2b;

  --cr-slate-900:     #0f172a;
  --cr-slate-700:     #334155;
  --cr-slate-600:     #475569;
  --cr-slate-400:     #94a3b8;
  --cr-slate-200:     #e2e8f0;
  --cr-slate-100:     #f1f5f9;

  --cr-bg:            #f8f9fa;
  --cr-white:         #ffffff;

  --cr-success:       #16a34a;
  --cr-warning:       #d97706;
  --cr-danger:        #dc2626;
  --cr-info:          #0891b2;

  --cr-radius:        12px;
  --cr-shadow:        0 1px 3px rgba(0,0,0,0.04), 0 4px 14px rgba(0,0,0,0.06);
  --cr-shadow-sm:     0 1px 2px rgba(0,0,0,0.05);
  --cr-font:          "Suisse Int'l", -apple-system, BlinkMacSystemFont,
                      "Segoe UI", Roboto, Oxygen, Ubuntu, Cantarell,
                      "Helvetica Neue", sans-serif;
}

/* ---- Base ---- */

*,
*::before,
*::after { box-sizing: border-box; }

body {
  background: var(--cr-bg);
  color: var(--cr-slate-900);
  font-family: var(--cr-font);
  font-size: 0.9375rem;
  line-height: 1.6;
  -webkit-font-smoothing: antialiased;
  -moz-osx-font-smoothing: grayscale;
}

.cr-page {
  max-width: 1280px;
  margin: 0 auto;
  padding: 24px 20px 48px;
}

/* ---- Brand bar (above hero) ---- */

.cr-hero-brand {
  display: flex;
  align-items: center;
  gap: 10px;
  margin-bottom: 16px;
}

.cr-hero-brand img {
  height: 36px;
  width: auto;
}

.cr-hero-brand span {
  font-size: 0.875rem;
  font-weight: 600;
  color: var(--cr-slate-600);
  letter-spacing: 0.04em;
  text-transform: uppercase;
}

/* ---- Hero header ---- */

.cr-hero {
  background: linear-gradient(135deg, #1a1a2e 0%, #2d2d3f 60%, #3a3a4a 100%);
  color: var(--cr-white);
  border-radius: var(--cr-radius);
  padding: 36px 40px 32px;
  margin-bottom: 28px;
}

.cr-hero h1 {
  font-size: 1.75rem;
  font-weight: 700;
  letter-spacing: -0.01em;
  margin: 0 0 6px;
  line-height: 1.25;
}

.cr-hero .cr-sub {
  font-size: 0.875rem;
  opacity: 0.72;
  margin: 0;
  font-weight: 400;
}

.cr-hero .cr-timestamp {
  font-size: 0.75rem;
  opacity: 0.5;
  margin-top: 10px;
  font-variant-numeric: tabular-nums;
}

/* ---- Stat cards ---- */

.cr-stats { margin-bottom: 28px; }

.cr-stat-card {
  background: var(--cr-white);
  border-radius: var(--cr-radius);
  padding: 22px 18px;
  text-align: center;
  box-shadow: var(--cr-shadow);
  height: 100%;
  border: 1px solid rgba(0,0,0,0.04);
  transition: box-shadow 0.15s ease;
}

.cr-stat-card:hover {
  box-shadow: 0 2px 6px rgba(0,0,0,0.06), 0 8px 24px rgba(0,0,0,0.08);
}

.cr-stat-value {
  font-size: 1.5rem;
  font-weight: 700;
  line-height: 1.15;
  font-variant-numeric: tabular-nums;
  letter-spacing: -0.02em;
}

.cr-stat-label {
  font-size: 0.75rem;
  text-transform: uppercase;
  letter-spacing: 0.05em;
  color: var(--cr-slate-400);
  margin-top: 6px;
  font-weight: 500;
}

/* ---- Section cards ---- */

.cr-section {
  background: var(--cr-white);
  border-radius: var(--cr-radius);
  box-shadow: var(--cr-shadow);
  margin-bottom: 24px;
  overflow: hidden;
  border: 1px solid rgba(0,0,0,0.04);
}

.cr-section-head {
  padding: 14px 24px;
  font-weight: 600;
  font-size: 1rem;
  color: var(--cr-slate-900);
  border-bottom: 1px solid var(--cr-slate-100);
  letter-spacing: -0.005em;
  display: flex;
  align-items: center;
  gap: 8px;
}

.cr-section-body {
  padding: 20px 24px;
}

.cr-section-body.cr-flush { padding: 0; }

.cr-section-body p:last-child { margin-bottom: 0; }

.cr-note {
  font-size: 0.8125rem;
  color: var(--cr-slate-400);
  margin-top: 10px;
  line-height: 1.65;
}

.cr-text {
  font-size: 0.875rem;
  color: var(--cr-slate-700);
  line-height: 1.6;
  margin-bottom: 14px;
}

/* ---- Data tables ---- */

.cr-table-wrap { overflow-x: auto; }

.cr-table {
  width: 100%;
  border-collapse: separate;
  border-spacing: 0;
  font-size: 0.875rem;
  font-variant-numeric: tabular-nums;
}

.cr-table thead th {
  font-size: 0.75rem;
  font-weight: 600;
  text-transform: uppercase;
  letter-spacing: 0.04em;
  color: var(--cr-slate-600);
  padding: 10px 14px;
  border-bottom: 2px solid var(--cr-slate-200);
  white-space: nowrap;
  position: sticky;
  top: 0;
  background: var(--cr-white);
  z-index: 1;
}

.cr-table.cr-sortable thead th {
  cursor: pointer;
  user-select: none;
  transition: color 0.12s;
}

.cr-table.cr-sortable thead th:hover {
  color: var(--cr-primary);
}

.cr-table.cr-sortable thead th .cr-sort-arrow {
  display: inline-block;
  margin-left: 4px;
  font-size: 0.75rem;
  opacity: 0.35;
  transition: opacity 0.12s;
}

.cr-table.cr-sortable thead th:hover .cr-sort-arrow { opacity: 0.7; }

.cr-table tbody td {
  padding: 9px 14px;
  border-bottom: 1px solid var(--cr-slate-100);
  vertical-align: middle;
  color: var(--cr-slate-700);
}

.cr-table tbody tr:last-child td { border-bottom: none; }

.cr-table tbody tr:hover td {
  background: var(--cr-primary-light);
}

.cr-table .cr-cell-name {
  font-weight: 500;
  color: var(--cr-slate-900);
  max-width: 220px;
  overflow: hidden;
  text-overflow: ellipsis;
  white-space: nowrap;
}

.cr-table .cr-cell-num {
  text-align: right;
  font-variant-numeric: tabular-nums;
}

.cr-table .cr-cell-good { color: var(--cr-success); font-weight: 600; }
.cr-table .cr-cell-warn { color: var(--cr-warning); font-weight: 600; }
.cr-table .cr-cell-bad  { color: var(--cr-danger);  font-weight: 600; }

/* ---- Plotly chart containers ---- */

.cr-chart {
  width: 100%;
  min-height: 320px;
}

/* ---- Badges / pills ---- */

.cr-badge {
  display: inline-block;
  padding: 2px 10px;
  border-radius: 10px;
  font-size: 0.75rem;
  font-weight: 600;
  letter-spacing: 0.03em;
  text-transform: uppercase;
  line-height: 1.6;
  vertical-align: middle;
}

.cr-badge-primary  { background: var(--cr-primary);     color: var(--cr-white); }
.cr-badge-brand    { background: var(--cr-brand);        color: var(--cr-white); }
.cr-badge-success  { background: var(--cr-success);      color: var(--cr-white); }
.cr-badge-warning  { background: var(--cr-warning);      color: var(--cr-white); }
.cr-badge-danger   { background: var(--cr-danger);       color: var(--cr-white); }
.cr-badge-info     { background: var(--cr-info);         color: var(--cr-white); }
.cr-badge-muted    { background: var(--cr-slate-400);    color: var(--cr-white); }

.cr-badge-outline {
  background: transparent;
  border: 1.5px solid currentColor;
}

/* ---- Mol* viewer ---- */

.cr-molstar-wrap {
  position: relative;
  width: 100%;
  height: 520px;
  border-radius: 8px;
  overflow: hidden;
  border: 1px solid var(--cr-slate-200);
  background: var(--cr-slate-100);
}

.cr-toggle-bar {
  display: flex;
  flex-wrap: wrap;
  gap: 6px;
  margin-bottom: 12px;
}

.cr-toggle-btn {
  background: var(--cr-slate-400);
  color: var(--cr-white);
  border: none;
  padding: 4px 14px;
  border-radius: 16px;
  font-size: 0.8125rem;
  font-weight: 500;
  cursor: pointer;
  transition: opacity 0.15s, background-color 0.15s;
}

.cr-toggle-btn:hover { filter: brightness(1.1); }
.cr-toggle-btn.cr-off { opacity: 0.4; }

/* ---- Methods section ---- */

.cr-methods h6 {
  font-size: 0.875rem;
  font-weight: 600;
  color: var(--cr-slate-900);
  margin-bottom: 10px;
  display: flex;
  align-items: center;
  gap: 8px;
}

.cr-methods .cr-table td:first-child {
  color: var(--cr-slate-400);
  white-space: nowrap;
  padding-right: 16px;
  width: 1%;
}

.cr-methods .cr-table td:last-child {
  color: var(--cr-slate-700);
}

.cr-methods code {
  font-size: 0.8em;
  background: var(--cr-slate-100);
  padding: 1px 5px;
  border-radius: 3px;
  color: var(--cr-primary);
}

/* ---- Side-by-side layout ---- */

.cr-side-by-side {
  display: grid;
  grid-template-columns: 1fr 1fr;
  gap: 20px;
}

@media (max-width: 768px) {
  .cr-side-by-side { grid-template-columns: 1fr; }
}

/* ---- Footer ---- */

.cr-footer {
  margin-top: 40px;
  padding: 24px 0 0;
  border-top: 1px solid var(--cr-slate-200);
  display: flex;
  align-items: center;
  justify-content: space-between;
  gap: 16px;
}

.cr-footer-brand {
  display: flex;
  align-items: center;
  gap: 10px;
}

.cr-footer-logo {
  height: 28px;
  width: auto;
  opacity: 0.7;
}

.cr-footer-name {
  font-size: 0.82rem;
  font-weight: 600;
  color: var(--cr-slate-600);
  letter-spacing: 0.04em;
  text-transform: uppercase;
}

.cr-footer-meta {
  font-size: 0.75rem;
  color: var(--cr-slate-400);
  text-align: right;
  font-variant-numeric: tabular-nums;
}

/* ---- Utility ---- */

.cr-color-primary { color: var(--cr-primary); }
.cr-color-brand   { color: var(--cr-brand); }
.cr-color-success { color: var(--cr-success); }
.cr-color-warning { color: var(--cr-warning); }
.cr-color-danger  { color: var(--cr-danger); }
.cr-color-info    { color: var(--cr-info); }
.cr-color-muted   { color: var(--cr-slate-400); }

/* ---- Print-friendly ---- */

@media print {
  body { background: white; }
  .cr-hero { break-inside: avoid; }
  .cr-section { break-inside: avoid; box-shadow: none; border: 1px solid #ddd; }
  .cr-stat-card { box-shadow: none; border: 1px solid #ddd; }
  .cr-footer { page-break-inside: avoid; }
}

/* ---- Responsive tweaks ---- */

@media (max-width: 576px) {
  .cr-hero { padding: 24px 20px 20px; }
  .cr-hero h1 { font-size: 1.3rem; }
  .cr-hero-logo { height: 32px; }
  .cr-section-body { padding: 16px; }
  .cr-stat-value { font-size: 1.5rem; }
}
"""

# ---------------------------------------------------------------------------
# Table-sorting JS (vanilla, minimal)
# ---------------------------------------------------------------------------
TABLE_SORT_JS = """
document.querySelectorAll('.cr-sortable').forEach(function(table) {
  var heads = table.querySelectorAll('thead th');
  heads.forEach(function(th, colIdx) {
    th.addEventListener('click', function() {
      var body = table.querySelector('tbody');
      var rows = Array.from(body.querySelectorAll('tr'));
      var asc = th.dataset.sortDir !== 'asc';
      heads.forEach(function(h) { h.dataset.sortDir = ''; });
      th.dataset.sortDir = asc ? 'asc' : 'desc';
      rows.sort(function(a, b) {
        var at = a.cells[colIdx].textContent.trim();
        var bt = b.cells[colIdx].textContent.trim();
        var an = parseFloat(at.replace(/[^\\d.\\-]/g, ''));
        var bn = parseFloat(bt.replace(/[^\\d.\\-]/g, ''));
        if (!isNaN(an) && !isNaN(bn)) return asc ? an - bn : bn - an;
        return asc ? at.localeCompare(bt) : bt.localeCompare(at);
      });
      rows.forEach(function(r) { body.appendChild(r); });
      heads.forEach(function(h) {
        var arrow = h.querySelector('.cr-sort-arrow');
        if (arrow) arrow.textContent = h === th ? (asc ? '\\u25B2' : '\\u25BC') : '\\u25B4';
      });
    });
  });
});
"""
