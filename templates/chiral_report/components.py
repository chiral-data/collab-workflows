"""Chiral Report — HTML component renderers (internal)."""

from __future__ import annotations

import html as _html
import json
import uuid
from typing import Any

from . import theme
from .assets import LOGO_HERO, LOGO_FOOTER
from .plotly_theme import apply_theme


# -------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------

def _esc(text: str) -> str:
    return _html.escape(str(text))


def _uid() -> str:
    return uuid.uuid4().hex[:8]


def _color_class(color: str | None) -> str:
    mapping = {
        "primary": "cr-color-primary",
        "brand": "cr-color-brand",
        "success": "cr-color-success",
        "warning": "cr-color-warning",
        "danger": "cr-color-danger",
        "info": "cr-color-info",
        "muted": "cr-color-muted",
    }
    if color in mapping:
        return mapping[color]
    return ""


def _color_style(color: str | None) -> str:
    mapping = {
        "primary": theme.PRIMARY,
        "brand": theme.BRAND,
        "success": theme.SUCCESS,
        "warning": theme.WARNING,
        "danger": theme.DANGER,
        "info": theme.INFO,
        "muted": theme.SLATE_400,
    }
    if color in mapping:
        return f' style="color:{mapping[color]}"'
    if color and color.startswith("#"):
        return f' style="color:{_esc(color)}"'
    return ""


# -------------------------------------------------------------------
# Page wrapper
# -------------------------------------------------------------------

def page_html(
    title: str,
    body: str,
    *,
    include_plotly: bool = True,
    include_molstar: bool = False,
    extra_css: str = "",
    extra_js: str = "",
) -> str:
    parts = [
        "<!DOCTYPE html>",
        '<html lang="en">',
        "<head>",
        '  <meta charset="UTF-8">',
        '  <meta name="viewport" content="width=device-width, initial-scale=1.0">',
        f"  <title>{_esc(title)}</title>",
        f'  <link href="{theme.BOOTSTRAP_CSS}" rel="stylesheet">',
    ]
    if include_molstar:
        parts.append(f'  <link rel="stylesheet" href="{theme.MOLSTAR_CSS}">')
    if include_plotly:
        parts.append(f'  <script src="{theme.PLOTLY_JS}"></script>')
    if include_molstar:
        parts.append(f'  <script src="{theme.MOLSTAR_JS}"></script>')
    parts += [
        "  <style>",
        theme.CSS,
        extra_css,
        "  </style>",
        "</head>",
        "<body>",
        '<div class="cr-page">',
        body,
        "</div>",  # .cr-page
    ]
    # Table sort JS
    parts += [
        "<script>",
        theme.TABLE_SORT_JS,
        extra_js,
        "</script>",
        "</body>",
        "</html>",
    ]
    return "\n".join(parts)


# -------------------------------------------------------------------
# Hero
# -------------------------------------------------------------------

def hero(
    title: str,
    subtitle: str = "",
    *,
    timestamp: str = "",
    show_logo: bool = True,
) -> str:
    brand_bar = ""
    if show_logo:
        brand_bar = f"""
<div class="cr-hero-brand">
  <img src="{LOGO_HERO}" alt="Chiral">
  <span>Produced by CHIRAL</span>
</div>
"""

    ts = ""
    if timestamp:
        ts = f'<div class="cr-timestamp">{_esc(timestamp)}</div>'

    sub = ""
    if subtitle:
        sub = f'<p class="cr-sub">{_esc(subtitle)}</p>'

    return f"""
{brand_bar}
<div class="cr-hero">
  <h1>{_esc(title)}</h1>
  {sub}
  {ts}
</div>
"""


# -------------------------------------------------------------------
# Stat cards
# -------------------------------------------------------------------

def stat_cards(cards: list[dict]) -> str:
    n = len(cards)
    col_class = "col-6 col-md-3" if n == 4 else f"col-6 col-md-{12 // min(n, 6)}"

    items = []
    for c in cards:
        val = _esc(str(c["value"]))
        label = _esc(str(c["label"]))
        color_attr = _color_style(c.get("color"))
        items.append(f"""
      <div class="{col_class}">
        <div class="cr-stat-card">
          <div class="cr-stat-value"{color_attr}>{val}</div>
          <div class="cr-stat-label">{label}</div>
        </div>
      </div>""")

    return f"""
<div class="cr-stats">
  <div class="row g-3">
    {"".join(items)}
  </div>
</div>
"""


# -------------------------------------------------------------------
# Section
# -------------------------------------------------------------------

def section(title: str, body: str, *, note: str = "", flush: bool = False) -> str:
    note_html = f'<p class="cr-note">{note}</p>' if note else ""
    flush_cls = " cr-flush" if flush else ""
    return f"""
<div class="cr-section">
  <div class="cr-section-head">{title}</div>
  <div class="cr-section-body{flush_cls}">
    {body}
    {note_html}
  </div>
</div>
"""


# -------------------------------------------------------------------
# Data table
# -------------------------------------------------------------------

def data_table(
    headers: list[str],
    rows: list[list],
    *,
    table_id: str = "",
    sortable: bool = False,
    color_rules: dict[int, dict] | None = None,
) -> str:
    tid = f' id="{_esc(table_id)}"' if table_id else ""
    sort_cls = " cr-sortable" if sortable else ""
    arrow = ' <span class="cr-sort-arrow">&#x25B4;</span>' if sortable else ""

    ths = "".join(
        f"<th>{_esc(h)}{arrow}</th>" for h in headers
    )

    row_htmls = []
    for row in rows:
        cells = []
        for i, val in enumerate(row):
            cls_parts = []
            style = ""
            val_str = _esc(str(val)) if val is not None else "&mdash;"

            # first col is name
            if i == 0:
                cls_parts.append("cr-cell-name")
            else:
                cls_parts.append("cr-cell-num")

            # color rules: {col_index: {threshold_good: float, threshold_warn: float, lower_is_better: bool}}
            if color_rules and i in color_rules:
                rule = color_rules[i]
                try:
                    fval = float(val)
                    low_better = rule.get("lower_is_better", True)
                    good = rule.get("good")
                    warn = rule.get("warn")
                    if good is not None and warn is not None:
                        if low_better:
                            if fval <= good:
                                cls_parts.append("cr-cell-good")
                            elif fval <= warn:
                                cls_parts.append("cr-cell-warn")
                            else:
                                cls_parts.append("cr-cell-bad")
                        else:
                            if fval >= good:
                                cls_parts.append("cr-cell-good")
                            elif fval >= warn:
                                cls_parts.append("cr-cell-warn")
                            else:
                                cls_parts.append("cr-cell-bad")
                except (ValueError, TypeError):
                    pass

            cls = f' class="{" ".join(cls_parts)}"' if cls_parts else ""
            cells.append(f"<td{cls}>{val_str}</td>")
        row_htmls.append(f"<tr>{''.join(cells)}</tr>")

    return f"""
<div class="cr-table-wrap">
  <table class="cr-table{sort_cls}"{tid}>
    <thead><tr>{ths}</tr></thead>
    <tbody>
      {"".join(row_htmls)}
    </tbody>
  </table>
</div>
"""


# -------------------------------------------------------------------
# Plotly chart
# -------------------------------------------------------------------

def plotly_chart(div_id: str, *, height: int = 380) -> str:
    return f'<div id="{_esc(div_id)}" class="cr-chart" style="height:{height}px;"></div>'


def plotly_script(charts: list[dict]) -> str:
    """Generate JS for multiple Plotly charts.

    Each dict: {"div_id": str, "traces": list, "layout": dict}
    Layout is auto-merged with the Chiral theme.
    """
    blocks = []
    for chart in charts:
        div_id = chart["div_id"]
        traces = json.dumps(chart["traces"], default=str)
        layout = json.dumps(apply_theme(chart.get("layout", {})), default=str)
        config = json.dumps({"responsive": True, "displayModeBar": False})
        blocks.append(f"Plotly.newPlot('{div_id}', {traces}, {layout}, {config});")
    return "\n".join(blocks)


# -------------------------------------------------------------------
# Mol* viewer
# -------------------------------------------------------------------

def _escape_for_js(s: str) -> str:
    return s.replace("\\", "\\\\").replace("`", "\\`").replace("${", "\\${")


def molstar_viewer(
    viewer_id: str,
    structures: list[dict],
    *,
    height: int = 520,
) -> str:
    """Mol* 3D viewer with toggle buttons.

    Each structure: {"label": str, "data": str, "format": str, "color": str}
    format: "pdb" or "sdf"
    """
    buttons = []
    for i, s in enumerate(structures):
        color = _esc(s.get("color", theme.SLATE_400))
        label = _esc(s["label"])
        buttons.append(
            f'<button id="toggle-struct-{i}" onclick="toggleStruct({i})" '
            f'class="cr-toggle-btn" style="background:{color};">{label}</button>'
        )

    structs_js = []
    for s in structures:
        data_escaped = _escape_for_js(s["data"])
        structs_js.append(
            f'{{label:`{_escape_for_js(s["label"])}`, '
            f'data:`{data_escaped}`, '
            f'format:"{s["format"]}"}}'
        )

    js = f"""
(function() {{
  var structs = [{",".join(structs_js)}];
  var visible = new Array(structs.length).fill(true);
  var molPlugin = null;

  window.toggleStruct = function(idx) {{
    if (!molPlugin) return;
    visible[idx] = !visible[idx];
    var show = visible[idx];
    var hierarchy = molPlugin.managers.structure.hierarchy.current.structures;
    if (!hierarchy[idx]) return;
    var isHidden = !!hierarchy[idx].cell.state.isHidden;
    if ((show && isHidden) || (!show && !isHidden)) {{
      try {{
        molPlugin.managers.structure.hierarchy.toggleVisibility([hierarchy[idx]]);
      }} catch(e) {{
        try {{
          molstar.PluginCommands.State.ToggleVisibility(molPlugin, {{
            state: molPlugin.state.data,
            ref: hierarchy[idx].cell.transform.ref
          }});
        }} catch(e2) {{ console.warn('toggle failed:', e2); }}
      }}
    }}
    var btn = document.getElementById('toggle-struct-' + idx);
    if (btn) btn.style.opacity = show ? '1' : '0.35';
  }};

  molstar.Viewer.create('{viewer_id}', {{
    layoutIsExpanded: false,
    layoutShowControls: false,
    layoutShowLeftPanel: true,
    layoutShowSequence: false,
    layoutShowLog: false,
    layoutShowRemoteState: false,
    viewportShowAnimation: false,
    viewportShowExpand: true,
    viewportShowSelectionMode: false
  }}).then(function(viewer) {{
    molPlugin = viewer.plugin;
    var p = molPlugin;
    var chain = Promise.resolve();
    structs.forEach(function(s) {{
      chain = chain
        .then(function() {{ return p.builders.data.rawData({{data: s.data, label: s.label}}); }})
        .then(function(d) {{ return p.builders.structure.parseTrajectory(d, s.format); }})
        .then(function(t) {{ return p.builders.structure.hierarchy.applyPreset(t, 'default'); }});
    }});
    chain.catch(function(e) {{
      document.getElementById('{viewer_id}').innerHTML =
        '<p style="color:red;padding:16px;">Mol* error: ' + e + '</p>';
    }});
  }}).catch(function(e) {{
    document.getElementById('{viewer_id}').innerHTML =
      '<p style="color:red;padding:16px;">Mol* failed to initialize: ' + e + '</p>';
  }});
}})();
"""

    return f"""
<div class="cr-toggle-bar">
  {"".join(buttons)}
</div>
<div id="{_esc(viewer_id)}" class="cr-molstar-wrap" style="height:{height}px;"></div>
<script>{js}</script>
"""


# -------------------------------------------------------------------
# Methods table
# -------------------------------------------------------------------

def methods_section(sections: list[dict]) -> str:
    """Render methods/parameters tables.

    Each dict: {"title": str, "badge": str (optional), "badge_color": str (optional),
                "rows": [(label, value), ...]}
    """
    parts = []
    for sec in sections:
        badge = ""
        if sec.get("badge"):
            color = sec.get("badge_color", "primary")
            if color.startswith("#"):
                badge = f'<span class="cr-badge me-2" style="background:{_esc(color)};color:#fff">{_esc(sec["badge"])}</span>'
            else:
                badge = f'<span class="cr-badge cr-badge-{_esc(color)} me-2">{_esc(sec["badge"])}</span>'

        rows = "".join(
            f"<tr><td>{_esc(label)}</td><td>{value}</td></tr>"
            for label, value in sec["rows"]
        )
        parts.append(f"""
      <div class="col-12 col-md-6">
        <h6>{badge}{_esc(sec["title"])}</h6>
        <table class="cr-table"><tbody>{rows}</tbody></table>
      </div>""")

    return f"""
<div class="cr-methods">
  <div class="row g-4">
    {"".join(parts)}
  </div>
</div>
"""


# -------------------------------------------------------------------
# Side-by-side
# -------------------------------------------------------------------

def side_by_side(left: str, right: str) -> str:
    return f"""
<div class="row g-3">
  <div class="col-12 col-md-6">{left}</div>
  <div class="col-12 col-md-6">{right}</div>
</div>
"""


# -------------------------------------------------------------------
# Badge
# -------------------------------------------------------------------

def badge(text: str, color: str = "primary") -> str:
    return f'<span class="cr-badge cr-badge-{_esc(color)}">{_esc(text)}</span>'


# -------------------------------------------------------------------
# Footer
# -------------------------------------------------------------------

def footer(timestamp: str = "") -> str:
    ts = ""
    if timestamp:
        ts = f'{_esc(timestamp)} · '
    return f"""
<div class="cr-footer">
  <div class="cr-footer-meta">
    {ts}Generated with chiral_report
  </div>
</div>
"""


# -------------------------------------------------------------------
# Auto-render helpers (for Report class)
# -------------------------------------------------------------------

def render_item(item: Any) -> str:
    """Auto-detect type and render to HTML."""
    # Plotly Figure
    if _is_plotly_figure(item):
        return _render_plotly(item)
    # pandas DataFrame
    if _is_dataframe(item):
        return _render_dataframe(item)
    # raw HTML (marked with __html__)
    if hasattr(item, "__html__"):
        return item.__html__()
    # string → paragraph
    if isinstance(item, str):
        return f'<p class="cr-text">{_esc(item)}</p>'
    return f'<p class="cr-text">{_esc(str(item))}</p>'


def _is_plotly_figure(obj: Any) -> bool:
    t = type(obj)
    return (
        t.__module__.startswith("plotly.") and "Figure" in t.__name__
    ) if hasattr(t, "__module__") else False


def _is_dataframe(obj: Any) -> bool:
    t = type(obj)
    return (
        t.__module__.startswith("pandas") and t.__name__ == "DataFrame"
    ) if hasattr(t, "__module__") else False


def _render_plotly(fig: Any) -> str:
    div_id = f"cr-fig-{_uid()}"
    traces = fig.to_dict()["data"]
    layout = fig.to_dict().get("layout", {})
    chart_html = plotly_chart(div_id)
    script = plotly_script([{"div_id": div_id, "traces": traces, "layout": layout}])
    return f"{chart_html}\n<script>{script}</script>"


def _render_dataframe(df: Any) -> str:
    headers = list(df.columns)
    rows = df.values.tolist()
    return data_table(headers, rows, sortable=True)
