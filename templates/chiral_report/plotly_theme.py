"""Chiral Report — unified Plotly layout and color palette."""

from . import theme

COLORWAY = theme.CHART_COLORS

LAYOUT_DEFAULTS = {
    "template": "plotly_white",
    "font": {
        "family": "-apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif",
        "size": 13,
        "color": theme.SLATE_700,
    },
    "colorway": COLORWAY,
    "margin": {"t": 32, "r": 20, "b": 50, "l": 60},
    "paper_bgcolor": "rgba(0,0,0,0)",
    "plot_bgcolor": "rgba(0,0,0,0)",
    "hoverlabel": {
        "bgcolor": theme.WHITE,
        "bordercolor": theme.SLATE_200,
        "font": {"size": 12, "color": theme.SLATE_900},
    },
    "xaxis": {
        "gridcolor": "#eef1f5",
        "zerolinecolor": "#dde1e7",
        "title": {"font": {"size": 12, "color": theme.SLATE_600}},
    },
    "yaxis": {
        "gridcolor": "#eef1f5",
        "zerolinecolor": "#dde1e7",
        "title": {"font": {"size": 12, "color": theme.SLATE_600}},
    },
}


def apply_theme(layout: dict) -> dict:
    """Deep-merge LAYOUT_DEFAULTS under *layout* (user values win)."""
    import copy
    merged = copy.deepcopy(LAYOUT_DEFAULTS)
    _deep_update(merged, layout)
    return merged


def _deep_update(base: dict, override: dict) -> None:
    for k, v in override.items():
        if isinstance(v, dict) and isinstance(base.get(k), dict):
            _deep_update(base[k], v)
        else:
            base[k] = v
