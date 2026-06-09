"""
chiral_report — Unified scientific report template for Chiral workflows.

Usage::

    from chiral_report import Report

    report = Report("My Analysis", subtitle="Target: EGFR")
    report.stat_cards([{"value": "42", "label": "Compounds"}])
    report.section("Results", fig, df, "Some interpretation text.")
    report.save("report.html")
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path
from typing import Any

from . import components as C
from . import theme


class RawHTML:
    """Wrapper to pass pre-built HTML through render_item unchanged."""

    def __init__(self, html: str):
        self._html = html

    def __html__(self) -> str:
        return self._html


class Report:
    """Build a self-contained HTML report with Chiral branding."""

    def __init__(
        self,
        title: str,
        *,
        subtitle: str = "",
        timestamp: str | None = None,
    ):
        self.title = title
        self.subtitle = subtitle
        self.timestamp = timestamp or datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        self._blocks: list[str] = []
        self._extra_js: list[str] = []
        self._need_molstar = False

        self._blocks.append(
            C.hero(title, subtitle, timestamp=self.timestamp)
        )

    # ---- Public API ----

    def stat_cards(self, cards: list[dict]) -> "Report":
        """Add a row of summary stat cards.

        Each dict: {"value": str, "label": str, "color": str (optional)}
        color: "primary", "success", "warning", "danger", "info", "muted", or hex
        """
        self._blocks.append(C.stat_cards(cards))
        return self

    def section(
        self,
        title: str,
        *items: Any,
        note: str = "",
        flush: bool = False,
        molstar_viewer: list[dict] | None = None,
    ) -> "Report":
        """Add a content section.

        *items: Plotly Figure, DataFrame, str, or RawHTML — auto-detected.
        molstar_viewer: list of structure dicts for Mol* 3D viewer.
        """
        body_parts = [C.render_item(item) for item in items]

        if molstar_viewer is not None:
            self._need_molstar = True
            viewer_id = f"cr-molstar-{C._uid()}"
            body_parts.append(
                C.molstar_viewer(viewer_id, molstar_viewer)
            )

        body = "\n".join(body_parts)
        self._blocks.append(C.section(title, body, note=note, flush=flush))
        return self

    def side_by_side_sections(
        self,
        left_title: str,
        left_body: str,
        right_title: str,
        right_body: str,
    ) -> "Report":
        """Add two section cards side by side."""
        left = C.section(left_title, left_body)
        right = C.section(right_title, right_body)
        self._blocks.append(C.side_by_side(left, right))
        return self

    def table(
        self,
        title: str,
        headers: list[str],
        rows: list[list],
        *,
        sortable: bool = True,
        note: str = "",
        color_rules: dict | None = None,
    ) -> "Report":
        """Add a standalone data table section."""
        tbl = C.data_table(
            headers, rows,
            sortable=sortable,
            color_rules=color_rules,
        )
        self._blocks.append(
            C.section(title, tbl, note=note, flush=True)
        )
        return self

    def methods(self, sections: list[dict]) -> "Report":
        """Add a methods/parameters section.

        Each dict: {"title": str, "badge": str, "badge_color": str,
                     "rows": [(label, value), ...]}
        """
        body = C.methods_section(sections)
        self._blocks.append(C.section("Methods", body))
        return self

    def chart(
        self,
        title: str,
        charts: list[dict],
        *,
        note: str = "",
        height: int = 380,
    ) -> "Report":
        """Add a section with Plotly charts (raw trace/layout dicts).

        charts: [{"div_id": str, "traces": list, "layout": dict}]
        """
        body_parts = [C.plotly_chart(c["div_id"], height=height) for c in charts]
        js = C.plotly_script(charts)
        self._extra_js.append(js)
        body = "\n".join(body_parts)
        self._blocks.append(C.section(title, body, note=note))
        return self

    def raw(self, html: str) -> "Report":
        """Inject raw HTML directly."""
        self._blocks.append(html)
        return self

    def save(self, path: str | Path) -> Path:
        """Write the report as a self-contained HTML file."""
        self._blocks.append(C.footer(self.timestamp))

        body = "\n".join(self._blocks)
        html = C.page_html(
            self.title,
            body,
            include_plotly=True,
            include_molstar=self._need_molstar,
            extra_js="\n".join(self._extra_js),
        )

        out = Path(path)
        out.write_text(html, encoding="utf-8")
        return out
