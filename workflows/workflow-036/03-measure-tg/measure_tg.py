"""Estimate Tg, thermal expansion, and dimensional stability from the
melt-quench cooling series (node 02).

Method: fit a bilinear (two-segment) line to specific volume vs. temperature.
Tg is the intersection of the rubbery (high-T) and glassy (low-T) fitted
lines — the standard specific-volume kink method. With the fixed 5-point
ladder (melt, 200, 150, 80, 25 C) there are only two candidate breakpoints
that leave >=2 points on each side; this is the minimum viable case for a
kink fit and is documented as such in the report.
"""
import json
import os
import pathlib


def linreg(x, y):
    """Return (slope, intercept, sse) from a simple linear regression."""
    n = len(x)
    sx, sy = sum(x), sum(y)
    sxx = sum(xi * xi for xi in x)
    sxy = sum(xi * yi for xi, yi in zip(x, y))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        intercept = sy / n
        return 0.0, intercept, sum((yi - intercept) ** 2 for yi in y)
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    sse = sum((yi - (slope * xi + intercept)) ** 2 for xi, yi in zip(x, y))
    return slope, intercept, sse


def bilinear_tg_fit(points):
    """points: list of (temp_c, specific_volume) sorted ascending by temp_c.

    Returns dict with tg_c, cte_glassy, cte_rubbery, and per-segment fit
    diagnostics, or None if fewer than 4 points (can't form 2 valid segments).
    """
    n = len(points)
    if n < 4:
        return None

    temps = [p[0] for p in points]
    vols  = [p[1] for p in points]

    best = None
    for split in range(2, n - 1):  # split points into [0:split] glassy, [split:] rubbery
        low_x, low_y   = temps[:split], vols[:split]
        high_x, high_y = temps[split:], vols[split:]

        slope_g, intercept_g, sse_g = linreg(low_x, low_y)
        slope_r, intercept_r, sse_r = linreg(high_x, high_y)
        total_sse = sse_g + sse_r

        if best is None or total_sse < best["total_sse"]:
            best = {
                "split": split,
                "slope_glassy": slope_g, "intercept_glassy": intercept_g,
                "slope_rubbery": slope_r, "intercept_rubbery": intercept_r,
                "n_glassy": len(low_x), "n_rubbery": len(high_x),
                "total_sse": total_sse,
            }

    slope_diff = best["slope_glassy"] - best["slope_rubbery"]
    if abs(slope_diff) < 1e-12:
        tg_c = None
    else:
        tg_c = (best["intercept_rubbery"] - best["intercept_glassy"]) / slope_diff

    # The two segments' slopes come from short, noisy default MD runs and can
    # end up nearly parallel — dividing by that tiny slope_diff then
    # extrapolates the "intersection" far outside the sampled temperature
    # range (e.g. -1000+ C). Flag rather than silently report a fit that
    # isn't really constrained by the data.
    margin = 0.5 * (max(temps) - min(temps))
    tg_in_range = (
        tg_c is not None
        and (min(temps) - margin) <= tg_c <= (max(temps) + margin)
    )

    mean_v_glassy  = sum(vols[:best["split"]]) / best["n_glassy"]
    mean_v_rubbery = sum(vols[best["split"]:]) / best["n_rubbery"]

    cte_glassy  = best["slope_glassy"]  / mean_v_glassy  if mean_v_glassy  else None
    cte_rubbery = best["slope_rubbery"] / mean_v_rubbery if mean_v_rubbery else None

    return {
        "tg_c":               round(tg_c, 1) if tg_c is not None else None,
        "tg_reliable":        tg_in_range,
        "cte_glassy_per_c":   round(cte_glassy, 6)  if cte_glassy  is not None else None,
        "cte_rubbery_per_c":  round(cte_rubbery, 6) if cte_rubbery is not None else None,
        "fit_split_n_glassy":  best["n_glassy"],
        "fit_split_n_rubbery": best["n_rubbery"],
        "fit_total_sse":       round(best["total_sse"], 8),
    }


def main():
    cooling = json.loads(pathlib.Path("inputs/cooling_series.json").read_text())
    build   = json.loads(pathlib.Path("inputs/build_report.json").read_text())

    selected_temp_c = float(os.environ.get("PARAM_TEMPERATURE", "150.0"))

    series = sorted(cooling["series"], key=lambda s: s["temp_c"])
    points = [(s["temp_c"], s["avg_specific_volume_cm3_g"]) for s in series]

    fit = bilinear_tg_fit(points)
    if fit is None:
        print("WARNING: fewer than 4 series points — cannot fit Tg", flush=True)
        fit = {"tg_c": None, "tg_reliable": False, "cte_glassy_per_c": None, "cte_rubbery_per_c": None,
               "fit_split_n_glassy": None, "fit_split_n_rubbery": None, "fit_total_sse": None}
    elif not fit["tg_reliable"]:
        print(f"WARNING: fitted Tg ({fit['tg_c']} C) falls far outside the sampled "
              f"temperature range — glassy/rubbery slopes are nearly parallel "
              f"(noisy/under-converged segments). Treat as unreliable.", flush=True)

    # ── dimensional stability: % volume change vs 25 C reference ────────────────
    by_temp = {round(s["temp_c"], 1): s for s in series}
    ref = by_temp.get(25.0)
    # pick the series point closest to the student's selected temperature
    selected = min(series, key=lambda s: abs(s["temp_c"] - selected_temp_c))

    if ref and ref["avg_specific_volume_cm3_g"]:
        volume_change_pct = (
            (selected["avg_specific_volume_cm3_g"] - ref["avg_specific_volume_cm3_g"])
            / ref["avg_specific_volume_cm3_g"] * 100.0
        )
    else:
        volume_change_pct = None

    report = {
        "resin_type":               build["resin_type"],
        "crystallinity":            build["crystallinity"],
        "selected_temperature_c":   selected_temp_c,
        "matched_series_temp_c":    selected["temp_c"],
        "literature_tg_c":          build.get("literature_tg_c"),
        "tg_c":                     fit["tg_c"],
        "tg_reliable":              fit["tg_reliable"],
        "cte_glassy_per_c":         fit["cte_glassy_per_c"],
        "cte_rubbery_per_c":        fit["cte_rubbery_per_c"],
        "fit_split_n_glassy":       fit["fit_split_n_glassy"],
        "fit_split_n_rubbery":      fit["fit_split_n_rubbery"],
        "fit_total_sse":            fit["fit_total_sse"],
        "volume_change_pct":        round(volume_change_pct, 3) if volume_change_pct is not None else None,
        "reference_temp_c":         25.0,
        "density_temp_series":      series,
        "fit_note":                 ("Bilinear kink fit over the 5-point melt-quench ladder — "
                                      "the minimum viable case for a two-segment fit. "
                                      "Increase stage_time_ps for a better-converged curve."
                                      if fit["tg_reliable"] else
                                      "Glassy/rubbery segment slopes are nearly parallel, so the "
                                      "fitted intersection falls far outside the sampled temperature "
                                      "range — not a physically meaningful Tg. Increase stage_time_ps "
                                      "for a better-converged, more separable curve."),
    }

    outdir = pathlib.Path("outputs")
    outdir.mkdir(exist_ok=True)
    (outdir / "tg_report.json").write_text(json.dumps(report, indent=2))

    print(f"  Resin           : {report['resin_type']} ({report['crystallinity']} crystallinity)")
    reliability = "" if report["tg_reliable"] else "  [UNRELIABLE — see fit_note]"
    print(f"  Tg (estimated)  : {report['tg_c']} C{reliability}   (literature: {report['literature_tg_c']} C)")
    print(f"  CTE glassy      : {report['cte_glassy_per_c']} /C")
    print(f"  CTE rubbery     : {report['cte_rubbery_per_c']} /C")
    print(f"  Volume change   : {report['volume_change_pct']}% at {report['matched_series_temp_c']} C vs 25 C")


if __name__ == "__main__":
    main()
