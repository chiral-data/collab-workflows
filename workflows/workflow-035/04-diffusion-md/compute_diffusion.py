"""
Compute diffusion coefficient D from MSD via the Einstein relation.
D = lim_{t→∞} MSD(t) / (6t)

Fits the last half of the MSD where slope ≈ 1 on log-log (diffusive regime).
"""
import argparse, json, math, pathlib, sys

parser = argparse.ArgumentParser()
parser.add_argument("--msd",     required=True)
parser.add_argument("--report",  required=True)
parser.add_argument("--penetrant", default="O2")
parser.add_argument("--resin",   default="?")
parser.add_argument("--temp",    type=float, default=23.0)
args = parser.parse_args()


def read_xvg(path):
    times, vals = [], []
    for line in pathlib.Path(path).read_text().splitlines():
        if not line or line.startswith(("#", "@")):
            continue
        parts = line.split()
        times.append(float(parts[0]))
        vals.append(float(parts[1]))
    return times, vals


def linreg(x, y):
    n = len(x)
    sx, sy = sum(x), sum(y)
    sxx = sum(xi*xi for xi in x)
    sxy = sum(xi*yi for xi, yi in zip(x, y))
    d = n*sxx - sx*sx
    if abs(d) < 1e-30:
        return 0.0, sy/n
    return (n*sxy - sx*sy)/d, (sy - (n*sxy - sx*sy)/d * sx)/n


times_ps, msd_nm2 = read_xvg(args.msd)

if len(times_ps) < 10:
    print("WARNING: fewer than 10 MSD points — simulation too short", file=sys.stderr)
    D_cm2_s = None
    slope    = None
    regime   = "insufficient data"
else:
    # Use last 50% of trajectory (diffusive regime)
    n_fit = max(len(times_ps) // 2, 5)
    t_fit   = times_ps[-n_fit:]
    msd_fit = msd_nm2[-n_fit:]

    # Log-log slope to check if diffusive (should be ≈ 1)
    log_t   = [math.log(t) for t in t_fit if t > 0]
    log_msd = [math.log(m) for m in msd_fit if m > 0]
    slope, _ = linreg(log_t, log_msd) if len(log_t) > 2 else (None, None)

    # D from linear fit of MSD vs t: MSD = 6D·t
    # t in ps → s (*1e-12), MSD in nm² → cm² (*1e-14)
    slope_lin, _ = linreg(t_fit, msd_fit)
    D_nm2_ps  = slope_lin / 6.0
    D_cm2_s   = D_nm2_ps * 1e-14 / 1e-12   # nm²/ps → cm²/s = ×1e-14/1e-12 = ×0.01...
    # Actually: 1 nm² = 1e-14 cm², 1 ps = 1e-12 s → nm²/ps = 1e-14/1e-12 cm²/s = 1e-2 cm²/s
    # That's wrong. Let me be careful:
    # D [nm²/ps] × (1e-7 cm/nm)² / (1e-12 s/ps) = D × 1e-14 / 1e-12 = D × 1e-2  cm²/s
    # For typical polymers D ~ 1e-7 nm²/ps → D ~ 1e-9 cm²/s ✓
    D_cm2_s = D_nm2_ps * 1e-2

    if slope is not None and 0.7 <= slope <= 1.3:
        regime = "diffusive (slope ≈ 1)"
    elif slope is not None and slope < 0.7:
        regime = "sub-diffusive — run longer for reliable D"
    else:
        regime = "super-diffusive / ballistic — check simulation"

# Literature reference values (cm²/s at 23°C) for context
LIT_D = {
    ("PET",  "O2"):  3.4e-10,
    ("PET",  "H2O"): 5.0e-12,
    ("LDPE", "O2"):  4.5e-7,
    ("LDPE", "H2O"): 1.5e-8,
    ("PP",   "O2"):  2.0e-8,
    ("PP",   "H2O"): 3.0e-10,
    ("EVOH", "O2"):  2.0e-13,
    ("EVOH", "H2O"): 1.5e-11,
    ("PA6",  "O2"):  1.5e-9,
    ("PA6",  "H2O"): 5.0e-11,
}
lit = LIT_D.get((args.resin, args.penetrant))

report = {
    "resin_type":           args.resin,
    "penetrant":            args.penetrant,
    "temperature_c":        args.temp,
    "D_cm2_s":              D_cm2_s,
    "D_cm2_s_sci":          f"{D_cm2_s:.2e}" if D_cm2_s else "N/A",
    "log_log_slope":        round(slope, 3) if slope else None,
    "diffusive_regime":     regime,
    "n_msd_points":         len(times_ps),
    "sim_time_ps":          times_ps[-1] if times_ps else 0,
    "D_literature_cm2_s":   lit,
    "D_lit_sci":            f"{lit:.2e}" if lit else "N/A",
}
pathlib.Path(args.report).write_text(json.dumps(report, indent=2))

print(f"  D (Einstein) : {report['D_cm2_s_sci']} cm²/s")
print(f"  D (lit. ref) : {report['D_lit_sci']} cm²/s")
print(f"  Log-log slope: {slope:.2f} ({regime})" if slope else f"  Regime: {regime}")
