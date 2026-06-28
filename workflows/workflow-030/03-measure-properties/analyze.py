"""Extract density, Young's modulus, and specific modulus from GROMACS outputs."""
import json, os, pathlib, sys


def read_xvg(path):
    """Return (times, values) from a two-column XVG file."""
    times, vals = [], []
    for line in pathlib.Path(path).read_text().splitlines():
        if not line or line.startswith(("#", "@")):
            continue
        parts = line.split()
        times.append(float(parts[0]))
        vals.append(float(parts[1]))
    return times, vals


def linreg(x, y):
    """Return (slope, intercept) from a simple linear regression."""
    n = len(x)
    sx, sy = sum(x), sum(y)
    sxx = sum(xi * xi for xi in x)
    sxy = sum(xi * yi for xi, yi in zip(x, y))
    denom = n * sxx - sx * sx
    if abs(denom) < 1e-30:
        return 0.0, sy / n
    slope = (n * sxy - sx * sy) / denom
    intercept = (sy - slope * sx) / n
    return slope, intercept


def last_third_avg(vals):
    n = max(len(vals) // 3, 1)
    return sum(vals[-n:]) / n


# ── Load metadata ─────────────────────────────────────────────────────────────
build = json.loads(pathlib.Path("inputs/build_report.json").read_text())
equil = json.loads(pathlib.Path("inputs/equil_report.json").read_text())

target_temp_k = float(os.environ.get("TARGET_TEMP_K", equil["target_temp_k"]))
l0_nm = float(os.environ["L0_NM"])           # initial box x-length (nm)
deform_rate = float(os.environ["DEFORM_RATE"])  # nm/ps
deform_time = float(os.environ["DEFORM_TIME_PS"])
n_chains = build["n_chains"]

# ── Density from production run ───────────────────────────────────────────────
_, density_vals = read_xvg("outputs/density_prod.xvg")
avg_density_kg_m3 = last_third_avg(density_vals)
avg_density_gcc = avg_density_kg_m3 / 1000.0

# ── Young's modulus from deformation run ─────────────────────────────────────
times_ps, pxx_bar = read_xvg("outputs/stress_strain.xvg")

# Engineering strain ε(t) = deform_rate * t / L0  (dimensionless)
# Engineering stress σ_xx (Pa) = −Pres_XX * 1e5  (bar → Pa, sign convention)
strains  = [deform_rate * t / l0_nm for t in times_ps]
stresses_pa = [-p * 1e5 for p in pxx_bar]

# Fit only the elastic regime: strain 0 – 2%
elastic_pairs = [(e, s) for e, s in zip(strains, stresses_pa) if e <= 0.02]
if len(elastic_pairs) >= 10:
    e_vals, s_vals = zip(*elastic_pairs)
    slope_pa, _ = linreg(list(e_vals), list(s_vals))
    youngs_gpa = slope_pa / 1e9
else:
    youngs_gpa = None
    print("WARNING: fewer than 10 points in elastic regime — E not computed",
          file=sys.stderr)

# ── Specific modulus (EV lightweighting metric) ───────────────────────────────
# specific modulus = E / ρ  [GPa / (kg/m³)] = [kN·m/kg]
if youngs_gpa and avg_density_kg_m3 > 0:
    specific_modulus_kNm_kg = (youngs_gpa * 1e9) / avg_density_kg_m3 / 1e3
else:
    specific_modulus_kNm_kg = None

# ── Write report ──────────────────────────────────────────────────────────────
report = {
    "resin_type":              build["resin_type"],
    "temperature_c":           target_temp_k - 273.15,
    "fiber_loading_wt_pct":    build["fiber_loading_wt_pct"],
    "avg_density_kg_m3":       round(avg_density_kg_m3, 1),
    "avg_density_gcc":         round(avg_density_gcc, 4),
    "youngs_modulus_gpa":      round(youngs_gpa, 3) if youngs_gpa else None,
    "specific_modulus_kNm_kg": round(specific_modulus_kNm_kg, 1) if specific_modulus_kNm_kg else None,
    "equil_density_kg_m3":     equil["avg_density_kg_m3"],
    "n_chains":                n_chains,
    "prod_time_ps":            float(os.environ["PROD_TIME_PS"]),
    "deform_time_ps":          deform_time,
    "l0_nm":                   round(l0_nm, 4),
}

pathlib.Path("outputs/properties.json").write_text(json.dumps(report, indent=2))

print(f"  Resin          : {report['resin_type']} @ {report['temperature_c']:.0f} °C")
print(f"  Density        : {avg_density_kg_m3:.1f} kg/m³  ({avg_density_gcc:.4f} g/cc)")
if youngs_gpa:
    print(f"  Young's modulus: {youngs_gpa:.3f} GPa")
if specific_modulus_kNm_kg:
    print(f"  Specific modulus: {specific_modulus_kNm_kg:.1f} kN·m/kg")
