"""
Halpin-Tsai short-fiber composite correction.

Scales neat-resin MD properties (E, density) to glass-fiber loaded values.
Reference: Halpin & Tsai (1969); random-orientation average by Christensen & Waals (1972).
"""
import json, math, pathlib

# ── Glass fiber constants (E-glass, injection-moulded short fiber) ─────────────
E_FIBER_GPA    = 72.0     # E-glass longitudinal modulus
RHO_FIBER      = 2540.0   # kg/m³
ASPECT_RATIO   = 20.0     # l_f / d_f — typical for IM short GF


def weight_to_volume_fraction(w_f: float, rho_m: float) -> float:
    """Convert fiber weight fraction to volume fraction."""
    if w_f <= 0.0:
        return 0.0
    w_m = 1.0 - w_f
    return (w_f / RHO_FIBER) / (w_f / RHO_FIBER + w_m / rho_m)


def halpin_tsai_random(E_m_gpa: float, V_f: float, ar: float = ASPECT_RATIO) -> float:
    """
    Random short-fiber composite modulus [GPa] via Halpin-Tsai + Christensen-Waals.

    E_m_gpa : neat-resin modulus (GPa)
    V_f     : fiber volume fraction
    ar      : fiber aspect ratio l/d
    """
    if V_f <= 0.0:
        return E_m_gpa

    ratio = E_FIBER_GPA / E_m_gpa

    xi_L  = 2.0 * ar   # longitudinal reinforcement parameter
    xi_T  = 2.0         # transverse reinforcement parameter

    eta_L = (ratio - 1.0) / (ratio + xi_L)
    eta_T = (ratio - 1.0) / (ratio + xi_T)

    E_11 = E_m_gpa * (1.0 + xi_L * eta_L * V_f) / (1.0 - eta_L * V_f)
    E_22 = E_m_gpa * (1.0 + xi_T * eta_T * V_f) / (1.0 - eta_T * V_f)

    # Isotropic in-plane average (random short-fiber orientation)
    return 3.0 / 8.0 * E_11 + 5.0 / 8.0 * E_22


def composite_density(rho_m: float, V_f: float) -> float:
    """Rule-of-mixtures density [kg/m³]."""
    return rho_m * (1.0 - V_f) + RHO_FIBER * V_f


# ── Load inputs ───────────────────────────────────────────────────────────────
props = json.loads(pathlib.Path("inputs/properties.json").read_text())
build = json.loads(pathlib.Path("inputs/build_report.json").read_text())

fiber_wt_frac = build["fiber_loading_wt_pct"] / 100.0

# Use literature density for composite mechanics (MD density is under-converged
# at the default short equilibration time)
rho_m_lit = build["density_gcc"] * 1000.0   # g/cc → kg/m³

V_f = weight_to_volume_fraction(fiber_wt_frac, rho_m_lit)

# ── Apply correction ──────────────────────────────────────────────────────────
E_m = props["youngs_modulus_gpa"]

if fiber_wt_frac <= 0.0 or E_m is None:
    E_corrected     = E_m
    rho_corrected   = props["avg_density_kg_m3"]
    gf_applied      = False
    note            = "fiber_loading=0 — neat resin properties"
else:
    # Halpin-Tsai needs a positive, physical E_m.
    # If MD E is unphysical (negative / very small), fall back to literature
    # rule-of-thumb for amorphous polymers: E ≈ 0.003 * density_MPa
    if E_m is None or E_m <= 0:
        E_m_eff = rho_m_lit * 0.003 / 1e3   # crude fallback in GPa
        note = (f"MD E unphysical ({E_m} GPa); fell back to estimate "
                f"{E_m_eff:.3f} GPa for Halpin-Tsai input")
    else:
        E_m_eff = E_m
        note = "Halpin-Tsai + Christensen-Waals random short-fiber"

    E_corrected   = halpin_tsai_random(E_m_eff, V_f)
    rho_corrected = composite_density(rho_m_lit, V_f)
    gf_applied    = True

spec_mod = (E_corrected * 1e9 / rho_corrected / 1e3
            if (E_corrected and rho_corrected) else None)

# ── Write output ──────────────────────────────────────────────────────────────
out = {**props}
out.update({
    "gf_correction_applied":        gf_applied,
    "fiber_volume_fraction":        round(V_f, 4),
    "E_fiber_gpa":                  E_FIBER_GPA,
    "aspect_ratio":                 ASPECT_RATIO,
    "corrected_youngs_modulus_gpa": round(E_corrected, 3) if E_corrected else None,
    "corrected_density_kg_m3":      round(rho_corrected, 1),
    "corrected_specific_modulus_kNm_kg": round(spec_mod, 1) if spec_mod else None,
    "note":                         note,
})

pathlib.Path("outputs/corrected_properties.json").write_text(
    json.dumps(out, indent=2)
)

print(f"  Resin/loading : {build['resin_type']} + {build['fiber_loading_wt_pct']} wt% GF")
print(f"  V_f           : {V_f:.3f}")
print(f"  E (neat MD)   : {props['youngs_modulus_gpa']} GPa")
print(f"  E (corrected) : {E_corrected:.3f} GPa" if E_corrected else "  E (corrected) : N/A")
print(f"  Density (corr): {rho_corrected:.1f} kg/m³")
print(f"  Spec. modulus : {spec_mod:.1f} kN·m/kg" if spec_mod else "  Spec. modulus : N/A")
print(f"  Note          : {note}")
