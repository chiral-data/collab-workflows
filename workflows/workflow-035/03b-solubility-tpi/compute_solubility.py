"""
Compute the Henry-regime solubility coefficient S from the TPI excess
chemical potential mu_ex (gmx tpi reports "<mu> = ... kJ/mol" in its log).

Derivation (ideal-gas / Henry's-law equilibrium between the polymer and gas
phases, standard for Widom-insertion solubility estimates):
    mu_gas(p)  = kT ln(p / kT * Lambda^3)
    mu_polymer = kT ln(rho_sol * Lambda^3) + mu_ex
    mu_gas = mu_polymer  =>  rho_sol = (p / kT) * exp(-mu_ex / kT)
    S = rho_sol / p = 1/(kT) * exp(-mu_ex / kT)
so S has units of (number density)/(pressure); per mole, S = 1/(RT)*exp(-mu_ex/RT).
"""
import argparse, math, json, pathlib, re, sys

R = 8.314e-3  # kJ/(mol K)

parser = argparse.ArgumentParser()
parser.add_argument("--log",       required=True, help="tpi.log from gmx mdrun")
parser.add_argument("--report",    required=True)
parser.add_argument("--penetrant", default="O2")
parser.add_argument("--resin",     default="?")
parser.add_argument("--temp",      type=float, default=23.0, help="Temperature in °C")
parser.add_argument("--n_frames",  type=int, default=0, help="Trajectory frames reran over")
parser.add_argument("--n_insertions_per_frame", type=int, default=0)
args = parser.parse_args()

log_text = pathlib.Path(args.log).read_text()
m = re.search(r"<mu>\s*=\s*([-\d.eE+]+)\s*kJ/mol", log_text)
if not m:
    sys.exit("Could not find '<mu> = ... kJ/mol' in tpi.log — TPI run likely failed or produced no insertions.")
mu_ex_kj_mol = float(m.group(1))

T_k = args.temp + 273.15
RT = R * T_k  # kJ/mol

# S in SI units: mol / (m^3 * Pa)
S_si = (1.0 / (RT * 1000.0)) * math.exp(-mu_ex_kj_mol / RT)

# S in the conventional membrane-science unit: cm^3(STP) gas / (cm^3 polymer * cmHg)
# 1 mol ideal gas at STP = 22414 cm^3; 1 m^3 = 1e6 cm^3; 1 Pa = 1/1333.22 cmHg
STP_CM3_PER_MOL = 22414.0
PA_PER_CMHG = 1333.22
S_cmHg = S_si * STP_CM3_PER_MOL * 1e-6 * PA_PER_CMHG

qualitative = args.penetrant == "H2O"

report = {
    "resin_type":        args.resin,
    "penetrant":          args.penetrant,
    "temperature_c":      args.temp,
    "mu_ex_kj_mol":       round(mu_ex_kj_mol, 4),
    "S_mol_m3_pa":        S_si,
    "S_cm3stp_cm3_cmhg":  S_cmHg,
    "n_frames":           args.n_frames,
    "n_insertions_per_frame": args.n_insertions_per_frame,
    "qualitative":        qualitative,
    "water_model":        "SPC/E" if args.penetrant == "H2O" else None,
    "note": (
        "H2O solubility is Henry-regime (infinite dilution) and reported as a "
        "qualitative trend only — real water sorption in polar resins (EVOH, PA6) "
        "deviates from Henry's law via clustering/swelling."
    ) if qualitative else None,
}
pathlib.Path(args.report).write_text(json.dumps(report, indent=2))

print(f"  mu_ex        : {mu_ex_kj_mol:.4f} kJ/mol")
print(f"  S (SI)        : {S_si:.4e} mol/(m3*Pa)")
print(f"  S (cmHg conv.): {S_cmHg:.4e} cm3(STP)/(cm3*cmHg)")
if qualitative:
    print("  NOTE: H2O — qualitative trend only (see report note)")
