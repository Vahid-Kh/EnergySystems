import os
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

# --- 1) REFPROP install path ---
rp_root = r"C:\Program Files (x86)\REFPROP"  # adjust if needed (match Python bitness)
RP = REFPROPFunctionLibrary(rp_root)
RP.SETPATHdll(rp_root)

# --- 2) Load fluid ---
resp = RP.SETUPdll(1, "R134A.FLD", "HMX.BNC", "DEF")
if resp.ierr != 0:
    raise RuntimeError(resp.herr)

# --- 3) Select unit system = Mass-based SI ---
try:
    units_enum = RP.GETENUMdll(1, "Mass Base SI")  # iFlag=1 → unit systems
    iUnits = units_enum.iEnum
except Exception:
    iUnits = 21  # fallback; typical code for Mass-based SI

# --- 4) State point ---
T = 29 + 273.15   # K
P = 9 * 101.325      # kPa
z = [1.0]

# --- 5) Properties request ---
req_props = "H;S;D;CP;W"
tokens = req_props.split(";")

# Fallback unit labels for Mass-based SI (aligned with the tokens)
fallback_units_map = {
    "H":  "kJ/kg",
    "S":  "kJ/kg-K",
    "D":  "kg/m^3",
    "CP": "kJ/kg-K",
    "W":  "m/s",
    # If you add more tokens later, extend this dict (e.g., "V": "Pa·s", "L": "W/m-K")
}

r = RP.REFPROPdll("R134A", "TP", req_props, iUnits, 0, 0, T, P, z)
if r.ierr != 0:
    raise RuntimeError(r.herr)

# --- 6) Safe unpack of numeric outputs ---
# Some wrapper builds return more items; only take as many as we requested
values = list(r.Output[:len(tokens)])

# --- 7) Units: prefer r.hUnits if present & aligned; else build from map ---
units_str = getattr(r, "hUnits", "") or ""
units_list = [u.strip() for u in units_str.split(";") if u.strip()]
if len(units_list) != len(tokens):
    # Build units from our fallback map
    units_list = [fallback_units_map.get(tok, "") for tok in tokens]

# --- 8) Pretty print ---
print(f"T = {T:.2f} K, P = {P:.3f} kPa")
for tok, val, unit in zip(tokens, values, units_list):
    if tok == "H":
        label = "Enthalpy     "
    elif tok == "S":
        label = "Entropy      "
    elif tok == "D":
        label = "Density      "
    elif tok == "CP":
        label = "Cp           "
    elif tok == "W":
        label = "Speed of snd."
    else:
        label = f"{tok:<12}"
    # Use generic formatting for safety (float-like)
    try:
        print(f"{label}: {float(val):.6f} {unit}")
    except Exception:
        print(f"{label}: {val} {unit}")

# --- 9) (Optional) Debug helpers ---
# print("len(Output) :", len(r.Output))
# print("hUnits raw  :", getattr(r, "hUnits", None))
# print("units_list  :", units_list)

import CoolProp.CoolProp as CP

T = 29 + 273.15 # Convert 29 degC to Kelvin
P = 10 * 100000 # Convert 10 bar to Pa (or 1.0 MPa) - check your wrapper's unit convention

# Get density in kg/m³
density_kg_per_m3 = CP.PropsSI('D', 'T', T, 'P', P, 'R134a')

# Convert to kg/L
density_kg_per_L = density_kg_per_m3 / 1000.0

print(f"Density in kg/L: {density_kg_per_L}")