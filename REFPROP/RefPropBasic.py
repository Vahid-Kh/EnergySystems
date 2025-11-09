"""
Compute fluid densities using CoolProp with both the native CoolProp backend
and the REFPROP backend (if installed), including pure fluids and mixtures.

Requirements:
- CoolProp installed (e.g., `pip install CoolProp`)
- REFPROP installed (optional) and accessible via RPPREFIX environment variable
  (on Windows, typically: C:\Program Files (x86)\REFPROP or C:\Program Files\REFPROP)

Notes:
- Pressure is in SI units (Pa), temperature in K, density in kg/m^3.
- For REFPROP backend usage with CoolProp, fluid names must be prefixed with "REFPROP::".
- Mixture composition is expressed in molar fractions.
- If you copied this script from a rendered HTML page, make sure to replace HTML entities:
  - Use ">" instead of "&gt;"
  - Use "&" instead of "&amp;"
"""

import os
import CoolProp.CoolProp as CP

# ---------------------------------------------------------------------------
# 1) Configure REFPROP integration (optional)
#    Set the REFPROP installation path via the RPPREFIX environment variable.
#    Update the path below to match your local installation if different.
# ---------------------------------------------------------------------------
os.environ['RPPREFIX'] = r'C:\Program Files (x86)\REFPROP'

# ---------------------------------------------------------------------------
# 2) Define state point (SI units)
#    Pressure (P): 6.9 MPa  ->  6.9e6 Pa
#    Temperature (T): 29 °C  ->  29 + 273.15 K
# ---------------------------------------------------------------------------
P = 69e5                 # Pa (6.9e6 Pa)
T = 29 + 273.15          # K

# ---------------------------------------------------------------------------
# 3) Example A: Pure CO2 density from CoolProp (native backend)
#    "CO2" without REFPROP prefix uses CoolProp's built-in EOS.
# ---------------------------------------------------------------------------
fluid = 'CO2'
rho_co2_cp = CP.PropsSI("D", "P", P, "T", T, fluid)
print("From COOLPROP  -->", rho_co2_cp, "kg/m^3")

# ---------------------------------------------------------------------------
# 4) Example B: Pure CO2 density from REFPROP
#    Prefix "REFPROP::" to route the query through the REFPROP backend.
#    Requires REFPROP to be properly installed and RPPREFIX set.
# ---------------------------------------------------------------------------
fluid = "REFPROP::CO2"
rho_co2_rp = CP.PropsSI("D", "P", P, "T", T, fluid)
print("From REFPROP   -->", rho_co2_rp, "kg/m^3")

# ---------------------------------------------------------------------------
# 5) Example C: Pure R1234ze(E) from REFPROP
#    Pay attention to the exact REFPROP naming and case.
# ---------------------------------------------------------------------------
fluid_refprop = "REFPROP::R1234ZE(E)"
rho_refprop = CP.PropsSI("D", "P", P, "T", T, fluid_refprop)  # kg/m^3
print("R1234ze(E) from REFPROP -->", rho_refprop, "kg/m^3")


# Given state: T = 45 °C, quality Q = 0.5 (two-phase)
T_C = 45.0
T_K = T_C + 273.15  # convert to Kelvin
Q = 0.7

# Speed of sound (A) in m/s
a = CP.PropsSI("A", "T", T_K, "Q", Q, fluid_refprop)

print(f"R1234ze(E) from REFPROP -> speed of sound at T={T_C} °C, Q={Q}: {a:.6g} m/s")


# ---------------------------------------------------------------------------
# 6) Example D: Native CoolProp backend for R1234ze(E) (for comparison)
# ---------------------------------------------------------------------------
fluid_coolprop = "R1234ze(E)"
rho_coolprop = CP.PropsSI("D", "P", P, "T", T, fluid_coolprop)
print("R1234ze(E) from COOLPROP -->", rho_coolprop, "kg/m^3")

# ---------------------------------------------------------------------------
# 7) Mixture examples (molar fractions)
#    IMPORTANT: Use '&' (ampersand) to separate components in CoolProp mixture strings.
#               Do NOT use HTML-escaped '&amp;'.
#
#    Example A: 60% R1234ze(E) + 40% R32
#    Example B: 50% R32 + 50% R125
# ---------------------------------------------------------------------------
mix = "REFPROP::R1234ZE(E)[0.6]&R32[0.4]"
mix2 = "REFPROP::R32[0.5]&R125[0.5]"

rho_mix = CP.PropsSI("D", "P", P, "T", T, mix)
rho_mix2 = CP.PropsSI("D", "P", P, "T", T, mix2)

print("Mixture (60% R1234ze(E) + 40% R32) from REFPROP -->", rho_mix, "kg/m^3")
print("Mixture (50% R32 + 50% R125) from REFPROP -->", rho_mix2, "kg/m^3")


"""
Density calculations via CoolProp/REFPROP with reusable helpers.

Run:
    python densities.py

Customize:
- Adjust RPPREFIX to your REFPROP installation path.
- Edit the STATE dict or the lists in main() for different fluids/mixtures.

Dependencies:
- CoolProp (pip install CoolProp)
- REFPROP (optional) with valid license/installation.
"""

import os
import sys
import CoolProp.CoolProp as CP

# Configure REFPROP path (update as needed)
os.environ.setdefault('RPPREFIX', r'C:\Program Files (x86)\REFPROP')

# State point (SI units)
STATE = {
    "P": 69e5,            # Pa
    "T": 29.0 + 273.15,   # K
}

def density(fluid: str, P: float, T: float) -> float:
    """
    Return density [kg/m^3] for a given fluid string, pressure [Pa], temperature [K].

    Fluid string formats:
    - Native CoolProp: "CO2", "R1234ze(E)"
    - REFPROP: "REFPROP::CO2", "REFPROP::R1234ZE(E)",
               mixtures like "REFPROP::R1234ZE(E)[0.6]&R32[0.4]"

    Raises:
        ValueError: If the property call fails (e.g., bad fluid name, out of range).
    """
    try:
        return CP.PropsSI("D", "P", P, "T", T, fluid)
    except Exception as exc:
        raise ValueError(f"Failed to compute density for fluid '{fluid}' at P={P} Pa, T={T} K: {exc}")

def print_density(label: str, fluid: str, P: float, T: float):
    """Helper to compute and print density with a label."""
    try:
        rho = density(fluid, P, T)
        print(f"{label:<45s} --> {rho:.6f} kg/m^3")
    except ValueError as e:
        print(f"{label:<45s} --> ERROR: {e}")

def main():
    P = STATE["P"]
    T = STATE["T"]

    print(f"Using REFPROP path: {os.environ.get('RPPREFIX')}")
    print(f"State point: P={P:.3e} Pa, T={T:.2f} K\n")

    # Pure fluids
    print_density("CO2 (CoolProp)", "CO2", P, T)
    print_density("CO2 (REFPROP)", "REFPROP::CO2", P, T)
    print_density("R1234ze(E) (CoolProp)", "R1234ze(E)", P, T)

    print_density("R1234ze(E) (REFPROP)", "REFPROP::R1234ZE(E)", P, T)

    # Mixtures (molar fractions)
    print()  # blank line
    mix_a = "REFPROP::R1234ZE(E)[0.6]&R32[0.4]"
    mix_b = "REFPROP::R32[0.5]&R125[0.5]"
    print_density("Mixture: 60% R1234ze(E) + 40% R32 (REFPROP)", mix_a, P, T)
    print_density("Mixture: 50% R32 + 50% R125 (REFPROP)", mix_b, P, T)

if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        sys.exit(130)

