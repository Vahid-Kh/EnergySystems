
"""MINIMALIST VERSION"""
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

# Path to REFPROP (adjust as needed)
rp_root = r"C:\Program Files (x86)\REFPROP"
RP = REFPROPFunctionLibrary(rp_root)
RP.SETPATHdll(rp_root)

# Load fluid (R134a)
RP.SETUPdll(1, "R134A.FLD", "HMX.BNC", "DEF")

# State point (T in K, P in kPa), pure fluid
T = 29.66 + 273.15
P = 9 * 101.325
z = [1.0]

# Flash calculation
D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr, herr = RP.TPFLSHdll(T, P, z)
if ierr != 0:
    raise RuntimeError(herr)

print(f"Density (mol/L): {D:.6f}")
print(f"Enthalpy (J/mol): {h:.2f}")
print(f"Cp (J/mol/K): {Cp:.2f}")
print(f"Speed of sound (m/s): {w:.2f}")

"""Full Version"""
import os
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

# 1) Point to your REFPROP installation folder
#    (Adjust the path to where you installed REFPROP)
rp_root = r"C:\Program Files (x86)\REFPROP"  # Windows example
RP = REFPROPFunctionLibrary(rp_root)
RP.SETPATHdll(rp_root)  # ensure the DLL finds /Fluids and /Mixtures

# 2) Load fluid(s) – here pure water
#    SETUPdll(nc, "FLUID.FLD[|...]", "HMX.BNC", "DEF")

"""Default units"""
RP.SETUPdll(1, "R134A.FLD", "HMX.BNC", "DEF")

# 3) Flash calculation at T (K) and P (kPa); composition z is mole fractions
T = 29.66+ 273.15         # K
P = 9*101.325        # kPa
z = [1.0]          # 100% composition

# TPFLSHdll returns: D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr, herr (default units noted below)
out = RP.TPFLSHdll(T, P, z)
D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr, herr = out

if ierr != 0:
    raise RuntimeError(herr)

print(f"Density (mol/L): {D:.6f}")
print(f"Enthalpy (J/mol): {h:.2f}")
print(f"Cp (J/mol/K): {Cp:.2f}")
print(f"Speed of sound (m/s): {w:.2f}")


""" Using P& T and P& Q"""
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary

# --- Setup ---
RP = REFPROPFunctionLibrary(r"C:\Program Files (x86)\REFPROP")
RP.SETPATHdll(r"C:\Program Files (x86)\REFPROP")
RP.SETUPdll(1, "R1234ZEE.FLD", "HMX.BNC", "DEF")  # pure R1234ZEE

def _wm(z=[1.0]):
    """Return molar mass in kg/kmol; robust to different return shapes."""
    res = RP.WMOLdll(z)
    return res[0] if isinstance(res, (list, tuple)) else float(res)

def h_PT(T_K, P_kPa, mass_basis=False):
    """Enthalpy from T [K], P [kPa]. Returns J/mol (default) or J/kg."""
    D, Dl, Dv, x, y, q, e, h, s, Cv, Cp, w, ierr, herr = RP.TPFLSHdll(T_K, P_kPa, [1.0])
    if ierr: raise RuntimeError(herr)
    return h if not mass_basis else (h * 1000.0 / _wm())

def h_Pq(P_kPa, q, mass_basis=False, kq=2):
    """Enthalpy from P [kPa], quality q in [0..1]; kq=2 => vapor quality (molar)."""
    T, D, Dl, Dv, x, y, e, h, s, Cv, Cp, w, ierr, herr = RP.PQFLSHdll(P_kPa, q, [1.0], kq)
    if ierr: raise RuntimeError(herr)
    return h if not mass_basis else (h * 1000.0 / _wm())

# --- Quick checks ---
print("h(T,P) J/mol:", h_PT(300.0, 300.0))
print("h(P,q=1) J/kg:", h_Pq(300.0, 1.0, mass_basis=True))


"""Coolprop Interface"""

from CoolProp.CoolProp import PropsSI

# --- Fluid selection via REFPROP backend ---
# Using REFPROP backend explicitly; CoolProp will call REFPROP if available.
FLUID = "REFPROP::R1234ZEE"  # pure R1234ZEE

def _wm():
    """
    Return molar mass in kg/mol (CoolProp/REFPROP basis).
    """
    M = PropsSI("M", "T", 300.0, "P", 101325.0, FLUID)  # any valid state is fine
    return float(M)

def h_PT(T_K, P_kPa, mass_basis=False):
    """
    Enthalpy from T [K], P [kPa].
    Returns J/mol (default) or J/kg if mass_basis=True.
    """
    P_Pa = P_kPa * 1000.0
    # CoolProp returns mass-specific enthalpy in J/kg
    h_mass = PropsSI("H", "T", T_K, "P", P_Pa, FLUID)
    if mass_basis:
        return h_mass
    else:
        # convert J/kg -> J/mol using M [kg/mol]
        return h_mass * _wm()

def h_Pq(P_kPa, q, mass_basis=False, kq=2):
    """
    Enthalpy from P [kPa], quality q in [0..1].
    In CoolProp, 'Q' is quality (mass-based vapor quality).
    'kq' is kept for compatibility but is ignored (must be vapor quality).
    Returns J/mol (default) or J/kg if mass_basis=True.
    """
    P_Pa = P_kPa * 1000.0
    # Mass-based enthalpy at given pressure and quality
    h_mass = PropsSI("H", "P", P_Pa, "Q", q, FLUID)
    if mass_basis:
        return h_mass
    else:
        return h_mass * _wm()

# --- Quick checks (similar to your originals) ---
if __name__ == "__main__":
    print("h(T,P) J/mol:", h_PT(300.0, 300.0))
    print("h(P,q=1) J/kg:", h_Pq(300.0, 1.0, mass_basis=True))
