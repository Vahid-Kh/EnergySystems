#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Pump power for 1 MW cooling across a 4,000 m pipe:
- Includes Darcy–Weisbach (Churchill friction), material roughness, viscosity, density.
- Mass flow is computed for each medium to deliver 1 MW.
- Line diameter is sized to respect a max liquid velocity (default 2 m/s).
"""

import math
from dataclasses import dataclass
from typing import Optional, Dict, Tuple, List
from tabulate import tabulate

# ----------------------------- CONFIG ---------------------------------

LENGTH_M = 4000.0            # fixed total equivalent length (straight + fittings if you like)
PUMP_EFF = 0.70              # pump efficiency (shaft/hydraulic)
MINOR_LOSS_K = 0.0           # lumped minor-loss coefficient (add if elbows/valves, etc.)
T_REF_C = 29.0               # reference temperature for properties
DT_SENSIBLE = 10.0           # ΔT for single-phase fluids (K)
V_MAX_LIQ = 2.0              # max liquid-line velocity (m/s) used for diameter sizing

# Choose pipe material (sets roughness); default SS304
PIPE_MATERIAL = "SS304"

# Optional: per-fluid override of velocity limits (leave empty to use V_MAX_LIQ for all)
PER_FLUID_VMAX: Dict[str, float] = {
    # "CO2": 2.5,    # example of override
}

# Roughness map [m] for common materials (typical clean, drawn tubing)
ROUGHNESS_MAP = {
    "SS304": 1.5e-5,      # stainless steel (drawn/annealed ~1–2e-5 m)
    "SS316": 1.5e-5,
    "CarbonSteel": 4.5e-5,# commercial steel
    "Copper": 1.5e-5,
    "Aluminum": 1.5e-5,
    "Titanium": 1.5e-5,
    "PVC": 1.5e-6,
}

# Refrigerants to evaluate
# phase: "1Ø" (single-phase liquid) or "2Ø" (two-phase refrigeration; pump on liquid line)
REFRIGERANTS = [
    {"name": "Water",            "phase": "1Ø"},
    {"name": "INCOMP::MEG-25%",  "phase": "1Ø"},
    {"name": "CO2",              "phase": "2Ø"},
    {"name": "Ammonia",          "phase": "2Ø"},  # R717
    {"name": "R290",             "phase": "2Ø"},  # Propane
    {"name": "R1234ze(E)",       "phase": "2Ø"},
    {"name": "R1233zd(E)",       "phase": "2Ø"},
    {"name": "R134a",            "phase": "2Ø"},
    {"name": "R32",              "phase": "2Ø"},
]

# ------------------------ CoolProp (optional) --------------------------

try:
    from CoolProp.CoolProp import PropsSI
    HAVE_COOLPROP = True
except Exception:
    HAVE_COOLPROP = False
    PropsSI = None

# -------------------------- Utilities ---------------------------------

def to_K(t_c: float) -> float:
    return t_c + 273.15

def friction_factor_churchill(Re: float, rel_eps: float) -> float:
    """
    Churchill correlation for Darcy friction factor f (smooth-laminar to fully rough-turbulent).
    """
    Re = max(Re, 1e-12)
    if Re < 2300.0:
        return max(64.0/Re, 1e-8)
    # Churchill 1977
    term = (7.0/Re)**0.9 + 0.27*rel_eps
    term = max(term, 1e-12)
    A = (2.457*math.log(1.0/term))**16.0
    B = (37530.0/Re)**16.0
    f = 8.0 * (( (8.0/Re)**12.0 ) + 1.0/((A + B)**1.5) )**(1.0/12.0)
    return max(f, 1e-8)

def darcy_delta_p(rho: float, mu: float, L: float, D: float, v: float,
                  roughness_m: float, K_minor: float = 0.0) -> Tuple[float, float, float]:
    """
    Compute Δp [Pa] with Darcy–Weisbach (major + minor), returning (Δp_total, f, Re)
    """
    Re = rho * v * D / max(mu, 1e-12)
    f = friction_factor_churchill(Re, roughness_m / max(D, 1e-12))
    dyn_head = 0.5 * rho * v * v
    dp_major = f * (L / max(D, 1e-12)) * dyn_head
    dp_minor = K_minor * dyn_head
    return dp_major + dp_minor, f, Re

def pick_vmax(name: str) -> float:
    return PER_FLUID_VMAX.get(name, V_MAX_LIQ)

# ---------------------- Property helpers ------------------------------

@dataclass
class FluidProps:
    rho: float         # kg/m3
    mu: float          # Pa·s
    cp: Optional[float] = None   # J/kg/K (single-phase)
    h_fg: Optional[float] = None # J/kg (two-phase latent at T_ref)
    note: str = ""

def get_single_phase_props(name: str, T_K: float) -> FluidProps:
    """
    Single-phase liquid properties at (T_K, P=1 bar) when possible.
    Fallbacks used if CoolProp not available.
    """
    if HAVE_COOLPROP:
        rho = PropsSI("D", "T", T_K, "P", 1e5, name)
        mu  = PropsSI("V", "T", T_K, "P", 1e5, name)
        cp  = PropsSI("C", "T", T_K + DT_SENSIBLE/2.0, "P", 1e5, name)
        return FluidProps(rho=rho, mu=mu, cp=cp)
    # Fallbacks (approximate at ~25–30°C)
    if name.lower().startswith("incomp::meg"):
        return FluidProps(rho=1030.0, mu=0.0030, cp=3300.0, note="fallback MEG-25% @~25°C")
    if name.lower() == "water":
        return FluidProps(rho=997.0, mu=0.00089, cp=4180.0, note="fallback water @25°C")
    # Generic fallback
    return FluidProps(rho=1000.0, mu=0.0015, cp=3500.0, note="generic fallback liquid")

def get_two_phase_liquid_line_props(name: str, T_K: float) -> FluidProps:
    """
    For two-phase refrigerants, we size/pump the LIQUID line. Use saturated liquid ρ, μ at T_K.
    Also compute latent heat h_fg for mass flow.
    """
    if HAVE_COOLPROP:
        rho_L = PropsSI("D", "T", T_K, "Q", 0.0, name)
        mu_L  = PropsSI("V", "T", T_K, "Q", 0.0, name)
        h_v   = PropsSI("H", "T", T_K, "Q", 1.0, name)
        h_l   = PropsSI("H", "T", T_K, "Q", 0.0, name)
        h_fg  = max(h_v - h_l, 1e-6)
        return FluidProps(rho=rho_L, mu=mu_L, h_fg=h_fg)
    # Fallback latent heats (very rough, order-of-magnitude at ~25–35°C)
    fallback_hfg = {
        "co2": 2.0e5,       # J/kg (varies strongly with T, near critical)
        "ammonia": 1.15e6,  # J/kg
        "r717": 1.15e6,
        "r290": 3.5e5,
        "r1234ze(e)": 1.5e5,
        "r1233zd(e)": 1.2e5,
        "r134a": 2.0e5,
        "r32":  2.9e5,
    }
    key = name.lower()
    hfg = fallback_hfg.get(key, 2.0e5)
    # Very rough liquid properties (only to let the math run if CoolProp missing)
    rho_L = 1100.0 if key in ("ammonia", "r717") else 1000.0
    mu_L  = 0.00025
    return FluidProps(rho=rho_L, mu=mu_L, h_fg=hfg, note="fallback two-phase")

# ------------------------- Core calculation ---------------------------

@dataclass
class ResultRow:
    name: str
    phase: str
    material: str
    roughness_m: float
    rho: float
    mu: float
    mdot_kg_s: float
    v_m_s: float
    D_mm: float
    Re: float
    f: float
    dp_bar: float
    P_hyd_kW: float
    P_shaft_kW: float
    note: str

def size_and_compute(name: str, phase: str,
                     T_C: float,
                     L_m: float,
                     material: str) -> Optional[ResultRow]:
    T_K = to_K(T_C)
    eps = ROUGHNESS_MAP.get(material, 4.5e-5)

    try:
        if phase == "1Ø":
            fp = get_single_phase_props(name, T_K)
            # Mass flow for 1 MW (sensible)
            mdot = 1.0e6 / (max(fp.cp, 1e-9) * DT_SENSIBLE)
        elif phase == "2Ø":
            fp = get_two_phase_liquid_line_props(name, T_K)
            # Mass flow for 1 MW (latent)
            mdot = 1.0e6 / max(fp.h_fg, 1e-9)
        else:
            raise ValueError(f"Unknown phase '{phase}'")

        vmax = pick_vmax(name)
        # Size diameter from mdot = rho * v * A  =>  A = mdot/(rho*v)
        A = mdot / (max(fp.rho, 1e-12) * max(vmax, 1e-9))
        D = math.sqrt(4.0 * A / math.pi)
        v = mdot / (fp.rho * (math.pi * D * D / 4.0))  # should be ~vmax

        # Pressure loss & power
        dp_Pa, f, Re = darcy_delta_p(fp.rho, fp.mu, L_m, D, v, eps, K_minor=MINOR_LOSS_K)
        Q_vol = mdot / max(fp.rho, 1e-12)
        P_hyd = dp_Pa * Q_vol
        P_shaft = P_hyd / max(min(PUMP_EFF, 0.9999), 1e-6)

        return ResultRow(
            name=name, phase=phase, material=material, roughness_m=eps,
            rho=fp.rho, mu=fp.mu, mdot_kg_s=mdot, v_m_s=v,
            D_mm=D*1000.0, Re=Re, f=f, dp_bar=dp_Pa/1e5,
            P_hyd_kW=P_hyd/1000.0, P_shaft_kW=P_shaft/1000.0, note=fp.note
        )
    except Exception as e:
        print(f"[WARN] Skipping {name} ({phase}): {e}")
        return None

# ------------------------------- Run ----------------------------------

def main():
    rows: List[ResultRow] = []
    for r in REFRIGERANTS:
        res = size_and_compute(
            name=r["name"],
            phase=r["phase"],
            T_C=T_REF_C,
            L_m=LENGTH_M,
            material=PIPE_MATERIAL
        )
        if res:
            rows.append(res)

    headers = [
        "Medium", "Phase", "Material(ε m)", "ρ [kg/m³]", "μ [Pa·s]",
        "ṁ [kg/s]", "v [m/s]", "D [mm]", "Re", "f", "Δp [bar]",
        "P_hyd [kW]", "P_shaft [kW]", "Note"
    ]

    table = []
    for r in rows:
        table.append([
            r.name, r.phase, f"{r.material} ({r.roughness_m:.2e})",
            f"{r.rho:.1f}", f"{r.mu:.3e}",
            f"{r.mdot_kg_s:.2f}", f"{r.v_m_s:.3f}", f"{r.D_mm:.1f}",
            f"{r.Re:.0f}", f"{r.f:.4f}", f"{r.dp_bar:.3f}",
            f"{r.P_hyd_kW:.2f}", f"{r.P_shaft_kW:.2f}", r.note
        ])

    print("\nAssumptions:")
    print(f"- Length = {LENGTH_M:.0f} m, Pump η = {PUMP_EFF:.2f}, Minor K = {MINOR_LOSS_K:.2f}")
    print(f"- T_ref = {T_REF_C:.1f} °C, ΔT (single-phase) = {DT_SENSIBLE:.1f} K")
    print(f"- Diameter sized for v_max (liquid) = {V_MAX_LIQ:.2f} m/s (overrides per fluid allowed)")
    print(f"- Pipe material roughness = {PIPE_MATERIAL} → ε = {ROUGHNESS_MAP.get(PIPE_MATERIAL, 4.5e-5):.2e} m")
    print("- Two‑phase entries assume mass flow set by latent heat at T_ref and pumping the liquid line.\n")

    print(tabulate(table, headers=headers, tablefmt="grid", floatfmt=".3f"))

if __name__ == "__main__":
    main()
