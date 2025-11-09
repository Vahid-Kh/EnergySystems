import math
import CoolProp.CoolProp as CP

def _mix_rho_from_quality(x, rho_l, rho_v):
    """
    Two-phase mixture density using the lever rule on specific volumes:
        1/rho = (1-x)/rho_l + x/rho_v
    where x is mass quality (kg_vapor / kg_total).
    """
    return 1.0 / ((1.0 - x) / rho_l + x / rho_v)

def _sat_props_T(fluid, T):
    """
    Saturated properties at temperature T [K].
    Returns dict with p, rho_l, rho_v, s_l, s_v (SI units).
    """
    p = CP.PropsSI("P", "T", T, "Q", 0.5, fluid)  # saturation pressure
    rho_l = CP.PropsSI("D", "T", T, "Q", 0.0, fluid)
    rho_v = CP.PropsSI("D", "T", T, "Q", 1.0, fluid)
    s_l = CP.PropsSI("S", "T", T, "Q", 0.0, fluid)
    s_v = CP.PropsSI("S", "T", T, "Q", 1.0, fluid)
    return {"p": p, "rho_l": rho_l, "rho_v": rho_v, "s_l": s_l, "s_v": s_v}

def _sat_props_P(fluid, p):
    """
    Saturated properties at pressure p [Pa].
    Returns dict with T, rho_l, rho_v, s_l, s_v (SI units).
    """
    T = CP.PropsSI("T", "P", p, "Q", 0.5, fluid)  # saturation temperature
    rho_l = CP.PropsSI("D", "P", p, "Q", 0.0, fluid)
    rho_v = CP.PropsSI("D", "P", p, "Q", 1.0, fluid)
    s_l = CP.PropsSI("S", "P", p, "Q", 0.0, fluid)
    s_v = CP.PropsSI("S", "P", p, "Q", 1.0, fluid)
    return {"T": T, "rho_l": rho_l, "rho_v": rho_v, "s_l": s_l, "s_v": s_v}

def hem_speed_of_sound_TQ(fluid, T, Q, dp=None, max_halvings=10):
    """
    Homogeneous-Equilibrium-Model (HEM) speed of sound for a two-phase saturated state.
    State is specified by (T [K], Q mass quality in (0,1)).
    Returns speed of sound [m/s].

    Method: numerical isentropic derivative across saturation: a^2 = (∂p/∂ρ)_s.
    """
    if not (0.0 < Q < 1.0):
        raise ValueError("Quality Q must be strictly between 0 and 1 for two-phase HEM.")

    # Base saturated properties at (T,Q)
    sat0 = _sat_props_T(fluid, T)
    p0, rho_l0, rho_v0, s_l0, s_v0 = sat0["p"], sat0["rho_l"], sat0["rho_v"], sat0["s_l"], sat0["s_v"]

    # Mixture entropy at the nominal state (mass-weighted; equilibrium)
    s_mix0 = (1.0 - Q) * s_l0 + Q * s_v0

    # Choose a small pressure step (Pa). Scale with p0, but clamp to a practical range.
    if dp is None:
        dp = max(50.0, 1e-6 * p0)  # default ~1 ppm of p or at least 50 Pa

    # Helper to compute rho_mix at p*, keeping s_mix constant (HEM constraint)
    def rho_at_p(p_target, dp_try):
        # Adapt dp so we remain inside the dome (i.e., x in (0,1)).
        p_left = p_target
        for _ in range(max_halvings):
            # Saturation at the new pressure
            sat = _sat_props_P(fluid, p_left)
            s_l, s_v = sat["s_l"], sat["s_v"]

            denom = (s_v - s_l)
            if abs(denom) < 1e-9:
                # Near critical: entropy difference is tiny -> ill-conditioned
                raise RuntimeError("Entropy difference between phases is too small (near critical?).")

            x = (s_mix0 - s_l) / denom  # new equilibrium quality along s = const
            if 0.0 < x < 1.0:
                # Valid two-phase state
                rho_l, rho_v = sat["rho_l"], sat["rho_v"]
                return _mix_rho_from_quality(x, rho_l, rho_v), x, p_left
            # Else, halve the step away from the boundary
            p_left = p0 + (p_left - p0) * 0.5

        # If we reach here, could not stay inside the dome
        raise RuntimeError("Failed to keep the perturbed state inside the two-phase region. Try a smaller dp or move away from the dome edges.")

    # Compute central difference for (∂ρ/∂p)_s
    rho_plus, x_plus, p_used_plus = rho_at_p(p0 + dp, dp)
    rho_minus, x_minus, p_used_minus = rho_at_p(p0 - dp, dp)

    d_rho_dp = (rho_plus - rho_minus) / (p_used_plus - p_used_minus)  # central slope
    if d_rho_dp <= 0.0:
        raise RuntimeError(f"Non-positive dρ/dp encountered ({d_rho_dp:.3e}). Check state (near critical?) or reduce dp.")

    a2 = 1.0 / d_rho_dp
    if a2 <= 0.0 or not math.isfinite(a2):
        raise RuntimeError("Invalid a^2 computed. Try adjusting dp or avoiding near-critical states.")

    return math.sqrt(a2)

def hem_speed_of_sound_PQ(fluid, P, Q, dp=None):
    """
    Same HEM method but the base state is specified by (P [Pa], Q).
    Internally, it maps to the saturation temperature at (P,Q), then calls TQ version.
    """
    if not (0.0 < Q < 1.0):
        raise ValueError("Quality Q must be strictly between 0 and 1 for two-phase HEM.")
    T = CP.PropsSI("T", "P", P, "Q", 0.5, fluid)
    return hem_speed_of_sound_TQ(fluid, T, Q, dp=dp)

def frozen_speed_of_sound_TQ(fluid, T, Q, dp=None):
    """
    'Frozen' (no phase change during the wave) homogeneous model using Wood's formula:
        1 / (rho_mix * a^2) = α_l / (rho_l * a_l^2) + α_v / (rho_v * a_v^2)
    where α_i are volume fractions, and a_l, a_v are *single-phase* speeds of sound
    just off saturation (nudged by ±dp).

    This is provided for comparison; it generally yields a much *higher* speed
    than HEM because phase change is 'frozen' here.
    """
    if not (0.0 < Q < 1.0):
        raise ValueError("Quality Q must be strictly between 0 and 1 for two-phase (frozen) model.")

    sat0 = _sat_props_T(fluid, T)
    p0, rho_l, rho_v = sat0["p"], sat0["rho_l"], sat0["rho_v"]
    # volume fractions
    alpha_v = (Q / rho_v) / (Q / rho_v + (1.0 - Q) / rho_l)
    alpha_l = 1.0 - alpha_v

    if dp is None:
        dp = max(50.0, 1e-6 * p0)

    # Get single-phase speeds close to saturation
    a_l = CP.PropsSI("A", "T", T, "P", p0 + dp, fluid)  # subcooled liquid side
    a_v = CP.PropsSI("A", "T", T, "P", p0 - dp, fluid)  # superheated vapor side

    rho_mix = _mix_rho_from_quality(Q, rho_l, rho_v)

    inv_rho_a2 = alpha_l / (rho_l * a_l**2) + alpha_v / (rho_v * a_v**2)
    a2 = 1.0 / (rho_mix * inv_rho_a2)
    return math.sqrt(a2)



fluid = "REFPROP::R1234ZE(E)"
T_C = 45.0
T = T_C + 273.15
Q = 0.9  # two-phase mid-quality

a_hem = hem_speed_of_sound_TQ(fluid, T, Q)          # HEM (equilibrium) speed of sound
a_frozen = frozen_speed_of_sound_TQ(fluid, T, Q)    # Wood's (frozen) for comparison

print(f"HEM speed of sound at T={T_C} °C, Q={Q}:    {a_hem:.4g} m/s")
print(f"Frozen speed of sound at T={T_C} °C, Q={Q}: {a_frozen:.4g} m/s")



import math
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

# -------------------------------
# HEM and helper implementations
# -------------------------------

def _mix_rho_from_quality(x, rho_l, rho_v):
    """Two-phase mixture density using lever rule on specific volumes."""
    return 1.0 / ((1.0 - x) / rho_l + x / rho_v)

def _sat_props_T(fluid, T):
    """Saturation properties at temperature T [K]."""
    p = CP.PropsSI("P", "T", T, "Q", 0.5, fluid)  # saturation pressure
    rho_l = CP.PropsSI("D", "T", T, "Q", 0.0, fluid)
    rho_v = CP.PropsSI("D", "T", T, "Q", 1.0, fluid)
    s_l = CP.PropsSI("S", "T", T, "Q", 0.0, fluid)
    s_v = CP.PropsSI("S", "T", T, "Q", 1.0, fluid)
    return {"p": p, "rho_l": rho_l, "rho_v": rho_v, "s_l": s_l, "s_v": s_v}

def _sat_props_P(fluid, p):
    """Saturation properties at pressure p [Pa]."""
    T = CP.PropsSI("T", "P", p, "Q", 0.5, fluid)  # saturation temperature
    rho_l = CP.PropsSI("D", "P", p, "Q", 0.0, fluid)
    rho_v = CP.PropsSI("D", "P", p, "Q", 1.0, fluid)
    s_l = CP.PropsSI("S", "P", p, "Q", 0.0, fluid)
    s_v = CP.PropsSI("S", "P", p, "Q", 1.0, fluid)
    return {"T": T, "rho_l": rho_l, "rho_v": rho_v, "s_l": s_l, "s_v": s_v}

def hem_speed_of_sound_TQ(fluid, T, Q, dp=None, max_halvings=10):
    """
    Homogeneous-Equilibrium-Model (HEM) speed of sound [m/s] at saturated state (T,Q).
    Uses numerical isentropic derivative: a^2 = (∂p/∂ρ)_s via central difference in pressure.
    """
    if not (0.0 < Q < 1.0):
        raise ValueError("Quality Q must be strictly between 0 and 1 for two-phase HEM.")

    sat0 = _sat_props_T(fluid, T)
    p0, rho_l0, rho_v0, s_l0, s_v0 = sat0["p"], sat0["rho_l"], sat0["rho_v"], sat0["s_l"], sat0["s_v"]

    # Mixture entropy at base state
    s_mix0 = (1.0 - Q) * s_l0 + Q * s_v0
    if dp is None:
        dp = max(50.0, 1e-6 * p0)  # ~1 ppm of p or at least 50 Pa

    def rho_at_p(p_target):
        # Try to keep the perturbed state within the two-phase region by step halving
        p_try = p_target
        for _ in range(max_halvings):
            sat = _sat_props_P(fluid, p_try)
            s_l, s_v = sat["s_l"], sat["s_v"]
            denom = (s_v - s_l)
            if abs(denom) < 1e-9:
                raise RuntimeError("Entropy difference between phases too small (near critical?).")

            x = (s_mix0 - s_l) / denom
            if 0.0 < x < 1.0:
                rho_l, rho_v = sat["rho_l"], sat["rho_v"]
                rho_mix = _mix_rho_from_quality(x, rho_l, rho_v)
                return rho_mix, x, p_try
            # Halve the distance towards p0 to stay inside dome
            p_try = p0 + 0.5 * (p_try - p0)
        raise RuntimeError("Failed to keep perturbed state inside two-phase region. Reduce dp or move away from dome edges.")

    rho_plus, x_plus, p_plus = rho_at_p(p0 + dp)
    rho_minus, x_minus, p_minus = rho_at_p(p0 - dp)

    d_rho_dp = (rho_plus - rho_minus) / (p_plus - p_minus)
    if d_rho_dp <= 0.0 or not math.isfinite(d_rho_dp):
        raise RuntimeError(f"Invalid dρ/dp={d_rho_dp:.3e}. Try smaller dp or avoid near-critical states.")

    a2 = 1.0 / d_rho_dp
    if a2 <= 0.0 or not math.isfinite(a2):
        raise RuntimeError("Invalid a^2 computed.")
    return math.sqrt(a2)

def frozen_speed_of_sound_TQ(fluid, T, Q, dp=None):
    """
    'Frozen' model using Wood's formula (no phase change during the wave).
    Useful as an upper-bound comparison to HEM.
    """
    if not (0.0 < Q < 1.0):
        raise ValueError("Quality Q must be strictly between 0 and 1 for two-phase (frozen) model.")
    sat0 = _sat_props_T(fluid, T)
    p0, rho_l, rho_v = sat0["p"], sat0["rho_l"], sat0["rho_v"]

    # Volume fractions
    alpha_v = (Q / rho_v) / (Q / rho_v + (1.0 - Q) / rho_l)
    alpha_l = 1.0 - alpha_v

    if dp is None:
        dp = max(50.0, 1e-6 * p0)

    # Single-phase speeds just off saturation
    a_l = CP.PropsSI("A", "T", T, "P", p0 + dp, fluid)  # liquid side
    a_v = CP.PropsSI("A", "T", T, "P", p0 - dp, fluid)  # vapor side

    rho_mix = _mix_rho_from_quality(Q, rho_l, rho_v)
    inv_rho_a2 = alpha_l / (rho_l * a_l**2) + alpha_v / (rho_v * a_v**2)
    a2 = 1.0 / (rho_mix * inv_rho_a2)
    return math.sqrt(a2)

# -------------------------------
# User Inputs
# -------------------------------
fluid = "REFPROP::R1234ZE(E)"  # Use "HEOS::R1234ZE(E)" if you don't have REFPROP
T_evap_C = 45.0                # Evaporating saturation temperature [°C] (adjust as needed)
subcool_K = 5.0                # extent of subcooling below Tsat [K]
superheat_K = 10.0             # extent of superheat above Tsat [K]

# -------------------------------
# Setup baseline saturation state
# -------------------------------
T_sat = T_evap_C + 273.15
p_evap = CP.PropsSI("P", "T", T_sat, "Q", 0.5, fluid)  # evaporating pressure at T_evap
dp_side = max(50.0, 1e-6 * p_evap)                     # tiny offset to ensure single-phase calls

# -------------------------------
# Build the composite x-axis
#   x in [-subcool_K, 0)  : subcooled (ΔT wrt Tsat, negative)
#   x in [0, 1]           : two-phase (quality Q)
#   x in (1, 1+superheat] : superheated (1 + ΔT wrt Tsat, positive)
# -------------------------------
# Subcooling segment (exclude exactly 0 to avoid saturation)
x_sub = np.linspace(-subcool_K, -1e-6, 40) if subcool_K > 0 else np.array([])
T_sub = T_sat + x_sub  # since x_sub is negative K
# Two-phase quality segment
Q_vals = np.linspace(0.01, 0.99, 60)  # avoid Q=0 and Q=1 to keep well inside dome
x_tp = np.linspace(0.0, 1.0, Q_vals.size)
# Superheat segment
x_sup = 1.0 + np.linspace(1e-6, superheat_K, 60) if superheat_K > 0 else np.array([])
T_sup = T_sat + (x_sup - 1.0)

# -------------------------------
# Compute speeds of sound
# -------------------------------
a_sub = []
for Ti in T_sub:
    # Compressed liquid at same evaporating pressure
    a_liq = CP.PropsSI("A", "T", Ti, "P", p_evap + dp_side, fluid)
    a_sub.append(a_liq)
a_sub = np.array(a_sub)

a_tp_hem = []
a_tp_frozen = []
for Qi in Q_vals:
    a_hem = hem_speed_of_sound_TQ(fluid, T_sat, Qi, dp=dp_side)
    a_tp_hem.append(a_hem)
    # Optional: frozen model for comparison
    a_fr = frozen_speed_of_sound_TQ(fluid, T_sat, Qi, dp=dp_side)
    a_tp_frozen.append(a_fr)
a_tp_hem = np.array(a_tp_hem)
a_tp_frozen = np.array(a_tp_frozen)

a_sup = []
for Ti in T_sup:
    # Superheated vapor at same evaporating pressure
    a_vap = CP.PropsSI("A", "T", Ti, "P", p_evap - dp_side, fluid)
    a_sup.append(a_vap)
a_sup = np.array(a_sup)

# -------------------------------
# Plot
# -------------------------------
plt.figure(figsize=(9.5, 5.8))

# Plot segments
if x_sub.size:
    plt.plot(x_sub, a_sub, label=f"Subcooled liquid (ΔT={subcool_K:.0f} K)", color="#1f77b4", lw=2)
plt.plot(x_tp, a_tp_hem, label="Two-phase (HEM)", color="#d62728", lw=2)
plt.plot(x_tp, a_tp_frozen, label="Two-phase (Frozen / Wood)", color="#9467bd", lw=2, ls="--")
if x_sup.size:
    plt.plot(x_sup, a_sup, label=f"Superheated vapor (ΔT={superheat_K:.0f} K)", color="#2ca02c", lw=2)

# Vertical markers at the boundaries
plt.axvline(0.0, color="k", lw=1, alpha=0.6)
plt.axvline(1.0, color="k", lw=1, alpha=0.6)

# Axis labels and annotations
plt.title(f"Speed of Sound across Subcool → Two-phase → Superheat at T_evap={T_evap_C:.1f} °C ({fluid})")
plt.ylabel("Speed of sound [m/s]")

# Custom x-axis ticks to clarify segments
ticks = []
labels = []

# Subcooling ticks (every 1 K)
if subcool_K > 0:
    ticks += list(np.arange(-subcool_K, 0.1, 1.0))
    labels += [f"{int(tick):d}" for tick in np.arange(-int(subcool_K), 1, 1)]

# Two-phase region ticks (Q)
ticks += [0.0, 0.25, 0.5, 0.75, 1.0]
labels += ["Q=0", "0.25", "0.5", "0.75", "1.0"]

# Superheat ticks (every 2 K)
if superheat_K > 0:
    sup_ticks = 1.0 + np.arange(0, superheat_K + 0.1, 2.0)
    ticks += list(sup_ticks)
    labels += [f"{int(t-1):d}" for t in sup_ticks]

plt.xticks(ticks, labels)
plt.xlabel("Left: Subcooling ΔT [K]   |   Middle: Quality Q   |   Right: Superheat ΔT [K]")

plt.grid(True, ls=":", alpha=0.6)
plt.legend()
plt.tight_layout()

# After plt.tight_layout()
output_file = "speed_of_sound_plot.svg"
plt.savefig(output_file, format="svg")
print(f"SVG plot saved as: {output_file}")

plt.show()