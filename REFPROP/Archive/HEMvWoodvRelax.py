# -*- coding: utf-8 -*-
"""
Speed of sound: Subcool -> Two-phase -> Superheat
Two-phase region visually extended; single-phase regions visually compressed.

Models: Improved HEM (equilibrium), Frozen (Wood), Relaxation (frequency-dependent)
Saves SVG plot.
"""

import math
import numpy as np
import matplotlib.pyplot as plt
import CoolProp.CoolProp as CP

# ==========================
# Config (edit as needed)
# ==========================
fluid = "REFPROP::R1234ZE(E)"   # or "HEOS::R1234ZE(E)"
T_evap_C = 45.0                 # evaporating saturation temperature [°C]

# Thermal extents (actual thermodynamic ranges)
subcool_K = 2.0                 # was 5 K -> reduced to 2 K
superheat_K = 5.0               # was 10 K -> reduced to 5 K
# Visual widths (axis space allocation for each segment)
L_sub = 0.5                     # width of subcooled segment on x-axis (compressed)
L_tp  = 3.0                     # width of two-phase segment on x-axis (extended)
L_sup = 0.5                     # width of superheated segment on x-axis (compressed)

# Model params
dT_sat = 0.05                   # central-diff step along saturation [K] for improved HEM
freq_hz = 200.0                 # frequency for relaxation model [Hz]
tau_s = 0.005                   # relaxation time [s] (e.g., 5 ms)
svg_file = "speed_of_sound_HEM_Frozen_Relaxed_TP.svg"

# ==========================
# Helpers
# ==========================
def _sat_props_T(fluid, T):
    p = CP.PropsSI("P", "T", T, "Q", 0.5, fluid)
    rho_l = CP.PropsSI("D", "T", T, "Q", 0.0, fluid)
    rho_v = CP.PropsSI("D", "T", T, "Q", 1.0, fluid)
    v_l = 1.0 / rho_l
    v_v = 1.0 / rho_v
    s_l = CP.PropsSI("S", "T", T, "Q", 0.0, fluid)
    s_v = CP.PropsSI("S", "T", T, "Q", 1.0, fluid)
    return {"p": p, "v_l": v_l, "v_v": v_v, "s_l": s_l, "s_v": s_v}

def _central_diff(f_plus, f_minus, dT):
    return (f_plus - f_minus) / (2.0 * dT)

# ==========================
# Improved HEM (equilibrium)
# ==========================
def hem_speed_of_sound_eq_TQ(fluid, T, Q, dT=0.05):
    if not (0.0 < Q < 1.0):
        raise ValueError("Quality Q must be strictly between 0 and 1 for two-phase HEM.")
    if dT <= 0.0:
        raise ValueError("dT must be positive.")

    base = _sat_props_T(fluid, T)
    p0, v_l0, v_v0, s_l0, s_v0 = base["p"], base["v_l"], base["v_v"], base["s_l"], base["s_v"]

    ds = (s_v0 - s_l0)
    if abs(ds) < 1e-6:
        raise RuntimeError("Entropy difference between phases is too small (near critical).")

    plus = _sat_props_T(fluid, T + dT)
    minus = _sat_props_T(fluid, T - dT)

    dpdT  = _central_diff(plus["p"],   minus["p"],   dT)
    dvldT = _central_diff(plus["v_l"], minus["v_l"], dT)
    dvvdT = _central_diff(plus["v_v"], minus["v_v"], dT)
    dsldT = _central_diff(plus["s_l"], minus["s_l"], dT)
    dsvdT = _central_diff(plus["s_v"], minus["s_v"], dT)

    numerator = (1.0 - Q) * dsldT + Q * dsvdT
    dQdT = - numerator / ds

    v_mix = (1.0 - Q) * v_l0 + Q * v_v0
    dvdT  = (1.0 - Q) * dvldT + Q * dvvdT + (v_v0 - v_l0) * dQdT

    drhodT = - dvdT / (v_mix**2)
    if abs(dpdT) < 1e-12 or not math.isfinite(dpdT):
        raise RuntimeError("dp/dT along saturation is too small or invalid (near critical?).")
    drhodp_s = drhodT / dpdT

    if drhodp_s <= 0 or not math.isfinite(drhodp_s):
        raise RuntimeError(f"Non-physical or ill-conditioned dρ/dp|_s = {drhodp_s:.3e}.")
    a2 = 1.0 / drhodp_s
    if a2 <= 0 or not math.isfinite(a2):
        raise RuntimeError("Computed a^2 is invalid.")
    return math.sqrt(a2)

# ==========================
# Frozen (Wood)
# ==========================
def frozen_speed_of_sound_TQ(fluid, T, Q, dp=100.0):
    if not (0.0 < Q < 1.0):
        raise ValueError("Quality Q must be in (0,1).")
    p_sat = CP.PropsSI("P", "T", T, "Q", 0.5, fluid)
    rho_l = CP.PropsSI("D", "T", T, "Q", 0.0, fluid)
    rho_v = CP.PropsSI("D", "T", T, "Q", 1.0, fluid)

    alpha_v = (Q / rho_v) / (Q / rho_v + (1.0 - Q) / rho_l)
    alpha_l = 1.0 - alpha_v

    a_l = CP.PropsSI("A", "T", T, "P", p_sat + abs(dp), fluid)  # liquid side
    a_v = CP.PropsSI("A", "T", T, "P", p_sat - abs(dp), fluid)  # vapor side

    v_mix = (1.0 - Q) / rho_l + Q / rho_v
    rho_mix = 1.0 / v_mix

    inv_rho_a2 = alpha_l / (rho_l * a_l**2) + alpha_v / (rho_v * a_v**2)
    a2 = 1.0 / (rho_mix * inv_rho_a2)
    return math.sqrt(a2)

# ==========================
# Relaxation (frequency-dependent)
# ==========================
def relaxed_speed_of_sound_TQ(fluid, T, Q, freq_hz, tau_s, dT=0.05, dp_for_frozen=100.0):
    if not (0.0 < Q < 1.0):
        raise ValueError("Quality Q must be in (0,1).")

    rho_l = CP.PropsSI("D", "T", T, "Q", 0.0, fluid)
    rho_v = CP.PropsSI("D", "T", T, "Q", 1.0, fluid)
    v_mix = (1.0 - Q) / rho_l + Q / rho_v
    rho_mix = 1.0 / v_mix

    a_eq = hem_speed_of_sound_eq_TQ(fluid, T, Q, dT=dT)
    a_fr = frozen_speed_of_sound_TQ(fluid, T, Q, dp=dp_for_frozen)

    K_eq = rho_mix * a_eq**2
    K_fr = rho_mix * a_fr**2

    omega = 2.0 * math.pi * freq_hz
    jwt = 1j * omega * tau_s
    K_omega = K_eq + (K_fr - K_eq) * (jwt / (1.0 + jwt))
    a2_omega = K_omega / rho_mix
    return a2_omega**0.5  # complex

# ==========================
# 1) PRINT COMPARISON
# ==========================
T_sat = T_evap_C + 273.15
print("=== Comparison at T_sat = {:.2f} K ({:.2f} °C) ===".format(T_sat, T_evap_C))
print("Model params: freq = {:.1f} Hz, tau = {:.3e} s, dT_saturation = {:.3f} K".format(freq_hz, tau_s, dT_sat))
print("{:>6s}  {:>12s}  {:>12s}  {:>16s}".format("Q", "HEM [m/s]", "Frozen [m/s]", "Relax |a| [m/s]"))

for Q in [0.1, 0.3, 0.5, 0.7, 0.9]:
    a_eq = hem_speed_of_sound_eq_TQ(fluid, T_sat, Q, dT=dT_sat)
    a_fr = frozen_speed_of_sound_TQ(fluid, T_sat, Q, dp=100.0)
    a_rel_c = relaxed_speed_of_sound_TQ(fluid, T_sat, Q, freq_hz=freq_hz, tau_s=tau_s, dT=dT_sat)
    print("{:6.2f}  {:12.4f}  {:12.4f}  {:16.4f}".format(Q, a_eq, a_fr, abs(a_rel_c)))

# ==========================
# 2) PLOT with extended two-phase visual width
# ==========================
# Compute evaporating pressure and safe nudge
p_evap = CP.PropsSI("P", "T", T_sat, "Q", 0.5, fluid)
dp_side = max(50.0, 1e-6 * p_evap)

# ---- Build composite x-axis with controllable segment widths ----
# We'll map each physical segment to a visual segment of width L_sub, L_tp, L_sup.

# Subcooled segment (physical: ΔT ∈ [-subcool_K, 0], visual x ∈ [0, L_sub])
N_sub = max(6, int(6 * max(subcool_K, 1))) if subcool_K > 0 else 0
if N_sub > 0:
    t_sub = np.linspace(-subcool_K, 0.0, N_sub, endpoint=False)  # avoid 0 to not overlap with TP
    x_sub = np.linspace(0.0, L_sub, N_sub, endpoint=False)
    T_sub = T_sat + t_sub
else:
    x_sub = np.array([])
    T_sub = np.array([])

# Two-phase segment (physical: Q ∈ [0,1], visual x ∈ [L_sub, L_sub + L_tp])
N_tp = 100  # dense sampling for two-phase
Q_vals = np.linspace(0.02, 0.98, N_tp)  # avoid edges
x_tp = L_sub + (Q_vals - 0.0) * (L_tp / 1.0)

# Superheated segment (physical: ΔT ∈ (0, superheat_K], visual x ∈ (L_sub + L_tp, L_sub + L_tp + L_sup])
N_sup = max(6, int(6 * max(superheat_K, 1))) if superheat_K > 0 else 0
if N_sup > 0:
    t_sup = np.linspace(0.0, superheat_K, N_sup, endpoint=True)
    # Avoid duplicating the boundary point visually
    x_sup = L_sub + L_tp + np.linspace(0.0, L_sup, N_sup)
    T_sup = T_sat + t_sup
else:
    x_sup = np.array([])
    T_sup = np.array([])

# ---- Compute speeds ----
# Single-phase (subcooled liquid): call 'A' just off saturation pressure
if x_sub.size:
    a_sub = np.array([CP.PropsSI("A", "T", Ti, "P", p_evap + dp_side, fluid) for Ti in T_sub])
else:
    a_sub = np.array([])

# Two-phase models
a_tp_hem = []
a_tp_frozen = []
a_tp_relax_mag = []
for Qi in Q_vals:
    a_tp_hem.append(hem_speed_of_sound_eq_TQ(fluid, T_sat, Qi, dT=dT_sat))
    a_tp_frozen.append(frozen_speed_of_sound_TQ(fluid, T_sat, Qi, dp=100.0))
    a_tp_relax_mag.append(abs(relaxed_speed_of_sound_TQ(fluid, T_sat, Qi, freq_hz=freq_hz, tau_s=tau_s, dT=dT_sat)))
a_tp_hem = np.array(a_tp_hem)
a_tp_frozen = np.array(a_tp_frozen)
a_tp_relax_mag = np.array(a_tp_relax_mag)

# Single-phase (superheated vapor)
if x_sup.size:
    a_sup = np.array([CP.PropsSI("A", "T", Ti, "P", p_evap - dp_side, fluid) for Ti in T_sup])
else:
    a_sup = np.array([])

# ---- Plot ----
plt.figure(figsize=(11, 6.2))

if x_sub.size:
    plt.plot(x_sub, a_sub, label=f"Subcooled liquid (ΔT={subcool_K:.0f} K)", color="#1f77b4", lw=2)
plt.plot(x_tp, a_tp_hem, label="Two-phase (HEM - equilibrium)", color="#d62728", lw=2)
plt.plot(x_tp, a_tp_frozen, label="Two-phase (Frozen / Wood)", color="#9467bd", lw=2, ls="--")
plt.plot(x_tp, a_tp_relax_mag, label=f"Two-phase (Relaxation |a|, f={freq_hz:.0f} Hz, τ={tau_s*1e3:.0f} ms)", color="#ff7f0e", lw=2, ls="-.")

if x_sup.size:
    plt.plot(x_sup, a_sup, label=f"Superheated vapor (ΔT={superheat_K:.0f} K)", color="#2ca02c", lw=2)

# Boundaries
plt.axvline(L_sub, color="k", lw=1, alpha=0.6)
plt.axvline(L_sub + L_tp, color="k", lw=1, alpha=0.6)

# Labels
plt.title(f"Speed of Sound at T_evap={T_evap_C:.1f} °C — Extended Two-phase Region\n{fluid}")
plt.ylabel("Speed of sound [m/s]")

# Build custom x-ticks:
ticks = []
labels = []

# Subcool ticks: relative to T_sat (negative means subcooling)
if x_sub.size:
    ticks += list(np.linspace(0.0, L_sub, min(5, N_sub)))
    # map ticks to approximate ΔT labels
    tvals = np.linspace(-subcool_K, 0.0, min(5, N_sub), endpoint=False)
    labels += [f"{int(tv):d}" for tv in tvals]

# Two-phase ticks: show Q
tp_ticks = np.linspace(L_sub, L_sub + L_tp, 5)
ticks += list(tp_ticks)
labels += ["Q=0", "0.25", "0.5", "0.75", "1.0"]

# Superheat ticks: ΔT labels
if x_sup.size:
    ticks += list(np.linspace(L_sub + L_tp, L_sub + L_tp + L_sup, min(5, N_sup)))
    tvals = np.linspace(0.0, superheat_K, min(5, N_sup))
    labels += [f"{int(tv):d}" for tv in tvals]

plt.xticks(ticks, labels)
plt.xlabel(
    "Left: Subcooling ΔT [K]      |       Middle: Quality Q      |      Right: Superheat ΔT [K] "
)

plt.grid(True, ls=":", alpha=0.6)
plt.legend(loc="best")
plt.tight_layout()
plt.savefig(svg_file, format="svg")
print(f"\nSVG plot saved to: {svg_file}")
plt.show()
