from tabulate import tabulate
import math

# Try to import CoolProp
try:
    from CoolProp.CoolProp import PropsSI
except ImportError:
    PropsSI = None

# Yield strengths and prices (USD/m³)

MATERIAL_PROPERTIES = {
    "SS304": {"yield_strength": 215, "price_usd_m3": 5250},  # avg of 3500–7000
    "SS316": {"yield_strength": 205, "price_usd_m3": 6500},  # avg of 4500–8500
    "Aluminum 6061": {"yield_strength": 276, "price_usd_m3": 6000},
    "Copper": {"yield_strength": 70, "price_usd_m3": 31000},
    "Polyethylene": {"yield_strength": 20, "price_usd_m3": 1400},
    "Polypropylene": {"yield_strength": 25, "price_usd_m3": 1350},
    "Nylon 6": {"yield_strength": 45, "price_usd_m3": 2350},
    "PVC": {"yield_strength": 52, "price_usd_m3": 1700},
    "Titanium": {"yield_strength": 240, "price_usd_m3": 42000},  # Grade 1[1](https://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MTU010)
    "Silver": {"yield_strength": 45, "price_usd_m3": 310000},  # avg price ~$800/kg, density ~10490 kg/m³[2](https://material-properties.org/silver-properties-applications-price-production/)
    "Gold": {"yield_strength": 205, "price_usd_m3": 961000},  # price ~$49900/kg, density ~19300 kg/m³[3](https://material-properties.org/gold-properties-applications-price-production/)
}

Piping_Costing = {
    "material": {"as % of cost":62.5},
    "installation": {"as % of cost":17.5},
    "flush&balance":{"as % of cost":17},
    "other":{"as % of cost":3}
    }



# --- Pump sizing constants (NEW) ---
PIPE_TOTAL_LENGTH_M = 1000.0     # <<< total loop length for pump sizing (1 km)
MAX_VELOCITY_M_S    = 2.0        # <<< max allowable fluid velocity
# PipeLenghtPerMW = ((15*2+16+5)*2)/2 # Pipe length forward and return per MW, in drawing we have 2.2 MW roughly
"""UNIVERSAL VARIABLES"""
T = 25 + 273.15
T_amb = 25 + 273.15
P_amb = 1e5
dT = 10
P_amb = 1e5
area = 1
PipeLenghtPerMW = 1
v_limit_liq = 2
Q_liq = 0.0 # to assign density for liquid refrigerant supply

""" Toggle here :  """
Q = 0.5
# Q = 0.8
"""____________________________________"""

""" Toggle here :  """
# Direction = "Supply"
Direction = "Return"
"""____________________________________"""


""" Toggle here :  """
# ref2Plot = "CO2"
# ref2Plot = "R1234ze(E)"
# ref2Plot = "Ammonia"
ref2Plot = "R1233zd(E)"
"""____________________________________"""


""" Toggle here :  """
# Max_V_Overwrite = False
Max_V_Overwrite, Max_V_ref  = True, 5
"""____________________________________"""

PickMaterial = "SS304"
# PickMaterial = "SS316"
# PickMaterial = "Aluminum 6061"
# PickMaterial = "Copper"
# PickMaterial = "Polyethylene"
# PickMaterial = "Polypropylene"
# PickMaterial = "Nylon 6"
# PickMaterial = "PVC"
# PickMaterial = "Titanium"
# PickMaterial = "Silver"
# PickMaterial = "Gold"


SheetName =[str(PickMaterial) + ' Q=' + str(Q) +" " + str(Direction), 'PartCost Q=' + str(Q) + str(Direction)]


def fetch_yield_strength(material):
    return MATERIAL_PROPERTIES.get(material, {}).get("yield_strength")

def fetch_material_price(material):
    return MATERIAL_PROPERTIES.get(material, {}).get("price_usd_m3")

def calculate_required_thickness(material, pressure_bar, phase, outer_diameter_mm):
    allowable_stress = fetch_yield_strength(material)
    if allowable_stress is None:
        return "Material not found"
    safety_factor = 3 if phase == '1Ø' else 5
    pressure_mpa = pressure_bar / 10
    required_thickness_mm = (pressure_mpa * outer_diameter_mm) / (2 * (allowable_stress / safety_factor))
    return required_thickness_mm


def FloatingGCref(Tamb, dTgc):
    A0 = 4.81864e7
    A1 = -4.21537e5
    A2 = 9.44047e2
    alpha = 0.1
    Pref = 67e5
    T100 = 39.0 + 273.15
    dTsub = 3
    HRMODE = 0
    if HRMODE < 2:
        T = min(50, max(7, Tamb)) + 273.15
        Tref = ((A1 * A1 - 4.0 * (A0 - Pref) * A2) ** 0.5 - A1) / (2.0 * A2)
        tmp1 = T100 + dTsub
        tmp2 = T100 - Tref
        beta = (1e7 - A0 - (A1 + A2 * tmp1) * tmp1) / (tmp2 + (tmp2 * tmp2 + alpha * alpha) ** 0.5)
        b0 = A0 + (A1 + A2 * dTsub) * dTsub - beta * Tref
        b1 = A1 + 2 * A2 * dTsub + beta
        b2 = A2
        delta = T - Tref
        root = (delta * delta + alpha * alpha) ** 0.5
        P = b0 + (b1 + b2 * T) * T + beta * root
        Sgc_ref = T + dTgc - 273.15
        Pgc_ref = P * 1e-5
        return Sgc_ref, Pgc_ref
    else:
        return None, None


_, P_CO2 = FloatingGCref(T, 1)

class Refrigerant:
    def __init__(self, name, gwp, phase, pfas, Rclass, Vlimit,RCL, material=PickMaterial, pipe_length_m=PipeLenghtPerMW):
        self.name = name
        self.gwp = gwp
        self.phase = phase
        self.pfas = pfas
        self.Rclass = Rclass

        if Max_V_Overwrite and self.phase=="2Ø":
            self.Vlimit = Max_V_ref
        else:
            self.Vlimit = Vlimit

        self.RCL = RCL
        self.material = material
        self.pipe_length_m = pipe_length_m
        self.psat_25C = self.calculate_psat()
        self.pressure_used_bar = self.determine_pressure_used()
        self.pipe_schedule = self.pipe_schedule()
        self.kW_per_m2 = self.calc_kW_per_m2()
        self.area_for_1MW_m2 = self.calculate_area_for_1MW()
        self.pipe_diameter_mm = self.calculate_pipe_diameter()
        self.pipe_thickness_mm = self.calculate_pipe_thickness()
        self.pipe_cost_usd = self.calculate_pipe_cost()
        self.pump_power_kW = self.calculate_pump_power()
        self.mu = self.calc_mu()
        self.rho = self.calc_rho()
        self.sigma = self.calc_surface_tension()
        self.CHF_Zuber = self.calc_chf_zuber()
        self.rho_l_v_rat = self.calc_rho_l_v_rat()



    def calc_mu(self):
        if self.phase == '1Ø':
            mu_amb = PropsSI('VISCOSITY', "T", T_amb, "P", P_amb, self.name)

            return mu_amb
        elif self.phase == '2Ø':
            mu_amb = PropsSI('VISCOSITY', 'T', T_amb, 'Q', 0.5, self.name)

            return mu_amb
        else:
            raise ValueError("Phase must be '1Ø' or '2Ø'")

    def calc_rho(self):
        if self.phase == '1Ø':
            rho_amb = PropsSI("D", "T", T_amb, "P", P_amb, self.name)
            return rho_amb

        elif self.phase == '2Ø':
            rho_amb = PropsSI('D', 'T', T_amb, 'Q', 0.5, self.name)
            return rho_amb
        else:
            raise ValueError("Phase must be '1Ø' or '2Ø'")

    def calc_rho_l_v_rat(self):
        if self.phase == '1Ø':

            return 1

        elif self.phase == '2Ø':
            if self.name=="CO2" and T>PropsSI('Tcrit', self.name):
                rho_amb = PropsSI('D', 'P', P_CO2, 'T', T, self.name) / PropsSI('D', 'P', P_CO2, 'T', T+dT, self.name)
            else:
                rho_amb = PropsSI('D', 'T', T, 'Q', 0.0, self.name) / PropsSI('D', 'T', T, 'Q', 1, self.name)
            return rho_amb
        else:
            raise ValueError("Phase must be '1Ø' or '2Ø'")

    def enthalpy_evap_per_kg( self, T_min_C: float = 15.0, T_max_C: float = 50.0, Q_input: float = None, T_step_C: float = 1.0, ):
        """
        Δh(T) = H(T, Q) - H(T, 0.01) if T_triple < T < T_crit;
                else Δh = cp(T,P) * dT  (equivalent sensible cooling from T to T+dT)
        """
        # Pull Q and dT from globals (or instance)
        if Q_input is None:
            try:
                Q_val = Q
            except NameError:
                Q_val = getattr(self, "Q", None)
        else:
            Q_val = Q_input
        if Q_val is None:
            raise ValueError("Q not provided. Pass Q_input or define global Q or self.Q.")
        if not (0.0 <= Q_val <= 1.0):
            raise ValueError(f"Q must be in [0,1]. Got {Q_val}")

        try:
            dT_step = float(dT)
        except NameError:
            dT_step = 10.0

        from CoolProp.CoolProp import PropsSI

        Tc = PropsSI('Tcrit', self.name)
        Ttriple = PropsSI('Ttriple', self.name)

        # Pressure for sensible cp (above Tc or below Ttriple)
        try:
            P_used_bar = float(self.pressure_used_bar)
            P_used_Pa = P_used_bar * 1e5
        except Exception:
            try:
                P_used_Pa = float(P_amb)
            except NameError:
                P_used_Pa = 1e5

        if T_min_C > T_max_C:
            T_min_C, T_max_C = T_max_C, T_min_C
        if T_step_C <= 0:
            raise ValueError("T_step_C must be positive.")

        results = []
        Q_base = 0.01

        T_C = T_min_C
        while T_C <= T_max_C + 1e-12:
            T_K = T_C + 273.15
            in_two_phase = (T_K > Ttriple) and (T_K < Tc)

            if in_two_phase:
                if Q_val <= Q_base:
                    delta_h_kJ_per_kg = 0.0
                    mode = 'latent'
                else:
                    try:
                        H_Q = PropsSI('H', 'T', T_K, 'Q', float(Q_val), self.name)
                        H_base = PropsSI('H', 'T', T_K, 'Q', Q_base, self.name)
                        delta_h_kJ_per_kg = max(0.0, (H_Q - H_base) / 1000.0)
                        mode = 'latent'
                    except Exception:
                        # Fallback to sensible if CoolProp fails at edge
                        cp = PropsSI('Cpmass', 'T', T_K, 'P', P_used_Pa, self.name)
                        delta_h_kJ_per_kg = (cp * dT_step) / 1000.0
                        mode = 'sensible'
                valid_two_phase = True
            else:
                cp = PropsSI('Cpmass', 'T', T_K, 'P', P_used_Pa, self.name)
                delta_h_kJ_per_kg = (cp * dT_step) / 1000.0
                mode = 'sensible'
                valid_two_phase = False

            results.append({
                "T_C": float(T_C),
                "h_0_to_Q_kJ_per_kg": float(delta_h_kJ_per_kg),
                "valid_two_phase": bool(valid_two_phase),
                "mode": mode,
                "Q_base_used": Q_base if mode == 'latent' else None,
            })
            T_C += T_step_C

        return results

    def calc_surface_tension(self):

        if "INCOMP" in self.name:

            return None
        elif self.phase == '1Ø':
            return -1
        elif self.phase == '2Ø':
            sigma = PropsSI('SURFACE_TENSION', 'T', T_amb, 'Q', 0.5, self.name)

            return sigma
        else:
            raise ValueError("Phase must be '1Ø' or '2Ø'")

    def calc_chf_zuber(self):
        # Zuber pool-boiling CHF with optional density-ratio adjustment
        # -----------------------------------------------------------------------------
        # This expression estimates the **critical heat flux (CHF)** at which saturated
        # pool boiling on a horizontal surface transitions to film boiling, based on
        # Zuber’s hydrodynamic instability model. The base correlation (without the last
        # factor) was derived by Novak Zuber (1959):
        #     q_chf = C * h_fg * (rho_v)**0.5 * [sigma * g * (rho_l - rho_v)]**0.25  [1](https://www.osti.gov/servlets/purl/4175511)
        #
        # Here we include an additional *engineering adjustment* multiplying factor
        # (1 + rho_v/rho_l):
        #
        #     q_chf = C * h_fg * (rho_v ** 0.5) * ((sigma * g * (rho_l - rho_v)) ** 0.25) * (1 + rho_v/rho_l)
        #
        # Variables (SI units recommended):
        #   C       : empirical constant (~0.131 for water in Zuber’s original work)  [1](https://www.osti.gov/servlets/purl/4175511)
        #   h_fg    : latent heat of vaporization, J/kg
        #   rho_v   : saturated vapor density, kg/m^3
        #   rho_l   : saturated liquid density, kg/m^3
        #   sigma   : surface tension at saturation, N/m
        #   g       : gravitational acceleration, m/s^2
        #
        # Assumptions & applicability:
        #   • Saturated **pool boiling** on a relatively large, horizontal surface.
        #   • Hydrodynamic (Rayleigh–Taylor/Helmholtz) instability limit for CHF
        #     as in Zuber’s derivation.  [1](https://www.osti.gov/servlets/purl/4175511)
        #
        # Note on the (1 + rho_v/rho_l) factor:
        #   • This factor is **not** part of Zuber’s original correlation; it is sometimes
        #     used in practice to slightly increase CHF when vapor density is non-negligible
        #     (e.g., near critical conditions or for certain refrigerants). For water at
        #     low pressures, rho_v << rho_l, so the factor ≈ 1 and has negligible impact.
        #     Broad reviews of CHF models and common modifications are provided by
        #     Liang & Mudawar (2017).  [2](https://engineering.purdue.edu/mudawar/files/articles-all/2018/2018-01.pdf)
        #
        # References:
        #   • Zuber, N. (1959), *Hydrodynamic Aspects of Boiling Heat Transfer*, AECU‑4439
        #     (U.S. AEC/OSTI).  [1](https://www.osti.gov/servlets/purl/4175511)
        #   • Liang, G., & Mudawar, I. (2017), “Pool boiling critical heat flux (CHF) –
        #     Part 1: Review of mechanisms, models, and correlations,” *Int. J. Heat Mass
        #     Transfer*, 115, 1353–1370.  [2](https://engineering.purdue.edu/mudawar/files/articles-all/2018/2018-01.pdf)
        #
        # Implementation example:
        #   q_chf = C * h_fg * (rho_v ** 0.5) * ((sigma * g * (rho_l - rho_v)) ** 0.25) * (1.0 + rho_v/rho_l)

        if self.phase != '2Ø':

            return -1

        try:
            h_fg = PropsSI('H', 'T', T_amb, 'Q', 1, self.name) - PropsSI('H', 'T', T_amb, 'Q', 0, self.name)
            rho_l = PropsSI('D', 'T', T_amb, 'Q', 0, self.name)
            rho_v = PropsSI('D', 'T', T_amb, 'Q', 1, self.name)
            sigma = PropsSI('SURFACE_TENSION', 'T', T_amb, 'Q', 0.5, self.name)
            g = 9.81
            C = 0.131

            q_chf = C * h_fg * (rho_v ** 0.5) * ((sigma * g * (rho_l - rho_v)) ** 0.25)*(1+rho_v/rho_l)

            return q_chf
        except Exception as e:

            return None

    def calculate_psat(self):
        if PropsSI is None:
            return "CoolProp not available"
        try:
            pressure_pa = PropsSI("P", "T", T, "Q", 0, self.name)
            return pressure_pa / 1e5
        except:
            return "N/A"


    def plot_enthalpy_evap_svg(
            self,
            T_min_C: float = 15.0,
            T_max_C: float = 50.0,
            T_step_C: float = 1.0,
            Qs=None,
            show_tc: bool = True,
            shade_two_phase: bool = True,
            figsize=(8, 5),
            savepath: str = None,
            show: bool = True,
            transparent: bool = True,
            keep_text_as_text: bool = True,
            # --- NEW: PG25 (25% propylene glycol) flat ΔT line options ---
            show_pg25_dt10: bool = True,  # draw the reference line by default
            cp_pg25_kJ_per_kgK: float = 3.7,  # typical Cp for PG25 in HVAC range
            dt_pg25_C: float = 10.0,  # constant ΔT in °C
            pg25_style: dict | None = None,  # custom matplotlib style
            pg25_autoscale: bool = True,  # ensure the line is visible in y-limits
    ):
        """
        Plot Δh per kg vs temperature for the refrigerant (evaporation enthalpy increment),
        and optionally overlay a constant sensible-cooling reference line for PG25 with ΔT=10°C.

        The PG25 line uses Δh = c_p · ΔT with a constant c_p (default 3.7 kJ/(kg·K)).
        """
        import matplotlib.pyplot as plt
        from matplotlib import rc_context

        # Resolve Qs default
        if Qs is None:
            try:
                q_default = float(Q)  # global
            except NameError:
                q_default = getattr(self, "Q", None)
                if q_default is None:
                    raise ValueError("Provide Qs, or define global Q or self.Q.")
            Qs = [q_default]
        elif isinstance(Qs, (int, float)):
            Qs = [float(Qs)]
        else:
            Qs = [float(q) for q in Qs]

        # SVG-friendly settings
        rc_updates = {
            "svg.fonttype": "none" if keep_text_as_text else "path",
            "path.simplify": True,
            "path.simplify_threshold": 0.0,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
        }

        with rc_context(rc_updates):
            fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

            # Draw refrigerant curves
            for q in Qs:
                results = self.enthalpy_evap_per_kg(
                    T_min_C=T_min_C, T_max_C=T_max_C, Q_input=q, T_step_C=T_step_C
                )
                T_C = [r["T_C"] for r in results]
                dh = [r["h_0_to_Q_kJ_per_kg"] for r in results]
                ax.plot(T_C, dh, lw=2, label=f"Q = {q:g}")

            # Context (critical, two-phase)
            try:
                from CoolProp.CoolProp import PropsSI
                Tc = PropsSI('Tcrit', self.name) - 273.15
                Ttr = PropsSI('Ttriple', self.name) - 273.15
            except Exception:
                Tc = None
                Ttr = None

            if shade_two_phase and (Tc is not None and Ttr is not None):
                left = max(T_min_C, Ttr)
                right = min(T_max_C, Tc)
                if right > left:
                    ax.axvspan(left, right, color="tab:gray", alpha=0.10, label="Two-phase region")

            if show_tc and (Tc is not None) and (T_min_C <= Tc <= T_max_C):
                ax.axvline(Tc, color="k", ls="--", lw=1.2, alpha=0.7, label="Tcrit")

            # --- NEW: Add PG25 ΔT flat line ---
            if show_pg25_dt10:
                dh_pg25 = cp_pg25_kJ_per_kgK * dt_pg25_C  # kJ/kg
                style = dict(color="tab:green", ls="--", lw=2.0, zorder=5)
                if pg25_style:
                    style.update(pg25_style)

                ax.plot(
                    [T_min_C, T_max_C], [dh_pg25, dh_pg25],
                    label=f"PG25 ΔT={dt_pg25_C:g}°C (≈{dh_pg25:.1f} kJ/kg)",
                    **style
                )

                # Ensure the horizontal line is visible within y-limits
                if pg25_autoscale:
                    ymin, ymax = ax.get_ylim()
                    pad = max(0.05 * (ymax - ymin), 0.5)  # small padding
                    new_ymin = min(ymin, dh_pg25 - pad)
                    new_ymax = max(ymax, dh_pg25 + pad)
                    ax.set_ylim(new_ymin, new_ymax)

            # Labels and style
            ax.set_title(f"{self.name}: Δh from Q≈0.01 to Q over Temperature")
            ax.set_xlabel("Temperature [°C]")
            ax.set_ylabel("Δh [kJ/kg]")
            ax.grid(True, which="both", ls=":", alpha=0.5)
            ax.legend(loc="best", frameon=True)

            # Transparent backgrounds, optional
            if transparent:
                fig.patch.set_alpha(0.0)
                ax.set_facecolor((1, 1, 1, 0))

            # Save as SVG (if requested)
            if savepath:
                if not savepath.lower().endswith(".svg"):
                    savepath += ".svg"

                # Python 3.13+ safe, with fallback for older versions
                try:
                    from datetime import datetime, UTC
                except ImportError:
                    from datetime import datetime, timezone
                    UTC = timezone.utc  # type: ignore

                metadata = {
                    "Title": f"{self.name} Δh vs T",
                    "Description": f"Δh per kg from Q≈0.01 to Q over {T_min_C}–{T_max_C} °C",
                    "Creator": "Vahid Khorshidi",
                    "Date": datetime.now(UTC).isoformat(),
                }

                fig.savefig(
                    savepath,
                    format="svg",
                    bbox_inches="tight",
                    metadata=metadata,
                    transparent=transparent,
                )

            if show:
                plt.show()

            return fig, ax

    def plot_enthalpy_evap_2_svg(
            self,
            T_min_C: float = 15.0,
            T_max_C: float = 50.0,
            T_step_C: float = 1.0,
            Qs=None,
            show_tc: bool = True,
            shade_two_phase: bool = True,
            figsize=(8, 5),
            savepath: str | None = None,
            show: bool = True,
            transparent: bool = True,
            keep_text_as_text: bool = True,

            # --- NEW (kW/m² model inputs; default values so your call site remains unchanged) ---
            dT_C: float = 5.0,
            area_m2: float = 1.0,
            Direction: str = "Return",  # "Supply" or "Return"
            Q_liq: float = 0.01,
            P_amb_Pa: float = 101_325.0,
            P_CO2_bar: float = 80.0,
            v_limit_liq: float | None = None,
            Vlimit_m_s: float | None = None,  # if None, uses self.Vlimit

            # --- Optional: PG25 (25% PG) constant-ΔT reference line in kW/m² ---
            show_pg25_dt10: bool = True,  # draw the reference line by default
            cp_pg25_kJ_per_kgK: float = 3.7,  # constant Cp for reference line
            rho_pg25_kg_per_m3: float = 1030.0,
            dt_pg25_C: float = 10.0,  # ΔT for reference line
            pg25_style: dict | None = None,
            pg25_autoscale: bool = True,
            v_limit_l: float = v_limit_liq,
    ):
        """
        Plot cooling capacity per area (kW/m²) vs temperature using your kW/m² model (calc_kW_per_m2).
        - Handles 1Ø for Water/MEG-25% and 2Ø including transcritical per your logic.
        - Draws multiple curves for Q in Qs (for 2Ø).
        - Keeps SVG-friendly output, two-phase shading, critical line, and a PG25 flat reference in kW/m².

        NOTE: The original function name kept for drop-in compatibility, but the Y-axis is now kW/m².
        """

        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib import rc_context

        # ---- Resolve Qs default (qualities for 2Ø curves)
        if Qs is None:
            q_default = getattr(self, "Q", 0.5)
            Qs = [float(q_default)]
        elif isinstance(Qs, (int, float)):
            Qs = [float(Qs)]
        else:
            Qs = [float(q) for q in Qs]

        # ---- Resolve Vlimit
        if Vlimit_m_s is None:
            Vlimit_m_s = getattr(self, "Vlimit", 0.5)  # reasonable default

        # ---- SVG-friendly settings
        rc_updates = {
            "svg.fonttype": "none" if keep_text_as_text else "path",
            "path.simplify": True,
            "path.simplify_threshold": 0.0,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "legend.fontsize": 10,
        }

        # ---- Helpers
        def frange(start, stop, step):
            vals = []
            x = start
            while x <= stop + 1e-9:
                vals.append(x)
                x += step
            return vals

        def _explicit_kw_per_m2(
                T_C: float,
                Q: float,
        ) -> float | None:
            """Fallback explicit calculation (mirrors your calc_kW_per_m2 logic) -> kW/m²."""
            try:
                from CoolProp.CoolProp import PropsSI
            except Exception:
                raise RuntimeError("CoolProp is required for this calculation.")

            name = self.name
            phase = getattr(self, "phase", "2Ø")
            T = T_C + 273.15
            dT = dT_C
            P_CO2 = P_CO2_bar
            Vlimit = Vlimit_m_s
            if Vlimit is None:
                raise ValueError("Vlimit not defined on self; pass Vlimit_m_s or set self.Vlimit.")

            if phase == "1Ø" and name in ["INCOMP::MEG-25%", "Water"]:
                cp = PropsSI("C", "T", T + dT / 2.0, "P", P_amb_Pa, name)  # J/(kg·K)
                rho = PropsSI("D", "T", T, "P", P_amb_Pa, name)  # kg/m³
                cooling_power_W = area_m2 * Vlimit * rho * cp * dT  # W

                return cooling_power_W / 1000.0 / area_m2

            if phase == "2Ø":
                try:
                    Tc = PropsSI('Tcrit', name)
                except Exception:
                    Tc = None

                if Tc is None or T < Tc - 1:  # subcritical two-phase region
                    if Direction == "Supply":
                        rho = PropsSI("D", "T", T, "Q", Q_liq, name)
                        volumetric_flow = (v_limit_liq if v_limit_liq is not None else Vlimit) * area_m2
                    else:
                        rho = PropsSI("D", "T", T, "Q", Q, name)
                        volumetric_flow = Vlimit * area_m2

                    mass_flow = rho * volumetric_flow
                    H_Q = PropsSI("H", "T", T, "Q", Q, name)
                    H_001 = PropsSI("H", "T", T, "Q", 0.01, name)
                    delta_h = H_Q - H_001  # J/kg
                    cooling_power_W = mass_flow * delta_h
                else:
                    # Transcritical: hold pressure at P_CO2 * 1e5 Pa
                    P_CO2_Pa = P_CO2 * 1e5
                    if Direction == "Supply":
                        rho = PropsSI("D", "T", T, "P", P_CO2_Pa, name)
                        volumetric_flow = Vlimit * area_m2
                    else:
                        rho = PropsSI("D", "T", T + dT, "P", P_CO2_Pa, name)
                        volumetric_flow = Vlimit * area_m2

                    mass_flow = rho * volumetric_flow
                    H_out = PropsSI("H", "T", T + dT, "P", P_CO2_Pa, name)
                    H_in = PropsSI("H", "T", T, "P", P_CO2_Pa, name)
                    delta_h = H_out - H_in
                    cooling_power_W = mass_flow * delta_h

                return cooling_power_W / 1000.0 / area_m2

            return None

        def _call_user_calc_kw_per_m2(
                T_C: float, Q: float
        ) -> float | None:
            """
            Try to call user's self.calc_kW_per_m2() by injecting the free variables it expects
            into its global namespace. Falls back to explicit model on error.
            """
            if not hasattr(self, "calc_kW_per_m2"):
                return _explicit_kw_per_m2(T_C, Q)

            fn = self.calc_kW_per_m2
            g = fn.__globals__
            # Stash existing globals we will override
            keys = ["T", "dT", "area", "Direction", "Q", "Q_liq", "P_amb", "P_CO2", "v_limit_liq"]
            stash = {k: g[k] for k in keys if k in g}

            try:
                g["T"] = T_C + 273.15
                g["dT"] = dT_C
                g["area"] = area_m2
                g["Direction"] = Direction
                g["Q"] = Q
                g["Q_liq"] = Q_liq
                g["P_amb"] = P_amb_Pa
                g["P_CO2"] = P_CO2_bar
                g["v_limit_liq"] = v_limit_liq
                # Ensure self.Vlimit exists or pass via attribute for your method
                if getattr(self, "Vlimit", None) is None and Vlimit_m_s is not None:
                    setattr(self, "Vlimit", Vlimit_m_s)

                val = fn()  # returns kW (in your method) or kW/m²? — your code returns kW; normalize:
                # Your method returns kW total (cooling_power/1000). Normalize to kW/m²:
                if val is None:
                    return None
                return float(val) / float(area_m2)

            except Exception:
                # On any issue, use the explicit local model
                return _explicit_kw_per_m2(T_C, Q)
            finally:
                # Restore globals
                for k in keys:
                    if k in stash:
                        g[k] = stash[k]
                    elif k in g:
                        del g[k]

        # ---------- Plotting ----------
        with rc_context(rc_updates):
            fig, ax = plt.subplots(figsize=figsize, constrained_layout=True)

            Ts = frange(T_min_C, T_max_C, T_step_C)

            for Q in Qs:
                ys = []
                for T_C in Ts:
                    y = _call_user_calc_kw_per_m2(T_C, Q)  # kW/m²
                    ys.append(y if y is not None else np.nan)
                ax.plot(Ts, ys, lw=2, label=f"Q={Q:g}, ΔT={dT_C:g}°C, {Direction}")

            # ---- Context: critical temp & two-phase shading
            try:
                from CoolProp.CoolProp import PropsSI
                Tc = PropsSI('Tcrit', self.name) - 273.15
                Ttr = PropsSI('Ttriple', self.name) - 273.15
            except Exception:
                Tc = None
                Ttr = None

            if shade_two_phase and (Tc is not None and Ttr is not None):
                left = max(T_min_C, Ttr)
                right = min(T_max_C, Tc)
                if right > left:
                    ax.axvspan(left, right, color="tab:gray", alpha=0.10, label="Two-phase region")

            if show_tc and (Tc is not None) and (T_min_C <= Tc <= T_max_C):
                ax.axvline(Tc, color="k", ls="--", lw=1.2, alpha=0.7, label="Tcrit")

            # ---- PG25 constant-ΔT reference line (flat, kW/m²)
            if show_pg25_dt10:
                # q'' [kW/m²] = rho * Vlimit * cp_kJ/kg-K * ΔT
                V_for_pg = v_limit_l
                kwm2_pg25 = rho_pg25_kg_per_m3 * V_for_pg * cp_pg25_kJ_per_kgK * dt_pg25_C
                # print(kwm2_pg25, rho_pg25_kg_per_m3 , V_for_pg ,cp_pg25_kJ_per_kgK , dt_pg25_C)
                style = dict(color="tab:green", ls="--", lw=2.0, zorder=5)
                if pg25_style:
                    style.update(pg25_style)

                ax.plot(
                    [T_min_C, T_max_C], [kwm2_pg25, kwm2_pg25],
                    label=f"PG25 ΔT={dt_pg25_C:g}°C (≈{kwm2_pg25:.1f} kW/m²)",
                    **style
                )
                if pg25_autoscale:
                    ymin, ymax = ax.get_ylim()
                    pad = max(0.05 * (ymax - ymin), 0.5)
                    ax.set_ylim(min(ymin, kwm2_pg25 - pad), max(ymax, kwm2_pg25 + pad))

            # ---- Labels & styling
            ax.set_title(f"{self.name}: Cooling capacity per area vs Temperature")
            ax.set_xlabel("Temperature [°C]")
            ax.set_ylabel("Cooling capacity [kW/m²]")
            ax.grid(True, which="both", ls=":", alpha=0.5)
            ax.legend(loc="best", frameon=True)

            # ---- Transparent background, for SVG overlay workflows
            if transparent:
                fig.patch.set_alpha(0.0)
                ax.set_facecolor((1, 1, 1, 0))

            # ---- Save as SVG (if requested)
            if savepath:
                if not savepath.lower().endswith(".svg"):
                    savepath += ".svg"
                # Python 3.13+ UTC compatibility
                try:
                    from datetime import datetime, UTC
                except ImportError:
                    from datetime import datetime, timezone
                    UTC = timezone.utc  # type: ignore

                metadata = {
                    "Title": f"{self.name} kW/m² vs T",
                    "Description": f"Cooling capacity per area over {T_min_C}–{T_max_C} °C",
                    "Creator": "Vahid Khorshidi",
                    "Date": datetime.now(UTC).isoformat(),
                }

                fig.savefig(
                    savepath,
                    format="svg",
                    bbox_inches="tight",
                    metadata=metadata,
                    transparent=transparent,
                )

            if show:
                plt.show()

            return fig, ax

    def determine_pressure_used(self):
        if "water" in self.name.lower():
            # print(self.name, self.phase)
            return 3.0
        elif self.phase == "1Ø":
            # print(self.name, self.phase)
            return 3.0
        else:
            # print(self.name, self.phase)
            return self.psat_25C if isinstance(self.psat_25C, (int, float)) else P_CO2

    def pipe_schedule(self):
        if self.pressure_used_bar < 10:
            return "40"
        else:
            return "80"
        # Define the datasets

    def calc_kW_per_m2(self):


        if self.phase == "1Ø" and self.name in ["INCOMP::MEG-25%", "Water"]:
            cp = PropsSI("C", "T", T + dT / 2, "P", P_amb, self.name)
            rho = PropsSI("D", "T", T, "P", P_amb, self.name)
            cooling_power = area * self.Vlimit * rho * cp * dT
            return cooling_power / 1000

        if self.phase == "2Ø":
            if T < PropsSI('Tcrit', self.name)-1:
                if Direction == "Supply":
                    rho = PropsSI("D", "T", T, "Q", Q_liq, self.name)
                    volumetric_flow = v_limit_liq * area
                else:
                    rho = PropsSI("D", "T", T, "Q", Q, self.name)
                    volumetric_flow = self.Vlimit * area


                mass_flow = rho * volumetric_flow
                delta_h = PropsSI("H", "T", T, "Q", Q, self.name) - PropsSI("H", "T", T, "Q", 0.01, self.name)
                # print(self.name, delta_h)
                cooling_power = mass_flow * delta_h

            else:

                if Direction == "Supply":
                    rho = PropsSI("D", "T", T, "P", P_CO2*1e5, self.name)
                    volumetric_flow =  self.Vlimit * area
                else:
                    rho = PropsSI("D", "T", T+dT, "P", P_CO2*1e5, self.name)
                    volumetric_flow = self.Vlimit * area

                mass_flow = rho * volumetric_flow
                # Function to calculate P_gc
                delta_h = PropsSI("H", "T", T+dT, "P", P_CO2*1e5, self.name) - PropsSI("H", "T", T, "P", P_CO2*1e5, self.name)

                cooling_power = mass_flow * delta_h


            return cooling_power / 1000
        return None

    def calculate_area_for_1MW(self):
        try:
            if isinstance(self.kW_per_m2, (int, float)) and self.kW_per_m2 > 0:
                # print(round(1000 / self.kW_per_m2, 6))
                return 1000 / self.kW_per_m2
            else:
                return "N/A"
        except:
            return "N/A"

    def calculate_pipe_diameter(self):
        try:
            if isinstance(self.area_for_1MW_m2, (int, float)):
                diameter_m = math.sqrt(4 * self.area_for_1MW_m2 / math.pi)
                # print(diameter_m)
                return diameter_m * 1000
            else:
                return "N/A"
        except:
            return "N/A"

    def calculate_pipe_thickness(self):
        try:
            if isinstance(self.pipe_diameter_mm, (int, float)):
                return calculate_required_thickness(self.material, self.pressure_used_bar, self.phase, self.pipe_diameter_mm)
            else:
                return "N/A"
        except Exception as e:
            return f"Error: {e}"

    def interpolateDNShedule(self):
        schedule_40 = [
            (0.005588, 2.2), (0.00762, 2.5), (0.010668, 3.230769231), (0.01397, 3.666666667),
            (0.018796, 4.933333333), (0.024384, 5.333333333), (0.032512, 6.736842105),
            (0.0381, 7.5), (0.049276, 8.818181818), (0.058928, 8.285714286),
            (0.07366, 9.666666667), (0.085344, 10.5), (0.097282, 11.26470588),
            (0.122174, 12.65789474), (0.146304, 13.39534884), (0.193802, 15.26),
            (0.242824, 16.20338983), (0.289052, 16.49275362), (0.3175, 16.66666667),
            (0.363474, 17.03571429), (0.409702, 17.15957447), (0.455676, 17.41747573),
            (0.547624, 17.67213115)
        ]

        schedule_80 = [
            (0.006858, 3.857142857), (0.009144, 4), (0.012446, 5.444444444),
            (0.015748, 5.636363636), (0.020828, 7.454545455), (0.02667, 8.076923077),
            (0.035052, 9.857142857), (0.040894, 10.73333333), (0.052578, 13.8),
            (0.062738, 12.35), (0.077978, 13.95454545), (0.09017, 15.43478261),
            (0.102362, 16.79166667), (0.12827, 19.42307692), (0.154178, 21.67857143),
            (0.202692, 24.9375), (0.254508, 27.08108108), (0.303276, 29.12195122),
            (0.333502, 29.84090909), (0.381, 30), (0.428752, 30.14285714),
            (0.477774, 31.88135593), (0.574802, 32.79710145)
        ]

        dataset = schedule_40 if self.pipe_schedule == '40' else schedule_80 if self.pipe_schedule == '80' else None
        if dataset is None:
            raise ValueError("Invalid schedule type. Use '40' or '80'.")

        diameter = float(self.pipe_diameter_mm)/1000

        # Exact match
        for x, y in dataset:
            if diameter == x:
                return y

        # Interpolation
        for i in range(len(dataset) - 1):
            x0, y0 = dataset[i]
            x1, y1 = dataset[i + 1]
            if x0 <= diameter <= x1:
                return y0 + (diameter - x0) * (y1 - y0) / (x1 - x0)

        # Extrapolation
        if diameter < dataset[0][0]:
            x0, y0 = dataset[0]
            x1, y1 = dataset[1]
        elif diameter > dataset[-1][0]:
            x0, y0 = dataset[-2]
            x1, y1 = dataset[-1]
        else:
            raise ValueError("Unexpected error during extrapolation.")

        return y0 + (diameter - x0) * (y1 - y0) / (x1 - x0)

    def calculate_pipe_cost(self):

        try:

            price_per_m3 = fetch_material_price(self.material)
            if price_per_m3 is None or not isinstance(self.pipe_diameter_mm, (int, float)) or not isinstance(self.interpolateDNShedule(), (int, float)):
                return "N/A"

            outer_radius_m = self.pipe_diameter_mm / 2000
            inner_radius_m = outer_radius_m - (self.interpolateDNShedule() / 1000)
            cross_section_area_m2 = math.pi * (outer_radius_m**2 - inner_radius_m**2)
            volume_m3 = cross_section_area_m2 * self.pipe_length_m
            cost = volume_m3 * price_per_m3

            return cost*2
        except:
            return "N/A"

        # --- Friction/roughness helpers ------------------------------------------------

    def _roughness_from_material(self, fallback: float = 4.5e-5) -> float:
        """Map material → absolute roughness ε [m] (fallback to commercial steel)."""
        mat = str(getattr(self, "material", "")).strip().lower()
        map_exact = {
            "ss304": 1.5e-5, "ss316": 1.5e-5, "copper": 1.5e-5, "aluminum 6061": 1.5e-5,
            "polyethylene": 5.0e-6, "polypropylene": 5.0e-6, "nylon 6": 5.0e-6,
            "pvc": 1.5e-6, "titanium": 1.5e-5
        }
        map_keywords = {"stainless": 1.5e-5, "steel": 4.5e-5, "hdpe": 5.0e-6, "pe": 5.0e-6,
                        "cpvc": 2.0e-6, "smooth": 1.0e-6}
        if mat in map_exact: return map_exact[mat]
        for k, v in map_keywords.items():
            if k in mat: return v
        return fallback

    def _friction_factor(self, Re: float, rel_eps: float, method: str = "churchill") -> float:
        """Churchill default; Haaland optional. Handles laminar cap."""
        import math
        def f_churchill(Re, rel_eps):
            term = (7.0 / Re) ** 0.9 + 0.27 * rel_eps
            term = max(term, 1e-12)
            A = (2.457 * math.log(1.0 / term)) ** 16
            B = (37530.0 / Re) ** 16
            return 8.0 * (((8.0 / Re) ** 12) + 1.0 / ((A + B) ** 1.5)) ** (1.0 / 12.0)

        def f_haaland(Re, rel_eps):
            inner = (rel_eps / 3.7) ** 1.11 + 6.9 / Re
            inner = max(inner, 1e-12)
            return 1.0 / (-1.8 * math.log10(inner)) ** 2

        if Re < 2300:  # laminar
            return max(64.0 / max(Re, 1e-12), 1e-8)
        return max((f_haaland if method.lower() == "haaland" else f_churchill)(Re, rel_eps), 1e-8)

    def equivalent_length_for_delta_p( self, delta_p_bar: float = 4.0, mu_Pa_s: float = None, rho_kg_m3: float = None, velocity_m_per_s: float = 2.0, diameter_mm: float = 126.0, roughness_m: float = None, minor_loss_K: float = 0.0, method: str = "churchill",) -> float:
        """
        Compute the pipe length L_eq that produces the target Δp (bar) at given v and D.
        Default is the calibration condition: 4 bar, water-like properties, v=2 m/s, D=126 mm.
        """
        import math
        # If not provided, infer properties from CoolProp for water at (T, P_amb)
        if (mu_Pa_s is None) or (rho_kg_m3 is None):
            if PropsSI is not None:
                try:
                    rho_kg_m3 = PropsSI("D", "T", T, "P", P_amb, "Water") if rho_kg_m3 is None else rho_kg_m3
                    mu_Pa_s = PropsSI("V", "T", T, "P", P_amb, "Water") if mu_Pa_s is None else mu_Pa_s
                except Exception:
                    pass
        # Fallback typical values @ ~25°C if still None
        rho_kg_m3 = 997.0 if rho_kg_m3 is None else rho_kg_m3
        mu_Pa_s = 0.00089 if mu_Pa_s is None else mu_Pa_s

        D = self.pipe_diameter_mm / 1000.0
        roughness_m = self._roughness_from_material() if roughness_m is None else roughness_m
        area_m2 = math.pi * (D ** 2) / 4.0
        v = velocity_m_per_s
        Re = rho_kg_m3 * v * D / mu_Pa_s
        f = self._friction_factor(Re, roughness_m / D, method)
        dynamic_head = 0.5 * rho_kg_m3 * v * v
        delta_p_Pa = delta_p_bar * 1e5
        numerator = max(delta_p_Pa - minor_loss_K * dynamic_head, 0.0)
        L_eq = (D / f) * (numerator / dynamic_head) if dynamic_head > 0 else 0.0
        return L_eq

    # --- Pump power API ------------------------------------------------------------
    def _infer_mu_rho_for_hydraulics(self, use_liquid_for_2phase: bool = True):
        if PropsSI is None:
            raise ValueError("CoolProp not available to infer μ and ρ.")
        if self.phase == "1Ø":
            rho = PropsSI("D", "T", T, "P", P_amb, self.name)
            mu = PropsSI("V", "T", T, "P", P_amb, self.name)  # Pa·s
            return mu, rho
        if self.phase == "2Ø":
            if use_liquid_for_2phase:
                rho = PropsSI("D", "T", T, "Q", 0.0, self.name)
                mu = PropsSI("V", "T", T, "Q", 0.0, self.name)
                return mu, rho
            else:
                q_use = Q if Direction != "Supply" else Q_liq
                q_use = min(max(q_use, 0.0), 1.0)
                rho = PropsSI("D", "T", T, "Q", q_use, self.name)
                mu = PropsSI("V", "T", T, "Q", q_use, self.name)
                return mu, rho
        raise ValueError(f"Unsupported phase '{self.phase}' for μ/ρ inference.")

    # UPDATED: route to viscosity-based method with length=4 km
    def calculate_pump_power( self, delta_p_bar: float = None,   velocity_m_per_s: float = MAX_VELOCITY_M_S,   diameter_mm: float = None, pump_efficiency: float = 0.70, return_components: bool = False):
        """
        Pump power computed via Darcy–Weisbach over a fixed 4 km SS304 pipe,
        using the refrigerant's own 1 MW diameter and v = 2 m/s.
        """
        if diameter_mm is None:
            diameter_mm = self.pipe_diameter_mm if isinstance(self.pipe_diameter_mm, (int, float)) else 126.0

        return self.calculate_pump_power_from_viscosity(
            mu_Pa_s=None, rho_kg_m3=None,  # infer via CoolProp
            velocity_m_per_s=velocity_m_per_s,
            diameter_mm=diameter_mm,
            length_m=PIPE_TOTAL_LENGTH_M,  # <<< 4 km
            roughness_m=None,  # uses SS304 via self.material
            minor_loss_K=0.0,
            pump_efficiency=pump_efficiency,
            use_liquid_for_2phase=True,
            method="churchill",
            return_components=return_components,
            calibrate_to_4bar=False  # <<< do NOT calibrate to 4 bar
        )

    def _roughness_from_material(self, fallback: float = 4.5e-5) -> float:
        """
        Map your material labels to an absolute roughness ε [m].
        Defaults to commercial steel if unknown.
        """
        mat = str(getattr(self, "material", "")).strip().lower()
        map_exact = {
            "ss304": 1.5e-5,  # stainless steel (drawn/annealed ~1–2e-5 m)
            "ss316": 1.5e-5,
            "copper": 1.5e-5,
            "aluminum 6061": 1.5e-5,
            "polyethylene": 5.0e-6,  # HDPE/PE
            "polypropylene": 5.0e-6,
            "nylon 6": 5.0e-6,
            "pvc": 1.5e-6,
            "titanium": 1.5e-5
        }
        map_keywords = {
            "stainless": 1.5e-5,
            "steel": 4.5e-5,  # commercial steel
            "hdpe": 5.0e-6,
            "pe": 5.0e-6,
            "cpvc": 2.0e-6,
            "smooth": 1.0e-6,
        }
        if mat in map_exact:
            return map_exact[mat]
        for k, v in map_keywords.items():
            if k in mat:
                return v
        return fallback

    def calculate_pump_power_from_viscosity( self, mu_Pa_s: float = None, rho_kg_m3: float = None, velocity_m_per_s: float = MAX_VELOCITY_M_S,   diameter_mm: float = None,   length_m: float = PIPE_TOTAL_LENGTH_M,  roughness_m: float = None, minor_loss_K: float = 0.0, pump_efficiency: float = 0.70, use_liquid_for_2phase: bool = True, method: str = "churchill", return_components: bool = False, calibrate_to_4bar: bool = False,   calib_delta_p_bar: float = 4.0, calib_mu_Pa_s: float = None, calib_rho_kg_m3: float = None, calib_velocity_m_per_s: float = None, calib_diameter_mm: float = None,):
        """
        μ-based pump power using Darcy–Weisbach over a real length.
        """
        import math

        # Geometry
        if diameter_mm is None:
            diameter_mm = self.pipe_diameter_mm if isinstance(self.pipe_diameter_mm, (int, float)) else 126.0  # <<<
        D = diameter_mm / 1000.0
        roughness_m = self._roughness_from_material() if roughness_m is None else roughness_m
        area_m2 = math.pi * (D ** 2) / 4.0
        v = velocity_m_per_s
        Q = v * area_m2

        # Infer μ, ρ for the target fluid if missing
        if (mu_Pa_s is None) or (rho_kg_m3 is None):
            if PropsSI is not None:
                try:
                    if self.phase == "1Ø":
                        rho_kg_m3 = PropsSI("D", "T", T, "P", P_amb, self.name) if rho_kg_m3 is None else rho_kg_m3
                        mu_Pa_s = PropsSI("V", "T", T, "P", P_amb, self.name) if mu_Pa_s is None else mu_Pa_s

                    elif self.phase == "2Ø":
                        q_use = 0.0 if use_liquid_for_2phase else max(min(Q, 1.0), 0.0)
                        rho_kg_m3 = PropsSI("D", "T", T, "Q", q_use, self.name) if rho_kg_m3 is None else rho_kg_m3
                        mu_Pa_s = PropsSI("V", "T", T, "Q", q_use, self.name) if mu_Pa_s is None else mu_Pa_s
                except Exception:
                    pass
        if rho_kg_m3 is None or mu_Pa_s is None:
            raise ValueError("mu_Pa_s and rho_kg_m3 must be provided or inferable via CoolProp.")

        # (Optional) calibration to 4 bar (disabled by default now)
        if calibrate_to_4bar and (length_m is None):
            calib_velocity_m_per_s = v if calib_velocity_m_per_s is None else calib_velocity_m_per_s
            calib_diameter_mm = (diameter_mm if calib_diameter_mm is None else calib_diameter_mm)
            if (calib_mu_Pa_s is None) or (calib_rho_kg_m3 is None):
                if PropsSI is not None:
                    try:
                        calib_rho_kg_m3 = PropsSI("D", "T", T, "P", P_amb,
                                                  "Water") if calib_rho_kg_m3 is None else calib_rho_kg_m3
                        calib_mu_Pa_s = PropsSI("V", "T", T, "P", P_amb,
                                                "Water") if calib_mu_Pa_s is None else calib_mu_Pa_s
                    except Exception:
                        pass
            if calib_rho_kg_m3 is None: calib_rho_kg_m3 = 997.0
            if calib_mu_Pa_s is None: calib_mu_Pa_s = 0.00089
            length_m = self.equivalent_length_for_delta_p(
                delta_p_bar=calib_delta_p_bar,
                mu_Pa_s=calib_mu_Pa_s,
                rho_kg_m3=calib_rho_kg_m3,
                velocity_m_per_s=calib_velocity_m_per_s,
                diameter_mm=calib_diameter_mm,
                roughness_m=roughness_m,
                minor_loss_K=minor_loss_K,
                method=method,
            )
        # Pressure loss for the target fluid
        Re = rho_kg_m3 * v * D / mu_Pa_s
        f = self._friction_factor(Re, roughness_m / D, method)
        dynamic_head = 0.5 * rho_kg_m3 * v * v
        delta_p_major = f * (length_m / D) * dynamic_head
        delta_p_minor = minor_loss_K * dynamic_head
        delta_p_total = delta_p_major + delta_p_minor

        # Power
        P_hyd_W = delta_p_total * Q
        eta = min(max(pump_efficiency, 1e-6), 0.9999)
        P_shaft_W = P_hyd_W / eta

        result = {
            "mu_Pa_s": mu_Pa_s, "rho_kg_m3": rho_kg_m3,
            "Re": Re, "friction_factor": f,
            "velocity_m_per_s": v, "volumetric_flow_m3_s": Q,
            "delta_p_major_Pa": delta_p_major, "delta_p_minor_Pa": delta_p_minor,
            "delta_p_total_Pa": delta_p_total, "length_m": length_m,
            "hydraulic_power_kW": P_hyd_W / 1000.0, "shaft_power_kW": P_shaft_W / 1000.0,
            "diameter_m": D, "roughness_m": roughness_m,
        }
        return result if return_components else result["shaft_power_kW"]

    def calculate_pump_power( self, pump_efficiency: float = 0.70, return_components: bool = False, diameter_mm: float = None, velocity_m_per_s: float = None):
        """
        Calculates pump shaft power for delivering 1 MW cooling over 1 km of pipe.
        Accepts old parameters `diameter_mm` and `velocity_m_per_s` for backward compatibility, but ignores them.
        """
        import math

        # Always use calculated diameter for 1 MW
        if not isinstance(self.pipe_diameter_mm, (int, float)):
            raise ValueError("Pipe diameter for 1 MW not calculated.")
        D = self.pipe_diameter_mm / 1000.0
        area_m2 = math.pi * (D ** 2) / 4.0

        # Flow for 1 MW
        if not isinstance(self.kW_per_m2, (int, float)) or self.kW_per_m2 <= 0:
            raise ValueError("kW_per_m2 not calculated or invalid.")
        velocity_calc = (1000.0 / self.kW_per_m2) / self.area_for_1MW_m2
        volumetric_flow_m3_s = velocity_calc * area_m2

        # Fluid properties
        mu_Pa_s, rho_kg_m3 = self._infer_mu_rho_for_hydraulics()

        # Friction factor
        roughness_m = self._roughness_from_material()
        Re = rho_kg_m3 * velocity_calc * D / mu_Pa_s
        f = self._friction_factor(Re, roughness_m / D, "churchill")

        # Pressure drop for 1 km
        length_m = 1000.0
        dynamic_head = 0.5 * rho_kg_m3 * velocity_calc ** 2
        delta_p_major = f * (length_m / D) * dynamic_head
        delta_p_minor = 0.0
        delta_p_total = delta_p_major + delta_p_minor

        # Pump power
        P_hyd_W = delta_p_total * volumetric_flow_m3_s
        P_shaft_W = P_hyd_W / max(min(pump_efficiency, 0.9999), 1e-6)

        result = {
            "mu_Pa_s": mu_Pa_s,
            "rho_kg_m3": rho_kg_m3,
            "Re": Re,
            "friction_factor": f,
            "velocity_m_per_s": velocity_calc,
            "volumetric_flow_m3_s": volumetric_flow_m3_s,
            "delta_p_major_Pa": delta_p_major,
            "delta_p_minor_Pa": delta_p_minor,
            "delta_p_total_Pa": delta_p_total,
            "length_m": length_m,
            "hydraulic_power_kW": P_hyd_W / 1000.0,
            "shaft_power_kW": P_shaft_W / 1000.0,
            "diameter_m": D,
            "roughness_m": roughness_m,
        }
        return result if return_components else result["shaft_power_kW"]
    #

# Create refrigerant objects
refrigerants = [
# 1Ø
    Refrigerant('INCOMP::MEG-25%', 0, "1Ø", "No", "A1", v_limit_liq,"N/A"),
    Refrigerant("Water", 0, "1Ø", "No", "A1", v_limit_liq, "N/A"),

# 2Ø
    # Natural
    Refrigerant("CO2", 1, "2Ø", "No", "A1", 10.0,35.8),
    Refrigerant("Ammonia", 1.37, "2Ø", "No", "B2L", 25.0,0),
    Refrigerant("R290", 0.02, "2Ø", "No", "A3", 18.0,8),
    Refrigerant("Water", 0, "2Ø", "No", "A1", 30.0,"N/A"),

    # Low GWP
    Refrigerant("R1234ze(E)", 1.37, "2Ø", "Yes", "A2L", 20.0,76),
    Refrigerant("R1233zd(E)", 1, "2Ø", "No", "B1", 20.0,75),
    Refrigerant("R152a", 124, "2Ø", "Yes", "A2", 20.0,280),

    # Medium GWP
    Refrigerant("R32", 675, "2Ø", "Yes", "A2L", 20.0,77),
    Refrigerant("R134a", 1430, "2Ø", "Yes", "A1", 20.0,210),
    # Refrigerant("R515b", 1430, "2Ø", "Yes", "A1", 20.0),
    # Refrigerant("R245fa", 1430, "2Ø", "Yes", "A1", 20.0),
    # Refrigerant("R515b", 1430, "2Ø", "Yes", "A1", 20.0),
]



normalized_pipe_cost_water = 0
for r in refrigerants:
    if r.name=="Water" and r.phase=="1Ø":
        normalized_pipe_cost_water = r.pipe_cost_usd

table_data = []
error_refrigerants = []
# ---------- 1) Extend headers ----------
# ---------- extend headers (add both pump columns if not already done) ----------

if Direction == "Supply":
    headers = [
    "Media", "t CO2 eq.", "2Ø / 1Ø", "PFAS", "Ref. Class","RCL [gr./m3]", "Psat@25°C","Visc@25°C","Dens@25°C","SurfTens@25°C","CritHF@25°C",
    "V_lim liq."," Liq. rho", "P_@calc(bar)", "Pipe Schedule","Working fluid", "Cool. cap. / d_pipe [MW/m²] Liq",
    "Pipe Diameter (mm)", "Pipe Thickness Theoritical (mm)", "Pipe Thickness ANSI (mm)",
    "Pipe Tot. Cost for *"+str(PickMaterial)+"* (USD, 1m)", "Normalized cost",

        "Pump kW/1MW_cooling (Eta = 0.7, μ-based)" # keep if you implemented μ-based
]

else:
    headers = [
    "Media", "t CO2 eq.", "2Ø / 1Ø", "PFAS", "Ref. Class","RCL [gr./m3]", "Psat@25°C","Visc@25°C","Dens@25°C","SurfTens@25°C","CritHF@25°C",
    "V_lim Vap (L, 2Ø)","L/V rho", "P_@calc(bar)", "Pipe Schedule","Working fluid", "Cool. cap. / d_pipe [MW/m²] V+L",
    "Pipe Diameter (mm)", "Pipe Thickness Theoritical (mm)", "Pipe Thickness ANSI (mm)",
    "Pipe Tot. Cost for *"+str(PickMaterial)+"* (USD, 1m)", "Normalized cost"  # keep if you implemented μ-based
]


table_data = []
error_refrigerants = []
import matplotlib as mpl
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["font.sans-serif"] = ["Arial"]  # or "Calibri" if installed
path_svg = r"C:\Users\U375297\Git\Python\EnergySystems\Data Center\Plots"
for r in refrigerants:

    try:
        # Use each refrigerant's own LIQUID-line diameter for the 2 m/s baseline.
        diam_mm = r.pipe_diameter_mm if isinstance(r.pipe_diameter_mm, (int, float)) else 126.0
        if r.phase =="2Ø" and r.name==ref2Plot:
            r.plot_enthalpy_evap_svg(
                T_min_C=15, T_max_C=50, Qs=[0.5, 0.8],
                savepath=path_svg+"/enthalpy_evap_"+str(r.name)+"_V="+str(r.Vlimit)+".svg",
                transparent=True,
                keep_text_as_text=True
            )

            diam_mm = r.pipe_diameter_mm if isinstance(r.pipe_diameter_mm, (int, float)) else 126.0

            if r.phase == "2Ø" and r.name == ref2Plot:
                r.plot_enthalpy_evap_2_svg(
                    T_min_C=15, T_max_C=50, Qs=[0.5, 0.8],
                    savepath=path_svg+"/enthalpy_evap_2_" + str(r.name) +"_V="+str(r.Vlimit)+ ".svg",
                    transparent=True,
                    keep_text_as_text=True
                )



        # (A) Baseline: fixed Δp = 4 bar, v = 2 m/s, diameter = media-specific
        baseline_kW = r.calculate_pump_power(
            diameter_mm=diam_mm,
            velocity_m_per_s=2.0,
            pump_efficiency=0.70,
            return_components=False
        )

        # (B) μ-based (optional): same geometry as baseline; calibrate to 4 bar with water
        try:

            visc_kW = r.calculate_pump_power_from_viscosity(
                mu_Pa_s=None, rho_kg_m3=None,  # let CoolProp infer
                diameter_mm=diam_mm, velocity_m_per_s=2.0,
                length_m=4000,  # <— real length
                minor_loss_K=0.0,
                pump_efficiency=0.70,
                calibrate_to_4bar=False,  # <— turn OFF calibration
                return_components=False
            )

        except Exception:
            visc_kW = None

        if Direction == "Supply":
            row = [
                r.name, r.gwp, r.phase, r.pfas, r.Rclass,r.RCL,
                r.psat_25C,r.calc_mu(), r.calc_rho(),
                r.calc_surface_tension(),r.calc_chf_zuber(),
                f"{r.Vlimit} ({v_limit_liq}, {round(r.Vlimit / 3, 0)})", r.calc_rho_l_v_rat(),
                r.pressure_used_bar, r.pipe_schedule,r.name,
                r.kW_per_m2 / 1000,
                r.pipe_diameter_mm,
                r.pipe_thickness_mm, r.interpolateDNShedule(),
                r.pipe_cost_usd / (Piping_Costing.get("material", {}).get("as % of cost") / 100),
                r.pipe_cost_usd / normalized_pipe_cost_water,
                r.pump_power_kW,
            ]
        else:
            row = [
                r.name, r.gwp, r.phase, r.pfas, r.Rclass, r.RCL, r.psat_25C,r.calc_mu(), r.calc_rho(),
                r.calc_surface_tension(),r.calc_chf_zuber(),
                f"{r.Vlimit} ({v_limit_liq}, {round(r.Vlimit / 3, 0)})",r.calc_rho_l_v_rat(),
                r.pressure_used_bar, r.pipe_schedule,r.name,r.kW_per_m2 / 1000,
                 r.pipe_diameter_mm,
                r.pipe_thickness_mm, r.interpolateDNShedule(),
                r.pipe_cost_usd / (Piping_Costing.get("material", {}).get("as % of cost") / 100),
                r.pipe_cost_usd / normalized_pipe_cost_water
            ]

        table_data.append(row)

    except (TypeError, AttributeError) as e:
        error_name = getattr(r, 'name', 'Unknown Refrigerant')
        print(f"Error with refrigerant '{error_name}': {e}")
        error_refrigerants.append(error_name)
        continue

# Print the table

print(tabulate(table_data, headers=headers, tablefmt="grid", colalign=("center",)*len(headers)))



import os
import pandas as pd

df = pd.DataFrame(table_data, columns=headers)

# print(df)
# Create output folder if it doesn't exist
output_folder = r"C:\Users\U375297\Danfoss\AST System Solutions & Technology - Documents\Technology Projects\@ System & Efficiency\@ Data Center Thermal Management"

output_folder = r"C:\Users\U375297\Git\Python\EnergySystems\Data Center\Output"
os.makedirs(output_folder, exist_ok=True)

# Save the Excel file
excel_path = os.path.join(output_folder, "3 - 1Øvs2Ø.xlsx")
# df.to_excel(excel_path, index=False)

print(f"\nExcel file saved to: {excel_path}")


# Write to the "Updated" sheet without touching other sheets
with pd.ExcelWriter(excel_path, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
    df.to_excel(writer, sheet_name = SheetName[0], index=False)

from tabulate import tabulate

def calculate_component_costs():
    # Component data: (name, unit price EUR, quantity or length)
    components = [
        ("DN100 Pipe", 25, 35),
        ("DN50 Hose", 5, 64),
        ("T-Joints", 50, 32),
        ("Elbows DN100", 12, 4),
        ("Main Isolation Valves DN100", 500, 8),
        ("Rack Isolation Valves DN50", 400, 32),
        ("Weld Neck Flanges DN100", 24, 8)
    ]

    # Prepare table data
    table_data = []
    total_cost = 0

    for name, unit_price, quantity in components:
        cost = unit_price * quantity
        total_cost += cost
        table_data.append([name, unit_price, quantity, cost])

    # Add total cost row
    table_data.append(["Total", "", "", total_cost])



# Run the function
calculate_component_costs()

from tabulate import tabulate

# Component data: (name, unit price EUR, quantity or length)
components = [
    ("DN100 Pipe", 40, 35),
    ("DN50 Hose", 200, 64),
    ("T-Joints", 40, 32),
    ("Elbows DN100", 30, 4),
    ("Main Isolation Valves DN100", 1000, 8),
    ("Rack Isolation Valves DN50", 150, 32),
    ("Weld Neck Flanges DN100", 150, 8)
]

# Determine normalized pipe cost for water
normalized_pipe_cost_water = next((r.pipe_cost_usd for r in refrigerants if r.name == "Water"), 1)

# Prepare table data
table_data = []
total_cost = 0

for name, unit_price, quantity in components:
    cost = unit_price * quantity
    total_cost += cost
    row = [name, unit_price, quantity, cost]
    # Add normalized cost columns for each refrigerant
    for r in refrigerants:
        normalized_factor = r.pipe_cost_usd / normalized_pipe_cost_water
        row.append(cost * normalized_factor)
    table_data.append(row)

# Add total cost row
total_row = ["Total based on *SS304*", "", "", total_cost]
for r in refrigerants:
    normalized_factor = r.pipe_cost_usd / normalized_pipe_cost_water/(Piping_Costing.get("material", {}).get("as % of cost") / 100)
    total_row.append(total_cost * normalized_factor)
table_data.append(total_row)

# Define headers
headers = ["Component", "Unit Price (EUR)", "Quantity/Length", "Total(EUR)"]
headers += [f"({r.name})" for r in refrigerants]

# Print the table
print(tabulate(table_data, headers=headers, tablefmt="grid"))
df = pd.DataFrame(table_data, columns=headers)


os.makedirs(output_folder, exist_ok=True)

# Save the Excel file
excel_path = os.path.join(output_folder, "3 - 1Øvs2Ø.xlsx")

print(f"\nExcel file saved to: {excel_path}")

# Write to the "Updated" sheet without touching other sheets
with pd.ExcelWriter(excel_path, engine='openpyxl', mode='a', if_sheet_exists='replace') as writer:
    df.to_excel(writer, sheet_name=SheetName[1], index=False)


