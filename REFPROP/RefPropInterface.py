"""
RefPropInterface.py
A robust, pythonic interface to REFPROP via CoolProp's AbstractState.

Highlights
----------
- Pure fluids and mixtures (by list, dict, or string)
- Set state by: PT, PH, PS, HS, TQ, PQ (mass-based)
- Safe accessors: will NEVER raise from state() / all_props(); return None if undefined
- Properties: P, T, rho, rhomolar, h/s/u, cp/cv (mass & molar), a, mu, k, sigma, Pr,
              Q, M, Ru, R_specific, Z, gamma, g, f, nu, alpha_th, extras (if available)
- Saturation helpers, limits, backend metadata

Units (SI)
----------
P [Pa], T [K], rho [kg/m^3], h/s/u [J/kg], cp/cv [J/kg/K], a [m/s],
mu [Pa·s], k [W/m/K], sigma [N/m], Pr [-], M [kg/mol], etc.

Setup
-----
- Install CoolProp: `pip install CoolProp`
- (Optional) Install REFPROP and set its folder via environment variable RPPREFIX,
  for example on Windows:

      import os
      os.environ["RPPREFIX"] = r"C:\Program Files (x86)\REFPROP"

  or pass it directly when constructing the interface:

      rp = RefpropFluid("CO2", backend="REFPROP",
                        rpprefix=r"C:\Program Files (x86)\REFPROP")

Notes
-----
- Use raw strings (prefix with 'r') or escaped backslashes for Windows paths to avoid
  invalid escape sequence warnings in Python 3.12+.
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Union

import CoolProp.CoolProp as CP


# ----------------------------- Utilities ------------------------------------ #

def _strip_refprop_prefix(name: str) -> str:
    """Normalize 'REFPROP::CO2' to 'CO2' so caller can pass either form."""
    prefix = "REFPROP::"
    return name[len(prefix):] if name.upper().startswith(prefix) else name


def _ensure_sum_to_one(zs: Sequence[float], tol: float = 1e-10) -> List[float]:
    s = float(sum(zs))
    if s <= 0.0:
        raise ValueError("Mixture composition must sum to > 0.")
    if abs(s - 1.0) > tol:
        zs = [z / s for z in zs]
    return list(zs)


def _safe_call(fn, *args, **kwargs) -> Optional[float]:
    """
    Call a function and return its float value, or None on any exception.
    """
    try:
        val = fn(*args, **kwargs)
        return float(val)
    except Exception:
        return None


# --------------------------- Main Interface --------------------------------- #

@dataclass
class BackendInfo:
    backend: str
    rpprefix: Optional[str]
    refprop_path: Optional[str]
    refprop_version: Optional[str]


class RefpropFluid:
    """
    High-level interface to CoolProp's AbstractState with REFPROP or HEOS backend.

    Examples
    --------
    rp = RefpropFluid("CO2", backend="REFPROP", rpprefix=r"C:\Program Files (x86)\REFPROP")
    rp.PT(6.9e6, 302.15)
    print(rp.state())         # safe snapshot (never raises)

    # Mixture
    rp = RefpropFluid(["R1234ZE(E)", "R32"], composition=[0.6, 0.4], backend="REFPROP")
    rp.PT(6.9e6, 302.15)
    print(rp.all_props())     # extended snapshot, still safe
    """

    def __init__(
        self,
        fluid: Union[str, Sequence[str]],
        composition: Optional[Union[Sequence[float], Dict[str, float]]] = None,
        *,
        backend: str = "REFPROP",
        rpprefix: Optional[str] = None,
    ) -> None:
        """
        Parameters
        ----------
        fluid:
            - Pure fluid name: "CO2", "R1234ZE(E)" (case must match REFPROP if backend='REFPROP')
            - Sequence of component names for mixtures: ["R32", "R125"]
            - 'REFPROP::' prefix allowed; ignored internally.
        composition:
            - For mixtures: list of mole fractions [x1, x2, ...] or dict {"R32":0.5, "R125":0.5}
            - For pure fluids: None or [1.0]
        backend:
            - "REFPROP" (default) or "HEOS" (CoolProp native EOS)
        rpprefix:
            - Optional absolute path to REFPROP folder; sets os.environ['RPPREFIX']
        """
        self.backend = backend.upper()
        if rpprefix:
            os.environ["RPPREFIX"] = rpprefix  # must be set before creating state

        # Normalize input names and composition
        if isinstance(fluid, str):
            name = _strip_refprop_prefix(fluid)
            names = [name]
        else:
            names = [_strip_refprop_prefix(n) for n in fluid]

        if composition is None:
            zs = [1.0] if len(names) == 1 else None
        elif isinstance(composition, dict):
            comp_map = { _strip_refprop_prefix(k): float(v) for k, v in composition.items() }
            try:
                zs = [comp_map[n] for n in names]
            except KeyError as ke:
                raise ValueError(f"Composition dict missing key for component '{ke.args[0]}'") from ke
        else:
            zs = list(map(float, composition))

        if len(names) > 1:
            if zs is None:
                raise ValueError("Mixture specified but no composition provided.")
            zs = _ensure_sum_to_one(zs)

        self._names: List[str] = names
        self._zs: Optional[List[float]] = zs

        # Construct AbstractState
        fluid_string = "&".join(self._names)
        try:
            self._as = CP.AbstractState(self.backend, fluid_string)
        except Exception as exc:
            msg = (
                f"Failed to initialize AbstractState('{self.backend}', '{fluid_string}').\n"
                f"Details: {exc}\n\n"
                "Tips:\n"
                "- Ensure CoolProp is installed (pip install CoolProp).\n"
                "- If using REFPROP: verify REFPROP is installed and RPPREFIX points to its folder.\n"
                "- Verify fluid/component names match the backend ('R1234ZE(E)' vs 'R1234ze(E)')."
            )
            raise RuntimeError(msg) from exc

        if self.is_mixture:
            self._as.set_mole_fractions(self._zs)  # type: ignore[arg-type]

    # ---------------------------- Introspection ---------------------------- #

    @property
    def is_mixture(self) -> bool:
        return len(self._names) > 1

    def backend_info(self) -> BackendInfo:
        """Return backend/path/version info (best-effort)."""
        def _safe_param(key: str) -> Optional[str]:
            try:
                return CP.get_global_param_string(key)
            except Exception:
                return None

        return BackendInfo(
            backend=self.backend,
            rpprefix=os.environ.get("RPPREFIX"),
            refprop_path=_safe_param("refprop_path"),
            refprop_version=_safe_param("REFPROP_version"),
        )

    # ----------------------------- Updaters -------------------------------- #

    def PT(self, P: float, T: float) -> "RefpropFluid":
        """Set state by Pressure [Pa] and Temperature [K]."""
        return self._update(CP.PT_INPUTS, P, T)

    def PH(self, P: float, h: float) -> "RefpropFluid":
        """Set state by Pressure [Pa] and Mass Enthalpy [J/kg]."""
        return self._update(CP.HmassP_INPUTS, h, P)

    def PS(self, P: float, s: float) -> "RefpropFluid":
        """Set state by Pressure [Pa] and Mass Entropy [J/kg/K]."""
        return self._update(CP.SmassP_INPUTS, s, P)

    def HS(self, h: float, s: float) -> "RefpropFluid":
        """Set state by Mass Enthalpy [J/kg] and Mass Entropy [J/kg/K]."""
        return self._update(CP.HmassSmass_INPUTS, h, s)

    def TQ(self, T: float, Q: float) -> "RefpropFluid":
        """Set saturation state by Temperature [K] and Quality [0..1]."""
        return self._update(CP.QT_INPUTS, Q, T)

    def PQ(self, P: float, Q: float) -> "RefpropFluid":
        """Set saturation state by Pressure [Pa] and Quality [0..1]."""
        return self._update(CP.QP_INPUTS, Q, P)

    def _update(self, inputs: int, v1: float, v2: float) -> "RefpropFluid":
        try:
            self._as.update(inputs, float(v1), float(v2))
        except Exception as exc:
            raise ValueError(f"State update failed (inputs={inputs}, v1={v1}, v2={v2}): {exc}")
        return self

    # ------------------------------ Safe Getters --------------------------- #
    # All getters below MUST NEVER raise; they return None on error/undefined.

    @property
    def T(self) -> Optional[float]: return _safe_call(self._as.T)
    @property
    def p(self) -> Optional[float]: return _safe_call(self._as.p)
    @property
    def rho(self) -> Optional[float]: return _safe_call(self._as.rhomass)
    @property
    def rhomolar(self) -> Optional[float]: return _safe_call(self._as.rhomolar)

    @property
    def Q(self) -> Optional[float]: return _safe_call(self._as.Q)  # -1 single, 0..1 two-phase

    # Energetics (mass)
    @property
    def h(self) -> Optional[float]: return _safe_call(self._as.hmass)
    @property
    def s(self) -> Optional[float]: return _safe_call(self._as.smass)
    @property
    def u(self) -> Optional[float]: return _safe_call(self._as.umass)
    @property
    def cp(self) -> Optional[float]: return _safe_call(self._as.cpmass)
    @property
    def cv(self) -> Optional[float]: return _safe_call(self._as.cvmass)

    # Energetics (molar)
    @property
    def hmolar(self) -> Optional[float]: return _safe_call(self._as.hmolar)
    @property
    def smolar(self) -> Optional[float]: return _safe_call(self._as.smolar)
    @property
    def umolar(self) -> Optional[float]: return _safe_call(self._as.umolar)
    @property
    def cpmolar(self) -> Optional[float]: return _safe_call(self._as.cpmolar)
    @property
    def cvmolar(self) -> Optional[float]: return _safe_call(self._as.cvmolar)

    # Transport / acoustic
    @property
    def a(self) -> Optional[float]: return _safe_call(self._as.speed_sound)
    @property
    def mu(self) -> Optional[float]: return _safe_call(self._as.viscosity)
    @property
    def k(self) -> Optional[float]: return _safe_call(self._as.conductivity)

    @property
    def sigma(self) -> Optional[float]:
        """
        Surface tension [N/m]; only defined in two-phase.
        Returns None if not in two-phase or not supported by backend.
        """
        Q = self.Q
        if Q is None or not (0.0 <= Q <= 1.0):
            return None
        return _safe_call(self._as.surface_tension)

    @property
    def Pr(self) -> Optional[float]: return _safe_call(self._as.Prandtl)

    # Composition & constants
    @property
    def molar_mass(self) -> Optional[float]: return _safe_call(self._as.molar_mass)

    @property
    def gas_constant_molar(self) -> Optional[float]: return _safe_call(self._as.gas_constant)

    # Derived properties (safe)
    @property
    def gamma(self) -> Optional[float]:
        cp, cv = self.cp, self.cv
        try:
            if cp is None or cv in (None, 0.0):
                return None
            return cp / cv
        except Exception:
            return None

    @property
    def R_specific(self) -> Optional[float]:
        Ru, M = self.gas_constant_molar, self.molar_mass
        try:
            if Ru is None or M in (None, 0.0):
                return None
            return Ru / M
        except Exception:
            return None

    @property
    def Z(self) -> Optional[float]:
        """
        Compressibility factor Z = p / (rho * R_specific * T).
        Returns None in two-phase or if inputs are invalid.
        """
        try:
            Q = self.Q
            if Q is not None and 0.0 <= Q <= 1.0:
                return None
        except Exception:
            pass

        p, rho, Rs, T = self.p, self.rho, self.R_specific, self.T
        try:
            denom = (rho or 0.0) * (Rs or 0.0) * (T or 0.0)
            if p is None or denom == 0.0:
                return None
            return p / denom
        except Exception:
            return None

    @property
    def g(self) -> Optional[float]:
        T, s, h = self.T, self.s, self.h
        try:
            if None in (T, s, h):
                return None
            return h - T * s  # J/kg
        except Exception:
            return None

    @property
    def f(self) -> Optional[float]:
        T, s, u = self.T, self.s, self.u
        try:
            if None in (T, s, u):
                return None
            return u - T * s  # J/kg
        except Exception:
            return None

    @property
    def nu(self) -> Optional[float]:
        mu, rho = self.mu, self.rho
        try:
            if mu is None or rho in (None, 0.0):
                return None
            return mu / rho  # m^2/s
        except Exception:
            return None

    @property
    def alpha_th(self) -> Optional[float]:
        k, rho, cp = self.k, self.rho, self.cp
        try:
            denom = (rho or 0.0) * (cp or 0.0)
            if k is None or denom == 0.0:
                return None
            return k / denom  # m^2/s
        except Exception:
            return None

    def _try_method(self, name: str) -> Optional[float]:
        """
        Attempt to call a method on AbstractState if present. Returns None if unsupported.
        """
        try:
            meth = getattr(self._as, name)
        except AttributeError:
            return None
        return _safe_call(meth)

    # ------------------------------ Snapshots ------------------------------- #

    def state(self) -> Dict[str, Optional[float]]:
        """
        Curated state snapshot (SI units). Guaranteed not to raise.
        """
        return {
            "P_Pa": self.p,
            "T_K": self.T,
            "rho_kg_per_m3": self.rho,
            "rhomolar_mol_per_m3": self.rhomolar,
            "h_J_per_kg": self.h,
            "s_J_per_kgK": self.s,
            "u_J_per_kg": self.u,
            "cp_J_per_kgK": self.cp,
            "cv_J_per_kgK": self.cv,
            "hmolar_J_per_mol": self.hmolar,
            "smolar_J_per_molK": self.smolar,
            "umolar_J_per_mol": self.umolar,
            "cpmolar_J_per_molK": self.cpmolar,
            "cvmolar_J_per_molK": self.cvmolar,
            "a_m_per_s": self.a,
            "mu_Pa_s": self.mu,
            "k_W_per_mK": self.k,
            "sigma_N_per_m": self.sigma,  # None unless two-phase
            "Pr": self.Pr,
            "Q": self.Q,
            "M_kg_per_mol": self.molar_mass,
            "R_specific_J_per_kgK": self.R_specific,
            "Z": self.Z,
            "gamma_cp_over_cv": self.gamma,
            "g_J_per_kg": self.g,
            "f_J_per_kg": self.f,
            "nu_m2_per_s": self.nu,
            "alpha_th_m2_per_s": self.alpha_th,
            # Best-effort extras (if available)
            "isothermal_compressibility_1_Pa": self._try_method("isothermal_compressibility"),
            "isobaric_expansion_coeff_1_K": self._try_method("isobaric_expansion_coefficient"),
            "isentropic_compressibility_1_Pa": self._try_method("isentropic_compressibility"),
            "fundamental_derivative_of_gas_dynamics": self._try_method("fundamental_derivative_of_gas_dynamics"),
        }

    def all_props(self) -> Dict[str, Optional[float]]:
        """
        Extended snapshot including metadata. Guaranteed not to raise.
        """
        data = self.state().copy()

        # Friendly duplicates
        data.update({
            "P_bar": (data["P_Pa"] / 1e5) if data["P_Pa"] is not None else None,
            "T_C": (data["T_K"] - 273.15) if data["T_K"] is not None else None,
            "rho_mol_per_l": (data["rhomolar_mol_per_m3"] / 1000.0) if data["rhomolar_mol_per_m3"] is not None else None,
        })

        # Phase index
        try:
            data["phase_index"] = int(self._as.phase())
        except Exception:
            data["phase_index"] = None

        # Limits / critical (safe)
        def _add_limit(key: str, fn_name: str):
            data[key] = _safe_call(getattr(self._as, fn_name))

        _add_limit("T_min_K", "Tmin")
        _add_limit("T_max_K", "Tmax")
        _add_limit("p_max_Pa", "pmax")
        _add_limit("T_critical_K", "T_critical")
        _add_limit("p_critical_Pa", "p_critical")
        _add_limit("rho_critical_kg_per_m3", "rhomass_critical")

        return data

    # ------------------------------ Saturation ------------------------------ #

    def saturation_at_T(self, T: float) -> Dict[str, Dict[str, Optional[float]]]:
        """
        Return saturation properties at a given temperature for Q=0 (liq) and Q=1 (vap).
        Guaranteed not to raise.
        """
        out: Dict[str, Dict[str, Optional[float]]] = {}
        for Q, label in [(0.0, "sat_liq"), (1.0, "sat_vap")]:
            try:
                self._as.update(CP.QT_INPUTS, Q, float(T))
                out[label] = self.state().copy()
            except Exception:
                out[label] = {"error": None}  # safe placeholder
        return out

    def saturation_at_P(self, P: float) -> Dict[str, Dict[str, Optional[float]]]:
        """
        Return saturation properties at a given pressure for Q=0 (liq) and Q=1 (vap).
        Guaranteed not to raise.
        """
        out: Dict[str, Dict[str, Optional[float]]] = {}
        for Q, label in [(0.0, "sat_liq"), (1.0, "sat_vap")]:
            try:
                self._as.update(CP.QP_INPUTS, Q, float(P))
                out[label] = self.state().copy()
            except Exception:
                out[label] = {"error": None}
        return out

    # ------------------------------ Composition ----------------------------- #

    def set_composition(self, composition: Union[Sequence[float], Dict[str, float]]) -> None:
        """Update mixture composition (mole fractions)."""
        if not self.is_mixture:
            raise ValueError("set_composition only allowed for mixtures.")
        if isinstance(composition, dict):
            comp_map = { _strip_refprop_prefix(k): float(v) for k, v in composition.items() }
            try:
                zs = [comp_map[n] for n in self._names]
            except KeyError as ke:
                raise ValueError(f"Composition dict missing key for component '{ke.args[0]}'") from ke
        else:
            zs = list(map(float, composition))
        zs = _ensure_sum_to_one(zs)
        self._as.set_mole_fractions(zs)
        self._zs = zs


# ------------------------- One-shot Convenience API ------------------------- #

def props(
    fluid: Union[str, Sequence[str]],
    *,
    backend: str = "REFPROP",
    rpprefix: Optional[str] = None,
    composition: Optional[Union[Sequence[float], Dict[str, float]]] = None,
    PT: Optional[Tuple[float, float]] = None,
    PH: Optional[Tuple[float, float]] = None,
    PS: Optional[Tuple[float, float]] = None,
    HS: Optional[Tuple[float, float]] = None,
    TQ: Optional[Tuple[float, float]] = None,
    PQ: Optional[Tuple[float, float]] = None,
) -> Dict[str, Optional[float]]:
    """
    One-call interface returning an extended properties dict (SI units).
    Guaranteed not to raise for retrieval.
    """
    rp = RefpropFluid(fluid, composition=composition, backend=backend, rpprefix=rpprefix)
    set_args = [x is not None for x in (PT, PH, PS, HS, TQ, PQ)]
    if sum(set_args) != 1:
        raise ValueError("Provide exactly one of: PT, PH, PS, HS, TQ, PQ")

    if PT: rp.PT(*PT)
    elif PH: rp.PH(*PH)
    elif PS: rp.PS(*PS)
    elif HS: rp.HS(*HS)
    elif TQ: rp.TQ(*TQ)
    elif PQ: rp.PQ(*PQ)

    return rp.all_props()


# Int_Test.py
from RefPropInterface import RefpropFluid

RPPATH = r"C:\Program Files (x86)\REFPROP"  # raw string avoids SyntaxWarning

def print_props(title: str, d: dict):
    print("\n" + title)
    print("-" * len(title))
    for k, v in d.items():
        print(f"{k:40s}: {v}")

def main():
    # Example 1: Pure CO2 (REFPROP), PT
    rp1 = RefpropFluid("CO2", backend="REFPROP", rpprefix=RPPATH)
    rp1.PT(6.9e6, 29.0 + 273.15)
    print_props("CO2 @ 6.9 MPa, 302.15 K (REFPROP)", rp1.all_props())

    # Example 2: Pure CO2 (HEOS), PT — safe even though sigma is undefined
    rp2 = RefpropFluid("CO2", backend="HEOS")
    rp2.PT(6.9e6, 29.0 + 273.15)
    print_props("CO2 @ 6.9 MPa, 302.15 K (HEOS)", rp2.all_props())

    # Example 3: Mixture 60% R1234ZE(E) + 40% R32 (REFPROP), PT
    rp3 = RefpropFluid(["R1234ZE(E)", "R32"], composition=[0.6, 0.4],
                       backend="REFPROP", rpprefix=RPPATH)
    rp3.PT(6.9e6, 302.15)
    print_props("Mix 60% R1234ZE(E) + 40% R32 (REFPROP)", rp3.all_props())

    # Example 4: Saturation — surface tension defined here (Q=0 or 1)
    rp4 = RefpropFluid("R134A", backend="REFPROP", rpprefix=RPPATH)
    sat = rp4.saturation_at_T(278.15)
    print_props("R134a Saturated Liquid @ 278.15 K", sat["sat_liq"])
    print_props("R134a Saturated Vapor  @ 278.15 K", sat["sat_vap"])

if __name__ == "__main__":
    main()