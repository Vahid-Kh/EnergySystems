"""
refprop_interface.py
A high-level, pythonic interface to REFPROP via CoolProp's AbstractState.

Features
--------
- Pure fluids and mixtures (by list, dict, or string)
- Set state by common input pairs: PT, PH, PS, HS, TQ, PQ (mass-based properties)
- Accessors for: P, T, rho, h, s, u, cp, cv, a (speed of sound), mu (viscosity),
  k (thermal conductivity), Pr (Prandtl), Q (quality), Z (compressibility)
- Saturation helpers at given T or P
- Limits/critical properties
- Helpful errors and guidance for REFPROP setup

Units (SI)
----------
P [Pa], T [K], rho [kg/m^3], h/s/u [J/kg], a [m/s], mu [Pa·s], k [W/m/K]

Requirements
------------
pip install CoolProp
REFPROP optional (licensed), set via RPPREFIX or pass rpprefix=... in RefpropFluid()
"""

from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple, Union

import CoolProp.CoolProp as CP
# ----------------------------- Utilities ------------------------------------ #

def _strip_refprop_prefix(name: str) -> str:
    """Allow users to pass names like 'REFPROP::CO2' and normalize to 'CO2'."""
    prefix = "REFPROP::"
    return name[len(prefix):] if name.upper().startswith(prefix) else name
def _ensure_sum_to_one(zs: Sequence[float], tol: float = 1e-10) -> List[float]:
    s = float(sum(zs))
    if s <= 0.0:
        raise ValueError("Mixture composition must sum to > 0.")
    if abs(s - 1.0) > tol:
        # Normalize gently
        zs = [z / s for z in zs]
    return list(zs)


def _as_dict(items: Sequence[Tuple[str, float]]) -> Dict[str, float]:
    return {k: float(v) for k, v in items}


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
    Pure CO2, PT:
        rp = RefpropFluid("CO2", rpprefix=r"C:\\Program Files (x86)\\REFPROP")
        rp.PT(6.9e6, 29 + 273.15)
        print(rp.state())

    R1234ze(E)/R32 mixture, PT:
        rp = RefpropFluid(["R1234ZE(E)", "R32"], composition=[0.6, 0.4])
        rp.PT(6.9e6, 302.15)
        print(rp.rho, rp.a, rp.mu)

    Quality at saturation (TQ):
        rp = RefpropFluid("R134a")
        rp.TQ(278.15, Q=0.0)   # saturated liquid
        print(rp.p, rp.rho)
        rp.TQ(278.15, Q=1.0)   # saturated vapor
        print(rp.p, rp.rho)
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
            - You may also pass 'REFPROP::CO2' etc.; prefix will be ignored internally.
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
            # Ensure order matches names; if user passed dict, order by 'names'
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

        # Apply composition (if mixture)
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

    # ------------------------------ Getters -------------------------------- #

    # Base state
    @property
    def T(self) -> float: return self._as.T()

    @property
    def p(self) -> float: return self._as.p()

    @property
    def rho(self) -> float: return self._as.rhomass()

    @property
    def Q(self) -> float: return self._as.Q()  # -1 for single phase; 0..1 for two-phase

    # Energetics
    @property
    def h(self) -> float: return self._as.hmass()

    @property
    def s(self) -> float: return self._as.smass()

    @property
    def u(self) -> float: return self._as.umass()

    # Heat capacities
    @property
    def cp(self) -> float: return self._as.cpmass()

    @property
    def cv(self) -> float: return self._as.cvmass()

    # Transport / acoustic
    @property
    def a(self) -> float: return self._as.speed_sound()

    @property
    def mu(self) -> float: return self._as.viscosity()

    @property
    def k(self) -> float: return self._as.conductivity()

    @property
    def Pr(self) -> float: return self._as.Prandtl()

    # Other handy derived properties
    @property
    def Z(self) -> Optional[float]:
        """
        Compressibility factor Z = p / (rho * R_specific * T).

        Note:
            In two-phase (0 <= Q <= 1), Z is not well-defined as an ideal-gas-like factor,
            so this method returns None in that region.
        """
        try:
            Q = self.Q
            if 0.0 <= Q <= 1.0:
                return None  # not meaningful in two-phase mixture
        except Exception:
            # If Q evaluation failed, proceed to compute Z anyway for single-phase cases
            pass

        # Robust against CoolProp version differences:
        # Use AbstractState's gas_constant() [J/mol/K] and molar_mass() [kg/mol]
        R_u = self._as.gas_constant()        # universal gas constant [J/mol/K]
        M   = self._as.molar_mass()          # molar mass [kg/mol]
        if M <= 0.0:
            return None
        R_specific = R_u / M                 # [J/kg/K]
        denom = self.rho * R_specific * self.T
        if denom == 0.0:
            return None
        return self.p / denom

    # ------------------------------ Helpers -------------------------------- #

    def state(self) -> Dict[str, Optional[float]]:
        """Return a dict snapshot of the current state (SI)."""
        out: Dict[str, Optional[float]] = {
            "P_Pa": self.p,
            "T_K": self.T,
            "rho_kg_per_m3": self.rho,
            "h_J_per_kg": self.h,
            "s_J_per_kgK": self.s,
            "u_J_per_kg": self.u,
            "cp_J_per_kgK": self.cp,
            "cv_J_per_kgK": self.cv,
            "a_m_per_s": self.a,
            "mu_Pa_s": self.mu,
            "k_W_per_mK": self.k,
            "Pr": self.Pr,
            "Q": self.Q,
            "Z": self.Z,  # may be None in two-phase or degenerate cases
        }
        return out

    def limits(self) -> Dict[str, float]:
        """Return temperature/pressure/critical limits for this fluid/mixture."""
        return {
            "T_min_K": self._as.Tmin(),
            "T_max_K": self._as.Tmax(),
            "p_max_Pa": self._as.pmax(),
            "T_critical_K": self._as.T_critical(),
            "p_critical_Pa": self._as.p_critical(),
            "rho_critical_kg_per_m3": self._as.rhomass_critical(),
        }

    def saturation_at_T(self, T: float) -> Dict[str, Dict[str, Optional[float]]]:
        """
        Return saturation properties at a given temperature for Q=0 (liq) and Q=1 (vap).
        """
        out: Dict[str, Dict[str, Optional[float]]] = {}
        for Q, label in [(0.0, "sat_liq"), (1.0, "sat_vap")]:
            self._as.update(CP.QT_INPUTS, Q, float(T))
            out[label] = self.state().copy()
        return out
    def saturation_at_P(self, P: float) -> Dict[str, Dict[str, Optional[float]]]:
        """
        Return saturation properties at a given pressure for Q=0 (liq) and Q=1 (vap).
        """
        out: Dict[str, Dict[str, Optional[float]]] = {}
        for Q, label in [(0.0, "sat_liq"), (1.0, "sat_vap")]:
            self._as.update(CP.QP_INPUTS, Q, float(P))
            out[label] = self.state().copy()
        return out

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
    One-call interface returning a full state dict (SI units).

    Examples
    --------
    props("CO2", PT=(6.9e6, 302.15))
    props(["R32", "R125"], composition=[0.5, 0.5], PT=(1.2e6, 300.0))
    props("R134a", TQ=(278.15, 0.0))  # saturated liquid
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

    return rp.state()