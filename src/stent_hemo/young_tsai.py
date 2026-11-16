"""Stenosis pressure-drop correlation (Young-Tsai form).

Steady pressure drop across a stenosis modeled as a viscous Poiseuille
term through the throat plus an inertial (Borda-Carnot-type) expansion
loss from throat area back to nominal area. Young & Tsai (1973) measured
this two-term form experimentally for axisymmetric stenoses; the module
implements it in the standard form

    dp = 8 mu L_throat Q / (pi R_throat^4) + K * rho/2 * (Q/A_throat - Q/A_0)^2

with K an empirical loss coefficient calibrated against the Young-Tsai
measurements (K grows with severity; we use K = 1.52 as reported for the
high-severity regime and document the assumption).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .core import BloodProperties, VesselGeometry

K_YOUNG_TSAI = 1.52  # inertial loss coefficient (Young & Tsai 1973, high severity)


@dataclass
class StenosisDrop:
    """Pressure drop across the stenosis, split by mechanism."""

    dp_viscous: float  # Pa (Poiseuille through the throat)
    dp_inertial: float  # Pa (expansion loss)
    dp_total: float  # Pa
    v_throat: float  # throat mean velocity (m/s)


def stenosis_pressure_drop(
    vessel: VesselGeometry,
    blood: BloodProperties,
    q: float,
    k_loss: float = K_YOUNG_TSAI,
    throat_length: float | None = None,
) -> StenosisDrop:
    """Steady pressure drop (Pa) across the stenosis at flow rate q (m^3/s).

    Args:
        vessel: Vessel geometry (severity sets the throat).
        blood: Blood properties (Newtonian viscosity used).
        q: Flow rate (m^3/s).
        k_loss: Inertial loss coefficient (-).
        throat_length: Effective throat length (m); defaults to 0.4x the
            stenosis length (the region of near-throat radius).
    """
    r_t = vessel.r_throat
    a_t = vessel.a_throat
    a0 = vessel.a0
    l_t = throat_length if throat_length is not None else 0.4 * vessel.stenosis_length

    dp_visc = 8.0 * blood.mu_newtonian * q * l_t / (math.pi * r_t**4)
    v_t = q / a_t
    v_0 = q / a0
    dp_inert = k_loss * blood.rho / 2.0 * (v_t - v_0) ** 2
    return StenosisDrop(
        dp_viscous=dp_visc,
        dp_inertial=dp_inert,
        dp_total=dp_visc + dp_inert,
        v_throat=v_t,
    )
