"""Hand calculations: Reynolds, Womersley, Poiseuille WSS, pressure gradient.

Closed-form engineering estimates that anchor the whole study; every
simulator result is compared against these numbers.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from .core import BloodProperties, FlowWaveform, VesselGeometry


@dataclass
class HandCalc:
    """All hand-calculated hemodynamics numbers at one operating point."""

    u_mean: float  # mean axial velocity (m/s)
    reynolds: float  # Re = rho U D / mu
    reynolds_peak: float  # Re at peak systolic flow
    womersley: float  # alpha = R sqrt(omega rho / mu)
    tau_poiseuille: float  # Poiseuille wall shear stress (Pa)
    dp_per_l: float  # Poiseuille pressure gradient (Pa/m)
    dp_vessel: float  # dp/L * straight-section length (Pa)
    tau_peak: float  # Poiseuille WSS at peak flow (Pa)
    q_peak: float  # peak flow rate (m^3/s)
    q_min: float  # minimum flow rate (m^3/s)


def flow_hand_calc(
    vessel: VesselGeometry,
    blood: BloodProperties,
    wf: FlowWaveform,
) -> HandCalc:
    """Hand-calculate the hemodynamic numbers for the nominal lumen."""
    q_mean = wf.mean_flow()
    a0 = vessel.a0
    u_mean = q_mean / a0
    re = blood.rho * u_mean * vessel.d0 / blood.mu_newtonian

    # peak/minimum of the waveform over one cycle
    n = 600
    dt = wf.period / n
    q_peak = q_min = wf.flow(0.0)
    for i in range(1, n):
        q = wf.flow(i * dt)
        q_peak = max(q_peak, q)
        q_min = min(q_min, q)

    u_peak = q_peak / a0
    re_peak = blood.rho * u_peak * vessel.d0 / blood.mu_newtonian
    alpha = vessel.r0 * math.sqrt(wf.omega * blood.rho / blood.mu_newtonian)

    # Poiseuille: tau = 4 mu U / R ; dp/dL = 8 mu Q / (pi R^4)
    tau = 4.0 * blood.mu_newtonian * u_mean / vessel.r0
    tau_peak = 4.0 * blood.mu_newtonian * u_peak / vessel.r0
    dp_per_l = 8.0 * blood.mu_newtonian * q_mean / (math.pi * vessel.r0**4)
    dp_vessel = dp_per_l * vessel.length

    return HandCalc(
        u_mean=u_mean,
        reynolds=re,
        reynolds_peak=re_peak,
        womersley=alpha,
        tau_poiseuille=tau,
        dp_per_l=dp_per_l,
        dp_vessel=dp_vessel,
        tau_peak=tau_peak,
        q_peak=q_peak,
        q_min=q_min,
    )
