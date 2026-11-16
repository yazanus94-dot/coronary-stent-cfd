"""Womersley analytic solution for pulsatile flow in a straight tube.

For a harmonic pressure gradient dp/dz = G * exp(i omega t), the velocity
profile is (Womersley 1955):

    u(r, t) = Re{ G / (i rho omega) * [1 - J0(alpha * i^(3/2) * r/R)
                                           / J0(alpha * i^(3/2))] * e^(i omega t) }

with alpha = R sqrt(omega rho / mu) and J0 the Bessel function of the
first kind (scipy.special.jv(0, z) accepts complex z).

Properties used for validation:
- alpha -> 0 recovers the parabolic Poiseuille profile.
- The flow rate driven by a real pressure gradient G has the classical
  complex admittance form; this module also provides the inverse map
  (velocity-profile formulation given a prescribed flow waveform).
"""

from __future__ import annotations

import cmath
import math

import numpy as np
from scipy.special import jv

from .core import BloodProperties, FlowWaveform, VesselGeometry

_I32 = 1j ** 1.5  # i^(3/2)


def _j0(z: complex) -> complex:
    """J0 Bessel function for a complex argument."""
    return complex(jv(0, z))


def womersley_profile(
    r_over_r: np.ndarray | float,
    t: float,
    alpha: float,
    rho: float,
    mu: float,
    omega: float,
    grad_amp: float,
    grad_phase: float = 0.0,
) -> np.ndarray | float:
    """Womersley axial velocity u(r, t) for one harmonic pressure gradient.

    Args:
        r_over_r: Normalized radius r/R (array or scalar), 0 center, 1 wall.
        t: Time (s).
        alpha: Womersley number (-).
        rho: Density (kg/m^3).
        mu: Viscosity (Pa s).
        omega: Angular frequency (rad/s).
        grad_amp: Pressure-gradient amplitude dp/dz (Pa/m).
        grad_phase: Pressure-gradient phase (rad).

    Returns:
        Velocity (m/s) with the same shape as r_over_r.
    """
    lam = alpha * _I32
    j_wall = _j0(lam)
    rr = np.asarray(r_over_r, dtype=float)
    shape_factor = 1.0 - np.array(
        [_j0(lam * r) / j_wall for r in np.atleast_1d(rr).astype(complex)]
    )
    u_complex = (
        grad_amp
        / (1j * rho * omega)
        * shape_factor
        * np.exp(1j * (omega * t + grad_phase))
    )
    out = u_complex.real
    if np.isscalar(r_over_r) or np.asarray(r_over_r).ndim == 0:
        return float(out[0])
    return out


def womersley_wss(
    alpha: float, rho: float, mu: float, omega: float, grad_amp: float
) -> float:
    """Wall shear stress amplitude (Pa) for one harmonic of the solution.

    tau_w = mu * du/dr |_{r=R}; from the Womersley form this evaluates to
    |G| * [J1(alpha i^{3/2}) / J0(alpha i^{3/2})] * (1 / (alpha i^{3/2})) in
    complex amplitude; the returned value is its magnitude.
    """
    lam = alpha * _I32
    ratio = _jv1(lam) / _j0(lam)
    # tau_w complex amplitude = -G * J1(lam) / (lam * J0(lam))
    tau_w = grad_amp * ratio / lam
    return abs(tau_w)


def _jv1(z: complex) -> complex:
    """J1 Bessel function for a complex argument."""
    return complex(jv(1, z))


def harmonic_pressure_gradient(
    wf: FlowWaveform, blood: BloodProperties, vessel: VesselGeometry
) -> list[tuple[int, float, float]]:
    """Fourier decomposition of the flow waveform -> pressure gradient.

    Uses the Womersley flow-rate admittance: a harmonic flow rate q_k
    e^(i omega_k t) is driven by the pressure gradient

        G_k = q_k / C_k,  C_k = pi R^2 / (i rho omega_k) * [1 - 2 J1(lam)/(lam J0(lam))]

    Returns a list of (harmonic k, amplitude Pa/m, phase rad) including the
    k=0 (steady Poiseuille) component first.
    """
    r0 = vessel.r0
    a0 = vessel.a0
    mu = blood.mu_newtonian
    rho = blood.rho

    # Numerical Fourier decomposition of Q(t) over one cycle.
    n = 512
    dt = wf.period / n
    ts = np.arange(n) * dt
    qs = np.array([wf.flow(t) for t in ts])
    q0 = float(qs.mean())

    grads: list[tuple[int, float, float]] = []
    # k = 0: steady Poiseuille, G0 = 8 mu Q / (pi R^4)
    grads.append((0, 8.0 * mu * q0 / (math.pi * r0**4), 0.0))

    for k, _ in enumerate(wf.coeffs, start=1):
        wk = k * wf.omega
        # complex Fourier coefficient of Q at harmonic k
        q_k = (2.0 / n) * np.sum(qs * np.exp(-1j * wk * ts)) * 1.0
        # q_k = mean(Q e^{-i w t}) * 2 -> matches cos decomposition
        lam = math.sqrt(wk * rho / mu) * r0 * _I32
        j0l = _j0(lam)
        j1l = _jv1(lam)
        admittance = (1j * rho * wk) / (
            1.0 - 2.0 * j1l / (lam * j0l)
        ) * (1.0 / (math.pi * r0**2))
        g_k = q_k * admittance
        grads.append((k, float(abs(g_k)), float(cmath.phase(g_k))))
    return grads


def womersley_case(
    vessel: VesselGeometry,
    blood: BloodProperties,
    wf: FlowWaveform,
    n_r: int = 25,
    n_t: int = 40,
):
    """Full Womersley velocity field over the vessel for the given waveform.

    Returns:
        (r_over_r array, t array, u array [n_t, n_r]) reconstructed from
        all waveform harmonics, plus the k=0 Poiseuille parabola.
    """
    rr = np.linspace(0.0, 1.0, n_r)
    ts = np.linspace(0.0, wf.period, n_t, endpoint=False)
    grads = harmonic_pressure_gradient(wf, blood, vessel)
    u = np.zeros((n_t, n_r))
    for k, amp, ph in grads:
        if k == 0:
            q0 = wf.mean_flow()
            g0 = 8.0 * blood.mu_newtonian * q0 / (math.pi * vessel.r0**4)
            # u = G0 R^2 (1 - rr^2) / (4 mu)  (parabolic, centerline = 2 U_mean)
            u += g0 * vessel.r0**2 * (1.0 - rr**2) / (4.0 * blood.mu_newtonian)
            continue
        alpha_k = vessel.r0 * math.sqrt(k * wf.omega * blood.rho / blood.mu_newtonian)
        for it, t in enumerate(ts):
            u[it, :] += np.asarray(
                womersley_profile(rr, t, alpha_k, blood.rho, blood.mu_newtonian, k * wf.omega, amp, ph)
            )
    return rr, ts, u
