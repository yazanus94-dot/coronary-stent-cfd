"""Core parameter definitions: blood, vessel, waveform, stent.

All SI units. These dataclasses are the single source of truth for every
downstream module (hand calc, Womersley, correlation, geometry, meshing,
solver case writers) so the whole study stays consistent.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class BloodProperties:
    """Blood rheology.

    Newtonian baseline plus the Carreau (Cho & Kensey 1991) model:
        mu(g) = mu_inf + (mu0 - mu_inf) * (1 + (lambda * g)^2)^((n - 1) / 2)
    with g the shear rate (1/s).

    Attributes:
        rho: Density (kg/m^3).
        mu_newtonian: Newtonian viscosity (Pa s).
        mu0: Zero-shear viscosity (Pa s).
        mu_inf: Infinite-shear viscosity (Pa s).
        lam: Carreau time constant lambda (s).
        n: Carreau power-law exponent (-).
    """

    rho: float = 1056.0
    mu_newtonian: float = 0.0035
    mu0: float = 0.056
    mu_inf: float = 0.00345
    lam: float = 3.313
    n: float = 0.3568

    def carreau_viscosity(self, shear_rate: float) -> float:
        """Carreau apparent viscosity (Pa s) at a shear rate (1/s)."""
        if shear_rate <= 0.0:
            return self.mu0
        return self.mu_inf + (self.mu0 - self.mu_inf) * (
            1.0 + (self.lam * shear_rate) ** 2
        ) ** ((self.n - 1.0) / 2.0)


@dataclass(frozen=True)
class VesselGeometry:
    """Idealized straight vessel with a cosine-shaped stenosis throat.

    The lumen radius follows r(x) = R0 - a * 0.5 * (1 + cos(pi * (x - xc) / Ls))
    inside the stenosis section, smoothly clamped elsewhere, so the throat
    radius is R0 - a. Severity is defined on area:
        severity = 1 - (R_throat / R0)^2   (fraction, e.g. 0.5 or 0.7)

    Attributes:
        d0: Nominal inner diameter (m), coronary scale ~3 mm.
        length: Total vessel length (m) including inlet/outlet extensions.
        stenosis_length: Axial length of the stenosis section (m).
        severity: Area-reduction severity in [0, 1).
        xc_frac: Axial center of the stenosis as a fraction of length.
    """

    d0: float = 3.0e-3
    length: float = 90.0e-3  # ~10 D upstream + 15 D downstream at 3 mm
    stenosis_length: float = 10.0e-3
    severity: float = 0.5
    xc_frac: float = 0.45

    @property
    def r0(self) -> float:
        """Nominal radius (m)."""
        return self.d0 / 2.0

    @property
    def r_throat(self) -> float:
        """Minimum lumen radius at the throat (m)."""
        return self.r0 * math.sqrt(1.0 - self.severity)

    @property
    def xc(self) -> float:
        """Axial position of the stenosis center (m)."""
        return self.xc_frac * self.length

    @property
    def a0(self) -> float:
        """Nominal lumen area (m^2)."""
        return math.pi * self.r0**2

    @property
    def a_throat(self) -> float:
        """Throat area (m^2)."""
        return math.pi * self.r_throat**2

    def radius_at(self, x):
        """Lumen radius (m) at axial position x (m); vectorized-safe."""
        import numpy as np

        half = self.stenosis_length / 2.0
        if isinstance(x, np.ndarray) and x.ndim > 0:
            bump = np.clip(
                0.5 * (1.0 + np.cos(np.pi * (x - self.xc) / half)),
                0.0,
                None,
            )
            bump = np.where(
                np.abs(x - self.xc) >= half, 0.0, bump
            )
            return self.r0 - (self.r0 - self.r_throat) * bump
        if abs(x - self.xc) >= half:
            return self.r0
        # cosine bump: 0 at edges, 1 at the throat
        bump = 0.5 * (1.0 + math.cos(math.pi * (x - self.xc) / half))
        return self.r0 - (self.r0 - self.r_throat) * bump

    def inlet_extension(self) -> float:
        """Length of the upstream extension section (m)."""
        return self.xc - self.stenosis_length / 2.0

    def outlet_extension(self) -> float:
        """Length of the downstream recovery section (m)."""
        return self.length - (self.xc + self.stenosis_length / 2.0)


@dataclass(frozen=True)
class FlowWaveform:
    """Pulsatile flow waveform (coronary-like, normalized Fourier form).

    Q(t) = Q_mean * [1 + sum_k c_k * cos(2 pi k t / T + phi_k)]

    The default coefficients give a coronary-shaped cycle: low diastolic
    baseline with a systolic dip (coronary flow is diastolic-dominant).
    Attributes:
        q_mean: Mean flow rate (m^3/s).
        heart_rate: Beats per minute.
        coeffs: (amplitude, phase_rad) per harmonic k = 1..n.
    """

    q_mean: float = 1.2e-6  # ~1.2 mL/s for a coronary artery
    heart_rate: float = 75.0
    coeffs: tuple[tuple[float, float], ...] = (
        (0.55, math.pi),        # systolic reduction (coronary is diastolic)
        (0.25, 2.0 * math.pi * 0.35),
        (0.12, 2.0 * math.pi * 0.15),
    )

    @property
    def period(self) -> float:
        """Cardiac cycle duration T (s)."""
        return 60.0 / self.heart_rate

    @property
    def omega(self) -> float:
        """Fundamental angular frequency (rad/s)."""
        return 2.0 * math.pi / self.period

    def flow(self, t: float) -> float:
        """Instantaneous flow rate (m^3/s) at time t (s)."""
        q = 1.0
        for k, (amp, phi) in enumerate(self.coeffs, start=1):
            q += amp * math.cos(k * self.omega * t + phi)
        return self.q_mean * q

    def mean_flow(self, n_samples: int = 400) -> float:
        """Numerically integrated mean over one cycle (m^3/s)."""
        n = max(2, n_samples)
        dt = self.period / n
        total = 0.0
        for i in range(n):
            total += self.flow(i * dt) * dt
        return total / self.period


@dataclass(frozen=True)
class StentParams:
    """Parametric lattice stent.

    The lattice consists of sinusoidal struts arranged in rings along the
    vessel axis plus axial connectors between adjacent rings. Struts are
    partly embedded in the wall (embed_frac of the strut thickness) to
    avoid sub-millimetre gaps that break CFD meshes.

    Attributes:
        strut_thickness: Strut thickness (m), typical 80-120 um.
        strut_width: Strut width in the axial direction (m).
        crowns: Crowns (strut half-waves) per ring.
        n_rings: Number of strut rings.
        ring_pitch: Axial spacing between rings (m).
        embed_frac: Fraction of strut thickness embedded in the wall.
        cover_frac: Axial vessel length covered by the stent (fraction).
    """

    strut_thickness: float = 100.0e-6
    strut_width: float = 120.0e-6
    crowns: int = 7
    n_rings: int = 6
    ring_pitch: float = 4.0e-3
    embed_frac: float = 0.35
    cover_frac: float = 0.28

    def strut_count(self) -> int:
        """Total number of sinusoidal strut bars."""
        return self.n_rings * self.crowns
