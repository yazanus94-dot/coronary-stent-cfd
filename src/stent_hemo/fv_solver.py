"""Axisymmetric finite-volume SIMPLE solver for stenosed-artery flow.

Solves steady or transient incompressible Navier-Stokes in (z, r) polar
coordinates on the structured grid (uniform axial spacing; radial cells
uniform in the local-radius fraction). STAGGERED arrangement (classic
Patankar SIMPLE): axial velocity u on z-faces, radial velocity v on
r-faces, pressure at cell centers - checkerboard-free and robust.

First-order upwind convection, central diffusion. Newtonian or Carreau
(outer viscosity iteration) rheology. Boundary conditions: developed
(Poiseuille or Womersley-harmonic) inlet profile, no-slip wall, axis
symmetry, fixed-pressure outlet.

This solver provides the real, locally-runnable CFD numbers of the
study: Poiseuille and Womersley validation, stenosis pressure drops,
and the grid-convergence study. Full 3-D pulsatile stented flow runs in
OpenFOAM (Docker/CI) via the generated cases.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from .core import BloodProperties, FlowWaveform, VesselGeometry
from .struct_grid import StructGrid
from .womersley import harmonic_pressure_gradient


@dataclass
class FlowSolution:
    """Result of one solver run (cell-centered fields)."""

    u: np.ndarray  # axial velocity (nz, nr) (m/s)
    v: np.ndarray  # radial velocity (nz, nr) (m/s)
    p: np.ndarray  # pressure (nz, nr) (Pa)
    mu: np.ndarray  # apparent viscosity field (nz, nr) (Pa s)
    iterations: int
    mass_residual: float
    dp: float  # area-averaged static pressure drop inlet-to-outlet (Pa)
    grid: StructGrid
    u_history: np.ndarray | None = None  # transient: (n_samples, nz, nr)
    continuity: dict[str, float] | None = None  # continuity diagnostics (m^3/s, 1/s)


def harmonic_inlet_profile(
    vessel: VesselGeometry,
    blood: BloodProperties,
    wf: FlowWaveform,
    r_frac: np.ndarray,
    t: float,
) -> np.ndarray:
    """Womersley inlet velocity (m/s) at radius fractions and time t."""
    grads = harmonic_pressure_gradient(wf, blood, vessel)
    from scipy.special import jv

    rr = np.asarray(r_frac, dtype=float)
    u = np.zeros_like(rr)
    for k, amp, ph in grads:
        if k == 0:
            g0 = 8.0 * blood.mu_newtonian * wf.mean_flow() / (math.pi * vessel.r0**4)
            u += g0 * vessel.r0**2 * (1.0 - rr**2) / (4.0 * blood.mu_newtonian)
            continue
        alpha_k = vessel.r0 * math.sqrt(k * wf.omega * blood.rho / blood.mu_newtonian)
        lam = alpha_k * (1 + 1j) ** 1.5
        shape = 1.0 - jv(0, lam * rr.astype(complex)) / jv(0, lam)
        u += np.real(
            amp
            / (1j * blood.rho * k * wf.omega)
            * shape
            * np.exp(1j * (k * wf.omega * t + ph))
        )
    return u


class AxisymSolver:
    """Staggered SIMPLE solver on a StructGrid (axisymmetric)."""

    def __init__(
        self,
        grid: StructGrid,
        blood: BloodProperties,
        use_carreau: bool = False,
    ):
        self.grid = grid
        self.blood = blood
        self.use_carreau = use_carreau
        self.nz, self.nr = grid.shape
        self.dz = grid.vessel.length / self.nz
        self._geom()

    # ------------------------------------------------------------------
    def _geom(self) -> None:
        g = self.grid
        self.rc = g.r_center  # (nz,) column radii at cell centers
        self.rf = g.r_face  # (nz+1,) radii at axial faces
        f = g.r_frac  # (nr+1,)
        self.f = f
        # cell-centered p: ring areas and volumes (per radian)
        self.az = g.axial_face_area  # (nz+1, nr) rings at z faces
        self.ar = g.radial_face_area()  # (nz, nr+1) cylindrical faces
        self.vol = g.cell_volume  # (nz, nr)
        # u-CV (face j between columns j-1, j): mean column radius
        self.ru = 0.5 * (self.rc[:-1] + self.rc[1:])  # (nz-1,)
        self.vu = self.dz * 0.5 * (
            (self.ru[:, None] * f[None, 1:]) ** 2
            - (self.ru[:, None] * f[None, :-1]) ** 2
        )  # (nz-1, nr) u-CV volumes
        self.aru_in = self.dz * self.ru[:, None] * f[None, :-1]  # inner cyl area
        self.aru_out = self.dz * self.ru[:, None] * f[None, 1:]
        # v-CV (face k between p-cells k-1, k): radial extent between the
        # HALF-INTEGER ring centers (p-cell ring midpoints)
        fi = 0.5 * (f[None, :-2] + f[None, 1:-1])  # (f_{k-1}+f_k)/2
        fo = 0.5 * (f[None, 1:-1] + f[None, 2:])  # (f_k+f_{k+1})/2
        rc_f = self.rc[:, None]
        self.vv = self.dz * 0.5 * ((rc_f * fo) ** 2 - (rc_f * fi) ** 2)  # (nz, nr-1)
        self.arv_in = self.dz * rc_f * fi  # (nz, nr-1) inner radial face
        self.arv_out = self.dz * rc_f * fo
        # axial annulus areas for the v-CV at z faces j and j+1
        self.azv_lo = 0.5 * self.rf[:-1, None] ** 2 * (fo**2 - fi**2)
        self.azv_hi = 0.5 * self.rf[1:, None] ** 2 * (fo**2 - fi**2)
        self.dru = self.ru[:, None] / self.nr * np.ones((1, self.nr))
        self.drv = self.rc[:, None] / self.nr * np.ones((1, self.nr - 1))

    # ------------------------------------------------------------------
    def _face_mu(self, mu: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Viscosity at u faces (nz+1, nr) and v faces (nz, nr+1)."""
        mu_u = np.empty((self.nz + 1, self.nr))
        mu_u[1:-1] = 0.5 * (mu[:-1] + mu[1:])
        mu_u[0] = mu[0]
        mu_u[-1] = mu[-1]
        mu_v = np.empty((self.nz, self.nr + 1))
        mu_v[:, 1:-1] = 0.5 * (mu[:, :-1] + mu[:, 1:])
        mu_v[:, 0] = mu[:, 0]
        mu_v[:, -1] = mu[:, -1]
        return mu_u, mu_v

    # ------------------------------------------------------------------
    @staticmethod
    def _pair(
        rows, cols, vals, diag, b,
        ip: int, io: int | None, f: float, df: float, ub: float | None,
    ) -> None:
        """Upwind convection + diffusion assembly for one face (own side).

        Adds ONLY the calling cell's (ip's) equation terms: each interior
        face is visited once from each side, so the neighbour contributes
        its own terms when its loop visits the same face.

        f: mass flux from ip toward io (kg/s); df: diffusion coeff;
        ub: face value for a Dirichlet boundary (io is None), or None.
        """
        diag[ip] += df + max(f, 0.0)
        if io is not None:
            rows.append(ip)
            cols.append(io)
            vals.append(-df + min(f, 0.0))
        else:
            # boundary: inflow (f<0) brings ub; diffusion sees ub too
            b[ip] += (max(-f, 0.0) + df) * (ub if ub is not None else 0.0)

    # ------------------------------------------------------------------
    def _u_momentum(
        self,
        u: np.ndarray,
        v: np.ndarray,
        p: np.ndarray,
        mu: np.ndarray,
        dt: float | None,
        u_old: np.ndarray | None,
        relax: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solve u-momentum for interior z-faces; returns (u_new, a_P)."""
        nz, nr = self.nz, self.nr
        nu = nz - 1  # unknown faces j = 1..nz-1
        n = nu * nr
        mu_u, mu_v = self._face_mu(mu)
        rows, cols, vals = [], [], []
        diag = np.zeros(n)
        b = np.zeros(n)
        pair = self._pair

        def sub(j: int, k: int) -> int:
            return (j - 1) * nr + k

        for j in range(1, nz):
            for k in range(nr):
                ip = sub(j, k)
                # axial east face (io = j+1 or outlet boundary)
                if j <= nz - 2:
                    uf = 0.5 * (u[j, k] + u[j + 1, k])
                    pair(rows, cols, vals, diag, b, ip, sub(j + 1, k),
                         self.blood.rho * uf * self.az[j + 1][k],
                         mu_u[j + 1][k] * self.az[j + 1][k] / self.dz, None)
                else:  # outlet: zero-gradient outflow, convection only
                    uf = u[j, k]
                    ff = self.blood.rho * uf * self.az[nz][k]
                    diag[ip] += max(ff, 0.0)
                # axial west face (io = j-1 or inlet boundary)
                if j >= 2:
                    uf = 0.5 * (u[j - 1, k] + u[j, k])
                    pair(rows, cols, vals, diag, b, ip, sub(j - 1, k),
                         -self.blood.rho * uf * self.az[j][k],
                         mu_u[j - 1][k] * self.az[j][k] / self.dz, None)
                else:  # j == 1: inlet Dirichlet (inflow, f measured ip->io = -z)
                    uf = 0.5 * (u[0, k] + u[1, k])
                    f = -self.blood.rho * uf * self.az[0][k]  # < 0: inflow
                    # face at zc_0, unknown at z=dz: distance dz/2
                    df = 2.0 * mu_u[0][k] * self.az[0][k] / self.dz
                    diag[ip] += df  # upwind value is the boundary -> no diag conv
                    b[ip] += (df + max(-f, 0.0)) * u[0, k]
                # radial faces: inner (io = k-1) and outer (io = k+1)
                if k <= nr - 2:
                    vbar = 0.5 * (v[j - 1, k + 1] + v[j, k + 1])
                    pair(rows, cols, vals, diag, b, ip, sub(j, k + 1),
                         self.blood.rho * vbar * self.aru_out[j - 1][k],
                         0.5 * (mu_v[j - 1, k + 1] + mu_v[j, k + 1])
                         * self.aru_out[j - 1][k] / self.dru[j - 1][k], None)
                else:  # wall (k+1 = nr): v = 0, diffusion Dirichlet u = 0
                    a_r = self.aru_out[j - 1][k]
                    # unknown at ring center, wall at f=1: distance dru/2
                    df = (mu_v[j - 1, nr] + mu_v[j, nr]) * a_r / self.dru[j - 1][k]
                    diag[ip] += df
                if k >= 1:
                    vbar = 0.5 * (v[j - 1, k] + v[j, k])
                    pair(rows, cols, vals, diag, b, ip, sub(j, k - 1),
                         -self.blood.rho * vbar * self.aru_in[j - 1][k],
                         0.5 * (mu_v[j - 1, k] + mu_v[j, k])
                         * self.aru_in[j - 1][k] / self.dru[j - 1][k], None)
                # else k == 0: axis face, zero area -> nothing
                # pressure driving
                b[ip] += -(p[j, k] - p[j - 1, k]) / self.dz * self.vu[j - 1, k]
        # temporal term
        if dt is not None and u_old is not None:
            coeff = self.blood.rho * self.vu / dt
            diag += coeff.ravel()
            b += coeff.ravel() * u_old[1:nz].ravel()
        # under-relaxation (Patankar)
        u_guess = u[1:nz].ravel()
        b = b + (1.0 - relax) / relax * diag * u_guess
        diag = diag / relax
        rows.extend(range(n))
        cols.extend(range(n))
        vals.extend(diag)
        a_mat = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
        sol = spla.spsolve(a_mat, b).reshape(nu, nr)
        unew = u.copy()
        unew[1:nz] = sol
        return unew, diag.reshape(nu, nr)

    # ------------------------------------------------------------------
    def _v_momentum(
        self,
        u: np.ndarray,
        v: np.ndarray,
        p: np.ndarray,
        mu: np.ndarray,
        dt: float | None,
        v_old: np.ndarray | None,
        relax: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        """Solve v-momentum for interior r-faces; returns (v_new, a_P)."""
        nz, nr = self.nz, self.nr
        nv = nr - 1  # unknown faces k = 1..nr-1 (all nz columns)
        n = nz * nv
        mu_u, mu_v = self._face_mu(mu)
        rows, cols, vals = [], [], []
        diag = np.zeros(n)
        b = np.zeros(n)
        pair = self._pair

        def sub(j: int, k: int) -> int:
            return j * nv + (k - 1)

        for j in range(nz):
            for k in range(1, nr):
                ip = sub(j, k)
                # radial outer face (io = k+1 or wall)
                if k <= nr - 2:
                    vf = 0.5 * (v[j, k] + v[j, k + 1])
                    pair(rows, cols, vals, diag, b, ip, sub(j, k + 1),
                         self.blood.rho * vf * self.arv_out[j][k - 1],
                         mu_v[j][k + 1] * self.arv_out[j][k - 1]
                         / self.drv[j][k - 1], None)
                else:  # wall (k+1 = nr): v = 0, diffusion Dirichlet
                    a_r = self.arv_out[j][k - 1]
                    df = mu_v[j][nr] * a_r / self.drv[j][k - 1]
                    diag[ip] += df
                # radial inner face (io = k-1 or axis)
                if k >= 2:
                    vf = 0.5 * (v[j, k - 1] + v[j, k])
                    pair(rows, cols, vals, diag, b, ip, sub(j, k - 1),
                         -self.blood.rho * vf * self.arv_in[j][k - 1],
                         mu_v[j][k - 1] * self.arv_in[j][k - 1]
                         / self.drv[j][k - 1], None)
                # k == 1: axis face has zero area -> nothing
                # axial east face (io = j+1 or outlet boundary)
                if j <= nz - 2:
                    uf = 0.5 * (u[j + 1, k - 1] + u[j + 1, k])
                    pair(rows, cols, vals, diag, b, ip, sub(j + 1, k),
                         self.blood.rho * uf * self.azv_hi[j][k - 1],
                         0.5 * (mu_u[j + 1, k - 1] + mu_u[j + 1, k])
                         * self.azv_hi[j][k - 1] / self.dz, None)
                else:  # outlet: zero-gradient outflow
                    uf = 0.5 * (u[nz, k - 1] + u[nz, k])
                    ff = self.blood.rho * uf * self.azv_hi[j][k - 1]
                    diag[ip] += max(ff, 0.0)
                # axial west face (io = j-1 or inlet boundary)
                if j >= 1:
                    uf = 0.5 * (u[j, k - 1] + u[j, k])
                    pair(rows, cols, vals, diag, b, ip, sub(j - 1, k),
                         -self.blood.rho * uf * self.azv_lo[j][k - 1],
                         0.5 * (mu_u[j, k - 1] + mu_u[j, k])
                         * self.azv_lo[j][k - 1] / self.dz, None)
                else:  # j == 0: inlet, convected v = 0 (Dirichlet)
                    # face at z=0, unknown at zc_0: distance dz/2
                    df = (mu_u[0, k - 1] + mu_u[0, k]) * self.azv_lo[0][k - 1] / self.dz
                    diag[ip] += df  # upwind value is the boundary -> no diag conv
                # pressure driving
                b[ip] += (
                    -(p[j, k] - p[j, k - 1]) / self.drv[j][k - 1] * self.vv[j][k - 1]
                )
        if dt is not None and v_old is not None:
            coeff = self.blood.rho * self.vv / dt
            diag += coeff.ravel()
            b += coeff.ravel() * v_old[:, 1:nr].ravel()
        v_guess = v[:, 1:nr].ravel()
        b = b + (1.0 - relax) / relax * diag * v_guess
        diag = diag / relax
        rows.extend(range(n))
        cols.extend(range(n))
        vals.extend(diag)
        a_mat = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
        sol = spla.spsolve(a_mat, b).reshape(nz, nv)
        vnew = v.copy()
        vnew[:, 1:nr] = sol
        return vnew, diag.reshape(nz, nv)

    # ------------------------------------------------------------------
    def _poisson(
        self,
        u: np.ndarray,
        v: np.ndarray,
        du: np.ndarray,
        dv: np.ndarray,
    ) -> np.ndarray:
        """Pressure-correction Poisson from star-face fluxes.

        du: a_P-inverse for u (nu, nr); dv: for v (nz, nv).
        """
        nz, nr, nv = self.nz, self.nr, self.nr - 1
        n = nz * nr
        rows, cols, vals = [], [], []
        diag = np.zeros(n)
        b = np.zeros(n)
        rho = self.blood.rho

        def pidx(j: int, k: int) -> int:
            return j * nr + k

        def uidx(j: int, k: int) -> int:
            return (j - 1) * nr + k

        def vidx(j: int, k: int) -> int:
            return j * nv + (k - 1)

        # axial u faces j = 1..nz-1 (interior)
        for j in range(1, nz):
            for k in range(nr):
                ip, io = pidx(j - 1, k), pidx(j, k)
                a_z = self.az[j][k]
                d_f = rho * a_z * du[j - 1, k] / self.dz
                diag[ip] += d_f
                diag[io] += d_f
                rows.extend((ip, io))
                cols.extend((io, ip))
                vals.extend((-d_f, -d_f))
                flux = rho * u[j, k] * a_z
                b[ip] -= flux  # outflow of ip
                b[io] += flux  # inflow of io
        # inlet faces j = 0 (fixed inflow)
        for k in range(nr):
            ip = pidx(0, k)
            b[ip] += rho * u[0, k] * self.az[0][k]
        # outlet faces j = nz (p' Dirichlet 0)
        for k in range(nr):
            ip = pidx(nz - 1, k)
            a_z = self.az[nz][k]
            d_f = rho * a_z * du[nz - 2, k] / self.dz if nz > 2 else rho * a_z * du[0, k] / self.dz
            diag[ip] += d_f
            b[ip] -= rho * u[nz, k] * a_z
        # radial v faces k = 1..nr-1 (interior)
        for j in range(nz):
            for k in range(1, nr):
                ip, io = pidx(j, k - 1), pidx(j, k)
                a_r = self.ar[j][k]
                d_f = rho * a_r * dv[j, k - 1] / self.drv[j][k - 1]
                diag[ip] += d_f
                diag[io] += d_f
                rows.extend((ip, io))
                cols.extend((io, ip))
                vals.extend((-d_f, -d_f))
                flux = rho * v[j, k] * a_r
                b[ip] -= flux
                b[io] += flux
        rows.extend(range(n))
        cols.extend(range(n))
        vals.extend(diag)
        a_mat = sp.csc_matrix((vals, (rows, cols)), shape=(n, n))
        pp = spla.spsolve(a_mat, b)
        return pp.reshape(nz, nr)

    # ------------------------------------------------------------------
    def _carreau_field(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Apparent viscosity (nz, nr) from the local strain rate.

        ``u`` is the axial velocity on the staggered z-faces (nz+1, nr)
        and ``v`` the radial velocity on the staggered r-faces
        (nz, nr+1). Interior face values are interpolated to cell centers
        (face average); the shear-rate magnitude follows the axisymmetric
        no-swirl tensor magnitude

            gamma^2 = dudz^2 + dvdr^2 + 2 (v_r / r)^2
                      + 0.5 (dudr + dvdz)^2

        where dudz = du_z/dz, dvdr = du_r/dr, dvdz = du_r/dz and dudr =
        du_z/dr, all evaluated on the cell centers with the cylindrical
        spacing (rc/nr) between adjacent cell centers.
        """
        r_in, r_out = self.grid.cell_r_bounds()
        r_c = 0.5 * (r_in + r_out)  # (nz, nr) cell-center mean radii
        u_c = 0.5 * (u[:-1] + u[1:])  # (nz, nr) cell-centered axial velocity
        v_c = 0.5 * (v[:, :-1] + v[:, 1:])  # (nz, nr) cell-centered radial velocity
        dz = self.dz
        dr = self.rc / self.nr  # (nz,) radial spacing, uniform per column

        def rad_grad(field: np.ndarray) -> np.ndarray:
            """Cylindrical gradient of a cell-centered field (1/s)."""
            out = np.empty_like(field)
            out[:, 1:-1] = (field[:, 2:] - field[:, :-2]) / (2.0 * dr[:, None])
            out[:, 0] = (field[:, 1] - field[:, 0]) / dr
            out[:, -1] = (field[:, -1] - field[:, -2]) / dr
            return out

        dudz = np.gradient(u_c, dz, axis=0)  # du_z/dz
        dvdz = np.gradient(v_c, dz, axis=0)  # du_r/dz
        dvdr = rad_grad(v_c)  # du_r/dr
        dudr = rad_grad(u_c)  # du_z/dr
        vt = v_c / np.maximum(r_c, 1e-12)  # v_r / r at the cell center
        rate = np.sqrt(
            dudz**2 + dvdr**2 + 2.0 * vt**2 + 0.5 * (dudr + dvdz) ** 2
        )
        vec_carreau = np.vectorize(self.blood.carreau_viscosity)
        return vec_carreau(rate)

    # ------------------------------------------------------------------
    def solve(
        self,
        wf: FlowWaveform,
        steady: bool = True,
        dt: float = 2.0e-3,
        n_cycles: int = 3,
        max_outer: int = 400,
        tol: float = 1e-4,
        relax: float = 0.7,
        verbose: bool = False,
    ) -> FlowSolution:
        """Run the solver (steady Poiseuille or transient Womersley inlet)."""
        nz, nr = self.nz, self.nr
        vessel: VesselGeometry = self.grid.vessel
        q_scale = wf.mean_flow()
        r_in = self.rf[0]
        f0 = self.f
        # ring-FLUX-AVERAGED parabolic inlet: u_ring = U (1 - (f_k^2+f_{k+1}^2)/2)
        # makes sum(u[0]*az[0]) = q exactly
        u = np.zeros((nz + 1, nr))
        u[0] = 2.0 * q_scale / (math.pi * r_in**2) * (
            1.0 - 0.5 * (f0[None, :-1] ** 2 + f0[None, 1:] ** 2)
        )
        vel = np.zeros((nz, nr + 1))
        p = np.zeros((nz, nr))
        mu = np.full((nz, nr), self.blood.mu_newtonian)
        u_old = u.copy()
        v_old = vel.copy()

        if steady:
            u[1:] = u[0][None, :]
            it_final, res_final = 0, float("inf")
            for outer in range(1, max_outer + 1):
                if self.use_carreau:
                    mu = self._carreau_field(u, vel)
                u, a_u = self._u_momentum(u, vel, p, mu, None, None, relax)
                v, a_v = self._v_momentum(u, vel, p, mu, None, None, relax)
                u[nz] = u[nz - 1]  # outlet extrapolation
                d_u = self.vu / np.maximum(a_u, 1e-300)  # V/a_P (nz-1, nr)
                d_v = self.vv / np.maximum(a_v, 1e-300)  # (nz, nr-1)
                # Poisson star-face fluxes from the PREDICTED fields: the
                # predicted axial velocity (from _u_momentum) and the
                # predicted radial velocity (from _v_momentum) are what the
                # pressure correction integrates; the correction below
                # updates the predicted field itself, so the radial
                # momentum predictor is part of the SIMPLE sequence.
                pp = self._poisson(u, v, d_u, d_v)
                rp = 0.5  # pressure under-relaxation
                p = p + rp * pp
                # apply corrections from the CORRECTION field pp only
                for j in range(1, nz):
                    u[j] -= d_u[j - 1] * rp * (pp[j] - pp[j - 1]) / self.dz
                for j in range(nz):
                    for k in range(1, nr):
                        v[j, k] -= d_v[j, k - 1] * rp * (
                            pp[j, k] - pp[j, k - 1]
                        ) / self.drv[j][k - 1]
                v[:, 0] = 0.0
                v[:, nr] = 0.0
                # the corrected predicted field becomes the working field
                vel = v
                u[nz] = u[nz - 1]
                res = self._mass_res(u, vel)
                it_final, res_final = outer, res
                if verbose and outer % 25 == 0:
                    print(f"  iter {outer}: res {res:.3e}")
                if res < tol * abs(q_scale):
                    break
            uc = 0.5 * (u[:-1] + u[1:])
            dp = float(np.average(p[0], weights=self.az[0]) - np.average(p[-1], weights=self.az[-1]))
            continuity = self._continuity_stats(u, vel)
            return FlowSolution(uc, vel[:, :-1], p, mu, it_final, res_final, dp, self.grid,
                                continuity=continuity)

        # ---- transient (Womersley pulsatile inlet) -----------------------
        n_per = int(round(wf.period / dt))
        n_steps = n_per * n_cycles
        # store the last TWO consecutive cycles so periodicity can be
        # validated between identical phases of consecutive cycles
        keep = 2 * n_per if n_cycles >= 2 else n_per
        u_hist = []
        it_final, res_final = 0, 0.0
        for step in range(1, n_steps + 1):
            t = step * dt
            u_old, v_old = u.copy(), vel.copy()
            u[0] = harmonic_inlet_profile(
                vessel, self.blood, wf, 0.5 * (f0[:-1] + f0[1:]), t % wf.period
            )
            u[nz] = u[0]
            mu_loc = mu
            for _ in range(3):
                if self.use_carreau:
                    mu_loc = self._carreau_field(u, vel)
                u, a_u = self._u_momentum(u, vel, p, mu_loc, dt, u_old, 0.8)
                v, a_v = self._v_momentum(u, vel, p, mu_loc, dt, v_old, 0.8)
                u[nz] = u[nz - 1]
                d_u = self.vu / np.maximum(a_u, 1e-300)
                d_v = self.vv / np.maximum(a_v, 1e-300)
                # Poisson star-face fluxes from the PREDICTED fields (same
                # integration of the radial momentum predictor as steady).
                pp = self._poisson(u, v, d_u, d_v)
                rp = 0.8  # pressure under-relaxation
                p = p + rp * pp
                for j in range(1, nz):
                    u[j] -= d_u[j - 1] * rp * (pp[j] - pp[j - 1]) / self.dz
                for j in range(nz):
                    for k in range(1, nr):
                        v[j, k] -= d_v[j, k - 1] * rp * (
                            pp[j, k] - pp[j, k - 1]
                        ) / self.drv[j][k - 1]
                v[:, 0] = 0.0
                v[:, nr] = 0.0
                u[nz] = u[nz - 1]
                vel = v
            res = self._mass_res(u, v)
            it_final, res_final = step, res
            if step > n_steps - keep:
                u_hist.append(0.5 * (u[:-1] + u[1:]).copy())
        uc = 0.5 * (u[:-1] + u[1:])
        dp = float(np.average(p[0], weights=self.az[0]) - np.average(p[-1], weights=self.az[-1]))
        return FlowSolution(
            uc, v[:, :-1], p, mu, it_final, res_final, dp, self.grid,
            u_history=np.array(u_hist) if u_hist else None,
            continuity=self._continuity_stats(u, v),
        )

    # ------------------------------------------------------------------
    def _mass_res(self, u: np.ndarray, v: np.ndarray) -> float:
        """Net volumetric outflow imbalance over the whole domain (m^3/s).

        Sum of the per-cell FV fluxes through the staggered faces: the
        axial terms telescope to ``u[nz]*az[nz] - u[0]*az[0]``; the radial
        terms telescope to the zero wall/axis fluxes. The density is NOT
        multiplied in -- the residual is a volumetric flow imbalance in
        SI units, matching the ``tol * |q|`` convergence tolerance. The
        magnitude is returned (non-negative), as in the original
        implementation, so the ``res < tol*|q|`` check stays meaningful.
        """
        return abs(float(np.sum(self._cell_continuity(u, v))))

    # ------------------------------------------------------------------
    def _cell_continuity(self, u: np.ndarray, v: np.ndarray) -> np.ndarray:
        """Per-cell net volumetric outflow (m^3/s) of the FV face fluxes.

        ``u`` is the axial velocity on z-faces (nz+1, nr); ``v`` the radial
        velocity on r-faces (nz, nr+1); cells are (nz, nr).
        """
        dz_flux = u[1:, :] * self.az[1:, :] - u[:-1, :] * self.az[:-1, :]
        r_flux = v[:, 1:] * self.ar[:, 1:] - v[:, :-1] * self.ar[:, :-1]
        return dz_flux + r_flux

    # ------------------------------------------------------------------
    def _continuity_stats(self, u: np.ndarray, v: np.ndarray) -> dict[str, float]:
        """Continuity diagnostics on the corrected face fields.

        ``net_outflow`` is the domain-level imbalance (m^3/s, the same
        quantity as ``mass_residual``); ``cell_max_abs`` and ``cell_rms``
        are the normalized per-cell divergence residuals (1/s) of the FV
        fluxes divided by the cell volumes; ``cell_rel_max`` is the
        dimensionless per-cell continuity residual: the magnitude of the
        unbalanced flux as a fraction of the cell's gross face fluxes.
        """
        flux = self._cell_continuity(u, v)
        divergence = flux / np.maximum(self.vol, 1e-300)
        # per-face flux magnitudes: out-faces and in-faces of each cell
        dz_out = np.abs(u[1:, :] * self.az[1:, :])
        dz_in = np.abs(u[:-1, :] * self.az[:-1, :])
        r_out = np.abs(v[:, 1:] * self.ar[:, 1:])
        r_in = np.abs(v[:, :-1] * self.ar[:, :-1])
        gross = dz_out + dz_in + r_out + r_in
        cell_rel = np.abs(flux) / np.maximum(gross, 1e-300)
        return {
            "net_outflow": float(np.sum(flux)),
            "cell_max_abs": float(np.max(np.abs(divergence))),
            "cell_rms": float(np.sqrt(np.mean(divergence**2))),
            "cell_rel_max": float(np.max(cell_rel)),
        }
