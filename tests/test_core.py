"""Unit tests: geometry, analytic helpers, mesh writers, solver basics.

The heavy V&V runs live in stent_hemo.validation and are exercised on
small grids here; the full suite is run by `stent-hemo validate` / CI.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from stent_hemo.core import (
    BloodProperties,
    FlowWaveform,
    StentParams,
    VesselGeometry,
)
from stent_hemo.fv_solver import AxisymSolver, harmonic_inlet_profile
from stent_hemo.gci import grid_convergence_index
from stent_hemo.hand_calc import flow_hand_calc
from stent_hemo.mesh_io import write_gmsh, write_su2
from stent_hemo.struct_grid import StructGrid
from stent_hemo.young_tsai import stenosis_pressure_drop


# ---- geometry ---------------------------------------------------------
def test_radius_profile_bounds():
    for sev in (0.0, 0.5, 0.7):
        v = VesselGeometry(severity=sev)
        z = np.linspace(0.0, v.length, 200)
        r = v.radius_at(z)
        assert r.min() > 0.0
        assert r.max() <= v.r0 + 1e-12
        if sev > 0.0:
            # smooth cosine bump narrows but does not reach the nominal
            # throat radius exactly (documented in design-decisions)
            assert v.r0 * (1.0 - sev) <= r.min() < v.r0 - 1e-9


# ---- analytic helpers -------------------------------------------------
def test_hand_calc_poiseuille_dp_matches_formula():
    blood = BloodProperties()
    vessel = VesselGeometry(severity=0.0, length=30e-3)
    wf = FlowWaveform()
    hc = flow_hand_calc(vessel, blood, wf)
    exact = 8 * blood.mu_newtonian * wf.mean_flow() * vessel.length / (
        math.pi * vessel.r0**4
    )
    assert hc.dp_vessel == pytest.approx(exact, rel=1e-9)


def test_reynolds_order_of_magnitude():
    blood = BloodProperties()
    vessel = VesselGeometry(severity=0.0, length=90e-3)
    hc = flow_hand_calc(vessel, blood, FlowWaveform())
    assert 50 < hc.reynolds < 500  # laminar coronary regime


def test_young_tsai_increases_with_severity():
    blood = BloodProperties()
    q = FlowWaveform().mean_flow()
    drops = [
        stenosis_pressure_drop(VesselGeometry(severity=s, length=90e-3), blood, q).dp_total
        for s in (0.0, 0.5, 0.7)
    ]
    assert drops[0] < drops[1] < drops[2]
    assert drops[0] < 10.0  # unstenosed throat -> negligible excess drop


# ---- GCI --------------------------------------------------------------
def test_gci_on_exact_quadratic():
    # f(h) = f_exact + c*h^2 with h = 1, 1/2, 1/4 -> p == 2
    f_exact, c = 10.0, 3.0
    res = grid_convergence_index(
        f_exact + c / 16, f_exact + c / 4, f_exact + c, 2.0
    )
    assert res.p == pytest.approx(2.0, rel=1e-6)
    assert res.f_exact == pytest.approx(f_exact, rel=1e-6)


# ---- mesh writers -----------------------------------------------------
@pytest.mark.parametrize("writer", [write_gmsh, write_su2])
def test_mesh_writers_smoke(tmp_path, writer):
    grid = StructGrid(VesselGeometry(severity=0.5, length=30e-3), 8, 4)
    path = writer(grid, tmp_path / ("out." + ("msh" if writer is write_gmsh else "su2")))
    text = path.read_text()
    assert text.strip()
    if writer is write_su2:
        assert "NDIME= 2" in text
        for marker in ("INLET", "OUTLET", "WALL"):
            assert marker in text


# ---- solver -----------------------------------------------------------
def test_steady_poiseuille_small_grid():
    """Coarse straight-tube run: dp within 15% of analytic, no NaN."""
    blood = BloodProperties()
    vessel = VesselGeometry(severity=0.0, length=30e-3)
    grid = StructGrid(vessel, 20, 8)
    solver = AxisymSolver(grid, blood)
    wf = FlowWaveform()
    sol = solver.solve(wf, steady=True, max_outer=400, tol=1e-6)
    u = np.asarray(sol.u)
    assert np.isfinite(u).all()
    assert np.isfinite(sol.p).all()
    dp_exact = 8 * blood.mu_newtonian * wf.mean_flow() * vessel.length / (
        math.pi * vessel.r0**4
    )
    assert abs(sol.dp - dp_exact) / dp_exact < 0.15
    assert sol.u.min() > -1e-9  # no reversed axial flow in a straight tube


def test_transient_inlet_harmonic_periodicity():
    """Inlet profile at t and t+T must match (shape, values)."""
    blood = BloodProperties()
    vessel = VesselGeometry(severity=0.0, length=30e-3)
    grid = StructGrid(vessel, 10, 4)
    wf = FlowWaveform()
    f = grid.r_frac
    rmid = 0.5 * (f[:-1] + f[1:])
    u0 = harmonic_inlet_profile(vessel, blood, wf, rmid, 0.0)
    u1 = harmonic_inlet_profile(vessel, blood, wf, rmid, wf.period)
    np.testing.assert_allclose(u0, u1, rtol=1e-9)


def test_transient_no_divergence_short():
    """Three coarse-cycle transient run must stay finite (regression for
    the Carreau shear-thinning feedback loop blow-up)."""
    blood = BloodProperties()
    vessel = VesselGeometry(severity=0.5, length=90e-3)
    grid = StructGrid(vessel, 20, 6)
    solver = AxisymSolver(grid, blood)
    wf = FlowWaveform()
    sol = solver.solve(wf, steady=False, dt=2e-3, n_cycles=2)
    hist = sol.u_history
    assert hist is not None and np.isfinite(hist).all()


def test_stent_params_consistency():
    sp = StentParams()
    assert 0.0 < sp.embed_frac < 1.0
    assert 0.0 < sp.cover_frac <= 1.0
    assert sp.crowns > 0 and sp.n_rings > 0
