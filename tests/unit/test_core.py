"""Unit tests: geometry, analytic helpers, mesh writers."""

from __future__ import annotations

import math

import numpy as np
import pytest

from stent_hemo.core import StentParams, VesselGeometry
from stent_hemo.fv_solver import harmonic_inlet_profile
from stent_hemo.hand_calc import flow_hand_calc
from stent_hemo.mesh_io import write_gmsh, write_su2
from stent_hemo.struct_grid import StructGrid
from stent_hemo.young_tsai import stenosis_pressure_drop

# ---- geometry ---------------------------------------------------------


def test_radius_profile_bounds(vessel_straight, vessel_stenosis_50, vessel_stenosis_70):
    """Lumen radius stays positive and within nominal bounds."""
    for v in (vessel_straight, vessel_stenosis_50, vessel_stenosis_70):
        z = np.linspace(0.0, v.length, 200)
        r = v.radius_at(z)
        assert r.min() > 0.0
        assert r.max() <= v.r0 + 1e-12
        if v.severity > 0.0:
            # smooth cosine bump narrows but does not reach nominal throat
            assert v.r0 * (1.0 - v.severity) <= r.min() < v.r0 - 1e-9


# ---- analytic helpers -------------------------------------------------


def test_hand_calc_poiseuille_dp_matches_formula(vessel_straight, wf, blood):
    """Hand-calculated Poiseuille dp matches the analytic formula."""
    hc = flow_hand_calc(vessel_straight, blood, wf)
    exact = (
        8
        * blood.mu_newtonian
        * wf.mean_flow()
        * vessel_straight.length
        / (math.pi * vessel_straight.r0**4)
    )
    assert hc.dp_vessel == pytest.approx(exact, rel=1e-9)


def test_reynolds_order_of_magnitude(vessel_straight, wf, blood):
    """Reynolds number stays in the laminar coronary regime."""
    hc = flow_hand_calc(vessel_straight, blood, wf)
    assert 50 < hc.reynolds < 500


def test_young_tsai_increases_with_severity(
    vessel_straight, vessel_stenosis_50, vessel_stenosis_70, wf, blood
):
    """Young-Tsai pressure drop increases monotonically with severity."""
    q = wf.mean_flow()
    drops = [
        stenosis_pressure_drop(v, blood, q).dp_total
        for v in (vessel_straight, vessel_stenosis_50, vessel_stenosis_70)
    ]
    assert drops[0] < drops[1] < drops[2]
    assert drops[0] < 10.0


# ---- GCI --------------------------------------------------------------


def test_gci_on_exact_quadratic():
    """GCI recovers p=2 for f(h) = f_exact + c*h^2."""
    from stent_hemo.gci import grid_convergence_index

    f_exact, c = 10.0, 3.0
    res = grid_convergence_index(f_exact + c / 16, f_exact + c / 4, f_exact + c, 2.0)
    assert res.p == pytest.approx(2.0, rel=1e-6)
    assert res.f_exact == pytest.approx(f_exact, rel=1e-6)


# ---- mesh writers -----------------------------------------------------


@pytest.mark.parametrize("writer", [write_gmsh, write_su2])
def test_mesh_writers_smoke(tmp_path, writer):
    """Gmsh and SU2 mesh writers produce non-empty output with markers."""
    grid = StructGrid(VesselGeometry(severity=0.5, length=30e-3), 8, 4)
    ext = "msh" if writer is write_gmsh else "su2"
    path = writer(grid, tmp_path / f"out.{ext}")
    text = path.read_text()
    assert text.strip()
    if writer is write_su2:
        assert "NDIME= 2" in text
        for marker in ("INLET", "OUTLET", "WALL"):
            assert marker in text


# ---- solver basics -----------------------------------------------------


def test_steady_poiseuille_small_grid(grid_coarse, wf, blood):
    """Coarse straight-tube run: dp within tolerance of analytic, no NaN."""
    solver = __import__("stent_hemo.fv_solver", fromlist=["AxisymSolver"]).AxisymSolver(
        grid_coarse, blood
    )
    sol = solver.solve(wf, steady=True, max_outer=400, tol=1e-6)
    u = np.asarray(sol.u)
    assert np.isfinite(u).all()
    assert np.isfinite(sol.p).all()
    dp_exact = (
        8
        * blood.mu_newtonian
        * wf.mean_flow()
        * grid_coarse.vessel.length
        / (math.pi * grid_coarse.vessel.r0**4)
    )
    assert abs(sol.dp - dp_exact) / dp_exact < 0.15
    assert sol.u.min() > -1e-9


def test_transient_inlet_harmonic_periodicity(grid_coarse, wf, blood):
    """Inlet profile at t and t+T must match (shape, values)."""

    f = grid_coarse.r_frac
    rmid = 0.5 * (f[:-1] + f[1:])
    u0 = harmonic_inlet_profile(grid_coarse.vessel, blood, wf, rmid, 0.0)
    u1 = harmonic_inlet_profile(
        grid_coarse.vessel, blood, wf, rmid, wf.period
    )
    np.testing.assert_allclose(u0, u1, rtol=1e-9)


def test_stent_params_consistency():
    """Stent parameters are internally consistent."""
    sp = StentParams()
    assert 0.0 < sp.embed_frac < 1.0
    assert 0.0 < sp.cover_frac <= 1.0
    assert sp.crowns > 0 and sp.n_rings > 0
