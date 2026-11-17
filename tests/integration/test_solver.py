"""Integration tests: solver end-to-end, pipeline stages, validation suite."""

from __future__ import annotations

import numpy as np

from stent_hemo.fv_solver import AxisymSolver
from stent_hemo.struct_grid import StructGrid
from stent_hemo.validation import run_validation

# ---- solver integration ------------------------------------------------


def test_steady_stenosis_solve_finite(vessel_stenosis_50, wf, blood):
    """50% stenosis steady solve produces finite fields and reasonable dp."""
    grid = StructGrid(vessel_stenosis_50, 60, 12)
    solver = AxisymSolver(grid, blood)
    sol = solver.solve(wf, steady=True, max_outer=600, tol=1e-6)
    assert np.isfinite(sol.u).all()
    assert np.isfinite(sol.p).all()
    assert sol.dp > 0  # pressure drops across stenosis


def test_transient_stays_finite(vessel_stenosis_50, wf, blood):
    """Transient run on stenosis geometry stays finite (Courant-stable)."""
    grid = StructGrid(vessel_stenosis_50, 20, 6)
    solver = AxisymSolver(grid, blood)
    sol = solver.solve(wf, steady=False, dt=2e-3, n_cycles=2)
    hist = sol.u_history
    assert hist is not None
    assert np.isfinite(hist).all()


# ---- validation integration --------------------------------------------


def test_validation_quick_accepts():
    """Quick V&V suite passes acceptance criteria."""
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp:
        results_dir = Path(tmp) / "results"
        out = run_validation(results_dir, quick=True)
    assert out["accept"] is True


def test_validation_full_grid_converges():
    """Full validation grid converges: fine-grid Poiseuille dp error < 5%."""
    import tempfile
    from pathlib import Path

    with tempfile.TemporaryDirectory() as tmp:
        results_dir = Path(tmp) / "results"
        out = run_validation(results_dir, quick=False)
    assert out["accept"] is True
    pois = out["poiseuille"]
    assert pois[-1]["dp_error_pct"] < 5.0
