"""Validation & verification suite.

Runs the in-house solver against analytic/correlation references and writes
``results/validation.json`` + a markdown table fragment:

1. Steady Poiseuille in a straight tube (exact Hagen-Poiseuille).
2. Transient Womersley inlet (periodicity + qualitative check).
3. Stenosis pressure drop vs. the Young & Tsai (1973) correlation.
4. Grid Convergence Index (Roache) on three refinement levels.

Exit code 0 if every acceptance criterion holds (see ``accept``).
"""

from __future__ import annotations

import json
import math
import time
from pathlib import Path
from typing import Any

from .core import BloodProperties, FlowWaveform, VesselGeometry
from .fv_solver import AxisymSolver
from .gci import grid_convergence_index
from .struct_grid import StructGrid
from .young_tsai import stenosis_pressure_drop


def run_validation(
    results_dir: str | Path = "results",
    quick: bool = False,
) -> dict[str, Any]:
    """Run the full V&V suite; returns the results dict and writes JSON/MD."""
    results_dir = Path(results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    blood = BloodProperties()
    wf = FlowWaveform()
    out: dict[str, Any] = {"poiseuille": [], "womersley": {}, "stenosis": [], "gci": {}}
    t_start = time.perf_counter()

    # ---- 1. Poiseuille ------------------------------------------------
    vessel = VesselGeometry(severity=0.0, length=30e-3)
    q = wf.mean_flow()
    mu = blood.mu_newtonian
    r = vessel.r0
    dp_exact = 8 * mu * q * vessel.length / (math.pi * r**4)
    uc_exact = 2 * q / (math.pi * r**2)
    levels = [(20, 8), (40, 14), (80, 28)] if not quick else [(20, 8)]
    for nz, nr in levels:
        grid = StructGrid(vessel, nz, nr)
        solver = AxisymSolver(grid, blood)
        t0 = time.perf_counter()
        sol = solver.solve(wf, steady=True, max_outer=400, tol=1e-6)
        dp_err = 100 * abs(sol.dp - dp_exact) / dp_exact
        uc_err = 100 * abs(float(sol.u[nz // 2].max()) - uc_exact) / uc_exact
        out["poiseuille"].append(
            {
                "grid": f"{nz}x{nr}",
                "dp": sol.dp,
                "dp_error_pct": dp_err,
                "u_center": float(sol.u[nz // 2].max()),
                "u_center_error_pct": uc_err,
                "iterations": sol.iterations,
                "mass_residual": sol.mass_residual,
                "wall_s": time.perf_counter() - t0,
            }
        )

    # ---- 2. Womersley transient ---------------------------------------
    grid = StructGrid(vessel, 40, 14)
    solver = AxisymSolver(grid, blood)
    t0 = time.perf_counter()
    sol = solver.solve(wf, steady=False, dt=2e-3, n_cycles=3)
    hist = sol.u_history
    assert hist is not None and len(hist) >= 2
    n_per = hist.shape[0]
    periodicity = float(np_maxabs(hist[0], hist[-1]))
    out["womersley"] = {
        "grid": "40x14",
        "dt_s": 2e-3,
        "cycles": 3,
        "steps": sol.iterations,
        "u_center_range": [float(hist.min()), float(hist.max())],
        "periodicity_last_cycle_endpoints": periodicity,
        "flow_reversal": bool(hist.min() < 0.0),
        "wall_s": time.perf_counter() - t0,
    }

    # ---- 3. Stenosis dp vs Young-Tsai ---------------------------------
    # The Young-Tsai correlation gives the EXCESS drop across the stenosis
    # section; compare with dp_solver minus the analytic straight-tube
    # baseline for the same full vessel length.
    for sev in (0.5, 0.7):
        vessel = VesselGeometry(severity=sev, length=90e-3)
        dp_straight = 8 * mu * q * vessel.length / (math.pi * vessel.r0**4)
        grid = StructGrid(vessel, 120, 24)
        solver = AxisymSolver(grid, blood)
        t0 = time.perf_counter()
        sol = solver.solve(wf, steady=True, max_outer=600, tol=1e-6)
        vt = stenosis_pressure_drop(vessel, blood, q)
        excess = sol.dp - dp_straight
        out["stenosis"].append(
            {
                "severity": sev,
                "dp_solver": sol.dp,
                "dp_straight_baseline": dp_straight,
                "dp_excess_solver": excess,
                "dp_young_tsai": vt.dp_total,
                "ratio_to_correlation": excess / vt.dp_total,
                "iterations": sol.iterations,
                "wall_s": time.perf_counter() - t0,
            }
        )

    # ---- 4. GCI on the Poiseuille dp ----------------------------------
    if len(out["poiseuille"]) == 3:
        f_coarse, f_med, f_fine = (row["dp"] for row in out["poiseuille"])
        res = grid_convergence_index(f_fine, f_med, f_coarse, 2.0)
        out["gci"] = {
            "quantity": "dp (Pa), straight tube",
            "values": [f_coarse, f_med, f_fine],
            "r": 2.0,
            "p": res.p,
            "f_exact": res.f_exact,
            "gci_fine_pct": res.gci_fine_pct,
            "gci_coarse_pct": res.gci_coarse_pct,
            "monotonic": res.monotonic,
        }

    out["total_wall_s"] = time.perf_counter() - t_start
    ok = accept(out)
    out["accept"] = ok

    (results_dir / "validation.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    (results_dir / "validation_tables.md").write_text(
        render_tables(out), encoding="utf-8"
    )
    return out


def np_maxabs(a, b) -> float:
    import numpy as np

    return float(np.abs(np.asarray(a) - np.asarray(b)).max())


def accept(out: dict[str, Any]) -> bool:
    """Project acceptance criteria (documented in docs/)."""
    pois = out["poiseuille"]
    if not pois:
        return False
    finest = pois[-1]
    if finest["dp_error_pct"] > 5.0:  # straight-tube dp on the fine grid
        return False
    if out["womersley"].get("periodicity_last_cycle_endpoints", 1e9) > 0.05:
        return False
    for row in out["stenosis"]:
        if not (0.5 <= row["ratio_to_correlation"] <= 2.0):
            return False
    return True


def render_tables(out: dict[str, Any]) -> str:
    lines = ["# Validation tables (auto-generated)", ""]
    lines += [
        "| grid | dp (Pa) | dp err % | u_center | uc err % | iters |",
        "|---|---|---|---|---|---|",
    ]
    for r in out["poiseuille"]:
        lines.append(
            f"| {r['grid']} | {r['dp']:.3f} | {r['dp_error_pct']:.2f} "
            f"| {r['u_center']:.5f} | {r['u_center_error_pct']:.2f} "
            f"| {r['iterations']} |"
        )
    lines.append("")
    w = out["womersley"]
    lines += [
        f"Womersley: steps={w['steps']}, uc range "
        f"{w['u_center_range'][0]:.4f}..{w['u_center_range'][1]:.4f}, "
        f"flow reversal={w['flow_reversal']}",
        "",
        "| severity | dp solver (Pa) | dp excess (Pa) | dp Young-Tsai (Pa) | ratio |",
        "|---|---|---|---|---|",
    ]
    for r in out["stenosis"]:
        lines.append(
            f"| {r['severity']:.1f} | {r['dp_solver']:.2f} "
            f"| {r['dp_excess_solver']:.2f} "
            f"| {r['dp_young_tsai']:.2f} | {r['ratio_to_correlation']:.3f} |"
        )
    lines.append("")
    g = out.get("gci") or {}
    if g:
        lines += [f"GCI ({g['quantity']}): {g.get('gci_fine_pct', 'n/a')}"]
    return "\n".join(lines) + "\n"
