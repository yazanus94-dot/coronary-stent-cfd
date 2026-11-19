"""Post-processing: WSS metrics, plots, and results.md generation.

Produces the study deliverables for the two CFD cases (50% and 70% stenosis)
and the stent-surface comparison:

* ``results/results.md`` - the headline numbers and tables.
* ``results/figures/*.png`` - velocity/waveform/WSS metric plots.

All quantities are computed from the in-house axisymmetric solver; OpenFOAM
cross-check numbers are inserted from ``openfoam/results.csv`` when present
(see README - Reproducibility), otherwise marked as pending (TBD).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from .core import BloodProperties, FlowWaveform, StentParams, VesselGeometry
from .fv_solver import AxisymSolver
from .mesh_io import face_shear_metrics
from .struct_grid import StructGrid
from .young_tsai import stenosis_pressure_drop


def _centerline_u(u_cell: np.ndarray) -> np.ndarray:
    """Axis-adjacent ring of the cell-centered field (nz, nr)."""
    return u_cell[:, 0]


def _wall_wss_steady(
    grid: StructGrid, u_cell: np.ndarray, mu: float
) -> np.ndarray:
    """Steady WSS magnitude (Pa) per axial cell from the wall-adjacent ring.

    ``u_cell`` is the cell-centered velocity (nz, nr); the wall-adjacent
    value sits at fraction (nr-0.5)/nr of the local lumen radius.
    """
    dr_wall = 0.5 * grid.r_center / grid.nr
    uw = np.abs(u_cell[:, grid.nr - 1])
    return mu * uw / dr_wall


def run_postprocess(
    results_dir: str | Path = "results",
    fig_dir: str | Path = "results/figures",
    fine_transient: bool = False,
) -> dict[str, Any]:
    """Run both stenosis cases + stented variant and write all outputs.

    Steady fields use the production grid (120x24). Transient WSS metrics
    default to a 60x12 grid with dt = 2 ms (Courant-stable); 120x24 needs
    dt <= 1 ms (verified stable, ~7x slower) - enable via ``fine_transient``.
    """
    results_dir = Path(results_dir)
    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)
    blood = BloodProperties()
    wf = FlowWaveform()
    out: dict[str, Any] = {"cases": []}

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    for sev in (0.5, 0.7):
        vessel = VesselGeometry(severity=sev, length=90e-3)
        grid = StructGrid(vessel, 120, 24)
        solver = AxisymSolver(grid, blood)
        sol_steady = solver.solve(wf, steady=True, max_outer=600, tol=1e-6)
        uc = sol_steady.u
        wss = _wall_wss_steady(grid, uc, blood.mu_newtonian)
        z_mm = grid.z_centers * 1e3
        r_profile = grid.r_center * 1e3

        sol_trans_grid = grid if fine_transient else StructGrid(vessel, 60, 12)
        trans_dt = 1e-3 if fine_transient else 2e-3
        solver_trans = AxisymSolver(sol_trans_grid, blood)
        sol_trans = solver_trans.solve(wf, steady=False, dt=trans_dt, n_cycles=3)
        hist = sol_trans.u_history
        assert hist is not None
        metrics = face_shear_metrics(
            sol_trans_grid, hist, trans_dt, blood.mu_newtonian
        )
        if np.isnan(metrics["tawss"]).any() or np.isnan(metrics["osi"]).any():
            raise RuntimeError(
                f"transient WSS metrics diverged (sev={sev}); "
                "reduce dt or refine the grid (Courant limit)"
            )

        case = {
            "severity": sev,
            "dp_steady": sol_steady.dp,
            "dp_young_tsai": stenosis_pressure_drop(vessel, blood, wf.mean_flow()).dp_total,
            "u_center_max_steady": float(uc.max()),
            "wss_steady_max": float(wss.max()),
            "wss_steady_throat": float(wss[np.argmin(grid.r_center)]),
            "tawss_max": float(metrics["tawss"].max()),
            "osi_max": float(metrics["osi"].max()),
            "rrt_max": float(metrics["rrt"].max()),
            "transient_grid": f"{sol_trans_grid.nz}x{sol_trans_grid.nr}",
            "transient_dt_s": trans_dt,
        }
        out["cases"].append(case)

        # ---- figures -----------------------------------------------
        fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
        ax = axes[0]
        for frac, lbl in ((0.05, "near wall"), (0.5, "mid radius"), (1.0, "centerline")):
            k = min(int(frac * grid.nr), grid.nr - 1)
            ax.plot(z_mm, uc[:, k], label=lbl)
        ax.axvline(vessel.xc * 1e3, color="gray", ls="--", lw=0.8)
        ax.set_xlabel("z (mm)")
        ax.set_ylabel("u (m/s)")
        ax.set_title(f"Steady velocity, {int(sev*100)}% stenosis")
        ax.legend(fontsize=8)
        ax = axes[1]
        ax.plot(z_mm, r_profile, color="tab:red")
        ax.set_xlabel("z (mm)")
        ax.set_ylabel("lumen radius (mm)")
        ax.set_title("Lumen profile")
        ax = axes[2]
        ax.plot(z_mm, wss, color="tab:purple")
        ax.set_xlabel("z (mm)")
        ax.set_ylabel("WSS (Pa)")
        ax.set_title("Steady wall shear stress")
        fig.tight_layout()
        fig.savefig(fig_dir / f"steady_sev{int(sev*100)}.png", dpi=150)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6, 3.6))
        z = np.arange(hist.shape[0]) * trans_dt * wf.heart_rate / 60.0
        ax.plot(z, hist[:, :, 0].max(axis=1), label="solver uc (t)")
        ts = np.linspace(0, z[-1], 200)
        qt = np.array([wf.flow(t) for t in ts])
        ax.plot(ts, qt / qt.max() * hist.max(), "--", label="Q(t) shape")
        ax.set_xlabel("time (periods)")
        ax.set_ylabel("u_center (m/s)")
        ax.set_title(f"Transient inlet response, {int(sev*100)}% stenosis")
        ax.legend(fontsize=8)
        fig.tight_layout()
        fig.savefig(fig_dir / f"transient_sev{int(sev*100)}.png", dpi=150)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6, 3.6))
        zt_mm = sol_trans_grid.z_centers * 1e3
        ax.plot(zt_mm, metrics["tawss"], label="TAWSS")
        ax.set_xlabel("z (mm)")
        ax.set_ylabel("TAWSS (Pa)")
        ax.set_title(f"Time-averaged WSS, {int(sev*100)}% stenosis")
        fig.tight_layout()
        fig.savefig(fig_dir / f"tawss_sev{int(sev*100)}.png", dpi=150)
        plt.close(fig)

    # ---- stent effect (geometry-derived protrusion estimate) -----------
    # The stented-WSS estimate is built from the ACTUAL stent geometry
    # (strut thickness, embed fraction, stent span) and the ACTUAL lumen
    # radii of the 70% stenosis geometry: each strut protrudes by
    # h = t*(1 - embed) into the lumen, and the nominal wall shear is
    # scaled by the protrusion ratio (1 + 2h/Dr) with Dr the local lumen
    # radius. This is an engineering estimate, NOT a strut-resolved CFD
    # run -- flow around individual struts requires the 3-D wedge
    # simulation or an FSI coupling (see docs/design-decisions.md, D5).
    stent = StentParams()
    vessel70 = VesselGeometry(severity=0.7, length=90e-3)
    grid = StructGrid(vessel70, 120, 24)
    solver = AxisymSolver(grid, blood)
    sol = solver.solve(wf, steady=True, max_outer=600, tol=1e-6)
    uc = sol.u
    wss_clean = _wall_wss_steady(grid, uc, blood.mu_newtonian)
    rc = grid.r_center
    zh = np.linspace(0.0, vessel70.length, len(rc))
    in_stent = (np.abs(zh - vessel70.xc) <= stent_span_len(stent, vessel70) / 2)
    h = stent.strut_thickness * (1.0 - stent.embed_frac)  # protrusion height
    wss_multiplier = 1.0 + 2.0 * h / np.maximum(rc, 1e-9)  # nominal ratio
    wss_stent = np.where(in_stent, wss_clean * wss_multiplier, wss_clean)
    out["stent_effect"] = {
        "model": (
            "geometry-derived protrusion estimate: nominal multiplier "
            "1 + 2h/Dr with h = strut_thickness*(1 - embed_frac) and Dr = "
            "local lumen radius; not a strut-resolved CFD"
        ),
        "strut_protrusion_h_m": float(h),
        "stent_span_m": float(stent_span_len(stent, vessel70)),
        "wss_max_clean": float(wss_clean.max()),
        "wss_max_stented": float(wss_stent.max()),
        "tawss_ratio_throat": float(
            wss_stent[np.argmin(rc)] / max(wss_clean[np.argmin(rc)], 1e-12)
        ),
        "limitation": (
            "No strut-resolved flow simulation: flow around individual "
            "struts requires the 3-D wedge OpenFOAM run or an FSI "
            "coupling; see docs/design-decisions.md (D5) and the "
            "Limitations section of the README."
        ),
    }

    fig, ax = plt.subplots(figsize=(6, 3.6))
    ax.plot(zh * 1e3, wss_clean, label="clean lumen")
    ax.plot(zh * 1e3, wss_stent, label="stented (est.)")
    ax.set_xlabel("z (mm)")
    ax.set_ylabel("WSS (Pa)")
    ax.set_title("70% stenosis: clean vs stented WSS (geometry-derived estimate)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(fig_dir / "stent_wss_70.png", dpi=150)
    plt.close(fig)

    out["openfoam_cross_check"] = _read_openfoam(results_dir)
    out["figures"] = sorted(p.name for p in fig_dir.glob("*.png"))

    (results_dir / "postprocess.json").write_text(
        json.dumps(out, indent=2), encoding="utf-8"
    )
    _write_results_md(out, results_dir)
    return out


def stent_span_len(stent: StentParams, vessel: VesselGeometry) -> float:
    return stent.cover_frac * vessel.length


def _read_openfoam(results_dir: Path) -> dict[str, Any]:
    """Insert OpenFOAM cross-check numbers if the Docker run produced them."""
    csv = results_dir / "openfoam_results.csv"
    if csv.exists():
        rows = csv.read_text().strip().splitlines()
        return {"source": str(csv), "rows": rows}
    return {
        "status": "TBD",
        "note": "Run the OpenFOAM case in Docker (docker compose run openfoam) "
        "and drop the extracted dp / WSS into results/openfoam_results.csv.",
    }


def _write_results_md(out: dict[str, Any], results_dir: Path) -> None:
    lines = [
        "# Results (auto-generated by stent_hemo.postprocess)",
        "",
        "All numbers from the in-house structured axisymmetric FV solver",
        "(SIMPLE, staggered grid). See docs/design-decisions.md for the",
        "discretization choices and results/validation.json for the V&V run.",
        "",
        "| case | dp steady (Pa) | Young-Tsai excess (Pa) | u_center max (m/s) | WSS max (Pa) | TAWSS max (Pa) | OSI max |",
        "|---|---|---|---|---|---|---|",
    ]
    for c in out["cases"]:
        lines.append(
            f"| {int(c['severity']*100)}% stenosis | {c['dp_steady']:.2f} "
            f"| {c['dp_young_tsai']:.2f} | {c['u_center_max_steady']:.4f} "
            f"| {c['wss_steady_max']:.3f} | {c['tawss_max']:.3f} "
            f"| {c['osi_max']:.3f} |"
        )
    st = out["stent_effect"]
    lines += [
        "",
        f"Stent effect ({st['model']}): WSS max {st['wss_max_clean']:.3f} Pa (clean)",
        f"vs {st['wss_max_stented']:.3f} Pa (stented estimate).",
        "",
        "## OpenFOAM cross-check",
        "",
        json.dumps(out["openfoam_cross_check"], indent=2),
        "",
        "## Figures",
        "",
    ]
    lines += [f"- results/figures/{name}" for name in out["figures"]]
    (results_dir / "results.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
