"""One-command study pipeline.

``python -m stent_hemo.pipeline [all|mesh|cases|validate|postprocess]``

Runs the deliverables of the study in dependency order:

1. ``mesh``       - Gmsh/SU2 meshes for the 50% and 70% cases.
2. ``cases``      - OpenFOAM / SU2 / CalculiX decks for cross-check runs.
3. ``validate``   - the in-house V&V suite (Poiseuille, Womersley,
                    Young-Tsai, GCI).
4. ``postprocess``- headline metrics, figures and results.md.

``all`` runs everything (about 5 minutes on a laptop CPU).
"""

from __future__ import annotations

import sys
from pathlib import Path

from .core import BloodProperties, FlowWaveform, StentParams, VesselGeometry
from .external_cases import write_calculix_inp, write_openfoam_case, write_su2_cfg
from .mesh_io import write_gmsh, write_su2
from .struct_grid import StructGrid


def build_meshes(mesh_dir: str | Path = "mesh") -> list[Path]:
    """Write Gmsh + SU2 meshes for both stenosis severities."""
    mesh_dir = Path(mesh_dir)
    out = []
    for sev in (0.5, 0.7):
        grid = StructGrid(VesselGeometry(severity=sev, length=90e-3), 120, 24)
        out.append(write_gmsh(grid, mesh_dir / f"lumen_sev{int(sev*100)}.msh"))
        out.append(write_su2(grid, mesh_dir / f"lumen_sev{int(sev*100)}.su2"))
    return out


def build_cases(
    openfoam_dir: str | Path = "openfoam",
    su2_dir: str | Path = "su2",
    calculix_dir: str | Path = "calculix",
) -> list[Path]:
    """Write external-solver decks (OpenFOAM wedge, SU2, CalculiX strut)."""
    blood = BloodProperties()
    wf = FlowWaveform()
    out = []
    for sev in (0.5, 0.7):
        grid = StructGrid(VesselGeometry(severity=sev, length=90e-3), 120, 24)
        out.append(
            write_openfoam_case(
                grid, blood, wf, Path(openfoam_dir) / f"stenosis_{int(sev*100)}"
            )
        )
        out.append(
            write_su2_cfg(
                grid, blood, wf, Path(su2_dir) / f"stenosis_{int(sev*100)}.cfg"
            )
        )
    out.append(write_calculix_inp(StentParams(), Path(calculix_dir) / "stent_strut.inp"))
    return out


def main(argv: list[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    step = argv[0] if argv else "all"
    quick = "--quick" in argv

    if step in ("mesh", "all"):
        paths = build_meshes()
        print(f"[mesh] wrote {len(paths)} files")
    if step in ("cases", "all"):
        paths = build_cases()
        print(f"[cases] wrote {len(paths)} decks")
    if step in ("validate", "all"):
        from .validation import run_validation

        out = run_validation("results", quick=quick)
        status = "PASS" if out["accept"] else "FAIL"
        print(f"[validate] {status} ({out['total_wall_s']:.0f}s)")
        if not out["accept"]:
            return 1
    if step in ("postprocess", "all"):
        from .postprocess import run_postprocess

        out = run_postprocess("results", "results/figures")
        for c in out["cases"]:
            print(
                f"[post] {int(c['severity']*100)}%: dp={c['dp_steady']:.1f} Pa "
                f"WSSmax={c['wss_steady_max']:.2f} Pa OSI={c['osi_max']:.3f}"
            )
        print("[post] results/results.md + results/figures written")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
