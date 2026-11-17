"""Gmsh mesher: wall-refined tetrahedral CFD meshes.

Imports the CadQuery STEP solids, identifies inlet/outlet/wall/stent
surfaces, tags physical groups, and builds tetrahedral meshes with a
dual distance-field size control:

- wall refinement: first cell size ~1 % of the radius at the wall,
  growing to the bulk size (replaces boundary-layer prisms, which this
  gmsh build's BoundaryLayer field does not support for 3D faces; the
  OpenFOAM stage adds real prism layers with snappyHexMesh addLayers);
- throat refinement: 0.5x bulk size within 6 R of the stenosis.

Writes Gmsh 2.2 meshes at three levels (coarse / medium / fine) for the
convergence study (GCI), and an SU2-format mesh for the cross-check.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from .core import VesselGeometry


@dataclass(frozen=True)
class MeshLevel:
    """One refinement level."""

    name: str
    char_length: float  # m
    wall_first: float  # m, first-cell size at the wall (~1 % of R)
    throat_ratio: float = 0.5


MESH_LEVELS: tuple[MeshLevel, ...] = (
    MeshLevel("coarse", 150.0e-6, 60.0e-6),
    MeshLevel("medium", 100.0e-6, 40.0e-6),
    MeshLevel("fine", 70.0e-6, 28.0e-6),
)

# physical-group ids used consistently across meshing and the OF case
PHYS_ID = {"inlet": 11, "outlet": 12, "wall": 13, "stent": 14}


def _classify_surfaces(model, vessel: VesselGeometry) -> dict[str, list[int]]:
    """Split surfaces into inlet / outlet / wall / stent by geometry."""
    inlet, outlet, wall, stent = [], [], [], []
    for dim, tag in model.getEntities(2):
        com = model.occ.getCenterOfMass(dim, tag)
        if abs(com[0]) < 1e-6:
            inlet.append(tag)
        elif abs(com[0] - vessel.length) < 1e-6:
            outlet.append(tag)
        else:
            # strut faces carry a radius beyond the smooth lumen profile
            r_com = math.hypot(com[1], com[2])
            if r_com > vessel.radius_at(com[0]) + 0.25 * vessel.r0:
                stent.append(tag)
            else:
                wall.append(tag)
    return {"inlet": inlet, "outlet": outlet, "wall": wall, "stent": stent}


def mesh_vessel(
    step_path: Path,
    out_msh: Path,
    vessel: VesselGeometry,
    level: MeshLevel,
    with_stent: bool,
    also_su2: bool = False,
) -> str:
    """Generate one wall-refined tetrahedral mesh.

    Returns a status message; raises ImportError when gmsh is missing.
    """
    try:
        import gmsh
    except ImportError as exc:
        raise ImportError(
            "gmsh NOT INSTALLED - CFD meshes are TBD (install: pip install gmsh)"
        ) from exc

    if not step_path.exists():
        return (
            f"SKIP {out_msh.name}: {step_path.name} missing "
            "(run the geometry stage first)"
        )

    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.occ.importShapes(str(step_path))
        gmsh.model.occ.synchronize()

        groups = _classify_surfaces(gmsh.model, vessel)
        for name, tags in groups.items():
            if tags:
                gid = PHYS_ID[name]
                gmsh.model.addPhysicalGroup(2, tags, gid)
                gmsh.model.setPhysicalName(2, gid, name)

        gmsh.option.setNumber("Mesh.MeshSizeMin", level.wall_first * 0.5)
        gmsh.option.setNumber("Mesh.MeshSizeMax", level.char_length)
        gmsh.option.setNumber("Mesh.Algorithm3D", 10)  # HXT (robust Delaunay)
        gmsh.option.setNumber("Mesh.RandomFactor", 1e-5)

        # --- size field 1: wall refinement (distance to all wetted faces)
        wetted = groups["wall"] + (groups["stent"] if with_stent else [])
        f_wall = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f_wall, "SurfacesList", wetted)
        f_wall_t = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f_wall_t, "IField", f_wall)
        gmsh.model.mesh.field.setNumber(f_wall_t, "LcMin", level.wall_first)
        gmsh.model.mesh.field.setNumber(f_wall_t, "LcMax", level.char_length)
        gmsh.model.mesh.field.setNumber(f_wall_t, "DistMin", 0.0)
        gmsh.model.mesh.field.setNumber(f_wall_t, "DistMax", 2.0 * vessel.r0)

        # --- size field 2: stenosis/stent region refinement
        f_dist2 = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumbers(f_dist2, "SurfacesList", wetted)
        f_throat = gmsh.model.mesh.field.add("Threshold")
        gmsh.model.mesh.field.setNumber(f_throat, "IField", f_dist2)
        gmsh.model.mesh.field.setNumber(f_throat, "LcMin", level.throat_ratio * level.char_length)
        gmsh.model.mesh.field.setNumber(f_throat, "LcMax", level.char_length)
        gmsh.model.mesh.field.setNumber(f_throat, "DistMin", 1.0 * vessel.r0)
        gmsh.model.mesh.field.setNumber(f_throat, "DistMax", 6.0 * vessel.r0)

        f_min = gmsh.model.mesh.field.add("Min")
        gmsh.model.mesh.field.setNumbers(f_min, "FieldsList", [f_wall_t, f_throat])
        gmsh.model.mesh.field.setAsBackgroundMesh(f_min)

        gmsh.model.mesh.generate(3)
        n_elem = sum(
            len(tags) for tags in gmsh.model.mesh.getElements(3)[1]
        )
        out_msh.parent.mkdir(parents=True, exist_ok=True)
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.write(str(out_msh))
        if also_su2:
            gmsh.write(str(out_msh.with_suffix(".su2")))
    finally:
        gmsh.finalize()

    extra = " + SU2" if also_su2 else ""
    return (
        f"Wrote {out_msh}{extra}: {n_elem} elements, "
        f"wall first cell {level.wall_first * 1e6:.0f} um, "
        f"bulk {level.char_length * 1e6:.0f} um"
    )
