"""Mesh writers for the structured axisymmetric grid.

Three formats are produced from one :class:`~stent_hemo.struct_grid.StructGrid`:

* Gmsh ``.msh`` (2.2 ASCII, quads, physical groups) - general interchange.
* SU2 ``.su2`` (2-D axisymmetric-ready quads with INLET/OUTLET/WALL markers).
* OpenFOAM ``constant/polyMesh`` (5-degree wedge sector, hexes, wedge patches)
  for ``foamRun -solver incompressibleFluid`` (OpenFOAM Foundation v11).

All writers are pure Python (no CAD kernel): the geometry is the structured
(z, r) grid itself, which is also what the in-house FV solver uses. The Gmsh
tetrahedral path via cadquery/OCC remains available in ``meshing.py`` but is
documented as blocked (see docs/design-decisions.md).
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np

from .struct_grid import StructGrid

# physical/marker ids shared by the writers
INLET_ID, OUTLET_ID, WALL_ID = 11, 12, 13


def _node_index(nz: int, nr: int) -> np.ndarray:
    """Node id of (i, k) for the (nz+1, nr+1) structured node grid."""
    idx = np.arange((nz + 1) * (nr + 1)).reshape(nz + 1, nr + 1)
    return idx


def _nodes(grid: StructGrid) -> np.ndarray:
    """Node coordinates (nz+1, nr+1, 2) as (z, r)."""
    zf = grid.z_faces
    rf = grid.r_face[:, None] * grid.r_frac[None, :]
    z = np.broadcast_to(zf[:, None], rf.shape)
    return np.stack([z, rf], axis=-1)


def write_gmsh(grid: StructGrid, path: str | Path) -> Path:
    """Write a 2-D structured quad mesh in Gmsh 2.2 ASCII format."""
    nz, nr = grid.shape
    pts = _nodes(grid).reshape(-1, 2)
    idx = _node_index(nz, nr)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["$MeshFormat", "2.2 0 8", "$EndMeshFormat"]
    lines.append(f"$Nodes\n{len(pts)}")
    for i, (z, r) in enumerate(pts, 1):
        lines.append(f"{i} {z:.9e} {r:.9e} 0")
    lines.append("$EndNodes")
    quads = []
    for j in range(nz):
        for k in range(nr):
            n = (idx[j, k], idx[j + 1, k], idx[j + 1, k + 1], idx[j, k + 1])
            quads.append((j, k, [x + 1 for x in n]))
    lines.append(f"$Elements\n{len(quads) + 2 * (nz + nr)}")
    eid = 1
    # volume (surface in 2D) quads: physical 1
    for _j, _k, n in quads:
        lines.append(f"{eid} 3 1 1 {' '.join(map(str, n))}")
        eid += 1
    # inlet (j=0), outlet (j=nz), wall (k=nr) line elements
    for k in range(nr):
        lines.append(
            f"{eid} 2 {INLET_ID} {INLET_ID} {idx[0, k] + 1} {idx[0, k + 1] + 1}"
        )
        eid += 1
    for k in range(nr):
        lines.append(
            f"{eid} 2 {OUTLET_ID} {OUTLET_ID} {idx[nz, k] + 1} {idx[nz, k + 1] + 1}"
        )
        eid += 1
    for j in range(nz):
        lines.append(
            f"{eid} 2 {WALL_ID} {WALL_ID} {idx[j, nr] + 1} {idx[j + 1, nr] + 1}"
        )
        eid += 1
    lines += ["$EndElements", ""]
    path.write_text("\n".join(lines), encoding="utf-8")
    return path


def write_su2(grid: StructGrid, path: str | Path) -> Path:
    """Write an SU2 ASCII mesh (2-D, quads, 3 boundary markers).

    Use ``AXISYMMETRIC=YES`` in the .cfg to reuse the (z, r) mesh as an
    axisymmetric sector of unit angle.
    """
    nz, nr = grid.shape
    pts = _nodes(grid).reshape(-1, 2)
    idx = _node_index(nz, nr)
    quads = []
    for j in range(nz):
        for k in range(nr):
            quads.append(
                [idx[j, k], idx[j + 1, k], idx[j + 1, k + 1], idx[j, k + 1]]
            )
    inlet = [idx[0, k] for k in range(nr)] + [idx[0, nr]]
    outlet = [idx[nz, k] for k in range(nr)] + [idx[nz, nr]]
    wall = [idx[j, nr] for j in range(nz)] + [idx[nz, nr]]
    n_elems = len(quads) + (nr + 1) + (nr + 1) + (nz + 1)
    lines = [
        "NDIME= 2",
        f"NELEM= {n_elems}",
        "NEARBF= 0",
        f"NPOIN= {len(pts)}",
        "NPOIN2D= %d" % len(pts),
    ]
    for a in quads:
        lines.append("%d\t%d\t%d\t%d\t9" % (a[0], a[1], a[2], a[3]))
    for k in range(nr):  # inlet quads (2-node lines as degenerate quads)
        a, b = inlet[k], inlet[k + 1]
        lines.append("%d\t%d\t-1\t-1\t3" % (a, b))
    for k in range(nr):
        a, b = outlet[k], outlet[k + 1]
        lines.append("%d\t%d\t-1\t-1\t3" % (a, b))
    for j in range(nz):
        a, b = wall[j], wall[j + 1]
        lines.append("%d\t%d\t-1\t-1\t3" % (a, b))
    lines.append("")
    for i, (z, r) in enumerate(pts):
        lines.append(f"{z:.9e}\t{r:.9e}\t{i}")
    lines += [
        "NMARK= 3",
        "MARKER_TAG= INLET",
        f"MARKER_ELEMS= {nr}",
        "",
        *[f"{inlet[k]}\t{inlet[k + 1]}" for k in range(nr)],
        "MARKER_TAG= OUTLET",
        f"MARKER_ELEMS= {nr}",
        "",
        *[f"{outlet[k]}\t{outlet[k + 1]}" for k in range(nr)],
        "MARKER_TAG= WALL",
        f"MARKER_ELEMS= {nz}",
        "",
        *[f"{wall[j]}\t{wall[j + 1]}" for j in range(nz)],
    ]
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_openfoam_polymesh(grid: StructGrid, case_dir: str | Path) -> Path:
    """Write a 5-degree wedge polyMesh into ``case_dir/constant/polyMesh``.

    The (z, r) grid is extruded into a narrow sector spanning
    ``[-wedge/2, +wedge/2]`` around the x-axis; OpenFOAM's ``wedge`` boundary
    condition turns it into the axisymmetric equivalent.
    """
    nz, nr = grid.shape
    half = math.radians(grid.wedge_deg) / 2.0
    pts = _nodes(grid).reshape(-1, 2)  # (z, r)
    npz = len(pts)
    idx = _node_index(nz, nr)

    # points: two theta planes (back at -half on axis x, front at +half)
    def xyz(z: float, r: float, th: float) -> tuple[float, float, float]:
        return z, r * math.cos(th), r * math.sin(th)

    points: list[tuple[float, float, float]] = []
    for th in (-half, half):  # plane 0: back, plane 1: front
        for z, r in pts:
            points.append(xyz(z, r, th))
    def pid(i: int, k: int, plane: int) -> int:
        return plane * npz + idx[i, k]

    # internal faces: hexahedra cells (j, k) with 4 back + 4 front nodes
    # ordering: ccw when viewed from inside (OpenFOAM owner/neighbour)
    cells = {}
    for j in range(nz):
        for k in range(nr):
            cells[(j, k)] = (
                pid(j, k, 0), pid(j + 1, k, 0), pid(j + 1, k + 1, 0),
                pid(j, k + 1, 0),
                pid(j, k, 1), pid(j + 1, k, 1), pid(j + 1, k + 1, 1),
                pid(j, k + 1, 1),
            )
    faces: list[tuple[int, ...]] = []
    owners: list[int] = []
    neighbours: list[int] = []

    def add_face(nodes4: tuple[int, ...], owner: int, neigh: int | None) -> None:
        faces.append(nodes4)
        owners.append(owner)
        neighbours.append(-1 if neigh is None else neigh)

    cell_id = lambda j, k: j * nr + k  # noqa: E731

    # axial interior faces (between (j-1,k) and (j,k)), j = 1..nz-1
    for j in range(1, nz):
        for k in range(nr):
            add_face(
                (pid(j, k, 0), pid(j, k + 1, 0), pid(j, k + 1, 1), pid(j, k, 1)),
                cell_id(j - 1, k), cell_id(j, k),
            )
    # radial interior faces (between (j,k-1) and (j,k)), k = 1..nr-1
    for j in range(nz):
        for k in range(1, nr):
            add_face(
                (pid(j, k, 0), pid(j + 1, k, 0), pid(j + 1, k, 1), pid(j, k, 1)),
                cell_id(j, k - 1), cell_id(j, k),
            )
    n_internal = len(faces)
    # boundary faces: inlet (j=0), outlet (j=nz), wall (k=nr), wedge planes
    inlet_f, outlet_f, wall_f = [], [], []
    for k in range(nr):
        f = (pid(0, k, 0), pid(0, k + 1, 0), pid(0, k + 1, 1), pid(0, k, 1))
        add_face(f, cell_id(0, k), None)
        inlet_f.append(len(faces) - 1)
    for k in range(nr):
        f = (pid(nz, k, 0), pid(nz, k + 1, 0), pid(nz, k + 1, 1), pid(nz, k, 1))
        add_face(f, cell_id(nz - 1, k), None)
        outlet_f.append(len(faces) - 1)
    for j in range(nz):
        f = (pid(j, nr, 0), pid(j + 1, nr, 0), pid(j + 1, nr, 1), pid(j, nr, 1))
        add_face(f, cell_id(j, nr - 1), None)
        wall_f.append(len(faces) - 1)
    # wedge planes: generated plane by plane so the boundary patches
    # wedge_back (theta=-half) and wedge_front (+half) are contiguous and
    # in increasing startFace order, as OpenFOAM requires.
    wedge_back: list[int] = []
    wedge_front: list[int] = []
    for j in range(nz):
        for k in range(nr):
            f = (
                pid(j, k, 0), pid(j + 1, k, 0), pid(j + 1, k + 1, 0),
                pid(j, k + 1, 0),
            )
            add_face(f, cell_id(j, k), None)
            wedge_back.append(len(faces) - 1)
    for j in range(nz):
        for k in range(nr):
            f = (
                pid(j, k, 1), pid(j, k + 1, 1), pid(j + 1, k + 1, 1),
                pid(j + 1, k, 1),
            )
            add_face(f, cell_id(j, k), None)
            wedge_front.append(len(faces) - 1)

    pm = Path(case_dir) / "constant" / "polyMesh"
    pm.mkdir(parents=True, exist_ok=True)
    with open(pm / "points", "w", encoding="utf-8") as fh:
        fh.write(f"/* points: {len(points)} */\nFoamFile\n{{\n format ascii;\n class vectorField;\n object points;\n}}\n")
        fh.write(f"{len(points)}\n(\n")
        for p in points:
            fh.write("(%.9e %.9e %.9e)\n" % p)
        fh.write(")\n")
    with open(pm / "faces", "w", encoding="utf-8") as fh:
        fh.write("FoamFile\n{\n format ascii;\n class faceList;\n object faces;\n}\n")
        fh.write(f"{len(faces)}\n(\n")
        for f in faces:
            fh.write("4(%d %d %d %d)\n" % f)
        fh.write(")\n")
    with open(pm / "owner", "w", encoding="utf-8") as fh:
        fh.write("FoamFile\n{\n format ascii;\n class labelList;\n object owner;\n}\n")
        fh.write(f"{len(owners)}\n(\n")
        fh.write("\n".join(map(str, owners)))
        fh.write("\n)\n")
    internal_neighbours = neighbours[:n_internal]
    with open(pm / "neighbour", "w", encoding="utf-8") as fh:
        fh.write("FoamFile\n{\n format ascii;\n class labelList;\n object neighbour;\n}\n")
        fh.write(f"{len(internal_neighbours)}\n(\n")
        fh.write("\n".join(map(str, internal_neighbours)))
        fh.write("\n)\n")
    boundary_patches = [
        ("inlet", "patch", inlet_f),
        ("outlet", "patch", outlet_f),
        ("wall", "wall", wall_f),
        ("wedge_back", "wedge", wedge_back),
        ("wedge_front", "wedge", wedge_front),
    ]
    with open(pm / "boundary", "w", encoding="utf-8") as fh:
        fh.write("FoamFile\n{\n format ascii;\n class polyBoundaryMesh;\n object boundary;\n}\n")
        fh.write(f"{len(boundary_patches)}\n(\n")
        for name, typ, flist in boundary_patches:
            fh.write(
                f"{name}\n{{\n type {typ};\n nFaces {len(flist)};\n startFace {flist[0]};\n}}\n"
            )
        fh.write(")\n")
    return pm


def face_shear_metrics(
    grid: StructGrid,
    u_hist: np.ndarray,
    dt: float,
    mu: float,
) -> dict[str, np.ndarray]:
    """Wall shear stress metrics over the last cardiac cycle.

    u_hist: (n_steps, nz, nr) cell-centered axial velocity history.
    Returns dict with tws (TAWSS), osi, rrt arrays of shape (nz,).
    WSS is estimated from the wall-adjacent cell: tau = mu * du/dr at
    r = r_wall - dr/2 (two-point difference to the no-slip wall).
    """
    nz, nr = grid.shape
    rc = grid.r_center
    dr_wall = rc / nr * 0.5  # cell center to wall distance
    # u at wall-adjacent ring center
    uw = u_hist[:, :, nr - 1]
    # du/dr at wall ~ (0 - uw) / (dr/2)
    wss = mu * np.abs(uw) / dr_wall[None, :]  # (n_steps, nz), magnitude
    tawss = wss.mean(axis=0)
    # OSI needs the signed WSS direction; the axisymmetric model only has the
    # axial component, so OSI is computed on the signed axial shear proxy
    dirw = mu * uw / dr_wall[None, :]
    mean = dirw.mean(axis=0)
    osi = 0.5 * (1.0 - np.abs(mean) / np.maximum(tawss, 1e-12))
    rrt = 1.0 / np.maximum((1.0 - 2.0 * osi) * tawss, 1e-12)
    return {"tawss": tawss, "osi": osi, "rrt": rrt, "wss_hist": wss}
