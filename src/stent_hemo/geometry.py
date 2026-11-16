"""Parametric geometry: vessel, wall shell, lattice stent (CadQuery).

Builds four solids (all in metres, X = axial direction):
- lumen: the blood domain of the unstented vessel (cosine stenosis)
- wall: the vessel wall shell (thickness wall_thickness around the lumen)
- stent: a sinusoidal-strut lattice stent that props the stenosis open
- stented lumen: lumen restored to nominal radius over the stent span,
  minus the protruding part of the struts (struts partly embedded)

Exports STEP files and a rendered SVG projection.
"""

from __future__ import annotations

import math
from pathlib import Path

from .core import StentParams, VesselGeometry


def _radius_samples(v: VesselGeometry, n: int = 220) -> list[tuple[float, float]]:
    """Sampled (x, r) profile of the lumen for the revolve operation."""
    xs = [v.length * i / (n - 1) for i in range(n)]
    return [(x, v.radius_at(x)) for x in xs]


def _lumen_profile_workplane(pts: list[tuple[float, float]]):
    """Build the revolve profile workplane from sampled (x, r) points.

    The lumen arc is an OCC interpolating spline (not a polyline):
    hundreds of nearly-collinear polyline segments create degenerate
    triangles in Gmsh's 3D Delaunay recovery ('vertex lies in a
    segment'). The spline adds ~1.2 % volume versus the analytic
    integral; documented in docs/design-decisions.md.
    """
    import cadquery as cq

    l_end = pts[-1][0]
    top = list(reversed(pts))  # (L, r(L)) ... (0, r0)
    return (
        cq.Workplane("XY")
        .moveTo(0.0, 0.0)
        .lineTo(l_end, 0.0)
        .lineTo(l_end, top[0][1])
        .spline(top[1:], includeCurrent=True)
        .close()
    )


def build_lumen(v: VesselGeometry):
    """Revolved lumen solid (blood domain of the unstented vessel)."""
    return _lumen_profile_workplane(_radius_samples(v)).revolve(
        360, (0, 0, 0), (1, 0, 0)
    )


def build_wall(v: VesselGeometry, wall_thickness: float = 0.5e-3):
    """Wall shell: outer body minus the lumen (annular tube)."""

    outer_pts = [(x, v.radius_at(x) + wall_thickness) for x, _ in _radius_samples(v)]
    outer = _lumen_profile_workplane(outer_pts).revolve(360, (0, 0, 0), (1, 0, 0))
    return outer.cut(build_lumen(v))


def stent_span(v: VesselGeometry, s: StentParams) -> tuple[float, float]:
    """Axial (start, end) of the stent, centred on the stenosis."""
    cover = s.cover_frac * v.length
    return v.xc - cover / 2.0, v.xc + cover / 2.0


def build_stent(v: VesselGeometry, s: StentParams):
    """Straight-strut zigzag lattice stent with connector bands.

    Each ring is a set of ``crowns`` straight struts alternately tilted
    +/- one crown arc around the circumference (classic zigzag ring);
    rings are joined by thin circumferential connector bands. Struts
    occupy a radial band [R - embed*t, R + (1-embed)*t] so they are
    partly embedded in the wall and partly protrude into the lumen.
    One oriented box per strut, fused in a single boolean operation.
    """
    import cadquery as cq

    x0, x1 = stent_span(v, s)
    r_surf = v.r0  # stent props the vessel open at the nominal radius
    t = s.strut_thickness
    w = s.strut_width
    r_inner = r_surf - s.embed_frac * t
    r_outer = r_surf + (1.0 - s.embed_frac) * t
    r_mid = 0.5 * (r_inner + r_outer)
    ring_pitch = (x1 - x0) / s.n_rings
    arc = 2.0 * math.pi / s.crowns  # crown arc around the circumference

    shapes: list[cq.Shape] = []
    for ring in range(s.n_rings):
        x_start = x0 + ring * ring_pitch
        direction = 1.0 if ring % 2 == 0 else -1.0
        for crown in range(s.crowns):
            th0 = 2.0 * math.pi * crown / s.crowns
            th1 = th0 + direction * arc
            p0 = (x_start, r_mid * math.cos(th0), r_mid * math.sin(th0))
            p1 = (
                x_start + ring_pitch,
                r_mid * math.cos(th1),
                r_mid * math.sin(th1),
            )
            shapes.append(_segment_box(p0, p1, w, t))

    # connector bands between adjacent rings
    for ring in range(s.n_rings - 1):
        joint_x = x0 + (ring + 1) * ring_pitch
        band = (
            cq.Workplane("XZ")
            .workplane(offset=joint_x)
            .center(r_mid, 0.0)
            .rect(w, t)
            .revolve(360, (joint_x, 0, 0), (joint_x, 1, 0))
        )
        shapes.append(band.val())

    fused = shapes[0].fuse(*shapes[1:])
    return cq.Workplane("XY").newObject([fused])


def _segment_box(
    p0: tuple[float, float, float],
    p1: tuple[float, float, float],
    w: float,
    t: float,
):
    """A thin box spanning p0->p1 with cross-section w (tangential) x t
    (radial), built at the origin, rotated onto the segment, translated
    to its midpoint (slightly lengthened so neighbours overlap cleanly).
    """
    import cadquery as cq

    dx, dy, dz = p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]
    length = math.sqrt(dx * dx + dy * dy + dz * dz)
    seg_len = length + 1.6 * t  # overlap for a watertight fuse
    box = cq.Solid.makeBox(seg_len, w, t)
    box = box.translate(
        cq.Vector(-seg_len / 2.0, -w / 2.0, -t / 2.0)
    )
    d = cq.Vector(dx, dy, dz).normalized()
    x_axis = cq.Vector(1.0, 0.0, 0.0)
    cos_a = max(-1.0, min(1.0, x_axis.dot(d)))
    ang = math.degrees(math.acos(cos_a))
    if ang > 1e-9:
        axis = x_axis.cross(d)
        if axis.Length < 1e-9:  # antiparallel
            axis = cq.Vector(0.0, 1.0, 0.0)
        box = box.rotate(cq.Vector(0, 0, 0), axis, ang)
    mid = cq.Vector((p0[0] + p1[0]) / 2.0, (p0[1] + p1[1]) / 2.0, (p0[2] + p1[2]) / 2.0)
    return box.translate(mid)


def build_stented_lumen(v: VesselGeometry, s: StentParams):
    """Lumen restored to nominal radius over the stent span, minus struts.

    This models stenting: the stent holds the stenosis open at the nominal
    diameter; the protruding strut parts locally reduce the area.
    """

    x0, x1 = stent_span(v, s)
    pts = _radius_samples(v)
    restored = [
        (x, max(r, v.r0) if x0 <= x <= x1 else r) for x, r in pts
    ]
    profile = _lumen_profile_workplane(restored)
    return profile.revolve(360, (0, 0, 0), (1, 0, 0)).cut(build_stent(v, s))


def export_geometry(
    out_dir: Path,
    v: VesselGeometry,
    s: StentParams,
) -> dict[str, str]:
    """Export all solids as STEP plus an SVG render of the vessel profile."""
    try:
        import cadquery as cq
    except ImportError as exc:
        raise ImportError(
            "CadQuery NOT INSTALLED - geometry STEP files are TBD "
            "(install: pip install cadquery)"
        ) from exc

    out_dir.mkdir(parents=True, exist_ok=True)
    msgs: dict[str, str] = {}
    solids = {
        "vessel_lumen": build_lumen(v),
        "vessel_wall": build_wall(v),
        "stent": build_stent(v, s),
        "stented_lumen": build_stented_lumen(v, s),
    }
    for name, solid in solids.items():
        path = out_dir / f"{name}.stp"
        cq.exporters.export(solid, str(path), exportType="STEP")
        msgs[name] = f"Wrote {path}"

    svg = out_dir / "vessel_profile.svg"
    cq.exporters.export(solids["vessel_lumen"], str(svg), opt={
        "width": 800, "height": 240, "marginLeft": 20, "marginTop": 20,
    })
    msgs["render"] = f"Wrote {svg}"
    return msgs
