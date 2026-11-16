"""Structured axisymmetric grid shared by the FV solver and mesh writers.

A cylindrical (z, r) grid: uniform axial spacing, radial cells uniform
in the FRACTION of the local radius (so radial size follows the local
stenosis). The wedge mesh for OpenFOAM adds a 5-degree theta sector;
the SU2 writer emits the plain 2D (z, r) mesh for AXISYMMETRIC=YES.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import VesselGeometry


@dataclass
class StructGrid:
    """Structured axisymmetric grid for a stenosed lumen."""

    vessel: VesselGeometry
    nz: int
    nr: int
    wedge_deg: float = 5.0  # sector angle for the OpenFOAM wedge mesh

    # ---- axial -------------------------------------------------------
    @property
    def z_faces(self) -> np.ndarray:
        return np.linspace(0.0, self.vessel.length, self.nz + 1)

    @property
    def z_centers(self) -> np.ndarray:
        z = self.z_faces
        return 0.5 * (z[:-1] + z[1:])

    # ---- radial (per axial column) ------------------------------------
    def radius_at_faces(self, z: np.ndarray | float) -> np.ndarray:
        """Lumen radius at arbitrary axial positions."""
        return self.vessel.radius_at(np.atleast_1d(z))

    @property
    def r_frac(self) -> np.ndarray:
        """Radial cell-boundary fractions i/Nr, shape (Nr+1,)."""
        return np.linspace(0.0, 1.0, self.nr + 1)

    @property
    def r_center(self) -> np.ndarray:
        """Lumen radius at cell-center axial positions, shape (nz,)."""
        return self.radius_at_faces(self.z_centers)

    @property
    def r_face(self) -> np.ndarray:
        """Lumen radius at axial face positions, shape (nz+1,)."""
        return self.radius_at_faces(self.z_faces)

    def cell_r_bounds(self) -> tuple[np.ndarray, np.ndarray]:
        """Inner/outer radial bounds of each cell, shape (nz, nr).

        Cell k spans local radius fractions [k/Nr, (k+1)/Nr] of the
        column radius (radius at the cell-center axial position).
        """
        f = self.r_frac
        rc = self.r_center
        r_in = rc[:, None] * f[None, :-1]
        r_out = rc[:, None] * f[None, 1:]
        return r_in, r_out

    # ---- FV geometry ---------------------------------------------------
    @property
    def cell_volume(self) -> np.ndarray:
        """Cell volumes (per radian), shape (nz, nr)."""
        r_in, r_out = self.cell_r_bounds()
        dz = self.vessel.length / self.nz
        return 0.5 * dz * (r_out**2 - r_in**2)

    @property
    def axial_face_area(self) -> np.ndarray:
        """Ring areas of axial faces (per radian), shape (nz+1, nr)."""
        f = self.r_frac
        rf = self.r_face
        return 0.5 * rf[:, None] ** 2 * (f[None, 1:] ** 2 - f[None, :-1] ** 2)

    def radial_face_area(self) -> np.ndarray:
        """Cylindrical areas of radial faces (per radian), shape (nz, nr+1)."""
        dz = self.vessel.length / self.nz
        return dz * self.r_center[:, None] * self.r_frac[None, :]

    @property
    def shape(self) -> tuple[int, int]:
        return (self.nz, self.nr)

    # ---- mesh metrics ---------------------------------------------------
    @property
    def n_cells(self) -> int:
        return self.nz * self.nr

    @property
    def cell_size(self) -> tuple[float, float]:
        """(mean axial, mean radial at inlet) cell size in metres."""
        return self.vessel.length / self.nz, self.vessel.r0 / self.nr
