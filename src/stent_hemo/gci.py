"""Grid Convergence Index (Roache 1994) for mesh-convergence studies.

Given three grid levels of a quantity f (fine, medium, coarse) with grid
refinement ratios, compute the observed order of accuracy p, the
Richardson-extrapolated value, and the GCI for the fine grid:

    p        = ln(|(f3 - f2)/(f2 - f1)|) / ln(r)
    f_exact  = f1 + (f1 - f2) / (r^p - 1)
    GCI_fine = Fs * |(f2 - f1) / f1| / (r^p - 1) * 100 %

with Fs = 1.25 the standard safety factor.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

FS = 1.25  # Roache's safety factor


@dataclass
class GCIResult:
    """Outcome of a three-grid convergence study."""

    p: float  # observed order of accuracy
    f_fine: float
    f_exact: float  # Richardson-extrapolated value
    gci_fine_pct: float  # GCI of the fine grid (%)
    gci_coarse_pct: float
    monotonic: bool


def grid_convergence_index(
    f_fine: float, f_medium: float, f_coarse: float, r: float = 1.3
) -> GCIResult:
    """Compute the GCI for three uniform refinement levels.

    Args:
        f_fine: Quantity on the fine grid.
        f_medium: Quantity on the medium grid.
        f_coarse: Quantity on the coarse grid.
        r: Constant grid refinement ratio (h_coarse / h_medium = h_medium /
            h_fine = r).

    Raises:
        ValueError: If r <= 1 or the differences are degenerate (no
            observable convergence trend).
    """
    if r <= 1.0:
        raise ValueError("refinement ratio must be > 1")
    d21 = f_medium - f_fine
    d32 = f_coarse - f_medium
    if abs(d21) < 1e-14 or abs(d32) < 1e-14:
        raise ValueError("grid differences are degenerate; no convergence trend")
    p = math.log(abs(d32 / d21)) / math.log(r)
    denom = r**p - 1.0
    if denom <= 0:
        raise ValueError("non-convergent refinement (r^p <= 1)")
    f_exact = f_fine - d21 / denom
    gci_f = FS * abs(d21 / f_fine) / denom * 100.0
    gci_c = FS * abs(d32 / f_fine) / denom * 100.0
    monotonic = (d21 * d32) > 0
    return GCIResult(
        p=p,
        f_fine=f_fine,
        f_exact=f_exact,
        gci_fine_pct=gci_f,
        gci_coarse_pct=gci_c,
        monotonic=monotonic,
    )
