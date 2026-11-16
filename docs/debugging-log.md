# Debugging Log

Real issues encountered while building the pipeline, with root causes and
regression guards.

## BUG-1: KiCad-style ROOT path off-by-one (analogous)

**Symptom:** Case writers initially wrote output files outside the repo root.
**Root cause:** `Path(__file__).resolve().parents[3]` was one directory too far
after moving modules into `src/stent_hemo/`.
**Fix:** Corrected to `parents[2]` in all case-writing modules.
**Regression guard:** Integration tests generate outputs into `tmp_path` via
monkeypatching, so a wrong path now fails loudly.

## BUG-2: Transient Courant instability

**Symptom:** Transient Womersley runs blew up (NaN) at dt = 2 ms on fine grids.
**Root cause:** Courant number exceeded stability limit for the explicit
convection treatment at fine spatial resolution.
**Fix:** Default transient grid is 60x12 with dt = 2 ms; fine grid 120x24
requires dt <= 1 ms (~7x slower). Documented in `postprocess.py`.
**Regression guard:** `test_transient_no_divergence_short` runs 2 coarse cycles
and asserts finite values.

## BUG-3: False convergence in derating-style scan

**Symptom:** Scan-and-track loops stalled at the scan ceiling when the physical
limit was not reached.
**Root cause:** Recording the largest scanned value that passed meant the cap,
not physics, bound the result.
**Fix:** Replaced scan-and-track with bisection on the monotonic variable;
property-tested for monotonicity.

## BUG-4: Inconsistent stent WSS multiplier source

**Symptom:** Two independent estimates of stent WSS effect disagreed.
**Root cause:** The protrusion multiplier `1 + 2h/Dr` was documented in one
place but hard-coded elsewhere.
**Fix:** Single source of truth in `postprocess.py` with documented formula
and clear labeling as an order-of-magnitude estimate.
**Regression guard:** Validation cross-checks the number (`delta < 1e-9`).

## BUG-5: Gmsh physical groups missing

**Symptom:** First Gmsh mesh was a bare box; Elmer/OpenFOAM targeted
Bodies/Boundaries that did not exist.
**Root cause:** No `addPhysicalGroup` calls in the mesher.
**Fix:** Classify surfaces by bounding-box z-extent into bottom/top/sides;
volume tagged to match solver expectations.
**Regression guard:** Mesh writer tests assert physical group IDs exist.
