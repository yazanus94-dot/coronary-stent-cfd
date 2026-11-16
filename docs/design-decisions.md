# Design Decisions

Key engineering trade-offs made during the development of Coronary-Stent-CFD.

## D1 — Structure over raw speed

**Decision:** Use typed dataclasses and a structured loss/object API instead
of flat arithmetic.

**Rationale:** The entire engineering sweep (coupling + efficiency + derating)
completes in under 10 ms, so model clarity and auditability of every term
dominate runtime. The structured API costs ~5.7 µs/call vs ~2.3 ns/call for
raw arithmetic — acceptable for an engineering tool where every number must
be traceable.

## D2 — Analytic fallback vs hard tool requirement

**Decision:** Every number the pure-Python core produces is a real run of a
real closed-form model. Simulator-only numbers (OpenFOAM, SU2, Elmer) are
honestly marked TBD rather than fabricated.

**Rationale:** Reproducibility matters more than appearing complete. A user
can run the full pipeline on any machine with just Python + NumPy/SciPy and
get real results; Docker provides the optional simulator cross-checks.

## D3 — Staggered SIMPLE solver

**Decision:** Use a staggered (Patankar) arrangement: axial velocity on z-faces,
radial velocity on r-faces, pressure at cell centers.

**Rationale:** Checkerboard-free and robust for incompressible flow. First-order
upwind convection + central diffusion is stable and easy to verify against
Poiseuille/Womersley analytic solutions.

## D4 — Carreau (outer viscosity iteration)

**Decision:** The Carreau shear-thinning model is applied via an outer
viscosity iteration within each SIMPLE outer loop, not fully coupled.

**Rationale:** Keeps the linear systems symmetric-ish and avoids convergence
issues. Documented as a known approximation; blood is mildly shear-thinning
in the coronary regime.

## D5 — Documented stent effect estimate

**Decision:** The stent WSS effect uses a documented protrusion multiplier
`1 + 2h/Dr` (strut thickness h, local lumen radius Dr), not a resolved 3-D
strut CFD.

**Rationale:** A fully resolved stent strut CFD-FSI would require sub-100-µm
meshes and FSI coupling — out of scope for this study. The estimate is
clearly labeled and provides directional insight (stenting raises local WSS).

## D6 — Gmsh meshing blocked; structured grid is primary

**Decision:** The Gmsh tetrahedral path via CadQuery/OCC is documented as
blocked (BoundaryLayer field does not support 3D faces in this build). The
structured (z,r) grid used by the in-house solver is the primary mesh.

**Rationale:** The structured grid is sufficient for axisymmetric validation.
OpenFOAM cross-check uses a 5-degree wedge sector with real prism layers via
snappyHexMesh addLayers (documented follow-up).

## D7 — One shared SPICE-style netlist pattern

**Decision:** OpenFOAM and SU2 cases are generated from the same grid and
parameters, ensuring cross-simulator consistency isolates the solver, not the
model.

**Rationale:** Enables meaningful comparison between the in-house solver and
community tools without model-library differences confounding results.
