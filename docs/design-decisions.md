# Design Decisions

Key engineering trade-offs made during the development of Coronary-Stent-CFD.

## D1 — Structure over raw speed

**Decision:** Use typed dataclasses and a structured solver API instead of
flat arithmetic.

**Rationale:** The production steady solve completes in ~5 s and the full
V&V suite in ~80 s, so model clarity and auditability of every term dominate
runtime. The structured API costs little and keeps every number traceable.

## D2 — Analytic fallback vs hard tool requirement

**Decision:** Every number the pure-Python core produces is a real run of a
real closed-form or numerically solved model. Simulator-only numbers
(OpenFOAM, SU2) are honestly marked TBD rather than fabricated.

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

## D5 — Geometry-derived protrusion estimate (not strut-resolved CFD)

**Decision:** The stent WSS effect uses an engineering protrusion estimate
built from the actual stent geometry: protrusion height `h = t*(1 - embed_frac)`
of each strut into the lumen and a nominal WSS multiplier `1 + 2h/Dr` with
`Dr` the local lumen radius. It is explicitly NOT a strut-resolved CFD.

**Rationale:** A fully resolved stent strut CFD-FSI would require sub-100-µm
meshes around the struts and FSI coupling — out of scope for this study. The
estimate uses actual strut parameters and the actual lumen radii, is clearly
labeled, and provides directional insight (stenting raises local WSS). The
documented limitation is that flow around individual struts is not resolved;
resolving it is a documented follow-up (pending the OpenFOAM 3-D wedge
cross-check or an FSI coupling).

## D6 — Gmsh meshing blocked; structured grid is primary

**Decision:** The Gmsh tetrahedral path via CadQuery/OCC is documented as
blocked (BoundaryLayer field does not support 3D faces in this build). The
structured (z,r) grid used by the in-house solver is the primary mesh.

**Rationale:** The structured grid is sufficient for axisymmetric validation.
OpenFOAM cross-check uses a 5-degree wedge sector with real prism layers via
snappyHexMesh addLayers (documented follow-up).

## D7 — One shared grid for all external solvers

**Decision:** OpenFOAM and SU2 cases are generated from the same grid and
parameters, ensuring cross-simulator consistency isolates the solver, not the
model.

**Rationale:** Enables meaningful comparison between the in-house solver and
community tools without model-library differences confounding results.
