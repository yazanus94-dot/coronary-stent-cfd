# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/) and the
project adheres to [Semantic Versioning](https://semver.org/).

## [1.0.2] - 2026-11-19

### Fixed
- OpenFOAM wedge cases (`external_cases.write_openfoam_case`, `mesh_io.write_openfoam_polymesh`):
  - `write_openfoam_polymesh`: `neighbour` file now has `n_internal` entries
    (previously wrote the total face count as the count but only the
    internal values — `foamRun` rejected the mesh with
    "expected int32_t, found ')'").
  - `write_openfoam_polymesh`: wedge planes are generated plane-by-plane so
    the `wedge_back`/`wedge_front` boundary patches are contiguous and in
    increasing startFace order (previously interleaved, breaking the
    boundary tiling).
  - `write_openfoam_polymesh`: boundary file uses the generated patch lists.
  - `write_openfoam_case`: pulsatile coronary inlet implemented via the
    standard `timeVaryingMappedFixedValue` BC + generated
    `constant/boundaryData/inlet` data files (points + `<time>/U` files;
    uniform waveform values since the axisymmetric inlet is uniform;
    non-collinear sample points required — the planar interpolation plane
    degenerates on the collinear inlet face centres).
  - `write_openfoam_case`: `viscosityModel Newtonian` (v11 `incompressibleFluid`
    supports only Newtonian/constant; BirdCarreau is available to
    `pimpleFluid` via MomentumTransportModels, used by the in-house FV
    solver).
  - `write_boundary_data`: raw pointField/field lists (no FoamFile header —
    the readers use `List::operator>>` which rejects headers).

### Tests
- Structural validation of both cases (mesh consistency: owner/neighbour/
  faces/boundary counts, boundary patch names and contiguity, boundaryData
  waveform data) — see the local validation script.

### OpenFOAM execution status
- Cases (`openfoam/stenosis_50`, `openfoam/stenosis_70`) generated and
  structurally validated; `foamRun` (OpenFOAM Foundation v11, Docker)
  starts the time loop (mesh parses, solver selects, first iteration) and
  then aborts in the U-momentum solve — the wedge mesh has degenerate
  cells (`checkMesh` reports zero/negative cell volumes, ~all 2880 cells;
  also "possible hole in boundary description", non-orthogonality, face
  pyramids). The wedge-mesh geometry/orientation issues are documented as
  a limitation (see README §14). Not executed to completion in this
  environment; CI runs the cases.

## [1.0.1] - 2026-11-19

### Fixed
- SIMPLE pressure-correction integration (`fv_solver.py`): the Poisson
  star-face fluxes now use the predicted radial velocity from
  `_v_momentum`, and the pressure correction updates the predicted
  field (the old `vel` field was passed in and corrected instead), so
  the radial momentum predictor is part of the SIMPLE sequence.
- `Carreau` viscosity (`_carreau_field`): rewritten for the staggered
  arrangement — interior face values are interpolated to cell centers,
  the shear-rate magnitude is the axisymmetric no-swirl tensor
  magnitude with cylindrical cell-center spacing (rc/nr), and the
  array dimensions are consistent (axial (nz-1, nr) and radial
  (nz, nr-1) inputs map to cell-centered (nz, nr) fields).
- Mass conservation (`_mass_res`): the residual is now the volumetric
  outflow imbalance in m^3/s (density factor removed), consistent with
  the `tol * |Q|` convergence tolerance; added per-cell continuity
  residuals (`_cell_continuity`, `_continuity_stats`: net_outflow,
  cell_max_abs, cell_rms, cell_rel_max).
- Transient periodicity validation: the solver stores two consecutive
  cycles and `validation.py` compares identical phases of consecutive
  cycles instead of the previous start-vs-end comparison within the
  final cycle.
- GCI methodology: consistent refinement ratios (nz 20 → 40 → 80);
  full-mode acceptance now requires a valid GCI result with monotonic
  convergence; quick mode unchanged.
- `postprocess.py` stent effect: geometry-derived protrusion estimate
  built from actual strut parameters (h = strut_thickness*(1 - embed_frac));
  limitation documented (not strut-resolved CFD).
- Docs: README numbers regenerated to match the committed
  `results/validation.json` / `results/results.md`; removed double-encoded
  characters; `docs/design-decisions.md` D5 updated.
- CI: `.github/workflows/ci.yml` now triggers on the default branch
  `master` (was `main`).

### Tests
- Added: Newtonian + Carreau mode execution tests (steady + transient)
  with continuity residual assertions.

### Results
- Poiseuille three-grid study present: 20x8 dp = 58.867 (7.12%),
  40x14 dp = 61.243 (3.37%), 80x28 dp = 62.394 (1.55%); r = 2.0,
  GCI fine = 2.16 %, monotonic convergence.
- Womersley: phase-matched periodicity 2e-15 (two stored consecutive
  cycles agree exactly at identical phases).
- Stenosis: 50 % Δp = 240.90 Pa (excess 50.77, ratio 0.89), 70 % Δp =
  401.33 Pa (excess 211.19, ratio 0.96).
- Stent effect: WSS max 20.07 → 23.21 Pa at the 70 % throat.
- Validation acceptance: PASS (full three-grid).
- OpenFOAM: not executed in this environment (Docker daemon
  unavailable); case files structurally generated and validated.

## [1.0.0] - 2026-11-16

### Added
- In-house axisymmetric SIMPLE FV solver (`fv_solver.py`) — steady + transient
  pulsatile flow with Carreau blood rheology
- Analytic core: Poiseuille, Womersley, Young-Tsai pressure-drop correlation,
  hand calculations, Roache GCI
- CadQuery parametric geometry (`geometry.py`) — lumen, wall, stent,
  stented lumen → STEP + SVG
- Mesh writers (`mesh_io.py`) — Gmsh, SU2, OpenFOAM 5-degree wedge polyMesh
- Gmsh wall-refined tetrahedral mesher (`meshing.py`)
- External case writers (`external_cases.py`) — OpenFOAM v11, SU2, CalculiX
- V&V suite (`validation.py`) — Poiseuille, Womersley, Young-Tsai, GCI
- Post-processing (`postprocess.py`) — WSS/TAWSS/OSI/RRT, results.md, figures
- One-command pipeline (`pipeline.py`) — mesh / cases / validate / postprocess / all
- Structured axisymmetric grid (`struct_grid.py`)
- Docker / docker-compose with OpenFOAM v11
- GitHub Actions CI: Ruff, pytest, quick V&V, OpenFOAM cross-check
- 40+ automated tests: unit + Hypothesis property-based + integration
- Professional repository format: CONTRIBUTING.md, requirements.txt,
  .dockerignore, scripts/, docs/, benchmarks/, reorganized tests/

### Results
- 50 % stenosis: Δp = 240.42 Pa, WSS max = 7.135 Pa, TAWSS max = 7.828 Pa, OSI max = 0.038
- 70 % stenosis: Δp = 399.30 Pa, WSS max = 18.688 Pa, TAWSS max = 19.010 Pa, OSI max = 0.300
- Stent effect (estimate): WSS max 18.688 Pa (clean) vs 21.612 Pa (stented)
- GCI fine grid: 1.93 %, monotonic convergence
- Validation acceptance: PASS

### OpenFOAM cross-check
- TBD — run `docker compose run openfoam` to generate `results/openfoam_results.csv`
