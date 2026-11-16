# Changelog

All notable changes to this project are documented here.
The format follows [Keep a Changelog](https://keepachangelog.com/) and the
project adheres to [Semantic Versioning](https://semver.org/).

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
