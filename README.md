# Coronary-Stent-CFD

**A professional Python engineering project demonstrating hemodynamics of a stenosed coronary artery and the effect of a stent.**

![Python](https://img.shields.io/badge/Python-3.10+-blue)
![Docker](https://img.shields.io/badge/Docker-Supported-blue)
![Testing](https://img.shields.io/badge/Tests-Pytest-green)
![License](https://img.shields.io/badge/License-MIT-green)

Pulsatile, shear-thinning blood flow in a stenosed artery (50 % and 70 % area stenosis), a parametric lattice stent, and the change in wall shear stress, pressure drop and recirculation after stenting. An in-house axisymmetric finite-volume solver produces the reference numbers; OpenFOAM, SU2 and CalculiX decks are generated for cross-check runs in Docker.

## 1. Project Overview

**Problem solved:** Blood flow through a stenosed coronary artery creates complex hemodynamic patterns â€” elevated wall shear stress at the throat, flow separation downstream, and oscillatory shear that correlates with plaque progression and stent failure. A stent props the vessel open but introduces strut-induced flow disturbance. Quantifying these effects requires coupled fluid-structure analysis across multiple simulation tools.

**Why it is challenging:** The answer requires coupled physics (pulsatile Navier-Stokes + non-Newtonian Carreau blood rheology), grid-convergence verification (GCI), cross-validation between independent solvers (in-house FV vs OpenFOAM vs SU2), a parametric stent geometry (CadQuery STEP), and a reproducible Docker environment.

**Engineering solution:** A one-command reproducible pipeline: analytic hand calculations as ground truth â†’ axisymmetric SIMPLE FV solver â†’ Poiseuille/Womersley/Young-Tsai/GCI validation â†’ CadQuery STEP geometry â†’ Gmsh/SU2/OpenFOAM mesh writers â†’ post-processing with WSS/TAWSS/OSI metrics â†’ results and figures. Every external tool is auto-detected and gracefully skipped when not installed.

## 2. Key Results

| Case | dp steady (Pa) | Young-Tsai excess (Pa) | WSS max (Pa) | OSI max |
|------|----------------|------------------------|--------------|---------|
| 50 % stenosis | 240.4 | 56.9 (solver excess 50.3, ratio 0.88) | 7.13 | 0.04 |
| 70 % stenosis | 399.3 | 219.8 (solver excess 209.2, ratio 0.95) | 18.69 | 0.30 |

- Straight-tube verification: dp error 1.4 % on the fine grid; centerline
  velocity error 0.1 %; GCI fine-grid 1.9 % (monotonic).
- Transient Womersley flow (3 cycles, dt = 2 ms): capture of flow reversal
  (u_center range âˆ’0.027 to 0.514 m/s), periodic within 1 %.
- Stent effect (documented protrusion estimate): WSS max 18.69 â†’ 21.61 Pa at
  the 70 % throat.

## 3. Key Features

- In-house axisymmetric SIMPLE FV solver (steady + transient pulsatile flow)
- Non-Newtonian Carreau blood rheology with outer viscosity iteration
- Analytic anchors: Poiseuille, Womersley, Young-Tsai pressure-drop correlation
- Grid Convergence Index (Roache) verification with monotonic convergence
- CadQuery parametric stent geometry (lumen, wall, stent, stented lumen â†’ STEP + SVG)
- Gmsh, SU2, and OpenFOAM 5-degree wedge polyMesh writers
- OpenFOAM v11 cross-check cases in Docker
- CalculiX stent-strut structural beam model
- WSS / TAWSS / OSI / RRT post-processing with figures
- One-command pipeline: `python -m stent_hemo.pipeline all`
- 40+ automated tests: unit + Hypothesis property-based + integration
- Docker / docker-compose environment with OpenFOAM v11
- GitHub Actions CI: Ruff, pytest, quick V&V suite, OpenFOAM cross-check

## 4. System Architecture

![System Architecture](docs/images/architecture.png)

Data flows from `core.py` (blood, vessel, waveform, stent parameters) into the
SIMPLE solver (`fv_solver.py`) and out to post-processing (`postprocess.py`).
External simulators (OpenFOAM, SU2, CalculiX) consume the same grid via
`mesh_io.py` and `external_cases.py`. The validation suite (`validation.py`)
cross-checks against analytic references. The pipeline (`pipeline.py`) orchestrates
all stages.

```mermaid
flowchart TD
    Core["core.py\nBlood, Vessel, Waveform, Stent params"]:::core
    Geometry["geometry.py\nCadQuery STEP/SVG"]:::geo
    Mesh["mesh_io.py\nGmsh / SU2 / OpenFOAM writers"]:::mesh
    Solver["fv_solver.py\nAxisymmetric SIMPLE FV"]:::solver
    Validation["validation.py\nPoiseuille/Womersley/Young-Tsai/GCI"]:::val
    Post["postprocess.py\nWSS/TAWSS/OSI + figures"]:::post
    Pipeline["pipeline.py\nmesh / cases / validate / postprocess / all"]:::pipe
    External["external_cases.py\nOpenFOAM / SU2 / CalculiX"]:::ext

    Core --> Geometry
    Core --> Mesh
    Core --> Solver
    Core --> Validation
    Core --> Post
    Core --> External
    Mesh --> Pipeline
    Solver --> Pipeline
    Validation --> Pipeline
    Post --> Pipeline
    External --> Pipeline
    Pipeline --> Results["results/\nvalidation.json, results.md, figures/*.png"]

    classDef core fill:#e1f5fe,stroke:#01579b
    classDef geo fill:#f3e5f5,stroke:#4a148c
    classDef mesh fill:#e8f5e9,stroke:#1b5e20
    classDef solver fill:#fff3e0,stroke:#e65100
    classDef val fill:#fce4ec,stroke:#880e4f
    classDef post fill:#f1f8e9,stroke:#33691e
    classDef pipe fill:#e3f2fd,stroke:#0d47a1
    classDef ext fill:#f5f5f5,stroke:#424242
```

## 5. Technology Stack

| Technology | Purpose |
|------------|---------|
| Python 3.10+ (type hints) | Core implementation |
| NumPy / SciPy | Numerical linear algebra and special functions |
| matplotlib | Figures and plots |
| Docker / docker compose | Reproducible environment with OpenFOAM v11 |
| Pytest | Automated testing |
| Hypothesis | Property-based physics-invariant testing |
| Ruff | Lint and format |
| GitHub Actions | Continuous integration |
| OpenFOAM v11 | 3-D wedge cross-check (Docker) |
| SU2 | 2-D axisymmetric incompressible config |
| CalculiX | Stent-strut structural beam model |
| CadQuery | Parametric STEP/SVG geometry |
| Gmsh | Mesh generation (auto-detected) |

## 6. Repository Structure

```
Coronary-Stent-CFD/
â”œâ”€â”€ README.md, LICENSE, CONTRIBUTING.md, CHANGELOG.md
â”œâ”€â”€ pyproject.toml          # package + Ruff + pytest config
â”œâ”€â”€ requirements.txt        # pip dependencies
â”œâ”€â”€ Dockerfile, docker-compose.yml, .dockerignore
â”œâ”€â”€ scripts/setup.sh        # venv + install
â”œâ”€â”€ scripts/run_tests.sh    # lint + test + validate
â”œâ”€â”€ src/stent_hemo/         # the engineering library
â”‚   â”œâ”€â”€ core.py             # BloodProperties, VesselGeometry, FlowWaveform, StentParams
â”‚   â”œâ”€â”€ fv_solver.py        # Axisymmetric SIMPLE FV solver (steady + transient)
â”‚   â”œâ”€â”€ geometry.py         # CadQuery STEP/SVG geometry builder
â”‚   â”œâ”€â”€ mesh_io.py          # Gmsh / SU2 / OpenFOAM polyMesh writers
â”‚   â”œâ”€â”€ meshing.py          # Gmsh wall-refined tetrahedral mesher
â”‚   â”œâ”€â”€ external_cases.py   # OpenFOAM / SU2 / CalculiX case writers
â”‚   â”œâ”€â”€ validation.py       # Poiseuille / Womersley / Young-Tsai / GCI V&V suite
â”‚   â”œâ”€â”€ postprocess.py      # WSS/TAWSS/OSI metrics + results.md + figures
â”‚   â”œâ”€â”€ pipeline.py         # one-command: mesh / cases / validate / postprocess / all
â”‚   â”œâ”€â”€ struct_grid.py      # Structured axisymmetric (z, r) grid
â”‚   â”œâ”€â”€ young_tsai.py       # Stenosis pressure-drop correlation
â”‚   â”œâ”€â”€ hand_calc.py        # Closed-form hand calculations
â”‚   â”œâ”€â”€ womersley.py        # Womersley analytic pulsatile solution
â”‚   â””â”€â”€ gci.py              # Roache grid-convergence index
â”œâ”€â”€ tests/                  # unit / integration / property
â”‚   â”œâ”€â”€ conftest.py         # shared fixtures
â”‚   â”œâ”€â”€ unit/test_core.py
â”‚   â”œâ”€â”€ integration/test_solver.py
â”‚   â””â”€â”€ property/test_invariants.py
â”œâ”€â”€ benchmarks/             # performance benchmarks
â”‚   â”œâ”€â”€ benchmark.py
â”‚   â””â”€â”€ results.md
â”œâ”€â”€ docs/                   # architecture, decisions, debugging, security
â”‚   â”œâ”€â”€ architecture.md
â”‚   â”œâ”€â”€ design-decisions.md
â”‚   â”œâ”€â”€ debugging-log.md
â”‚   â””â”€â”€ security.md
â”œâ”€â”€ cad/                    # generated STEP solids + vessel profile SVG
â”œâ”€â”€ mesh/                   # generated Gmsh + SU2 meshes
â”œâ”€â”€ openfoam/               # generated OpenFOAM wedge cases (50 %, 70 %)
â”œâ”€â”€ su2/                    # SU2 configs
â”œâ”€â”€ calculix/               # CalculiX stent-strut .inp
â”œâ”€â”€ results/                # validation.json, results.md, figures/*.png
â””â”€â”€ .github/workflows/ci.yml
```

## 7. Installation and Setup

```bash
git clone https://github.com/yazanus94-dot/coronary-stent-cfd.git
cd coronary-stent-cfd
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install --upgrade pip
pip install -e .[dev]
```

Optional: `pip install cadquery gmsh` for the geometry and mesh stages
(auto-detected; gracefully skipped when absent).

Or with Docker (includes real OpenFOAM v11):

```bash
docker compose build
```

## 8. Execution

```bash
# One-command pipeline: mesh, external cases, validation, postprocessing
python -m stent_hemo.pipeline all

# Or individual stages:
python -m stent_hemo.pipeline mesh      # Gmsh + SU2 meshes
python -m stent_hemo.pipeline cases     # OpenFOAM / SU2 / CalculiX decks
python -m stent_hemo.pipeline validate  # V&V suite
python -m stent_hemo.pipeline postprocess  # results.md + figures

# Quick validation (coarse grids only, for CI):
stent-hemo validate --quick
```

**Expected outputs:** `results/validation.json`, `results/results.md`,
`results/validation_tables.md`, `results/postprocess.json`,
`results/figures/*.png`. Missing external tools print a clear message (e.g.
`Gmsh NOT INSTALLED`) and the pipeline still completes.

Or with Docker:

```bash
docker compose up          # runs the full pipeline
docker compose down
docker compose run openfoam   # runs the generated OpenFOAM wedge cases
```

## 9. Testing and Verification

```bash
pytest tests -q                    # all unit + integration + property tests
pytest --cov=stent_hemo --cov-report=term-missing
ruff check src tests               # lint
stent-hemo validate --quick        # quick V&V suite
```

**Actual results:** All tests pass. Validation acceptance: **PASS**.

| Test layer | What it proves |
|------------|----------------|
| Unit (`tests/unit`) | Radius bounds, Poiseuille formula match, Reynolds laminar range, Young-Tsai monotonicity, GCI quadratic recovery, mesh writer smoke, solver basics |
| Integration (`tests/integration`) | Stenosis solve finite, transient Courant stability, quick validation accepts, full grid convergence |
| Property (`tests/property`, Hypothesis) | Young-Tsai monotonic in severity, hand-calc positivity, radius positivity over random inputs |

**Validation acceptance criteria:**
- Straight-tube Poiseuille dp error < 5 % on fine grid
- Womersley periodicity error < 0.05
- Stenosis dp / Young-Tsai ratio in [0.5, 2.0]
- Grid Convergence Index: monotonic, fine-grid GCI < 5 %

**Grid-convergence index (GCI, fine grid): 1.93 %, monotonic convergence.**

## 10. Performance Benchmarking

`python benchmarks/benchmark.py` measures the solver and post-processing stages:

| Benchmark | Workload | Time | Min | Max |
|-----------|----------|------|-----|-----|
| solver_steady_20x8 | 20x8 steady solve | 0.28 s | 0.28 s | 0.28 s |
| solver_steady_40x14 | 40x14 steady solve | 1.73 s | 1.69 s | 1.83 s |
| solver_steady_80x28 | 80x28 steady solve | 2.55 s | 2.51 s | 2.65 s |
| solver_steady_120x24 | 120x24 steady solve | 4.86 s | 4.83 s | 4.93 s |
| postprocess_full | 50% + 70% + stent effect + figures | 181 s | 180 s | 182 s |
| validation_quick | Poiseuille + Womersley + stenosis | 75.6 s | â€” | â€” |

Full table: [`benchmarks/results.md`](benchmarks/results.md).

## 11. Debugging and Problem Resolution

Documented in full in [`docs/debugging-log.md`](docs/debugging-log.md). Highlights:

1. **Patankar relaxation formula error** â€” the source term applied
   `(1 - relax)/relax * diag * u_guess * relax`, silently shifting the fixed
   point of every momentum solve; converged fields were consistently wrong
   even with tiny residuals. Removing the stray `* relax` restored the exact
   Poiseuille profile on a straight tube.
2. **Double-assembled matrix coefficients** â€” the matrix-assembly helper was
   called once per neighbour *and* again when the loop reached the partner
   cell, doubling every interior coefficient (SciPy COO sums duplicates).
   Fixed by assembling only the calling cell's side per visit.
3. **Inlet/wall Dirichlet double count** â€” boundary convection went into both
   the diagonal and the RHS; only the RHS belongs there for a Dirichlet face.
4. **v-control-volume geometry** â€” the radial velocity CVs must span
   half-integer ring centers (`(f_k-1+f_k)/2`), not integer ones; the original
   widths made the mass operator and momentum operator inconsistent.
5. **Transient Courant instability** â€” dt=2 ms is stable on 60x12 but diverges
   on 120x24 through the Carreau shear-thinning feedback (throat jet ->
   viscosity drop -> faster jet); fixed with dt <= 1 ms on fine grids and a
   NaN guard in post-processing.
6. **Richardson-extrapolation sign** â€” `f_exact = f_fine + (f_medium -
   f_fine)/(r^p - 1)` had the wrong sign; the corrected form reproduces the
   analytic Poiseuille dp to 0.06 %.

## 12. Design Decisions

Full rationale in [`docs/design-decisions.md`](docs/design-decisions.md). Key trade-offs:

- **Structure over raw speed** â€” Typed dataclasses and structured API; clarity dominates when sweeps run in <10 ms.
- **Analytic fallback vs hard tool requirement** â€” Real closed-form numbers always; simulator-only results honestly marked TBD.
- **Staggered SIMPLE solver** â€” Checkerboard-free, robust for incompressible flow.
- **Carreau outer viscosity iteration** â€” Keeps linear systems manageable; documented approximation.
- **Documented stent effect estimate** â€” Protrusion multiplier `1 + 2h/Dr`, clearly labeled as order-of-magnitude.
- **Structured grid primary; Gmsh blocked** â€” Gmsh tetrahedral path documented as blocked; structured (z,r) grid is the primary mesh.
- **One shared grid for all external solvers** â€” Isolates solver differences from model differences.

## 13. Security Considerations

- No secrets in the repository (see [`docs/security.md`](docs/security.md)).
- The pipeline performs no network access (except optional Docker image pull).
- Executes no user-supplied code.
- External tools invoked with explicit argument lists (no shell interpolation); subprocess timeouts prevent runaway solves.
- Writes confined to `results/`, `cad/`, `mesh/`.

## 14. Limitations

- OpenFOAM cross-check results are **TBD** â€” run `docker compose run openfoam` to generate `results/openfoam_results.csv`.
- Stent effect uses a documented protrusion estimate (`1 + 2h/Dr`), not a resolved 3-D strut CFD-FSI.
- Gmsh tetrahedral meshing path is blocked; structured (z,r) grid is primary.
- Core model is Newtonian by default; Carreau shear-thinning via outer iteration.
- Single operating point (coronary flow ~1.2 mL/s, 75 bpm, 90 mm vessel).
- No transient thermal, no scaffold degradation, no polymer erosion.
- Research and education only â€” not a medical device.

## 15. Future Improvements

- Run the OpenFOAM cross-check automatically in CI (already wired in
  `.github/workflows/ci.yml`; needs the Docker runner to publish
  `results/openfoam_results.csv`).
- Resolve the blocked Gmsh 3-D tetrahedral path (see
  [`docs/design-decisions.md`](docs/design-decisions.md)) or adopt
  snappyHexMesh for a full 3-D stent-resolved lumen.
- Patient-specific inlet waveforms (measured Doppler tracings instead of the
  synthetic Fourier waveform).
- Fluid-structure interaction for the stent struts (couple CalculiX with the
  fluid solve).
- Realistic tapering, curvature and multi-branch coronary anatomy.
- Newtonian-vs-Carreau sensitivity study across the operating range.

## 16. License

MIT License â€” see [LICENSE](LICENSE).