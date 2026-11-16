# Coronary-Stent-CFD

**Hemodynamics of a Stenosed Coronary Artery and the Effect of a Stent**

Pulsatile, shear-thinning blood flow in a stenosed artery (50 % and 70 % area stenosis), a parametric lattice stent, and the change in wall shear stress, pressure drop and recirculation after stenting. An in-house axisymmetric finite-volume solver produces the reference numbers; OpenFOAM, SU2 and CalculiX decks are generated for cross-check runs in Docker.

| Component | Tool | Status |
|-----------|------|--------|
| In-house axisymmetric FV solver (SIMPLE, staggered) | Python + NumPy/SciPy | ✅ implemented + validated |
| Analytic anchors: Poiseuille, Womersley, Young-Tsai, GCI | Python | ✅ implemented + validated |
| CAD geometry (lumen, wall, stent, stented lumen) | CadQuery → STEP + SVG | ✅ implemented |
| Mesh writers: Gmsh, SU2, OpenFOAM wedge polyMesh | Python | ✅ implemented |
| OpenFOAM v11 cross-check cases | Docker (`openfoam/openfoam11-paraview510`) | ⏳ generated, TBD run |
| SU2 axisymmetric configs | SU2 | ✅ generated |
| CalculiX stent-strut structural model | CalculiX | ✅ generated |
| V&V suite (accept = true) | pytest | ✅ 4 tests pass |

## Headline results (in-house solver)

| Case | Δp steady (Pa) | WSS max (Pa) | TAWSS max (Pa) | OSI max |
|------|----------------|--------------|----------------|---------|
| 50 % stenosis | 240.42 | 7.135 | 7.828 | 0.038 |
| 70 % stenosis | 399.30 | 18.688 | 19.010 | 0.300 |

Stent effect (documented protrusion estimate `1 + 2h/Dr`): WSS max 18.688 Pa (clean lumen) vs 21.612 Pa (stented estimate) at 70 % stenosis.

Grid-convergence index (GCI, fine grid): **1.93 %**, monotonic convergence. Validation acceptance: **PASS**.

OpenFOAM cross-check: **TBD** — run `docker compose run openfoam` and drop `results/openfoam_results.csv`.

## Repository layout

```
Coronary-Stent-CFD/
├── pyproject.toml          # package metadata (stent-hemo 0.1.0)
├── Dockerfile, docker-compose.yml, docker/run_openfoam.sh
├── .github/workflows/ci.yml
├── src/stent_hemo/         # the Python package
│   ├── core.py             # blood, vessel, waveform, stent parameters
│   ├── fv_solver.py        # axisymmetric SIMPLE solver (steady + transient)
│   ├── geometry.py         # CadQuery STEP/SVG geometry builder
│   ├── mesh_io.py          # Gmsh / SU2 / OpenFOAM polyMesh writers
│   ├── meshing.py          # Gmsh wall-refined tetrahedral mesher
│   ├── external_cases.py   # OpenFOAM / SU2 / CalculiX case writers
│   ├── validation.py       # Poiseuille / Womersley / Young-Tsai / GCI
│   ├── postprocess.py      # WSS/TAWSS/OSI metrics + results.md + figures
│   ├── pipeline.py         # one-command: mesh / cases / validate / postprocess / all
│   ├── struct_grid.py      # structured axisymmetric (z, r) grid
│   ├── young_tsai.py       # stenosis pressure-drop correlation
│   ├── hand_calc.py        # closed-form hand calculations
│   ├── womersley.py        # Womersley analytic pulsatile solution
│   └── gci.py              # Roache grid-convergence index
├── tests/test_core.py      # unit tests
├── cad/                    # generated STEP solids + vessel profile SVG
├── mesh/                   # generated Gmsh + SU2 meshes
├── openfoam/               # generated OpenFOAM wedge cases (50 %, 70 %, wedge)
├── su2/                    # SU2 configs
├── calculix/               # CalculiX stent-strut .inp
├── results/                # validation.json, results.md, figures/*.png
└── tmp_of_check/           # scratch (not committed)
```

## Installation

```bash
cd Coronary-Stent-CFD
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install -e .[dev]
```

Optional: `pip install cadquery gmsh` for the geometry and mesh stages (auto-detected; gracefully skipped when absent).

## Running the pipeline

```bash
# Everything: mesh, external cases, validation, postprocessing (~5 min laptop)
python -m stent_hemo pipeline all

# Or individual stages:
python -m stent_hemo pipeline mesh      # Gmsh + SU2 meshes
python -m stent_hemo pipeline cases     # OpenFOAM / SU2 / CalculiX decks
python -m stent_hemo pipeline validate  # V&V suite (writes results/validation.json)
python -m stent_hemo pipeline postprocess  # results.md + figures

# Quick validation (coarse grid only, for CI):
stent-hemo validate --quick
```

## Docker (full reproducible environment, incl. OpenFOAM)

```bash
docker compose up        # builds image, runs the full pipeline
docker compose down
docker compose run openfoam   # runs the generated OpenFOAM wedge cases
```

OpenFOAM results are written to `results/openfoam_results.csv` and consumed by postprocessing.

## Testing

```bash
pytest tests -q                 # unit tests
stent-hemo validate --quick     # quick V&V
ruff check src tests            # lint
```

## Reproducibility

All numbers in `results/` are produced by the in-house solver from the committed source. The OpenFOAM cross-check requires the Docker image (`openfoam/openfoam11-paraview510`) and is marked TBD until run.

## License

MIT — see LICENSE (to be added).