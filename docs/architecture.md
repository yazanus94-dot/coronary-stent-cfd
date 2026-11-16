# Architecture

## Overview

Coronary-Stent-CFD is a hemodynamic simulation study of stenosed coronary
arteries with stent effect analysis. The core is an in-house axisymmetric
finite-volume (FV) SIMPLE solver that computes pressure drop, wall shear
stress (WSS), time-averaged WSS (TAWSS), and oscillatory shear index (OSI)
for 50 % and 70 % area stenosis cases. OpenFOAM, SU2, and CalculiX decks are
generated for cross-check runs in Docker.

## Module diagram

```mermaid
flowchart TD
    Core["core.py\nBlood, Vessel, Waveform, Stent params"]:::core
    Geometry["geometry.py\nCadQuery STEP/SVG"]:::geo
    Mesh["mesh_io.py\nGmsh / SU2 / OpenFOAM writers"]:::mesh
    Solver["fv_solver.py\nAxisymmetric SIMPLE FV"]:::solver
    Validation["validation.py\nPoiseuille/Womersley/Young-Tsai/GCI"]:::val
    Post["postprocess.py\nWSS/TAWSS/OSI + results.md + figures"]:::post
    Pipeline["pipeline.py\nmesh / cases / validate / postprocess / all"]:::pipe
    External["external_cases.py\nOpenFOAM / SU2 / CalculiX writers"]:::ext

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

## Data flow

1. **`core.py`** defines all physical parameters (blood rheology, vessel
   geometry, flow waveform, stent parameters) as frozen dataclasses — the
   single source of truth.
2. **`geometry.py`** builds CAD solids (lumen, wall, stent, stented lumen)
   with CadQuery and exports STEP + SVG.
3. **`mesh_io.py`** writes Gmsh, SU2, and OpenFOAM wedge polyMesh from the
   structured axisymmetric grid.
4. **`fv_solver.py`** solves steady or transient incompressible Navier-Stokes
   on a staggered grid (SIMPLE, first-order upwind, central diffusion).
5. **`validation.py`** runs the V&V suite against analytic/correlation
   references and writes `results/validation.json`.
6. **`postprocess.py`** computes WSS/TAWSS/OSI metrics, generates figures, and
   writes `results/results.md`.
7. **`pipeline.py`** orchestrates all stages with a single command.

## External solvers

| Solver | Purpose | Status |
|--------|---------|--------|
| OpenFOAM v11 | 3-D wedge cross-check (incompressibleFluid) | Docker image |
| SU2 | 2-D axisymmetric incompressible config | Generated |
| CalculiX | Stent-strut structural beam model | Generated |

External binaries are auto-detected; missing tools print a clear message and
are marked TBD rather than crashing the pipeline.
