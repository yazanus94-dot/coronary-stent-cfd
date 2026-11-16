# Security Considerations

## Threat model

Coronary-Stent-CFD is a local scientific computing pipeline. It performs no
network access, executes no user-supplied code, and writes only to
`results/`, `kicad/`, `cad/`, `mesh/`.

## Security boundaries

- **No secrets:** No API keys, tokens, or credentials are stored in the
  repository. The `LICENSE`, README, and source files contain no sensitive
  material.
- **No network access:** The pipeline runs entirely offline. Docker pulls
  `openfoam/openfoam11-paraview510` only when explicitly requested.
- **No user-supplied code execution:** All inputs are hard-coded parameters
  or files within the repository. No `eval`, `exec`, or shell interpolation
  of user input.
- **Subprocess safety:** External tools (Gmsh, OpenFOAM, SU2, CalculiX) are
  invoked with explicit argument lists (no shell interpolation). Subprocess
  timeouts (300–600 s) prevent runaway solves.
- **File writes confined:** Outputs are written only to `results/`,
  `kicad/`, `cad/`, `mesh/`. No writes to system directories or home.

## Known limitations

- STEP/CAD files (`.stp`) are binary formats; treat them as untrusted if
  obtained from external sources.
- OpenFOAM Docker image is pulled from Docker Hub; verify image integrity
  in production environments.
- No dependency pinning beyond minimum versions in `pyproject.toml`; use a
  lock file (`pip freeze` or `poetry.lock`) for reproducible builds.

## Reporting security issues

Open a GitHub issue or contact the repository owner. Do not disclose
potential vulnerabilities in public issues; use private disclosure when
appropriate.
