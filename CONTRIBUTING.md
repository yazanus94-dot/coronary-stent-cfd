# Contributing to Coronary-Stent-CFD

Thank you for considering a contribution to this project! This is a hemodynamic
simulation study of stenosed coronary arteries with stent effect analysis,
built with an in-house axisymmetric finite-volume solver and validated against
Poiseuille, Womersley, Young-Tsai, and grid-convergence references.

## How to contribute

### 1. Fork and clone

```bash
git clone https://github.com/yazanus94-dot/coronary-stent-cfd.git
cd coronary-stent-cfd
```

### 2. Set up the development environment

```bash
python -m venv .venv
# Windows: .venv\Scripts\activate
# Linux/macOS: source .venv/bin/activate
pip install -e .[dev]
```

### 3. Run the checks

```bash
# Lint
ruff check src tests

# Unit tests
pytest tests -q

# Quick validation (CI-style)
stent-hemo validate --quick
```

### 4. Submit a pull request

- Create a feature branch from `main`
- Make your changes
- Add or update tests as needed
- Ensure all existing tests pass
- Open a pull request with a clear description of what changed and why

## Code style

- Python 3.10+ with type hints
- Line length 88 (Ruff default)
- Docstrings on all public modules and functions
- Meaningful commit messages

## Reporting issues

Open a GitHub issue for bugs, feature requests, or questions about the
hemodynamics model. Please include your operating system, Python version,
and steps to reproduce.

## License

MIT Licensed. See [LICENSE](LICENSE).