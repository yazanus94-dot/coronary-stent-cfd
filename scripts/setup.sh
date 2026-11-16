#!/usr/bin/env bash
# Setup script for the Coronary-Stent-CFD development environment.
# Creates a virtual environment and installs the package with dev dependencies.

set -euo pipefail

echo "[setup] Creating virtual environment..."
python -m venv .venv

echo "[setup] Activating virtual environment..."
# shellcheck disable=SC1091
source .venv/bin/activate

echo "[setup] Installing package with dev dependencies..."
pip install --upgrade pip
pip install -e .[dev]

echo "[setup] Done. Run the pipeline with:"
echo "  python -m stent_hemo pipeline all"
echo "  pytest tests -q"