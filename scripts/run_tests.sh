#!/usr/bin/env bash
# Run the full test suite for Coronary-Stent-CFD.
# Usage: ./scripts/run_tests.sh

set -euo pipefail

echo "[tests] Running lint..."
ruff check src tests

echo "[tests] Running unit tests..."
pytest tests -q

echo "[tests] Running quick validation suite..."
stent-hemo validate --quick

echo "[tests] All checks passed."