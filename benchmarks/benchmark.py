"""Performance benchmarks for the Coronary-Stent-CFD pipeline.

Measures the in-house solver and post-processing stages on representative
workloads. Run with:

    python benchmarks/benchmark.py

Results are written to benchmarks/benchmark_results.csv.
"""

from __future__ import annotations

import csv
import time
from pathlib import Path

import numpy as np

from stent_hemo.core import BloodProperties, FlowWaveform, VesselGeometry
from stent_hemo.fv_solver import AxisymSolver
from stent_hemo.struct_grid import StructGrid


def bench_solver_steady(nz: int, nr: int, repeats: int = 5) -> dict[str, float]:
    """Benchmark the steady SIMPLE solver."""
    blood = BloodProperties()
    wf = FlowWaveform()
    vessel = VesselGeometry(severity=0.7, length=90e-3)
    grid = StructGrid(vessel, nz, nr)
    solver = AxisymSolver(grid, blood)

    # Warm-up
    solver.solve(wf, steady=True, max_outer=100, tol=1e-4)

    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        solver.solve(wf, steady=True, max_outer=600, tol=1e-6)
        times.append(time.perf_counter() - t0)

    return {
        "benchmark": f"solver_steady_{nz}x{nr}",
        "workload": f"{nz}x{nr} steady solve",
        "time_s": float(np.mean(times)),
        "min_s": float(np.min(times)),
        "max_s": float(np.max(times)),
    }


def bench_postprocess(repeats: int = 3) -> dict[str, float]:
    """Benchmark the full post-processing stage."""
    import json

    from stent_hemo.postprocess import run_postprocess

    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        results_dir = Path(tmp) / "results"
        fig_dir = Path(tmp) / "figures"
        times = []
        for _ in range(repeats):
            t0 = time.perf_counter()
            run_postprocess(results_dir, fig_dir)
            times.append(time.perf_counter() - t0)

    return {
        "benchmark": "postprocess_full",
        "workload": "50% + 70% stenosis + stent effect + figures",
        "time_s": float(np.mean(times)),
        "min_s": float(np.min(times)),
        "max_s": float(np.max(times)),
    }


def bench_validation_quick() -> dict[str, float]:
    """Benchmark the quick validation suite."""
    from stent_hemo.validation import run_validation

    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        results_dir = Path(tmp) / "results"
        t0 = time.perf_counter()
        run_validation(results_dir, quick=True)
        elapsed = time.perf_counter() - t0

    return {
        "benchmark": "validation_quick",
        "workload": "Poiseuille (20x8) + Womersley + stenosis (50%, 70%)",
        "time_s": elapsed,
        "min_s": elapsed,
        "max_s": elapsed,
    }


def main() -> None:
    results: list[dict[str, float]] = []

    print("== Coronary-Stent-CFD benchmarks ==")
    print()

    for nz, nr in ((20, 8), (40, 14), (80, 28), (120, 24)):
        r = bench_solver_steady(nz, nr)
        results.append(r)
        print(f"{r['benchmark']}: {r['time_s']:.3f}s "
              f"(min {r['min_s']:.3f}, max {r['max_s']:.3f})")

    r = bench_postprocess()
    results.append(r)
    print(f"{r['benchmark']}: {r['time_s']:.3f}s "
          f"(min {r['min_s']:.3f}, max {r['max_s']:.3f})")

    r = bench_validation_quick()
    results.append(r)
    print(f"{r['benchmark']}: {r['time_s']:.3f}s")

    # Write CSV
    out = Path(__file__).parent / "benchmark_results.csv"
    fields = ["benchmark", "workload", "time_s", "min_s", "max_s"]
    with open(out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    print(f"\nResults written to {out}")


if __name__ == "__main__":
    main()