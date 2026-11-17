# Benchmark Results

Real measurements from `python benchmarks/benchmark.py` (CPython 3.14,
Windows laptop, single core; `benchmarks/benchmark_results.csv` is
regenerated on every run).

| Benchmark | Workload | Time (s) | Min (s) | Max (s) |
|-----------|----------|----------|---------|---------|
| solver_steady_20x8 | 20x8 steady solve | 0.282 | 0.279 | 0.284 |
| solver_steady_40x14 | 40x14 steady solve | 1.728 | 1.692 | 1.829 |
| solver_steady_80x28 | 80x28 steady solve | 2.553 | 2.505 | 2.647 |
| solver_steady_120x24 | 120x24 steady solve | 4.857 | 4.829 | 4.929 |
| postprocess_full | 50% + 70% stenosis + stent effect + figures | 180.839 | 179.565 | 182.441 |
| validation_quick | Poiseuille (20x8) + Womersley + stenosis | 75.606 | 75.606 | 75.606 |

## Reproduce

```bash
python benchmarks/benchmark.py
```

## Interpretation

The steady solver scales close to linearly with cell count (20x8 -> 120x24 is
~24x the cells and ~17x the time; the coarser grids converge in fewer SIMPLE
iterations). The full post-processing stage runs both stenosis cases plus the
stent-effect estimate and figure generation and is dominated by the two
transient (3-cycle, 1200-step) solves. The quick validation suite (coarse
grids only, used in CI) completes in about 75 s, dominated by the transient
Womersley run.
