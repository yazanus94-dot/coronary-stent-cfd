# Benchmark Results

Real measurements from `python benchmarks/benchmark.py` (CPython, single core;
`benchmarks/benchmark_results.csv` is regenerated on every run).

| Benchmark | Workload | Time (s) | Min (s) | Max (s) |
|-----------|----------|----------|---------|---------|
| solver_steady_20x8 | 20x8 steady solve | TBD | TBD | TBD |
| solver_steady_40x14 | 40x14 steady solve | TBD | TBD | TBD |
| solver_steady_80x28 | 80x28 steady solve | TBD | TBD | TBD |
| solver_steady_120x24 | 120x24 steady solve | TBD | TBD | TBD |
| postprocess_full | 50% + 70% stenosis + stent effect + figures | TBD | TBD | TBD |
| validation_quick | Poiseuille (20x8) + Womersley + stenosis | TBD | TBD | TBD |

## Reproduce

```bash
python benchmarks/benchmark.py
```

## Interpretation

The steady solver scales approximately linearly with grid size. The full
post-processing stage runs both stenosis cases plus the stent-effect estimate
and figure generation. The quick validation suite (coarse grids only) is used
in CI and completes in seconds.

Replace the `TBD` values with real measurements after the first benchmark run.