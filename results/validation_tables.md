# Validation tables (auto-generated)

| grid | dp (Pa) | dp err % | u_center | uc err % | iters | net outflow (m^3/s) | cell div max (1/s) |
|---|---|---|---|---|---|---|---|
| 20x8 | 58.867 | 7.12 | 0.33468 | 1.43 | 17 | -4.035e-14 | 2.213e-02 |
| 40x14 | 61.243 | 3.37 | 0.33799 | 0.45 | 16 | 5.514e-13 | 2.365e-02 |
| 80x28 | 62.394 | 1.55 | 0.33922 | 0.09 | 19 | 4.765e-13 | 5.541e-03 |

Womersley: steps=1200, uc range -0.0275..0.5167, flow reversal=True, cycles stored=2, phase-matched periodicity=0.000000

| severity | dp solver (Pa) | dp excess (Pa) | dp Young-Tsai (Pa) | ratio |
|---|---|---|---|---|
| 0.5 | 240.90 | 50.77 | 56.93 | 0.892 |
| 0.7 | 401.33 | 211.19 | 219.82 | 0.961 |

GCI (dp (Pa), straight tube): p=1.047, f_exact=63.472 Pa, gci_fine=2.161% (r=2.0), gci_coarse=4.465%, monotonic=True
