#!/usr/bin/env bash
# Run every generated OpenFOAM case (wedge axisymmetric) and extract
# the area-averaged inlet-to-outlet pressure drop into
# results/openfoam_results.csv for the postprocess cross-check table.
set -euo pipefail
cd /work

mkdir -p results
out=results/openfoam_results.csv
echo "case,dp_pa" > "$out"

source /opt/openfoam11/etc/bashrc 2>/dev/null || true

for case in openfoam/stenosis_*; do
    [ -d "$case" ] || continue
    (cd "$case" && bash Allrun)
    dp=$(python3 - "$case" <<'PY'
import re, sys
from pathlib import Path
case = sys.argv[1]
def last(field):
    data = (Path(case) / "postProcessing" / "surfaces" / "1"
            / "surface.dat").read_text()
    rows = [l.split() for l in data.splitlines()
            if l and not l.startswith("#")]
    return float(rows[-1][1]) if rows else float("nan")
print(f"{abs(last('p')):.2f}")
PY
)
    echo "$(basename "$case"),$dp" >> "$out"
    echo "[of] $(basename "$case") dp = $dp Pa"
done
