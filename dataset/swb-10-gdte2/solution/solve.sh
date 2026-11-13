#!/bin/bash
# Oracle solution: regenerate every deliverable from the raw data in /app/data.
# Nothing is hard-coded; the reference pipeline recomputes all vectors, the
# boundary curve, the bias-anomaly scores, the Fourier feature, the point-spectrum
# edges and the forward/retrace validation straight from the NPZ/CSV inputs, then
# build_summary.py maps those results into the agent-facing summary.json schema.
set -euo pipefail

APP_DIR="${APP_DIR:-/app}"
DATA_DIR="${TASK_DATA_DIR:-$APP_DIR/data}"
OUT_DIR="${TASK_OUTPUT_DIR:-$APP_DIR/output}"
SOLUTION_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PY="$(command -v python3 || command -v python)"

mkdir -p "$OUT_DIR"

# 1) Run the deterministic reference pipeline. Emits six *_reference.png figures
#    plus reference_results.json / robustness_summary.json / reference_results.csv.
"$PY" "$SOLUTION_DIR/reference_analysis.py" \
    --data-dir "$DATA_DIR" \
    --output-dir "$OUT_DIR"

# 2) Map reference_results.json -> summary.json (agent schema) and rename the
#    six figures to the required deliverable names.
"$PY" "$SOLUTION_DIR/build_summary.py" \
    --reference-json "$OUT_DIR/reference_results.json" \
    --output-dir "$OUT_DIR"

# 3) The scientific report and the reproducible script are required deliverables.
cp "$SOLUTION_DIR/analysis_report.md" "$OUT_DIR/analysis_report.md"
cp "$SOLUTION_DIR/reference_analysis.py" "$OUT_DIR/analysis.py"

echo "Solution outputs written to $OUT_DIR:"
ls -l "$OUT_DIR"
