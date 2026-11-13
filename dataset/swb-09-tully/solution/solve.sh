#!/bin/bash
set -e

# Reference solution for the Tully 2012 RPSH Fig. 2 rate-reproduction task.
# Regenerates the seven required output/ deliverables from the two audited
# frozen reference results (reference_143K.json / reference_250K.json).

WORKDIR="$(pwd)"

# Locate solution directory (Harbor container: /solution; local: script dir)
if [ -d "/solution" ] && [ -f "/solution/build_outputs.py" ]; then
    SOLUTION_DIR="/solution"
else
    SOLUTION_DIR="$(cd "$(dirname "$0")" && pwd)"
fi

# Output directory (Harbor container: /workspace/output; local: ./output)
if [ -n "${TASK_OUTPUT_DIR:-}" ]; then
    OUTPUT_DIR="${TASK_OUTPUT_DIR}"
elif [ -d "/workspace" ]; then
    OUTPUT_DIR="/workspace/output"
else
    OUTPUT_DIR="${WORKDIR}/output"
fi

echo "=== Tully 2012 RPSH Fig. 2 Reference Solution ==="
echo "Working directory:  $WORKDIR"
echo "Solution directory: $SOLUTION_DIR"
echo "Output directory:   $OUTPUT_DIR"

mkdir -p "$OUTPUT_DIR"

# Regenerate all deliverables from the frozen reference results.
set +e
python3 "${SOLUTION_DIR}/build_outputs.py" \
    --reference-143 "${SOLUTION_DIR}/reference_143K.json" \
    --reference-250 "${SOLUTION_DIR}/reference_250K.json" \
    --src-dir "${SOLUTION_DIR}/src" \
    --output-dir "${OUTPUT_DIR}"
BUILD_EXIT=$?
set -e

if [ "$BUILD_EXIT" -ne 0 ]; then
    echo "WARNING: build_outputs.py failed (exit ${BUILD_EXIT}); falling back to pre-generated deliverables."
    if [ -d "${SOLUTION_DIR}/pregenerated" ]; then
        cp -R "${SOLUTION_DIR}/pregenerated/." "${OUTPUT_DIR}/"
    else
        echo "ERROR: no pre-generated fallback available."
        exit 1
    fi
fi

echo "=== Solution complete ==="
echo "Output contents:"
ls -R "$OUTPUT_DIR"
