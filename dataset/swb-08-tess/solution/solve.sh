#!/bin/bash
set -e

# Solution script for AST-TESS-001 (AU Mic flare-aware periodicity) task.
# Runs the reference analysis pipeline, which generates the complete ten-file
# submission (including a copy of itself) into output/ from public data only.

WORKDIR="$(pwd)"

# Locate data directory (Harbor container: /workspace/data; local: environment/data)
if [ -d "${WORKDIR}/data" ]; then
    DATA_DIR="${WORKDIR}/data"
elif [ -d "${WORKDIR}/environment/data" ]; then
    DATA_DIR="${WORKDIR}/environment/data"
else
    echo "ERROR: data directory not found"
    exit 1
fi

OUTPUT_DIR="${WORKDIR}/output"
SOLUTION_DIR="/solution"

# Fallback for local execution (solve.sh run outside the container)
if [ ! -d "$SOLUTION_DIR" ] || [ ! -f "$SOLUTION_DIR/analysis_pipeline.py" ]; then
    SOLUTION_DIR="$(cd "$(dirname "$0")" && pwd)"
fi

echo "=== Running AST-TESS-001 reference solution ==="
echo "Working directory:  $WORKDIR"
echo "Data directory:     $DATA_DIR"
echo "Solution directory: $SOLUTION_DIR"
echo "Output directory:   $OUTPUT_DIR"

mkdir -p "$OUTPUT_DIR"

# The pipeline is self-contained: it regenerates all ten required output files
# from the public data and copies itself into the output directory.
# MPLCONFIGDIR avoids matplotlib cache-permission issues in read-only homes.
export MPLCONFIGDIR="${OUTPUT_DIR}/.mplconfig"
export PYTHONDONTWRITEBYTECODE=1

python3 "$SOLUTION_DIR/analysis_pipeline.py" \
    --data "$DATA_DIR" \
    --output "$OUTPUT_DIR" \
    --bootstrap-replicates 300

# Remove the transient matplotlib config dir so output/ holds exactly 10 files
rm -rf "${OUTPUT_DIR}/.mplconfig"

echo "=== Generated submission files ==="
ls -1 "$OUTPUT_DIR"

echo "=== Solution complete ==="
