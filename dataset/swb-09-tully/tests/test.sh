#!/bin/bash
set -e

# Verification script for the Tully 2012 RPSH Fig. 2 rate-reproduction task.
# Runs the deterministic 8-checkpoint grader against the agent submission and
# writes Harbor's scalar reward and a separate detailed grading report.

WORKDIR="$(pwd)"

# Submission directory (agent output)
if [ -n "${TASK_OUTPUT_DIR:-}" ] && [ -d "${TASK_OUTPUT_DIR}" ]; then
    SUBMISSION_DIR="${TASK_OUTPUT_DIR}"
elif [ -d "${WORKDIR}/output" ]; then
    SUBMISSION_DIR="${WORKDIR}/output"
elif [ -d "/workspace/output" ]; then
    SUBMISSION_DIR="/workspace/output"
else
    SUBMISSION_DIR="${WORKDIR}/output"
fi

# Log directory (Harbor container: /logs/verifier; local: ./logs/verifier)
if [ -d "/logs" ]; then
    LOG_DIR="/logs/verifier"
else
    LOG_DIR="${WORKDIR}/logs/verifier"
fi
mkdir -p "$LOG_DIR"

# Tests directory
if [ -d "/tests" ] && [ -f "/tests/grader.py" ]; then
    TESTS_DIR="/tests"
else
    TESTS_DIR="$(cd "$(dirname "$0")" && pwd)"
fi

echo "=== Tully 2012 RPSH Fig. 2 Verification ==="
echo "Submission: ${SUBMISSION_DIR}"
echo "Tests:      ${TESTS_DIR}"
echo "Log:        ${LOG_DIR}"
echo ""

if [ ! -d "${SUBMISSION_DIR}" ]; then
    echo "FAIL: submission directory not found"
    cat > "${LOG_DIR}/grader_report.json" <<EOF
{"reward": 0.0, "score": 0, "max_score": 100, "message": "submission directory missing"}
EOF
    printf '0\n' > "$LOG_DIR/reward.txt"
    exit 0
fi

set +e
python3 "${TESTS_DIR}/grader.py" \
    --submission-dir "${SUBMISSION_DIR}" \
    --output-json "${LOG_DIR}/grader_report.json" 2>&1 | tee "${LOG_DIR}/grader_output.txt"
GRADER_EXIT=${PIPESTATUS[0]}
set -e

echo ""
echo "Grader exit code: ${GRADER_EXIT}"

if [ ! -f "${LOG_DIR}/grader_report.json" ]; then
    cat > "${LOG_DIR}/grader_report.json" <<EOF
{"reward": 0.0, "score": 0, "max_score": 100, "message": "grader produced no report"}
EOF
fi

echo ""
echo "=== grader_report.json ==="
cat "${LOG_DIR}/grader_report.json"
python3 - "$LOG_DIR" <<'PY'
import json
from pathlib import Path
import sys

logs = Path(sys.argv[1])
reward = float(json.loads((logs / "grader_report.json").read_text())["reward"])
if not 0 <= reward <= 1:
    raise ValueError(f"Invalid Tully reward: {reward}")
(logs / "reward.txt").write_text(f"{reward}\n")
PY
exit 0
