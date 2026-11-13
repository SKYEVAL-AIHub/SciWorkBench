#!/bin/bash
# Harbor verifier entrypoint. Runs from the working directory (typically /app).
# Produces the reward file(s) under /logs/verifier/ as required by Harbor.
set -uo pipefail

PY="$(command -v python3 || command -v python)"
SUBMISSION_DIR="${SUBMISSION_DIR:-/app/output}"
LOGS_DIR="${LOGS_DIR:-/logs/verifier}"
GRADER="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/grader.py"

mkdir -p "$LOGS_DIR"

"$PY" "$GRADER" --submission-dir "$SUBMISSION_DIR" --logs-dir "$LOGS_DIR"
rc=$?

if [ "$rc" -ne 0 ] || [ ! -f "$LOGS_DIR/reward.txt" ]; then
  echo "[test.sh] grader exited with code $rc; writing zero reward." >&2
  echo 0 > "$LOGS_DIR/reward.txt"
  echo '{"reward": 0.0}' > "$LOGS_DIR/reward.json"
fi

exit 0
