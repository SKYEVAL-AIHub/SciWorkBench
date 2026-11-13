#!/bin/bash
set -euo pipefail
umask 077
logs=/logs/verifier
mkdir -p "$logs"
chmod 700 "$logs"
if [ ! -d /workspace/output ]; then
    printf '0\n' > "$logs/reward.txt"
    exit 0
fi
python /tests/evaluate.py --data /workspace/data --submission /workspace/output \
    --rubric /tests/rubric.json --json-out "$logs/grader-report.json"
python - "$logs" <<'PY'
import json
from pathlib import Path
import sys

logs = Path(sys.argv[1])
report = json.loads((logs / "grader-report.json").read_text())
reward = float(report["score_final"])
if not 0 <= reward <= 1:
    raise ValueError(f"Invalid TESS reward: {reward}")
(logs / "reward.txt").write_text(f"{reward}\n")
PY
