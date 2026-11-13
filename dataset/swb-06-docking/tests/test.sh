#!/usr/bin/env bash
set -u
TEST_DIR="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
mkdir -p /logs/verifier
python3 "$TEST_DIR/grader.py" --output /app/output --oracle "$TEST_DIR/oracle" --report /logs/verifier/report.json --reward /logs/verifier/reward.txt
