#!/usr/bin/env bash
set -euo pipefail
mkdir -p /logs/verifier
python /tests/grade.py /app/output "/tests/reference" > /logs/verifier/reward.txt
