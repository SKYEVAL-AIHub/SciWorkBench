#!/usr/bin/env bash
set -euo pipefail
mkdir -p /logs/verifier
python /tests/grade.py /app/output /app/data > /logs/verifier/reward.txt
