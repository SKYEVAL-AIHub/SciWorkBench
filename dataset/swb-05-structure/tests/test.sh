#!/usr/bin/env bash
set -euo pipefail
mkdir -p /logs/verifier
if python /tests/reference/verify_submission.py /app --reward-file /logs/verifier/reward.txt; then
    exit 0
fi
test -s /logs/verifier/reward.txt
