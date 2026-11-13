#!/bin/bash
set -euo pipefail
exec python /tests/verify.py /app/output /logs/verifier
