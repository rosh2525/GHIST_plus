#!/bin/sh
set -eu
GHIST_REPRO_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
GHIST_REPRO_ENV="${GHIST_REPRO_ENV:-$GHIST_REPRO_ROOT/../work/ghist_plus_reproduction/.venv}"
GHIST_REPRO_PYTHON="${GHIST_REPRO_PYTHON:-python3}"
"$GHIST_REPRO_PYTHON" -m venv "$GHIST_REPRO_ENV"
"$GHIST_REPRO_ENV/bin/python" -m pip install -r "$GHIST_REPRO_ROOT/reproduction/requirements-analysis.txt"
printf 'Environment ready: %s\n' "$GHIST_REPRO_ENV"
