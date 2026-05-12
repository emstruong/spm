#!/usr/bin/env bash
# One-button parity check:
#   1. Julia generates input grids -> data/inputs_<fn>.mat
#   2. Octave runs original spm_<fn>.m on those inputs -> data/reference_<fn>.mat
#   3. Julia runs translations and diffs against reference.
#
# Exit code 0 iff all functions pass within their declared tolerance.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
JULIA_DIR="$(dirname "$HERE")"
PARITY_DIR="$JULIA_DIR/test/parity"
export PATH="/root/.juliaup/bin:$PATH"

cd "$PARITY_DIR"

echo "== [1/3] generating inputs via Julia =="
julia --project="$JULIA_DIR" gen_inputs.jl

echo "== [2/3] running reference via Octave =="
for fn in Npdf Ncdf Tcdf; do
    octave-cli --no-gui --quiet --path "$PARITY_DIR" \
        --eval "run_reference('$fn')"
done

echo "== [3/3] checking parity via Julia =="
julia --project="$JULIA_DIR" check_parity.jl
