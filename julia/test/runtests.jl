# Entry point for `julia --project=. -e 'using Pkg; Pkg.test()'` or
# `julia --project=. test/runtests.jl`.
#
# Runs only the in-Julia property tests (fast, no Octave needed).
# For full parity vs the original MATLAB SPM code, run scripts/parity.sh
# which additionally invokes Octave on the unmodified spm_*.m files.

include("properties/test_distributions.jl")
