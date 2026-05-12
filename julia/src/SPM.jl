module SPM

# SPM.jl — Julia port of selected SPM (Statistical Parametric Mapping) functions.
#
# Scope (this commit): pilot set of pure-numerical helpers translated from the
# original MATLAB sources in the parent SPM repository. File layout mirrors the
# MATLAB layout one-function-per-file with matching filenames.
#
# Out of scope: MEX kernels (src/), the matlabbatch/ system, the spm.m GUI,
# spm_orthviews/, external/ (vendored libraries), compat/, and the @class
# OOP folders. See julia/README.md for the migration-status table.

using SpecialFunctions: erf, beta_inc

include("spm_Npdf.jl")
include("spm_Ncdf.jl")
include("spm_Tcdf.jl")

export spm_Npdf, spm_Ncdf, spm_Tcdf

end # module
