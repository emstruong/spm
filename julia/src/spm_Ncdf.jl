"""
    spm_Ncdf(x, u=0, v=1)

Cumulative Distribution Function (CDF) for univariate Normal distributions.

    FORMAT F = spm_Ncdf(x, u, v)

    x - ordinates
    u - mean              [Defaults to 0]
    v - variance  (v>0)   [Defaults to 1]
    F - CDF of N(u,v) at x (lower-tail probability)

Algorithm: Phi(z) = 0.5 + erf(z/sqrt(2))/2  (Abramowitz & Stegun, 26.2.29)

For v <= 0 the result is NaN at those entries.

Translated from spm_Ncdf.m
(Andrew Holmes; (C) 1995-2022 Wellcome Centre for Human Neuroimaging).
"""
function spm_Ncdf(x, u=0, v=1)
    _check_compatible_sizes(x, u, v)
    return _ncdf.(x, u, v)
end

function _ncdf(x, u, v)
    v > 0 || return oftype(float(x - u), NaN)
    return 0.5 + 0.5 * erf((x - u) / sqrt(2 * v))
end
