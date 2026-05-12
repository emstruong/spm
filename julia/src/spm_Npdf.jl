"""
    spm_Npdf(x, u=0, v=1)

Probability Density Function (PDF) of univariate Normal distribution.

    FORMAT f = spm_Npdf(x, u, v)

    x - ordinates
    u - mean              [Defaults to 0]
    v - variance  (v>0)   [Defaults to 1]
    f - pdf of N(u,v) at x

Let random variable X have a Normal distribution with mean u and variance v,
then X~N(u,v). The PDF is

    f(x) = (1 / sqrt(v*2*pi)) * exp( -(x-u)^2 / (2v) )

For v <= 0 the result is NaN at those entries.

Reference: Evans, Hastings, Peacock (1993) "Statistical Distributions", Ch29;
Abramowitz & Stegun (1964). Translated from spm_Npdf.m
(Andrew Holmes; (C) 1994-2022 Wellcome Centre for Human Neuroimaging).
"""
function spm_Npdf(x, u=0, v=1)
    _check_compatible_sizes(x, u, v)
    return _npdf.(x, u, v)
end

# scalar kernel; broadcast does the rest
function _npdf(x, u, v)
    v > 0 || return oftype(float(x - u), NaN)
    z = x - u
    return exp(-(z*z) / (2*v)) / sqrt(2*pi*v)
end

"""
    _check_compatible_sizes(args...)

Mirror the SPM size-check: among non-scalar arguments, all sizes must match
exactly. (SPM does not use NumPy-style singleton broadcasting between
non-scalars; scalars do broadcast.)
"""
function _check_compatible_sizes(args...)
    sizes = Tuple{Vararg{Int}}[size(a) for a in args if !(a isa Number) && length(a) > 1]
    isempty(sizes) && return nothing
    ref = sizes[1]
    for s in sizes
        s == ref || error("non-scalar args must match in size")
    end
    return nothing
end
