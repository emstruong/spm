"""
    spm_Tcdf(x, v)

Cumulative Distribution Function (CDF) of Student's t-distribution.

    FORMAT F = spm_Tcdf(x, v)

    x - t-variate (range -Inf..Inf)
    v - degrees of freedom (v>0; non-integer accepted)
    F - CDF of Student's t with v df at points x

Algorithm (Abramowitz & Stegun 26.5.27, 26.7.1):
    Pr(|X| < x) = betainc(v/(v+x^2), v/2, 1/2)
    F(x) = (x>0) - sign(x) * 0.5 * betainc(v/(v+x^2), v/2, 1/2)
With special cases:
    F(0)  = 0.5
    v = 1 = Cauchy: F(x) = 0.5 + atan(x)/pi
    v <= 0 -> NaN

Note: MATLAB's `betainc(X, A, B)` is the lower regularized incomplete beta
I_X(A, B). In Julia we use `SpecialFunctions.beta_inc(A, B, X)[1]`.

Translated from spm_Tcdf.m
(Andrew Holmes; (C) 1992-2022 Wellcome Centre for Human Neuroimaging).
"""
function spm_Tcdf(x, v)
    _check_compatible_sizes(x, v)
    return _tcdf.(x, v)
end

function _tcdf(x, v)
    v > 0 || return oftype(float(x), NaN)
    x == 0 && return oftype(float(x), 0.5)
    v == 1 && return 0.5 + atan(x) / pi
    # Iₓ(a,b) with x = v/(v+x^2), a = v/2, b = 1/2
    z = v / (v + x*x)
    Ix = beta_inc(v/2, 1/2, z)[1]
    return x > 0 ? 1 - 0.5 * Ix : 0.5 * Ix
end
