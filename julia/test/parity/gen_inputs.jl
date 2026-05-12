# Generate a fixed input grid for parity testing. Writes one .mat per fn:
#   data/inputs_<fn>.mat with variables: name, args (cell {x[,u,v]})
#
# Inputs cover: random (uniform & log-uniform), adversarial scalars (0, ±Inf,
# NaN, domain boundaries), scalar-vs-vector broadcasting cases, and the
# algorithm-switch region of betainc for spm_Tcdf.

using MAT
using Random
using Printf

const HERE = @__DIR__
const DATA = joinpath(HERE, "data")
mkpath(DATA)

Random.seed!(20260512)

# ---- helpers --------------------------------------------------------------

# Concatenate a list of arrays/scalars into a single column vector of Float64,
# treating scalars as 1-element.
flat(xs) = collect(Float64, Iterators.flatten(((x isa AbstractArray) ? vec(x) : (x,)) for x in xs))

function adversarial_scalars()
    return [0.0, -0.0, 1.0, -1.0, 1e-300, -1e-300, 1e300, -1e300,
            Inf, -Inf, NaN, eps(), -eps(), prevfloat(1.0), nextfloat(0.0)]
end

# ---- input grids ----------------------------------------------------------

# spm_Npdf(x, u, v): x ∈ ℝ, u ∈ ℝ, v > 0  (v ≤ 0 → NaN, also tested)
function inputs_Npdf()
    xs  = flat([randn(2000), range(-10, 10; length=201), adversarial_scalars()])
    us  = zeros(length(xs))
    vs  = ones(length(xs))
    cases = Tuple[]
    push!(cases, (xs, us, vs))
    # Non-unit variance, non-zero mean
    push!(cases, (flat([randn(500) .* 3 .+ 2]), fill(2.0, 500), fill(9.0, 500)))
    # v <= 0 path (expect NaN)
    push!(cases, ([1.0, 2.0, 3.0], [0.0, 0.0, 0.0], [0.0, -1.0, -1e-12]))
    # Scalar broadcasting (scalar u, scalar v)
    push!(cases, (collect(range(-5, 5; length=11)), 0.0, 1.0))
    return cases
end

# spm_Ncdf(x, u, v): same shape as Npdf
function inputs_Ncdf()
    xs  = flat([randn(2000), range(-10, 10; length=201), adversarial_scalars()])
    cases = Tuple[]
    push!(cases, (xs, zeros(length(xs)), ones(length(xs))))
    push!(cases, (flat([randn(500) .* 3 .+ 2]), fill(2.0, 500), fill(9.0, 500)))
    push!(cases, ([1.0, 2.0, 3.0], [0.0, 0.0, 0.0], [0.0, -1.0, -1e-12]))
    push!(cases, (collect(range(-5, 5; length=11)), 0.0, 1.0))
    # Far-tail values (Worsley note in docstring): |z| > 6
    push!(cases, ([-10.0, -8.0, -7.0, 7.0, 8.0, 10.0, 37.0, -37.0],
                  zeros(8), ones(8)))
    return cases
end

# spm_Tcdf(x, v): x ∈ ℝ, v > 0
function inputs_Tcdf()
    cases = Tuple[]
    # Bulk: random t-variates across many df
    rng_x = flat([randn(2000) .* 3, range(-20, 20; length=401), adversarial_scalars()])
    rng_v = vcat(fill(1.0, length(rng_x)÷5),      # Cauchy special case
                 fill(2.0, length(rng_x)÷5),
                 fill(30.0, length(rng_x)÷5),
                 fill(0.5, length(rng_x)÷5),      # non-integer df
                 fill(1e6, length(rng_x) - 4*(length(rng_x)÷5)))  # ~normal
    @assert length(rng_x) == length(rng_v)
    push!(cases, (rng_x, rng_v))
    # x = 0 special case (expect exactly 0.5)
    push!(cases, ([0.0, 0.0, 0.0, 0.0], [1.0, 2.0, 10.0, 100.0]))
    # v <= 0 (expect NaN)
    push!(cases, ([1.0, -1.0, 0.0], [0.0, -1.0, -1e-6]))
    # Algorithm-switch region for betainc: z = v/(v+x^2) near 0.5
    # corresponds to x^2 ≈ v, e.g. x = sqrt(v)
    vs = [2.0, 3.0, 5.0, 10.0, 30.0, 100.0]
    xs = vcat([sqrt(v) for v in vs], [-sqrt(v) for v in vs],
              [sqrt(v) * 0.99 for v in vs], [sqrt(v) * 1.01 for v in vs])
    push!(cases, (xs, repeat(vs, 4)))
    return cases
end

# ---- write -----------------------------------------------------------------

function write_cases(name::AbstractString, cases)
    file = joinpath(DATA, "inputs_$(name).mat")
    # Flat schema to avoid MAT.jl/Octave cell-array interop issues:
    # for case i (1-indexed) and arg position k (1..nargs), variable is
    #   "case<i>_arg<k>"  : a column vector of Float64
    # plus:
    #   "ncases"          : scalar
    #   "nargs"           : scalar (assumed constant across cases for one fn)
    #   "name"            : string
    d = Dict{String,Any}()
    d["name"]   = name
    d["ncases"] = Float64(length(cases))
    d["nargs"]  = Float64(length(cases[1]))
    for (i, c) in enumerate(cases)
        for (k, a) in enumerate(c)
            v = a isa AbstractArray ? collect(Float64, vec(a)) : Float64[a]
            d["case$(i)_arg$(k)"] = v
        end
    end
    matwrite(file, d)
    @printf("wrote %s with %d cases (nargs=%d)\n", file, length(cases), length(cases[1]))
end

write_cases("Npdf", inputs_Npdf())
write_cases("Ncdf", inputs_Ncdf())
write_cases("Tcdf", inputs_Tcdf())
