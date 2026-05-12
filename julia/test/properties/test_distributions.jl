# Property tests for translated distributions. These catch bugs that the
# Octave-comparison harness can miss when both implementations are wrong the
# same way (e.g. a shared off-by-half-ulp normalization).

using Test
push!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "..", "src")))
using SPM

@testset "spm_Npdf properties" begin
    # symmetry: f(x|0,1) = f(-x|0,1)
    xs = range(-6, 6; length=121)
    @test all(SPM.spm_Npdf.(xs) .== SPM.spm_Npdf.(-xs))
    # peak at mean
    @test SPM.spm_Npdf(0.0) ≈ 1 / sqrt(2π)
    # default args
    @test SPM.spm_Npdf(0.0) == SPM.spm_Npdf(0.0, 0.0, 1.0)
    # NaN for v <= 0
    @test isnan(SPM.spm_Npdf(0.0, 0.0, 0.0))
    @test isnan(SPM.spm_Npdf(0.0, 0.0, -1.0))
end

@testset "spm_Ncdf properties" begin
    # monotone non-decreasing
    xs = sort!(randn(500))
    Fs = SPM.spm_Ncdf.(xs)
    @test all(diff(Fs) .>= -eps())
    # complement at zero: F(0) = 0.5
    @test SPM.spm_Ncdf(0.0) == 0.5
    # symmetry: F(-x) = 1 - F(x)
    xs = range(-5, 5; length=51)
    @test all(SPM.spm_Ncdf.(xs) .+ SPM.spm_Ncdf.(-xs) .≈ 1.0)
    # limits
    @test SPM.spm_Ncdf(-30.0) == 0.0
    @test SPM.spm_Ncdf(30.0) == 1.0
    # NaN for v <= 0
    @test isnan(SPM.spm_Ncdf(0.0, 0.0, -1.0))
end

@testset "spm_Tcdf properties" begin
    # monotone non-decreasing in x (for fixed v)
    xs = sort!(randn(500) .* 3)
    Fs = SPM.spm_Tcdf.(xs, 5.0)
    @test all(diff(Fs) .>= -eps())
    # F(0|v) = 0.5 exactly
    for v in (0.5, 1.0, 2.0, 30.0, 1e6)
        @test SPM.spm_Tcdf(0.0, v) == 0.5
    end
    # symmetry: F(-x|v) = 1 - F(x|v)
    xs = range(-5, 5; length=51)
    for v in (1.5, 4.0, 30.0)
        @test all(SPM.spm_Tcdf.(xs, v) .+ SPM.spm_Tcdf.(-xs, v) .≈ 1.0)
    end
    # v=1 Cauchy closed form
    @test SPM.spm_Tcdf(1.0, 1.0) ≈ 0.75
    @test SPM.spm_Tcdf(-1.0, 1.0) ≈ 0.25
    # large v converges to Ncdf
    for x in (-2.0, -0.5, 0.5, 2.0)
        @test isapprox(SPM.spm_Tcdf(x, 1e7), SPM.spm_Ncdf(x); atol=1e-6)
    end
    # NaN for v <= 0
    @test isnan(SPM.spm_Tcdf(1.0, 0.0))
    @test isnan(SPM.spm_Tcdf(1.0, -1.0))
end
