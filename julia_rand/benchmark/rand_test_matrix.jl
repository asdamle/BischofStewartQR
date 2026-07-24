# Port of matlab_rand/benchmark/rand_test_matrix.m: k-by-n matrices with
# orthonormal rows (M*M' = I_k) that stress the randomized column-selection
# scheme in different ways. The sampler's cost is governed by how many
# remaining columns are "acceptable" at each step; the leverage scores
# ell_j = ||M[:,j]||^2 (summing to k) are exactly the column quantities the
# sampler sees first, and the families below span uniform to highly
# concentrated leverage. Same family definitions as the MATLAB generator, but
# a different RNG -- matrices are reproducible per (family, k, n, seed) within
# Julia, not across languages.
#
# Included by both the test suite and the benchmark runner (the Julia analogue
# of matlab_rand tests addpath'ing matlab_rand/benchmark).

using LinearAlgebra
using Random

"""
    rand_test_matrix(family, k, n, seed = 0) -> Matrix{Float64}

`k×n` matrix with orthonormal rows. Families: `:gaussian` (near-uniform
leverage, benign), `:graded_leverage` (smooth decay), `:spiked_leverage` (a
few dominant columns), `:coherent` (columns cluster into a few groups),
`:collinear_cluster` (one high-leverage near-collinear cluster — an
R11-conditioning stress case), `:needle` (~k useful columns among near-null
ones; hardest for uniform sampling), `:chebyshev` (polynomial design matrix;
leverage concentrates near the ends).
"""
function rand_test_matrix(family::Symbol, k::Int, n::Int, seed::Integer = 0)
    rng = Xoshiro(seed)
    W = if family === :gaussian
        randn(rng, n, k)
    elseif family === :graded_leverage
        # Smoothly decaying row scaling -> graded leverage across columns.
        d = exp10.(range(0.0, -3.0; length = n))     # 1 ... 1e-3
        d .* randn(rng, n, k)
    elseif family === :spiked_leverage
        # A few dominant columns carry most of the leverage.
        nspike = max(1, ceil(Int, 1.5 * k))
        d = ones(n)
        d[randperm(rng, n)[1:min(nspike, n)]] .= 100.0
        d .* randn(rng, n, k)
    elseif family === :coherent
        # Rows of the factor cluster into r groups -> most columns of M are
        # near-duplicates of a few representatives.
        r = max(2, ceil(Int, k / 4))
        centers = randn(rng, r, k)
        assign = rand(rng, 1:r, n)
        centers[assign, :] .+ 0.05 .* randn(rng, n, k)
    elseif family === :collinear_cluster
        # One high-leverage, near-collinear cluster (all close to a single
        # direction v0) plus well-spread background columns: selecting more
        # than one cluster column makes the k-by-k R11 near-singular.
        ncl = max(2, round(Int, 0.20 * n))
        v0 = randn(rng, 1, k)
        Wc = randn(rng, n, k)
        Wc[1:ncl, :] .= 30.0 .* (v0 .+ 1e-3 .* randn(rng, ncl, k))
        Wc
    elseif family === :needle
        # ~k useful columns; the rest are near-null (tiny leverage).
        ng = min(n, k + ceil(Int, 0.25 * k))
        d = fill(1e-6, n)
        d[randperm(rng, n)[1:ng]] .= 1.0
        d .* randn(rng, n, k)
    elseif family === :chebyshev
        # Degree-(k-1) Chebyshev design matrix on n random points in [-1,1];
        # leverage (the Christoffel function) concentrates near the ends.
        x = sort!(2.0 .* rand(rng, n) .- 1.0)
        B = ones(k, n)
        if k >= 2
            B[2, :] .= x
        end
        for i in 3:k
            B[i, :] .= 2.0 .* x .* B[i - 1, :] .- B[i - 2, :]
        end
        Matrix(B')
    else
        throw(ArgumentError("unknown rand_test_matrix family $(family)"))
    end
    # Householder QR always produces an exactly orthonormal thin Q, even for
    # (near-)rank-deficient W, so no orth()+pad dance is needed here (the
    # MATLAB generator pads for degenerate orth() output).
    Qm = Matrix(qr(W).Q)          # n×k, orthonormal columns
    return Matrix(Qm')            # k×n, orthonormal rows
end
