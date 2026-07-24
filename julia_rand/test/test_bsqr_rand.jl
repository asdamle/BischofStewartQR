# Correctness tests for the randomized BSQR variant — the Julia port of
# matlab_rand/tests/test_bsqr_rand.m (minus the MEX-vs-mfile backend items:
# Julia has a single kernel). Pivot sequences are RNG-dependent and need not
# match either MATLAB backend; the assertions here are the invariants both
# languages share (exact factorization, per-step bound, stats accounting).

# k-by-n matrix with orthonormal rows (M*M' = I_k), the GKS setting.
orthonormal_rows(k, n, seed) = Matrix(Matrix(qr(randn(Xoshiro(seed), n, k)).Q)')

function check_factorization(M::AbstractMatrix, F::BSQRRandPivoted)
    m, n = size(M)
    k = F.k
    p = perm(F)
    @test sort(p) == collect(1:n)
    @test size(F.R11) == (k, k)
    @test istriu(F.R11)
    Qm = Q(F)
    @test size(Qm) == (m, k)
    k == 0 && return nothing
    sel = p[1:k]
    scale = max(1.0, norm(M[:, sel]))
    @test norm(Qm' * M[:, sel] - F.R11) / scale < 1e-10   # Q'*A(:,sel) = R11
    @test norm(M[:, sel] - Qm * F.R11) / scale < 1e-10    # A(:,sel) = Q*R11
    @test norm(Qm' * Qm - I) < 1e-10                      # Q orthonormal
    return nothing
end

# Per-step running bound and the final Osinsky bound (orthonormal-row input).
function verify_bound(stats::BSQRRandStats, k::Int, n::Int)
    tol = 1e-8 * max(1.0, stats.Fhat[end])
    @test all(stats.f2 .<= stats.Fhat .+ tol)
    @test stats.f2[end] <= k * (n - k + 1) * (1 + 1e-8)
    return nothing
end

@testset "bsqr_rand reconstruction (orthonormal rows)" begin
    for (c, (k, n)) in enumerate([(8, 40), (16, 64), (32, 200), (10, 10)])
        M = orthonormal_rows(k, n, 100 + c)
        F = bsqr_rand(M; seed = c)
        check_factorization(M, F)
    end
end

@testset "bsqr_rand general matrix" begin
    # Theory assumes orthonormal rows, but the factorization must still be exact.
    M = randn(Xoshiro(7), 20, 90)
    F = bsqr_rand(M; k = 12, batched = false, seed = 3)
    check_factorization(M, F)
end

@testset "bsqr_rand return_r12" begin
    k, n = 12, 60
    M = orthonormal_rows(k, n, 11)
    F = bsqr_rand(M; seed = 5, return_r12 = true)
    check_factorization(M, F)
    @test size(F.R12) == (k, n - k)
    p = perm(F)
    @test norm(Q(F)' * M[:, p[(k + 1):n]] - F.R12) < 1e-10
    @test r12(F) == F.R12
    @test r12(F) !== F.R12          # accessor returns a copy

    # R12 is opt-in: without return_r12 the field is nothing.
    F2 = bsqr_rand(M; seed = 5)
    @test F2.R12 === nothing
    @test r12(F2) === nothing

    # Full selection (k = n): R12 is k-by-0.
    Msq = orthonormal_rows(10, 10, 21)
    Fsq = bsqr_rand(Msq; k = 10, seed = 2, return_r12 = true)
    check_factorization(Msq, Fsq)
    @test size(Fsq.R12) == (10, 0)

    # The single-select path builds the same compact-WY store; R12 must be
    # exact off it too.
    Fss = bsqr_rand(M; seed = 7, batched = false, return_r12 = true)
    check_factorization(M, Fss)
    pss = perm(Fss)
    @test norm(Q(Fss)' * M[:, pss[(k + 1):n]] - Fss.R12) < 1e-10

    # k = 0 with return_r12: the 0-by-n early-return branch.
    F0 = bsqr_rand(M; k = 0, return_r12 = true)
    @test size(F0.R12) == (0, n)
end

@testset "bsqr_rand k = 0 edge" begin
    M = orthonormal_rows(5, 20, 1)
    F = bsqr_rand(M; k = 0)
    @test sort(perm(F)) == collect(1:20)
    @test size(F.R11) == (0, 0)
    @test size(Q(F)) == (5, 0)
    @test isempty(selected(F))
    @test F.stats.total_tested == 0
    @test F.stats.blocks_sampled == 0
end

@testset "bsqr_rand k = full (square orthonormal)" begin
    M = Matrix(qr(randn(Xoshiro(21), 15, 15)).Q)
    F = bsqr_rand(M; k = 15, seed = 1)
    check_factorization(M, F)
end

@testset "bsqr_rand determinism with seed/rng" begin
    M = orthonormal_rows(16, 80, 9)
    F1 = bsqr_rand(M; seed = 42)
    F2 = bsqr_rand(M; seed = 42)
    @test perm(F1) == perm(F2)
    @test F1.R11 == F2.R11
    @test F1.stats.samples_tested == F2.stats.samples_tested

    # An explicit rng reproduces too (fresh state per call).
    F3 = bsqr_rand(M; rng = Xoshiro(7))
    F4 = bsqr_rand(M; rng = Xoshiro(7))
    @test perm(F3) == perm(F4)
    @test F3.R11 == F4.R11
end

@testset "bsqr_rand threshold/sampling/pick knobs" begin
    M = orthonormal_rows(20, 120, 13)
    combos = [
        (:running_mean, :uniform, :best_in_block),
        (:worstcase_allowance, :uniform, :first),
        (:running_mean, :normweighted, :best_in_block),
        (:worstcase_allowance, :normweighted, :first),
    ]
    for (c, (mode, samp, pk)) in enumerate(combos)
        F = bsqr_rand(M; seed = c, threshold_mode = mode, sampling = samp, pick = pk)
        check_factorization(M, F)
        verify_bound(F.stats, 20, 120)
        # pick/visiting order only matter on the single-select path (batched
        # always takes the in-block minimizer), so run each combination there
        # too -- this is what actually exercises pick = :first.
        F2 = bsqr_rand(M; seed = c, batched = false,
            threshold_mode = mode, sampling = samp, pick = pk)
        check_factorization(M, F2)
        verify_bound(F2.stats, 20, 120)
    end
end

@testset "bsqr_rand rank guard" begin
    # k beyond the exact rank: only k-1 nonzero columns, so at step k every
    # remaining residual is exactly zero and the guard must fire rather than
    # propagate Inf into f2/R11. Separate code paths on batched/single-select.
    k, n = 6, 40
    M = [randn(Xoshiro(31), k, k - 1) zeros(k, n - (k - 1))]
    @test_throws RankDeficientError bsqr_rand(M; k = k, seed = 1)
    @test_throws RankDeficientError bsqr_rand(M; k = k, batched = false, seed = 1)
end

@testset "bsqr_rand zero-tail householder" begin
    # A selected column whose below-pivot tail is exactly zero exercises the
    # xnorm == 0 path of the reflector build (tau = 0, beta = alpha).
    # block_size = n puts every column in the first block, so the largest-norm
    # zero-tail column is picked first regardless of the RNG.
    m, n = 8, 30
    M = [[5.0; zeros(m - 1)] 0.1 .* randn(Xoshiro(20260715), m, n - 1)]
    F = bsqr_rand(M; k = 4, block_size = n, seed = 1)
    @test perm(F)[1] == 1
    check_factorization(M, F)
end

@testset "bsqr_rand argument validation" begin
    W = orthonormal_rows(6, 20, 4)
    @test_throws ArgumentError bsqr_rand(W; k = 99)
    @test_throws ArgumentError bsqr_rand(W; k = -1)
    @test_throws ArgumentError bsqr_rand(W; threshold_mode = :bogus)
    @test_throws ArgumentError bsqr_rand(W; sampling = :bogus)
    @test_throws ArgumentError bsqr_rand(W; pick = :bogus)
    @test_throws ArgumentError bsqr_rand(W; slack = 0.5)
    @test_throws ArgumentError bsqr_rand(W; slack = Inf)
    @test_throws ArgumentError bsqr_rand(W; norm_recomp_tol = -0.1)
    @test_throws ArgumentError bsqr_rand(W; norm_recomp_tol = 1.5)
    @test_throws ArgumentError bsqr_rand(W; block_size = 0)
    # check_finite defaults to false, but the default norm-weighted sampling
    # still detects non-finite input for free from its precomputed weights ...
    @test_throws ArgumentError bsqr_rand([1.0 Inf; 2.0 3.0])
    # ... while uniform sampling needs the explicit O(mn) scan.
    @test_throws ArgumentError bsqr_rand([1.0 Inf; 2.0 3.0];
        sampling = :uniform, check_finite = true)
end

@testset "bsqr_rand input type conversion" begin
    W = orthonormal_rows(8, 40, 6)
    F = bsqr_rand(Float32.(W); seed = 1)
    check_factorization(Float64.(Float32.(W)), F)
    Ai = [3 1 4 1 5 9 2 6; 2 7 1 8 2 8 1 8]
    Fi = bsqr_rand(Ai; k = 2, seed = 1)
    check_factorization(Float64.(Ai), Fi)
end

@testset "bsqr_rand instrumentation counts all work" begin
    # samples_tested / rounds must count ALL sampling work, including blocks
    # that yield no selection. Two regimes:
    #   (a) no-failure: block >= n means the first block contains every
    #       column, so all k selections come from one apply -- exact counts,
    #       independent of the RNG.
    #   (b) heavy-failure: uniform sampling on needle with a tiny block fails
    #       often; the totals must reflect that work.
    W = orthonormal_rows(12, 150, 7)
    F = bsqr_rand(W; k = 12, seed = 1, block_size = 200)
    @test F.stats.total_tested == 150   # one all-column block counts each column once
    @test F.stats.blocks_sampled == 1   # ... and is one reflector apply
    @test F.stats.samples_tested[1] == 150
    @test all(F.stats.samples_tested[2:end] .== 0)

    Mneedle = rand_test_matrix(:needle, 16, 400, 23)
    Fn = bsqr_rand(Mneedle; k = 16, seed = 5, sampling = :uniform, block_size = 8)
    @test Fn.stats.total_tested > 300   # failed blocks count toward total_tested
    @test Fn.stats.blocks_sampled > 20  # rounds count failed blocks too
end

@testset "bsqr_rand adversarial corners" begin
    # Kernel corners: degenerate block sizes, the forced-fallback path, the
    # final-step worstcase_allowance division (k-i = 1), a bound-voiding
    # slack, and fully-underflowed sampling weights. Each must either factor
    # exactly with the per-step bound intact, or fail with the documented
    # error.
    k, n = 12, 150
    W = orthonormal_rows(k, n, 5)
    needle = rand_test_matrix(:needle, k, n, 11)

    function corner_check(M, kk; kwargs...)
        F = bsqr_rand(M; k = kk, kwargs...)
        check_factorization(M, F)
        @test all(F.stats.f2 .<= F.stats.Fhat .* (1 + 1e-8) .+ 1e-8)
        return F
    end

    corner_check(W, k; seed = 1, block_size = 1, batched = true)
    corner_check(W, k; seed = 1, block_size = 1, batched = false)
    corner_check(W, k; seed = 2, block_size = 10 * n)
    corner_check(needle, k; seed = 3, block_size = 1, sampling = :uniform)
    corner_check(W, 2; seed = 4, threshold_mode = :worstcase_allowance)

    # slack >= 1 voids the bound but must not break exactness.
    Fs = bsqr_rand(W; k = k, seed = 5, slack = 1e3)
    check_factorization(W, Fs)

    # All sampling weights underflow to zero: squared norms leave the
    # documented domain, and the rank guard fires (rho^2 underflows too).
    @test_throws RankDeficientError bsqr_rand(W .* 1e-170; k = k, seed = 6)
end

@testset "bsqr_rand norm_recomp_tol safeguard" begin
    # The batched in-block path downdates running residual norms
    # incrementally; the safeguard recomputes them exactly once they decay
    # past tol * (last exact value). Exercise the recompute-every-step branch
    # (tol = 1), the default, and the pure-downdate branch (tol = 0) on
    # cancellation-prone inputs -- near-duplicate / near-collinear / near-null
    # columns whose residuals collapse during in-block reduction.
    k, n = 20, 240
    bound = k * (n - k + 1)
    for fam in (:coherent, :collinear_cluster, :needle)
        M = rand_test_matrix(fam, k, n, 29)
        for tol in (0.0, sqrt(eps(Float64)), 1.0)
            F = bsqr_rand(M; k = k, seed = 8, batched = true, norm_recomp_tol = tol)
            check_factorization(M, F)
            tolF = 1e-8 * max(1.0, F.stats.Fhat[end])
            @test all(F.stats.f2 .<= F.stats.Fhat .+ tolF)
            # Realized ||R11^{-1}||_F^2 (inverse-free, via singular values).
            @test sum(abs2, 1.0 ./ svdvals(F.R11)) <= bound * (1 + 1e-6)
        end
    end
end

@testset "bsqr_rand pivot validity on concentrated norms" begin
    # Regression guard for the Fenwick weighted sampler: on general input
    # whose squared column norms span many orders of magnitude, remove/restore
    # traffic can drift the tree's running total past its prefix sums; every
    # draw must still resolve to a valid in-pool column (the MATLAB MEX bug
    # fixed in 20f562f produced pivot index n+1 here). Input family matches
    # matlab_rand/tests/test_bsqr_rand_pivot_range.m ('spectrum_needle').
    function spectrum_needle(m, n, seed)
        V = rand_test_matrix(:needle, m, n, seed)
        U = Matrix(qr(randn(Xoshiro(1000 + seed), m, m)).Q)
        return U * (exp10.(range(0.0, -5.0; length = m)) .* V)
    end
    for (m, n) in [(64, 500), (64, 1000)]
        for mat_seed in 0:4
            A = spectrum_needle(m, n, mat_seed)
            for s in 0:4
                F = bsqr_rand(A; k = m, seed = s)
                @test sort(perm(F)) == collect(1:n)
            end
        end
    end
end
