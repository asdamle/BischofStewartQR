# The headline guarantee: ||R11^{-1}||_F stays bounded (orthonormal-row
# input, the GKS setting). Port of matlab_rand/tests/test_bsqr_rand_bounds.m.

import BSPivotQR

@testset "bsqr_rand running_mean bound" begin
    for (c, (k, n)) in enumerate([(8, 50), (16, 120), (32, 300), (24, 400)])
        M = orthonormal_rows(k, n, 200 + c)
        F = bsqr_rand(M; seed = c, threshold_mode = :running_mean)
        verify_bound(F.stats, k, n)
    end
end

@testset "bsqr_rand worstcase_allowance bound" begin
    for (c, (k, n)) in enumerate([(8, 50), (16, 120), (32, 300), (24, 400)])
        M = orthonormal_rows(k, n, 300 + c)
        F = bsqr_rand(M; seed = c, threshold_mode = :worstcase_allowance)
        verify_bound(F.stats, k, n)
    end
end

@testset "bsqr_rand f2 matches explicit ||R11^{-1}||_F^2" begin
    k, n = 20, 140
    M = orthonormal_rows(k, n, 55)
    F = bsqr_rand(M; seed = 6)
    explicit = sum(abs2, 1.0 ./ svdvals(F.R11))
    @test abs(F.stats.f2[end] - explicit) / explicit < 1e-8
    @test F.stats.frob_inv ≈ sqrt(F.stats.f2[end])
    @test F.stats.osinsky_bound ≈ sqrt(k * (n - k + 1))
end

@testset "bsqr_rand samples reported" begin
    k, n = 16, 200
    M = orthonormal_rows(k, n, 71)
    # Single-select: every step samples at least one block (>= 1 column tested).
    F1 = bsqr_rand(M; batched = false, seed = 8)
    @test F1.stats.total_tested >= k
    @test all(F1.stats.samples_tested .>= 1)
    # Batched (default): the per-block apply is attributed to the block's
    # first pick, so total_tested (not every step) is the meaningful count.
    F2 = bsqr_rand(M; seed = 8)
    @test F2.stats.total_tested >= k
    @test F2.stats.blocks_sampled >= 1
end

@testset "bsqr_rand comparable to deterministic BSQR" begin
    # Conditioning should be in the same ballpark as the deterministic kernel
    # (not identical pivots). BSPivotQR is a test-only dependency here -- the
    # package kernels stay decoupled.
    k, n = 20, 160
    M = orthonormal_rows(k, n, 88)
    F = bsqr_rand(M; seed = 2)
    rand_frobinv = F.stats.frob_inv

    Fdet = BSPivotQR.bsqr(M; k = k)
    Rdet = BSPivotQR.R(Fdet)[1:k, 1:k]
    det_frobinv = sqrt(sum(abs2, 1.0 ./ svdvals(Rdet)))
    # Randomized should be within a small constant factor of deterministic.
    @test rand_frobinv < 3 * det_frobinv + 1e-8
    # Always: respects the absolute Osinsky guarantee.
    @test rand_frobinv <= sqrt(k * (n - k + 1)) * (1 + 1e-8)
end
