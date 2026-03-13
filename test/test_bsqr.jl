using LinearAlgebra
using Random
using Test

using BSPivotQR

function _residual_and_orthogonality(A::Matrix{Float64}, F::BSQRPivoted)
    Q = BSPivotQR._explicit_q(F)
    T = BSPivotQR._packed_to_qt(F)
    Aperm = A[:, perm(F)]
    resid = norm(Aperm - Q * T) / max(norm(A), eps(Float64))
    orth = norm(I - Q' * Q)
    return resid, orth
end

@testset "Bischof-Stewart CPQR correctness" begin
    rng = MersenneTwister(20260223)
    shapes = ((12, 8), (10, 10), (8, 12))

    for (m, n) in shapes
        A = randn(rng, m, n)
        F = bsqr(A)

        @test F.ksteps == min(m, n)
        @test sort(perm(F)) == collect(1:n)

        resid, orth = _residual_and_orthogonality(A, F)
        tol_resid = 2.5e2 * eps(Float64) * max(m, n)
        tol_orth = 2.5e2 * eps(Float64) * m

        @test resid <= tol_resid
        @test orth <= tol_orth

        Rt = R(F)
        @test norm(Rt - triu(Rt)) <= 1e-12

        Arec = reconstruct(F, A)
        @test norm(Arec - A) / max(norm(A), eps(Float64)) <= tol_resid
    end
end

@testset "Default full factorization and edge cases" begin
    rng = MersenneTwister(17)
    A = randn(rng, 9, 6)

    F0 = bsqr(A; k = 0)
    @test F0.ksteps == 0
    @test length(F0.tau) == 0

    F1 = bsqr(A; k = 1)
    @test F1.ksteps == 1

    F3 = bsqr(A; k = 3)
    @test F3.ksteps == 3

    B = [ones(6) ones(6) randn(rng, 6, 4)]
    Fr = bsqr(B)
    @test Fr.ksteps == min(size(B)...)
    @test all(isfinite, Fr.factors)

    Frstop = bsqr(B; rank_stop = true)
    @test Frstop.ksteps < min(size(B)...)

    C = [ones(5) ones(5) ones(5) ones(5)]
    Ftie = bsqr(C)
    @test perm(Ftie) == collect(1:size(C, 2))
end

@testset "Prescribed k stopping" begin
    rng = MersenneTwister(2027)
    A = randn(rng, 24, 16)
    k = 6
    F = bsqr(A; k = k)
    @test F.ksteps == k

    B = [ones(20) ones(20) randn(rng, 20, 10)]
    Fauto = bsqr(B)
    @test Fauto.ksteps == min(size(B)...)

    Ffixed = bsqr(B; k = 5, rank_stop = false)
    @test Ffixed.ksteps == 5
    @test all(isfinite, Ffixed.factors)
end

@testset "Rank-stop policy" begin
    rng = MersenneTwister(2028)
    B = [ones(30) ones(30) randn(rng, 30, 10)]
    kmax = min(size(B)...)

    Fstop = bsqr(B; k = kmax, rank_stop = true)
    @test Fstop.ksteps < kmax

    Ffull = bsqr(B; k = kmax, rank_stop = false)
    @test Ffull.ksteps == kmax
    @test all(isfinite, Ffull.factors)
    @test all(isfinite, Ffull.tau)
end

@testset "Criterion-consistent pivot sequence" begin
    rng = MersenneTwister(1234)
    A = randn(rng, 7, 6)
    m, n = size(A)
    k = min(m, n)

    Awork = copy(A)
    ws = BSPivotQR.BSWorkspace(m, n, k)
    tau = zeros(Float64, k)
    jpvt = collect(1:n)
    hist = Int[]
    BSPivotQR._bsqr_kernel!(Awork, tau, jpvt, ws, k; pivot_history = hist)

    M = copy(A)
    p = collect(1:n)
    expected = Int[]
    for i in 1:k
        best_j = i
        best_c = Inf

        if i == 1
            for j in i:n
                rho2 = norm(view(M, i:m, j))^2
                c = rho2 > 0.0 ? 1.0 / rho2 : Inf
                if c < best_c
                    best_c = c
                    best_j = j
                end
            end
        else
            R11 = UpperTriangular(view(M, 1:(i - 1), 1:(i - 1)))
            for j in i:n
                rho2 = norm(view(M, i:m, j))^2
                if rho2 <= 0.0
                    c = Inf
                else
                    wj = R11 \ view(M, 1:(i - 1), j)
                    c = (1.0 + dot(wj, wj)) / rho2
                end
                if c < best_c
                    best_c = c
                    best_j = j
                end
            end
        end

        p[i], p[best_j] = p[best_j], p[i]
        M[:, [i, best_j]] = M[:, [best_j, i]]
        push!(expected, p[i])

        if norm(view(M, i:m, i)) <= 1e-12
            break
        end

        tau_i, beta_i = BSPivotQR._householder!(view(M, i:m, i))
        if tau_i != 0.0
            M[i, i] = 1.0
            BSPivotQR._apply_householder_left!(M, i, tau_i, zeros(Float64, n))
        end
        M[i, i] = beta_i
    end

    @test hist == expected
end

@testset "Comparison with Julia qr(ColumnNorm())" begin
    rng = MersenneTwister(99)
    shapes = ((40, 24), (24, 40), (32, 32))

    for (m, n) in shapes
        A = randn(rng, m, n)

        Fbs = bsqr(A)
        r_bs, q_bs = _residual_and_orthogonality(A, Fbs)

        Fq = qr(copy(A), ColumnNorm())
        p = Vector(Fq.p)
        Qq = Matrix(Fq.Q)
        Rq = Matrix(Fq.R)
        r_q = norm(A[:, p] - Qq * Rq) / max(norm(A), eps(Float64))
        q_q = norm(I - Qq' * Qq)

        tol = 5.0e2 * eps(Float64) * max(m, n)
        @test r_bs <= tol
        @test q_bs <= tol
        @test r_q <= tol
        @test q_q <= tol
    end
end

@testset "Norm recompute tolerance knob" begin
    rng = MersenneTwister(2026)
    m, n = 48, 36
    r = min(m, n)
    U = Matrix(qr(randn(rng, m, r)).Q)
    V = Matrix(qr(randn(rng, n, r)).Q)
    s = exp.(range(0.0, stop = -log(1.0e10), length = r))
    A = U * Diagonal(s) * V'

    Fdefault = bsqr(A)
    Frelaxed = bsqr(A; norm_recomp_tol = 1.0e-12)

    r_def, q_def = _residual_and_orthogonality(A, Fdefault)
    r_rel, q_rel = _residual_and_orthogonality(A, Frelaxed)

    tol = 1.0e3 * eps(Float64) * max(size(A)...)
    @test r_def <= tol
    @test q_def <= tol
    @test r_rel <= tol
    @test q_rel <= tol

    @test_throws ArgumentError bsqr(A; norm_recomp_tol = -1.0)
    @test_throws ArgumentError bsqr(A; norm_recomp_tol = 2.0)
end

@testset "Shared validation path parity" begin
    rng = MersenneTwister(2029)
    A = randn(rng, 12, 8)
    m, n = size(A)
    k = min(m, n)
    tau = zeros(Float64, k)
    jpvt = collect(1:n)
    ws = BSPivotQR.BSWorkspace(m, n, k)

    Ainf = copy(A)
    Ainf[1, 1] = Inf

    @test_throws ArgumentError bsqr!(copy(A); norm_recomp_tol = -1.0)
    @test_throws ArgumentError bsqr!(copy(A), tau, jpvt, ws; norm_recomp_tol = -1.0)

    @test_throws ArgumentError bsqr!(copy(A); norm_recomp_tol = 2.0)
    @test_throws ArgumentError bsqr!(copy(A), tau, jpvt, ws; norm_recomp_tol = 2.0)

    @test_throws ArgumentError bsqr!(copy(A); k = k + 1)
    @test_throws ArgumentError bsqr!(copy(A), tau, jpvt, ws; k = k + 1)

    @test_throws ArgumentError bsqr!(copy(Ainf); check = true)
    @test_throws ArgumentError bsqr!(copy(Ainf), tau, jpvt, ws; check = true)
end

@testset "Optional R11^{-1}R12 return" begin
    rng = MersenneTwister(7)
    A = randn(rng, 18, 12)
    k = 8

    F = bsqr(A; k = k, return_rinv_r12 = true)
    @test F.rinv_r12 !== nothing
    @test rinv_r12(F) == F.rinv_r12

    T = BSPivotQR._packed_to_qt(F)
    R11 = UpperTriangular(Matrix(view(T, 1:F.ksteps, 1:F.ksteps)))
    R12 = Matrix(view(T, 1:F.ksteps, (F.ksteps + 1):size(T, 2)))
    expected = Matrix(R11 \ R12)

    got = F.rinv_r12
    @test size(got) == size(expected)
    @test norm(got - expected) / max(norm(expected), eps(Float64)) <= 5e2 * eps(Float64) * size(A, 2)

    Fnone = bsqr(A; k = k, return_rinv_r12 = false)
    @test Fnone.rinv_r12 === nothing

    B = randn(rng, 10, 10)
    Fsq = bsqr(B; k = min(size(B)...), return_rinv_r12 = true)
    @test Fsq.rinv_r12 !== nothing
    @test size(Fsq.rinv_r12) == (Fsq.ksteps, 0)

    Fzero = bsqr(B; k = 0, return_rinv_r12 = true)
    @test Fzero.rinv_r12 !== nothing
    @test size(Fzero.rinv_r12) == (0, size(B, 2))
end

@testset "Preallocated API and BLAS thread control" begin
    rng = MersenneTwister(314)
    A = randn(rng, 30, 18)
    m, n = size(A)
    k = 12
    kmax = min(m, n)

    Fref = bsqr(A; k = k, check = true, track_inverse_frob = true, return_rinv_r12 = true)

    Awork = copy(A)
    tau = zeros(Float64, k)
    jpvt = zeros(Int, n)
    ws = BSPivotQR.BSWorkspace(m, n, kmax)
    trace = Float64[]

    old_threads = BLAS.get_num_threads()
    ksteps = bsqr!(Awork, tau, jpvt, ws; k = k, check = true, frob_inv_trace = trace, blas_threads = 1)
    @test BLAS.get_num_threads() == old_threads

    @test ksteps == Fref.ksteps
    @test Awork ≈ Fref.factors
    @test tau ≈ Fref.tau
    @test jpvt == Fref.jpvt
    @test trace ≈ Fref.frob_inv_trace

    Ffrom_prealloc = BSQRPivoted(copy(Awork), copy(tau), copy(jpvt), ksteps, copy(trace), Matrix(view(ws.W, 1:ksteps, (ksteps + 1):n)))
    @test rinv_r12(Ffrom_prealloc) ≈ Fref.rinv_r12
end

@testset "Kernel helper invariants and fastpath knobs" begin
    ws = BSPivotQR.BSWorkspace(8, 5, 5)
    ws.s .= [4.0, 4.0, 1.0, 1.0, 0.0]
    ws.wnorm2 .= [0.0, 0.0, 0.0, 0.0, 0.0]

    j1, c1 = BSPivotQR._select_pivot_column!(ws, 1, 5, nothing)
    @test j1 == 1
    @test c1 == 0.25

    # Tie must keep first minimum index because the kernel uses strict '<'.
    j2, _ = BSPivotQR._select_pivot_column!(ws, 3, 5, nothing)
    @test j2 == 3

    tol_scale = eps(Float64) * 10.0
    ws.s[1] = 0.0
    ws.s_ref[1] = 1.0
    @test BSPivotQR._rank_stop_triggered(ws, 1, tol_scale)
    ws.s[1] = 10.0 * tol_scale
    @test !BSPivotQR._rank_stop_triggered(ws, 1, tol_scale)

    keys = [
        "BS_SHORT_WIDE_FASTPATH_ASPECT",
        "BS_SHORT_WIDE_FASTPATH_MMAX",
        "BS_SHORT_WIDE_FASTPATH_NMIN",
    ]
    old = Dict{String,Union{Nothing,String}}()
    for key in keys
        old[key] = haskey(ENV, key) ? ENV[key] : nothing
    end

    try
        haskey(ENV, "BS_SHORT_WIDE_FASTPATH_ASPECT") && delete!(ENV, "BS_SHORT_WIDE_FASTPATH_ASPECT")
        haskey(ENV, "BS_SHORT_WIDE_FASTPATH_MMAX") && delete!(ENV, "BS_SHORT_WIDE_FASTPATH_MMAX")
        haskey(ENV, "BS_SHORT_WIDE_FASTPATH_NMIN") && delete!(ENV, "BS_SHORT_WIDE_FASTPATH_NMIN")
        aspect, mmax, nmin = BSPivotQR._short_wide_fastpath_params()
        @test aspect == BSPivotQR._SHORT_WIDE_FASTPATH_ASPECT
        @test mmax == BSPivotQR._SHORT_WIDE_FASTPATH_MMAX
        @test nmin == BSPivotQR._SHORT_WIDE_FASTPATH_NMIN

        ENV["BS_SHORT_WIDE_FASTPATH_ASPECT"] = "6"
        ENV["BS_SHORT_WIDE_FASTPATH_MMAX"] = "300"
        ENV["BS_SHORT_WIDE_FASTPATH_NMIN"] = "700"
        @test BSPivotQR._use_short_wide_fastpath(128, 1024)
        @test !BSPivotQR._use_short_wide_fastpath(128, 640)

        ENV["BS_SHORT_WIDE_FASTPATH_ASPECT"] = "0"
        @test_throws ArgumentError BSPivotQR._short_wide_fastpath_params()
    finally
        for (key, val) in old
            if val === nothing
                haskey(ENV, key) && delete!(ENV, key)
            else
                ENV[key] = val
            end
        end
    end
end
