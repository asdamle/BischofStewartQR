using LinearAlgebra
using Random
using Test

function _bsqr_exact_score(M::Matrix{Float64}, i::Int, j::Int)
    m, _ = size(M)
    rho2 = norm(view(M, i:m, j))^2
    if rho2 <= 0.0
        return Inf, 0.0, 0.0
    elseif i == 1
        return 1.0 / rho2, 0.0, rho2
    else
        R11 = UpperTriangular(view(M, 1:(i - 1), 1:(i - 1)))
        wj = R11 \ view(M, 1:(i - 1), j)
        wnorm2 = dot(wj, wj)
        return (1.0 + wnorm2) / rho2, wnorm2, rho2
    end
end

function _bsqr_reference_pivots(A::Matrix{Float64})
    M = copy(A)
    m, n = size(M)
    k = min(m, n)
    p = collect(1:n)
    hist = Int[]

    for i in 1:k
        best_j = i
        best_c = Inf
        for j in i:n
            cj, _, _ = _bsqr_exact_score(M, i, j)
            if cj < best_c
                best_c = cj
                best_j = j
            end
        end

        p[i], p[best_j] = p[best_j], p[i]
        M[:, [i, best_j]] = M[:, [best_j, i]]
        push!(hist, p[i])

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

    return hist
end

function _bsqr_eager_kernel_pivots(A::Matrix{Float64}; rank_stop::Bool = true)
    m, n = size(A)
    k = min(m, n)
    M = copy(A)
    tau = zeros(Float64, k)
    jpvt = collect(1:n)
    ws = BSPivotQR.BSWorkspace(m, n, k)
    hist = Int[]
    BSPivotQR._bsqr_kernel!(M, tau, jpvt, ws, k; pivot_history = hist, rank_stop = rank_stop)
    return hist
end

function _bsqr_lazy_certificate_pivots(A::Matrix{Float64})
    M = copy(A)
    m, n = size(M)
    k = min(m, n)
    p = collect(1:n)
    hist = Int[]

    s_upper = [norm(view(M, :, j))^2 for j in 1:n]
    wnorm2_lb = zeros(Float64, n)
    exact_score = fill(Inf, n)
    exact_s = zeros(Float64, n)
    exact_wnorm2 = zeros(Float64, n)
    refreshed = falses(n)
    refresh_count = 0

    for i in 1:k
        while true
            best_j = 0
            best_c = Inf
            for j in i:n
                if refreshed[j] && exact_score[j] < best_c
                    best_c = exact_score[j]
                    best_j = j
                end
            end

            lb_j = 0
            lb_min = Inf
            for j in i:n
                if refreshed[j]
                    continue
                end
                lb = s_upper[j] > 0.0 ? (1.0 + wnorm2_lb[j]) / s_upper[j] : Inf
                if lb < lb_min
                    lb_min = lb
                    lb_j = j
                end
            end

            if best_j != 0 && (lb_j == 0 || best_c < lb_min)
                if best_j != i
                    p[i], p[best_j] = p[best_j], p[i]
                    M[:, [i, best_j]] = M[:, [best_j, i]]
                    s_upper[i], s_upper[best_j] = s_upper[best_j], s_upper[i]
                    wnorm2_lb[i], wnorm2_lb[best_j] = wnorm2_lb[best_j], wnorm2_lb[i]
                    exact_score[i], exact_score[best_j] = exact_score[best_j], exact_score[i]
                    exact_s[i], exact_s[best_j] = exact_s[best_j], exact_s[i]
                    exact_wnorm2[i], exact_wnorm2[best_j] = exact_wnorm2[best_j], exact_wnorm2[i]
                    refreshed[i], refreshed[best_j] = refreshed[best_j], refreshed[i]
                end
                break
            end

            @test lb_j != 0
            exact_score[lb_j], exact_wnorm2[lb_j], exact_s[lb_j] = _bsqr_exact_score(M, i, lb_j)
            refreshed[lb_j] = true
            refresh_count += 1
        end

        push!(hist, p[i])
        pivot_wnorm2 = exact_wnorm2[i]
        denom_factor = 1.0 + pivot_wnorm2

        if norm(view(M, i:m, i)) <= 1e-12
            break
        end

        tau_i, beta_i = BSPivotQR._householder!(view(M, i:m, i))
        if tau_i != 0.0
            M[i, i] = 1.0
            BSPivotQR._apply_householder_left!(M, i, tau_i, zeros(Float64, n))
        end
        M[i, i] = beta_i

        for j in (i + 1):n
            if refreshed[j]
                wnorm2_lb[j] = exact_wnorm2[j] / denom_factor
                s_upper[j] = exact_s[j]
            else
                wnorm2_lb[j] /= denom_factor
            end
            refreshed[j] = false
            exact_score[j] = Inf
            exact_s[j] = 0.0
            exact_wnorm2[j] = 0.0
        end
    end

    return hist, refresh_count
end

@testset "Lazy BSQR certificate assessment" begin
    rng = MersenneTwister(20260310)

    cases = Matrix{Float64}[]
    push!(cases, randn(rng, 8, 6))
    push!(cases, randn(rng, 6, 8))
    push!(cases, randn(rng, 7, 7))

    A_rank = randn(rng, 8, 6)
    A_rank[:, 5] .= A_rank[:, 2]
    A_rank[:, 6] .= A_rank[:, 2] .+ 1e-13 * randn(rng, 8)
    push!(cases, A_rank)

    A_tie = randn(rng, 8, 6)
    A_tie[:, 3] .= A_tie[:, 1] + 1e-14 * randn(rng, 8)
    A_tie[:, 4] .= A_tie[:, 2] + 1e-14 * randn(rng, 8)
    push!(cases, A_tie)

    r = 6
    U = Matrix(qr(randn(rng, 8, r)).Q)
    V = Matrix(qr(randn(rng, 6, r)).Q)
    s = exp.(range(0.0, stop = -log(1.0e10), length = r))
    push!(cases, U * Diagonal(s) * V')

    Qn = Matrix(qr(randn(rng, 24, 6)).Q)
    push!(cases, Matrix(transpose(Qn[:, 1:6])))

    for A in cases
        eager = _bsqr_reference_pivots(copy(A))
        lazy, _ = _bsqr_lazy_certificate_pivots(copy(A))
        @test lazy == eager
    end

    for _ in 1:12
        m = rand(rng, 5:9)
        n = rand(rng, 5:9)
        A = randn(rng, m, n)
        eager = _bsqr_reference_pivots(copy(A))
        lazy, _ = _bsqr_lazy_certificate_pivots(copy(A))
        @test lazy == eager
    end
end

@testset "Lazy BSQR certificate can skip refreshes on separated scores" begin
    A = Matrix(Diagonal([100.0, 10.0, 1.0, 0.1]))
    eager = _bsqr_reference_pivots(copy(A))
    lazy, refresh_count = _bsqr_lazy_certificate_pivots(copy(A))

    @test lazy == eager
    @test refresh_count < 10
end

@testset "BLAS-lazy prototype matches eager pivot history" begin
    rng = MersenneTwister(20260311)

    cases = Matrix{Float64}[]
    push!(cases, randn(rng, 8, 6))
    push!(cases, randn(rng, 6, 8))

    A_tie = Matrix(I, 5, 5)
    push!(cases, A_tie)

    A_rank = randn(rng, 8, 6)
    A_rank[:, 5] .= A_rank[:, 2]
    A_rank[:, 6] .= A_rank[:, 2] .+ 1e-13 * randn(rng, 8)
    push!(cases, A_rank)

    r = 6
    U = Matrix(qr(randn(rng, 8, r)).Q)
    V = Matrix(qr(randn(rng, 6, r)).Q)
    s = exp.(range(0.0, stop = -log(1.0e10), length = r))
    push!(cases, U * Diagonal(s) * V')

    Qn = Matrix(qr(randn(rng, 24, 6)).Q)
    push!(cases, Matrix(transpose(Qn[:, 1:6])))

    for A in cases
        eager = _bsqr_eager_kernel_pivots(copy(A); rank_stop = true)
        hist = Int[]
        stats = BSPivotQR.BSLazyBLASStats()
        F = BSPivotQR._bsqr_lazy_blas(
            copy(A);
            lazy_batch_size = 1,
            pivot_history = hist,
            selection_stats = stats,
        )
        @test hist == eager
        @test perm(F)[1:length(eager)] == eager
        @test stats.refresh_count >= length(eager)
    end
end

@testset "BLAS-lazy prototype matches eager factorization quality" begin
    rng = MersenneTwister(20260312)
    shapes = ((24, 16), (16, 24), (20, 20))

    for (m, n) in shapes
        A = randn(rng, m, n)
        Feager = bsqr(A; rank_stop = false)
        Flazy = BSPivotQR._bsqr_lazy_blas(copy(A); rank_stop = false, lazy_batch_fraction = 0.25)

        @test Flazy.ksteps == Feager.ksteps
        @test perm(Flazy) == perm(Feager)

        r_lazy, q_lazy = _residual_and_orthogonality(A, Flazy)
        tol = 5.0e2 * eps(Float64) * max(m, n)
        @test r_lazy <= tol
        @test q_lazy <= tol
    end
end

@testset "BLAS-lazy prototype skips and falls back as expected" begin
    A_sep = Matrix(Diagonal([100.0, 10.0, 1.0, 0.1]))
    stats_sep = BSPivotQR.BSLazyBLASStats()
    F_sep = BSPivotQR._bsqr_lazy_blas(
        copy(A_sep);
        rank_stop = false,
        lazy_batch_size = 1,
        selection_stats = stats_sep,
    )
    @test perm(F_sep) == _bsqr_reference_pivots(copy(A_sep))
    @test stats_sep.refresh_count < 10

    A_tie = Matrix{Float64}(I, 4, 4)
    stats_tie = BSPivotQR.BSLazyBLASStats()
    F_tie = BSPivotQR._bsqr_lazy_blas(
        copy(A_tie);
        rank_stop = false,
        lazy_batch_size = 1,
        selection_stats = stats_tie,
    )
    @test perm(F_tie) == _bsqr_reference_pivots(copy(A_tie))
    @test stats_tie.full_refresh_rounds > 0
end
