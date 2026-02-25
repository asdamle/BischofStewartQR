struct BSQRPivoted <: Factorization{Float64}
    factors::Matrix{Float64}
    tau::Vector{Float64}
    jpvt::Vector{Int}
    ksteps::Int
    frob_inv_trace::Union{Nothing,Vector{Float64}}
    rinv_r12::Union{Nothing,Matrix{Float64}}
end

Base.size(F::BSQRPivoted) = size(F.factors)

mutable struct BSKernelStats
    pivot_select_ns::Int
    householder_ns::Int
    apply_reflector_ns::Int
    w_update_ns::Int
    norm_downdate_ns::Int
    recompute_count::Int
end

BSKernelStats() = BSKernelStats(0, 0, 0, 0, 0, 0)

function _with_blas_threads(f::F, blas_threads::Union{Nothing,Integer}) where {F<:Function}
    if blas_threads === nothing
        return f()
    end
    nth = Int(blas_threads)
    nth >= 1 || throw(ArgumentError("blas_threads must be >= 1"))
    old = BLAS.get_num_threads()
    BLAS.set_num_threads(nth)
    try
        return f()
    finally
        BLAS.set_num_threads(old)
    end
end

function bsqr(
    A::AbstractMatrix{<:Real};
    k::Integer = min(size(A)...),
    check::Bool = true,
    track_inverse_frob::Bool = false,
    return_rinv_r12::Bool = false,
    rank_stop::Bool = true,
    norm_recomp_tol::Float64 = _DEFAULT_NORM_RECOMP_TOL,
    blas_threads::Union{Nothing,Integer} = nothing,
)
    Acopy = Matrix{Float64}(A)
    return bsqr!(
        Acopy;
        k = k,
        check = check,
        track_inverse_frob = track_inverse_frob,
        return_rinv_r12 = return_rinv_r12,
        rank_stop = rank_stop,
        norm_recomp_tol = norm_recomp_tol,
        blas_threads = blas_threads,
        workspace = nothing,
    )
end

"""
    bsqr!(A, tau, jpvt, workspace; kwargs...) -> ksteps

Allocation-minimal in-place Bischof-Stewart kernel for repeated calls.
`tau`, `jpvt`, and `workspace` are caller-provided scratch/state buffers.
"""
function bsqr!(
    A::StridedMatrix{Float64},
    tau::Vector{Float64},
    jpvt::Vector{Int},
    workspace::BSWorkspace;
    k::Integer = min(size(A)...),
    check::Bool = true,
    reset_pivots::Bool = true,
    frob_inv_trace::Union{Nothing,Vector{Float64}} = nothing,
    rank_stop::Bool = true,
    norm_recomp_tol::Float64 = _DEFAULT_NORM_RECOMP_TOL,
    blas_threads::Union{Nothing,Integer} = nothing,
)
    m, n = size(A)
    kmax = min(m, n)
    (0 <= k <= kmax) || throw(ArgumentError("k must satisfy 0 <= k <= min(m, n)"))
    kwork = Int(k)
    length(tau) >= k || throw(ArgumentError("tau has length $(length(tau)) but needs at least k=$k"))
    length(jpvt) >= n || throw(ArgumentError("jpvt has length $(length(jpvt)) but needs at least n=$n"))

    if check
        all(isfinite, A) || throw(ArgumentError("A contains non-finite values"))
    end
    (0.0 <= norm_recomp_tol <= 1.0) || throw(ArgumentError("norm_recomp_tol must satisfy 0 <= norm_recomp_tol <= 1"))

    ws = _require_workspace(workspace, m, n, kwork)
    if reset_pivots
        @inbounds for j in 1:n
            jpvt[j] = j
        end
    end
    fill!(view(tau, 1:k), 0.0)
    if frob_inv_trace !== nothing
        empty!(frob_inv_trace)
    end

    return _with_blas_threads(blas_threads) do
        _bsqr_kernel!(
            A,
            tau,
            jpvt,
            ws,
            k;
            frob_inv_trace = frob_inv_trace,
            rank_stop = rank_stop,
            norm_recomp_tol = norm_recomp_tol,
        )
    end
end

function bsqr!(
    A::StridedMatrix{Float64};
    k::Integer = min(size(A)...),
    check::Bool = true,
    track_inverse_frob::Bool = false,
    return_rinv_r12::Bool = false,
    rank_stop::Bool = true,
    norm_recomp_tol::Float64 = _DEFAULT_NORM_RECOMP_TOL,
    blas_threads::Union{Nothing,Integer} = nothing,
    workspace::Union{Nothing,BSWorkspace} = nothing,
)
    m, n = size(A)
    kmax = min(m, n)
    (0 <= k <= kmax) || throw(ArgumentError("k must satisfy 0 <= k <= min(m, n)"))
    kwork = Int(k)

    if check
        all(isfinite, A) || throw(ArgumentError("A contains non-finite values"))
    end
    (0.0 <= norm_recomp_tol <= 1.0) || throw(ArgumentError("norm_recomp_tol must satisfy 0 <= norm_recomp_tol <= 1"))

    ws = _require_workspace(workspace, m, n, kwork)
    tau = zeros(Float64, k)
    jpvt = collect(1:n)
    frob_trace = track_inverse_frob ? Float64[] : nothing

    ksteps = _with_blas_threads(blas_threads) do
        _bsqr_kernel!(
            A,
            tau,
            jpvt,
            ws,
            k;
            frob_inv_trace = frob_trace,
            rank_stop = rank_stop,
            norm_recomp_tol = norm_recomp_tol,
        )
    end
    rinv_r12_data = if return_rinv_r12
        if ksteps == 0
            zeros(Float64, 0, n)
        elseif ksteps < n
            Matrix(view(ws.W, 1:ksteps, (ksteps + 1):n))
        else
            zeros(Float64, ksteps, 0)
        end
    else
        nothing
    end

    factors = A isa Matrix{Float64} ? A : Matrix(A)
    return BSQRPivoted(factors, tau, jpvt, ksteps, frob_trace, rinv_r12_data)
end

function _bsqr_kernel!(
    A::StridedMatrix{Float64},
    tau::Vector{Float64},
    jpvt::Vector{Int},
    ws::BSWorkspace,
    k::Int;
    frob_inv_trace::Union{Nothing,Vector{Float64}} = nothing,
    pivot_history::Union{Nothing,Vector{Int}} = nothing,
    rank_stop::Bool = true,
    norm_recomp_tol::Float64 = _DEFAULT_NORM_RECOMP_TOL,
    norm_recomp_count::Union{Nothing,Base.RefValue{Int}} = nothing,
    kernel_stats::Union{Nothing,BSKernelStats} = nothing,
)
    m, n = size(A)
    k == 0 && return 0
    tol_scale = eps(Float64) * max(m, n)

    @inbounds for j in 1:n
        nj = BLAS.nrm2(view(A, :, j))
        sj = nj * nj
        ws.s[j] = sj
        ws.s_ref[j] = sj
        ws.wnorm2[j] = 0.0
    end

    ksteps = 0
    for i in 1:k
        best_j = i
        best_c = Inf
        t_pivot = kernel_stats === nothing ? 0 : time_ns()

        @inbounds for j in i:n
            sj = ws.s[j]
            cj = sj > 0.0 ? (1.0 + ws.wnorm2[j]) / sj : Inf
            if cj < best_c
                best_c = cj
                best_j = j
            end
        end
        if kernel_stats !== nothing
            kernel_stats.pivot_select_ns += time_ns() - t_pivot
        end

        if best_j != i
            _swap_columns!(A, i, best_j)
            if i > 1
                _swap_columns_prefix!(ws.W, i, best_j, i - 1)
            end
            _swap_entries!(ws.s, i, best_j)
            _swap_entries!(ws.s_ref, i, best_j)
            _swap_entries!(ws.wnorm2, i, best_j)
            _swap_entries!(jpvt, i, best_j)
        end

        if pivot_history !== nothing
            push!(pivot_history, jpvt[i])
        end

        pivot_norm2 = ws.s[i]
        atol = tol_scale * max(ws.s_ref[i], 1.0)
        if !(pivot_norm2 > atol)
            if rank_stop
                tau[i] = 0.0
                break
            end
        end

        t_hh = kernel_stats === nothing ? 0 : time_ns()
        tau_i, beta_i = _householder!(view(A, i:m, i))
        if kernel_stats !== nothing
            kernel_stats.householder_ns += time_ns() - t_hh
        end
        tau[i] = tau_i

        if tau_i != 0.0 && i < n
            A[i, i] = 1.0
            t_apply = kernel_stats === nothing ? 0 : time_ns()
            _apply_householder_left!(A, i, tau_i, ws.work)
            if kernel_stats !== nothing
                kernel_stats.apply_reflector_ns += time_ns() - t_apply
            end
        end

        A[i, i] = beta_i

        ws.s[i] = beta_i * beta_i
        ws.s_ref[i] = ws.s[i]

        nrem = n - i
        if nrem > 0
            t_update = kernel_stats === nothing ? 0 : time_ns()
            downdate_ns_local = 0
            alpha = view(A, i, (i + 1):n)
            beta_vec = view(ws.beta, 1:nrem)
            invdiag = beta_i != 0.0 ? (1.0 / beta_i) : 0.0
            @inbounds for t in 1:nrem
                beta_vec[t] = alpha[t] * invdiag
            end

            wnorm2_pivot = ws.wnorm2[i]
            dots = view(ws.dots, 1:nrem)
            if i > 1
                BLAS.gemv!(
                    'T',
                    1.0,
                    view(ws.W, 1:(i - 1), (i + 1):n),
                    view(ws.W, 1:(i - 1), i),
                    0.0,
                    dots,
                )
                BLAS.ger!(
                    -1.0,
                    view(ws.W, 1:(i - 1), i),
                    beta_vec,
                    view(ws.W, 1:(i - 1), (i + 1):n),
                )
            end

            copyto!(view(ws.W, i, (i + 1):n), beta_vec)

            wcoeff = wnorm2_pivot + 1.0
            if i > 1
                @inbounds for t in 1:nrem
                    j = i + t
                    b = beta_vec[t]
                    d = dots[t]

                    wn = ws.wnorm2[j] - 2.0 * b * d + b * b * wcoeff
                    ws.wnorm2[j] = wn > 0.0 ? wn : 0.0

                    if kernel_stats === nothing
                        recomputed = _downdate_colnorm2!(A, ws, i, j, alpha[t], m, norm_recomp_tol)
                    else
                        t_down = time_ns()
                        recomputed = _downdate_colnorm2!(A, ws, i, j, alpha[t], m, norm_recomp_tol)
                        downdate_ns_local += time_ns() - t_down
                    end
                    if norm_recomp_count !== nothing && recomputed
                        norm_recomp_count[] += 1
                    end
                    if kernel_stats !== nothing && recomputed
                        kernel_stats.recompute_count += 1
                    end
                end
            else
                @inbounds for t in 1:nrem
                    j = i + t
                    b = beta_vec[t]

                    wn = ws.wnorm2[j] + b * b * wcoeff
                    ws.wnorm2[j] = wn > 0.0 ? wn : 0.0

                    if kernel_stats === nothing
                        recomputed = _downdate_colnorm2!(A, ws, i, j, alpha[t], m, norm_recomp_tol)
                    else
                        t_down = time_ns()
                        recomputed = _downdate_colnorm2!(A, ws, i, j, alpha[t], m, norm_recomp_tol)
                        downdate_ns_local += time_ns() - t_down
                    end
                    if norm_recomp_count !== nothing && recomputed
                        norm_recomp_count[] += 1
                    end
                    if kernel_stats !== nothing && recomputed
                        kernel_stats.recompute_count += 1
                    end
                end
            end
            if kernel_stats !== nothing
                elapsed = time_ns() - t_update
                kernel_stats.norm_downdate_ns += downdate_ns_local
                wns = elapsed - downdate_ns_local
                kernel_stats.w_update_ns += wns > 0 ? wns : 0
            end
        end

        ksteps = i
        if frob_inv_trace !== nothing
            push!(frob_inv_trace, best_c)
        end
    end

    return ksteps
end

const _DEFAULT_NORM_RECOMP_TOL = sqrt(eps(Float64))

@inline function _downdate_colnorm2!(
    A::StridedMatrix{Float64},
    ws::BSWorkspace,
    i::Int,
    j::Int,
    alpha_ij::Float64,
    m::Int,
    norm_recomp_tol::Float64,
)
    old_s = ws.s[j]
    if !(old_s > 0.0)
        ws.s[j] = 0.0
        return false
    end

    sj = old_s - alpha_ij * alpha_ij
    sj = sj > 0.0 ? sj : 0.0

    # LAPACK-like partial norm safeguard: refresh only when norm estimate
    # has decayed significantly relative to its reference exact norm.
    if sj <= ws.s_ref[j] * norm_recomp_tol
        if i < m
            tail_norm = BLAS.nrm2(view(A, (i + 1):m, j))
            sj = tail_norm * tail_norm
        else
            sj = 0.0
        end
        ws.s_ref[j] = sj
        ws.s[j] = sj
        return true
    end

    ws.s[j] = sj
    return false
end

function _swap_entries!(v::Vector{T}, i::Int, j::Int) where {T}
    v[i], v[j] = v[j], v[i]
    return nothing
end

function _swap_columns!(A::StridedMatrix{Float64}, i::Int, j::Int)
    m = size(A, 1)
    @inbounds for r in 1:m
        A[r, i], A[r, j] = A[r, j], A[r, i]
    end
    return nothing
end

function _swap_columns_prefix!(A::StridedMatrix{Float64}, i::Int, j::Int, rmax::Int)
    @inbounds for r in 1:rmax
        A[r, i], A[r, j] = A[r, j], A[r, i]
    end
    return nothing
end

function _householder!(x::StridedVector{Float64})
    n = length(x)
    alpha = x[1]

    if n == 1
        return 0.0, alpha
    end

    xtail = view(x, 2:n)
    xnorm = BLAS.nrm2(xtail)
    if xnorm == 0.0
        return 0.0, alpha
    end

    beta = -copysign(hypot(alpha, xnorm), alpha)
    tau = (beta - alpha) / beta
    scale = 1.0 / (alpha - beta)
    BLAS.scal!(scale, xtail)
    x[1] = beta

    return tau, beta
end

function _apply_householder_left!(
    A::StridedMatrix{Float64},
    i::Int,
    tau::Float64,
    work::Vector{Float64},
)
    m, n = size(A)
    cols = n - i
    cols > 0 || return nothing
    w = view(work, 1:cols)

    copyto!(w, view(A, i, (i + 1):n))
    if i < m
        BLAS.gemv!(
            'T',
            1.0,
            view(A, (i + 1):m, (i + 1):n),
            view(A, (i + 1):m, i),
            1.0,
            w,
        )
    end
    BLAS.scal!(tau, w)

    BLAS.axpy!(-1.0, w, view(A, i, (i + 1):n))
    if i < m
        BLAS.ger!(-1.0, view(A, (i + 1):m, i), w, view(A, (i + 1):m, (i + 1):n))
    end

    return nothing
end
