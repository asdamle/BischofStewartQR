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

@inline function _validate_norm_recomp_tol(norm_recomp_tol::Float64)
    (0.0 <= norm_recomp_tol <= 1.0) ||
        throw(ArgumentError("norm_recomp_tol must satisfy 0 <= norm_recomp_tol <= 1"))
    return nothing
end

function _validate_bsqr_common_args(
    A::StridedMatrix{Float64},
    k::Integer,
    check::Bool,
    norm_recomp_tol::Float64,
)
    m, n = size(A)
    kmax = min(m, n)
    (0 <= k <= kmax) || throw(ArgumentError("k must satisfy 0 <= k <= min(m, n)"))
    if check
        all(isfinite, A) || throw(ArgumentError("A contains non-finite values"))
    end
    _validate_norm_recomp_tol(norm_recomp_tol)
    return m, n, Int(k)
end

@inline function _validate_prealloc_buffers!(
    tau::Vector{Float64},
    jpvt::Vector{Int},
    k::Int,
    n::Int,
)
    length(tau) >= k ||
        throw(ArgumentError("tau has length $(length(tau)) but needs at least k=$k"))
    length(jpvt) >= n ||
        throw(ArgumentError("jpvt has length $(length(jpvt)) but needs at least n=$n"))
    return nothing
end

@inline function _reset_pivots!(jpvt::Vector{Int}, n::Int)
    @inbounds for j in 1:n
        jpvt[j] = j
    end
    return nothing
end

@inline function _extract_rinv_r12(
    ws::BSWorkspace,
    ksteps::Int,
    n::Int,
)
    if ksteps == 0
        return zeros(Float64, 0, n)
    elseif ksteps < n
        return Matrix(view(ws.W, 1:ksteps, (ksteps + 1):n))
    else
        return zeros(Float64, ksteps, 0)
    end
end

"""
    bsqr(A; kwargs...) -> BSQRPivoted

Out-of-place Bischof-Stewart column-pivoted QR. Converts `A` to `Float64`,
factorizes in-place on a copy, and returns a `BSQRPivoted` wrapper.
"""
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
    m, n, kwork = _validate_bsqr_common_args(A, k, check, norm_recomp_tol)
    _validate_prealloc_buffers!(tau, jpvt, kwork, n)

    ws = _require_workspace(workspace, m, n, kwork)
    if reset_pivots
        _reset_pivots!(jpvt, n)
    end
    fill!(view(tau, 1:kwork), 0.0)
    if frob_inv_trace !== nothing
        empty!(frob_inv_trace)
    end

    return _with_blas_threads(blas_threads) do
        _bsqr_kernel!(
            A,
            tau,
            jpvt,
            ws,
            kwork;
            frob_inv_trace = frob_inv_trace,
            rank_stop = rank_stop,
            norm_recomp_tol = norm_recomp_tol,
        )
    end
end

"""
    bsqr!(A; kwargs...) -> BSQRPivoted

In-place Bischof-Stewart factorization on `A::StridedMatrix{Float64}` that
allocates internal scratch unless `workspace` is provided.
"""
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
    m, n, kwork = _validate_bsqr_common_args(A, k, check, norm_recomp_tol)

    ws = _require_workspace(workspace, m, n, kwork)
    tau = zeros(Float64, kwork)
    jpvt = collect(1:n)
    frob_trace = track_inverse_frob ? Float64[] : nothing

    ksteps = _with_blas_threads(blas_threads) do
        _bsqr_kernel!(
            A,
            tau,
            jpvt,
            ws,
            kwork;
            frob_inv_trace = frob_trace,
            rank_stop = rank_stop,
            norm_recomp_tol = norm_recomp_tol,
        )
    end
    rinv_r12_data = return_rinv_r12 ? _extract_rinv_r12(ws, ksteps, n) : nothing

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
    short_wide_fastpath = _use_short_wide_fastpath(m, n) && _short_wide_fastpath_enabled()
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
            # Bischof-Stewart pivot criterion: minimize (1 + ||w_j||^2) / ||a_j^(i)||^2.
            cj = sj > 0.0 ? (1.0 + ws.wnorm2[j]) / sj : Inf
            # Use strict '<' so ties keep the first minimum, matching reference behavior.
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
        # Rank-stop tolerance is scaled to matrix size/reference norm to avoid
        # stopping on benign floating-point noise.
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

            if short_wide_fastpath
                # Short-wide fast path: materialize W-row and beta in one pass to avoid
                # an extra strided copy over the full trailing width.
                @inbounds for t in 1:nrem
                    j = i + t
                    b = alpha[t] * invdiag
                    beta_vec[t] = b
                    ws.W[i, j] = b
                end
            else
                @inbounds for t in 1:nrem
                    beta_vec[t] = alpha[t] * invdiag
                end
            end

            wnorm2_pivot = ws.wnorm2[i]
            dots = view(ws.dots, 1:nrem)
            # ws.W stores R11^{-1}R12 rows built so far; ws.wnorm2 tracks per-column
            # ||w_j||^2 used directly by the pivot criterion.
            has_prefix = i > 1
            if has_prefix
                Wprefix = view(ws.W, 1:(i - 1), (i + 1):n)
                wpivot = view(ws.W, 1:(i - 1), i)
                BLAS.gemv!(
                    'T',
                    1.0,
                    Wprefix,
                    wpivot,
                    0.0,
                    dots,
                )
                BLAS.ger!(
                    -1.0,
                    wpivot,
                    beta_vec,
                    Wprefix,
                )
            end

            if !short_wide_fastpath
                copyto!(view(ws.W, i, (i + 1):n), beta_vec)
            end

            wcoeff = wnorm2_pivot + 1.0
            recompute_slots = view(ws.work, 1:nrem)
            nrecompute = 0
            if has_prefix
                @inbounds for t in 1:nrem
                    j = i + t
                    b = beta_vec[t]
                    d = dots[t]

                    wn = ws.wnorm2[j] - 2.0 * b * d + b * b * wcoeff
                    ws.wnorm2[j] = wn > 0.0 ? wn : 0.0

                    old_s = ws.s[j]
                    if !(old_s > 0.0)
                        ws.s[j] = 0.0
                        continue
                    end

                    sj = old_s - alpha[t] * alpha[t]
                    sj = sj > 0.0 ? sj : 0.0
                    ws.s[j] = sj

                    # LAPACK-style safeguard: refresh only when the running norm
                    # estimate has decayed enough relative to the last exact reference norm.
                    if sj <= ws.s_ref[j] * norm_recomp_tol
                        nrecompute += 1
                        recompute_slots[nrecompute] = t
                    end
                end
            else
                @inbounds for t in 1:nrem
                    j = i + t
                    b = beta_vec[t]

                    wn = ws.wnorm2[j] + b * b
                    ws.wnorm2[j] = wn > 0.0 ? wn : 0.0

                    old_s = ws.s[j]
                    if !(old_s > 0.0)
                        ws.s[j] = 0.0
                        continue
                    end

                    sj = old_s - alpha[t] * alpha[t]
                    sj = sj > 0.0 ? sj : 0.0
                    ws.s[j] = sj

                    if sj <= ws.s_ref[j] * norm_recomp_tol
                        nrecompute += 1
                        recompute_slots[nrecompute] = t
                    end
                end
            end
            downdate_ns_local += _refresh_recompute_colnorms!(
                A,
                ws,
                i,
                m,
                recompute_slots,
                nrecompute,
                norm_recomp_count,
                kernel_stats,
            )
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
const _SHORT_WIDE_FASTPATH_ASPECT = 4
const _SHORT_WIDE_FASTPATH_MMAX = 256

@inline _use_short_wide_fastpath(m::Int, n::Int) = (m <= _SHORT_WIDE_FASTPATH_MMAX) && (n >= _SHORT_WIDE_FASTPATH_ASPECT * m)
@inline _short_wide_fastpath_enabled() = strip(get(ENV, "BS_SHORT_WIDE_FASTPATH", "1")) != "0"
@inline _householder_lapack_larf_enabled() = strip(get(ENV, "BS_HOUSEHOLDER_LAPACK_LARF", "1")) != "0"

@inline function _refresh_recompute_colnorms!(
    A::StridedMatrix{Float64},
    ws::BSWorkspace,
    i::Int,
    m::Int,
    recompute_slots::AbstractVector{Float64},
    nrecompute::Int,
    norm_recomp_count::Union{Nothing,Base.RefValue{Int}},
    kernel_stats::Union{Nothing,BSKernelStats},
)
    nrecompute == 0 && return 0

    t_down = kernel_stats === nothing ? 0 : time_ns()
    @inbounds for q in 1:nrecompute
        t = Int(recompute_slots[q])
        j = i + t
        if i < m
            tail_norm = BLAS.nrm2(view(A, (i + 1):m, j))
            sj = tail_norm * tail_norm
        else
            sj = 0.0
        end
        ws.s_ref[j] = sj
        ws.s[j] = sj
    end

    if norm_recomp_count !== nothing
        norm_recomp_count[] += nrecompute
    end
    if kernel_stats !== nothing
        kernel_stats.recompute_count += nrecompute
        return time_ns() - t_down
    end
    return 0
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

    if _householder_lapack_larf_enabled()
        # Optional LAPACK path for A[i:m, i+1:n] -= tau * v * (v' * A[i:m, i+1:n]),
        # where v = A[i:m, i] with v[1] == 1 set by the caller.
        LAPACK.larf!(
            'L',
            view(A, i:m, i),
            tau,
            view(A, i:m, (i + 1):n),
            w,
        )
        return nothing
    end

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
