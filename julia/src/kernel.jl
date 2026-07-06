"""
    BSQRPivoted <: Factorization{Float64}

Result of [`bsqr`](@ref)/[`bsqr!`](@ref): a Bischof-Stewart column-pivoted QR.
Access the results with the exported accessors [`R`](@ref), [`Q`](@ref),
[`perm`](@ref), [`rinv_r12`](@ref), and [`reconstruct`](@ref); see `R` for the
factorization identities (including the early-stop case `ksteps < min(m, n)`).

Fields:
- `factors` — the reduced `m×n` matrix in LAPACK packed form: the upper
  trapezoid of the first `ksteps` rows holds `R`; the Householder reflector
  tails sit below the diagonal of the first `ksteps` columns. Together with
  `tau` they define `Q` implicitly (materialize with `Q(F)`).
- `tau` — the `ksteps` Householder coefficients.
- `jpvt` — the length-`n` column permutation (`jpvt[1:ksteps]` are the
  selected columns, in pivot order).
- `ksteps` — factorization steps performed: the requested `k` unless
  `rank_stop` fired earlier.
- `frob_inv_trace` — per-step pivot-criterion values, i.e. the exact
  increments to `‖R11⁻¹‖_F²` (their cumulative sum is the running squared
  inverse Frobenius norm); `nothing` unless `track_inverse_frob = true`.
- `rinv_r12` — the `ksteps×(n-ksteps)` matrix `R11⁻¹R12`; `nothing` unless
  `return_rinv_r12 = true`.
"""
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
    # The kernel's BLAS/LAPACK calls assume LAPACK layout: columns contiguous
    # in memory (unit first stride). A row-strided view type-checks as
    # StridedMatrix but silently factorizes incorrectly, so reject it here.
    stride(A, 1) == 1 || throw(ArgumentError(
        "bsqr requires unit column stride (stride(A, 1) == 1); pass a contiguous copy"))
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
factorizes in-place on a copy, and returns a [`BSQRPivoted`](@ref); access the
results with [`R`](@ref), [`Q`](@ref), [`perm`](@ref), [`rinv_r12`](@ref), and
[`reconstruct`](@ref).

# Keyword arguments
- `k::Integer = min(size(A)...)`: number of factorization steps. Early stop
  (`k < min(m, n)`) selects and reduces only `k` columns; see [`R`](@ref) for
  the identities that hold in that case.
- `check::Bool = true`: validate that `A` is finite.
- `track_inverse_frob::Bool = false`: record the per-step pivot criterion
  (the exact increment to `‖R11⁻¹‖_F²`) in the result's `frob_inv_trace`.
- `return_rinv_r12::Bool = false`: store `R11⁻¹R12` in the result's
  `rinv_r12`, read directly from the kernel workspace (no extra triangular
  solve); retrieve it with [`rinv_r12`](@ref).
- `rank_stop::Bool = false`: stop early once the selected pivot column is
  numerically negligible (running squared norm at or below an
  `eps·max(m,n)`-scaled tolerance); the result's `ksteps` records the steps
  actually taken.
- `norm_recomp_tol::Float64 = sqrt(eps)`: running column-norm recompute
  safeguard in `[0, 1]` (Businger-Golub): a downdated squared norm that decays
  past this fraction of its last exact value is recomputed from scratch.
- `blas_threads::Union{Nothing,Integer} = nothing`: temporarily pin the BLAS
  thread count for the duration of the call (restored afterwards).

Numerical domain: the pivot criterion tracks *squared* column norms, so the
selection guarantee assumes column norms roughly within `[1e-150, 1e150]`
(squared norms representable in `Float64`). Outside that range the
factorization remains exact but pivot quality degrades silently.
"""
function bsqr(
    A::AbstractMatrix{<:Real};
    k::Integer = min(size(A)...),
    check::Bool = true,
    track_inverse_frob::Bool = false,
    return_rinv_r12::Bool = false,
    rank_stop::Bool = false,
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

Allocation-minimal in-place Bischof-Stewart kernel for repeated calls. The
unblocked kernel is allocation-free after warm-up; the default panel kernel
additionally allocates its panel-local scratch (a few `n·nb`-sized buffers)
per call — set `BS_PANEL_NB=0` to force the unblocked kernel when strict
allocation-freedom matters more than the panel kernel's speed. No
`BSQRPivoted` is built; the results live in the caller-provided buffers:
`A` is overwritten with the packed factorization (upper trapezoid of the
first `ksteps` rows = `R`, reflector tails below the diagonal), `tau`
(`length ≥ k` required) receives the Householder coefficients, `jpvt`
(`length ≥ n` required) the column permutation. Returns the number of steps
performed (< `k` only when `rank_stop` fires).

Keyword arguments as in [`bsqr`](@ref) (`k`, `check`, `rank_stop`,
`norm_recomp_tol`, `blas_threads`), except:
- `reset_pivots::Bool = true`: reinitialize `jpvt` to `1:n` before running.
- `frob_inv_trace::Union{Nothing,Vector{Float64}} = nothing`: caller-provided
  vector for the per-step criterion trace (emptied, then appended to at each
  step) — the preallocated counterpart of `track_inverse_frob`.
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
    rank_stop::Bool = false,
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
        _run_kernel!(
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

# Kernel dispatch: the panel kernel is the default for problems above
# the empirical crossover; BS_PANEL_NB=0|1 forces the unblocked kernel and
# BS_PANEL_MIN_KN overrides the crossover (see docs/P3_BLOCKED_BSQR.md).
function _run_kernel!(A, tau, jpvt, ws, k; kwargs...)
    nb = _panel_nb()
    if nb >= 2 && k * size(A, 2) >= _panel_min_kn()
        return _bsqr_kernel_panel!(A, tau, jpvt, ws, k, nb; kwargs...)
    end
    return _bsqr_kernel!(A, tau, jpvt, ws, k; kwargs...)
end

"""
    bsqr!(A; kwargs...) -> BSQRPivoted

In-place variant of [`bsqr`](@ref) on `A::StridedMatrix{Float64}`: `A` is
overwritten with the packed factorization and becomes the `factors` field of
the returned [`BSQRPivoted`](@ref). `A` must have unit first stride (columns
contiguous in memory, the LAPACK layout) — offset views are fine, row-strided
views are rejected. Accepts the same keyword arguments as `bsqr`, plus
`workspace::Union{Nothing,BSWorkspace} = nothing` to reuse preallocated
scratch across calls (allocated internally when `nothing`).
"""
function bsqr!(
    A::StridedMatrix{Float64};
    k::Integer = min(size(A)...),
    check::Bool = true,
    track_inverse_frob::Bool = false,
    return_rinv_r12::Bool = false,
    rank_stop::Bool = false,
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
        _run_kernel!(
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

@inline function _init_running_colnorms!(A::StridedMatrix{Float64}, ws::BSWorkspace, n::Int)
    @inbounds for j in 1:n
        nj = BLAS.nrm2(view(A, :, j))
        sj = nj * nj
        ws.s[j] = sj
        ws.s_ref[j] = sj
        ws.wnorm2[j] = 0.0
    end
    return nothing
end

@inline function _select_pivot_column!(
    ws::BSWorkspace,
    i::Int,
    n::Int,
    kernel_stats::Union{Nothing,BSKernelStats},
)
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
    return best_j, best_c
end

@inline function _swap_pivot_state!(
    A::StridedMatrix{Float64},
    ws::BSWorkspace,
    jpvt::Vector{Int},
    i::Int,
    best_j::Int,
)
    best_j == i && return nothing
    _swap_columns!(A, i, best_j)
    if i > 1
        _swap_columns_prefix!(ws.W, i, best_j, i - 1)
    end
    _swap_entries!(ws.s, i, best_j)
    _swap_entries!(ws.s_ref, i, best_j)
    _swap_entries!(ws.wnorm2, i, best_j)
    _swap_entries!(jpvt, i, best_j)
    return nothing
end

@inline function _rank_stop_triggered(
    ws::BSWorkspace,
    i::Int,
    tol_scale::Float64,
)
    pivot_norm2 = ws.s[i]
    atol = tol_scale * max(ws.s_ref[i], 1.0)
    # Rank-stop tolerance scales with matrix size/reference norm to avoid
    # stopping on benign floating-point noise.
    return !(pivot_norm2 > atol)
end

@inline function _householder_stage!(
    A::StridedMatrix{Float64},
    tau::Vector{Float64},
    i::Int,
    m::Int,
    n::Int,
    ws::BSWorkspace,
    kernel_stats::Union{Nothing,BSKernelStats},
)
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
    return beta_i
end

@inline function _materialize_beta_row!(
    A::StridedMatrix{Float64},
    ws::BSWorkspace,
    i::Int,
    n::Int,
    invdiag::Float64,
    short_wide_fastpath::Bool,
)
    nrem = n - i
    alpha = view(A, i, (i + 1):n)
    beta_vec = view(ws.beta, 1:nrem)

    if short_wide_fastpath
        # Materialize the W row while forming beta to avoid an extra strided copy.
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
    return alpha, beta_vec
end

@inline function _update_trailing_state!(
    A::StridedMatrix{Float64},
    ws::BSWorkspace,
    i::Int,
    m::Int,
    n::Int,
    beta_i::Float64,
    short_wide_fastpath::Bool,
    norm_recomp_tol::Float64,
    norm_recomp_count::Union{Nothing,Base.RefValue{Int}},
    kernel_stats::Union{Nothing,BSKernelStats},
)
    nrem = n - i
    nrem > 0 || return 0

    t_update = kernel_stats === nothing ? 0 : time_ns()
    invdiag = beta_i != 0.0 ? (1.0 / beta_i) : 0.0
    alpha, beta_vec = _materialize_beta_row!(A, ws, i, n, invdiag, short_wide_fastpath)

    dots = view(ws.dots, 1:nrem)
    has_prefix = i > 1
    if has_prefix
        Wprefix = view(ws.W, 1:(i - 1), (i + 1):n)
        wpivot = view(ws.W, 1:(i - 1), i)
        BLAS.gemv!('T', 1.0, Wprefix, wpivot, 0.0, dots)
        BLAS.ger!(-1.0, wpivot, beta_vec, Wprefix)
    end

    if !short_wide_fastpath
        copyto!(view(ws.W, i, (i + 1):n), beta_vec)
    end

    # ws.W stores rows of R11^{-1}R12 built so far. ws.wnorm2 tracks ||w_j||^2,
    # which appears directly in the pivot criterion.
    wcoeff = ws.wnorm2[i] + 1.0
    recompute_slots = view(ws.recompute_idx, 1:nrem)
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

            # LAPACK-style safeguard: refresh only when the running norm estimate
            # has decayed enough relative to the last exact reference norm.
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

    downdate_ns = _refresh_recompute_colnorms!(
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
        kernel_stats.norm_downdate_ns += downdate_ns
        wns = elapsed - downdate_ns
        kernel_stats.w_update_ns += wns > 0 ? wns : 0
    end
    return nrecompute
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
    recomp_history::Union{Nothing,Vector{Int}} = nothing,
    kernel_stats::Union{Nothing,BSKernelStats} = nothing,
)
    m, n = size(A)
    short_wide_fastpath = _use_short_wide_fastpath(m, n) && _short_wide_fastpath_enabled()
    k == 0 && return 0
    tol_scale = eps(Float64) * max(m, n)
    _init_running_colnorms!(A, ws, n)

    ksteps = 0
    for i in 1:k
        best_j, best_c = _select_pivot_column!(ws, i, n, kernel_stats)
        _swap_pivot_state!(A, ws, jpvt, i, best_j)
        if pivot_history !== nothing
            push!(pivot_history, jpvt[i])
        end

        if _rank_stop_triggered(ws, i, tol_scale) && rank_stop
            tau[i] = 0.0
            break
        end

        beta_i = _householder_stage!(A, tau, i, m, n, ws, kernel_stats)
        nrecomp = _update_trailing_state!(
            A,
            ws,
            i,
            m,
            n,
            beta_i,
            short_wide_fastpath,
            norm_recomp_tol,
            norm_recomp_count,
            kernel_stats,
        )
        if recomp_history !== nothing
            push!(recomp_history, nrecomp)
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
const _SHORT_WIDE_FASTPATH_NMIN = 0
const _RECOMP_NORM2_LOOP_THRESHOLD = 256

@inline function _env_or_default_int(key::String, default::Int; minval::Int = 0)
    raw = strip(get(ENV, key, ""))
    isempty(raw) && return default
    v = tryparse(Int, raw)
    v === nothing && throw(ArgumentError("$key must be an integer"))
    v >= minval || throw(ArgumentError("$key must be >= $minval"))
    return v
end

@inline function _short_wide_fastpath_params()
    aspect = _env_or_default_int(
        "BS_SHORT_WIDE_FASTPATH_ASPECT",
        _SHORT_WIDE_FASTPATH_ASPECT;
        minval = 1,
    )
    mmax = _env_or_default_int(
        "BS_SHORT_WIDE_FASTPATH_MMAX",
        _SHORT_WIDE_FASTPATH_MMAX;
        minval = 1,
    )
    nmin = _env_or_default_int(
        "BS_SHORT_WIDE_FASTPATH_NMIN",
        _SHORT_WIDE_FASTPATH_NMIN;
        minval = 0,
    )
    return aspect, mmax, nmin
end

@inline function _use_short_wide_fastpath(m::Int, n::Int)
    aspect, mmax, nmin = _short_wide_fastpath_params()
    nbound = max(nmin, aspect * m)
    return (m <= mmax) && (n >= nbound)
end

@inline _short_wide_fastpath_enabled() = strip(get(ENV, "BS_SHORT_WIDE_FASTPATH", "1")) != "0"

@inline function _tail_colnorm2(
    A::StridedMatrix{Float64},
    row_start::Int,
    row_end::Int,
    col::Int,
)
    len = row_end - row_start + 1
    len > 0 || return 0.0
    if len <= _RECOMP_NORM2_LOOP_THRESHOLD
        acc = 0.0
        @inbounds @simd for r in row_start:row_end
            v = A[r, col]
            acc = muladd(v, v, acc)
        end
        return acc
    end
    nrm = BLAS.nrm2(view(A, row_start:row_end, col))
    return nrm * nrm
end

@inline function _refresh_recompute_colnorms!(
    A::StridedMatrix{Float64},
    ws::BSWorkspace,
    i::Int,
    m::Int,
    recompute_slots::AbstractVector{Int},
    nrecompute::Int,
    norm_recomp_count::Union{Nothing,Base.RefValue{Int}},
    kernel_stats::Union{Nothing,BSKernelStats},
)
    nrecompute == 0 && return 0

    t_down = kernel_stats === nothing ? 0 : time_ns()
    @inbounds for q in 1:nrecompute
        t = recompute_slots[q]
        j = i + t
        sj = _tail_colnorm2(A, i + 1, m, j)
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

# Reflector construction goes through LAPACK dlarfg so Julia and the MEX
# backend share identical semantics, including dlarfg's rescaling
# safeguards near under/overflow (docs/VALIDATION.md, V5). The stdlib
# LAPACK.larfg! discards beta, hence the direct ccall.
function _householder!(x::StridedVector{Float64})
    n = LinearAlgebra.BlasInt(length(x))
    n >= 1 || return 0.0, 0.0

    alpha = Ref{Float64}(x[1])
    tau = Ref{Float64}(0.0)
    ccall(
        (BLAS.@blasfunc(dlarfg_), BLAS.libblastrampoline),
        Cvoid,
        (
            Ref{LinearAlgebra.BlasInt},
            Ref{Float64},
            Ptr{Float64},
            Ref{LinearAlgebra.BlasInt},
            Ref{Float64},
        ),
        n,
        alpha,
        pointer(x, 2),
        LinearAlgebra.BlasInt(1),
        tau,
    )
    beta = alpha[]
    x[1] = beta
    return tau[], beta
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

    # Apply A[i:m, i+1:n] -= tau * v * (v' * A[i:m, i+1:n]),
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
