mutable struct BSLazyBLASWorkspace
    s_upper::Vector{Float64}
    wnorm2_lb::Vector{Float64}
    exact_score::Vector{Float64}
    exact_s::Vector{Float64}
    exact_wnorm2::Vector{Float64}
    refreshed::BitVector
    candidates::Vector{Int}
    solve_block::Matrix{Float64}
    work::Vector{Float64}
end

mutable struct BSLazyBLASStats
    refresh_count::Int
    full_refresh_rounds::Int
end

BSLazyBLASStats() = BSLazyBLASStats(0, 0)

function BSLazyBLASWorkspace(_m::Int, n::Int, kmax::Int)
    return BSLazyBLASWorkspace(
        zeros(Float64, n),
        zeros(Float64, n),
        fill(Inf, n),
        zeros(Float64, n),
        zeros(Float64, n),
        falses(n),
        zeros(Int, n),
        zeros(Float64, max(kmax, 1), max(n, 1)),
        zeros(Float64, max(n, 1)),
    )
end

function _require_lazy_workspace(
    workspace::Union{Nothing,BSLazyBLASWorkspace},
    m::Int,
    n::Int,
    kmax::Int,
)
    if workspace === nothing
        return BSLazyBLASWorkspace(m, n, kmax)
    end

    ws = workspace
    length(ws.s_upper) >= n || throw(ArgumentError("lazy workspace.s_upper too small"))
    length(ws.wnorm2_lb) >= n || throw(ArgumentError("lazy workspace.wnorm2_lb too small"))
    length(ws.exact_score) >= n || throw(ArgumentError("lazy workspace.exact_score too small"))
    length(ws.exact_s) >= n || throw(ArgumentError("lazy workspace.exact_s too small"))
    length(ws.exact_wnorm2) >= n || throw(ArgumentError("lazy workspace.exact_wnorm2 too small"))
    length(ws.refreshed) >= n || throw(ArgumentError("lazy workspace.refreshed too small"))
    length(ws.candidates) >= n || throw(ArgumentError("lazy workspace.candidates too small"))
    size(ws.solve_block, 1) >= max(kmax, 1) ||
        throw(ArgumentError("lazy workspace.solve_block has too few rows"))
    size(ws.solve_block, 2) >= max(n, 1) ||
        throw(ArgumentError("lazy workspace.solve_block has too few columns"))
    length(ws.work) >= max(n, 1) || throw(ArgumentError("lazy workspace.work too small"))
    return ws
end

@inline function _lazy_score_lower_bound(s_upper::Float64, wnorm2_lb::Float64)
    return s_upper > 0.0 ? (1.0 + wnorm2_lb) / s_upper : Inf
end

function _lazy_refresh_batch_size(
    remaining::Int,
    batch_size::Union{Nothing,Int},
    batch_fraction::Float64,
    batch_min::Int,
    batch_max::Union{Nothing,Int},
)
    remaining > 0 || return 0
    if batch_size !== nothing
        return min(remaining, batch_size)
    end
    batch = ceil(Int, batch_fraction * remaining)
    batch = max(batch, batch_min, 1)
    if batch_max !== nothing
        batch = min(batch, batch_max)
    end
    return min(batch, remaining)
end

function _lazy_collect_candidates!(
    ws::BSLazyBLASWorkspace,
    i::Int,
    n::Int,
    batch_count::Int,
)
    total = 0
    @inbounds for j in i:n
        if !ws.refreshed[j]
            total += 1
            ws.candidates[total] = j
        end
    end
    sort!(
        view(ws.candidates, 1:total);
        by = j -> (_lazy_score_lower_bound(ws.s_upper[j], ws.wnorm2_lb[j]), j),
    )
    return min(batch_count, total), total
end

function _lazy_refresh_exact_batch!(
    A::StridedMatrix{Float64},
    ws::BSLazyBLASWorkspace,
    i::Int,
    batch_count::Int,
    stats::Union{Nothing,BSLazyBLASStats},
)
    m, _ = size(A)
    batch_count > 0 || return nothing

    if i > 1
        rows = i - 1
        block = view(ws.solve_block, 1:rows, 1:batch_count)
        @inbounds for t in 1:batch_count
            j = ws.candidates[t]
            copyto!(view(block, :, t), view(A, 1:rows, j))
        end
        ldiv!(UpperTriangular(view(A, 1:rows, 1:rows)), block)
        @inbounds for t in 1:batch_count
            j = ws.candidates[t]
            sj = BLAS.nrm2(view(A, i:m, j))
            sj *= sj
            wn = BLAS.nrm2(view(block, :, t))
            wn *= wn
            ws.exact_s[j] = sj
            ws.exact_wnorm2[j] = wn
            ws.exact_score[j] = sj > 0.0 ? (1.0 + wn) / sj : Inf
            ws.refreshed[j] = true
        end
    else
        @inbounds for t in 1:batch_count
            j = ws.candidates[t]
            sj = BLAS.nrm2(view(A, :, j))
            sj *= sj
            ws.exact_s[j] = sj
            ws.exact_wnorm2[j] = 0.0
            ws.exact_score[j] = sj > 0.0 ? 1.0 / sj : Inf
            ws.refreshed[j] = true
        end
    end

    if stats !== nothing
        stats.refresh_count += batch_count
    end
    return nothing
end

function _lazy_choose_pivot!(
    A::StridedMatrix{Float64},
    ws::BSLazyBLASWorkspace,
    i::Int,
    batch_size::Union{Nothing,Int},
    batch_fraction::Float64,
    batch_min::Int,
    batch_max::Union{Nothing,Int},
    stats::Union{Nothing,BSLazyBLASStats},
)
    _, n = size(A)

    while true
        best_j = 0
        best_c = Inf
        @inbounds for j in i:n
            if ws.refreshed[j]
                cj = ws.exact_score[j]
                if cj < best_c
                    best_c = cj
                    best_j = j
                end
            end
        end

        min_lb = Inf
        remaining = 0
        @inbounds for j in i:n
            if !ws.refreshed[j]
                remaining += 1
                lb = _lazy_score_lower_bound(ws.s_upper[j], ws.wnorm2_lb[j])
                if lb < min_lb
                    min_lb = lb
                end
            end
        end

        if best_j != 0 && (remaining == 0 || best_c < min_lb)
            if stats !== nothing && remaining == 0
                stats.full_refresh_rounds += 1
            end
            return best_j, best_c
        end

        remaining > 0 || error("lazy BSQR could not certify a pivot")
        want = _lazy_refresh_batch_size(remaining, batch_size, batch_fraction, batch_min, batch_max)
        batch_count, _ = _lazy_collect_candidates!(ws, i, n, want)
        _lazy_refresh_exact_batch!(A, ws, i, batch_count, stats)
    end
end

function _swap_lazy_state!(ws::BSLazyBLASWorkspace, i::Int, j::Int)
    _swap_entries!(ws.s_upper, i, j)
    _swap_entries!(ws.wnorm2_lb, i, j)
    _swap_entries!(ws.exact_score, i, j)
    _swap_entries!(ws.exact_s, i, j)
    _swap_entries!(ws.exact_wnorm2, i, j)
    ws.refreshed[i], ws.refreshed[j] = ws.refreshed[j], ws.refreshed[i]
    return nothing
end

function _validate_lazy_batch_args(
    batch_size::Union{Nothing,Integer},
    batch_fraction::Float64,
    batch_min::Integer,
    batch_max::Union{Nothing,Integer},
)
    if batch_size !== nothing
        Int(batch_size) >= 1 || throw(ArgumentError("lazy_batch_size must be >= 1"))
    end
    (0.0 <= batch_fraction <= 1.0) ||
        throw(ArgumentError("lazy_batch_fraction must satisfy 0 <= value <= 1"))
    Int(batch_min) >= 1 || throw(ArgumentError("lazy_batch_min must be >= 1"))
    if batch_max !== nothing
        Int(batch_max) >= 1 || throw(ArgumentError("lazy_batch_max must be >= 1"))
    end
    return Int(batch_size === nothing ? 0 : batch_size),
        batch_fraction,
        Int(batch_min),
        batch_max === nothing ? nothing : Int(batch_max)
end

function _bsqr_lazy_blas_kernel!(
    A::StridedMatrix{Float64},
    tau::Vector{Float64},
    jpvt::Vector{Int},
    ws::BSLazyBLASWorkspace,
    k::Int;
    pivot_history::Union{Nothing,Vector{Int}} = nothing,
    rank_stop::Bool = true,
    lazy_batch_size::Union{Nothing,Integer} = nothing,
    lazy_batch_fraction::Float64 = 0.125,
    lazy_batch_min::Integer = 4,
    lazy_batch_max::Union{Nothing,Integer} = nothing,
    selection_stats::Union{Nothing,BSLazyBLASStats} = nothing,
)
    m, n = size(A)
    batch_size_int, batch_fraction, batch_min, batch_max = _validate_lazy_batch_args(
        lazy_batch_size,
        lazy_batch_fraction,
        lazy_batch_min,
        lazy_batch_max,
    )
    batch_size = lazy_batch_size === nothing ? nothing : batch_size_int
    k == 0 && return 0
    tol_scale = eps(Float64) * max(m, n)

    @inbounds for j in 1:n
        sj = BLAS.nrm2(view(A, :, j))
        sj *= sj
        ws.s_upper[j] = sj
        ws.wnorm2_lb[j] = 0.0
        ws.exact_score[j] = Inf
        ws.exact_s[j] = 0.0
        ws.exact_wnorm2[j] = 0.0
        ws.refreshed[j] = false
    end

    ksteps = 0
    for i in 1:k
        best_j, best_c = _lazy_choose_pivot!(
            A,
            ws,
            i,
            batch_size,
            batch_fraction,
            batch_min,
            batch_max,
            selection_stats,
        )

        if best_j != i
            _swap_columns!(A, i, best_j)
            _swap_lazy_state!(ws, i, best_j)
            _swap_entries!(jpvt, i, best_j)
        end

        if pivot_history !== nothing
            push!(pivot_history, jpvt[i])
        end

        pivot_norm2 = ws.exact_s[i]
        atol = tol_scale * max(pivot_norm2, 1.0)
        if !(pivot_norm2 > atol)
            if rank_stop
                tau[i] = 0.0
                break
            end
        end

        tau_i, beta_i = _householder!(view(A, i:m, i))
        tau[i] = tau_i

        if tau_i != 0.0 && i < n
            A[i, i] = 1.0
            _apply_householder_left!(A, i, tau_i, ws.work)
        end
        A[i, i] = beta_i

        ws.s_upper[i] = beta_i * beta_i
        ws.wnorm2_lb[i] = ws.exact_wnorm2[i]
        ws.exact_score[i] = best_c

        denom_factor = 1.0 + ws.exact_wnorm2[i]
        @inbounds for j in (i + 1):n
            if ws.refreshed[j]
                ws.wnorm2_lb[j] = ws.exact_wnorm2[j] / denom_factor
                ws.s_upper[j] = ws.exact_s[j]
            else
                ws.wnorm2_lb[j] /= denom_factor
            end
            ws.refreshed[j] = false
            ws.exact_score[j] = Inf
            ws.exact_s[j] = 0.0
            ws.exact_wnorm2[j] = 0.0
        end

        ksteps = i
    end

    return ksteps
end

function _bsqr_lazy_blas!(
    A::StridedMatrix{Float64};
    k::Integer = min(size(A)...),
    check::Bool = true,
    rank_stop::Bool = true,
    blas_threads::Union{Nothing,Integer} = nothing,
    workspace::Union{Nothing,BSLazyBLASWorkspace} = nothing,
    lazy_batch_size::Union{Nothing,Integer} = nothing,
    lazy_batch_fraction::Float64 = 0.125,
    lazy_batch_min::Integer = 4,
    lazy_batch_max::Union{Nothing,Integer} = nothing,
    pivot_history::Union{Nothing,Vector{Int}} = nothing,
    selection_stats::Union{Nothing,BSLazyBLASStats} = nothing,
)
    m, n, kwork = _validate_bsqr_common_args(A, k, check, _DEFAULT_NORM_RECOMP_TOL)
    ws = _require_lazy_workspace(workspace, m, n, kwork)
    tau = zeros(Float64, kwork)
    jpvt = collect(1:n)

    ksteps = _with_blas_threads(blas_threads) do
        _bsqr_lazy_blas_kernel!(
            A,
            tau,
            jpvt,
            ws,
            kwork;
            pivot_history = pivot_history,
            rank_stop = rank_stop,
            lazy_batch_size = lazy_batch_size,
            lazy_batch_fraction = lazy_batch_fraction,
            lazy_batch_min = lazy_batch_min,
            lazy_batch_max = lazy_batch_max,
            selection_stats = selection_stats,
        )
    end

    factors = A isa Matrix{Float64} ? A : Matrix(A)
    return BSQRPivoted(factors, tau, jpvt, ksteps, nothing, nothing)
end

function _bsqr_lazy_blas(
    A::AbstractMatrix{<:Real};
    k::Integer = min(size(A)...),
    check::Bool = true,
    rank_stop::Bool = true,
    blas_threads::Union{Nothing,Integer} = nothing,
    workspace::Union{Nothing,BSLazyBLASWorkspace} = nothing,
    lazy_batch_size::Union{Nothing,Integer} = nothing,
    lazy_batch_fraction::Float64 = 0.125,
    lazy_batch_min::Integer = 4,
    lazy_batch_max::Union{Nothing,Integer} = nothing,
    pivot_history::Union{Nothing,Vector{Int}} = nothing,
    selection_stats::Union{Nothing,BSLazyBLASStats} = nothing,
)
    Acopy = Matrix{Float64}(A)
    return _bsqr_lazy_blas!(
        Acopy;
        k = k,
        check = check,
        rank_stop = rank_stop,
        blas_threads = blas_threads,
        workspace = workspace,
        lazy_batch_size = lazy_batch_size,
        lazy_batch_fraction = lazy_batch_fraction,
        lazy_batch_min = lazy_batch_min,
        lazy_batch_max = lazy_batch_max,
        pivot_history = pivot_history,
        selection_stats = selection_stats,
    )
end
