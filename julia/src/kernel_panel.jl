# Panel/blocked BSQR kernel (P3 prototype; see docs/P3_BLOCKED_BSQR.md).
#
# dlaqps-style reorganization of the unblocked kernel: within a panel of nb
# steps, the trailing block of A and the prefix rows of W receive no per-step
# rank-1 writes. Each step does two read-only gemvs (A'v for the deferred-
# update factors Fa, W0'omega for the criterion dots) plus O(nb x nrem)
# corrections; the deferred rank-nb updates are applied as two gemms at panel
# end. All selection quantities equal the unblocked ones in exact arithmetic
# (regrouped sums only); pivot parity is enforced by the fixture suite.
#
# The panel kernel is the default (nb = 8, the measured optimum; see
# docs/P3_BLOCKED_BSQR.md). BS_PANEL_NB overrides the width; 0 or 1
# selects the unblocked reference kernel.

struct BSPanelWorkspace
    V::Matrix{Float64}      # m x nb reflector vectors (panel-local rows)
    Fa::Matrix{Float64}     # n x nb deferred-update factors for A
    B::Matrix{Float64}      # nb x n panel rows of W (eagerly maintained)
    Bsel::Matrix{Float64}   # nb x n beta rows frozen at selection time
    Omega::Matrix{Float64}  # kmax x nb pivot prefix w-vectors
    h::Vector{Float64}      # nb scratch
    q::Vector{Float64}      # nb scratch
    g::Vector{Float64}      # n scratch (criterion dots)
    y::Vector{Float64}      # n scratch (Fa gemv)
    flagged::Vector{Bool}   # columns awaiting exact norm refresh at flush
    flush_list::Vector{Int}
end

function BSPanelWorkspace(m::Int, n::Int, kmax::Int, nb::Int)
    return BSPanelWorkspace(
        zeros(m, nb),
        zeros(n, nb),
        zeros(nb, n),
        zeros(nb, n),
        zeros(max(kmax, 1), nb),
        zeros(nb),
        zeros(nb),
        zeros(n),
        zeros(n),
        fill(false, n),
        Int[],
    )
end

const _DEFAULT_PANEL_NB = 8
# Empirical crossover: below this k*n the panel bookkeeping (per-call buffers,
# extra BLAS-call latency) outweighs the gemm benefit. Measured boundary on
# Apple Silicon/Accelerate: losers <= 16k (128x128, 32x256), winners >= 33k
# (128x256, 192x192, 64x640).
const _DEFAULT_PANEL_MIN_KN = 24576

@inline function _panel_nb()
    raw = strip(get(ENV, "BS_PANEL_NB", ""))
    isempty(raw) && return _DEFAULT_PANEL_NB
    v = tryparse(Int, raw)
    v === nothing && throw(ArgumentError("BS_PANEL_NB must be an integer"))
    return v
end

@inline function _panel_min_kn()
    raw = strip(get(ENV, "BS_PANEL_MIN_KN", ""))
    isempty(raw) && return _DEFAULT_PANEL_MIN_KN
    v = tryparse(Int, raw)
    v === nothing && throw(ArgumentError("BS_PANEL_MIN_KN must be an integer"))
    return v
end

# Swap every column-indexed structure for the pivot exchange i <-> best_j.
# W prefix rows are 1..s-1 (panel rows live in B); panel bookkeeping is
# indexed by absolute column position.
@inline function _panel_swap!(
    A::StridedMatrix{Float64},
    ws::BSWorkspace,
    pws::BSPanelWorkspace,
    jpvt::Vector{Int},
    s::Int,
    i::Int,
    best_j::Int,
    t::Int,
)
    best_j == i && return nothing
    _swap_columns!(A, i, best_j)
    if s > 1
        _swap_columns_prefix!(ws.W, i, best_j, s - 1)
    end
    _swap_entries!(ws.s, i, best_j)
    _swap_entries!(ws.s_ref, i, best_j)
    _swap_entries!(ws.wnorm2, i, best_j)
    _swap_entries!(jpvt, i, best_j)
    @inbounds for u in 1:(t - 1)
        pws.Fa[i, u], pws.Fa[best_j, u] = pws.Fa[best_j, u], pws.Fa[i, u]
        pws.B[u, i], pws.B[u, best_j] = pws.B[u, best_j], pws.B[u, i]
        pws.Bsel[u, i], pws.Bsel[u, best_j] = pws.Bsel[u, best_j], pws.Bsel[u, i]
    end
    return nothing
end

function _bsqr_kernel_panel!(
    A::StridedMatrix{Float64},
    tau::Vector{Float64},
    jpvt::Vector{Int},
    ws::BSWorkspace,
    k::Int,
    nb::Int;
    frob_inv_trace::Union{Nothing,Vector{Float64}} = nothing,
    pivot_history::Union{Nothing,Vector{Int}} = nothing,
    rank_stop::Bool = true,
    norm_recomp_tol::Float64 = _DEFAULT_NORM_RECOMP_TOL,
    norm_recomp_count::Union{Nothing,Base.RefValue{Int}} = nothing,
    recomp_history::Union{Nothing,Vector{Int}} = nothing,
)
    m, n = size(A)
    k == 0 && return 0
    nb = max(nb, 2)
    tol_scale = eps(Float64) * max(m, n)
    _init_running_colnorms!(A, ws, n)
    pws = BSPanelWorkspace(m, n, k, nb)

    ksteps = 0
    s = 1
    while s <= k
        tmax = min(nb, k - s + 1)
        tdone = 0
        stopped = false

        for t in 1:tmax
            i = s + t - 1
            best_j, best_c = _select_pivot_column!(ws, i, n, nothing)
            _panel_swap!(A, ws, pws, jpvt, s, i, best_j, t)
            if pivot_history !== nothing
                push!(pivot_history, jpvt[i])
            end

            if _rank_stop_triggered(ws, i, tol_scale) && rank_stop
                tau[i] = 0.0
                stopped = true
                break
            end

            # --- Catch up the pivot column of A and reflect ---
            # Rows s..i-1 of column i were already finalized by the per-step
            # row writes; only the still-stale rows i:m carry deferred updates.
            if t > 1
                BLAS.gemv!(
                    'N', -1.0,
                    view(pws.V, i:m, 1:(t - 1)),
                    view(pws.Fa, i, 1:(t - 1)),
                    1.0,
                    view(A, i:m, i),
                )
            end
            tau_i, beta_i = _householder!(view(A, i:m, i))
            tau[i] = tau_i
            fill!(view(pws.V, s:(i - 1), t), 0.0)
            pws.V[i, t] = 1.0
            if i < m
                copyto!(view(pws.V, (i + 1):m, t), view(A, (i + 1):m, i))
            end

            nrem = n - i
            if nrem > 0
                # --- Fa(:,t) = tau * (A_true' v) over stale rows i:m ---
                yv = view(pws.y, (i + 1):n)
                BLAS.gemv!(
                    'T', tau_i,
                    view(A, i:m, (i + 1):n),
                    view(pws.V, i:m, t),
                    0.0,
                    yv,
                )
                if t > 1
                    hv = view(pws.h, 1:(t - 1))
                    BLAS.gemv!(
                        'T', 1.0,
                        view(pws.V, i:m, 1:(t - 1)),
                        view(pws.V, i:m, t),
                        0.0,
                        hv,
                    )
                    BLAS.gemv!(
                        'N', -tau_i,
                        view(pws.Fa, (i + 1):n, 1:(t - 1)),
                        hv,
                        1.0,
                        yv,
                    )
                end
                copyto!(view(pws.Fa, (i + 1):n, t), yv)

                # --- Finalize row i of the trailing block ---
                alpha = view(A, i, (i + 1):n)
                BLAS.gemv!(
                    'N', -1.0,
                    view(pws.Fa, (i + 1):n, 1:t),
                    view(pws.V, i, 1:t),
                    1.0,
                    alpha,
                )
            end

            A[i, i] = beta_i
            ws.s[i] = beta_i * beta_i
            ws.s_ref[i] = ws.s[i]

            if nrem > 0
                # --- W-side scalars from the deferred representation ---
                invdiag = beta_i != 0.0 ? (1.0 / beta_i) : 0.0
                alpha = view(A, i, (i + 1):n)
                beta_vec = view(ws.beta, 1:nrem)
                @inbounds for r in 1:nrem
                    beta_vec[r] = alpha[r] * invdiag
                end

                # omega_t: true prefix w of the pivot (rows 1..s-1)
                omega = view(pws.Omega, 1:(s - 1), t)
                copyto!(omega, view(ws.W, 1:(s - 1), i))
                if t > 1 && s > 1
                    BLAS.gemv!(
                        'N', -1.0,
                        view(pws.Omega, 1:(s - 1), 1:(t - 1)),
                        view(pws.Bsel, 1:(t - 1), i),
                        1.0,
                        omega,
                    )
                end
                qv = view(pws.q, 1:(t - 1))
                copyto!(qv, view(pws.B, 1:(t - 1), i))

                wcoeff = ws.wnorm2[i] + 1.0
                dots = view(ws.dots, 1:nrem)
                # d = W0' omega - Bsel' (Omega' omega) + B' q
                if s > 1
                    BLAS.gemv!(
                        'T', 1.0,
                        view(ws.W, 1:(s - 1), (i + 1):n),
                        omega,
                        0.0,
                        dots,
                    )
                else
                    fill!(dots, 0.0)
                end
                if t > 1
                    hv = view(pws.h, 1:(t - 1))
                    if s > 1
                        BLAS.gemv!(
                            'T', 1.0,
                            view(pws.Omega, 1:(s - 1), 1:(t - 1)),
                            omega,
                            0.0,
                            hv,
                        )
                        BLAS.gemv!(
                            'T', -1.0,
                            view(pws.Bsel, 1:(t - 1), (i + 1):n),
                            hv,
                            1.0,
                            dots,
                        )
                    end
                    BLAS.gemv!(
                        'T', 1.0,
                        view(pws.B, 1:(t - 1), (i + 1):n),
                        qv,
                        1.0,
                        dots,
                    )
                end

                # --- Eager panel-row updates of W (rows s..i-1 live in B) ---
                copyto!(view(pws.Bsel, t, (i + 1):n), beta_vec)
                if t > 1
                    BLAS.ger!(-1.0, qv, beta_vec, view(pws.B, 1:(t - 1), (i + 1):n))
                end
                copyto!(view(pws.B, t, (i + 1):n), beta_vec)

                # --- wnorm2 recurrence and norm downdates (as unblocked) ---
                # Columns whose downdate trips the guard are batched and
                # refreshed exactly at panel flush: a flagged column's running
                # s is <= tol * s_ref, so its criterion is far too large to
                # win a selection within the remaining <= nb-1 panel steps,
                # and the refreshed value is the same mathematical quantity
                # either way. (Per-step early flushing degenerates panels to
                # width ~1 in recompute-heavy regimes; measured in
                # docs/P3_BLOCKED_BSQR.md.)
                nrecompute = 0
                @inbounds for r in 1:nrem
                    j = i + r
                    b = beta_vec[r]
                    wn = ws.wnorm2[j] - 2.0 * b * dots[r] + b * b * wcoeff
                    ws.wnorm2[j] = wn > 0.0 ? wn : 0.0

                    old_s = ws.s[j]
                    if !(old_s > 0.0)
                        ws.s[j] = 0.0
                        continue
                    end
                    sj = old_s - alpha[r] * alpha[r]
                    sj = sj > 0.0 ? sj : 0.0
                    ws.s[j] = sj
                    if sj <= ws.s_ref[j] * norm_recomp_tol && !pws.flagged[j]
                        pws.flagged[j] = true
                        push!(pws.flush_list, j)
                        nrecompute += 1
                    end
                end
                if recomp_history !== nothing
                    push!(recomp_history, nrecompute)
                end
                if nrecompute > 0 && norm_recomp_count !== nothing
                    norm_recomp_count[] += nrecompute
                end
            elseif recomp_history !== nothing
                push!(recomp_history, 0)
            end

            tdone = t
            ksteps = i
            if frob_inv_trace !== nothing
                push!(frob_inv_trace, best_c)
            end
        end

        # --- Panel flush: apply the deferred rank-tdone updates ---
        if tdone > 0
            ie = s + tdone - 1
            if ie < m && ie < n
                BLAS.gemm!(
                    'N', 'T', -1.0,
                    view(pws.V, (ie + 1):m, 1:tdone),
                    view(pws.Fa, (ie + 1):n, 1:tdone),
                    1.0,
                    view(A, (ie + 1):m, (ie + 1):n),
                )
            end
            if ie < n
                if s > 1
                    BLAS.gemm!(
                        'N', 'N', -1.0,
                        view(pws.Omega, 1:(s - 1), 1:tdone),
                        view(pws.Bsel, 1:tdone, (ie + 1):n),
                        1.0,
                        view(ws.W, 1:(s - 1), (ie + 1):n),
                    )
                end
                copyto!(view(ws.W, s:ie, (ie + 1):n), view(pws.B, 1:tdone, (ie + 1):n))
            end

            # Exact norm refresh for columns flagged mid-panel (now true).
            for j in pws.flush_list
                sj = _tail_colnorm2(A, ie + 1, m, j)
                ws.s_ref[j] = sj
                ws.s[j] = sj
                pws.flagged[j] = false
            end
            empty!(pws.flush_list)
        end

        stopped && break
        s += max(tdone, 1)
    end

    return ksteps
end
