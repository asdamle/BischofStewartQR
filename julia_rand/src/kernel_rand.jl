# The randomized selection kernel: a port of the MEX backend
# (matlab_rand/mex/src/bsqr_rand_mex.cpp), behaviourally in lockstep with the
# MATLAB implementations (same acceptance rule, threshold formulas, fallback
# semantics, and stats attribution) but with its own RNG -- pivot sequences
# need not match either MATLAB backend; only the guarantees are the
# cross-language contract.
#
# Two paths (batched, default true):
#   batched=true  - in-block BSQR: each sampled block is brought to the current
#                   frame ONCE (one compact-WY apply); BSQR then runs within
#                   it, selecting the in-block minimizer greedily -- f2 and the
#                   threshold updated after EACH pick -- while the per-step
#                   bound holds. One expensive apply yields many selections
#                   (O(k^3) overall vs the single-select O(k^4)). The in-block
#                   w-vectors and running residual norms are downdated
#                   incrementally (BLAS-2), with the Businger-Golub recompute
#                   safeguard (norm_recomp_tol) on the residual norms.
#   batched=false - single-select: one accepted column per sampled block;
#                   norms are recomputed from scratch per block (no safeguard
#                   needed).
# If sampling keeps missing acceptable columns (since_last >= remaining), the
# kernel scans ALL remaining columns for the true global minimizer -- which is
# guaranteed acceptable up to rounding ties -- rather than force-accepting a
# poor column a block happened to miss.

# Immutable bundle of kernel state; the mutable scalars (f2, pool sizes,
# counters) live as locals in the run functions.
struct _RandKernel
    A::Matrix{Float64}          # input; columns are gathered, never mutated
    m::Int
    n::Int
    k::Int
    block::Int
    threshold_mode::Symbol
    slack::Float64
    norm_recomp_tol::Float64
    weighted::Bool
    pick::Symbol
    # factorization state
    V::Matrix{Float64}          # m x k reflector store, unit-diagonal packed
    Twy::Matrix{Float64}        # k x k compact-WY factor
    R11::Matrix{Float64}        # k x k
    tau::Vector{Float64}
    selected::Vector{Int}
    taken::Vector{Bool}
    # per-step stats
    st_f2::Vector{Float64}
    st_crit::Vector{Float64}
    st_thr::Vector{Float64}
    st_Fhat::Vector{Float64}
    st_samples::Vector{Int}
    st_rounds::Vector{Int}
    st_fallback::Vector{Bool}
    # scratch (bmax = block columns)
    X::Matrix{Float64}          # m x bmax gathered/reduced candidates
    WB::Matrix{Float64}         # k x bmax compact-WY scratch
    Xtop::Matrix{Float64}       # k x bmax w-vector scratch (R11^{-1} X_top)
    cbuf::Vector{Float64}
    tcol::Vector{Float64}
    xcol::Vector{Float64}
    ids::Vector{Int}
    # sampling state
    g::Vector{Float64}          # original squared column norms (weighted only)
    bit::Fenwick
    drawn::Vector{Int}
    inpool::Vector{Bool}
    remaining::Vector{Int}      # uniform pool (valid front)
end

function _make_kernel(
    A::Matrix{Float64},
    k::Int,
    block::Int,
    threshold_mode::Symbol,
    slack::Float64,
    norm_recomp_tol::Float64,
    weighted::Bool,
    pick::Symbol,
)
    m, n = size(A)
    bmax = max(block, 1)
    g = Float64[]
    inpool = Bool[]
    remaining = Int[]
    if weighted
        g = Vector{Float64}(undef, n)
        @inbounds for j in 1:n
            nj = BLAS.nrm2(view(A, :, j))
            g[j] = nj * nj
        end
        # Free non-finite detection: any NaN/Inf in A makes its column's
        # squared norm non-finite, so this O(n) scan of the already-computed
        # sampling weights catches bad input even with check_finite=false (the
        # default). Uniform sampling has no such pass; there bad input needs
        # check_finite=true.
        all(isfinite, g) || throw(ArgumentError("A contains non-finite values"))
        inpool = fill(true, n)
    else
        remaining = collect(1:n)
    end
    bit = weighted ? Fenwick(g) : Fenwick(Float64[])
    return _RandKernel(
        A, m, n, k, block, threshold_mode, slack, norm_recomp_tol, weighted, pick,
        zeros(m, k), zeros(k, k), zeros(k, k), zeros(k), zeros(Int, k), fill(false, n),
        zeros(k), zeros(k), zeros(k), zeros(k), zeros(Int, k), zeros(Int, k),
        fill(false, k),
        Matrix{Float64}(undef, m, bmax), zeros(k, bmax), zeros(k, bmax),
        zeros(bmax), zeros(k), zeros(m), zeros(Int, bmax),
        g, bit, sizehint!(Int[], n), inpool, remaining,
    )
end

@inline _theta(st::_RandKernel, f2::Float64, nsel::Int) =
    _threshold(st.threshold_mode, f2, nsel, st.k, st.n) * st.slack

# Draw `bcount` distinct columns into X[:, 1:bcount], recording their original
# indices in `idbuf`. Weighted draws come off the Fenwick tree (removed until
# the block is resolved); uniform uses a partial Fisher-Yates over
# remaining[from:rcount]. Returns nothing; the weighted draws are logged in
# st.drawn for later restore.
function _draw_block!(
    st::_RandKernel,
    bcount::Int,
    idbuf::Vector{Int},
    rng::AbstractRNG,
    rcount::Int,
    from::Int,
)
    @inbounds for t in 1:bcount
        local id::Int
        if st.weighted
            id = _draw_weighted_one!(st.bit, st.g, st.inpool, rng)
            push!(st.drawn, id)
        else
            p = from + t - 1
            j = rand(rng, p:rcount)
            st.remaining[p], st.remaining[j] = st.remaining[j], st.remaining[p]
            id = st.remaining[p]
        end
        idbuf[t] = id
        copyto!(view(st.X, :, t), view(st.A, :, id))
    end
    return nothing
end

# Record a completed selection's stats at `step` (1-based).
@inline function _record_step!(
    st::_RandKernel,
    step::Int,
    f2::Float64,
    crit::Float64,
    theta::Float64,
    fallback::Bool,
)
    st.st_f2[step] = f2
    st.st_crit[step] = crit
    st.st_thr[step] = theta
    st.st_Fhat[step] = _fhat(step - 1, st.k, st.n)
    st.st_fallback[step] = fallback
    return nothing
end

# ===== Batched in-block BSQR (default kernel path) ==========================

function _run_batched!(st::_RandKernel, rng::AbstractRNG)
    m, n, k, block = st.m, st.n, st.k, st.block
    A, X, V, Twy, R11 = st.A, st.X, st.V, st.Twy, st.R11
    bmax = max(block, 1)
    # In-block state: w-vectors, running squared trailing norms (+ last exact
    # value, the recompute reference), ||w||^2, and downdate scratch.
    Wblk = zeros(k, bmax)
    sblk = zeros(bmax)
    sblk_ref = zeros(bmax)
    wn2blk = zeros(bmax)
    dots = zeros(bmax)      # w_a . w_pivot (pre-update)
    betav = zeros(bmax)     # alpha/diag per rest column
    hw = zeros(bmax)        # reflector-apply gemv scratch
    idblk = zeros(Int, bmax)

    f2 = 0.0
    nsel = 0
    rcount = n              # uniform pool size (valid front of `remaining`)
    since_last = 0          # columns sampled since the last selection
    rounds_since_last = 0   # blocks (reflector applies) since the last selection

    while nsel < k
        rem = n - nsel
        if since_last >= rem
            f2 = _global_min_fallback!(st, f2, nsel, since_last, rounds_since_last)
            if !st.weighted
                # Compact the uniform pool: drop every taken column.
                w = 0
                @inbounds for jj in 1:rcount
                    if !st.taken[st.remaining[jj]]
                        w += 1
                        st.remaining[w] = st.remaining[jj]
                    end
                end
                rcount = w
            end
            since_last = 0
            rounds_since_last = 0
            nsel += 1
            continue
        end
        bcount = min(block, rem)

        st.weighted && empty!(st.drawn)
        _draw_block!(st, bcount, idblk, rng, rcount, 1)
        since_last += bcount
        rounds_since_last += 1

        # Bring the block to the current frame, then seed the in-block state:
        # w-vectors with one trsm, running norms exactly.
        _wy_apply!(X, bcount, V, Twy, st.WB, nsel)
        if nsel > 0
            @inbounds for t in 1:bcount
                copyto!(view(Wblk, 1:nsel, t), view(X, 1:nsel, t))
            end
            ldiv!(UpperTriangular(view(R11, 1:nsel, 1:nsel)), view(Wblk, 1:nsel, 1:bcount))
        end
        @inbounds for t in 1:bcount
            s = _colnorm2(X, nsel + 1, m, t)
            sblk[t] = s
            sblk_ref[t] = s     # exact value seeds the recompute reference
            wn2blk[t] = nsel > 0 ? _colnorm2(Wblk, 1, nsel, t) : 0.0
        end

        nact = bcount
        block_start = nsel

        # --- in-block greedy BSQR ---
        while nsel < k && nact > 0
            pbest = 1
            cbest = Inf
            @inbounds for a in 1:nact
                c = sblk[a] > 0.0 ? (1.0 + wn2blk[a]) / sblk[a] : Inf
                if c < cbest
                    cbest = c
                    pbest = a
                end
            end
            theta = _theta(st, f2, nsel)
            (isfinite(cbest) && cbest <= theta) || break   # none acceptable -> resample

            # Swap the pivot to the end of the active region (rest stays contiguous).
            nr = nact - 1
            if pbest != nact
                @inbounds for r in 1:m
                    X[r, pbest], X[r, nact] = X[r, nact], X[r, pbest]
                end
                @inbounds for r in 1:k
                    Wblk[r, pbest], Wblk[r, nact] = Wblk[r, nact], Wblk[r, pbest]
                end
                sblk[pbest], sblk[nact] = sblk[nact], sblk[pbest]
                sblk_ref[pbest], sblk_ref[nact] = sblk_ref[nact], sblk_ref[pbest]
                wn2blk[pbest], wn2blk[nact] = wn2blk[nact], wn2blk[pbest]
                idblk[pbest], idblk[nact] = idblk[nact], idblk[pbest]
            end

            step = nsel + 1
            wcoeff = wn2blk[nact] + 1.0    # pivot ||w||^2 + 1, read before overwrite
            copyto!(st.xcol, view(X, :, nact))
            tau_i = _commit_reflector!(V, Twy, R11, st.tau, st.tcol, st.xcol, nsel, m)
            beta = st.xcol[step]
            invdiag = beta != 0.0 ? 1.0 / beta : 0.0

            # Apply the new reflector to the rest of the block (rows step:m --
            # it is zero above), then downdate their w-vectors / running norms
            # incrementally (BLAS-2 over nr contiguous columns).
            if nr > 0
                vtail = view(V, step:m, step)
                Xtail = view(X, step:m, 1:nr)
                hwv = view(hw, 1:nr)
                mul!(hwv, Xtail', vtail)
                BLAS.ger!(-tau_i, vtail, hwv, Xtail)
                @inbounds for a in 1:nr
                    betav[a] = X[step, a] * invdiag
                end
                if nsel > 0
                    wpiv = view(Wblk, 1:nsel, nact)
                    mul!(view(dots, 1:nr), view(Wblk, 1:nsel, 1:nr)', wpiv)
                    BLAS.ger!(-1.0, wpiv, view(betav, 1:nr), view(Wblk, 1:nsel, 1:nr))
                else
                    fill!(view(dots, 1:nr), 0.0)
                end
                @inbounds for a in 1:nr
                    bv = betav[a]
                    Wblk[step, a] = bv                                   # new w-row
                    wn2blk[a] += -2.0 * bv * dots[a] + bv * bv * wcoeff  # ||w||^2 downdate
                    alpha = X[step, a]
                    s = sblk[a] - alpha * alpha                          # running-norm downdate
                    sblk[a] = s > 0.0 ? s : 0.0

                    # Businger-Golub safeguard (as in the deterministic
                    # kernel): once the running residual norm has decayed past
                    # norm_recomp_tol * (last exact value), recompute it
                    # exactly from the just-updated block column (rows
                    # step+1:m, the post-elimination residual). Residual norm
                    # only -- wn2blk is refreshed exactly at every block
                    # boundary, so it is not separately guarded.
                    if sblk[a] <= sblk_ref[a] * st.norm_recomp_tol
                        s_exact = _colnorm2(X, step + 1, m, a)
                        sblk[a] = s_exact
                        sblk_ref[a] = s_exact
                    end
                end
            end

            # Commit the selection.
            f2 += cbest
            st.selected[step] = idblk[nact]
            st.taken[idblk[nact]] = true
            _record_step!(st, step, f2, cbest, theta, false)
            if nsel == block_start
                # ALL sampling work since the previous selection -- including
                # blocks that yielded no selection -- is attributed to a
                # block's first pick (0 for its later in-block picks), so
                # total_tested / blocks_sampled are true totals.
                st.st_samples[step] = since_last
                st.st_rounds[step] = rounds_since_last
            end
            since_last = 0
            rounds_since_last = 0
            nact = nr   # pivot removed
            nsel += 1
        end

        # Return the unselected drawn columns to the pool.
        if st.weighted
            @inbounds for id in st.drawn
                if !st.taken[id]
                    fenwick_add!(st.bit, id, st.g[id])
                    st.inpool[id] = true
                end
            end
        else
            # Only the drawn columns (remaining[1:bcount]) can have been taken
            # this block, so swap-pop the taken ones out of the pool in
            # O(bcount) rather than rescanning the whole pool. Refill each
            # hole from the tail and re-examine it -- the tail can itself be a
            # taken drawn column.
            t = 1
            @inbounds while t <= bcount && t <= rcount
                if st.taken[st.remaining[t]]
                    st.remaining[t] = st.remaining[rcount]
                    rcount -= 1
                else
                    t += 1
                end
            end
        end
    end
    return f2
end

# Global-min fallback for the batched path: random blocks keep missing the
# acceptable columns (e.g. uniform sampling on concentrated-leverage input),
# so scan ALL remaining columns in chunks for the true minimizer and take it.
# One selection; the caller resumes normal sampling. Errors when every
# remaining residual is ~zero (rank-deficient input). Returns the updated f2;
# the caller handles uniform-pool compaction and the nsel increment.
function _global_min_fallback!(
    st::_RandKernel,
    f2::Float64,
    nsel::Int,
    since_last::Int,
    rounds_since_last::Int,
)
    m, n, k, block = st.m, st.n, st.k, st.block
    A, X, V, Twy, R11 = st.A, st.X, st.V, st.Twy, st.R11
    theta = _theta(st, f2, nsel)
    gmin_c = Inf
    gmin_id = 0
    scanned = 0
    chunks = 0
    j0 = 1
    while j0 <= n
        bc = 0
        while j0 <= n && bc < block
            if !st.taken[j0]
                bc += 1
                copyto!(view(X, :, bc), view(A, :, j0))
                st.ids[bc] = j0
            end
            j0 += 1
        end
        bc == 0 && break
        chunks += 1
        scanned += bc
        if nsel > 0
            _wy_apply!(X, bc, V, Twy, st.WB, nsel)
            @inbounds for t in 1:bc
                copyto!(view(st.Xtop, 1:nsel, t), view(X, 1:nsel, t))
            end
            ldiv!(UpperTriangular(view(R11, 1:nsel, 1:nsel)), view(st.Xtop, 1:nsel, 1:bc))
        end
        @inbounds for t in 1:bc
            rho2 = _colnorm2(X, nsel + 1, m, t)
            wn2 = nsel > 0 ? _colnorm2(st.Xtop, 1, nsel, t) : 0.0
            c = rho2 > 0.0 ? (1.0 + wn2) / rho2 : Inf
            if c < gmin_c
                gmin_c = c
                gmin_id = st.ids[t]
                copyto!(st.xcol, view(X, :, t))
            end
        end
    end
    isfinite(gmin_c) || throw(RankDeficientError(nsel + 1, k))

    # Reduce the global minimizer (its reduced column is already in xcol).
    step = nsel + 1
    _commit_reflector!(V, Twy, R11, st.tau, st.tcol, st.xcol, nsel, m)
    f2 += gmin_c
    st.selected[step] = gmin_id
    st.taken[gmin_id] = true
    if st.weighted
        fenwick_add!(st.bit, gmin_id, -st.g[gmin_id])
        st.inpool[gmin_id] = false
    end
    _record_step!(st, step, f2, gmin_c, theta, gmin_c > theta)
    # Attribute ALL work since the previous selection: the failed random
    # blocks plus this exhaustive scan.
    st.st_samples[step] = since_last + scanned
    st.st_rounds[step] = rounds_since_last + chunks
    return f2
end

# ===== Single-select path (batched = false) =================================

function _run_single!(st::_RandKernel, rng::AbstractRNG)
    m, n, k, block = st.m, st.n, st.k, st.block
    A, X, V, Twy, R11 = st.A, st.X, st.V, st.Twy, st.R11
    f2 = 0.0
    rem_count = n

    for nsel in 0:(k - 1)
        step = nsel + 1
        theta = _theta(st, f2, nsel)

        st.weighted && empty!(st.drawn)
        best_c = Inf
        best_id = 0
        best_pos = 0            # position in `remaining` (uniform only)
        accepted = false
        accept_id = 0
        accept_local = 0        # column within the last examined block
        accept_pos = 0
        accept_c = Inf
        tested = 0
        rounds = 0

        pos = 1
        while pos <= rem_count && !accepted
            bcount = min(block, rem_count - pos + 1)
            # Uniform draws freeze positions pos:pos+bcount-1 of the pool
            # (partial Fisher-Yates), so blocks within a step never repeat a
            # column; weighted draws come off the tree until the step resolves.
            _draw_block!(st, bcount, st.ids, rng, rem_count, pos)

            _wy_apply!(X, bcount, V, Twy, st.WB, nsel)
            if nsel > 0
                @inbounds for t in 1:bcount
                    copyto!(view(st.Xtop, 1:nsel, t), view(X, 1:nsel, t))
                end
                ldiv!(
                    UpperTriangular(view(R11, 1:nsel, 1:nsel)),
                    view(st.Xtop, 1:nsel, 1:bcount),
                )
            end

            blk_best = 1
            blk_best_c = Inf
            @inbounds for t in 1:bcount
                rho2 = _colnorm2(X, nsel + 1, m, t)
                wn2 = nsel > 0 ? _colnorm2(st.Xtop, 1, nsel, t) : 0.0
                c = rho2 > 0.0 ? (1.0 + wn2) / rho2 : Inf
                st.cbuf[t] = c
                if c < blk_best_c
                    blk_best_c = c
                    blk_best = t
                end
            end
            tested += bcount
            rounds += 1
            if blk_best_c < best_c
                best_c = blk_best_c
                best_id = st.ids[blk_best]
                best_pos = pos + blk_best - 1
            end

            if st.pick === :first
                @inbounds for t in 1:bcount
                    if st.cbuf[t] <= theta
                        accept_id = st.ids[t]
                        accept_c = st.cbuf[t]
                        accept_local = t
                        accept_pos = pos + t - 1
                        accepted = true
                        break
                    end
                end
            elseif blk_best_c <= theta
                accept_id = st.ids[blk_best]
                accept_c = blk_best_c
                accept_local = blk_best
                accept_pos = pos + blk_best - 1
                accepted = true
            end
            pos += bcount
        end

        fallback = false
        if !accepted
            accept_id = best_id
            accept_c = best_c
            accept_pos = best_pos
            fallback = true
        end

        # Rank guard: a non-finite increment means every remaining column had
        # a ~zero residual (rho^2 = 0).
        isfinite(accept_c) || throw(RankDeficientError(step, k))

        # Reduce the accepted column. If it was accepted from the last block
        # its reduced form is still in X (column accept_local); otherwise
        # re-gather and re-apply.
        if accepted
            copyto!(st.xcol, view(X, :, accept_local))
        else
            copyto!(st.xcol, view(A, :, accept_id))
            _wy_apply_col!(st.xcol, V, Twy, st.tcol, nsel)
        end
        _commit_reflector!(V, Twy, R11, st.tau, st.tcol, st.xcol, nsel, m)

        f2 += accept_c
        st.selected[step] = accept_id
        st.taken[accept_id] = true

        # Remove the pivot from the pool; restore the other drawn weights.
        if st.weighted
            @inbounds for id in st.drawn
                if id != accept_id
                    fenwick_add!(st.bit, id, st.g[id])
                    st.inpool[id] = true
                end
            end
            # accept_id stays removed (taken).
        else
            st.remaining[accept_pos] = st.remaining[rem_count]   # O(1) swap-pop
        end
        rem_count -= 1

        _record_step!(st, step, f2, accept_c, theta, fallback)
        st.st_samples[step] = tested
        st.st_rounds[step] = rounds
    end
    return f2
end
