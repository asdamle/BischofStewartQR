# Kernel helpers for the randomized Bischof-Stewart column selection.
#
# These are deliberate own copies, not imports from BSPivotQR: the randomized
# variant is kept decoupled from the deterministic kernels (repo CLAUDE.md), so
# a change here can never silently alter the tested deterministic code paths
# (and vice versa). The ports mirror matlab_rand/mex/src/bsqr_rand_mex.cpp,
# the performance reference implementation.

const _DEFAULT_NORM_RECOMP_TOL = sqrt(eps(Float64))

"""
    RankDeficientError(step, k)

Thrown when, at selection `step`, every remaining column has a ~zero residual
(`rho^2 == 0`): the input is rank-deficient for the requested `k`. The
per-step bound cannot be maintained past the numerical rank, so the kernel
fails loudly rather than propagate `Inf` into `f2`/`R11` (randomized BSQR
targets the full-rank GKS setting; the deterministic kernel's `rank_stop` has
no analogue here).
"""
struct RankDeficientError <: Exception
    step::Int
    k::Int
end

function Base.showerror(io::IO, e::RankDeficientError)
    print(io, "RankDeficientError: at step ", e.step,
        " all remaining columns have ~zero residual; the input appears ",
        "rank-deficient for k=", e.k, ". Reduce k to the numerical rank.")
end

# Per-step acceptance threshold on the criterion increment
# c = (1 + ||w||^2) / rho^2 to f2 = ||R11^{-1}||_F^2 (port of
# matlab_rand/private/bsqr_rand_threshold.m; see docs/RANDOMIZED_BSQR_PLAN.md).
@inline function _threshold(mode::Symbol, f2::Float64, nsel::Int, k::Int, n::Int)
    if mode === :running_mean
        # rho^2-weighted mean of c over the remaining columns (orthonormal-row
        # input), so the minimizer always qualifies; per-singular-value control.
        return (f2 + n - 2 * nsel) / (k - nsel)
    end
    # :worstcase_allowance -- spend accumulated slack against the deterministic
    # worst case Fhat_{nsel+1}; more permissive, same final Frobenius bound.
    return (nsel + 1.0) * (n - nsel) / (k - nsel) - f2
end

@inline _fhat(nsel::Int, k::Int, n::Int) = (nsel + 1.0) * (n - nsel) / (k - nsel)

# Reflector construction through LAPACK dlarfg, so the reflector semantics
# (including rescaling safeguards near under/overflow) match the MEX backend
# exactly. Deliberate copy of BSPivotQR's _householder! (decoupling rule); the
# stdlib LAPACK.larfg! discards beta, hence the direct ccall.
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

# Squared norm of X[r1:r2, col] (accumulated in-order with fma, matching the
# MEX's scalar loops).
@inline function _colnorm2(X::AbstractMatrix{Float64}, r1::Int, r2::Int, col::Int)
    acc = 0.0
    @inbounds @simd for r in r1:r2
        v = X[r, col]
        acc = muladd(v, v, acc)
    end
    return acc
end

# --- compact-WY application of the accumulated reflectors -------------------
# Q_nsel' = I - V T' V' with V the unit-diagonal packed reflector store and T
# the upper-triangular compact-WY factor (dlarft forward recurrence). The
# block apply is one gemm + trmm + gemm -- BLAS-3, the reason the k<<n win
# materializes (see the MEX header comment).

# X[:, 1:bc] <- Q_nsel' * X[:, 1:bc]; WB is nsel x bc scratch.
function _wy_apply!(
    X::AbstractMatrix{Float64},
    bc::Int,
    V::Matrix{Float64},
    Twy::Matrix{Float64},
    WB::Matrix{Float64},
    nsel::Int,
)
    (nsel == 0 || bc == 0) && return nothing
    Vv = view(V, :, 1:nsel)
    Xv = view(X, :, 1:bc)
    Wv = view(WB, 1:nsel, 1:bc)
    mul!(Wv, Vv', Xv)
    lmul!(adjoint(UpperTriangular(view(Twy, 1:nsel, 1:nsel))), Wv)
    mul!(Xv, Vv, Wv, -1.0, 1.0)
    return nothing
end

# Single-column variant: x <- Q_nsel' * x; tcol is length-k scratch.
function _wy_apply_col!(
    x::Vector{Float64},
    V::Matrix{Float64},
    Twy::Matrix{Float64},
    tcol::Vector{Float64},
    nsel::Int,
)
    nsel == 0 && return nothing
    tv = view(tcol, 1:nsel)
    mul!(tv, view(V, :, 1:nsel)', x)
    lmul!(adjoint(UpperTriangular(view(Twy, 1:nsel, 1:nsel))), tv)
    mul!(x, view(V, :, 1:nsel), tv, -1.0, 1.0)
    return nothing
end

# Reduce the (already frame-reduced) column xcol into the factorization at
# step nsel+1: write the R11 column, build the new reflector (xcol is
# overwritten: xcol[step] = beta, tail = reflector tail), append it to the
# reflector store, and extend the compact-WY T factor (dlarft forward
# recurrence: T[1:nsel, step] = -tau * T * (V' v), using only rows step:m of V
# since v is zero above its diagonal). Returns tau_i.
function _commit_reflector!(
    V::Matrix{Float64},
    Twy::Matrix{Float64},
    R11::Matrix{Float64},
    tau::Vector{Float64},
    tcol::Vector{Float64},
    xcol::Vector{Float64},
    nsel::Int,
    m::Int,
)
    step = nsel + 1
    @inbounds for r in 1:nsel
        R11[r, step] = xcol[r]
    end
    tau_i, beta = _householder!(view(xcol, step:m))
    R11[step, step] = beta
    V[step, step] = 1.0
    @inbounds for r in (step + 1):m
        V[r, step] = xcol[r]
    end
    tau[step] = tau_i
    if nsel > 0
        tv = view(tcol, 1:nsel)
        mul!(tv, view(V, step:m, 1:nsel)', view(V, step:m, step))
        tv .*= -tau_i
        lmul!(UpperTriangular(view(Twy, 1:nsel, 1:nsel)), tv)
        @inbounds for r in 1:nsel
            Twy[r, step] = tv[r]
        end
    end
    Twy[step, step] = tau_i
    return tau_i
end

# --- Fenwick tree for norm-weighted sampling without replacement ------------
# Sample an index with probability proportional to its current weight in
# O(log n); remove/restore a weight in O(log n). A step thus costs
# O(tested * log n) instead of the O(n log n) sort a key-based scheme needs.
# Port of the MEX's Fenwick (post-fix 20f562f): the tree's node sums and the
# separately accumulated cur_total drift apart under floating point when
# weights spanning many orders of magnitude are removed and restored, so
# find() alone cannot be trusted to return an in-pool index -- every draw is
# validated against the exact `inpool` membership array by the caller.
mutable struct Fenwick
    tree::Vector{Float64}
    n::Int
    log2n::Int
    cur_total::Float64
end

function Fenwick(w::AbstractVector{Float64})
    n = length(w)
    tree = zeros(Float64, n)
    log2n = 0
    while (1 << (log2n + 1)) <= n
        log2n += 1
    end
    # O(n) bottom-up build: seed each node with its own leaf weight and fold it
    # into its Fenwick parent in a single increasing pass (children j < i have
    # already contributed tree[j] into tree[i] by the time i is processed).
    total = 0.0
    @inbounds for i in 1:n
        tree[i] += w[i]
        total += w[i]
        parent = i + (i & -i)
        if parent <= n
            tree[parent] += tree[i]
        end
    end
    return Fenwick(tree, n, log2n, total)
end

function fenwick_add!(b::Fenwick, i::Int, d::Float64)
    b.cur_total += d
    j = i
    @inbounds while j <= b.n
        b.tree[j] += d
        j += j & -j
    end
    return nothing
end

# Smallest 1-based index whose cumulative weight exceeds x (0 <= x < total).
# May return n+1 when floating-point drift pushes x past the tree-reachable
# prefix total -- callers must validate the result against `inpool`.
function fenwick_find(b::Fenwick, x::Float64)
    pos = 0
    @inbounds for pw in b.log2n:-1:0
        nx = pos + (1 << pw)
        if nx <= b.n && b.tree[nx] <= x
            pos = nx
            x -= b.tree[nx]
        end
    end
    return pos + 1
end

# Draw one column by current weight, remove it from the pool, and return its
# index. Degenerate or drifted draws (all remaining weights ~0, x beyond the
# tree total, or a removed column hit at an exact prefix boundary) fall back
# to the first column still in the pool; one always exists because draws per
# block never exceed the number of in-pool columns.
function _draw_weighted_one!(
    bit::Fenwick,
    g::Vector{Float64},
    inpool::Vector{Bool},
    rng::AbstractRNG,
)
    n = bit.n
    id = n + 1
    total = bit.cur_total
    if total > 0.0
        x = rand(rng) * total
        if x >= total
            x = prevfloat(total)
        end
        id = fenwick_find(bit, x)
    end
    if id > n || !inpool[id]
        idx = findfirst(inpool)
        idx === nothing &&
            throw(ArgumentError("internal error: weighted sampling pool is empty"))
        id = idx
    end
    fenwick_add!(bit, id, -g[id])
    inpool[id] = false
    return id
end
