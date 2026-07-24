"""
    BSQRRandStats

Instrumentation from one [`bsqr_rand`](@ref) run. Per-step vectors have length
`k`, indexed by selection step:
- `f2` — running `‖R11⁻¹‖_F²` after each step
- `crit` — criterion increment `(1 + ‖w‖²)/ρ²` of the accepted pivot
- `threshold` — acceptance threshold used at each step
- `Fhat` — deterministic worst-case bound `F̂ᵢ = i(n-i+1)/(k-i+1)`
- `samples_tested` — candidate columns evaluated per selection, INCLUDING
  blocks that yielded no selection (batched: all work since the previous
  selection is attributed to a block's first pick, 0 for its later in-block
  picks)
- `rounds` — sampled blocks / reflector applies per selection (same
  attribution as `samples_tested`)
- `fallback` — true where the exhaustive scan force-accepted the global
  minimizer above the threshold (rounding-tie safety net)

Scalars: `frob_inv = sqrt(f2[end])`, `osinsky_bound = sqrt(k(n-k+1))` (the
Frobenius guarantee target), `total_tested = sum(samples_tested)`, and
`blocks_sampled = sum(rounds)` (true totals of candidate evaluations and
reflector applies).
"""
struct BSQRRandStats
    f2::Vector{Float64}
    crit::Vector{Float64}
    threshold::Vector{Float64}
    Fhat::Vector{Float64}
    samples_tested::Vector{Int}
    rounds::Vector{Int}
    fallback::Vector{Bool}
    frob_inv::Float64
    osinsky_bound::Float64
    total_tested::Int
    blocks_sampled::Int
end

"""
    BSQRRandPivoted <: Factorization{Float64}

Result of [`bsqr_rand`](@ref): a randomized Bischof-Stewart column selection
with its partial QR of the selected columns. Access the results with
[`perm`](@ref), [`selected`](@ref), [`Q`](@ref), [`R11`](@ref), [`r12`](@ref);
the identities are `A[:, selected(F)] ≈ Q(F) * F.R11` (equivalently
`Q(F)' * A[:, selected(F)] ≈ F.R11`).

Fields:
- `V` — the `m×k` Householder reflector store in unit-diagonal packed form
  (together with `tau` it defines Q implicitly; materialize with `Q(F)`).
- `tau` — the `k` Householder coefficients.
- `R11` — the `k×k` upper-triangular factor of the selected columns.
- `p` — the length-`n` permutation: `p[1:k]` are the selected columns in
  pivot order, `p[k+1:n]` the unselected ones.
- `k` — number of selected columns.
- `stats` — a [`BSQRRandStats`](@ref).
- `R12` — the `k×(n-k)` coupling block `Q' * A[:, p[k+1:n]]`; `nothing`
  unless `return_r12 = true` (it costs an extra `O(nk²)` pass).
"""
struct BSQRRandPivoted <: Factorization{Float64}
    V::Matrix{Float64}
    tau::Vector{Float64}
    R11::Matrix{Float64}
    p::Vector{Int}
    k::Int
    stats::BSQRRandStats
    R12::Union{Nothing,Matrix{Float64}}
end

Base.size(F::BSQRRandPivoted) = (size(F.V, 1), length(F.p))

"""
    bsqr_rand(A; kwargs...) -> BSQRRandPivoted

Randomized Bischof-Stewart column selection: select `k` columns of `A` with
the same `‖R11⁻¹‖_F ≤ sqrt(k(n-k+1))` guarantee as the deterministic kernel
(on orthonormal-row input, the GKS setting), but WITHOUT maintaining
`R11⁻¹R12` / column norms for every column at every step. The kernel tracks
only the running squared inverse Frobenius norm `f2 = ‖R11⁻¹‖_F²` and samples
candidate columns in blocks, keeping those that hold `f2` under the per-step
threshold. See `docs/RANDOMIZED_BSQR_PLAN.md` and §5 of the manuscript
(arXiv:2607.13532, Thm. 5.1).

This is the Julia port of `matlab_rand/bsqr_rand` (kernel behaviour in
lockstep with the MEX backend); it uses its own RNG, so pivot sequences need
not match either MATLAB backend — only the guarantees are the cross-language
contract.

# Keyword arguments
- `k::Integer = min(size(A)...)`: number of columns to select.
- `batched::Bool = true`: in-block BSQR — each sampled block is brought to the
  current frame once, then BSQR runs within it, taking as many columns as the
  bound allows before resampling (amortizes the per-block reflector apply:
  `O(k³)` vs the single-select `O(k⁴)`). `false` = one selection per block
  (tighter realized conditioning, more applies). Same bound either way.
- `block_size::Union{Nothing,Integer} = nothing`: candidates per sampled block
  (`nothing` = auto: `k` when batched, else `ceil(k/2)` clamped to `[16, 64]`).
  Larger improves realized quality at higher cost; the bound is unchanged.
- `threshold_mode::Symbol = :running_mean`: the per-step acceptance threshold,
  strict at every step (per-singular-value control) — or
  `:worstcase_allowance` (more permissive: spends accumulated slack, fewer
  samples, same final bound).
- `slack::Real = 1.0`: `≥ 1` multiplier loosening the threshold (values `> 1`
  no longer guarantee the Osinsky bound — experimental).
- `norm_recomp_tol::Real = sqrt(eps())`: running-norm recompute safeguard in
  `[0, 1]`, matching the deterministic kernel; used by the batched path's
  incremental in-block residual-norm downdates.
- `sampling::Symbol = :normweighted`: by starting squared column norms (robust
  across leverage profiles; adds `O(mn)`, drawn via a Fenwick tree in
  `O(log n)` per draw), or `:uniform`.
- `pick::Symbol = :best_in_block`: single-select only (`batched = false`): or
  `:first`. Batched always takes the in-block minimizer.
- `seed::Union{Nothing,Integer} = nothing`: convenience RNG seed (uses
  `Xoshiro(seed)`); overrides `rng`.
- `rng::AbstractRNG = Random.default_rng()`: the RNG to draw from.
- `return_r12::Bool = false`: compute the coupling block `R12` (an extra
  `O(nk²)` pass), stored in the result's `R12` field.
- `check_finite::Bool = false`: validate `A` is finite with a full `O(mn)`
  scan. Norm-weighted sampling (the default) still detects non-finite input
  for free from its precomputed weights; only `:uniform` needs `true`.

Numerical domain: the kernel tracks SQUARED column norms (`g`, `ρ²`) and the
squared inverse Frobenius norm `f2`, so column norms should stay roughly
within `[1e-150, 1e150]`; outside that range the acceptance rule degrades
silently. Throws [`RankDeficientError`](@ref) when `k` exceeds the numerical
rank.
"""
function bsqr_rand(
    A::AbstractMatrix{<:Real};
    k::Integer = min(size(A)...),
    batched::Bool = true,
    block_size::Union{Nothing,Integer} = nothing,
    threshold_mode::Symbol = :running_mean,
    slack::Real = 1.0,
    norm_recomp_tol::Real = _DEFAULT_NORM_RECOMP_TOL,
    sampling::Symbol = :normweighted,
    pick::Symbol = :best_in_block,
    seed::Union{Nothing,Integer} = nothing,
    rng::AbstractRNG = Random.default_rng(),
    return_r12::Bool = false,
    check_finite::Bool = false,
)
    m, n = size(A)
    kmax = min(m, n)
    (0 <= k <= kmax) || throw(ArgumentError("k must satisfy 0 <= k <= min(m, n)"))
    kw = Int(k)
    threshold_mode in (:running_mean, :worstcase_allowance) || throw(ArgumentError(
        "threshold_mode must be :running_mean or :worstcase_allowance"))
    sampling in (:normweighted, :uniform) ||
        throw(ArgumentError("sampling must be :normweighted or :uniform"))
    pick in (:best_in_block, :first) ||
        throw(ArgumentError("pick must be :best_in_block or :first"))
    (isfinite(slack) && slack >= 1.0) ||
        throw(ArgumentError("slack must be a finite scalar >= 1"))
    tol = Float64(norm_recomp_tol)
    (0.0 <= tol <= 1.0) ||
        throw(ArgumentError("norm_recomp_tol must satisfy 0 <= norm_recomp_tol <= 1"))
    if block_size === nothing
        # Auto block size (mirrors matlab_rand: batched favours block = k;
        # single-select uses rand_default_block, ceil(k/2) clamped to [16,64]).
        block = batched ? max(kw, 1) : clamp(cld(kw, 2), 16, 64)
    else
        block_size >= 1 || throw(ArgumentError("block_size must be a positive integer"))
        block = Int(block_size)
    end
    block = min(block, max(n, 1))
    if seed !== nothing
        rng = Xoshiro(seed)
    end

    Awork = A isa Matrix{Float64} ? A : Matrix{Float64}(A)
    if check_finite
        all(isfinite, Awork) || throw(ArgumentError("A contains non-finite values"))
    end

    if kw == 0
        stats = BSQRRandStats(
            Float64[], Float64[], Float64[], Float64[], Int[], Int[], Bool[],
            0.0, 0.0, 0, 0,
        )
        R12 = return_r12 ? zeros(0, n) : nothing
        return BSQRRandPivoted(
            zeros(m, 0), Float64[], zeros(0, 0), collect(1:n), 0, stats, R12)
    end

    st = _make_kernel(
        Awork, kw, block, threshold_mode, Float64(slack), tol,
        sampling === :normweighted, pick,
    )
    f2 = batched ? _run_batched!(st, rng) : _run_single!(st, rng)

    p = Vector{Int}(undef, n)
    copyto!(p, st.selected)
    w = kw
    @inbounds for j in 1:n
        if !st.taken[j]
            w += 1
            p[w] = j
        end
    end

    stats = BSQRRandStats(
        st.st_f2, st.st_crit, st.st_thr, st.st_Fhat,
        st.st_samples, st.st_rounds, st.st_fallback,
        sqrt(max(f2, 0.0)), sqrt(kw * (n - kw + 1.0)),
        sum(st.st_samples), sum(st.st_rounds),
    )

    R12 = nothing
    if return_r12
        R12 = _form_r12(st, p)
    end
    return BSQRRandPivoted(st.V, st.tau, triu!(st.R11), p, kw, stats, R12)
end

# R12 = Q' * A[:, p[k+1:n]] (top k rows), one compact-WY apply over the
# leftover columns.
function _form_r12(st::_RandKernel, p::Vector{Int})
    m, n, k = st.m, st.n, st.k
    n12 = n - k
    n12 > 0 || return zeros(k, 0)
    Xr = Matrix{Float64}(undef, m, n12)
    @inbounds for j in 1:n12
        copyto!(view(Xr, :, j), view(st.A, :, p[k + j]))
    end
    WBr = Matrix{Float64}(undef, k, n12)
    Vv = view(st.V, :, 1:k)
    mul!(WBr, Vv', Xr)
    lmul!(adjoint(UpperTriangular(view(st.Twy, 1:k, 1:k))), WBr)
    mul!(Xr, Vv, WBr, -1.0, 1.0)
    return Matrix(view(Xr, 1:k, :))
end

# --- accessors --------------------------------------------------------------

"""
    Q(F::BSQRRandPivoted) -> Matrix{Float64}

The economy orthogonal factor: `m×k` with orthonormal columns, satisfying
`A[:, selected(F)] ≈ Q(F) * F.R11`. `F` stores Q only implicitly (packed
reflectors in `V` + `tau`); each call materializes it with LAPACK `orgqr`,
one `O(mk²)` pass — callers that want only the subset never pay for it.
"""
function Q(F::BSQRRandPivoted)
    m = size(F.V, 1)
    k = F.k
    k == 0 && return Matrix{Float64}(undef, m, 0)
    Qm = Matrix{Float64}(undef, m, k)
    copyto!(Qm, F.V)
    return LAPACK.orgqr!(Qm, F.tau, k)
end

"""
    R11(F::BSQRRandPivoted) -> Matrix{Float64}

The `k×k` upper-triangular factor of the selected columns (a copy of `F.R11`).
"""
R11(F::BSQRRandPivoted) = copy(F.R11)

"""
    perm(F::BSQRRandPivoted) -> Vector{Int}

The length-`n` column permutation (a copy of `F.p`): `perm(F)[1:F.k]` are the
selected columns in pivot order, the rest the unselected columns.
"""
perm(F::BSQRRandPivoted) = copy(F.p)

"""
    selected(F::BSQRRandPivoted) -> Vector{Int}

The selected column subset `perm(F)[1:F.k]`, in pivot order.
"""
selected(F::BSQRRandPivoted) = F.p[1:F.k]

"""
    r12(F::BSQRRandPivoted) -> Union{Nothing, Matrix{Float64}}

The `k×(n-k)` coupling block `Q' * A[:, perm(F)[k+1:n]]` (a copy) — or
`nothing` unless the factorization was computed with `return_r12 = true`.
"""
r12(F::BSQRRandPivoted) = F.R12 === nothing ? nothing : copy(F.R12)
