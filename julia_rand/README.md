# julia_rand — randomized Bischof–Stewart column selection (Julia port)

`BSRandPivotQR`: the Julia implementation of the randomized BSQR variant in
`matlab_rand/`. It selects `k` columns of a `k×n` matrix with the same
theoretical guarantee on `||R11^{-1}||_F` as the deterministic kernel, but
**without** maintaining `R11^{-1}R12` / column norms for every column at every
step: it tracks only the running squared inverse Frobenius norm and samples
candidate columns in blocks, keeping those that hold the running value under
the per-step bound. By default it runs **in-block** (`batched`): each sampled
block is brought to the current frame once (one compact-WY BLAS-3 apply), then
BSQR is run within it to take as many columns as the bound allows before
resampling (`O(k^3)` overall vs the single-select `O(k^4)`). See
`docs/RANDOMIZED_BSQR_PLAN.md` and §5 of the manuscript
([arXiv:2607.13532](https://arxiv.org/abs/2607.13532), Thm. 5.1).

The kernel is behaviourally in lockstep with `matlab_rand` (same acceptance
rule, threshold formulas, fallback semantics, and stats attribution) and is a
port of the MEX backend's performance design: compact-WY block applies,
incremental in-block norm downdates with the Businger–Golub recompute
safeguard, a Fenwick tree for `O(log n)` norm-weighted draws (with exact
pool-membership validation of every draw), and a partial Fisher–Yates pool for
uniform draws. It uses its own RNG, so **pivot sequences need not match either
MATLAB backend** — only the guarantees are the cross-language contract; there
are no parity fixtures for the randomized variant.

Like `matlab_rand/`, this package is **decoupled from the deterministic
implementations**: it depends only on stdlib (`LinearAlgebra`, `Random`) and
keeps its own copies of kernel helpers. `BSPivotQR` appears only as a
test-target dependency (a conditioning-ballpark comparison) and in the
separate `benchmark/` environment.

## Quick start

```julia
julia --project=julia_rand
```

```julia
using BSRandPivotQR, LinearAlgebra

M = Matrix(Matrix(qr(randn(2000, 32)).Q)')   # 32×2000, orthonormal rows (GKS setting)
F = bsqr_rand(M; seed = 1)                   # BSQRRandPivoted
selected(F)                                  # the 32 selected column indices
perm(F)                                      # full permutation; selected first
F.R11                                        # 32×32 upper-triangular factor
Q(F)                                         # economy Q, formed lazily (orgqr) on request
F.stats.frob_inv / F.stats.osinsky_bound     # realized quality vs the guarantee
```

`Q` is formed only when requested — callers that want only the subset never
pay for it. `R12` is opt-in (an extra `O(nk^2)` pass, off by default):

```julia
F = bsqr_rand(M; return_r12 = true); F.R12   # k×(n-k) coupling block
```

Options (`k`, `batched`, `block_size`, `threshold_mode`, `slack`,
`norm_recomp_tol`, `sampling`, `pick`, `seed`/`rng`, `return_r12`,
`check_finite`) mirror `matlab_rand/bsqr_rand.m` one-for-one — see `?bsqr_rand`
and the option table in `matlab_rand/README.md`. The `stats` field carries the
same instrumentation (`f2`, `Fhat`, `crit`, `threshold`, `samples_tested`,
`rounds`, `fallback`, plus the scalar totals), with the same batched
attribution rule. The differences from MATLAB: there is a single kernel (no
`backend` option — Julia needs no MEX split), errors are `ArgumentError` /
`RankDeficientError` instead of MATLAB identifiers, and enums are `Symbol`s.

## Tests

```bash
julia --project=julia_rand -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=julia_rand -e 'using Pkg; Pkg.test()'
```

`test_bsqr_rand.jl` checks the factorization is exact across the option grid
and kernel corners (including the Fenwick-drift pivot-validity regression on
concentrated-norm input); `test_bsqr_rand_bounds.jl` checks the headline
guarantee (`||R11^{-1}||_F` under `sqrt(k(n-k+1))`, per-step) and the
ballpark comparison against deterministic `BSPivotQR.bsqr`. The opt-in
statistical sweep (seeds × paths × sampling × threshold modes × families):

```bash
julia --project=julia_rand julia_rand/test/stress_bsqr_rand_bounds.jl 100
```

## Benchmarks

Implementation-validation only (not publication material — no committed
snapshot, no plots, no external-method comparisons; those live in
`matlab_rand/benchmark/`):

```bash
julia --project=julia_rand/benchmark -e 'using Pkg; Pkg.instantiate()'   # once
julia --project=julia_rand/benchmark julia_rand/benchmark/run_rand_benchmarks.jl
```

Compares the randomized selection (`p`, lazy `Q`, `R11` — no `R12`, plus a
`+R12` row) against (a) deterministic `BSPivotQR.bsqr` and (b) built-in
`qr(M, ColumnNorm())` (LAPACK dgeqp3), all sides materializing `Q`, `R`, `p`,
on the `gaussian` and `needle` families. Reports speedups, the conditioning
ratio `||R11^{-1}||_F(rand) / (det)`, and candidates tested per selection;
writes `benchmark/results/rand_timings.csv` (git-ignored) and a
`metadata.txt`. Configure via `BS_RAND_SIZES`, `BS_RAND_FAMILIES`,
`BS_RAND_SEED`, `BS_RAND_OUTDIR`. On macOS the run uses Apple Accelerate (via
`julia/benchmark/setup_accelerate.jl`).
