# Julia BSQR

Canonical Julia implementation of Bischof-Stewart column-pivoted QR
(`BSPivotQR`). The package depends only on `LinearAlgebra`; the heavier
benchmark/dev dependencies live in the separate `julia/benchmark/` environment.

## Quick start

```julia
julia> include("startup.jl")        # from the repo root: activate + instantiate, then `using BSPivotQR`
julia> F = bsqr(A)                  # deterministic BS pivoted QR -> BSQRPivoted
julia> Q(F), R(F), perm(F)          # economy Q, R factor, column permutation
```

Equivalently: `julia --project=julia` then `using BSPivotQR`.

## Results and accessors

`bsqr(A; ...)` / `bsqr!(A; ...)` return a `BSQRPivoted` holding the packed
factorization (`k` steps, default `k = min(m,n)`); the exported accessors
materialize the pieces:

- `R(F)`: `k×n` upper-trapezoidal factor of the *permuted* columns;
  `R(F)[1:k,1:k]` is `R11`, `R(F)[:,k+1:n]` is `R12`. For the default full
  `k`, `A[:, perm(F)] ≈ Q(F)*R(F)`; with early stop (`k < min(m,n)`) that
  holds only for the selected block, `A[:, perm(F)[1:k]] ≈ Q(F)*R(F)[:,1:k]`,
  and `R12 = Q(F)'*A[:, perm(F)[k+1:n]]` is the unselected columns'
  projection onto `span(Q)`.
- `Q(F)`: `m×k` economy factor with orthonormal columns, materialized on each
  call from the packed reflectors (LAPACK `orgqr`, one `O(mk²)` pass).
- `perm(F)`: the length-`n` column permutation; `perm(F)[1:k]` are the
  selected columns in pivot order.
- `rinv_r12(F)`: `k×(n−k)` matrix `R11⁻¹R12`, captured from the kernel
  workspace with no extra solve — `nothing` unless the factorization was run
  with `return_rinv_r12 = true`.
- `reconstruct(F, A)`: rebuilds the original matrix from the factorization
  (validation/tests).

Kernel keyword arguments (`k`, `check`, `track_inverse_frob`,
`return_rinv_r12`, `rank_stop`, `norm_recomp_tol`, `blas_threads`, and the
preallocated `bsqr!(A, tau, jpvt, workspace; ...)` form) are documented in
the docstrings — `?bsqr`, `?bsqr!`, `?BSQRPivoted`.

## Setup

```bash
julia --project=julia           -e 'using Pkg; Pkg.instantiate()'   # core package env
julia --project=julia/benchmark -e 'using Pkg; Pkg.instantiate()'   # benchmark env
```

## Tests

Test dependencies (`Test`, `Random`, `DelimitedFiles`) live in the package's
`[targets].test`, so run them through `Pkg.test`:

```bash
julia --project=julia -e 'using Pkg; Pkg.test()'
```

## Benchmarks

The benchmark pipeline uses the `julia/benchmark` environment.

Publication benchmark:

```bash
julia --project=julia/benchmark julia/benchmark/run_publication_benchmarks.jl
```

Smoke benchmark:

```bash
julia --project=julia/benchmark julia/benchmark/run_publication_smoke_benchmark.jl
```

Plot from CSV:

```bash
julia --project=julia/benchmark julia/benchmark/plot_publication_results.jl
```

Publication plots default to PNG. Set `BS_PUB_FIG_FORMATS=png,pdf,eps`
to request vector outputs as well. Derived performance plots and tables report
relative time (`BSQR / baseline`), where `1.0` is parity.

Performance gate:

```bash
julia --project=julia/benchmark julia/benchmark/check_publication_perf_gate.jl
```

Cross-language timing comparison helper:

```bash
python3 julia/benchmark/compare_matlab_julia_timings.py \
  --matlab-csv matlab/benchmark/results/publication/publication_timings.csv \
  --julia-csv julia/benchmark/results/publication/publication_timings.csv \
  --method bsqr_full --julia-threads 1
```

## Output Defaults

Publication artifacts default to:

- `julia/benchmark/results/publication/publication_timings.csv`
- `julia/benchmark/results/publication/publication_summary.md`
- `julia/benchmark/results/publication/metadata.txt`
- `julia/benchmark/results/publication/plots/...`
- `julia/benchmark/results/publication/tables/...`

Detailed benchmark/test notes: `julia/docs/TESTS_AND_BENCHMARKS.md`.
