# Julia BSQR

Canonical Julia implementation of Bischof-Stewart column-pivoted QR
(`BSPivotQR`). The package depends only on `LinearAlgebra`; the heavier
benchmark/dev dependencies live in the separate `julia/benchmark/` environment.

## Quick start

```julia
julia> include("startup.jl")        # from the repo root: activate + instantiate, then `using BSPivotQR`
julia> F = bsqr(A)                  # deterministic BS pivoted QR -> BSQRPivoted
julia> R(F), perm(F)               # R factor, column permutation
```

Equivalently: `julia --project=julia` then `using BSPivotQR`.

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
