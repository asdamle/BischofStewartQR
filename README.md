# BSPivotQR

Julia implementation of Bischof-Stewart column-pivoted QR with an optional early-stop parameter `k`.

Detailed validation/performance documentation:

- `/Users/anildamle/Dropbox/Trystuff/BSpivoting/docs/TESTS_AND_BENCHMARKS.md`

## Requirements

- macOS 13.4+ (for Apple Accelerate integration)
- Julia 1.10+
- `AppleAccelerate.jl` (for Accelerate backend switching)

## Install dependencies

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## API

- `bsqr(A; k=min(size(A)...), check=true, track_inverse_frob=false, return_rinv_r12=false, rank_stop=true, norm_recomp_tol=sqrt(eps(Float64)), blas_threads=nothing)`
- `bsqr!(A::StridedMatrix{Float64}; k=min(size(A)...), check=true, track_inverse_frob=false, return_rinv_r12=false, rank_stop=true, norm_recomp_tol=sqrt(eps(Float64)), blas_threads=nothing, workspace=nothing)`
- `bsqr!(A, tau, jpvt, workspace; k=min(size(A)...), check=true, reset_pivots=true, frob_inv_trace=nothing, rank_stop=true, norm_recomp_tol=sqrt(eps(Float64)), blas_threads=nothing) -> ksteps`
- `R(F::BSQRPivoted)`
- `perm(F::BSQRPivoted)`
- `rinv_r12(F::BSQRPivoted)`
- `reconstruct(F, Aorig)`

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Apple Accelerate setup check

```bash
julia --project=. benchmark/setup_accelerate.jl
```

## Benchmarks

```bash
BS_USE_ACCELERATE=1 BS_REQUIRE_ACCELERATE=0 julia --project=. benchmark/bench_cpqr.jl
```

Set BLAS threads explicitly (recommended for reproducible timing):

```bash
BS_BLAS_THREADS=8 BS_USE_ACCELERATE=1 BS_REQUIRE_ACCELERATE=0 julia --project=. benchmark/bench_cpqr.jl
```

Benchmark fairness note: benchmark scripts compare `bsqr!` and `qr(..., ColumnNorm())` on equal footing using fresh `copy(A)` per trial, disable bsqr-only extras (`workspace` reuse, `track_inverse_frob`, `return_rinv_r12`), and force full-step runs (`k=min(m,n)`, `rank_stop=false`).
Use `BS_NORM_RECOMP_TOL` to tune the partial-norm recomputation threshold when exploring robustness/performance tradeoffs.

Quick smoke run:

```bash
BS_QUICK=1 BS_USE_ACCELERATE=1 BS_REQUIRE_ACCELERATE=0 julia --project=. benchmark/bench_cpqr.jl
```

Systematic size/aspect sweep benchmark:

```bash
BS_SWEEP_QUICK=1 julia --project=. benchmark/bench_cpqr_sweep.jl
```

Regime scaling benchmark (fixed `n` with increasing `m`, and fixed `m` with increasing `n`):

```bash
BS_REGIME_QUICK=1 julia --project=. benchmark/bench_cpqr_regimes.jl
```

Note: for the orthonormal-row family, regime benchmark only evaluates fixed-`m`/increasing-`n` (short-wide growth).

Kernel breakdown profile (pinpoint bsqr runtime origin by phase):

```bash
BS_PROFILE_QUICK=1 julia --project=. benchmark/profile_bsqr_breakdown.jl
```

## Validation workflow

Tier 1 (quick gate):

```bash
julia --project=. benchmark/run_validation_tier1.jl
```

Tier 2 (full gate):

```bash
julia --project=. benchmark/run_validation_tier2.jl
```

Baseline snapshot (for candidate-vs-baseline comparisons):

```bash
julia --project=. benchmark/freeze_baseline.jl
```

Compare current outputs against a baseline snapshot:

```bash
julia --project=. benchmark/compare_results.jl benchmark/results/baseline_YYYYMMDD_HHMMSS
```

Outputs are written to:

- `/Users/anildamle/Dropbox/Trystuff/BSpivoting/benchmark/results/timings.csv`
- `/Users/anildamle/Dropbox/Trystuff/BSpivoting/benchmark/results/summary.md`
- `/Users/anildamle/Dropbox/Trystuff/BSpivoting/benchmark/results/sweep_timings.csv`
- `/Users/anildamle/Dropbox/Trystuff/BSpivoting/benchmark/results/sweep_summary.md`
- `/Users/anildamle/Dropbox/Trystuff/BSpivoting/benchmark/results/regime_timings.csv`
- `/Users/anildamle/Dropbox/Trystuff/BSpivoting/benchmark/results/regime_summary.md`
- `/Users/anildamle/Dropbox/Trystuff/BSpivoting/benchmark/results/profile_breakdown.csv`
- `/Users/anildamle/Dropbox/Trystuff/BSpivoting/benchmark/results/profile_breakdown_summary.md`
- `/Users/anildamle/Dropbox/Trystuff/BSpivoting/benchmark/results/validation_report.md`
- `/Users/anildamle/Dropbox/Trystuff/BSpivoting/benchmark/results/guardrail_failures.csv`

## Plotting benchmark results

Install plotting dependency:

```bash
python3 -m pip install matplotlib
```

Generate plots from benchmark CSV:

```bash
julia --project=. benchmark/plot_results.jl
```

Optional input/output override:

```bash
julia --project=. benchmark/plot_results.jl benchmark/results/timings.csv benchmark/results/plots
```

Generate regime scaling plots:

```bash
julia --project=. benchmark/plot_regime_results.jl
```
