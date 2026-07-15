# Julia Tests And Publication Benchmarks

## Test Suite

Test deps (`Test`, `Random`, `DelimitedFiles`) are in the package's
`[targets].test`, so run the suite via `Pkg.test` (a direct `runtests.jl` call
would not see them):

```bash
julia --project=julia -e 'using Pkg; Pkg.test()'
```

## Publication Benchmark Workflow

The benchmark pipeline uses the separate `julia/benchmark` environment (it
carries BenchmarkTools, AppleAccelerate, etc.; the package itself depends only
on `LinearAlgebra`). Instantiate it once with
`julia --project=julia/benchmark -e 'using Pkg; Pkg.instantiate()'`.

Runner:

```bash
julia --project=julia/benchmark julia/benchmark/run_publication_benchmarks.jl
```

Smoke runner:

```bash
julia --project=julia/benchmark julia/benchmark/run_publication_smoke_benchmark.jl
```

Plot generation:

```bash
julia --project=julia/benchmark julia/benchmark/plot_publication_results.jl
```

Plot generation defaults to PNG. Set `BS_PUB_FIG_FORMATS=png,pdf,eps` to
also emit publication vector formats. Derived comparison artifacts report
relative time (`BSQR / baseline`); `1.0` is parity.

Perf gate:

```bash
julia --project=julia/benchmark julia/benchmark/check_publication_perf_gate.jl
```

## Environment Knobs

Publication knobs:

- `BS_PUB_OUTDIR` (default: `julia/benchmark/results/publication`)
- `BS_PUB_THREADS` (default: `1,4`)
- `BS_PUB_SEEDS` (default: `20260310,20260311`)
- `BS_PUB_FAMILIES` (default: `gaussian,ill_conditioned,orthonormal_rows`; constructions
  are defined in `docs/PUBLICATION_FIGURES_PLAN.md` and implemented in
  `julia/benchmark/matrix_generators.jl`)
- `BS_PUB_SQUARE_MS` (default: `64,128,256,384,512`)
- `BS_PUB_SHORT_MS` (default: `32,64,128,256,512`)
- `BS_PUB_SHORT_ASPECTS` (default: `2.0,4.0,8.0,10.0`)
- `BS_PUB_WARMUP` (default: `1`)
- `BS_PUB_SAMPLES` (default: `30`)
- `BS_NORM_RECOMP_TOL` (default: `sqrt(eps(Float64))`)
- `BS_PANEL_NB` (default: `8`; the panel/blocked kernel width — `0` or `1`
  selects the unblocked reference kernel; see `docs/P3_BLOCKED_BSQR.md`)
- `BS_PANEL_MIN_KN` (default: `24576`; minimum `k*n` for the panel kernel —
  below the crossover the unblocked kernel is faster and runs instead)
- `BS_PUB_FIG_FORMATS` (default: `png`; supported: `png,pdf,eps`)

Quality/CI knobs:

- `BS_PUB_CI_WARN_FRAC` (default: `0.5`)
- `BS_PUB_CI_FAIL_FRAC` (default: `10.0`)
- `BS_PUB_CI_MIN_TMED` (default: `1.0e-4`)
- `BS_PUB_CI_ENFORCE` (default: `0`)
- `BS_PUB_RESID_FACTOR` (default: `2500`)
- `BS_PUB_ORTH_FACTOR` (default: `2500`)

Environment/backend knobs (read by `bench_common.jl`, included by the runner):

- `BS_USE_ACCELERATE` (default: `1`; macOS only — load Apple Accelerate as the
  BLAS/LAPACK backend)
- `BS_REQUIRE_ACCELERATE` (default: `0`; macOS only — error if Accelerate did
  not activate instead of silently continuing on OpenBLAS)
- `BS_BLAS_THREADS` (default: unset; overrides the BLAS thread count for the
  whole run — the per-case `BS_PUB_THREADS` sweep still applies when unset)

Kernel fastpath knobs (read by the kernel itself, so they affect any
benchmarked factorization; defaults in `julia/src/kernel.jl`):

- `BS_SHORT_WIDE_FASTPATH` (default: `1`; `0` disables the short-wide fastpath)
- `BS_SHORT_WIDE_FASTPATH_ASPECT`, `BS_SHORT_WIDE_FASTPATH_MMAX`,
  `BS_SHORT_WIDE_FASTPATH_NMIN` (crossover tuning; see the kernel comments)

## Output Artifacts

- `julia/benchmark/results/publication/publication_timings.csv`
- `julia/benchmark/results/publication/publication_summary.md`
- `julia/benchmark/results/publication/metadata.txt`
- `julia/benchmark/results/publication/plots/plain/*`
- `julia/benchmark/results/publication/plots/rinv/*`
- `julia/benchmark/results/publication/tables/plain/*`
- `julia/benchmark/results/publication/tables/rinv/*`
