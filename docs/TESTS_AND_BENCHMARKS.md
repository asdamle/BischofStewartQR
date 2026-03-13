# Tests And Publication Benchmarks

This repository now maintains a single benchmark workflow focused on publication artifacts.
Legacy benchmark/validation pipelines and their generated outputs were removed.

## Test Suite

Run all tests:

```bash
julia --project=. test/runtests.jl
```

Or via package test entrypoint:

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Publication Benchmark Workflow

Primary benchmark runner:

```bash
julia --project=. benchmark/bench_cpqr_publication.jl
```

Plot + table generation:

```bash
julia --project=. benchmark/plot_publication_results.jl
```

Internal performance headroom gate (stage timings + variant comparison):

```bash
julia --project=. benchmark/perf_headroom_gate.jl
```

Publication timing compares:
- Plain pair: `bsqr_full` vs `dgeqp3`.
- `R11^{-1}R12` pair: `bsqr_rinv` vs `dgeqp3_trsm` (DGEQP3 + triangular solve).

### Important Environment Knobs

- `BS_PUB_OUTDIR`
  - output directory for publication artifacts.
  - default: `benchmark/results/publication`.
- `BS_PUB_THREADS`
  - comma-separated BLAS thread settings.
  - default: `1,4`.
- `BS_PUB_SEEDS`
  - comma-separated RNG seeds.
  - default: `20260310,20260311`.
- `BS_PUB_FAMILIES`
  - comma-separated family subset.
  - allowed: `gaussian`, `ill_conditioned`, `orthonormal_rows`.
  - default: all three.
- `BS_PUB_SQUARE_MS`
  - square case sizes.
  - default: `64,128,256,384,512`.
- `BS_PUB_SHORT_MS`
  - short-wide row sizes.
  - default: `32,64,128,256,512`.
- `BS_PUB_SHORT_ASPECTS`
  - short-wide aspect ratios `n/m`.
  - default: `2.0,4.0,8.0,10.0`.
- `BS_PUB_WARMUP`
  - warmup evaluations per timing point.
  - default: `1`.
- `BS_PUB_SAMPLES`
  - timing samples per method/case.
  - default: `30`.
- `BS_NORM_RECOMP_TOL`
  - partial-norm recomputation tolerance.
  - default: `sqrt(eps(Float64))`.
- `BS_PUB_CI_WARN_FRAC`
  - warning threshold for CI spread `(tci_high - tci_low) / tmed`.
  - default: `0.5`.
- `BS_PUB_CI_FAIL_FRAC`
  - hard-fail threshold for CI spread.
  - default: `10.0`.
- `BS_PUB_CI_MIN_TMED`
  - CI warn/fail checks apply only when `tmed >= BS_PUB_CI_MIN_TMED`.
  - default: `1.0e-4`.
- `BS_PUB_CI_ENFORCE`
  - set to `1` to hard-fail when CI spread exceeds `BS_PUB_CI_FAIL_FRAC`.
  - default: `0` (warn-only).
- `BS_PUB_RESID_FACTOR`
  - residual tolerance factor in `factor * eps(Float64) * max(m,n)`.
  - default: `2500`.
- `BS_PUB_ORTH_FACTOR`
  - orthogonality tolerance factor in `factor * eps(Float64) * m`.
  - default: `2500`.

Kernel fastpath tuning knobs (used by `bsqr` and perf gate):

- `BS_SHORT_WIDE_FASTPATH` (`1`/`0`, default `1`)
- `BS_SHORT_WIDE_FASTPATH_ASPECT` (default `4`)
- `BS_SHORT_WIDE_FASTPATH_MMAX` (default `256`)
- `BS_SHORT_WIDE_FASTPATH_NMIN` (default `0`)

Headroom gate knobs (`benchmark/perf_headroom_gate.jl`):

- `BS_HEADROOM_BLAS_THREADS` (default `1`)
- `BS_HEADROOM_WARMUP` (default `1`)
- `BS_HEADROOM_SAMPLES` (default `8`)
- `BS_HEADROOM_SEED` (default `20260312`)
- `BS_HEADROOM_VARIANTS` (default `baseline,fastpath_off,aspect6`)
- `BS_HEADROOM_FAMILIES` (default `gaussian,ill_conditioned,orthonormal_rows`)
- `BS_HEADROOM_GATE` (default `1.15`)
- `BS_HEADROOM_ENFORCE` (`1` to hard-fail when no variant reaches the gate; default `0`)

### Example Commands

Standard publication run:

```bash
julia --project=. benchmark/bench_cpqr_publication.jl
julia --project=. benchmark/plot_publication_results.jl
```

Explicit knobs:

```bash
BS_PUB_THREADS=1,4 BS_PUB_SEEDS=20260310,20260311 BS_PUB_WARMUP=1 BS_PUB_SAMPLES=30 julia --project=. benchmark/bench_cpqr_publication.jl
julia --project=. benchmark/plot_publication_results.jl benchmark/results/publication/publication_timings.csv benchmark/results/publication/plots benchmark/results/publication/tables
```

Smoke run:

```bash
BS_PUB_OUTDIR=/tmp/bs_pub_smoke BS_PUB_FAMILIES=gaussian BS_PUB_THREADS=1 BS_PUB_SEEDS=20260310 BS_PUB_SQUARE_MS=64 BS_PUB_SHORT_MS=32 BS_PUB_SHORT_ASPECTS=2 BS_PUB_WARMUP=0 BS_PUB_SAMPLES=2 julia --project=. benchmark/bench_cpqr_publication.jl
julia --project=. benchmark/plot_publication_results.jl /tmp/bs_pub_smoke/publication_timings.csv /tmp/bs_pub_smoke/plots /tmp/bs_pub_smoke/tables
```

## Publication Checklist

- Run tests: `julia --project=. test/runtests.jl`.
- Run publication benchmark with fixed seeds/threads.
- Confirm `publication_timings.csv` row count equals:
  - `|threads| * |seeds| * |families| * |cases| * 4 methods`.
- Confirm no CI/quality gate failures in benchmark runner output.
- Check `metadata.txt` contains:
  - `schema_version`, `git_sha`, `git_branch`, `git_dirty`, `blas_config`, `expected_rows`, `observed_rows`.
- Run plot pipeline and ensure CSV validation succeeds before figure generation.
- Verify expected artifact tree exists:
  - plots and tables for both `plain` and `rinv`.
- (Optional) Run `benchmark/perf_headroom_gate.jl` to decide whether performance tuning changes are substantial enough to keep.

## Output Artifacts

Publication benchmark runner writes:

- `benchmark/results/publication/publication_timings.csv`
- `benchmark/results/publication/publication_summary.md`
- `benchmark/results/publication/metadata.txt`

Publication plotting writes:

- `benchmark/results/publication/plots/plain/figure1_square_runtime.{png,pdf}`
- `benchmark/results/publication/plots/plain/figure2_shortwide_runtime.{png,pdf}`
- `benchmark/results/publication/plots/plain/figure3_shortwide_speedup_heatmap.{png,pdf}`
- `benchmark/results/publication/plots/plain/figure4_quality.{png,pdf}`
- `benchmark/results/publication/plots/plain/figure5_aggregate_speedup.{png,pdf}`
- `benchmark/results/publication/plots/rinv/figure1_square_runtime.{png,pdf}`
- `benchmark/results/publication/plots/rinv/figure2_shortwide_runtime.{png,pdf}`
- `benchmark/results/publication/plots/rinv/figure3_shortwide_speedup_heatmap.{png,pdf}`
- `benchmark/results/publication/plots/rinv/figure4_quality.{png,pdf}`
- `benchmark/results/publication/plots/rinv/figure5_aggregate_speedup.{png,pdf}`
- `benchmark/results/publication/tables/plain/table_square_speedup.csv`
- `benchmark/results/publication/tables/plain/table_shortwide_speedup.csv`
- `benchmark/results/publication/tables/plain/table_quality.csv`
- `benchmark/results/publication/tables/rinv/table_square_speedup.csv`
- `benchmark/results/publication/tables/rinv/table_shortwide_speedup.csv`
- `benchmark/results/publication/tables/rinv/table_quality.csv`

Headroom gate writes:

- no files; prints markdown tables + gate decision to stdout.
