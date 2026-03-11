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
- `BS_HOUSEHOLDER_LAPACK_LARF`
  - reflector-apply implementation toggle.
  - `1` (default): `LAPACK.larf!` path.
  - `0`: manual BLAS (`gemv`/`ger`) fallback.

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

## Output Artifacts

Publication benchmark runner writes:

- `benchmark/results/publication/publication_timings.csv`
- `benchmark/results/publication/publication_summary.md`
- `benchmark/results/publication/metadata.txt`

Publication plotting writes:

- `benchmark/results/publication/plots/figure1_square_runtime.{png,pdf}`
- `benchmark/results/publication/plots/figure2_shortwide_runtime.{png,pdf}`
- `benchmark/results/publication/plots/figure3_shortwide_speedup_heatmap.{png,pdf}`
- `benchmark/results/publication/plots/figure4_quality.{png,pdf}`
- `benchmark/results/publication/plots/figure5_aggregate_speedup.{png,pdf}`
- `benchmark/results/publication/tables/table_square_speedup.csv`
- `benchmark/results/publication/tables/table_shortwide_speedup.csv`
- `benchmark/results/publication/tables/table_quality.csv`
