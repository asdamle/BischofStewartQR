# BSPivotQR

Julia implementation of Bischof-Stewart column-pivoted QR with optional early stopping (`rank_stop=true`).

Detailed validation/performance documentation:

- `docs/TESTS_AND_BENCHMARKS.md`

## Requirements

- macOS 13.4+ (for Apple Accelerate integration)
- Julia 1.10+
- `AppleAccelerate.jl` (for Accelerate backend switching)

## Install dependencies

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

For validation pipelines and plotting scripts, also install matplotlib:

```bash
python3 -m pip install matplotlib
```

## API

- `bsqr(A; k=min(size(A)...), check=true, track_inverse_frob=false, return_rinv_r12=false, rank_stop=false, norm_recomp_tol=sqrt(eps(Float64)), blas_threads=nothing)`
- `bsqr!(A::StridedMatrix{Float64}; k=min(size(A)...), check=true, track_inverse_frob=false, return_rinv_r12=false, rank_stop=false, norm_recomp_tol=sqrt(eps(Float64)), blas_threads=nothing, workspace=nothing)`
- `bsqr!(A, tau, jpvt, workspace; k=min(size(A)...), check=true, reset_pivots=true, frob_inv_trace=nothing, rank_stop=false, norm_recomp_tol=sqrt(eps(Float64)), blas_threads=nothing) -> ksteps`
- `R(F::BSQRPivoted)`
- `perm(F::BSQRPivoted)`
- `rinv_r12(F::BSQRPivoted)`
- `reconstruct(F, Aorig)`

`rank_stop=false` is the default so `bsqr` performs full CPQR by default (QR drop-in behavior). Set `rank_stop=true` to enable early stopping.

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

## Apple Accelerate setup check

```bash
julia --project=. benchmark/setup_accelerate.jl
```

## Publication benchmark workflow

Primary benchmark script:

```bash
julia --project=. benchmark/bench_cpqr_publication.jl
```

Generate publication figures and caption-ready tables for both comparison sets:

- `plain`: `bsqr_full` vs `dgeqp3`
- `rinv`: `bsqr_rinv` vs `dgeqp3_trsm`

```bash
julia --project=. benchmark/plot_publication_results.jl
```

Internal performance headroom gate (stage breakdown + variant geomean speedups):

```bash
julia --project=. benchmark/perf_headroom_gate.jl
```

Useful knobs:

- `BS_PUB_OUTDIR` (default: `benchmark/results/publication`)
- `BS_PUB_THREADS` (default: `1,4`)
- `BS_PUB_SEEDS` (default: `20260310,20260311`)
- `BS_PUB_FAMILIES` (default: `gaussian,ill_conditioned,orthonormal_rows`)
- `BS_PUB_SQUARE_MS` (default: `64,128,256,384,512`)
- `BS_PUB_SHORT_MS` (default: `32,64,128,256,512`)
- `BS_PUB_SHORT_ASPECTS` (default: `2.0,4.0,8.0,10.0`)
- `BS_PUB_WARMUP` (default: `1`)
- `BS_PUB_SAMPLES` (default: `30`)
- `BS_NORM_RECOMP_TOL` (default: `sqrt(eps(Float64))`)
- `BS_PUB_CI_WARN_FRAC` (default: `0.5`)
- `BS_PUB_CI_FAIL_FRAC` (default: `10.0`)
- `BS_PUB_CI_MIN_TMED` (default: `1.0e-4`)
- `BS_PUB_CI_ENFORCE` (`1` to hard-fail high-CI rows, default: `0`)
- `BS_PUB_RESID_FACTOR` (default: `2500`)
- `BS_PUB_ORTH_FACTOR` (default: `2500`)
- `BS_SHORT_WIDE_FASTPATH` / `BS_SHORT_WIDE_FASTPATH_ASPECT` / `BS_SHORT_WIDE_FASTPATH_MMAX` / `BS_SHORT_WIDE_FASTPATH_NMIN`

Example run with explicit knobs:

```bash
BS_PUB_THREADS=1,4 BS_PUB_SEEDS=20260310,20260311 BS_PUB_WARMUP=1 BS_PUB_SAMPLES=30 julia --project=. benchmark/bench_cpqr_publication.jl
julia --project=. benchmark/plot_publication_results.jl benchmark/results/publication/publication_timings.csv benchmark/results/publication/plots benchmark/results/publication/tables
```

Outputs:

- `benchmark/results/publication/publication_timings.csv`
- `benchmark/results/publication/publication_summary.md`
- `benchmark/results/publication/metadata.txt`
- `benchmark/results/publication/plots/plain/*.png`
- `benchmark/results/publication/plots/plain/*.pdf`
- `benchmark/results/publication/plots/rinv/*.png`
- `benchmark/results/publication/plots/rinv/*.pdf`
- `benchmark/results/publication/tables/plain/*.csv`
- `benchmark/results/publication/tables/rinv/*.csv`
