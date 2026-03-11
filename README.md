# BSPivotQR

Julia implementation of Bischof-Stewart column-pivoted QR with an optional early-stop parameter `k`.

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

## Publication benchmark workflow

Primary benchmark script:

```bash
julia --project=. benchmark/bench_cpqr_publication.jl
```

Generate publication figures and caption-ready tables:

```bash
julia --project=. benchmark/plot_publication_results.jl
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
- `BS_HOUSEHOLDER_LAPACK_LARF` (`1` by default; set `0` to force the manual BLAS path)

Example run with explicit knobs:

```bash
BS_PUB_THREADS=1,4 BS_PUB_SEEDS=20260310,20260311 BS_PUB_WARMUP=1 BS_PUB_SAMPLES=30 julia --project=. benchmark/bench_cpqr_publication.jl
julia --project=. benchmark/plot_publication_results.jl benchmark/results/publication/publication_timings.csv benchmark/results/publication/plots benchmark/results/publication/tables
```

Outputs:

- `benchmark/results/publication/publication_timings.csv`
- `benchmark/results/publication/publication_summary.md`
- `benchmark/results/publication/metadata.txt`
- `benchmark/results/publication/plots/*.png`
- `benchmark/results/publication/plots/*.pdf`
- `benchmark/results/publication/tables/*.csv`
