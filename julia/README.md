# Julia BSQR

Canonical Julia implementation of Bischof-Stewart column-pivoted QR.

## Setup

```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate()'
```

## Tests

```bash
julia --project=julia julia/test/runtests.jl
```

## Benchmarks

Publication benchmark:

```bash
julia --project=julia julia/benchmark/run_publication_benchmarks.jl
```

Smoke benchmark:

```bash
julia --project=julia julia/benchmark/run_publication_smoke_benchmark.jl
```

Plot from CSV:

```bash
julia --project=julia julia/benchmark/plot_publication_results.jl
```

Performance gate:

```bash
julia --project=julia julia/benchmark/check_publication_perf_gate.jl
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
