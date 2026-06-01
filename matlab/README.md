# MATLAB BSQR

This directory provides a MATLAB implementation of Bischof-Stewart pivoted QR.

## API

```matlab
R = bsqr(A)
[Q,R] = bsqr(A)
[Q,R,E_or_p] = bsqr(A, 'pivot_format', 'matrix' | 'vector')
[Q,R,E_or_p,R11invR12] = bsqr(A, 'return_rinv_r12', true)
```

Name-value options:

- `'k'`: early-stop step count (default `min(size(A))`)
- `'return_rinv_r12'`: default `false`
- `'pivot_format'`: `'matrix'` (default) or `'vector'`
- `'backend'`: `'auto'` (default), `'mfile'`, `'mex'`
- `'norm_recomp_tol'`: running norm refresh tolerance (default `sqrt(eps)`)
- `'check_finite'`: validate finite input (default `true`)

`backend='auto'` uses `bsqr_mex` when available, otherwise the pure MATLAB kernel.

## Tests

```matlab
addpath('matlab')
run_tests
```

or from shell:

```bash
matlab -batch "addpath('matlab'); run_tests"
```

## Benchmarks

Run publication benchmark + artifact generation:

```matlab
addpath('matlab')
addpath('matlab/benchmark')
run_publication_benchmarks
```

Run a fast smoke benchmark (small cases, warmup enabled):

```matlab
addpath('matlab')
addpath('matlab/benchmark')
run_publication_smoke_benchmark
```

Benchmark contract:

- Built-in baseline uses economy pivoted QR with vector permutation (`qr(...,'econ','vector')`, with compatibility fallback to `qr(A,0,'vector')`).
- BSQR baseline uses the MEX backend directly (`bsqr_mex`) with vector permutation.
- Publication runner is mex-only and does not accept `cfg.bsqr_backend`.
- Timed publication path materializes `Q,R,p` for both BSQR and built-in QR.

Regenerate plots/tables from CSV:

```matlab
addpath('matlab/benchmark')
plot_publication_results
```

Publication plots default to PNG. Set `BS_MATLAB_PUB_FIG_FORMATS=png,pdf,eps`
to request vector outputs as well. Derived performance plots and tables report
relative time (`BSQR / baseline`), where `1.0` is parity.

Compare candidate timings against a baseline with slowdown gating:

```matlab
addpath('matlab/benchmark')
check_publication_perf_gate('baseline.csv', 'candidate.csv', 0.05)
```

Environment knobs (optional):

- `BS_MATLAB_PUB_OUTDIR`
- `BS_MATLAB_PUB_SEEDS`
- `BS_MATLAB_PUB_FAMILIES`
- `BS_MATLAB_PUB_SQUARE_MS`
- `BS_MATLAB_PUB_SHORT_MS`
- `BS_MATLAB_PUB_SHORT_ASPECTS`
- `BS_MATLAB_PUB_WARMUP`
- `BS_MATLAB_PUB_SAMPLES`
- `BS_MATLAB_NORM_RECOMP_TOL`
- `BS_MATLAB_PUB_FIG_FORMATS` (`png` default; supports `png,pdf,eps`)
- `BS_MATLAB_ALLOW_SHARED_OUTDIR` (`0` default). Keep at `0` to prevent MATLAB
  benchmark artifacts from writing into Julia benchmark roots
  (`benchmark/results/*` and `julia/benchmark/results/*`).

## Optional MEX

```matlab
addpath('matlab')
build_bsqr_mex
```

Add C/C++ sources under `matlab/mex/src/`.

- Publication benchmarks require `bsqr_mex`.
- The pure `.m` implementation remains available for correctness testing and development.
