# MATLAB BSQR

This directory provides a MATLAB implementation of Bischof-Stewart pivoted QR.

## API

```matlab
R = bsqr(A)
[Q,R] = bsqr(A)
[Q,R,E_or_p] = bsqr(A, 'pivot_format', 'matrix' | 'vector')
[Q,R,E_or_p,R11invR12] = bsqr(A, 'return_rinv_r12', true)
```

Outputs, for `A` `m×n` and `k` factorization steps (default `k = min(m,n)`):

- `R`: `k×n` upper-trapezoidal factor of the *permuted* columns; `R(1:k,1:k)` is
  `R11`, `R(:,k+1:n)` is `R12`. For the default `k = min(m,n)`,
  `A(:,p) = Q*R` (equivalently `A*E = Q*R`); with early stop (`k < min(m,n)`)
  that holds only for the selected block, `A(:,p(1:k)) = Q*R(:,1:k)`, and
  `R12 = Q'*A(:,p(k+1:n))` is the unselected columns' projection onto
  `span(Q)`. With a single output the MEX backend skips materializing `Q`.
- `Q`: `m×k` economy factor with orthonormal columns.
- `E_or_p`: the permutation — an `n×n` permutation matrix `E` with
  `A*E = A(:,p)` (`'matrix'`, the default, matching `qr(A)`'s three-output
  form), or the `1×n` index row `p` itself (`'vector'`).
- `R11invR12`: `k×(n−k)` matrix `R11^{-1}R12`, read directly from the kernel
  workspace (no extra triangular solve); `[]` unless `'return_rinv_r12'` is true.
- `trace` (5th output, requires `'trace', true`): per-step validation trace,
  see the option below.

Name-value options:

- `'k'`: early-stop step count (default `min(size(A))`)
- `'return_rinv_r12'`: default `false`
- `'pivot_format'`: `'matrix'` (default) or `'vector'`
- `'backend'`: `'auto'` (default), `'mfile'`, `'mex'`
- `'norm_recomp_tol'`: running norm refresh tolerance (default `sqrt(eps)`)
- `'check_finite'`: validate finite input (default `true`)
- `'trace'`: enable the per-step validation trace as a fifth output
  (`[Q,R,p,rinv,trace] = bsqr(A,'trace',true)`; struct with fields `crit`
  and `nrecomp`; see `docs/VALIDATION.md` V4)

`backend='auto'` uses `bsqr_mex` when available, otherwise the pure MATLAB kernel.

## Tests

```matlab
addpath('matlab/tests')
run_tests
```

or from shell:

```bash
matlab -batch "addpath('matlab/tests'); run_tests"
```

Cross-language parity fixtures live in `parity/` (oracle outputs consumed by both the
MATLAB and Julia test suites). Regenerate them only when `matlab/tests/parity_zoo.m`
changes, then rerun both suites:

```bash
matlab -batch "addpath('matlab/tests'); generate_parity_fixtures"
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
- `BS_MATLAB_PUB_FAMILIES` (constructions defined in `docs/PUBLICATION_FIGURES_PLAN.md`,
  identical to the Julia generators)
- `BS_MATLAB_PUB_SQUARE_MS`
- `BS_MATLAB_PUB_SHORT_MS`
- `BS_MATLAB_PUB_SHORT_ASPECTS`
- `BS_MATLAB_PUB_WARMUP`
- `BS_MATLAB_PUB_SAMPLES`
- `BS_MATLAB_NORM_RECOMP_TOL`
- `BS_MATLAB_PUB_FIG_FORMATS` (`png` default; supports `png,pdf,eps`)
- `BS_PANEL_NB` (default `8`; panel/blocked MEX kernel width, `0`/`1` selects
  the unblocked kernel — kept in lockstep with the Julia kernel; see
  `docs/P3_BLOCKED_BSQR.md`)
- `BS_PANEL_MIN_KN` (default `24576`; minimum `k*n` for the panel kernel —
  below the crossover the unblocked kernel runs instead)
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
