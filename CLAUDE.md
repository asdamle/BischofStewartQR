# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

Two parallel implementations of Bischof-Stewart column-pivoted QR (BSQR), built for a publication: `julia/` (canonical, package `BSPivotQR`) and `matlab/` (pure `.m` kernel plus a C++ MEX backend). The original Bischof and Stewart papers (the ground-truth description of the algorithm) and the algorithm writeup (`bischof_stewart_pivoting.tex`) live in the local-only `notes/` folder at the repo root — it is gitignored, so consult it on disk, not in git history. Everything is Float64-only.

The two implementations are intentionally kept algorithmically identical: same pivot criterion `(1 + ||w_j||^2) / ||a_j||^2` with strict `<` tie-breaking (first minimum wins), same running column-norm downdates with a recompute tolerance (default `sqrt(eps)`), same workspace layout. Changes to kernel behavior in one language generally need a matching change in the other.

## Commands

### Julia

```bash
julia --project=julia -e 'using Pkg; Pkg.instantiate()'   # setup
julia --project=julia julia/test/runtests.jl              # tests
```

All tests are `@testset` blocks in `julia/test/test_bsqr.jl`; there is no per-test runner — comment out or run the file via `include` in a REPL to iterate on one testset.

Benchmarks (publication pipeline):

```bash
julia --project=julia julia/benchmark/run_publication_benchmarks.jl     # full run
julia --project=julia julia/benchmark/run_publication_smoke_benchmark.jl  # fast smoke
julia --project=julia julia/benchmark/plot_publication_results.jl       # plots/tables from CSV
julia --project=julia julia/benchmark/check_publication_perf_gate.jl    # perf gate
```

### MATLAB

```bash
matlab -batch "addpath('matlab/tests'); run_tests"            # tests
matlab -batch "addpath('matlab'); build_bsqr_mex"             # build MEX (sources in matlab/mex/src/)
matlab -batch "addpath('matlab'); addpath('matlab/benchmark'); run_publication_benchmarks"
matlab -batch "addpath('matlab'); addpath('matlab/benchmark'); run_publication_smoke_benchmark"
matlab -batch "addpath('matlab/benchmark'); plot_publication_results"
```

Tests are in `matlab/tests/test_bsqr.m` (script-based `runtests`); run a single file with `runtests('matlab/tests/test_bsqr.m')` after adding paths.

### Cross-language comparison

```bash
python3 julia/benchmark/compare_matlab_julia_timings.py \
  --matlab-csv matlab/benchmark/results/publication/publication_timings.csv \
  --julia-csv julia/benchmark/results/publication/publication_timings.csv \
  --method bsqr_full --julia-threads 1
```

## Architecture

### Julia package (`julia/src/`)

- `workspace.jl` — `BSWorkspace`: preallocated scratch buffers (`W`, running norms `s`/`s_ref`, `wnorm2`, etc.) so repeated factorizations are allocation-minimal.
- `kernel.jl` — the core. `BSQRPivoted <: Factorization` result type; three entry points: `bsqr(A)` (copies), `bsqr!(A; ...)` (in-place, allocates scratch), and `bsqr!(A, tau, jpvt, workspace; ...)` (fully preallocated, returns `ksteps`). `_bsqr_kernel!` is the main loop. The kernel incrementally maintains `W = R11^{-1} R12` in `ws.W` — that is both the pivot-criterion state and why `return_rinv_r12` is cheap (extracted directly from the workspace, no solve).
- `interface.jl` — accessors on `BSQRPivoted`: `R`, `perm`, `rinv_r12`, `reconstruct`. Q is stored implicitly as Householder vectors in `factors` + `tau` (LAPACK `ormqr` to materialize).

Kernel knobs that tests rely on: `k` (early stop), `rank_stop`, `norm_recomp_tol`, `blas_threads` (temporarily pins BLAS thread count), `track_inverse_frob`.

### MATLAB (`matlab/`)

- `bsqr.m` — user-facing dispatcher. `backend='auto'` (default) uses `bsqr_mex` when built, else falls back to `private/bsqr_mfile.m`. The m-file kernel mirrors the Julia kernel; the MEX (`mex/src/bsqr_mex.cpp`, BLAS/LAPACK calls, persistent workspace) mirrors it too.
- The pure `.m` backend exists for correctness testing and development; publication benchmarks are MEX-only and error without `bsqr_mex`.
- MEX promotion rule (from `matlab/mex/README.md`): MEX stays opt-in as default until benchmarks show ≥20% median speedup on short-wide focus cases with no residual/orthogonality regression.

### Publication benchmark pipeline

Both languages follow the same shape: runner → `results/publication/{publication_timings.csv, publication_summary.md, metadata.txt}` → plot script → `plots/{plain,rinv}/` and `tables/{plain,rinv}/` (`plain` = standard factorization, `rinv` = the `R11^{-1}R12`-returning variant). Latest publication artifacts are committed to git; regenerating them is part of a benchmark change.

Conventions:
- Derived plots/tables report **relative time** (`BSQR / baseline`), where `1.0` is parity.
- Baselines: built-in economy pivoted QR (`qr(...,'econ','vector')` in MATLAB, `qr(A, ColumnNorm())` in Julia); the timed path materializes `Q, R, p` for both sides.
- Configuration is via env vars: `BS_PUB_*` for Julia (see `julia/docs/TESTS_AND_BENCHMARKS.md` for the full list), `BS_MATLAB_PUB_*` for MATLAB (see `matlab/README.md`). Figures default to PNG; set `BS_PUB_FIG_FORMATS=png,pdf,eps` (or the `BS_MATLAB_` equivalent) for vector output.
- `BS_MATLAB_ALLOW_SHARED_OUTDIR=0` (default) prevents MATLAB benchmark artifacts from writing into Julia result roots — keep it that way.
- On macOS, Julia benchmarks use Apple Accelerate as the BLAS/LAPACK backend (`julia/benchmark/setup_accelerate.jl` verifies it loaded).
- Perf gates (`check_publication_perf_gate.*`) compare a candidate timings CSV against a baseline with a slowdown threshold; run after performance-relevant kernel changes.
