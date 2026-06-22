# BSQR — Bischof–Stewart column-pivoted QR

Reference implementations of Bischof–Stewart (BS) column-pivoted QR, plus an
experimental randomized variant. This repository accompanies a forthcoming paper.

- **MATLAB** — deterministic `bsqr` (`matlab/`) and randomized `bsqr_rand`
  (`matlab_rand/`), each with a fast C++ MEX backend.
- **Julia** — deterministic `BSPivotQR` package (`julia/`).

## Quick start

### MATLAB

Launch MATLAB from this folder (or run `startup`); it puts the kernels on the
path and builds both MEX backends, then the default calls just work:

```matlab
startup                                  % add paths + build MEX (run once)
R = bsqr(A);                             % deterministic BS pivoted QR
[Q, R, p] = bsqr(A);
[p, reflectors, R11] = bsqr_rand(M);     % randomized variant (M has orthonormal rows)
```

Both default to the compiled MEX backend; pass `'backend','mfile'` for the
pure-MATLAB reference. See [matlab/README.md](matlab/README.md) and
[matlab_rand/README.md](matlab_rand/README.md).

### Julia

```julia
julia> include("startup.jl")             # activate + instantiate julia/, then `using BSPivotQR`
julia> F = bsqr(A);                      # deterministic BS pivoted QR -> BSQRPivoted
julia> R(F), perm(F)                     # R factor, column permutation
```

or `julia --project=julia` then `using BSPivotQR`. (The randomized variant is
MATLAB-only.) See [julia/README.md](julia/README.md).

## Layout

- `matlab/` — deterministic `bsqr` (m-file + MEX), tests, publication benchmarks.
- `matlab_rand/` — randomized `bsqr_rand` (m-file + MEX), tests, benchmarks.
- `julia/` — `BSPivotQR` package (depends only on `LinearAlgebra`);
  `julia/benchmark/` is a separate environment for the benchmark pipeline.
- `docs/` — design notes and validation.
- `startup.m` / `startup.jl` — one-shot setup for MATLAB / Julia.

See [CLAUDE.md](CLAUDE.md) for the full command set (tests, benchmarks, builds).
