# MATLAB MEX Backend

## Build

```matlab
addpath('matlab')
build_bsqr_mex
```

`build_bsqr_mex` looks for C/C++ files in `matlab/mex/src/` and builds `bsqr_mex`.

Publication benchmarks require `bsqr_mex`. The pure MATLAB backend remains
available for correctness checks and development workflows.

## Promotion Rule

`backend='auto'` (the dispatcher default) uses `bsqr_mex` whenever it has been
built and falls back to the pure-`.m` kernel otherwise; `startup` builds the
MEX, so it is the de-facto default backend while remaining optional. Keeping
that preference requires benchmark evidence: at least 20% median speedup on
short-wide focus cases with no residual/orthogonality regression. Status: to be
evaluated against the final publication rerun (Phase 6 of
`docs/PUBLICATION_READINESS_PLAN.md`).
