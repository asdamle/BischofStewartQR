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

Promote MEX to default only when benchmark evidence shows at least 20% median
speedup on short-wide focus cases with no residual/orthogonality regression.
