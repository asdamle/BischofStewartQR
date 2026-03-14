# MATLAB MEX Backend (Optional)

The MATLAB BSQR implementation defaults to the pure M-file backend.

## Build

```matlab
addpath('matlab')
build_bsqr_mex
```

`build_bsqr_mex` looks for C/C++ files in `matlab/mex/src/` and builds `bsqr_mex`.

## Promotion Rule

Promote MEX to default only when benchmark evidence shows at least 20% median
speedup on short-wide focus cases with no residual/orthogonality regression.
