# CPQR Benchmark Summary

Generated: 2026-02-24T23:54:41.745

BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)

| family | m | n | k | method | median (s) | min (s) | speedup vs dgeqp3 | residual | orthogonality |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| gaussian | 64 | 64 | 64 | bsqr_full | 8.233e-5 | 5.967e-5 | 0.3932 | 5.27e-16 | 6.286e-15 |
| gaussian | 64 | 64 | 64 | dgeqp3 | 3.238e-5 | 3.012e-5 | 1.0 | 5.898e-16 | 6.792e-15 |
| gaussian | 64 | 128 | 64 | bsqr_full | 0.0001254 | 0.0001233 | 0.5072 | 6.309e-16 | 6.507e-15 |
| gaussian | 64 | 128 | 64 | dgeqp3 | 6.358e-5 | 6.229e-5 | 1.0 | 6.971e-16 | 6.652e-15 |
| gaussian | 128 | 64 | 64 | bsqr_full | 8.767e-5 | 8.517e-5 | 0.759 | 4.841e-16 | 8.273e-15 |
| gaussian | 128 | 64 | 64 | dgeqp3 | 6.654e-5 | 6.521e-5 | 1.0 | 5.537e-16 | 4.942e-15 |
| ill_conditioned | 64 | 64 | 64 | bsqr_full | 6.058e-5 | 5.608e-5 | 0.5935 | 4.121e-16 | 6.391e-15 |
| ill_conditioned | 64 | 64 | 64 | dgeqp3 | 3.596e-5 | 3.421e-5 | 1.0 | 4.104e-16 | 6.336e-15 |
| ill_conditioned | 64 | 128 | 64 | bsqr_full | 0.0001233 | 0.0001192 | 0.5375 | 4.962e-16 | 5.902e-15 |
| ill_conditioned | 64 | 128 | 64 | dgeqp3 | 6.629e-5 | 6.221e-5 | 1.0 | 5.386e-16 | 6.324e-15 |
| ill_conditioned | 128 | 64 | 64 | bsqr_full | 8.575e-5 | 8.096e-5 | 0.81 | 4.702e-16 | 8.403e-15 |
| ill_conditioned | 128 | 64 | 64 | dgeqp3 | 6.946e-5 | 6.671e-5 | 1.0 | 4.517e-16 | 4.755e-15 |
| orthonormal_rows | 64 | 128 | 64 | bsqr_full | 0.0001279 | 0.0001245 | 0.5487 | 6.673e-16 | 5.997e-15 |
| orthonormal_rows | 64 | 128 | 64 | dgeqp3 | 7.017e-5 | 6.9e-5 | 1.0 | 7.931e-16 | 6.886e-15 |
| orthonormal_rows | 64 | 256 | 64 | bsqr_full | 0.0003205 | 0.0003128 | 0.4517 | 6.712e-16 | 5.99e-15 |
| orthonormal_rows | 64 | 256 | 64 | dgeqp3 | 0.0001448 | 0.0001378 | 1.0 | 7.585e-16 | 6.359e-15 |
| orthonormal_rows | 64 | 512 | 64 | bsqr_full | 0.0009648 | 0.0009501 | 0.4401 | 7.479e-16 | 6.463e-15 |
| orthonormal_rows | 64 | 512 | 64 | dgeqp3 | 0.0004246 | 0.0004096 | 1.0 | 7.817e-16 | 6.496e-15 |
