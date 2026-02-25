# CPQR Benchmark Summary

Generated: 2026-02-25T10:57:06.641

BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)

| family | m | n | k | method | median (s) | min (s) | speedup vs dgeqp3 | residual | orthogonality |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| gaussian | 64 | 64 | 64 | bsqr_full | 5.946e-5 | 5.05e-5 | 0.7358 | 5.27e-16 | 6.286e-15 |
| gaussian | 64 | 64 | 64 | dgeqp3 | 4.375e-5 | 3.088e-5 | 1.0 | 5.898e-16 | 6.792e-15 |
| gaussian | 64 | 128 | 64 | bsqr_full | 0.0001097 | 0.0001088 | 0.5752 | 6.309e-16 | 6.507e-15 |
| gaussian | 64 | 128 | 64 | dgeqp3 | 6.308e-5 | 6.267e-5 | 1.0 | 6.971e-16 | 6.652e-15 |
| gaussian | 128 | 64 | 64 | bsqr_full | 7.825e-5 | 7.475e-5 | 1.024 | 4.841e-16 | 8.273e-15 |
| gaussian | 128 | 64 | 64 | dgeqp3 | 8.013e-5 | 7.2e-5 | 1.0 | 5.537e-16 | 4.942e-15 |
| ill_conditioned | 64 | 64 | 64 | bsqr_full | 5.213e-5 | 4.796e-5 | 0.8473 | 4.121e-16 | 6.391e-15 |
| ill_conditioned | 64 | 64 | 64 | dgeqp3 | 4.417e-5 | 3.729e-5 | 1.0 | 4.104e-16 | 6.336e-15 |
| ill_conditioned | 64 | 128 | 64 | bsqr_full | 0.000182 | 0.0001347 | 0.3846 | 4.962e-16 | 5.902e-15 |
| ill_conditioned | 64 | 128 | 64 | dgeqp3 | 7.0e-5 | 6.929e-5 | 1.0 | 5.386e-16 | 6.324e-15 |
| ill_conditioned | 128 | 64 | 64 | bsqr_full | 7.858e-5 | 7.463e-5 | 0.9512 | 4.702e-16 | 8.403e-15 |
| ill_conditioned | 128 | 64 | 64 | dgeqp3 | 7.475e-5 | 7.238e-5 | 1.0 | 4.517e-16 | 4.755e-15 |
| orthonormal_rows | 64 | 128 | 64 | bsqr_full | 0.0001113 | 0.0001065 | 0.5489 | 6.673e-16 | 5.997e-15 |
| orthonormal_rows | 64 | 128 | 64 | dgeqp3 | 6.108e-5 | 5.612e-5 | 1.0 | 7.931e-16 | 6.886e-15 |
| orthonormal_rows | 64 | 256 | 64 | bsqr_full | 0.0002952 | 0.0002924 | 0.4236 | 6.712e-16 | 5.99e-15 |
| orthonormal_rows | 64 | 256 | 64 | dgeqp3 | 0.000125 | 0.0001199 | 1.0 | 7.585e-16 | 6.359e-15 |
| orthonormal_rows | 64 | 512 | 64 | bsqr_full | 0.0008662 | 0.0008589 | 0.4717 | 7.479e-16 | 6.463e-15 |
| orthonormal_rows | 64 | 512 | 64 | dgeqp3 | 0.0004086 | 0.0003946 | 1.0 | 7.817e-16 | 6.496e-15 |
