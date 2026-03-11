# CPQR Benchmark Summary

Generated: 2026-03-10T18:40:11.508

BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)

| family | m | n | k | method | median (s) | min (s) | speedup vs dgeqp3 | residual | orthogonality |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| gaussian | 64 | 64 | 64 | bsqr_full | 7.275e-5 | 5.442e-5 | 0.4502 | 5.27e-16 | 6.286e-15 |
| gaussian | 64 | 64 | 64 | bsqr_lazy_blas | 0.0003132 | 0.0002868 | 0.1045 | 5.27e-16 | 6.286e-15 |
| gaussian | 64 | 64 | 64 | dgeqp3 | 3.275e-5 | 3.054e-5 | 1.0 | 5.898e-16 | 6.792e-15 |
| gaussian | 64 | 128 | 64 | bsqr_full | 0.0001207 | 0.0001178 | 0.5319 | 6.309e-16 | 6.507e-15 |
| gaussian | 64 | 128 | 64 | bsqr_lazy_blas | 0.0013 | 0.001207 | 0.0494 | 6.309e-16 | 6.507e-15 |
| gaussian | 64 | 128 | 64 | dgeqp3 | 6.421e-5 | 6.25e-5 | 1.0 | 6.971e-16 | 6.652e-15 |
| gaussian | 128 | 64 | 64 | bsqr_full | 8.596e-5 | 8.121e-5 | 0.8245 | 4.841e-16 | 8.273e-15 |
| gaussian | 128 | 64 | 64 | bsqr_lazy_blas | 0.0003006 | 0.0002704 | 0.2358 | 4.841e-16 | 8.273e-15 |
| gaussian | 128 | 64 | 64 | dgeqp3 | 7.088e-5 | 6.788e-5 | 1.0 | 5.537e-16 | 4.942e-15 |
| ill_conditioned | 64 | 64 | 64 | bsqr_full | 5.662e-5 | 5.263e-5 | 0.5865 | 4.121e-16 | 6.391e-15 |
| ill_conditioned | 64 | 64 | 64 | bsqr_lazy_blas | 0.000407 | 0.0003607 | 0.08158 | 4.121e-16 | 6.391e-15 |
| ill_conditioned | 64 | 64 | 64 | dgeqp3 | 3.321e-5 | 3.067e-5 | 1.0 | 4.104e-16 | 6.336e-15 |
| ill_conditioned | 64 | 128 | 64 | bsqr_full | 0.0001384 | 0.000122 | 0.4863 | 4.962e-16 | 5.902e-15 |
| ill_conditioned | 64 | 128 | 64 | bsqr_lazy_blas | 0.001925 | 0.001852 | 0.03495 | 4.962e-16 | 5.902e-15 |
| ill_conditioned | 64 | 128 | 64 | dgeqp3 | 6.729e-5 | 6.475e-5 | 1.0 | 5.386e-16 | 6.324e-15 |
| ill_conditioned | 128 | 64 | 64 | bsqr_full | 9.383e-5 | 8.817e-5 | 1.03 | 4.702e-16 | 8.403e-15 |
| ill_conditioned | 128 | 64 | 64 | bsqr_lazy_blas | 0.0006 | 0.0004681 | 0.1611 | 4.702e-16 | 8.403e-15 |
| ill_conditioned | 128 | 64 | 64 | dgeqp3 | 9.667e-5 | 9.1e-5 | 1.0 | 4.517e-16 | 4.755e-15 |
| orthonormal_rows | 64 | 128 | 64 | bsqr_full | 0.0001172 | 0.0001138 | 0.6034 | 6.673e-16 | 5.997e-15 |
| orthonormal_rows | 64 | 128 | 64 | bsqr_lazy_blas | 0.0009664 | 0.0008935 | 0.07321 | 6.673e-16 | 5.997e-15 |
| orthonormal_rows | 64 | 128 | 64 | dgeqp3 | 7.075e-5 | 6.846e-5 | 1.0 | 7.931e-16 | 6.886e-15 |
| orthonormal_rows | 64 | 256 | 64 | bsqr_full | 0.0003289 | 0.0003068 | 0.4189 | 6.712e-16 | 5.99e-15 |
| orthonormal_rows | 64 | 256 | 64 | bsqr_lazy_blas | 0.002386 | 0.002286 | 0.05772 | 6.712e-16 | 5.99e-15 |
| orthonormal_rows | 64 | 256 | 64 | dgeqp3 | 0.0001378 | 0.0001295 | 1.0 | 7.585e-16 | 6.359e-15 |
| orthonormal_rows | 64 | 512 | 64 | bsqr_full | 0.0008839 | 0.0008803 | 0.6076 | 7.479e-16 | 6.463e-15 |
| orthonormal_rows | 64 | 512 | 64 | bsqr_lazy_blas | 0.004574 | 0.004455 | 0.1174 | 7.479e-16 | 6.463e-15 |
| orthonormal_rows | 64 | 512 | 64 | dgeqp3 | 0.000537 | 0.0005222 | 1.0 | 7.817e-16 | 6.496e-15 |
