# CPQR Benchmark Summary

Generated: 2026-02-25T12:53:57.736

BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)

| family | m | n | k | method | median (s) | min (s) | speedup vs dgeqp3 | residual | orthogonality |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| gaussian | 64 | 64 | 64 | bsqr_full | 6.912e-5 | 5.292e-5 | 0.437 | 5.27e-16 | 6.286e-15 |
| gaussian | 64 | 64 | 64 | dgeqp3 | 3.021e-5 | 2.804e-5 | 1.0 | 5.898e-16 | 6.792e-15 |
| gaussian | 64 | 128 | 64 | bsqr_full | 0.0001251 | 0.0001215 | 0.6179 | 6.309e-16 | 6.507e-15 |
| gaussian | 64 | 128 | 64 | dgeqp3 | 7.729e-5 | 6.483e-5 | 1.0 | 6.971e-16 | 6.652e-15 |
| gaussian | 128 | 64 | 64 | bsqr_full | 8.288e-5 | 7.938e-5 | 0.8034 | 4.841e-16 | 8.273e-15 |
| gaussian | 128 | 64 | 64 | dgeqp3 | 6.658e-5 | 6.492e-5 | 1.0 | 5.537e-16 | 4.942e-15 |
| ill_conditioned | 64 | 64 | 64 | bsqr_full | 5.629e-5 | 5.367e-5 | 0.5833 | 4.121e-16 | 6.391e-15 |
| ill_conditioned | 64 | 64 | 64 | dgeqp3 | 3.283e-5 | 2.942e-5 | 1.0 | 4.104e-16 | 6.336e-15 |
| ill_conditioned | 64 | 128 | 64 | bsqr_full | 0.0001128 | 0.0001093 | 0.7136 | 4.962e-16 | 5.902e-15 |
| ill_conditioned | 64 | 128 | 64 | dgeqp3 | 8.046e-5 | 6.846e-5 | 1.0 | 5.386e-16 | 6.324e-15 |
| ill_conditioned | 128 | 64 | 64 | bsqr_full | 0.0001076 | 8.8e-5 | 0.6286 | 4.702e-16 | 8.403e-15 |
| ill_conditioned | 128 | 64 | 64 | dgeqp3 | 6.762e-5 | 6.583e-5 | 1.0 | 4.517e-16 | 4.755e-15 |
| orthonormal_rows | 64 | 128 | 64 | bsqr_full | 0.0001165 | 0.0001142 | 0.5395 | 6.673e-16 | 5.997e-15 |
| orthonormal_rows | 64 | 128 | 64 | dgeqp3 | 6.283e-5 | 6.062e-5 | 1.0 | 7.931e-16 | 6.886e-15 |
| orthonormal_rows | 64 | 256 | 64 | bsqr_full | 0.0003078 | 0.0003066 | 0.4181 | 6.712e-16 | 5.99e-15 |
| orthonormal_rows | 64 | 256 | 64 | dgeqp3 | 0.0001287 | 0.0001239 | 1.0 | 7.585e-16 | 6.359e-15 |
| orthonormal_rows | 64 | 512 | 64 | bsqr_full | 0.000867 | 0.0008527 | 0.4891 | 7.479e-16 | 6.463e-15 |
| orthonormal_rows | 64 | 512 | 64 | dgeqp3 | 0.0004241 | 0.0004048 | 1.0 | 7.817e-16 | 6.496e-15 |
