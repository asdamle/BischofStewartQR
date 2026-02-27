# CPQR Benchmark Summary

Generated: 2026-02-25T13:18:13.202

BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)

| family | m | n | k | method | median (s) | min (s) | speedup vs dgeqp3 | residual | orthogonality |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|
| gaussian | 64 | 64 | 64 | bsqr_full | 6.529e-5 | 5.567e-5 | 0.4978 | 5.27e-16 | 6.286e-15 |
| gaussian | 64 | 64 | 64 | dgeqp3 | 3.25e-5 | 3.017e-5 | 1.0 | 5.898e-16 | 6.792e-15 |
| gaussian | 64 | 128 | 64 | bsqr_full | 0.0001205 | 0.0001182 | 0.5244 | 6.309e-16 | 6.507e-15 |
| gaussian | 64 | 128 | 64 | dgeqp3 | 6.317e-5 | 6.296e-5 | 1.0 | 6.971e-16 | 6.652e-15 |
| gaussian | 128 | 64 | 64 | bsqr_full | 8.588e-5 | 8.108e-5 | 0.7889 | 4.841e-16 | 8.273e-15 |
| gaussian | 128 | 64 | 64 | dgeqp3 | 6.775e-5 | 6.696e-5 | 1.0 | 5.537e-16 | 4.942e-15 |
| ill_conditioned | 64 | 64 | 64 | bsqr_full | 5.125e-5 | 4.758e-5 | 0.6349 | 4.121e-16 | 6.391e-15 |
| ill_conditioned | 64 | 64 | 64 | dgeqp3 | 3.254e-5 | 3.021e-5 | 1.0 | 4.104e-16 | 6.336e-15 |
| ill_conditioned | 64 | 128 | 64 | bsqr_full | 0.0001204 | 0.00012 | 0.5742 | 4.962e-16 | 5.902e-15 |
| ill_conditioned | 64 | 128 | 64 | dgeqp3 | 6.912e-5 | 6.504e-5 | 1.0 | 5.386e-16 | 6.324e-15 |
| ill_conditioned | 128 | 64 | 64 | bsqr_full | 7.925e-5 | 7.483e-5 | 0.8607 | 4.702e-16 | 8.403e-15 |
| ill_conditioned | 128 | 64 | 64 | dgeqp3 | 6.821e-5 | 6.579e-5 | 1.0 | 4.517e-16 | 4.755e-15 |
| orthonormal_rows | 64 | 128 | 64 | bsqr_full | 0.0001219 | 0.0001201 | 0.5121 | 6.673e-16 | 5.997e-15 |
| orthonormal_rows | 64 | 128 | 64 | dgeqp3 | 6.242e-5 | 6.042e-5 | 1.0 | 7.931e-16 | 6.886e-15 |
| orthonormal_rows | 64 | 256 | 64 | bsqr_full | 0.0003046 | 0.0002982 | 0.422 | 6.712e-16 | 5.99e-15 |
| orthonormal_rows | 64 | 256 | 64 | dgeqp3 | 0.0001285 | 0.0001265 | 1.0 | 7.585e-16 | 6.359e-15 |
| orthonormal_rows | 64 | 512 | 64 | bsqr_full | 0.0009596 | 0.0008896 | 0.4194 | 7.479e-16 | 6.463e-15 |
| orthonormal_rows | 64 | 512 | 64 | dgeqp3 | 0.0004025 | 0.0003989 | 1.0 | 7.817e-16 | 6.496e-15 |
