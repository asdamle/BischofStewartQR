# Publication Benchmark Summary

- Run ID: `20260316_151728`
- Generated: 2026-03-16T15:37:49.941
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1, 4
- Seeds: 20260310, 20260311

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr_full) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.54542 |
| gaussian | short_wide | 4 | 0.54231 |
| gaussian | square | 1 | 0.62455 |
| gaussian | square | 4 | 0.61167 |
| ill_conditioned | short_wide | 1 | 0.58259 |
| ill_conditioned | short_wide | 4 | 0.58489 |
| ill_conditioned | square | 1 | 0.6545 |
| ill_conditioned | square | 4 | 0.65024 |
| orthonormal_rows | short_wide | 1 | 0.53751 |
| orthonormal_rows | short_wide | 4 | 0.54678 |
| orthonormal_rows | square | 1 | 0.61391 |
| orthonormal_rows | square | 4 | 0.62869 |

| family | regime | blas_threads | geomean speedup (dgeqp3_trsm/bsqr_rinv) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.63235 |
| gaussian | short_wide | 4 | 0.63381 |
| gaussian | square | 1 | 0.62276 |
| gaussian | square | 4 | 0.63308 |
| ill_conditioned | short_wide | 1 | 0.68982 |
| ill_conditioned | short_wide | 4 | 0.67749 |
| ill_conditioned | square | 1 | 0.64296 |
| ill_conditioned | square | 4 | 0.65069 |
| orthonormal_rows | short_wide | 1 | 0.63204 |
| orthonormal_rows | short_wide | 4 | 0.62719 |
| orthonormal_rows | square | 1 | 0.60032 |
| orthonormal_rows | square | 4 | 0.62395 |
