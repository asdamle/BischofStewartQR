# Publication Benchmark Summary

- Run ID: `20260316_121419`
- Generated: 2026-03-16T12:32:58.713
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Bench Surface: factor_only
- Threads: 1, 4
- Seeds: 20260310, 20260311

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr_full) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.49784 |
| gaussian | short_wide | 4 | 0.47967 |
| gaussian | square | 1 | 0.50265 |
| gaussian | square | 4 | 0.50165 |
| ill_conditioned | short_wide | 1 | 0.52873 |
| ill_conditioned | short_wide | 4 | 0.54352 |
| ill_conditioned | square | 1 | 0.56367 |
| ill_conditioned | square | 4 | 0.54399 |
| orthonormal_rows | short_wide | 1 | 0.486 |
| orthonormal_rows | short_wide | 4 | 0.48179 |
| orthonormal_rows | square | 1 | 0.49996 |
| orthonormal_rows | square | 4 | 0.51005 |

| family | regime | blas_threads | geomean speedup (dgeqp3_trsm/bsqr_rinv) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.6075 |
| gaussian | short_wide | 4 | 0.593 |
| gaussian | square | 1 | 0.51406 |
| gaussian | square | 4 | 0.53752 |
| ill_conditioned | short_wide | 1 | 0.66326 |
| ill_conditioned | short_wide | 4 | 0.64246 |
| ill_conditioned | square | 1 | 0.51683 |
| ill_conditioned | square | 4 | 0.56676 |
| orthonormal_rows | short_wide | 1 | 0.60155 |
| orthonormal_rows | short_wide | 4 | 0.57738 |
| orthonormal_rows | square | 1 | 0.52262 |
| orthonormal_rows | square | 4 | 0.54041 |
