# Publication Benchmark Summary

- Run ID: `20260312_115647`
- Generated: 2026-03-12T11:58:37.228
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1
- Seeds: 20260310

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr_full) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.4489 |
| gaussian | square | 1 | 0.47834 |
| ill_conditioned | short_wide | 1 | 0.54558 |
| ill_conditioned | square | 1 | 0.60594 |
| orthonormal_rows | short_wide | 1 | 0.48046 |
| orthonormal_rows | square | 1 | 0.51195 |

| family | regime | blas_threads | geomean speedup (dgeqp3_trsm/bsqr_rinv) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.56607 |
| gaussian | square | 1 | 0.51979 |
| ill_conditioned | short_wide | 1 | 0.65843 |
| ill_conditioned | square | 1 | 0.61333 |
| orthonormal_rows | short_wide | 1 | 0.59318 |
| orthonormal_rows | square | 1 | 0.51443 |
