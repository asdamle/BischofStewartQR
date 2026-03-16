# Publication Benchmark Summary

- Run ID: `20260316_120942`
- Generated: 2026-03-16T12:09:45.458
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Bench Surface: factor_only
- Threads: 1
- Seeds: 20260310

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr_full) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.68521 |
| gaussian | square | 1 | 0.4786 |

| family | regime | blas_threads | geomean speedup (dgeqp3_trsm/bsqr_rinv) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.75905 |
| gaussian | square | 1 | 0.64919 |
