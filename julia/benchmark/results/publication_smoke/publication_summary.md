# Publication Benchmark Summary

- Run ID: `20260314_160752`
- Generated: 2026-03-14T16:07:55.186
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1
- Seeds: 20260310

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr_full) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.49795 |
| gaussian | square | 1 | 0.74556 |

| family | regime | blas_threads | geomean speedup (dgeqp3_trsm/bsqr_rinv) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.95653 |
| gaussian | square | 1 | 0.65215 |
