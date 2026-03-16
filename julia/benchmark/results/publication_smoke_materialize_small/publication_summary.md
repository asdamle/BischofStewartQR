# Publication Benchmark Summary

- Run ID: `20260316_120835`
- Generated: 2026-03-16T12:08:37.472
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Bench Surface: materialize_qrp
- Threads: 1
- Seeds: 20260310

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr_full) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.70973 |
| gaussian | square | 1 | 0.57642 |

| family | regime | blas_threads | geomean speedup (dgeqp3_trsm/bsqr_rinv) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.92455 |
| gaussian | square | 1 | 0.80262 |
