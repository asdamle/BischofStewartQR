# Publication Benchmark Summary

- Run ID: `20260316_123334`
- Generated: 2026-03-16T12:54:35.572
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Bench Surface: materialize_qrp
- Threads: 1, 4
- Seeds: 20260310, 20260311

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr_full) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.56701 |
| gaussian | short_wide | 4 | 0.53309 |
| gaussian | square | 1 | 0.61411 |
| gaussian | square | 4 | 0.62747 |
| ill_conditioned | short_wide | 1 | 0.60532 |
| ill_conditioned | short_wide | 4 | 0.57094 |
| ill_conditioned | square | 1 | 0.65034 |
| ill_conditioned | square | 4 | 0.63209 |
| orthonormal_rows | short_wide | 1 | 0.53092 |
| orthonormal_rows | short_wide | 4 | 0.51489 |
| orthonormal_rows | square | 1 | 0.59613 |
| orthonormal_rows | square | 4 | 0.58427 |

| family | regime | blas_threads | geomean speedup (dgeqp3_trsm/bsqr_rinv) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.63023 |
| gaussian | short_wide | 4 | 0.61089 |
| gaussian | square | 1 | 0.60405 |
| gaussian | square | 4 | 0.57985 |
| ill_conditioned | short_wide | 1 | 0.68027 |
| ill_conditioned | short_wide | 4 | 0.66243 |
| ill_conditioned | square | 1 | 0.65729 |
| ill_conditioned | square | 4 | 0.61651 |
| orthonormal_rows | short_wide | 1 | 0.61871 |
| orthonormal_rows | short_wide | 4 | 0.61958 |
| orthonormal_rows | square | 1 | 0.59462 |
| orthonormal_rows | square | 4 | 0.60141 |
