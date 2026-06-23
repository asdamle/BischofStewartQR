# Publication Benchmark Summary

- Run ID: `20260622_161419`
- Generated: 2026-06-22T17:00:47.390
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1, 4
- Seeds: 20260310, 20260311, 20260312, 20260313, 20260314


Relative time is BSQR median time divided by baseline median time; `1.0` is parity.

| family | regime | blas_threads | geomean relative time (bsqr_full/dgeqp3) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 1.492 |
| gaussian | short_wide | 4 | 1.5123 |
| gaussian | square | 1 | 1.2879 |
| gaussian | square | 4 | 1.2576 |
| ill_conditioned | short_wide | 1 | 1.394 |
| ill_conditioned | short_wide | 4 | 1.3901 |
| ill_conditioned | square | 1 | 1.1983 |
| ill_conditioned | square | 4 | 1.1861 |
| orthonormal_rows | short_wide | 1 | 1.507 |
| orthonormal_rows | short_wide | 4 | 1.5067 |
| orthonormal_rows | square | 1 | 1.2628 |
| orthonormal_rows | square | 4 | 1.2586 |

| family | regime | blas_threads | geomean relative time (bsqr_rinv/dgeqp3_trsm) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 1.2938 |
| gaussian | short_wide | 4 | 1.2985 |
| gaussian | square | 1 | 1.2155 |
| gaussian | square | 4 | 1.2539 |
| ill_conditioned | short_wide | 1 | 1.1874 |
| ill_conditioned | short_wide | 4 | 1.2086 |
| ill_conditioned | square | 1 | 1.1847 |
| ill_conditioned | square | 4 | 1.2029 |
| orthonormal_rows | short_wide | 1 | 1.3077 |
| orthonormal_rows | short_wide | 4 | 1.2984 |
| orthonormal_rows | square | 1 | 1.2713 |
| orthonormal_rows | square | 4 | 1.2555 |
