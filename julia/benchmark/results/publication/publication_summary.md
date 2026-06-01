# Publication Benchmark Summary

- Run ID: `20260316_151728`
- Generated: 2026-03-16T15:37:49.941
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1, 4
- Seeds: 20260310, 20260311

Relative time is BSQR median time divided by baseline median time; `1.0` is parity.

| family | regime | blas_threads | geomean relative time (bsqr_full/dgeqp3) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 1.8335 |
| gaussian | short_wide | 4 | 1.844 |
| gaussian | square | 1 | 1.6012 |
| gaussian | square | 4 | 1.6349 |
| ill_conditioned | short_wide | 1 | 1.7165 |
| ill_conditioned | short_wide | 4 | 1.7097 |
| ill_conditioned | square | 1 | 1.5279 |
| ill_conditioned | square | 4 | 1.5379 |
| orthonormal_rows | short_wide | 1 | 1.8604 |
| orthonormal_rows | short_wide | 4 | 1.8289 |
| orthonormal_rows | square | 1 | 1.6289 |
| orthonormal_rows | square | 4 | 1.5906 |

| family | regime | blas_threads | geomean relative time (bsqr_rinv/dgeqp3_trsm) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 1.5814 |
| gaussian | short_wide | 4 | 1.5778 |
| gaussian | square | 1 | 1.6057 |
| gaussian | square | 4 | 1.5796 |
| ill_conditioned | short_wide | 1 | 1.4497 |
| ill_conditioned | short_wide | 4 | 1.476 |
| ill_conditioned | square | 1 | 1.5553 |
| ill_conditioned | square | 4 | 1.5368 |
| orthonormal_rows | short_wide | 1 | 1.5822 |
| orthonormal_rows | short_wide | 4 | 1.5944 |
| orthonormal_rows | square | 1 | 1.6658 |
| orthonormal_rows | square | 4 | 1.6027 |
