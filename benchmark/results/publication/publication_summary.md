# Publication Benchmark Summary

- Run ID: `20260312_122004`
- Generated: 2026-03-12T12:38:08.146
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1, 4
- Seeds: 20260310, 20260311

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr_full) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.49188 |
| gaussian | short_wide | 4 | 0.4889 |
| gaussian | square | 1 | 0.49446 |
| gaussian | square | 4 | 0.50761 |
| ill_conditioned | short_wide | 1 | 0.54222 |
| ill_conditioned | short_wide | 4 | 0.52648 |
| ill_conditioned | square | 1 | 0.55365 |
| ill_conditioned | square | 4 | 0.5136 |
| orthonormal_rows | short_wide | 1 | 0.48882 |
| orthonormal_rows | short_wide | 4 | 0.47794 |
| orthonormal_rows | square | 1 | 0.51638 |
| orthonormal_rows | square | 4 | 0.51687 |

| family | regime | blas_threads | geomean speedup (dgeqp3_trsm/bsqr_rinv) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.60378 |
| gaussian | short_wide | 4 | 0.58605 |
| gaussian | square | 1 | 0.5153 |
| gaussian | square | 4 | 0.51308 |
| ill_conditioned | short_wide | 1 | 0.6566 |
| ill_conditioned | short_wide | 4 | 0.64053 |
| ill_conditioned | square | 1 | 0.55661 |
| ill_conditioned | square | 4 | 0.55371 |
| orthonormal_rows | short_wide | 1 | 0.603 |
| orthonormal_rows | short_wide | 4 | 0.58647 |
| orthonormal_rows | square | 1 | 0.51914 |
| orthonormal_rows | square | 4 | 0.51684 |
