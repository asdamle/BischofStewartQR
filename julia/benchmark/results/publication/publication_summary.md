# Publication Benchmark Summary

- Run ID: `20260611_204425`
- Generated: 2026-06-11T21:02:57.548
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1, 4
- Seeds: 20260310, 20260311


Relative time is BSQR median time divided by baseline median time; `1.0` is parity.

| family | regime | blas_threads | geomean relative time (bsqr_full/dgeqp3) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 1.5668 |
| gaussian | short_wide | 4 | 1.5228 |
| gaussian | square | 1 | 1.3237 |
| gaussian | square | 4 | 1.2658 |
| ill_conditioned | short_wide | 1 | 1.4128 |
| ill_conditioned | short_wide | 4 | 1.4297 |
| ill_conditioned | square | 1 | 1.2153 |
| ill_conditioned | square | 4 | 1.2989 |
| orthonormal_rows | short_wide | 1 | 1.5156 |
| orthonormal_rows | short_wide | 4 | 1.5494 |
| orthonormal_rows | square | 1 | 1.499 |
| orthonormal_rows | square | 4 | 1.36 |

| family | regime | blas_threads | geomean relative time (bsqr_rinv/dgeqp3_trsm) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 1.3383 |
| gaussian | short_wide | 4 | 1.3362 |
| gaussian | square | 1 | 1.2968 |
| gaussian | square | 4 | 1.2705 |
| ill_conditioned | short_wide | 1 | 1.2146 |
| ill_conditioned | short_wide | 4 | 1.2437 |
| ill_conditioned | square | 1 | 1.2246 |
| ill_conditioned | square | 4 | 1.2461 |
| orthonormal_rows | short_wide | 1 | 1.3406 |
| orthonormal_rows | short_wide | 4 | 1.325 |
| orthonormal_rows | square | 1 | 1.2668 |
| orthonormal_rows | square | 4 | 1.2921 |
