# Publication Benchmark Summary

- Run ID: `20260712_215448`
- Generated: 2026-07-12T22:12:26.631
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1, 4
- Seeds: 20260310, 20260311


Relative time is BSQR median time divided by baseline median time; `1.0` is parity.

| family | regime | blas_threads | geomean relative time (bsqr_full/dgeqp3) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 1.5408 |
| gaussian | short_wide | 4 | 1.5314 |
| gaussian | square | 1 | 1.2692 |
| gaussian | square | 4 | 1.2656 |
| ill_conditioned | short_wide | 1 | 1.4045 |
| ill_conditioned | short_wide | 4 | 1.4124 |
| ill_conditioned | square | 1 | 1.2054 |
| ill_conditioned | square | 4 | 1.2073 |
| orthonormal_rows | short_wide | 1 | 1.5337 |
| orthonormal_rows | short_wide | 4 | 1.5341 |
| orthonormal_rows | square | 1 | 1.2545 |
| orthonormal_rows | square | 4 | 1.2681 |

| family | regime | blas_threads | geomean relative time (bsqr_rinv/dgeqp3_trsm) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 1.3151 |
| gaussian | short_wide | 4 | 1.3049 |
| gaussian | square | 1 | 1.2462 |
| gaussian | square | 4 | 1.2688 |
| ill_conditioned | short_wide | 1 | 1.2358 |
| ill_conditioned | short_wide | 4 | 1.2303 |
| ill_conditioned | square | 1 | 1.2087 |
| ill_conditioned | square | 4 | 1.1977 |
| orthonormal_rows | short_wide | 1 | 1.3236 |
| orthonormal_rows | short_wide | 4 | 1.3182 |
| orthonormal_rows | square | 1 | 1.2726 |
| orthonormal_rows | square | 4 | 1.2611 |
