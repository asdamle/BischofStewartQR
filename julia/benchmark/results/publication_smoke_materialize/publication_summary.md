# Publication Benchmark Summary

- Run ID: `20260316_112940`
- Generated: 2026-03-16T11:50:53.255
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Bench Surface: materialize_qrp
- Threads: 1, 4
- Seeds: 20260310, 20260311

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr_full) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.55152 |
| gaussian | short_wide | 4 | 0.51837 |
| gaussian | square | 1 | 0.62784 |
| gaussian | square | 4 | 0.63964 |
| ill_conditioned | short_wide | 1 | 0.58436 |
| ill_conditioned | short_wide | 4 | 0.57961 |
| ill_conditioned | square | 1 | 0.62432 |
| ill_conditioned | square | 4 | 0.61535 |
| orthonormal_rows | short_wide | 1 | 0.5515 |
| orthonormal_rows | short_wide | 4 | 0.54234 |
| orthonormal_rows | square | 1 | 0.63936 |
| orthonormal_rows | square | 4 | 0.63374 |

| family | regime | blas_threads | geomean speedup (dgeqp3_trsm/bsqr_rinv) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.64951 |
| gaussian | short_wide | 4 | 0.61657 |
| gaussian | square | 1 | 0.60772 |
| gaussian | square | 4 | 0.60603 |
| ill_conditioned | short_wide | 1 | 0.69737 |
| ill_conditioned | short_wide | 4 | 0.68308 |
| ill_conditioned | square | 1 | 0.64545 |
| ill_conditioned | square | 4 | 0.65435 |
| orthonormal_rows | short_wide | 1 | 0.64719 |
| orthonormal_rows | short_wide | 4 | 0.64102 |
| orthonormal_rows | square | 1 | 0.63471 |
| orthonormal_rows | square | 4 | 0.62134 |
