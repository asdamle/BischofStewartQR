# Publication Benchmark Summary

- Run ID: `20260310_213519`
- Generated: 2026-03-10T21:40:46.238
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1, 8
- Seeds: 20260310, 20260311

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.4482 |
| gaussian | short_wide | 8 | 0.44592 |
| gaussian | square | 1 | 0.48594 |
| gaussian | square | 8 | 0.48724 |
| ill_conditioned | short_wide | 1 | 0.48954 |
| ill_conditioned | short_wide | 8 | 0.49425 |
| ill_conditioned | square | 1 | 0.53084 |
| ill_conditioned | square | 8 | 0.51388 |
| orthonormal_rows | short_wide | 1 | 0.44864 |
| orthonormal_rows | short_wide | 8 | 0.4446 |
| orthonormal_rows | square | 1 | 0.48854 |
| orthonormal_rows | square | 8 | 0.49282 |
