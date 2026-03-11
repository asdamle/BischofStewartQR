# Publication Benchmark Summary

- Run ID: `20260310_215816`
- Generated: 2026-03-10T22:03:36.442
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1, 4
- Seeds: 20260310, 20260311

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.44524 |
| gaussian | short_wide | 4 | 0.45573 |
| gaussian | square | 1 | 0.47511 |
| gaussian | square | 4 | 0.4853 |
| ill_conditioned | short_wide | 1 | 0.50319 |
| ill_conditioned | short_wide | 4 | 0.47308 |
| ill_conditioned | square | 1 | 0.52385 |
| ill_conditioned | square | 4 | 0.51934 |
| orthonormal_rows | short_wide | 1 | 0.44275 |
| orthonormal_rows | short_wide | 4 | 0.45531 |
| orthonormal_rows | square | 1 | 0.49064 |
| orthonormal_rows | square | 4 | 0.53373 |
