# Publication Benchmark Summary

- Run ID: `20260310_225439`
- Generated: 2026-03-10T23:02:40.570
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1, 4
- Seeds: 20260310, 20260311

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.48397 |
| gaussian | short_wide | 4 | 0.47454 |
| gaussian | square | 1 | 0.49914 |
| gaussian | square | 4 | 0.47625 |
| ill_conditioned | short_wide | 1 | 0.53258 |
| ill_conditioned | short_wide | 4 | 0.51939 |
| ill_conditioned | square | 1 | 0.5413 |
| ill_conditioned | square | 4 | 0.50614 |
| orthonormal_rows | short_wide | 1 | 0.48903 |
| orthonormal_rows | short_wide | 4 | 0.47115 |
| orthonormal_rows | square | 1 | 0.49613 |
| orthonormal_rows | square | 4 | 0.48595 |
