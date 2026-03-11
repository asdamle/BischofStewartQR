# Publication Benchmark Summary

- Run ID: `20260310_235245`
- Generated: 2026-03-11T00:00:33.688
- BLAS: LBTConfig([ILP64] libopenblas64_.dylib, [LP64] Accelerate, [ILP64] Accelerate)
- Threads: 1, 4
- Seeds: 20260310, 20260311

| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr) |
|---|---|---:|---:|
| gaussian | short_wide | 1 | 0.49556 |
| gaussian | short_wide | 4 | 0.49747 |
| gaussian | square | 1 | 0.52359 |
| gaussian | square | 4 | 0.51731 |
| ill_conditioned | short_wide | 1 | 0.55004 |
| ill_conditioned | short_wide | 4 | 0.54669 |
| ill_conditioned | square | 1 | 0.56699 |
| ill_conditioned | square | 4 | 0.56422 |
| orthonormal_rows | short_wide | 1 | 0.4964 |
| orthonormal_rows | short_wide | 4 | 0.49445 |
| orthonormal_rows | square | 1 | 0.53 |
| orthonormal_rows | square | 4 | 0.51984 |
