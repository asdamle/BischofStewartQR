# MATLAB Publication Benchmark Summary

- Run ID: `20260316_143619`
- Generated: 16-Mar-2026 14:40:34
- Seeds: 20260310, 20260311
- Families: gaussian, ill_conditioned, orthonormal_rows

- Bench Surface: materialize_qrp

| family | regime | geomean speedup (qr/bsqr_full) |
|---|---|---:|
| gaussian | short_wide | 0.47707 |
| gaussian | square | 0.60206 |
| ill_conditioned | short_wide | 0.54572 |
| ill_conditioned | square | 0.62556 |
| orthonormal_rows | short_wide | 0.49221 |
| orthonormal_rows | square | 0.58464 |

| family | regime | geomean speedup (qr_trsm/bsqr_rinv) |
|---|---|---:|
| gaussian | short_wide | 0.58358 |
| gaussian | square | 0.61804 |
| ill_conditioned | short_wide | 0.63057 |
| ill_conditioned | square | 0.62151 |
| orthonormal_rows | short_wide | 0.56468 |
| orthonormal_rows | square | 0.62716 |
