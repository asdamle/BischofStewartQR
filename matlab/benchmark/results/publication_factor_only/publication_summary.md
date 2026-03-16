# MATLAB Publication Benchmark Summary

- Run ID: `20260316_130208`
- Generated: 16-Mar-2026 13:06:55
- Seeds: 20260310, 20260311
- Families: gaussian, ill_conditioned, orthonormal_rows

- Bench Surface: factor_only

| family | regime | geomean speedup (qr/bsqr_full) |
|---|---|---:|
| gaussian | short_wide | 0.47684 |
| gaussian | square | 0.73597 |
| ill_conditioned | short_wide | 0.51857 |
| ill_conditioned | square | 0.77642 |
| orthonormal_rows | short_wide | 0.46763 |
| orthonormal_rows | square | 0.71015 |

| family | regime | geomean speedup (qr_trsm/bsqr_rinv) |
|---|---|---:|
| gaussian | short_wide | 0.59077 |
| gaussian | square | 0.75215 |
| ill_conditioned | short_wide | 0.59626 |
| ill_conditioned | square | 0.77925 |
| orthonormal_rows | short_wide | 0.54024 |
| orthonormal_rows | square | 0.71333 |
