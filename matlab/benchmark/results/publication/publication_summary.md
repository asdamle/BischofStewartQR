# MATLAB Publication Benchmark Summary

- Run ID: `20260316_125458`
- Generated: 16-Mar-2026 13:01:23
- Seeds: 20260310, 20260311
- Families: gaussian, ill_conditioned, orthonormal_rows

- Bench Surface: materialize_qrp

| family | regime | geomean speedup (qr/bsqr_full) |
|---|---|---:|
| gaussian | short_wide | 0.45199 |
| gaussian | square | 0.57128 |
| ill_conditioned | short_wide | 0.49469 |
| ill_conditioned | square | 0.56083 |
| orthonormal_rows | short_wide | 0.44606 |
| orthonormal_rows | square | 0.55647 |

| family | regime | geomean speedup (qr_trsm/bsqr_rinv) |
|---|---|---:|
| gaussian | short_wide | 0.52883 |
| gaussian | square | 0.58589 |
| ill_conditioned | short_wide | 0.57731 |
| ill_conditioned | square | 0.58751 |
| orthonormal_rows | short_wide | 0.51647 |
| orthonormal_rows | square | 0.51963 |
