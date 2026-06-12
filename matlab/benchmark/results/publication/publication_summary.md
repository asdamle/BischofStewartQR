# MATLAB Publication Benchmark Summary

- Run ID: `20260611_210948`
- Generated: 11-Jun-2026 21:13:19
- Seeds: 20260310, 20260311
- Families: gaussian, ill_conditioned, orthonormal_rows

Relative time is BSQR median time divided by baseline median time; `1.0` is parity.

| family | regime | geomean relative time (bsqr_full/qr) |
|---|---|---:|
| gaussian | short_wide | 1.67439 |
| gaussian | square | 1.30416 |
| ill_conditioned | short_wide | 1.52907 |
| ill_conditioned | square | 1.26500 |
| orthonormal_rows | short_wide | 1.67220 |
| orthonormal_rows | square | 1.34112 |

| family | regime | geomean relative time (bsqr_rinv/qr_trsm) |
|---|---|---:|
| gaussian | short_wide | 1.44983 |
| gaussian | square | 1.30896 |
| ill_conditioned | short_wide | 1.31115 |
| ill_conditioned | square | 1.24851 |
| orthonormal_rows | short_wide | 1.41874 |
| orthonormal_rows | square | 1.31048 |
