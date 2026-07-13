# MATLAB Publication Benchmark Summary

- Run ID: `20260712_215106`
- Generated: 12-Jul-2026 21:54:13
- Seeds: 20260310, 20260311
- Families: gaussian, ill_conditioned, orthonormal_rows

Relative time is BSQR median time divided by baseline median time; `1.0` is parity.

| family | regime | geomean relative time (bsqr_full/qr) |
|---|---|---:|
| gaussian | short_wide | 1.68965 |
| gaussian | square | 1.30674 |
| ill_conditioned | short_wide | 1.56493 |
| ill_conditioned | square | 1.24993 |
| orthonormal_rows | short_wide | 1.71833 |
| orthonormal_rows | square | 1.31795 |

| family | regime | geomean relative time (bsqr_rinv/qr_trsm) |
|---|---|---:|
| gaussian | short_wide | 1.44577 |
| gaussian | square | 1.24702 |
| ill_conditioned | short_wide | 1.35135 |
| ill_conditioned | square | 1.26257 |
| orthonormal_rows | short_wide | 1.47844 |
| orthonormal_rows | square | 1.32600 |
