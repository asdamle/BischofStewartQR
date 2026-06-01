# MATLAB Publication Benchmark Summary

- Run ID: `20260316_153811`
- Generated: 16-Mar-2026 15:44:11
- Seeds: 20260310, 20260311
- Families: gaussian, ill_conditioned, orthonormal_rows

Relative time is BSQR median time divided by baseline median time; `1.0` is parity.

| family | regime | geomean relative time (bsqr_full/qr) |
|---|---|---:|
| gaussian | short_wide | 2.1926 |
| gaussian | square | 1.6547 |
| ill_conditioned | short_wide | 1.8998 |
| ill_conditioned | square | 1.7577 |
| orthonormal_rows | short_wide | 2.0766 |
| orthonormal_rows | square | 1.8892 |

| family | regime | geomean relative time (bsqr_rinv/qr_trsm) |
|---|---|---:|
| gaussian | short_wide | 1.7483 |
| gaussian | square | 1.6393 |
| ill_conditioned | short_wide | 1.6902 |
| ill_conditioned | square | 1.6778 |
| orthonormal_rows | short_wide | 1.7811 |
| orthonormal_rows | square | 1.7719 |
