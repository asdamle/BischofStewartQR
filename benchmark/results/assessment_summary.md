# Assessment Summary

- Generated: 2026-02-25T12:58:27.992
- Raw dataset: `benchmark/results/assessment_raw.csv`
- Pair count: 200

## Overall Speedup Distribution (`dgeqp3/bsqr`)
- mean: 0.57793
- median: 0.51975
- faster cases (>1): 15/200
- significant wins: 2
- significant losses: 171

## Case-Level Aggregates

| family | regime | tol | count | faster | slower | sig wins | sig losses | overlap | mean speedup | median speedup | mean lower-bound speedup |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gaussian | fixed_m_vary_n | 0.0 | 6 | 0 | 6 | 0 | 6 | 0 | 0.47701 | 0.46854 | 0.36664 |
| gaussian | fixed_m_vary_n | 1.0e-12 | 6 | 0 | 6 | 0 | 6 | 0 | 0.46714 | 0.46893 | 0.37665 |
| gaussian | fixed_m_vary_n | 1.0e-10 | 6 | 0 | 6 | 0 | 6 | 0 | 0.48136 | 0.48385 | 0.37697 |
| gaussian | fixed_m_vary_n | 1.4901161193847656e-8 | 6 | 0 | 6 | 0 | 5 | 1 | 0.49592 | 0.51965 | 0.36339 |
| gaussian | fixed_n_vary_m | 0.0 | 6 | 0 | 6 | 0 | 5 | 1 | 0.59482 | 0.54575 | 0.50984 |
| gaussian | fixed_n_vary_m | 1.0e-12 | 6 | 0 | 6 | 0 | 5 | 1 | 0.5985 | 0.5399 | 0.47328 |
| gaussian | fixed_n_vary_m | 1.0e-10 | 6 | 0 | 6 | 0 | 5 | 1 | 0.61189 | 0.54653 | 0.50217 |
| gaussian | fixed_n_vary_m | 1.4901161193847656e-8 | 6 | 0 | 6 | 0 | 5 | 1 | 0.59472 | 0.53599 | 0.51107 |
| gaussian | short_wide | 0.0 | 2 | 0 | 2 | 0 | 2 | 0 | 0.45791 | 0.45791 | 0.39571 |
| gaussian | short_wide | 1.0e-12 | 2 | 0 | 2 | 0 | 2 | 0 | 0.46346 | 0.46346 | 0.37132 |
| gaussian | short_wide | 1.0e-10 | 2 | 0 | 2 | 0 | 2 | 0 | 0.47052 | 0.47052 | 0.39462 |
| gaussian | short_wide | 1.4901161193847656e-8 | 2 | 0 | 2 | 0 | 2 | 0 | 0.46679 | 0.46679 | 0.39544 |
| gaussian | square | 0.0 | 2 | 0 | 2 | 0 | 2 | 0 | 0.60989 | 0.60989 | 0.47774 |
| gaussian | square | 1.0e-12 | 2 | 0 | 2 | 0 | 1 | 1 | 0.62077 | 0.62077 | 0.49884 |
| gaussian | square | 1.0e-10 | 2 | 0 | 2 | 0 | 2 | 0 | 0.60161 | 0.60161 | 0.47982 |
| gaussian | square | 1.4901161193847656e-8 | 2 | 0 | 2 | 0 | 1 | 1 | 0.58975 | 0.58975 | 0.47619 |
| gaussian | tall_skinny | 0.0 | 2 | 2 | 0 | 1 | 0 | 1 | 1.0532 | 1.0532 | 0.94108 |
| gaussian | tall_skinny | 1.0e-12 | 2 | 2 | 0 | 0 | 0 | 2 | 1.0595 | 1.0595 | 0.94651 |
| gaussian | tall_skinny | 1.0e-10 | 2 | 2 | 0 | 0 | 0 | 2 | 1.0428 | 1.0428 | 0.93067 |
| gaussian | tall_skinny | 1.4901161193847656e-8 | 2 | 1 | 1 | 0 | 0 | 2 | 1.0434 | 1.0434 | 0.91507 |
| ill_conditioned | fixed_m_vary_n | 0.0 | 6 | 0 | 6 | 0 | 6 | 0 | 0.50726 | 0.49859 | 0.4313 |
| ill_conditioned | fixed_m_vary_n | 1.0e-12 | 6 | 0 | 6 | 0 | 6 | 0 | 0.50437 | 0.49244 | 0.42683 |
| ill_conditioned | fixed_m_vary_n | 1.0e-10 | 6 | 0 | 6 | 0 | 6 | 0 | 0.49837 | 0.49856 | 0.37455 |
| ill_conditioned | fixed_m_vary_n | 1.4901161193847656e-8 | 6 | 0 | 6 | 0 | 6 | 0 | 0.48902 | 0.48342 | 0.38427 |
| ill_conditioned | fixed_n_vary_m | 0.0 | 6 | 0 | 6 | 0 | 5 | 1 | 0.64287 | 0.57439 | 0.55118 |
| ill_conditioned | fixed_n_vary_m | 1.0e-12 | 6 | 0 | 6 | 0 | 5 | 1 | 0.61966 | 0.55163 | 0.51124 |
| ill_conditioned | fixed_n_vary_m | 1.0e-10 | 6 | 0 | 6 | 0 | 6 | 0 | 0.62408 | 0.55925 | 0.52715 |
| ill_conditioned | fixed_n_vary_m | 1.4901161193847656e-8 | 6 | 0 | 6 | 0 | 5 | 1 | 0.62297 | 0.56557 | 0.53085 |
| ill_conditioned | short_wide | 0.0 | 2 | 0 | 2 | 0 | 2 | 0 | 0.48185 | 0.48185 | 0.39107 |
| ill_conditioned | short_wide | 1.0e-12 | 2 | 0 | 2 | 0 | 2 | 0 | 0.49306 | 0.49306 | 0.39115 |
| ill_conditioned | short_wide | 1.0e-10 | 2 | 0 | 2 | 0 | 2 | 0 | 0.47456 | 0.47456 | 0.36253 |
| ill_conditioned | short_wide | 1.4901161193847656e-8 | 2 | 0 | 2 | 0 | 2 | 0 | 0.46614 | 0.46614 | 0.38245 |
| ill_conditioned | square | 0.0 | 2 | 0 | 2 | 0 | 2 | 0 | 0.59132 | 0.59132 | 0.48611 |
| ill_conditioned | square | 1.0e-12 | 2 | 0 | 2 | 0 | 1 | 1 | 0.63274 | 0.63274 | 0.49184 |
| ill_conditioned | square | 1.0e-10 | 2 | 0 | 2 | 0 | 2 | 0 | 0.59801 | 0.59801 | 0.46407 |
| ill_conditioned | square | 1.4901161193847656e-8 | 2 | 0 | 2 | 0 | 2 | 0 | 0.58281 | 0.58281 | 0.46652 |
| ill_conditioned | tall_skinny | 0.0 | 2 | 2 | 0 | 1 | 0 | 1 | 1.0809 | 1.0809 | 0.91166 |
| ill_conditioned | tall_skinny | 1.0e-12 | 2 | 2 | 0 | 0 | 0 | 2 | 1.0615 | 1.0615 | 0.76268 |
| ill_conditioned | tall_skinny | 1.0e-10 | 2 | 2 | 0 | 0 | 0 | 2 | 1.0494 | 1.0494 | 0.87188 |
| ill_conditioned | tall_skinny | 1.4901161193847656e-8 | 2 | 2 | 0 | 0 | 0 | 2 | 1.0584 | 1.0584 | 0.91755 |
| orthonormal_rows | fixed_m_vary_n | 0.0 | 6 | 0 | 6 | 0 | 6 | 0 | 0.47809 | 0.46539 | 0.36623 |
| orthonormal_rows | fixed_m_vary_n | 1.0e-12 | 6 | 0 | 6 | 0 | 6 | 0 | 0.47162 | 0.46422 | 0.38477 |
| orthonormal_rows | fixed_m_vary_n | 1.0e-10 | 6 | 0 | 6 | 0 | 5 | 1 | 0.48964 | 0.47792 | 0.39076 |
| orthonormal_rows | fixed_m_vary_n | 1.4901161193847656e-8 | 6 | 0 | 6 | 0 | 6 | 0 | 0.47975 | 0.46991 | 0.37826 |
| orthonormal_rows | fixed_n_vary_m | 0.0 | 4 | 0 | 4 | 0 | 4 | 0 | 0.5193 | 0.51055 | 0.42609 |
| orthonormal_rows | fixed_n_vary_m | 1.0e-12 | 4 | 0 | 4 | 0 | 4 | 0 | 0.53055 | 0.52052 | 0.44847 |
| orthonormal_rows | fixed_n_vary_m | 1.0e-10 | 4 | 0 | 4 | 0 | 4 | 0 | 0.52843 | 0.52247 | 0.4418 |
| orthonormal_rows | fixed_n_vary_m | 1.4901161193847656e-8 | 4 | 0 | 4 | 0 | 4 | 0 | 0.53602 | 0.52383 | 0.42393 |
| orthonormal_rows | short_wide | 0.0 | 2 | 0 | 2 | 0 | 2 | 0 | 0.45971 | 0.45971 | 0.39507 |
| orthonormal_rows | short_wide | 1.0e-12 | 2 | 0 | 2 | 0 | 2 | 0 | 0.46479 | 0.46479 | 0.40169 |
| orthonormal_rows | short_wide | 1.0e-10 | 2 | 0 | 2 | 0 | 2 | 0 | 0.45848 | 0.45848 | 0.38364 |
| orthonormal_rows | short_wide | 1.4901161193847656e-8 | 2 | 0 | 2 | 0 | 2 | 0 | 0.46799 | 0.46799 | 0.39812 |
| orthonormal_rows | square | 0.0 | 2 | 0 | 2 | 0 | 2 | 0 | 0.58819 | 0.58819 | 0.47306 |
| orthonormal_rows | square | 1.0e-12 | 2 | 0 | 2 | 0 | 2 | 0 | 0.60608 | 0.60608 | 0.48825 |
| orthonormal_rows | square | 1.0e-10 | 2 | 0 | 2 | 0 | 2 | 0 | 0.59036 | 0.59036 | 0.47636 |
| orthonormal_rows | square | 1.4901161193847656e-8 | 2 | 0 | 2 | 0 | 1 | 1 | 0.6316 | 0.6316 | 0.48845 |

## Quality-Risk Map by Tolerance

| tol | count | median residual | max residual | median orthogonality | max orthogonality | nonfinite |
|---:|---:|---:|---:|---:|---:|---:|
| 0.0 | 50 | 5.866e-16 | 8.7901e-16 | 6.1305e-15 | 3.3058e-14 | 0 |
| 1.0e-12 | 50 | 5.8772e-16 | 8.7901e-16 | 6.156e-15 | 3.3914e-14 | 0 |
| 1.0e-10 | 50 | 5.8772e-16 | 8.7901e-16 | 6.156e-15 | 3.3914e-14 | 0 |
| 1.4901161193847656e-8 | 50 | 5.8772e-16 | 8.7901e-16 | 6.156e-15 | 3.3914e-14 | 0 |

## Hotspot Correlation
- mean phase shares: apply=33.129%, w_update=39.583%, downdate=21.894%, pivot=1.1618%
- corr(slowdown, apply_pct)=-0.84789
- corr(slowdown, w_update_pct)=0.91979
- corr(slowdown, downdate_pct)=0.77305
- corr(slowdown, pivot_pct)=0.73421
- weighted slowdown contributions: apply=0.22273, w_update=0.37408, downdate=0.20944, pivot=0.011042

## Top 10 Best Cases
- gaussian tall_skinny m=512 n=64 k=64 tol=1.4901161193847656e-8: speedup=1.1008 significance=overlap
- ill_conditioned tall_skinny m=512 n=64 k=64 tol=1.0e-12: speedup=1.0997 significance=overlap
- ill_conditioned tall_skinny m=512 n=64 k=64 tol=1.4901161193847656e-8: speedup=1.0992 significance=overlap
- ill_conditioned tall_skinny m=512 n=64 k=64 tol=1.0e-10: speedup=1.0971 significance=overlap
- ill_conditioned tall_skinny m=512 n=64 k=64 tol=0.0: speedup=1.0931 significance=sig_win
- gaussian tall_skinny m=512 n=64 k=64 tol=0.0: speedup=1.092 significance=sig_win
- gaussian tall_skinny m=512 n=64 k=64 tol=1.0e-12: speedup=1.0806 significance=overlap
- ill_conditioned tall_skinny m=1024 n=64 k=64 tol=0.0: speedup=1.0687 significance=overlap
- gaussian tall_skinny m=512 n=64 k=64 tol=1.0e-10: speedup=1.0495 significance=overlap
- gaussian tall_skinny m=1024 n=64 k=64 tol=1.0e-12: speedup=1.0385 significance=overlap

## Top 10 Worst Cases
- gaussian fixed_m_vary_n m=32 n=512 k=32 tol=1.0e-10: speedup=0.38406 significance=sig_loss
- gaussian fixed_m_vary_n m=64 n=2048 k=64 tol=1.4901161193847656e-8: speedup=0.38931 significance=sig_loss
- gaussian fixed_m_vary_n m=32 n=512 k=32 tol=1.0e-12: speedup=0.40985 significance=sig_loss
- orthonormal_rows fixed_m_vary_n m=32 n=512 k=32 tol=1.0e-10: speedup=0.41104 significance=sig_loss
- orthonormal_rows fixed_m_vary_n m=32 n=2048 k=32 tol=0.0: speedup=0.42984 significance=sig_loss
- orthonormal_rows fixed_m_vary_n m=32 n=2048 k=32 tol=1.0e-10: speedup=0.43112 significance=sig_loss
- orthonormal_rows short_wide m=64 n=1024 k=64 tol=0.0: speedup=0.43299 significance=sig_loss
- orthonormal_rows fixed_m_vary_n m=32 n=128 k=32 tol=0.0: speedup=0.43467 significance=sig_loss
- orthonormal_rows short_wide m=64 n=1024 k=64 tol=1.0e-12: speedup=0.4349 significance=sig_loss
- orthonormal_rows fixed_m_vary_n m=32 n=512 k=32 tol=1.0e-12: speedup=0.43528 significance=sig_loss
