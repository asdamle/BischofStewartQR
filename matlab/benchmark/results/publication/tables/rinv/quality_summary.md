# Numerical quality summary

Across all 150 benchmark cases (3 matrix families, both regimes, all sizes and seeds):

| method | median rel. residual | max rel. residual | median ||I - Q^T Q||_F | max ||I - Q^T Q||_F |
|---|---:|---:|---:|---:|
| BSQR + $R_{11}^{-1}R_{12}$ | 8.26e-16 | 1.61e-15 | 1.02e-14 | 3.88e-14 |
| CPQR + $R_{11}^{-1}R_{12}$ | 8.57e-16 | 1.63e-15 | 1.03e-14 | 3.88e-14 |

The relative residual ||A Pi - QR||_F / ||A||_F of BSQR + $R_{11}^{-1}R_{12}$ never exceeded 1.61e-15, and its deviation from orthogonality ||I - Q^T Q||_F never exceeded 3.88e-14; both match the built-in baseline to within a small factor.
