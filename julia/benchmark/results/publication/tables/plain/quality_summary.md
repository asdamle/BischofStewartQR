# Numerical quality summary

Across all 750 benchmark cases (3 matrix families, both regimes, all sizes, seeds, and thread settings):

| method | median rel. residual | max rel. residual | median ‖I − QᵀQ‖_F | max ‖I − QᵀQ‖_F |
|---|---:|---:|---:|---:|
| BSQR | 8.36e-16 | 1.61e-15 | 1.02e-14 | 3.89e-14 |
| CPQR (built-in) | 8.83e-16 | 1.64e-15 | 1.12e-14 | 3.90e-14 |

The relative residual ‖AΠ − QR‖_F/‖A‖_F of BSQR never exceeded 1.61e-15, and its deviation from orthogonality ‖I − QᵀQ‖_F never exceeded 3.89e-14; both match the built-in baseline to within a small factor.
