# Suggested captions (rinv comparison)

Run: 20260611_210948; seeds: 2; MATLAB default threading. In all figures, BSQR + $R_{11}^{-1}R_{12}$ is compared against CPQR + $R_{11}^{-1}R_{12}$; both timed paths materialize Q, R, and the permutation.

**Test matrices.** Three families, regenerated per seed: *Gaussian* - i.i.d. standard normal entries; *ill-conditioned* - A = U S V^T with U, V orthonormal factors of Gaussian matrices and S geometrically graded from 1 down to 1e-10 (kappa = 1e10); *orthonormal rows* - A = Q^T with Q an orthonormal basis of a Gaussian n-by-m matrix (m <= n, so A A^T = I - the GKS column-selection setting). Square cases use m = n; short-wide cases sweep aspect ratios n/m in {2, 4, 8, 10}.

**fig_square_runtime.** Median runtime versus matrix size for square matrices (log-log). Lines are geometric means over 2 seeds; faint bands show the range across seeds. One panel per matrix family.

**fig_shortwide_runtime.** Median runtime versus column count n for short-wide matrices, one color per row count m (log-log). Lines/bands as above; marker and line style distinguish the methods.

**fig_relative_time.** Geometric-mean relative time ((BSQR + $R_{11}^{-1}R_{12}$) / (CPQR + $R_{11}^{-1}R_{12}$); 1 = parity, dashed line) per family and regime. Points are geomeans of per-seed geomeans; whiskers show the per-seed range.

**fig_relative_time_composite** (top-level plots/). The single timing figure: the relative time without (BSQR / CPQR) and with the interpolation matrix (both methods also form R11^{-1}R12; labelled (BSQR + R11^{-1}R12) / (CPQR + R11^{-1}R12)), overlaid on the shared rows and distinguished by colour. Square rows are m = n, short-wide rows m < n. Points/whiskers as for fig_relative_time.

Numerical quality is summarized in tables/quality_summary.md.
