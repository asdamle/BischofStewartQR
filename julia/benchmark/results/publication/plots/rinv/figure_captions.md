# Suggested captions (rinv comparison)

Run: 20260611_204425; seeds: 2; BLAS threads: 1, 4. In all figures, BSQR (W returned) is compared against CPQR + solve; both timed paths materialize Q, R, and the permutation.

**Test matrices.** Three families, regenerated per seed: *Gaussian* — i.i.d. standard normal entries; *ill-conditioned* — A = UΣVᵀ with U, V orthonormal factors of Gaussian matrices and Σ geometrically graded from 1 down to 1e-10 (κ = 1e10); *orthonormal rows* — A = Qᵀ with Q an orthonormal basis of a Gaussian n×m matrix (m ≤ n, so AAᵀ = I — the GKS column-selection setting). Square cases use m = n; short-wide cases sweep aspect ratios n/m ∈ {2, 4, 8, 10}.

**fig_square_runtime.** Median runtime versus matrix size for square matrices (log–log). Lines are geometric means over 2 seeds; faint bands show the range across seeds. Panels: matrix family × BLAS threads.

**fig_shortwide_runtime.** Median runtime versus column count n for short-wide matrices, one color per row count m (log–log). Lines/bands as above; marker and line style distinguish the methods.

**fig_relative_time.** Geometric-mean relative time (BSQR+W / CPQR+solve; 1 = parity, dashed line) per family and regime. Points are geomeans of per-seed geomeans; whiskers show the per-seed range.

Numerical quality is summarized in tables/quality_summary.md.
