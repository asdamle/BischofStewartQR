# `matlab_general/` — BSQR column selection for general short-wide A

The Bischof–Stewart selection guarantee (deterministic `bsqr` in `matlab/` and
randomized `bsqr_rand` in `matlab_rand/`) is a statement about short-wide
matrices with **orthonormal rows**. This folder provides `bsqr_general`, which
transplants that guarantee to a **general** full-row-rank A (m×n, m ≤ n) via a
row-space reduction:

1. Reduced QR of the transpose: `A' = Qt*Rt` (Qt is n×m with orthonormal
   columns, Rt is m×m upper triangular).
2. Run the selector (`bsqr` or `bsqr_rand`) on `G = Qt'` — m×n with
   orthonormal rows, exactly the setting its theory covers.
3. Reassemble a genuine economy QR of the permuted A: since `A = Rt'*G`,
   `A(:,p) = (Rt'*Qg)*Rg`, and the small QR `Rt'*Qg = Qs*Rs` (O(m²k)) gives
   `A(:,p) = Qs*(Rs*Rg)` with `Rs*Rg` upper trapezoidal.

## The guarantee (k = m)

Selecting all m columns S = p(1:m), the chain

```
sigma_min(A(:,S)) >= sigma_min(Rt') * sigma_min(G(:,S))
                  =  sigma_min(A)   * sigma_min(G(:,S))
                  >= sigma_min(A) / ||R11(G)^{-1}||_F
                  >= sigma_min(A) / sqrt(m*(n-m+1))
```

holds for both selectors (the last step is the Bischof–Stewart / Osinsky-type
Frobenius bound each selector maintains on the orthonormal-row basis). So the
smallest singular value of the selected m columns is controlled **relative to
sigma_min(A)** — no direct pivoting method on A itself carries such a bound.

The guarantee is a k = m statement. Early stopping (`'k', k` with k < m) is
supported mechanically — the factorization contract below still holds — but
carries no comparable bound.

## Usage

```matlab
startup                          % repo root: paths + MEX backends
[Q, R, E] = bsqr_general(A);                             % deterministic selector
[Q, R, p] = bsqr_general(A, 'pivot_format', 'vector');   % index-vector pivot
[Q, R, p, stats] = bsqr_general(A, 'selector', 'bsqr_rand', 'seed', 1, ...
    'pivot_format', 'vector');                           % randomized selector
```

Output contract mirrors `bsqr` (for A m×n, k steps, default k = m):

| output   | meaning |
|----------|---------|
| `R`      | k×n upper-trapezoidal factor of the permuted columns; `A(:,p) = Q*R` for k = m. With early stop, `A(:,p(1:k)) = Q*R(:,1:k)` and `R(:,k+1:n) = Q'*A(:,p(k+1:n))` (bsqr's semantics). |
| `Q`      | m×k economy factor, orthonormal columns. |
| `E_OR_P` | n×n permutation matrix `E` with `A*E = A(:,p)` (default), or the 1×n index row `p` (`'pivot_format','vector'`). |
| `STATS`  | phase timings (`t_qr_At`, `t_select`, `t_assemble`, `t_total`), `rt_diag_ratio` (min/max `|diag(Rt)|`), and the selector's own stats struct for `bsqr_rand`. Replaces `bsqr`'s `R11INV_R12` slot (the initial `qr(A','econ')` already costs as much as that product would). |

Options owned by the wrapper: `'selector'` (`'bsqr'` default, `'bsqr_rand'`),
`'k'` (default m), `'pivot_format'` (`'matrix'` default, `'vector'`),
`'check_finite'` (default true; the wrapper scans once and forwards `false` to
the selector), `'rank_tol'` (default `max(m,n)*eps`; below this relative
`diag(Rt)` ratio a `bsqr_general:IllConditioned` warning is raised — the
factorization stays valid, but the guarantee degrades as sigma_min(A) → 0).

**Every other name-value pair is forwarded verbatim to the selected selector**
(`'backend'`, `'norm_recomp_tol'`; for `bsqr_rand` also `'seed'`, `'sampling'`,
`'block_size'`, `'batched'`, `'threshold_mode'`, `'slack'`, `'pick'`), and
unknown names are rejected by the selector's own parser. The selectors'
output-shape options (`'return_rinv_r12'`, `'return_r12'`, `'trace'`) are
blocked — the wrapper owns its output contract.

## Relationship to the sibling folders

Unlike `matlab_rand/` (which deliberately keeps its own copies of helpers and
never touches the deterministic kernels), this layer is a **thin wrapper by
design**: it calls `matlab/bsqr.m` and `matlab_rand/bsqr_rand.m` directly and
contains no kernel code of its own. It is pure MATLAB — the heavy work happens
inside the built-in `qr` and the selectors' MEX backends. It never modifies
either sibling.

## Tests

```bash
matlab -batch "addpath('matlab_general/tests'); run_general_tests"
```

`test_bsqr_general.m` covers the output contract (reconstruction, both pivot
formats, early stop, option routing/blocking, warnings); the guarantee itself
is asserted in `test_bsqr_general_bounds.m` at k = m across adversarial
families from `benchmark/general_test_matrix.m`, for both selectors.

## Benchmark

```bash
matlab -batch "addpath('matlab_general/benchmark'); run_general_comparison; plot_general_comparison"
```

Five methods per family/size/seed, each timed while materializing `Q`, `R`,
and the index-vector permutation: built-in `qr(A,'econ','vector')` (column
pivoting), `bsqr` and `bsqr_rand` applied directly to A (outside their
theory), and `bsqr_general` with each selector. Metrics: wall time (plus the
wrapper's qr/select phase split), `sigma_min(A(:,S))` and its ratio to
`sigma_min(A)` against the `1/sqrt(m(n-m+1))` guarantee, and
`||A(:,S)^{-1}||_F`. Families (see `general_test_matrix.m`): graded column
norms, a loud near-collinear cluster, graded spectra with needle/coherent
right factors, and a Kahan-block decoy — cases where norm-driven pivoting on
A is misled while the reduction is scale-invariant by construction.

CSVs land in the git-ignored `benchmark/results/`; figures in
`benchmark/plots/`. Expectation to keep in mind when reading the timing
figure: the wrapper pays an extra O(nm²) for `qr(A','econ')` up front, so it
will not beat direct BSQR on time — the story is guaranteed subset quality at
a bounded constant-factor cost.

There is one timing regime where the wrapper wins outright — against
**built-in column-pivoted QR** at large n:

```bash
matlab -batch "addpath('matlab_general/benchmark'); run_general_largen; plot_general_largen"
```

runs the two-method head-to-head (default `graded_cols`, m = k = 256, n up to
128k, `fig_general_largen`). LAPACK `dgeqp3`'s per-column pivot/norm-update
pass is BLAS-2-bound and jumps super-linearly once the trailing matrix falls
out of cache, while the wrapper's unpivoted `qr(A','econ')` is blocked BLAS-3
and its randomized selection stays a small fraction of the budget: on the
benchmark machine `bsqr_general(A,'selector','bsqr_rand')` overtakes
`qr(A,'econ','vector')` at n ≈ 9k and holds a 1.2–1.6× lead through n = 128k,
while staying an order of magnitude above the sigma_min guarantee.
