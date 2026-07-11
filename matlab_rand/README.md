# matlab_rand — randomized Bischof–Stewart column selection (experimental)

A randomized variant of BSQR that selects `k` columns of a `k×n` matrix with the
same theoretical guarantees on `||R11^{-1}||_F` as the deterministic kernel, but
**without** maintaining `R11^{-1}R12` / column norms for every column at every
step. Instead it tracks only the running squared inverse Frobenius norm and
samples candidate columns in blocks, keeping those that hold the running value
under the per-step bound. By default it runs **in-block** (`batched`): each
sampled block is brought to the current frame once, then BSQR is run within it to
take as many columns as the bound allows before resampling -- amortizing the
per-block reflector apply over many selections (`O(k^3)` overall vs the
single-select `O(k^4)`). See `docs/RANDOMIZED_BSQR_PLAN.md` for the math and §5
of the manuscript (`notes/GKSevolved_draft.tex`, local-only until final;
Thm. 5.1 `thm:randBSpivot`) for the guarantee it relies on.

This is **separate from and does not modify** the deterministic implementations
in `matlab/` and `julia/`.

## Quick start

```matlab
addpath('matlab_rand'); addpath('matlab_rand/mex');   % mex auto-added by bsqr_rand

M = orth(randn(2000, 32))';            % 32×2000, orthonormal rows (GKS setting)
[p, Q, R11] = bsqr_rand(M);            % default product: subset + economy Q + R11
sel = p(1:32);                         % selected column indices
```

`Q` is the `m×k` economy factor with `Q'*M(:,p(1:k)) = R11`, formed lazily from
the kernel's accumulated reflectors (LAPACK `dorgqr` in the MEX) only when the
output is requested — `p = bsqr_rand(M)` never pays for it.

`R12` is opt-in (it costs an extra `O(n k^2)` pass and is off by default):

```matlab
[p, Q, R11, stats, R12] = bsqr_rand(M, 'return_r12', true);
```

## Build the MEX backend

```matlab
matlab -batch "addpath('matlab_rand'); build_bsqr_rand_mex"
```

`backend='auto'` (default) uses `bsqr_rand_mex` when built, else the pure-`.m`
reference `private/bsqr_rand_mfile.m`. The m-file is the correctness reference;
benchmarks use the MEX.

## Tests

```matlab
matlab -batch "addpath('matlab_rand'); run('matlab_rand/tests/run_rand_tests.m')"
```

`test_bsqr_rand.m` checks the factorization is exact; `test_bsqr_rand_bounds.m`
checks the headline guarantee (`||R11^{-1}||_F` stays under `sqrt(k(n-k+1))`).
The MEX uses its own RNG, so its pivots need not match the m-file — only the
guarantees must hold.

## Benchmarks

```matlab
matlab -batch "addpath('matlab_rand'); addpath('matlab'); addpath('matlab_rand/benchmark'); run_rand_benchmarks"
```

Compares the randomized selection (`[p, Q, R11]`, no `R12`) against the
deterministic factor path and reports speedup, the conditioning ratio
`||R11^{-1}||_F / sqrt(k(n-k+1))`, and candidate columns tested per pivot. Writes
`benchmark/results/rand_timings.csv`.

The full publication study (scaling, block size, sampling — with seed
min/max bands) is generated per `k` and written to `k`-tagged files:

```matlab
% one k value
run_rand_experiments('k', 128); plot_rand_experiments('k', 128);
% the committed sweep
for K = [64 128 256]; run_rand_experiments('k', K); plot_rand_experiments('k', K); end
```

This writes `benchmark/results/exp_*_k<K>.csv` and `benchmark/plots/fig_*_k<K>.png`.

## External comparison: `rejection_rpqr` (Adaptive Randomized Pivoting)

A standalone comparison against the `rejection_rpqr` selector from Epperly et
al.'s Adaptive Randomized Pivoting. The third-party code is **not** redistributed
here (`ext_comparisons/` is git-ignored); download it yourself:

```bash
# from the repo root
mkdir -p ext_comparisons && cd ext_comparisons
curl -L https://github.com/eepperly/Adaptive-Randomized-Pivoting/archive/refs/heads/main.tar.gz | tar xz
# -> ext_comparisons/Adaptive-Randomized-Pivoting-main/
```

It ships a compiled `rejection_helper` MEX for Apple Silicon; on other platforms
run its `code/compile_script.m`. Then:

```matlab
matlab -batch "addpath('matlab_rand/benchmark'); run_rpqr_comparison('k', 64); plot_rpqr_comparison('k', 64)"
```

This writes `exp_rpqr_k<K>.csv` and the three comparison figures
`fig_rpqr_{time,speedup,quality}_k<K>.png` (same metrics as the main suite). The
comparison is run fairly:

- **BSQR uses norm-weighted sampling**, matching `rejection_rpqr`'s
  squared-column-norm (leverage) sampling — the apples-to-apples choice.
- `rejection_rpqr` is the published `.m` implementation; `bsqr_rand` is our MEX,
  so timing reflects implementation language as well as algorithm.
- The two optimize **different objectives**: BSQR directly minimizes the growth of
  `||R11^{-1}||_F`, while ARP targets a volume/DPP criterion — so the
  `||R11^{-1}||_F` metric measures BSQR's objective and favors it by construction.

### Large-n scaling (sampler comparison)

A fixed-`k` (=64), large-`n` runtime study isolating the **sampling scheme**:
randomized BSQR with norm-weighted vs uniform sampling, and `rejection_rpqr`, on
two leverage regimes — `gaussian` (benign; uniform is fine) and `needle`
(~`k` high-leverage columns among many near-null ones; uniform keeps missing them
and must resample ~`n/k` times to find each).

```matlab
matlab -batch "addpath('matlab_rand'); addpath('matlab_rand/benchmark'); run_largen_scaling; plot_largen_scaling"
```

Writes `exp_largen.csv` and `fig_largen_scaling.{png,pdf}` (one log-log time-vs-`n`
panel per family; line = seed median, band = seed min/max; 20 seeds, `timeit`,
swept to `n = 1e6`). Same fair-timing discipline as `run_rpqr_comparison` (direct
`bsqr_rand_mex`, `check_finite=false`). The takeaway: norm-weighted BSQR and
`rejection_rpqr` both carry the common `O(mn)` norm work and scale alike (BSQR
~constant-factor faster, **no crossover** even at `n = 1e6`); uniform BSQR is
fastest on `gaussian` (it skips the norm precompute) but a constant factor (~20×)
slower on `needle` — it resamples to hit the rare high-leverage columns, so it
stays linear in `n` but well above the two norm-weighted methods.

### Matrix-approximation comparison

A second, application-facing comparison that scores the **low-rank approximation
error** of each method's column subset (the metric a practitioner cares about),
in both the Frobenius and spectral norms — a fairer counterpart to the
`||R11^{-1}||_F` comparison above, which is BSQR's own objective. It follows the
standard CSSP-via-leading-singular-vectors pipeline (exactly how ARP's `arp.m`
uses `rejection_rpqr`):

1. Build a full application matrix `A` (`approx_test_matrix` for synthetic
   families, `approx_real_matrix` for real data).
2. Compute accurate leading **right** singular vectors with `svds` (Lanczos):
   `W = V(:,1:k)'` is `k×n` with orthonormal rows. **Both selectors get the same
   `W`**, so the comparison isolates column selection, not subspace estimation —
   which is shared and accurate, and explicitly not the point. (ARP normally feeds
   a *sketched* `Q'`; exact `V_k'` to both is even more apples-to-apples.)
3. Each method picks `k` of the `n` columns; the subset is scored by the
   interpolation-free orthogonal-projection error `||A − P_S A||` (identical
   scoring for both). Reference: the optimal rank-`k` error.

The matrix is fixed per family (matching the ARP accuracy study's design); trials
vary only the selectors' RNG. `bsqr_rand` runs with its public defaults (batched,
norm-weighted). Run and plot:

```matlab
matlab -batch "addpath('matlab_rand'); addpath('matlab_rand/benchmark'); run_approx_comparison; plot_approx_comparison"
```

Synthetic families (always available, fully reproducible): `gmm_kernel` (RBF
kernel on a Gaussian-mixture cloud — Nyström landmark selection),
`integral_skeleton` (`1/dist` kernel between separated clouds — skeletonization /
interpolative decomposition), `snapshots` (parametric model-reduction snapshots —
empirical interpolation / DEIM). This writes `results/exp_approx.csv` and
`plots/fig_approx_quality.{png,pdf}`.

**Synthetic-spectrum companion.** To see how much the `||R11^{-1}||_F` differences
from the comparison above translate into *approximation* error, a companion run
uses matrices `A = U diag(s) V'` with a prescribed, deliberately interesting
spectrum `s` (a few large values, decay, a flatter section, more decay) and right
singular vectors `V` taken from the **same leverage families as the rpqr
comparison** (`gaussian`, `spiked_leverage`, `needle` — see `approx_synth_matrix`).
The three families share one spectrum, so the optimal-error curve is identical
across panels and only the leverage structure varies:

```matlab
matlab -batch "addpath('matlab_rand'); addpath('matlab_rand/benchmark'); run_approx_synth_comparison; plot_approx_comparison('tag','_synth')"
```

This writes `results/exp_approx_synth.csv` and `plots/fig_approx_synth_quality.{png,pdf}`.
The figure pairs the accuracy rows (Frobenius / spectral) with an
interpolation-coefficient row `max|R11^{-1}R12|` on the *same* matrices, so it shows
directly that accuracy is ~identical while the coefficients (the basis quality
`||R11^{-1}||` controls) differ between methods. (The accuracy-only application run
gets the same coefficient row.)

**Conditioning companion — where `||R11^{-1}||` matters.** Because the projection
error depends only on the *span* of the selected columns, `||R11^{-1}||` is
invisible to it (both methods come out near-optimal above). This companion measures
the basis-dependent quantities a CUR / interpolative-decomposition pipeline
actually pays for, over four rows: `||R11^{-1}||_F / bound`, the interpolation-
coefficient magnitude `max|R11^{-1}R12|`, the **rank-k ID reconstruction error**
(row 3, with the oblique leading-`k` coefficients `T` — above the projection optimum
by a factor growing with `||R11^{-1}||`), and the **orthogonal-projection error**
(row 4, the conditioning-blind best fit in `span(A(:,S))`, near-identical per method).
Both error rows include the black-dotted **best-possible rank-k error (SVD
truncation)**, the absolute lower bound below any column selection — so the story
reads as: projection (row 4) ≈ SVD optimum for both methods, but the oblique ID
coefficients (row 3) lift `rejection_rpqr` well above BSQR. (Lines are seed medians;
shaded bands are seed min/max. The runner also records a noisy-ID variant, not
plotted.) Families are the leverage profiles that make candidate columns near-collinear:
`gaussian` (control), `spiked_leverage`, `collinear_cluster`. Swept over **every
`k` up to 80** (linear axis); the prescribed spectrum has a sharp cliff at `k~50`
so the ID-error rows show a clear knee.

```matlab
matlab -batch "addpath('matlab_rand'); addpath('matlab_rand/benchmark'); run_approx_cond_comparison; plot_approx_cond_comparison"
```

Writes `results/exp_approx_cond.csv` and `plots/fig_approx_cond_quality.{png,pdf}`. The
CSV records every conditioning/error metric in both the Frobenius and spectral
(2-) norms; `plot_approx_cond_comparison('norm','2')` draws the spectral version
(`||R11^{-1}||_2` and the 2-norm approximation errors; `max|R11^{-1}R12|` stays
max-norm) as `plots/fig_approx_cond_quality_spec.{png,pdf}`.

The figure plots the cheap leading-k-frame coefficients `T = R11^{-1}R12`; the
standard least-squares/projection coefficients `T_proj = A(:,S)^+ A(:,rest)`
are **recorded in the CSV** (columns `maxTproj`, `noisy_id_*_proj`) but not
plotted. The takeaway from the recorded data: BSQR's `||R11^{-1}||` is
guaranteed `<=` the Osinsky bound while `rejection_rpqr`'s is 2–3× larger and
can exceed it, which directly inflates the **V_k-frame** coefficients and
their ID error. But `max|T_proj|` is actually *smaller* than `max|T|` and
nearly equal across methods, so with projection coefficients the
reconstruction gap largely closes — the conditioning guarantee is real, but it
shows up downstream when you use the cheap oblique coefficients rather than a
least-squares solve.

**Real data** is never committed. To include a real matrix, drop a `.mat` holding
a 2-D double variable `A` into `ext_comparisons/data/` (git-ignored) and pass its
name as a family:

```bash
mkdir -p ext_comparisons/data
# e.g. a matrix from the SuiteSparse Matrix Collection (https://sparse.tamu.edu/):
#   download <name>.mat (its 'Problem.A' field), save the matrix as variable A in
#   ext_comparisons/data/<name>.mat. (The ARP 'genetics' processed_data.mat also
#   works if you have it.)
```

```matlab
run_approx_comparison('families', {'gmm_kernel', '<name>'});   % missing real names are skipped with a note
```

## Outputs (`[p, Q, R11, stats, R12] = bsqr_rand(A, ...)`)

For `A` `m×n` and `k` selected columns (default `k = min(m,n)`):

| output | shape | meaning |
|---|---|---|
| `p` | `1×n` | permutation; `p(1:k)` is the selected subset (`A(:,p(1:k)) = Q*R11`), `p(k+1:n)` the unselected columns |
| `Q` | `m×k` | economy orthogonal factor; formed lazily — only when the output is requested — from the kernel's accumulated reflectors (LAPACK `dorgqr` in the MEX), so `p = bsqr_rand(A)` never pays for it |
| `R11` | `k×k` | upper-triangular factor of `A(:,p(1:k))` |
| `stats` | struct | instrumentation (see the `stats` section below) |
| `R12` | `k×(n−k)` | coupling block `Q'*A(:,p(k+1:n))`; **only** with `'return_r12', true` (an extra `O(nk^2)` pass) |

## Options (`bsqr_rand(A, 'name', value, ...)`)

| option | default | meaning |
|---|---|---|
| `k` | `min(m,n)` | columns to select |
| `batched` | `true` | in-block BSQR: many selections per sampled block, amortizing the per-block reflector apply (`O(k^3)` vs the single-select `O(k^4)`). `false` = one selection per block (tighter realized conditioning, more applies). Same bound either way |
| `block_size` | `k` (batched) / `ceil(k/2)` in `[16,64]` (single) | candidates evaluated per sampled block; larger improves realized quality at higher cost (`rand_default_block`) |
| `threshold_mode` | `running_mean` | the per-step acceptance threshold — strict at every step against the actual running `f2`, no slack carried between steps (per-singular-value control) — or `worstcase_allowance` (more permissive: spends accumulated slack, fewer samples, same final bound) |
| `slack` | `1.0` | `>=1` multiplier loosening the threshold |
| `norm_recomp_tol` | `sqrt(eps)` | running-norm recompute safeguard in `[0,1]`, matching the deterministic kernel. Used only by the MEX batched path (which downdates norms incrementally); the m-file recomputes exactly, so it is a no-op there |
| `sampling` | `normweighted` | by starting squared column norms (robust across leverage profiles; adds `O(mn)`), or `uniform` |
| `pick` | `best_in_block` | single-select only (`batched=false`): or `first`. Batched always takes the in-block minimizer |
| `seed` | `[]` | RNG seed for reproducibility |
| `return_r12` | `false` | compute `R12` as a 5th output |
| `backend` | `auto` | `auto` / `mfile` / `mex` |
| `check_finite` | `true` | validate inputs |

## `stats` (4th output)

Per-step `1×k` arrays: `f2` (running `||R11^{-1}||_F^2`), `Fhat` (per-step
worst-case bound), `crit`, `threshold`, `samples_tested`, `rounds`, `fallback`.
Scalars: `frob_inv`, `osinsky_bound`, `total_tested`, `blocks_sampled`. In
batched mode a block's `samples_tested`/`rounds` are attributed to its first
selection (`0` for the rest), so `total_tested` and `blocks_sampled` are the
meaningful aggregates.
