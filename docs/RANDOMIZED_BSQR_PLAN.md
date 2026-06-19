# Randomized Bischof–Stewart Column Selection

Design note and status for the experimental randomized variant in `matlab_rand/`.
The deterministic Julia/MATLAB kernels are **not** touched by any of this.

## 1. Idea

The deterministic kernel selects, at every step, the exact argmin of the
Bischof–Stewart criterion `c_j = (1 + ||w_j||^2) / rho_j^2`, which requires
maintaining `R11^{-1}R12` and the running column norms for **all** trailing
columns at every step — the `O(n k^2)` work that dominates the short-wide
regime. The guarantees, however, do not need the argmin: they only need each
step's increment to keep the running `||R11^{-1}||_F^2` under a per-step bound.
So we can instead **sample** a few candidate columns, bring them up to date with
the accumulated reflectors, and accept the first/best one that satisfies the
bound — never maintaining state for the columns we don't look at.

## 2. The quantity tracked: `f` is the *squared* inverse Frobenius norm

Throughout, **`f_i = ||R11^{(i)}^{-1}||_F^2`** (squared). The recurrence below
only closes on the squared quantity; the code names it `f2` to keep this
explicit. Seed `f_0 = 0`.

The per-step increment when column `j` is appended is exactly (writeup
eq. frob-growth)

```
c_j = (1 + ||w_j||^2) / rho_j^2,   w_j = R11^{-1} R12(:,j),   rho_j = ||tilde m_j||.
```

## 3. The acceptance rule

With `i` columns already selected (0-based), accept a sampled candidate `j` iff

```
c_j <= theta_i,     then set  f <- f + c_j.
```

Two threshold modes are implemented:

- **`running_mean`** (default — gives per-singular-value control):
  ```
  theta_i = (f_i + n - 2i) / (k - i)
  ```
  This is exactly the per-column bound the user specified:
  `f_{i+1} <= ((k-i+1)/(k-i)) f_i + (n-2i)/(k-i)`. For orthonormal-row input it
  equals the `rho^2`-weighted **mean** of `c_j` over the remaining columns
  (because `sum rho_j^2 = k-i` and `sum ||w_j||^2 = f_i - i`), so:
    - the minimizer always satisfies it (min <= mean) ⇒ an exhaustive pass can
      never fail to find an acceptable pivot (termination/feasibility);
    - maintaining it preserves `f_i <= Fhat_i = i(n-i+1)/(k-i+1)`, hence
      Osinsky's `||R11^{-1}||_F^2 <= k(n-k+1)` **and** the per-singular-value
      hierarchy (writeup Thm. 3.4 / Cor. 3.2).

- **`worstcase_allowance`** (more permissive, fewer samples):
  ```
  theta_i = Fhat_{i+1} - f_i,    Fhat_m = m(n-m+1)/(k-m+1).
  ```
  Since the algorithm tracks the *real* `f_i` (usually below the worst case
  `Fhat_i`), the gap `theta_wc - theta_running = (k-i+1)(Fhat_i - f_i)/(k-i) >= 0`
  is accumulated slack that can be spent. This admits more columns ⇒ accepts
  with fewer samples, and still provably keeps `f_i <= Fhat_i` (so Osinsky's
  bound still holds), at the cost of a looser final `||R11^{-1}||_F`.

A `slack` multiplier (`>= 1`) loosens either threshold further (experimental;
breaks the exact bound proportionally).

## 4. Why it is faster

Per accepted pivot, with `S_i` candidates tested at step `i`:

- apply the accumulated reflectors to the sampled block — done as one BLAS-3
  compact-WY block (`Q_i' = I - V T' V'`: gemm + trmm + gemm) plus one `dtrsm`
  for the block's `w`-solve;
- reduce the one chosen column — reused from the test block (no re-apply).

The MEX touches only the columns it actually tests: there is **no** `O(n)`
per-step pass. Total `~ O(sum_i S_i k i) + O(k^3) + O(n)` (the `O(n)` is a single
pool init; norm-weighted sampling adds a one-time `O(mn)` for column norms),
versus the deterministic `O(nk^2)`. With the threshold keeping `S_i = O(1)`, the
randomized kernel is essentially **n-independent**, so the speedup grows linearly
in `n` for `k << n`.

Measured (Apple Silicon, MEX, k=64, `gaussian` orthonormal rows, default block
`ceil(k/2)=32`, uniform sampling, `t_rand` is `[p,reflectors,R11]` only):

| size       | speedup vs det | speedup vs dgeqp3 |
|------------|---------------:|------------------:|
| 64×8000    | 25×            | 9×                |
| 64×32000   | 159×           | 36×               |
| 64×64000   | 345×           | 127×              |

`t_rand` is essentially flat across `n` while both baselines grow `O(nk^2)`, so
the speedup grows with `n`. (The default block was 16 in earlier revisions, which
gave larger raw speedups — e.g. 475× vs det at 64×64000 — but worse realized
conditioning; the larger default trades some of that speed for selection quality,
see §11.) On concentrated-
leverage inputs (where uniform sampling tests many columns) use `normweighted`
sampling to keep `tested/k ≈ block` and recover the same scaling (§11).

These numbers reflect three landed optimizations (earlier BLAS-2 / `O(n)`-shuffle
versions were ~3–6× only): see §10.

The speedup grows with `n`, exactly as the cost model predicts; the conditioning
stays far under the guarantee.

## 5. Outputs (`matlab_rand/bsqr_rand.m`)

Default product (the timed path): `[p, reflectors, R11]`.
- `p` — permutation row vector; `p(1:k)` is the selected column subset.
- `reflectors` — `struct('V', m×k unit-diagonal reflector store, 'tau', k×1, ...)`;
  materialize Q with `bsqr_rand_formQ`.
- `R11` — `k×k` upper-triangular factor of `A(:,p(1:k))`.

Optional: `[..., stats]` (instrumentation, always cheap) and `[..., R12]`
**only** with `'return_r12', true` (an extra `O(n k^2)` `Q'·A(:,unsel)` pass —
off and untimed by default).

## 6. Knobs

`k`, `block_size` (default 16), `threshold_mode` (default `running_mean`),
`slack` (default 1), `sampling` (`uniform` default | `normweighted` — by
*starting* squared column norms, adds an `O(mn)` precompute), `pick`
(`best_in_block` default | `first`), `seed`, `return_r12`, `backend`
(`auto`/`mfile`/`mex`), `check_finite`.

## 7. Instrumentation (the two requested metrics)

`stats` (per-step `1×k` arrays plus scalars):
- **bound on `||R11^{-1}||_F`**: `f2` (running squared), `Fhat` (per-step worst
  case), `frob_inv = sqrt(f2(end))`, `osinsky_bound = sqrt(k(n-k+1))`;
- **columns tested ("nested")**: `samples_tested`, `rounds`, `total_tested`;
- plus `crit` (accepted `c`), `threshold`, `fallback`.

## 8. Validation (`matlab_rand/tests/`)

- `test_bsqr_rand.m` — factorization is exact (`Q'·A(:,sel) = R11`,
  `A(:,sel) = Q·R11`, `Q` orthonormal), `p` is a permutation, optional `R12`
  correct, edge `k=0`/`k=n`, seed determinism, all knob combinations, and a
  MEX-vs-invariants check (the MEX has its own RNG, so pivots need not match the
  m-file — only the guarantees must).
- `test_bsqr_rand_bounds.m` — `f2_i <= Fhat_i` every step and final
  `||R11^{-1}||_F <= sqrt(k(n-k+1))` for both threshold modes; `f2` matches the
  explicit `||inv(R11)||_F^2`; conditioning within a small factor of the
  deterministic kernel.

Run: `matlab -batch "addpath('matlab_rand'); run('matlab_rand/tests/run_rand_tests.m')"`
(or `runtests('matlab_rand/tests')`).

## 9. Implementation notes

- The reference m-file and the MEX share the algorithm but **not** the RNG, so
  pivot sequences differ between backends (and from MATLAB's global `rng`). This
  is intentional — randomized variant, guarantees not pivots.
- `A` is never permuted or mutated; columns are gathered by original index and
  `remaining` is a shrinking index list reshuffled each step.
- Reflectors are stored in unit-diagonal packed form (`V(t,t)=1`, tail below,
  zeros above); `R11` is kept separately. `formQ` and the candidate-apply use
  the same convention.
- Feasibility caveat: the bound/feasibility derivation assumes orthonormal-row
  input (`MM^T = I`, the GKS setting `m = k`). General matrices run fine and the
  factorization stays exact, but only the orthonormal case is covered by theory;
  if `running_mean` ever yields a negative threshold (impossible for orthonormal
  input), the exhaustive global-min fallback fires and is flagged in `stats`.

## 10. Performance-measurement discipline

The timing harness (`benchmark/run_rand_benchmarks.m`, `run_rand_experiments.m`)
is built to time the *kernels* and nothing else, and to compare fairly:

- The compiled MEX kernels are called **directly** (`bsqr_mex`,
  `bsqr_rand_mex`) — the m-file dispatcher and its `inputParser` are never in the
  timed path.
- `'check_finite', false` on both methods — the `O(mn)` finiteness scan is
  identical overhead for both and is not part of either algorithm.
- Only the kernel call is inside the `timeit` thunk; matrix generation is
  outside. `timeit` supplies warm-up and a robust median.
- **Deterministic baseline:** `bsqr_mex(M,'k',k)` with one output (R only). This
  is the cheapest deterministic call that still performs the selection; it forms
  R12 as an unavoidable byproduct (the `O(nk^2)` work the randomized variant
  skips) but does **not** materialize Q.
- **Randomized:** `bsqr_rand_mex(...)` with three outputs `[p, reflectors, R11]`
  — the "R12 not needed" product. No Q, no R12.

Kernel optimizations landed (MEX), in order of impact for `k << n`:

1. **BLAS-3 compact-WY block apply.** The candidate block is reduced with the
   compact-WY form `Q_i' = I - V T' V'` (one `dgemm` + `dtrmm` + `dgemm`), with
   `T` maintained incrementally (`dlarft` forward recurrence). This replaces the
   previous loop of `nsel` tiny BLAS-2 rank-1 updates.
2. **O(tested) sampling, O(1) removal.** Uniform sampling uses a *partial*
   Fisher–Yates that only shuffles the columns actually drawn, and removes the
   pivot from the pool by swap-and-pop — eliminating the `O(n)` per-step shuffle
   that previously reintroduced an `O(nk)` term and capped the `k << n` speedup.
3. **Reuse the accepted column.** Its reduced form already sits in the test block,
   so it is not re-applied (only the rare exhaustive-fallback re-applies).

Together these moved `gaussian` k=64 from ~5× to **475×** at n=64000 (§4): the
kernel is now n-independent while the deterministic baseline grows `O(nk^2)`.

## 11. Experiment & plot suite (R12 not needed)

`run_rand_experiments.m` → CSVs in `benchmark/results/`; `plot_rand_experiments.m`
→ banded figures in `benchmark/plots/` (line = seed mean, shaded band = seed
min/max). Both are parameterized by `k` and the committed sweep covers
`k ∈ {64, 128, 256}`, written to `k`-tagged files (`exp_*_k<K>.csv`,
`fig_*_k<K>.png`). Stress families are in `rand_test_matrix.m`: `gaussian`
(benign, near-uniform leverage), `graded_leverage`, `spiked_leverage` (few
high-leverage columns), `coherent` (clustered/redundant columns), `needle`
(~k useful columns hidden among n near-null ones — hardest for uniform sampling).
The sampler's cost is set by how concentrated the leverage `ell_j = ||M(:,j)||^2`
is. Peak memory is modest (one `k×n` double is ~0.13 GB at k=256, n=64000; a
handful coexist → ~1 GB peak), so the full grid runs comfortably on 16 GB.

Baselines on every timing plot: the **deterministic BSQR** factor path *and* the
**built-in column-pivoted QR** (`qr(M,'econ','vector')` → LAPACK `dgeqp3`), a
different, vendor-tuned classical algorithm — the natural "can we beat the
library?" reference. (`dgeqp3` is itself faster than our deterministic BSQR, so
it is the tougher baseline.)

The plots use the **default block size** `ceil(k/2)` clamped to `[16,128]`
(`rand_default_block.m`; the MEX mirrors it) and `pick='best_in_block'` (accept
the block's minimum-criterion column when it meets the threshold). The
`worstcase_allowance` threshold stays a documented option (`bsqr_rand` help) but
is **omitted from the plots** to keep the narrative on the per-step-bounded
`running_mean` default.

1. **`fig_scaling_time` / `fig_scaling_speedup` / `fig_scaling_quality`** — the
   scaling study split into three single-metric figures (one former column
   each, so each has a clear legend): runtime vs `n`; randomized speedup over
   both baselines (and the R12-desired path vs `dgeqp3`); and selection quality
   `||R11^{-1}||_F / sqrt(k(n-k+1))`.
2. **`fig_blocksize`** — time, columns tested/`k`, and conditioning vs block
   size, uniform vs norm-weighted sampling (with a vendor-qr time reference).
3. **`fig_sampling`** — uniform vs norm-weighted across all families
   (columns tested/`k`, `||R11^{-1}||_F`, time; vendor-qr time overlaid).

### Findings (Apple Silicon, MEX, 5 seeds; numbers below are k=64 unless noted)

- **Speed & scaling.** With the optimized kernel the randomized method is
  essentially n-independent, so the speedup grows with `n`. On `gaussian` (k=64,
  default block 32) the speedup over the deterministic factor path is
  `25× → 159× → 345×` at `n = 8000 → 32000 → 64000`. On the stress families with
  *uniform* sampling the win is smaller because the sampler tests many columns —
  fixed by norm-weighted sampling (below).
- **Beats the vendor library, not just our own kernel.** Against built-in
  `dgeqp3` (gaussian, n=64000) randomized is `127× / 59× / 18×` faster at
  `k = 64 / 128 / 256` — and `dgeqp3` is the *tougher* baseline (it outpaces our
  deterministic BSQR). Even with `R12` formed it stays `~5–6×` faster than
  `dgeqp3` at n=64000. Selection quality is comparable and well under the bound
  (`||R11^{-1}||_F`/Osinsky ≈ 0.20–0.28 randomized vs ≈ 0.16–0.18 for `dgeqp3` on
  gaussian). The value proposition is speed at preserved guarantees. (Raw speedups
  scale inversely with the block size; the smaller block 16 used in earlier
  revisions roughly doubled–to–6×'d these numbers at the cost of conditioning.)
- **When R12 *is* needed** (the `rand_r12` series / `fig_scaling_*` purple line),
  the randomized method must apply the accumulated `Q'` to the `n-k` leftover
  columns — one `O(nk^2)` BLAS-3 pass, the same order as the deterministic
  kernel's total work. The dramatic n-scaling therefore collapses to a constant
  factor: `gaussian` `t_det/t_rand` settles at `~6× → 10× → 14×` over
  `n = 1000 → 8000 → 64000` (vs `7× → 34× → 517×` without R12), and the stress
  families to ~4–6×. It still wins because it does one clean `Q'`-apply instead of
  the deterministic kernel's per-step W-maintenance over every column — but this
  is the regime where the advantage is a modest constant, not orders of magnitude.
  Bottom line: the big win is specifically the *R12-not-needed* (subset +
  reflectors + R11) use case.
- **Conditioning.** `running_mean` keeps `||R11^{-1}||_F` far under the bound on
  *all* families (ratio ≈ 0.05–0.28 at the default block) and within ≈ 1.2–1.6× of
  the deterministic value. The bound is never violated — the stress matrices
  stress the *sample count*, not quality.
- **Block size trades compute for realized quality (now the main `k`-knob).**
  With `pick='best_in_block'` the accepted pivot is the block's minimum-criterion
  column, so a larger block minimizes over more candidates and drives the realized
  `||R11^{-1}||_F` toward the deterministic greedy — the *guaranteed* per-step
  bound is unchanged (it is enforced by the threshold at any block). On `gaussian`
  the ratio improves `0.30 → 0.27 → 0.24 → 0.20` at block `16 → 32 → 64 → 128`,
  roughly halving the gap to deterministic (≈ 0.17). Cost grows ~linearly in block
  (and on the stress families very small blocks are much slower — block=1 needs
  ~130–190 rounds). The default is now **`ceil(k/2)` clamped to `[16,128]`** (=
  32/64/128 for k = 64/128/256): a deliberate shift toward selection quality, at
  the price of some raw speed (e.g. k=256 vs `dgeqp3` drops from ~107× at block 16
  to ~18× at block 128). An adaptive block remains future work.
- *(`worstcase_allowance`, omitted from the plots, is a documented option that
  only guarantees the final Frobenius bound: on `spiked_leverage` it spends slack
  up to ~0.8 of the Osinsky bound vs ~0.1 for `running_mean`. Hence `running_mean`
  is the default and the sole plotted threshold.)*
- **Norm-weighted sampling — now performant, and decisive under stress.** Backed
  by a Fenwick tree (sample ∝ weight and remove in `O(log n)`; without-replacement
  via zero/restore), a step costs `O(tested · log n)` with no full sort. On
  concentrated-leverage families at n=32000 it cuts columns tested dramatically
  (spiked 544→32, needle 767→32 per k at the default block 32), improves
  conditioning, **and** is 8–11× faster than uniform in wall-clock (needle
  14.3 ms → 1.3 ms; spiked 9.6 ms → 1.3 ms) — ~65× over the deterministic
  baseline. Its only cost is
  the one-time `O(mn)` norm precompute, which makes it ~3× slower than uniform on
  *benign* families (where uniform already accepts in one block). **Rule of thumb:
  uniform for benign/unknown-uniform leverage, norm-weighted when leverage is (or
  may be) concentrated.**
- **Across k (k ∈ {64, 128, 256}, n=64000).** The story holds at every k; the
  numbers shift as the cost model predicts:
  - *No-R12 gaussian speedup* (vs deterministic, n=64000, default block) shrinks
    with k — `345× → 127× → 35×` at k = 64/128/256 — because the randomized apply
    grows with `k` *and* the default block `ceil(k/2)`, while the deterministic
    baseline grows `O(nk²)`. Versus `dgeqp3` it is `127× → 59× → 18×`.
  - *Realized conditioning improves with k* — `0.28 → 0.23 → 0.20` — since the
    default block `ceil(k/2)` minimizes over a larger fraction of candidates.
  - *Stress families* (uniform) still improve slightly with k (`needle`),
    because the useful set `ng ≈ 1.25k` grows with k so fewer samples are wasted.
  - *Conditioning* stays `~0.29–0.33` (gaussian) and lower for stress at all k —
    the Osinsky bound is respected throughout.
  - *Norm-weighted sampling* wins at every k: e.g. `needle` k=256 cuts tested/k
    `207 → 16` and time `105 ms → 10 ms` (~10×).

## 12. Review-pass outcomes and future work

Done in the review pass:
- **Rank guard** (both backends): a non-finite per-step increment (every remaining
  residual exactly zero) now errors `bsqr_rand:RankDeficient` instead of
  propagating `Inf` into `f2`/`R11`. *Caveat:* this only catches exact degeneracy;
  near-rank-deficiency (residual `~eps`, not exactly 0) is a documented
  precondition (`k ≤ numerical rank`), not a guarded case.
- **Test coverage** raised to 18: added the MEX `R12` (compact-WY) path, MEX
  `pick='first'`, and the rank-deficiency guard (both backends).
- Comment audit (m-file says BLAS-2/rank-1, MEX says BLAS-3 compact-WY — accurate)
  and the `slack` docstring caveat (>1 breaks the proven bound).

Still open:
- **Persistent MEX workspace.** `bsqr_rand_mex` allocates its scratch per call
  (unlike `bsqr_mex`'s static workspace), so under `timeit` the allocation is
  *inside* the timed region — the reported speedups are mildly conservative. A
  persistent workspace would remove the churn (and make the win look larger).
- **Adaptive block size.** Grow the block on consecutive misses to get small-block
  efficiency when columns are easy and large-block amortization when scarce.
- **Auto sampling selection.** Cheaply estimate leverage concentration and pick
  uniform vs norm-weighted automatically, so the `O(mn)` precompute is only paid
  when it helps.
- **Current-`rho^2` importance sampling / proper `rank_stop`** (graceful early
  stop returning `ksteps < k`) rather than erroring.
- A perf gate if this graduates from experiment.

(Done earlier, see §10: BLAS-3 compact-WY block apply, `O(tested)` sampling, and
a Fenwick-tree weighted sampler.)
