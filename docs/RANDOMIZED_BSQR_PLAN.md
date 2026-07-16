# Randomized Bischof–Stewart Column Selection

Design note and status for the experimental randomized variant in `matlab_rand/`.
The deterministic Julia/MATLAB kernels are **not** touched by any of this.

> **Update — in-block selection is now the default (`batched`).** When a block is
> sampled it is brought to the current frame once, then BSQR is run *within* the
> block, taking as many columns as the per-step bound allows before resampling
> (the running `f2`/threshold are updated after **each** in-block selection). One
> reflector apply yields many selections, cutting the dominant cost from
> `O(k^4)` to `O(k^3)`; the MEX maintains the block's `w`-vectors / running norms
> incrementally (BLAS-2). The in-block residual-norm downdate carries the same
> Businger–Golub recompute safeguard as the deterministic kernel
> (`norm_recomp_tol`, default `sqrt(eps)`): once a running norm decays past that
> fraction of its last exact value it is recomputed from scratch (inline, since
> the in-block apply has already updated the block column — unlike the deterministic
> *panel*, which defers to a flush). The single-select and global-min paths and the
> m-file kernel recompute norms exactly each step and need no guard. Measured: **1.5–4.6× faster than single-select** and
> **faster than `rejection_rpqr` at every tested point** (1.75–13.7× on the
> committed publication run, medians across all (family, n, k)) while keeping
> `||R11^{-1}||_F` far under the bound (rejection_rpqr can exceed it). `batched=false`
> recovers the single-select path. The default block is `k`; larger blocks trade
> speed for realized conditioning closer to single-select. The single-select
> sections below are retained as background; their measurements predate batching.

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

The per-step increment when column `j` is appended is exactly (manuscript
eq. (2.4))

```
c_j = (1 + ||w_j||^2) / rho_j^2,   w_j = R11^{-1} R12(:,j),   rho_j = ||tilde m_j||.
```

## 3. The acceptance rule

With `i` columns already selected (0-based), accept a sampled candidate `j` iff

```
c_j <= theta_i,     then set  f <- f + c_j.
```

Two threshold modes are implemented:

- **`running_mean`** (default — the per-step acceptance threshold; strict at
  every step against the actual running `f2`, no slack carried between steps;
  gives per-singular-value control):
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
      hierarchy (manuscript Cor. 3.2 `cor:osinsky` / Thm. 3.1 `thm:BSpivot`;
      the randomized statement is Thm. 5.1 `thm:randBSpivot`).

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
with `slack > 1` the Osinsky guarantee no longer applies).

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

Measured (2026-07-13 Phase 6 rerun on an idle machine: Apple Silicon, MEX,
k=64, `gaussian` orthonormal rows, batched default block `= k`, norm-weighted
sampling; `t_rand` times the default `[p, Q, R11]` product and the
deterministic baseline times `[Q, R, p]` with the vector permutation — every
timing via `timeit`; seed medians):

| size       | speedup vs det | speedup vs dgeqp3 |
|------------|---------------:|------------------:|
| 64×8000    | 53×            | 19×               |
| 64×32000   | 110×           | 25×               |
| 64×64000   | 148×           | 52×               |

`t_rand` is nearly flat across `n` while both baselines grow `O(nk^2)`, so the
speedup grows with `n` while the conditioning stays far under the guarantee.
The comparison uses norm-weighted (column-norm) sampling, which is robust
across leverage profiles, so the stress families land in the same range at the
top size (spiked_leverage 117×, needle 134× vs det); its `O(mn)` norm
precompute is the only `n`-growing term in `t_rand`. Uniform sampling skips
that precompute and is faster still on benign input, but is the right choice
only when leverage is known flat (see the largen study).

## 5. Outputs (`matlab_rand/bsqr_rand.m`)

Default product (the timed path): `[p, Q, R11]`.
- `p` — permutation row vector; `p(1:k)` is the selected column subset.
- `Q` — `m×k` economy orthogonal factor with `Q'*A(:,p(1:k)) = R11`. Formed
  lazily (only when the output is requested) from the kernel's accumulated
  Householder reflectors — LAPACK `dorgqr` in the MEX, one `O(m k²)`
  n-independent pass. `p = bsqr_rand(A)` never pays for it.
- `R11` — `k×k` upper-triangular factor of `A(:,p(1:k))`.

Optional: `[..., stats]` (instrumentation, always cheap) and `[..., R12]`
**only** with `'return_r12', true` (an extra `O(n k^2)` `Q'·A(:,unsel)` pass —
off and untimed by default).

## 6. Knobs

`k`, `batched` (default `true` — in-block BSQR; `false` = single-select),
`block_size` (default `k` when batched, else `ceil(k/2)` in `[16,64]`),
`threshold_mode` (default `running_mean`), `slack` (default 1),
`sampling` (`normweighted` default — by *starting* squared column norms, adds an
`O(mn)` precompute — | `uniform`), `pick` (single-select only: `best_in_block`
default | `first`), `seed`, `return_r12`, `backend` (`auto`/`mfile`/`mex`),
`check_finite` (default `false`; the O(mn) scan can rival the selection cost,
and norm-weighted sampling detects non-finite input for free from its
precomputed weights), `norm_recomp_tol` (default `sqrt(eps)` — the batched path's
inline norm-recompute safeguard; no-op in the m-file, which recomputes exactly).

## 7. Instrumentation (the two requested metrics)

`stats` (per-step `1×k` arrays plus scalars):
- **bound on `||R11^{-1}||_F`**: `f2` (running squared), `Fhat` (per-step worst
  case), `frob_inv = sqrt(f2(end))`, `osinsky_bound = sqrt(k(n-k+1))`;
- **columns tested**: `samples_tested`, `rounds`, `total_tested`, `blocks_sampled`
  (batched attributes a block's `samples_tested`/`rounds` to its first selection,
  so `total_tested` / `blocks_sampled` are the meaningful aggregates);
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
- **Deterministic baseline:** `bsqr_mex(M,'k',k)` with three outputs `[Q,R,p]`,
  matching the publication convention (Q, R, p materialized on every side) and
  the built-in baseline. Forming Q is one `O(k³)` `dorgqr` (m = k in the GKS
  setting), invisible next to the kernel's `O(nk²)` scan; R12 arrives as an
  unavoidable byproduct of that scan (the work the randomized variant skips).
- **Randomized:** `bsqr_rand_mex(...)` with three outputs `[p, Q, R11]` — the
  "R12 not needed" product. Q is the economy factor via `dorgqr` (`O(mk²)`,
  n-independent); no R12.

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

Together these took `gaussian` k=64 from ~5× to two orders of magnitude over the
deterministic baseline (~100× with the norm-weighted comparison default, several ×
more with uniform sampling on benign input; see §4): the kernel apply is now
n-independent while both baselines grow `O(nk^2)`.

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

The experiments run the **default in-block (`batched`) path with `block_size = k`**
(the apply-amortizing sweet spot). All cross-method comparisons (scaling vs
deterministic/`dgeqp3`, the quick benchmark, and the `rejection_rpqr` comparison
in §13) run BSQR with **norm-weighted (column-norm) sampling** — robust across
leverage profiles and the same sampling information the baselines use; the
`fig_sampling` / `fig_blocksize` studies are the exception, since they exist
precisely to compare uniform vs norm-weighted (and block size). The
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

### Findings (Apple Silicon, MEX; refreshed 2026-07-14 from the committed
### publication run — 20 seeds, batched in-block default with block = k;
### numbers are k=64 seed medians unless noted)

- **Speed & scaling.** With norm-weighted sampling the speedup over the
  deterministic factor path grows with `n` and is **consistent across families**.
  On `gaussian` (k=64) it is `53× → 110× → 148×` at `n = 8000 → 32000 → 64000`,
  and the stress families reach the same order at n=64000 (spiked_leverage 117×,
  needle 134×; norm-weighted keeps `tested/k` at a small constant everywhere).
  Uniform sampling is faster still on benign input (it skips the `O(mn)` norm
  precompute) but wastes samples on concentrated leverage — hence norm-weighted
  is the comparison default.
- **Beats the vendor library, not just our own kernel.** Against built-in
  `dgeqp3` (n=64000) randomized is `52× / 57× / 90×` faster at `k = 64 / 128 / 256`
  (gaussian; the stress families are similar) — and `dgeqp3` is the *tougher*
  baseline (it outpaces our deterministic BSQR). Even with `R12` formed it stays
  `~4.5×` faster than `dgeqp3` (k=64, n=64000). Selection quality is well under
  the bound: on the stress families the realized ratio is small for both
  (means ≈ 0.06 randomized vs ≈ 0.03 for `dgeqp3` and deterministic BSQR); on
  `gaussian` it is looser (≈ 0.37 vs ≈ 0.20) for the same speed advantage. The
  value proposition is speed at preserved guarantees.
- **When R12 *is* needed** (the `rand_r12` series / `fig_scaling_*` purple line),
  the randomized method must apply the accumulated `Q'` to the `n-k` leftover
  columns — one `O(nk^2)` BLAS-3 pass, the same order as the deterministic
  kernel's total work. The dramatic n-scaling therefore collapses to a constant
  factor: at k=64 `gaussian` `t_det/t_rand` settles at `~11× → 11× → 13×` over
  `n = 8000 → 32000 → 64000` (vs `53× → 110× → 148×` without R12). It still wins —
  one clean `Q'`-apply versus the deterministic kernel's per-step W-maintenance
  over every column — but the advantage here is a modest constant, not orders of
  magnitude. The big win is specifically the *R12-not-needed* (subset + Q
  + R11) use case.
- **Conditioning.** `running_mean` keeps `||R11^{-1}||_F` far under the bound on
  *all* families (mean ratio ≈ 0.06 on the stress families and ≈ 0.37 on
  `gaussian` at the batched default block = k; per-seed range 0.01–0.49) and
  within ≈ 2–3× of the deterministic value. The bound is never violated.
- **Block size trades compute for realized quality (now the main `k`-knob).**
  With `pick='best_in_block'` the accepted pivot is the block's minimum-criterion
  column, so a larger block minimizes over more candidates and drives the realized
  `||R11^{-1}||_F` toward the deterministic greedy — the *guaranteed* per-step
  bound is unchanged (it is enforced by the threshold at any block). On `gaussian`
  the ratio improves `0.30 → 0.27 → 0.24 → 0.20` at block `16 → 32 → 64 → 128`,
  roughly halving the gap to deterministic (≈ 0.17). Cost grows ~linearly in block
  (and on the stress families very small blocks are much slower — block=1 needs
  ~130–190 rounds). The default is now **`ceil(k/2)` clamped to `[16,64]`** (=
  32/64/64 for k = 64/128/256): a deliberate shift toward selection quality, at
  the price of some raw speed. The cap of 64 keeps that cost bounded for large k:
  at k=256 the default block 64 gives 26× vs `dgeqp3` (conditioning 0.22), whereas
  block 128 roughly halves the speed for only marginally better conditioning
  (0.20). *(These block numbers are for the single-select path, from the
  earlier 5-seed sweep — retained as the block-size intuition. The in-block
  `batched` default — see the note at the top — instead amortizes the apply over
  many selections, so its sweet spot is the larger `block = k`; it supersedes the
  "adaptive block" once listed here as future work.)*
- *(`worstcase_allowance`, omitted from the plots, is a documented option that
  only guarantees the final Frobenius bound: on `spiked_leverage` it spends slack
  up to ~0.8 of the Osinsky bound vs ~0.1 for `running_mean`. Hence `running_mean`
  is the default and the sole plotted threshold.)*
- **Norm-weighted sampling — now performant, and decisive under stress.** Backed
  by a Fenwick tree (sample ∝ weight and remove in `O(log n)`; without-replacement
  via zero/restore), a step costs `O(tested · log n)` with no full sort. On
  concentrated-leverage families at n=32000 it cuts columns tested dramatically
  (spiked 546→2, needle 793→1 tested per selected column at the batched default
  block = k), improves conditioning, **and** is 9–14× faster than uniform in
  wall-clock (needle 11.6 ms → 0.8 ms; spiked 7.7 ms → 0.8 ms). Its only cost is
  the one-time `O(mn)` norm precompute, which makes it ~5× slower than uniform on
  *benign* families (gaussian 0.86 ms vs 0.18 ms — where uniform already accepts
  in one block). **Rule of thumb:
  uniform for benign/unknown-uniform leverage, norm-weighted when leverage is (or
  may be) concentrated.**
- **Across k (k ∈ {64, 128, 256}, n=64000).** The story holds at every k:
  - *Gaussian speedup vs deterministic* stays in the ~100–160× range
    (`148× → 103× → 159×`); versus `dgeqp3` it *grows* with k:
    `52× → 57× → 90×`. The stress families track gaussian closely
    (norm-weighted sampling).
  - *Realized conditioning* is `0.37 → 0.34 → 0.30` (gaussian means) —
    comparable across k and well under the bound; on the stress families it is
    ≈ 0.06 vs ≈ 0.03 for the deterministic selection.

## 12. Review-pass outcomes and future work

Done:
- **In-block (`batched`) selection is now the default** (see the note at the top):
  one reflector apply per block yields many selections, `O(k^3)` vs `O(k^4)`;
  ~1.5–4.6× over single-select and faster than `rejection_rpqr`. The in-block
  applies are restricted to rows `nsel..m-1` (matching the deterministic kernel),
  a further ~10–20%. **Norm-weighted sampling is now the default** too.
- **Rank guard** (both backends): a non-finite per-step increment (every remaining
  residual exactly zero) now errors `bsqr_rand:RankDeficient` instead of
  propagating `Inf` into `f2`/`R11`. *Caveat:* this only catches exact degeneracy;
  near-rank-deficiency (residual `~eps`, not exactly 0) is a documented
  precondition (`k ≤ numerical rank`), not a guarded case.
- **Test coverage** at 28: the MEX `R12` (compact-WY) path, MEX `pick='first'`,
  the rank-deficiency guard (both backends), and the batched in-block path
  (`testBatchedInBlock`, both backends / samplings / blocks).
- Comment audit (m-file says BLAS-2/rank-1, MEX says BLAS-3 compact-WY — accurate)
  and the `slack` docstring caveat (>1 breaks the proven bound).

Still open:
- **Exact-threshold ties — the sampling-cost worst case, and an optional tolerance.**
  The bound is enforced by the strict test `c <= theta`. Generic inputs never tie,
  but a perfectly structured one can sit *exactly* on the threshold: a normalized
  Hadamard (k random rows of an N-by-N Hadamard, /sqrt(N) -> orthonormal rows,
  uniform leverage k/N) has `c_j = (1+0)/(k/N) = N/k = theta_0` for *every* column
  at step 0, and its quantized column overlaps keep the criterion pinned at theta.
  At certain `N/k` (empirically the odd powers of two) rounding then tips every
  remaining column a hair above theta at one step; no sampled block contains an
  acceptable column, the global-min fallback fires once, and the sampler degenerates
  to a single O(n) scan (`tested/k ~ n/k`) for that step — still correct and bounded,
  just not cheap. (It is intermittent: even powers of two sail through at
  `tested/k = 2`. Leverage is irrelevant here — gaussian, equally flat in leverage,
  never ties.) A defensible practical fix is to relax the test to
  `c <= theta * (1 + p*eps)` for a small integer `p`: rounding-level ties then pass
  normally and the bound is exceeded only at the rounding level. We don't ship a
  Hadamard test family (its conditioning just duplicates gaussian and the artifact
  flips with `n`), but the worst case is worth knowing.
- **Persistent MEX workspace.** `bsqr_rand_mex` allocates its scratch per call
  (unlike `bsqr_mex`'s static workspace), so under `timeit` the allocation is
  *inside* the timed region — the reported speedups are mildly conservative. A
  persistent workspace would remove the churn (and make the win look larger).
- **Deferred BLAS-3 in-block update.** The in-block reflector applies are BLAS-2;
  a dlaqps-style deferred update (as in the deterministic *panel* kernel) would
  make them BLAS-3 — a real but non-trivial lever.
- **Adaptive block size.** Within batched, grow the block on consecutive misses to
  get small-block efficiency when columns are easy and large-block amortization
  when scarce.
- **Auto sampling selection.** Cheaply estimate leverage concentration and pick
  uniform vs norm-weighted automatically, so the `O(mn)` precompute is only paid
  when it helps.
- **Current-`rho^2` importance sampling / proper `rank_stop`** (graceful early
  stop returning `ksteps < k`) rather than erroring.
- A perf gate if this graduates from experiment.

(Done earlier, see §10: BLAS-3 compact-WY block apply, `O(tested)` sampling, and
a Fenwick-tree weighted sampler.)

## 13. External comparison: `rejection_rpqr` (Adaptive Randomized Pivoting)

`run_rpqr_comparison.m` / `plot_rpqr_comparison.m` compare `bsqr_rand` against the
`rejection_rpqr` selector from Epperly's Adaptive Randomized Pivoting
(download into `ext_comparisons/`, git-ignored — see `matlab_rand/README.md`).
Same three metrics, separate figures `fig_rpqr_{time,speedup,quality}_k<K>.png`.

Both methods sample by squared column norms (leverage), so BSQR is run with
`normweighted` sampling — the apples-to-apples choice (and consistent with the
main suite, where all cross-method comparisons use column-norm sampling). Two
honest caveats: `rejection_rpqr` is the published `.m` implementation (vs our
MEX), and the two optimize *different objectives* — BSQR minimizes the growth of
`||R11^{-1}||_F`, while ARP targets a volume/DPP criterion — so the `||R11^{-1}||_F`
metric measures BSQR's objective.

`rejection_rpqr` is run at its default `l = k` (the proposal size — its analog of
our block; the ARP authors' own `arp.m` uses this default, and it is near its
speed sweet spot: it is *faster* with larger `l`, not smaller).

Findings (n=32000 unless noted; committed publication run, 20 seeds, batched
default block = k):
- **Conditioning (BSQR's objective) — guaranteed vs uncontrolled.** This is the
  headline. BSQR's `||R11^{-1}||_F / Osinsky` is tight and *always ≤ 1*
  (gaussian means ≈ 0.30–0.38 across k; smaller on the stress families).
  `rejection_rpqr`'s is much larger and **uncontrolled** — gaussian means
  ≈ 0.82–0.92 with high variance, and individual selections routinely *exceed*
  the Osinsky bound (up to ~2.9× at k=64, ~8.7× on spiked_leverage k=128, ~3.2×
  at k=256). ARP optimizes a volume/DPP criterion, not this one, so it has no
  guarantee here; on gaussian BSQR's mean ratio is ~2.5–3× smaller, and it is
  provably bounded where `rejection_rpqr` is not.
- **Runtime — faster at every tested point, with the variance on rejection's
  side.** BSQR's runtime is stable seed-to-seed; `rejection_rpqr`'s varies more
  (its round count is itself random — e.g. 3.8–6.1 ms across seeds at k=128)
  because of the rejection accept/reject. On the batched default BSQR is
  `2.3–2.5×` faster at k=64, `2.0–2.4×` at k=128, and `1.7–1.8×` at k=256
  (medians); across every tested `(family, n, k)` point the median ratio is
  1.75–13.7×, and the gap is asymptotically stable at ~2× in the large-n study
  (`run_largen_scaling`).
- **Bottom line:** at consistently faster (and steadier) runtime, BSQR delivers
  substantially better-conditioned subsets *with a guarantee* `rejection_rpqr`
  does not provide.

## 14. External comparison: low-rank approximation quality

`run_approx_comparison.m` / `plot_approx_comparison.m` add an application-facing
companion to §13 that scores the **approximation error** of each method's column
subset (Frobenius *and* spectral norms) rather than `||R11^{-1}||_F` — a fairer,
practitioner-facing metric (the conditioning metric is BSQR's own objective).

It uses the standard CSSP-via-leading-singular-vectors pipeline (as in ARP's own
`arp.m`): build a full application matrix `A`, compute accurate leading right
singular vectors `V_k` with `svds` (Lanczos), hand the *same* `W = V_k'` to both
selectors, and score the chosen `k` columns by the interpolation-free
orthogonal-projection error `||A − P_S A||` against the optimal rank-`k` error.
Sharing `V_k` isolates the column-selection step — subspace estimation is not the
point. Synthetic families (`gmm_kernel`, `integral_skeleton`, `snapshots`) are
reproducible and committed; real matrices live in `ext_comparisons/data/`
(git-ignored, with download instructions in `matlab_rand/README.md`). It sweeps
the rank `k` on a fixed matrix per family and writes `results/exp_approx.csv` plus
`plots/fig_approx_quality.{png,pdf}`. `bsqr_rand` runs with its public defaults
(batched, norm-weighted).

A **synthetic-spectrum companion** (`run_approx_synth_comparison.m` →
`approx_synth_matrix.m`; `plots/fig_approx_synth_quality.{png,pdf}`) asks how much
the §13 `||R11^{-1}||_F` differences translate into approximation error: it builds
`A = U diag(s) V'` with a prescribed interesting spectrum `s` (few large, decay,
flatter section, more decay) and right singular vectors `V` from the *same*
leverage families used in §13 (`gaussian`, `spiked_leverage`, `needle`, via
`rand_test_matrix`). The shared spectrum makes the optimal-error curve identical
across the three families, isolating the effect of leverage structure.

A **conditioning companion** (`run_approx_cond_comparison.m` /
`plot_approx_cond_comparison.m`; `plots/fig_approx_cond_quality.{png,pdf}`) closes
the loop: since the projection error depends only on the selected *span*,
`||R11^{-1}||` is invisible to it, so this run measures the basis-dependent
quantities a CUR/ID pipeline pays for — `||R11^{-1}||_F/bound`, the
interpolation-coefficient magnitude `max|R11^{-1}R12|`, and the rank-k ID
reconstruction error split into a *noiseless* part (oblique-coefficient penalty,
already amplified by `||R11^{-1}||` above the orthogonal-projection lower bound) and
a *noisy* part (measurement noise propagated through `T ~ ||R11^{-1}||`).
Families are near-collinear leverage profiles (`gaussian` control, `spiked_leverage`,
and `collinear_cluster` in `rand_test_matrix`), swept over a dense linear `k` grid
with a sharp-cliff spectrum (a clear knee in the ID-error rows). Result: projection accuracy is
identical, but BSQR's `||R11^{-1}||` stays `<=` the Osinsky bound (guaranteed)
while `rejection_rpqr`'s runs 2–3× larger and can exceed it, inflating its
interpolation coefficients (≈2–4×) and noise-amplified ID error correspondingly.
Computed with no `inv()` (conditioning via `svd`, coefficients via backslash).

The reconstruction rows compare two coefficient choices: the cheap leading-k-frame
`T = R11^{-1}R12` and the standard projection coefficients
`T_proj = A(:,S)^+ A(:,rest)` (columns `maxTproj`, `noisy_id_*_proj`; recorded
in the CSV, not plotted). `max|T_proj|` turns out *smaller* than `max|T|` and nearly method-agnostic,
so the conditioning gap inflates the V_k-frame ID but **largely closes under
projection coefficients** — the guarantee matters downstream with the cheap oblique
coefficients, not with a least-squares solve. Every metric is recorded in both
norms; the spectral figure is `fig_approx_cond_quality_spec.{png,pdf}`
(`plot_approx_cond_comparison('norm','2')`).
