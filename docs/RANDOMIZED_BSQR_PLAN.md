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

Measured (Apple Silicon, MEX, k=64, `gaussian` orthonormal rows, block 16,
uniform sampling, `t_rand` is `[p,reflectors,R11]` only):

| size       | t_rand (ms) | t_det (ms) | speedup | ||Rinv|| rand/det |
|------------|------------:|-----------:|--------:|------------------:|
| 64×8000    | 0.11        | 6.6        | 42×     | ~1.9              |
| 64×32000   | 0.39        | 95.8       | 238×    | ~1.7              |
| 64×64000   | 0.40        | 220        | 475×    | ~1.9              |

`t_rand` is flat across `n` while `t_det` grows `O(nk^2)`. On concentrated-
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
min/max). Stress families are in `rand_test_matrix.m`: `gaussian` (benign,
near-uniform leverage), `graded_leverage`, `spiked_leverage` (few high-leverage
columns), `coherent` (clustered/redundant columns), `needle` (~k useful columns
hidden among n near-null ones — hardest for uniform sampling). The sampler's cost
is set by how concentrated the leverage `ell_j = ||M(:,j)||^2` is.

1. **`fig_scaling`** — time, speedup `t_det/t_rand`, and conditioning
   `||R11^{-1}||_F / sqrt(k(n-k+1))` vs `n` (per family), deterministic vs both
   threshold modes.
2. **`fig_blocksize`** — time, columns tested/`k`, and conditioning vs block
   size, uniform vs norm-weighted sampling.
3. **`fig_sampling`** — uniform vs norm-weighted across all families
   (columns tested/`k`, `||R11^{-1}||_F`, time).

### Findings (Apple Silicon, MEX, k=64; 5 seeds)

- **Speed & scaling.** With the optimized kernel the randomized method is
  essentially n-independent. On `gaussian` (uniform sampling) the speedup over the
  deterministic factor path grows `9× → 42× → 238× → 475×` as `n` goes
  `1000 → 8000 → 32000 → 64000`. On the stress families with *uniform* sampling
  the win is smaller (≈6–10× at n=64000) because the sampler tests many columns —
  fixed by norm-weighted sampling (below).
- **Conditioning.** `running_mean` keeps `||R11^{-1}||_F` far under the bound on
  *all* families (ratio 0.06–0.32) and ~1.7–2× the deterministic value. The bound
  is never violated — the stress matrices stress the *sample count*, not quality.
- **`worstcase_allowance` is risky under stress.** Fine on benign inputs, but it
  spends slack aggressively: on `spiked_leverage` it rides to ~0.8 of the Osinsky
  bound (vs ~0.1 for `running_mean`). Keep `running_mean` the default; use
  `worstcase_allowance` only when the input is known benign.
- **Block size has a family-dependent sweet spot.** Benign (`gaussian`): small
  blocks (2–8) are fastest (best-in-block tests the whole block; conditioning
  still improves with block size). Stressed (`spiked`/`needle`) with *uniform*
  sampling: small blocks are much slower (block=1 → ~130–190 rounds); ~16–32
  amortizes the round overhead. Default 16; an adaptive block would dominate both.
- **Norm-weighted sampling — now performant, and decisive under stress.** Backed
  by a Fenwick tree (sample ∝ weight and remove in `O(log n)`; without-replacement
  via zero/restore), a step costs `O(tested · log n)` with no full sort. On
  concentrated-leverage families at n=32000 it cuts columns tested dramatically
  (spiked 519→16, needle 765→16 per k), improves conditioning 3–6×, **and** is now
  11–16× faster than uniform in wall-clock (needle: 17.7 ms → 1.1 ms; spiked:
  12.3 ms → 1.1 ms) — i.e. ~88× over the deterministic baseline. Its only cost is
  the one-time `O(mn)` norm precompute, which makes it ~3× slower than uniform on
  *benign* families (where uniform already accepts in one block). **Rule of thumb:
  uniform for benign/unknown-uniform leverage, norm-weighted when leverage is (or
  may be) concentrated.**

## 12. Future work / not yet done

- **Adaptive block size.** Start small and grow the block on consecutive misses
  to get the best of both regimes (small-block efficiency when columns are easy,
  large-block amortization when they are scarce) — would remove the need to pick
  uniform-vs-stress block sizes by hand.
- **Auto sampling selection.** Cheaply estimate leverage concentration (e.g. from
  a small sample of column norms) and pick uniform vs norm-weighted automatically,
  so the benign-case `O(mn)` precompute is only paid when it pays off.
- **Current-`rho^2` importance sampling.** `normweighted` uses *starting* column
  norms (the only weights available without touching all columns); the
  theory-clean weight is the current `rho_j^2`, unavailable cheaply.
- A perf gate if this graduates from experiment.

(Done, see §10: BLAS-3 compact-WY block apply, `O(tested)` sampling, and a
Fenwick-tree weighted sampler.)
