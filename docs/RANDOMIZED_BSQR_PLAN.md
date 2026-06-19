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

- apply the accumulated reflectors to the sampled block — `O(k i)` per candidate
  (BLAS-2 `dgemv`/`dger`, or one `dtrsm` for the whole block's `w`-solve);
- reduce the one chosen column — `O(k(k-i))` amortized.

Total `~ O(sum_i S_i k i) + O(k^3) + O(nk)` (the last to read columns / seed
weights), versus the deterministic `O(nk^2)`. With the threshold keeping
`S_i = O(1)`, that is `O(nk + k^3)` — a ~`k×` win when `k^2 << n`.

Measured (Apple Silicon, MEX, short-wide orthonormal rows, block_size 16):

| size       | t_rand (ms) | t_det (ms) | speedup | ratio (||Rinv||/Osinsky) | tested/k |
|------------|------------:|-----------:|--------:|-------------------------:|---------:|
| 64×4000    | 4.1         | 12.9       | 3.2×    | 0.29                     | 16       |
| 64×8000    | 5.9         | 33.2       | 5.6×    | 0.32                     | 16       |
| 128×16000  | 23.7        | 112.1      | 4.7×    | 0.29                     | 16       |

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

## 10. Future work / not yet done

- **Blocked WY apply.** The candidate-apply is currently a per-reflector BLAS-2
  loop; accumulating reflectors in a compact-WY block and applying a panel to a
  candidate block via `dgemm` (mirroring `docs/P3_BLOCKED_BSQR.md`) would cut the
  apply traffic further. The current speedups are already with the simple form.
- **Current-`rho^2` importance sampling.** `normweighted` uses *starting* column
  norms (the only weights available without touching all columns); the
  theory-clean weight is the current `rho_j^2`, which is unavailable cheaply.
  Worth measuring how much the starting-norm proxy helps `tested/k`.
- A perf gate / publication-style plot script if this graduates from experiment.
