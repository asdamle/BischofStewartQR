# matlab_rand — randomized Bischof–Stewart column selection (experimental)

A randomized variant of BSQR that selects `k` columns of a `k×n` matrix with the
same theoretical guarantees on `||R11^{-1}||_F` as the deterministic kernel, but
**without** maintaining `R11^{-1}R12` / column norms for every column at every
step. Instead it tracks only the running squared inverse Frobenius norm and, per
step, samples candidate columns and accepts the first/best one that keeps the
running value under the per-step bound. See `docs/RANDOMIZED_BSQR_PLAN.md` for
the math and `notes/bischof_stewart_pivoting.tex` §3.4 for the bound it relies on.

This is **separate from and does not modify** the deterministic implementations
in `matlab/` and `julia/`.

## Quick start

```matlab
addpath('matlab_rand'); addpath('matlab_rand/mex');   % mex auto-added by bsqr_rand

M = orth(randn(2000, 32))';            % 32×2000, orthonormal rows (GKS setting)
[p, reflectors, R11] = bsqr_rand(M);   % default product: subset + reflectors + R11
sel = p(1:32);                         % selected column indices
Q   = bsqr_rand_formQ(reflectors);     % materialize Q on demand
```

`R12` is opt-in (it costs an extra `O(n k^2)` pass and is off by default):

```matlab
[p, reflectors, R11, stats, R12] = bsqr_rand(M, 'return_r12', true);
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

Compares the randomized selection (`[p, reflectors, R11]`, no `R12`) against the
deterministic factor path and reports speedup, the conditioning ratio
`||R11^{-1}||_F / sqrt(k(n-k+1))`, and candidate columns tested per pivot. Writes
`benchmark/results/rand_timings.csv`.

## Options (`bsqr_rand(A, 'name', value, ...)`)

| option | default | meaning |
|---|---|---|
| `k` | `min(m,n)` | columns to select |
| `block_size` | `16` | candidates evaluated per sampling round |
| `threshold_mode` | `running_mean` | per-column bound (per-singular-value control) or `worstcase_allowance` (more permissive, fewer samples) |
| `slack` | `1.0` | `>=1` multiplier loosening the threshold |
| `sampling` | `uniform` | or `normweighted` (by starting squared column norms; adds `O(mn)`) |
| `pick` | `best_in_block` | or `first` |
| `seed` | `[]` | RNG seed for reproducibility |
| `return_r12` | `false` | compute `R12` as a 5th output |
| `backend` | `auto` | `auto` / `mfile` / `mex` |
| `check_finite` | `true` | validate inputs |

## `stats` (4th output)

Per-step `1×k` arrays: `f2` (running `||R11^{-1}||_F^2`), `Fhat` (per-step
worst-case bound), `crit`, `threshold`, `samples_tested`, `rounds`, `fallback`.
Scalars: `frob_inv`, `osinsky_bound`, `total_tested`.
