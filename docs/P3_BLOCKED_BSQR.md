# Panel/Blocked BSQR: Method, Analysis, and Experimental Results

Working notes behind Appendix A of the manuscript (`notes/GKSevolved_draft.tex`; P3 of
`docs/VALIDATION_AND_PERF_PLAN.md`; landed and default as of 2026-06-11). Notation follows
Algorithm A.1 of the manuscript, 1-indexed, keeping this file's legacy symbols `s_j`, `β_j`
for the manuscript's `ρ_j²`, `α_j` (see the symbol map in `docs/VALIDATION.md`). A panel
starts at global step `s` and runs local steps `t = 1..nb` (global step `i = s+t-1`);
"prefix" means rows `1..s-1`; `M = m-s+1`, `N = n-s+1` are trailing dimensions at panel start;
`rem` denotes the trailing columns `i+1..n` at the current step.

## 1. Motivation (measured, not assumed)

Phase-level profiling of the unblocked kernel (P0, `docs/PERF_P0_FINDINGS.md`) shows it is two
BLAS-2 walls and nothing else: the trailing Householder application takes 44–57% of kernel
time and the W update 38–53% (dominant in short-wide); pivot scan, norm downdates, and all
other bookkeeping total ≤3%. Both walls are gemv+ger pairs — four big-block memory passes per
step. The control experiment against LAPACK's *unblocked* pivoted QR (`dgeqpf`) splits the gap
to `dgeqp3` into a blocking gap (1.12–1.28 on square, ≈0 on short-wide) and an extra-work gap
matching the flop model `(HH + W + Qmat)/(HH + Qmat)`; short-wide BSQR was already *at* its
flop floor, so only reorganizing the W work to BLAS-3 could move that regime.

The reorganization below extends the `dlaqps` deferred-update scheme (the panel kernel inside
`dgeqp3`) to BSQR's additional state. Per step it leaves one *read-only* gemv pass over each
big block; the two rank-1 write streams accumulate into two rank-`nb` gemms per panel.

## 2. State and accumulators

Unblocked per-step state: `A` (reduced in place; reflectors below the diagonal), `W` (rows
`1..i-1`; `W(1:i-1,j) = R11^{-1}R12(:,j)`), running squared tail norms `s_j`, running
`wnorm2_j = ||W(1:i-1,j)||^2`.

Panel accumulators (workspace, reused across panels):

| symbol | shape | content |
|---|---|---|
| `V` | `M × nb` | reflector vectors `v_t`, stored at global rows with **zeros above row `i`** and unit diagonal; unit-lower-trapezoidal in panel rows |
| `Fa` | `N × nb` | deferred-update factors for A: true trailing block = `A_stored − V·Faᵀ` (rows below the current step) |
| `B` | `nb × N` | the β rows appended to W during the panel, **eagerly maintained** (they receive the panel's later rank-1 updates) |
| `Bsel` | `nb × N` | β rows frozen at selection time (`Bsel(u,j) = β_u[j]`) — the coefficients of the deferred prefix updates |
| `Ω` | `(s-1) × nb` | `Ω(:,t) = ω_t`, the prefix part (rows `1..s-1`) of the pivot's w-vector at its selection step |

The structural point that makes the W side work without ragged updates: the step-`u` rank-1
update `W(1:i_u-1, rem) -= w*_u·β_uᵀ` splits into a *prefix* part (rows `1..s-1`, uniform
height across the panel — deferred via `Ω`/`Bsel`) and a *panel* part (rows `s..i_u-1`, at
most `nb-1` rows — applied eagerly to `B` at `O(nb·N)` per step). The true W during the panel
is therefore

```
W_true(1:s-1, j) = W₀(1:s-1, j) − Ω(:,1:t)·Bsel(1:t, j)      (prefix, deferred)
W_true(s:i-1, j) = B(1:t-1, j)                                (panel rows, current)
```

## 3. Per-step recipe (local step t, global step i)

1. **Pivot selection** — unchanged from the unblocked kernel: `argmin_j (1+wnorm2_j)/s_j`
   over the running scalar arrays, which are kept current every step (see 5–6). First-minimum
   tie-breaking as everywhere.
2. **Swap** `i ↔ j*` in: stored `A` columns (full height), stored prefix `W₀` rows, `s`,
   `s_ref`, `wnorm2`, `jpvt`, and the panel bookkeeping indexed by absolute column position —
   `Fa` rows, `B` columns, `Bsel` columns (entries of the `t-1` prior steps). `V` is indexed
   by step and never swaps.
3. **Catch up the pivot column and reflect.** Only the still-stale rows `i:m` carry deferred
   updates (rows `s..i-1` of column `i` were finalized by earlier per-step row writes):
   `a_i(i:m) -= V(i:m, 1:t-1)·Fa(i, 1:t-1)ᵀ`, then `(τ_t, ρ_t) =` `dlarfg` on `a_i(i:m)`;
   store `v_t` into `V` (zeros above `i`, unit at `i`, tail from `A`).
4. **`Fa` column** (`dlaqps` recurrence, restricted to rows `i:m` — exact because `v_t` is
   zero above row `i`):
   `Fa(rem, t) = τ_t·(A_stored(i:m, rem)ᵀ v_t − Fa(rem, 1:t-1)·(V(i:m, 1:t-1)ᵀ v_t))`
   — one read-only gemv over the stored trailing block plus `O((M+N)·nb)` corrections.
5. **Finalize row i of the trailing block** (needed for β and the norm downdates):
   `α(rem) = A_stored(i, rem) − V(i, 1:t)·Fa(rem, 1:t)ᵀ`, written back. Then `β = α/ρ_t` and
   the usual clamped downdate `s_j -= α_j²`.
6. **W-side scalars from the deferred representation.** Catch up the single pivot column of
   the prefix: `ω_t = W₀(1:s-1, i) − Ω(:,1:t-1)·Bsel(1:t-1, i)` (store as `Ω(:,t)`), and read
   its panel part `q_t = B(1:t-1, i)` (already true). With `w*_t = [ω_t; q_t]`, the criterion
   dot products for all remaining `j` expand to
   `d_j = w_jᵀw*_t = W₀(:,j)ᵀω_t − Bsel(1:t-1, j)ᵀ(Ω(:,1:t-1)ᵀω_t) + B(1:t-1, j)ᵀq_t`
   — one read-only gemv `g = W₀ᵀω_t` over the prefix block (the shared small vector
   `Ωᵀω_t` is computed once) plus `O(nb·N)` corrections. The running-norm recurrence is then
   the unblocked one verbatim: `wnorm2_j ← wnorm2_j − 2β_j d_j + β_j²(1 + wnorm2_pivot)`,
   using the pivot's *running* `wnorm2` exactly as the unblocked kernel does.
7. **Eager panel-row updates of W**: record `Bsel(t, rem) = β`; apply
   `B(1:t-1, rem) -= q_t·βᵀ`; append `B(t, rem) = β`.

### 3.1 The O(1) running-norm recurrence: provenance and its role here

The Householder reduction itself costs `Θ(mnk)` flops regardless of pivoting policy
(`Θ(nk²)` in the manuscript's `k×n` GKS setting, where `m = k`). On top of that, the criterion
needs `||w_j||²` for every remaining column at every step, and the cost of *maintaining the
selection quantities* has a ladder worth making explicit in the appendix:

- **Naive**: recompute `w_j = R11^{-1}r_j` by triangular solve — `O(i²)` per column per
  step, `Θ(nk³)` for the maintenance over the whole factorization, which would dominate the
  `Θ(mnk)` reduction (the manuscript's Appendix A already rejects this).
- **Stewart (1990), Fig. 2.2**: maintain `S = R11^{-1}R12` incrementally (`O(i)` per column
  per step, `Θ(nk²)` total — at most the reduction's own cost, and equal to it when
  `m = k`), then recompute `σ_j = ||s_j||` *directly from the updated column* each step —
  asymptotically free (another `O(i)` pass per column) but operationally a second full sweep
  over the prefix block per step.
- **Bischof (1990)** sidesteps exact maintenance entirely: his incremental condition
  estimator carries an approximate singular vector updated by a 2×2 eigenproblem per step —
  `O(n)` per step, but it yields an *estimate* of `σ_min`, not the exact Frobenius quantity
  the selection rule here minimizes.
- **The manuscript's eq. (A.2)**, which both kernels implement: expand
  `||s_j − β_j s_k||² + β_j²` once and update
  `wnorm2_j ← wnorm2_j − 2β_j d_j + β_j²(1 + ||w*||²)` with `d_j = w_jᵀw*` — `O(1)` per
  column given the dot products, so the norm tracking rides along with the same `Θ(nk²)`
  S-maintenance term instead of adding a pass. This recurrence is **not in Stewart**; it is
  a derived, exact-arithmetic-equivalent refinement of his direct recomputation (recorded as
  deviation #9 in `docs/VALIDATION.md`, pinned by the oracle tests; the V2 stress analysis
  bounds its rounding drift, with the structural protection that the selection rule never
  picks a large-`||w*||` pivot).

Every exact variant therefore runs at `Θ(mnk + nk²) = Θ(mnk)` overall — the policies differ
in the constant and in memory passes, not in the exponent; only the naive policy breaks the
order.

In the *unblocked* kernel the refinement is roughly cost-neutral: the dots are one gemv over
the W prefix, replacing Stewart's strided column-norm sweep over the same block — one
optimized BLAS-2 pass either way. Its real payoff is that it is **what makes the panel kernel
possible**. The recurrence needs only the scalars `d_j`, and step 6 above shows `d_j` is
computable from the *deferred* representation (`W₀`, `Ω`, `Bsel`, `B`) at the cost of one
read-only gemv plus `O(nb·N)` corrections — the updated W columns never need to exist
mid-panel. Stewart's direct `σ_j = ||s_j||` recomputation, by contrast, requires the updated
column itself, which would force a flush at every step and defeat the deferral entirely. The
O(1) recurrence is therefore the hinge that decouples "selection scalars current" (every
step) from "W block current" (once per panel) — without it, the W side of BSQR has no blocked
formulation, and the short-wide regime would be stuck at the unblocked flop floor measured in
P0.

## 4. Panel flush (after nb steps, or fewer at `rank_stop`/end)

- `A(ie+1:m, rem) -= V(ie+1:m, 1:tdone)·Fa(rem, 1:tdone)ᵀ` — one gemm (rows `s..ie` were
  finalized per-step).
- `W₀(1:s-1, rem) -= Ω(:,1:tdone)·Bsel(1:tdone, rem)` — one gemm.
- Copy `B` rows into `W` rows `s..ie` (they join the prefix for the next panel).
- Exact-norm refresh of the flagged columns (next section).

## 5. The norm-refresh policy, and why the delay is safe

The recompute guard (`s_j ≤ s_ref_j·tol`, `tol = √eps`) needs the *true* tail of column `j`,
which is stale mid-panel. Two policies were implemented and measured:

1. *Early flush* (`dlaqps`'s `lsticc` mechanism): truncate the panel at the trip, flush,
   refresh, restart. Refresh timing is step-identical to the unblocked kernel, but
   recompute-heavy regimes trip nearly every step and panels degenerate to width ~1: the full
   publication run showed ill-conditioned short-wide **regressing** 1.72 → 1.97. (This also
   explains the P0 observation that `dgeqp3 ≈ dgeqpf` on short-wide: LAPACK's own blocking
   degenerates the same way there.)
2. *Batched refresh* (**adopted**): flag each tripping column once per panel and refresh all
   flagged columns exactly at the natural flush. Ill-conditioned short-wide improves to 1.41
   under this policy.

**Safety analysis.** The criterion is `(1+||w_j||²)/s_j`; the delay affects only the
denominator of flagged columns — the `wnorm2` recurrence runs for every column at every step
in both kernels, so the `R11^{-1}R12`-dependent numerator is never stale. For the denominator:

- At flag time `s_down ≤ tol·s_ref` with absolute error `≲ c(n)·eps·s_ref` — the trigger
  `tol = √eps` is calibrated so the running value still has ~`√eps` relative accuracy there.
- After the flag, both kernels subtract the same `α_j²` sequence (up to per-step `eps`
  rounding), so the gap between the panel's stale value and the unblocked kernel's refreshed
  value stays **frozen at the downdate-noise floor** `~c(n)·eps·s_ref`; it does not amplify
  during the window. Only its relative size grows as the true value shrinks.
- A selection can therefore deviate from the unblocked kernel's only when a competitor's
  criterion lies within a band of relative width `O(c·n·eps·s_ref/s_true)` of the flagged
  column's — near-ties at the noise scale — or when `s_true ≲ c·n·eps·s_ref` outright (true
  tail at the `√(n·eps)·||a_j||` scale), where pivot order is rounding-determined for *every*
  variant of the algorithm and below any rank-revealing resolution. Backward stability of the
  factorization is independent of pivot choice throughout; only selection quality among
  numerically indistinguishable candidates is at stake, and the manuscript's guarantees
  presuppose full-rank progress through the `k` steps in any case.
- The exposure window is ≤ `nb−1` steps, after which the exact refresh restores truth;
  `rank_stop` covers genuine rank-deficient continuation.

Empirically: identical pivot sequences on the full tie-free fixture suite including the
recompute-heavy members (hundreds of trips per factorization), criterion traces within
`rtol_crit`, and unchanged `||R11^{-1}||_F` on the cancellation-stress constructions.

A bookkeeping consequence: the per-step `nrecomp` trace counts one flag per column per panel,
unlike the unblocked kernel's per-trip counts; the cross-backend trace-parity test pins
`BS_PANEL_NB=0` accordingly.

## 6. Equivalence and cost accounting

*Exact arithmetic:* every quantity entering selection (`s_j`, `wnorm2_j`), every reflector,
and the final `R`/`W` are the unblocked kernel's values with sums regrouped; there is no
approximation anywhere. *Rounding:* differs, so the change is algorithm-visible under the
project's change taxonomy — Julia and MEX moved together, and pivot parity was revalidated on
the oracle fixtures in both languages and both dispatch modes.

Per step, big-block traffic drops from four passes (gemv+ger on the trailing block, gemv+ger
on the W prefix) to two read-only gemv passes; the two ger streams become two rank-`nb` gemms
per panel (amortized ~`1/nb` of a pass each, at gemm throughput). Correction terms cost
`O(nb·(M+N))` per step; the per-panel extras (pivot-column catch-ups, `Ω` maintenance) are
`O(nb²·(M+s))`. Workspace overhead beyond the unblocked kernel: `(M + 3N)·nb + (s-1)·nb`
doubles plus two `nb` and two `N` scratch vectors.

## 7. Parameters

| parameter | default | provenance |
|---|---|---|
| panel width `nb` (`BS_PANEL_NB`; `0`/`1` = unblocked kernel) | 8 | sweep below; best at every measured size in both languages |
| size crossover (`BS_PANEL_MIN_KN`), panel used iff `k·n ≥` this | 24576 | measured boundary below; panel bookkeeping dominates µs-scale factorizations |
| `norm_recomp_tol` | `√eps` (unchanged) | safeguard measured cost-free (≤3%, P0.4); also calibrates the safety analysis of §5 |

Panel-width sweep (Julia, fair materialized `Q,R,p` path vs `dgeqp3`, gaussian, 1 BLAS
thread, post-P1.1):

| case | unblocked | nb=4 | nb=8 | nb=16 | nb=32 |
|---|---:|---:|---:|---:|---:|
| square 256×256 | 1.547 | 1.283 | **1.146** | 1.221 | 1.348 |
| square 384×384 | 1.847 | 1.330 | **1.187** | 1.250 | 1.351 |
| square 512×512 | 1.932 | 1.642 | **1.290** | 1.299 | 1.386 |
| short_wide 64×640 | 1.830 | 1.725 | **1.549** | 1.639 | 1.870 |
| short_wide 128×1024 | 1.875 | 1.608 | **1.442** | 1.521 | 1.683 |
| short_wide 256×1024 | 2.031 | 1.728 | **1.446** | 1.462 | 1.599 |

Crossover measurement (panel time / unblocked time; gaussian, 1 thread): losers at
`k·n ≤ 16384` — 32×32 (2.07), 64×64 (1.58), 96×96 (1.24), 128×128 (1.12), 32×256 (1.57),
64×128 (1.28); winners at `k·n ≥ 32768` — 192×192 (0.85), 256×256 (0.74), 64×640 (0.85),
128×256 (0.86), 128×1024 (0.80). The default 24576 sits in the clean separation between.

## 8. Implementation notes and considerations

- **Row-range discipline is the correctness crux.** Per-step row finalization means rows
  `s..i-1` of the trailing block are *current* while rows `i:m` are *stale*; every deferred
  operation must restrict to the stale range. The pivot-column catch-up over `s:m` instead of
  `i:m` double-applies updates to finalized rows (the one bug hit during bring-up — caught
  immediately by the harness as "pivots/criterion right, `R` wrong"). The `Fa` gemv restricts
  itself for free because `v_t` is zero above row `i`.
- **Swap bookkeeping.** All panel structures indexed by absolute column position (`Fa` rows,
  `B`/`Bsel` columns) must participate in pivot swaps; `V` (indexed by step) must not. The
  prefix-`W` swap covers rows `1..s-1` only — the panel rows live in `B`.
- **`Bsel` vs `B`.** The deferred prefix update needs the β coefficients *as of selection
  time*, but `B` rows mutate (they receive later panel rank-1s) — hence the frozen copy.
- **Pivot `wnorm2` for the recurrence.** Use the running scalar (as the unblocked kernel
  does), not a recomputation from `ω_t`/`q_t` — keeps the recurrence self-consistent and the
  trace comparable.
- **Early termination.** `rank_stop` (and panel-truncating events generally) flush whatever
  `tdone` steps completed; flush handles `tdone < nb` uniformly.
- **Where the small-size overhead comes from.** Per-call panel-buffer allocation/zeroing and
  ~2 extra BLAS-call latencies per step dominate factorizations in the tens of µs — hence the
  `k·n` crossover rather than a `min(m,n)` one (64×128 loses while 64×640 wins at the same
  `k`).
- **Threading interaction.** The unblocked BLAS-2 kernel interacts badly with multithreaded
  BLAS (most visible under MATLAB's default threading: up to 3.8× vs `qr` at 384²); the panel
  gemms turn threads from a liability into an asset. This is why the MEX-side gains exceed
  the Julia 1-thread gains.

## 9. Validation methodology and outcomes

Every claim above is pinned by the Part-I harness rather than by inspection:

- **Oracle fixtures** (tie-screened inputs with expected outputs from the literal-transcription
  oracle): exact pivot-sequence equality, `R`/`rinv` tolerances, and per-step criterion traces
  within `rtol_crit` — run in both languages and in both dispatch modes (default and
  forced-panel), so both execution paths stay pinned.
- **Direct panel-vs-unblocked comparisons** on nine constructions (square/tall/short-wide/
  graded/scaled/early-stop/`nb>k`/`nb=2`): identical pivots, `R` agreement at 1e-16.
- **Cancellation-stress matrices** (graded spectra + range-orthogonal near-duplicates driving
  `||R11^{-1}||_F` to ~1e10): selection quality ratio 1.0, `rinv` at eps.
- A cautionary negative result worth keeping: comparing variants at `k` beyond the numerical
  rank flips pivots among noise-level columns — by design, not by bug; validation excludes
  that regime via gap screening (and it is where `rank_stop` applies).

## 10. Key experimental outcomes (publication grid, committed artifacts)

Geomean relative time vs the built-in pivoted QR baseline (old = unblocked kernel, new =
panel default with crossover; both languages on Apple Silicon/Accelerate):

| family/regime | Julia old → new | MATLAB old → new |
|---|---|---|
| gaussian square | 1.60 → **1.32** | 1.65 → **1.30** |
| gaussian short-wide | 1.83 → **1.57** | 2.19 → **1.67** |
| ill-conditioned square | 1.53 → **1.22** | 1.76 → **1.27** |
| ill-conditioned short-wide | 1.72 → **1.41** | 1.90 → **1.53** |
| orthonormal-rows square | 1.63 → **1.45** | 1.89 → **1.34** |
| orthonormal-rows short-wide | 1.86 → **1.55** | 2.08 → **1.67** |
| `rinv` variants | down to **1.30** | down to **1.40** |

Perf gates vs the pre-panel baselines: MATLAB 0/300 violations; Julia passes with a 1 ms
noise floor (all flagged cells were same-code noise at m ≤ 64; median change on gated cells
is a 14–15% *speedup*). Square sizes ≥256 run *below* the unblocked 1.25 flop-floor
prediction — the deferred gemms exceed `dgeqp3`'s average throughput at those sizes — and the
remaining short-wide deficit is the irreducible W-maintenance flops, now partially at gemm
throughput.

Context for the appendix narrative: P0's floor analysis showed short-wide BSQR already at its
flop floor *for an unblocked kernel*; the panel reorganization is what moved the floor itself.
