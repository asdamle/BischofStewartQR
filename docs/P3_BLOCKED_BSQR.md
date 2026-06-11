# P3: Panel/Blocked BSQR — Derivation

Phase P3 of `docs/VALIDATION_AND_PERF_PLAN.md`. Goal: reorganize the kernel so the two BLAS-2
walls measured in P0 (trailing reflector application, 44–57%; W update, 38–53%) become one
read-only gemv pass per step plus one rank-`nb` gemm per panel — the same restructuring
`dlaqps` applies to `dgeqp3`, extended to BSQR's W/wnorm2 state. Everything below is a
regrouping of the unblocked sums: all selection quantities are identical in exact arithmetic,
so this is an *algorithm-visible* change only through rounding (taxonomy: must land in Julia
and MEX together, parity suite green).

Notation follows Algorithm 1 of the writeup (1-indexed). A panel starts at global step `s` and
runs local steps `t = 1..nb` (global `i = s+t-1`). "Prefix" means rows `1..s-1`; `M = m-s+1`,
`N = n-s+1` are the trailing block dimensions at panel start.

## State and accumulators

Unblocked per-step state: `A` (reduced in place), `W` (rows `1..i-1`), running `s_j` (squared
tail norms), `wnorm2_j = ||W(1:i-1,j)||^2`.

Panel accumulators (allocated once in the workspace):

| symbol | shape | content |
|---|---|---|
| `V` | `M × nb` | panel reflector vectors `v_t` (unit lower-trapezoidal in panel rows) |
| `Fa` | `N × nb` | `dlaqps`-style deferred-update factors for A: the true trailing block is `A₀(s:m, s:n) − V(:,1:t)·Fa(:,1:t)ᵀ` |
| `B` | `nb × N` | the β rows appended to W during the panel, **eagerly maintained** (they receive the panel's later rank-1 updates) |
| `Ω` | `(s-1) × nb` | `Ω(:,t)` = prefix part (rows `1..s-1`) of the pivot's w-vector `w*_t` at its selection step |

Deferred W: the true prefix is `W₀(1:s-1, j) − Ω(:,1:t)·b̃(j,1:t)ᵀ`, where `b̃(j,u) = β_u[j]`
(the entries of `B` at selection time of step `u` — record them as `Bsel(u,j)` since `B` rows
mutate). The panel rows of W (`s..i-1`) are exactly the current rows of `B`.

Key structural point that makes this work without raggedness: the step-`u` rank-1 update
`W(1:i-1, rem) -= w*_u · β_uᵀ` splits into a *prefix* part (rows `1..s-1`, uniform height —
deferred via `Ω`/`Bsel`) and a *panel* part (rows `s..s+u-2`, at most `nb-1` rows — applied
eagerly to `B` at `O(nb·N)` per step).

## Per-step recipe (local step t, global step i = s+t-1)

1. **Pivot selection** — unchanged: `argmin (1+wnorm2_j)/s_j` over the running scalars. Both
   scalar arrays are kept current every step (see 5–6), so selection is identical to unblocked.
2. **Swap** columns `i ↔ j*` in: stored `A`, stored prefix `W₀`, `s`, `s_ref`, `wnorm2`, `jpvt`,
   plus the panel bookkeeping rows: columns of `Fa`, `B`, `Bsel` (indexed by trailing position).
3. **Catch up the pivot column of A**: `a_i := a_i − V(:,1:t-1)·Fa(p_i,1:t-1)ᵀ` — `O(M·nb)`.
   Build the reflector `(τ_t, ρ_t)` from the caught-up tail (`dlarfg`), store `v_t` in `V`.
4. **Fa column** (`dlaqps` recurrence): `Fa(:,t) = τ_t·(A₀(s:m,s:n)ᵀ v_t − Fa(:,1:t-1)·(V(:,1:t-1)ᵀ v_t))`
   — one *read-only* gemv over the stored trailing block plus `O((M+N)·nb)` corrections.
5. **Updated row i** (needed for β/downdates, again `dlaqps`-style):
   `α(rem) = A₀(i, rem) − V(i,1:t)·Fa(rem,1:t)ᵀ` — `O(nb·N)`. Then `β = α/ρ_t`;
   downdate `s_j -= α_j²` with the usual clamp.
6. **W-side scalars without touching the prefix block**, using the deferred form. Let
   `ω_t = Ω(:,t)` be the pivot's *true* prefix w (catch up the single column:
   `ω_t = W₀(1:s-1, p_i) − Ω(:,1:t-1)·Bsel(1:t-1, p_i)`, `O(s·nb)`), and let
   `q_t = B(1:t-1, p_i)` be its panel part (current `B` column — already true). Record
   `Bsel(t, :) = β`. Then for every remaining `j`:
   `d_j = w_jᵀ w*_t = (W₀(:,j))ᵀω_t − Bsel(1:t-1,j)ᵀ(Ω(:,1:t-1)ᵀω_t) + B(1:t-1,j)ᵀ q_t`
   — one *read-only* gemv `g = W₀ᵀ ω_t` over the prefix block plus `O(nb·N)` corrections.
   Update the running norms by the usual recurrence
   `wnorm2_j ← wnorm2_j − 2β_j d_j + β_j²(1 + ||w*_t||²)` with
   `||w*_t||² = ||ω_t||² + ||q_t||²` (equals the running `wnorm2` of the pivot).
7. **Eager panel-row updates of W**: `B(1:t-1, rem) -= q_t · βᵀ` (`O(nb·N)`), then append
   `B(t, rem) = β`.

## Panel flush (after nb steps, or early — see below)

- `A₀(i+1:m, rem) -= V(below,1:nb)·Fa(rem,1:nb)ᵀ` — one gemm. (Rows `s..i` of A are finalized
  per-step by 3/5 as in `dlaqps`.)
- `W₀(1:s-1, rem) -= Ω(:,1:nb)·Bsel(1:nb, rem)` — one gemm.
- Copy the `B` rows into `W` rows `s..s+nb-1`.

## Norm-recompute guard

The guard (`s_j ≤ s_ref_j·tol`) needs the *true* tail of column j, which is stale mid-panel.
Mirror LAPACK's `dlaqps`: when any downdate trips the guard, **end the panel early** at the
current step, flush the deferred updates, recompute the flagged norms exactly, and start a new
panel. This preserves the mathematical quantity (the recompute returns the same value either
way); it changes only rounding/timing of the refresh, which is within the sanctioned
deviation-#1 envelope. P0 measured 0 recomputes on square gaussian and ~600–900 on short-wide
(roughly one trip per step at n−i large) — if early termination degenerates panels to width ~1
in the short-wide regime, the guard policy needs refinement (e.g. batch the flagged columns and
recompute them against the *flushed* state at panel end, accepting one step of stale `s_ref`);
decide empirically in the prototype.

## Equivalence and cost

*Exact arithmetic*: every quantity entering selection (`s_j`, `wnorm2_j`), every reflector, and
the final `R`/`W` are the unblocked values with sums regrouped; no approximation anywhere.
*Rounding*: differs, so pivot sequences must be revalidated on the tie-free parity fixtures and
the crit traces must stay within `rtol_crit`.

Per-step big-block work drops from four passes (gemv+ger on A, gemv+ger on W-prefix) to two
read-only gemvs; the two ger streams become two rank-`nb` gemms per panel. Upper bound on the
kernel win if gemm time is negligible: ~2× on the 83–97% of kernel time the walls represent —
realistically (gemm is not free, corrections cost `O(nb·N)` per step) the target from the P0
model is: square toward its ~1.25 flop floor (from 1.55–1.99), short-wide toward
`(HH + W/2 + …)`-style ~1.4–1.5 (from 1.8–2.0). The P3 decision gate (≥20% median improvement
on the regimes that matter) is judged against these.

## Prototype status (2026-06-11): implemented, validated, **gate passed**

`julia/src/kernel_panel.jl` implements the derivation behind `BS_PANEL_NB` (off by default;
`bsqr!` dispatches when `nb ≥ 2`). Validation: the **full Julia suite passes with the panel
kernel enabled**, including exact pivot-sequence parity against the MATLAB-oracle fixtures, the
criterion-trace tolerances, and the cancellation-stress testset; a direct panel-vs-unblocked
comparison on nine constructions (square/tall/short-wide/graded/scaled/early-stop/nb>k/nb=2)
shows identical pivots and `R` agreement at 1e-16.

Benchmark (fair materialized path vs. `dgeqp3`, gaussian, 1 thread, post-P1.1):

| case | unblocked | nb=4 | nb=8 | nb=16 | nb=32 |
|---|---:|---:|---:|---:|---:|
| square 256x256 | 1.547 | 1.283 | **1.146** | 1.221 | 1.348 |
| square 384x384 | 1.847 | 1.330 | **1.187** | 1.250 | 1.351 |
| square 512x512 | 1.932 | 1.642 | **1.290** | 1.299 | 1.386 |
| short_wide 64x640 | 1.830 | 1.725 | **1.549** | 1.639 | 1.870 |
| short_wide 128x1024 | 1.875 | 1.608 | **1.442** | 1.521 | 1.683 |
| short_wide 256x1024 | 2.031 | 1.728 | **1.446** | 1.462 | 1.599 |

`nb = 8` wins everywhere: 26–36% on square (now *below* the 1.25 unblocked flop floor — the
deferred gemms run above `dgeqp3`'s average throughput at these sizes), 15–29% on short-wide
despite recompute-triggered early flushes. The ≥20% adoption gate is met with zero parity
regressions.

## Remaining to land (in order)

1. Port the identical panel algorithm to the MEX (algorithm-visible change: Julia and MEX move
   together per the change taxonomy; the m-file stays the unblocked reference).
2. Flip the default (`nb = 8`, env-overridable) in both, rerun the complete validation
   harness in both languages.
3. Regenerate publication benchmark artifacts and run the perf gates; revisit the
   `BS_PUB_THREADS` grid (gemm now benefits from threads).
4. Fold this derivation into the writeup as an appendix.
