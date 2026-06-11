# BSQR Validation Audit

Companion to `docs/VALIDATION_AND_PERF_PLAN.md` (Part I). This file is the living record of
(a) how the code maps onto the references, (b) every place an implementation deviates from the
literal Algorithm 1 of `notes/bischof_stewart_pivoting.tex`, and (c) which test pins each claim.

Ground truth: Bischof (1990), Stewart (1990) — both in `notes/` — distilled as Algorithm 1 in
`notes/bischof_stewart_pivoting.tex` ("the writeup"). The executable ground truth is
`matlab/tests/oracle_bsqr.m` (V1 oracle): a literal transcription that recomputes
`w_j = R11^{-1} R12(:,j)` by triangular solve and tail norms from scratch each step, sharing no
recurrence with the production kernels.

## Symbol map (writeup ↔ code)

| Writeup | Julia (`julia/src/`) | MATLAB m-file | MEX (`bsqr_mex.cpp`) |
|---|---|---|---|
| `W(1:i-1, j)` (`w_j = R11^{-1} R12(:,j)`) | `ws.W` | `W` | `ws.W` |
| `‖w_j‖²` (running) | `ws.wnorm2` | `wnorm2` | `wnorm2` |
| `s(j)` (squared tail norms) | `ws.s` / `ws.s_ref` | `s` / `s_ref` | `s` / `s_ref` |
| `ρ_i` (new diagonal) | `beta_i` from `_householder!` | `beta_i` | `beta` from `dlarfg` |
| `β_j = α_j/ρ_i` | `ws.beta` (`beta_vec`) | `beta_vec` implicit | `beta_vec` |
| `w_ℓᵀ w*` (for ‖w‖² update) | `ws.dots` via `gemv` | inside `update_trailing_state` | `dots` via `dgemv` |
| criterion `c_j = (1+‖w_j‖²)/s_j` | `_select_pivot_column!` | `select_pivot_column` | `select_pivot_column` |
| `π` | `jpvt` | `p` | `p` |

## V0 paper audit (completed 2026-06-11)

Sources: Bischof, *Incremental Condition Estimation*, SIAM J. Matrix Anal. Appl. 11(2), 1990
(`notes/Bischof.pdf`); Stewart, *Incremental Condition Calculation and Column Selection*,
UMIACS TR, 1990 (`notes/Stewart.pdf`).

**Attribution.** The algorithm we implement is Stewart's. Bischof's ICE maintains an
*approximate* singular vector to estimate `σ_min` of a growing triangular matrix (his eqs.
(1)–(10): a 2×2 eigenproblem per step, motivated as a secular-equation approximation), and his
§5 only *guards* classical column pivoting with that estimate — he does not propose a pivot
criterion. Stewart's condition calculator computes `‖R11⁻¹‖_F` *exactly* via the `S` matrix and
turns it into the selection rule we use. Writeup-wording nuance: the writeup's introduction says
Bischof "tracks the growth of `‖R⁻¹‖_F`"; precisely, Bischof tracks an estimate of `σ_min`,
and the exact Frobenius tracking is Stewart's contribution.

**Equation map** (Stewart → writeup → code):

| Stewart (1990) | Writeup | Code |
|---|---|---|
| §2: `S ≡ (s_{k+1},…,s_n) = R11⁻¹R12` | eq. (w-def), `W` | `ws.W` / `W` |
| eq. (2.2): `ν̃ = √(ν² + α⁻²(1+‖s‖²))` | eq. (frob-growth) | `frob_inv_trace` accumulation (Julia), `crit_best` (oracle) |
| eq. (2.3): `s_j ← (s_j − α⁻¹ρ_{k+1,j} s_{k+1};  α⁻¹ρ_{k+1,j})` | eq. (w-update), `β_j = α_j/ρ` | `ger`/`dger` rank-1 update + β row |
| Fig. 2.2: minimize `ω_j⁻²(1+σ_j²)` | eq. (criterion) | `_select_pivot_column!` / `select_pivot_column` |
| Fig. 2.2 init: `σ_j = 0`, `ω_j` = column norms | Remark "First step simplification" (= Businger–Golub at step 1) | `wnorm2 = 0`, `s` init |
| Fig. 2.2: `ω_j ← √(ω_j² − ρ_{k+1,j}²)` | line 19 downdate | `s[j] -= α²` |
| §2 end: Householder variant, `ω_i` = norms of `A22` columns | Algorithm 1 (the variant we implement) | tail norms `s` |
| Footnote 2: care in `ω` updates, cites LINPACK `sqrdc`; recompute `S` column *ab initio* on cancellation | Remark "Column norm downdating" (ω part only) | recompute guard on `s` only — see open items |

**Confirmations.** The criterion, the `S`/`w` update, the `ω` downdate with safeguard, the
first-step reduction to Businger–Golub, and the Householder-variant tail-norm denominators all
match Stewart exactly. The running-`ν` identity (writeup Remark "Tracking the full inverse
norm") is Stewart's (2.2) and is pinned by `testCriterionTraceMatchesInverseFrobNorm`.

**Differences found (all benign, now documented):**

1. *`σ` maintenance:* Stewart's Fig. 2.2 recomputes `σ_j = ‖s_j‖` directly from the updated
   `s_j` each step. The writeup's eq. (wnorm-update) is an O(1)-per-column recurrence obtained
   by expanding `‖s_j − β_j s_{k+1}‖² + β_j²` — algebraically identical, different rounding.
   The kernels implement the recurrence; the oracle uses neither (from-scratch solves), so
   parity tests pin the equivalence. Added as deviation #9 below.
2. *Tie-breaking:* Stewart only implies a convention ("breaks ties by doing nothing", i.e. the
   current column wins a tie) in his Kahan discussion. Our strict-`<` first-minimum scan keeps
   the current column on ties and is consistent with, and a total refinement of, Stewart's.
3. *Deviation #1 provenance:* the norm-downdate safeguard is not merely "LAPACK-style" — it is
   prescribed by Stewart himself (footnote 2, citing LINPACK `sqrdc`).

**Intentionally out of scope (papers contain, we do not implement):**

- Stewart §3 recovery procedure (swap-out heuristic that fixes the Kahan failure). The greedy
  forward selection — Stewart's, the writeup's, and ours — fails on Kahan's matrix by design;
  Osinsky's guarantees cover the orthonormal-rows regime. The Kahan matrix belongs in the test
  zoo as a documented known-failure of greedy selection, not as a bug.
- Stewart's plane-rotation variant for selecting within a *precomputed* triangular `R`; we
  implement his Householder variant for a general `m×n` matrix.
- Bischof's ICE estimator itself and his §5 guarded-pivoting scheme.

**Open item from the audit — RESOLVED (no guard needed):** Stewart's footnote 2 prescribes
recomputing a column of `S` *ab initio* if cancellation occurs in its update. Our kernels
safeguard only the `ω`/`s` downdates; the `W` columns and the `wnorm2` recurrence (which
contains a subtraction, `−2β_j·(w_jᵀw*)`) have no analogous safeguard. Stress experiments
(2026-06-11) engineered maximal one-step cancellation — near-duplicate columns `c·B + T` with
unit tails `T ⊥ range(B)` over graded spectra driving `‖R11⁻¹‖_F` to ~1e10 — and measured:
`wnorm2` recurrence vs. tracked-`W` drift ~1e-15, `W` vs. from-scratch solve ~1e-15,
selected-pivot criterion drift vs. oracle ≤ ~1e-7, zero pivot flips, selection-quality ratio
exactly 1, `rinv_r12` at eps level. The structural reason the failure mode cannot fire: a
catastrophic one-step collapse of a large `wnorm2` requires the *pivot's* `‖w*‖` to be large,
but the criterion `min (1+‖w‖²)/ρ²` never selects such a column, so the subtracted term stays
moderate. The selection rule itself provides the protection Stewart's remedy targets. Pinned by
`matlab/tests/test_cancellation_stress.m` (both backends: criterion drift, quality ratio,
`rinv_r12` accuracy) and the Julia testset "Cancellation stress: wnorm2/W recurrence drift"
(internal three-layer drift probe via the workspace).

## Deviation table

Sanctioned deviations from the literal Algorithm 1. "Pinned by" names the test that fails if the
equivalence claim is wrong.

| # | Deviation | Where | Why mathematically equivalent | Pinned by |
|---|---|---|---|---|
| 1 | Norm downdate `s_j -= α_j²` with recompute guard `s_j ≤ s_ref_j·tol` instead of from-scratch norms | all three kernels | recomputation returns the same mathematical quantity; guard is the writeup's Remark "Column norm downdating" | `test_oracle_parity` (oracle uses exact norms; pivot sequences must still match); Julia testset "Norm recompute tolerance knob" |
| 2 | Householder sign convention `ρ = -sign(x₁)‖x‖` (LAPACK), identity reflector when the tail is zero; writeup states `ρ = ‖x‖ ≥ 0` | all three kernels and the oracle (Julia and MEX both call `dlarfg` directly since V5) | row sign flips cancel in `(D·R11)⁻¹(D·R12)`, so every `w_j` and the criterion are unchanged; `Q` absorbs `D` (proof sketch in `oracle_bsqr.m` header) | direct (sign-normalization-free) Q/R comparison in `test_oracle_parity` |
| 3 | Clamp `wnorm2` and `s` at 0 after downdates | all three kernels | guards `-eps`-level negatives of nonnegative quantities | indirectly by parity; dedicated near-zero unit test pending |
| 4 | `s(i) = β_i²` refresh of the pivot column post-reflection | all three kernels | exact identity after the reflection | `test_oracle_parity` |
| 5 | BLAS-2 organization of the W update (`gemv` for dots, `ger` for the rank-1) | Julia, MEX | writeup Remark "Practical considerations"; same recurrence, regrouped | `test_oracle_parity` (oracle has no recurrence at all) |
| 6 | `rank_stop` early-exit option | Julia only, off by default in `bsqr` | extension, not in Algorithm 1; cross-language defaults agree (no early exit) | Julia testset "Rank-stop policy"; cross-language default parity pending (V3 fixtures) |
| 7 | Short-wide fastpath (inline W-row materialization) | Julia only | pure memory-layout reordering; identical arithmetic | Julia testset "Kernel helper invariants and fastpath knobs" |
| 8 | First-minimum tie-breaking via strict `<` | all three kernels and the oracle | writeup's argmin leaves ties unspecified; consistent with Stewart's "ties by doing nothing" | `testPivotTieStability` (MATLAB); Julia "Criterion-consistent pivot sequence" |
| 9 | `‖w_j‖²` maintained by the writeup's O(1) recurrence instead of Stewart's direct `σ_j = ‖s_j‖` recomputation | all three kernels | algebraic expansion of `‖s_j − β_j s_k‖² + β_j²`; identical in exact arithmetic | `test_oracle_parity` (oracle computes `‖w_j‖` from scratch) |

## Validation assets

| Asset | Status |
|---|---|
| V1 oracle (`matlab/tests/oracle_bsqr.m`) | done |
| Oracle/backend parity + invariants + Osinsky canaries (`matlab/tests/test_oracle_parity.m`) | done (mfile + MEX) |
| V0 paper → writeup audit | done (see above) |
| V3 cross-language fixtures (`parity/`, generated by `matlab/tests/generate_parity_fixtures.m`; consumed by `matlab/tests/test_parity_fixtures.m` and `julia/test/test_parity_fixtures.jl`) | done |
| V4 per-step trace instrumentation: `trace` option (crit + nrecomp per step) in mfile and MEX, `recomp_history`/`frob_inv_trace` hooks in the Julia kernel; oracle crit traces in `parity/` compared by all three implementations; mfile↔MEX `nrecomp` equality pinned | done |
| V5 Householder unification: Julia `_householder!` is a direct `dlarfg` ccall, sharing LAPACK's construction (and rescaling safeguards) with the MEX | done |
| V2 follow-up: `wnorm2`/W cancellation stress test (V0 open item) | done — no guard needed; see resolution above |

Notes on test design choices:

- Exact pivot-sequence equality across implementations is demanded only on inputs screened for
  criterion near-ties (`testZooIsTieFree` asserts every step's runner-up gap exceeds `1e-6`
  relative). With fixed seeds this is deterministic; a failing member gets a new seed, never a
  looser tolerance.
- Forward agreement of `Q`, `rinv_r12`, and `inv(R11)` is condition-sensitive even with
  identical pivots, so graded-spectrum zoo members check pivot sequence + `R` agreement only;
  their `Q`/solve quality is covered by backward-style contract tests instead.
