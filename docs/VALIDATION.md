# BSQR Validation Audit

Companion to `docs/VALIDATION_AND_PERF_PLAN.md` (Part I). This file is the living record of
(a) how the code maps onto the references, (b) every place an implementation deviates from the
literal Algorithm A.1 of the manuscript, and (c) which test pins each claim.

Ground truth: Bischof (1990), Stewart (1990) — both in `notes/` — distilled as Algorithm A.1
in Appendix A of the manuscript (`notes/GKSevolved_draft.tex`, compiled as
`notes/manuscript_draft.pdf`; local-only until the final version is added to the repository).
"The writeup" below refers to that appendix. The executable ground truth is
`matlab/tests/oracle_bsqr.m` (V1 oracle): a literal transcription that recomputes
`w_j = R11^{-1} R12(:,j)` by triangular solve and tail norms from scratch each step, sharing no
recurrence with the production kernels.

## Symbol map (manuscript ↔ code)

| Manuscript (Alg. A.1) | Julia (`julia/src/`) | MATLAB m-file | MEX (`bsqr_mex.cpp`) |
|---|---|---|---|
| `w_j` (`= R11^{-1} R12(:,j)`, maintained by line 9) | `ws.W` | `W` | `ws.W` |
| `‖w_j‖²` (running) | `ws.wnorm2` | `wnorm2` | `wnorm2` |
| `ρ_j²` (squared tail norms) | `ws.s` / `ws.s_ref` | `s` / `s_ref` | `s` / `s_ref` |
| `R(i,i)` (new diagonal after the reflection) | `beta_i` from `_householder!` | `beta_i` | `beta` from `dlarfg` |
| `α_j = R(i,j)/R(i,i)` | `ws.beta` (`beta_vec`) | `beta_vec` implicit | `beta_vec` |
| `w_jᵀw_i` (for the eq. (A.2) ‖w‖² update) | `ws.dots` via `gemv` | inside `update_trailing_state` | `dots` via `dgemv` |
| criterion `(1+‖w_j‖²)/ρ_j²` (line 4) | `_select_pivot_column!` | `select_pivot_column` | `select_pivot_column` |
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
turns it into the selection rule we use. The manuscript's attribution is precise (Appendix A
credits the selection rule to Stewart); Bischof tracks an estimate of `σ_min`, and the exact
Frobenius tracking is Stewart's contribution.

**Equation map** (Stewart → manuscript → code):

| Stewart (1990) | Manuscript | Code |
|---|---|---|
| §2: `S ≡ (s_{k+1},…,s_n) = R11⁻¹R12` | `w_j` state, Alg. A.1 line 9 | `ws.W` / `W` |
| eq. (2.2): `ν̃ = √(ν² + α⁻²(1+‖s‖²))` | eq. (2.4) | `frob_inv_trace` accumulation (Julia), `crit_best` (oracle) |
| eq. (2.3): `s_j ← (s_j − α⁻¹ρ_{k+1,j} s_{k+1};  α⁻¹ρ_{k+1,j})` | Alg. A.1 line 9 (`α_j = R(i,j)/R(i,i)`) | `ger`/`dger` rank-1 update + β row |
| Fig. 2.2: minimize `ω_j⁻²(1+σ_j²)` | eq. (2.3) / Alg. A.1 line 4 | `_select_pivot_column!` / `select_pivot_column` |
| Fig. 2.2 init: `σ_j = 0`, `ω_j` = column norms | the `‖w_j‖ = 0` convention on Alg. A.1 line 4 (= Businger–Golub at step 1) | `wnorm2 = 0`, `s` init |
| Fig. 2.2: `ω_j ← √(ω_j² − ρ_{k+1,j}²)` | Alg. A.1 line 8 downdate | `s[j] -= α²` |
| §2 end: Householder variant, `ω_i` = norms of `A22` columns | Algorithm A.1 (stated directly for general `m×n`) | tail norms `s` |
| Footnote 2: care in `ω` updates, cites LINPACK `sqrdc`; recompute `S` column *ab initio* on cancellation | Remark 7 (safeguarded recomputation; `ρ²` part) | recompute guard on `s` only — see open items |

**Confirmations.** The criterion, the `S`/`w` update, the `ω` downdate with safeguard, the
first-step reduction to Businger–Golub, and the Householder-variant tail-norm denominators all
match Stewart exactly. The running-`ν` identity (manuscript eq. (2.4)) is Stewart's (2.2) and
is pinned by `testCriterionTraceMatchesInverseFrobNorm`.

**Differences found (all benign, now documented):**

1. *`σ` maintenance:* Stewart's Fig. 2.2 recomputes `σ_j = ‖s_j‖` directly from the updated
   `s_j` each step. The manuscript's eq. (A.2) is an O(1)-per-column recurrence obtained
   by expanding `‖s_j − β_j s_{k+1}‖² + β_j²` — algebraically identical, different rounding.
   The kernels implement the recurrence; the oracle uses neither (from-scratch solves), so
   parity tests pin the equivalence. Added as deviation #9 below.
2. *Tie-breaking:* Stewart only implies a convention ("breaks ties by doing nothing", i.e. the
   current column wins a tie) in his Kahan discussion. Our strict-`<` first-minimum scan keeps
   the current column on ties and is consistent with, and a total refinement of, Stewart's.
3. *Deviation #1 provenance:* the norm-downdate safeguard is not merely "LAPACK-style" — it is
   prescribed by Stewart himself (footnote 2, citing LINPACK `sqrdc`) and now by the
   manuscript's Remark 7.

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

## Re-audit (2026-07-05, publication readiness Phase 1)

Re-diffed all four kernel paths — Julia unblocked (`_bsqr_kernel!`), Julia panel
(`_bsqr_kernel_panel!`), MATLAB m-file, and the MEX unblocked + panel paths — against
`oracle_bsqr.m` and the writeup. No undocumented kernel/oracle mismatch was found. Changes to
this file from the re-audit: added deviations #10–#11 (the panel kernel's deferred/regrouped
updates and its batched norm-recompute refresh, previously documented only in
`docs/P3_BLOCKED_BSQR.md`) and #12 (the zero-tail-norm/rank-deficiency extension shared with
the oracle's guard); widened #5's scope to include the m-file; replaced the stale "pending"
note in #6 (the V3 fixtures landed); added the dedicated Phase 1 tie tests to #8's pins; and
scoped #7 to the unblocked kernel.

## Deviation table

Sanctioned deviations from the literal Algorithm 1. "Pinned by" names the test that fails if the
equivalence claim is wrong.

"All three kernels" means the Julia, MATLAB m-file, and MEX implementations. Since P3
(2026-06-11) Julia and the MEX also have a panel/blocked kernel
(`julia/src/kernel_panel.jl`, the panel path in `bsqr_mex.cpp`), which is the **default**
dispatch when `k·n ≥ 24576` (`BS_PANEL_NB=0` forces the unblocked kernel). The panel kernel
inherits deviations #1–#5, #8, #9 unchanged and adds #10–#11; its state symbols (`V`, `Fa`,
`B`, `Bsel`, `Ω`) are mapped in `docs/P3_BLOCKED_BSQR.md` §2. The manuscript states
Algorithm A.1 directly for a general `m×n` matrix with a step-count `k` (Stewart's Householder
variant, per the V0 audit above); the §3 selection guarantees are stated for the
orthonormal-rows (GKS) case.

| # | Deviation | Where | Why mathematically equivalent | Pinned by |
|---|---|---|---|---|
| 1 | Norm downdate `s_j -= α_j²` with recompute guard `s_j ≤ s_ref_j·tol` instead of from-scratch norms | all three kernels | recomputation returns the same mathematical quantity; the guard is prescribed by the manuscript's Remark 7 | `test_oracle_parity` (oracle uses exact norms; pivot sequences must still match); Julia testset "Norm recompute tolerance knob" |
| 2 | Householder sign convention `ρ = -sign(x₁)‖x‖` (LAPACK), identity reflector when the tail is zero; Alg. A.1 leaves the sign free (`±‖x‖e₁`), so this row pins the shared choice rather than a deviation | all three kernels and the oracle (Julia and MEX both call `dlarfg` directly since V5; the m-file and the oracle implement the same convention in plain MATLAB — hypot-based, without `dlarfg`'s subnormal rescaling loop, which is inert within the documented column-norm domain) | row sign flips cancel in `(D·R11)⁻¹(D·R12)`, so every `w_j` and the criterion are unchanged; `Q` absorbs `D` (proof sketch in `oracle_bsqr.m` header) | direct (sign-normalization-free) Q/R comparison in `test_oracle_parity` |
| 3 | Clamp `wnorm2` and `s` at 0 after downdates | all three kernels | guards `-eps`-level negatives of nonnegative quantities | indirectly by parity; dedicated near-zero unit test pending |
| 4 | `s(i) = β_i²` refresh of the pivot column post-reflection | all three kernels | exact identity after the reflection | `test_oracle_parity` |
| 5 | BLAS-2 organization of the W update: all dots `w_jᵀw*` computed from the pre-update `W` prefix (`gemv`), then the rank-1 update (`ger`) | all three kernels (the m-file via the equivalent vectorized matrix ops) | regrouping of Alg. A.1 line 9's per-column updates (Stewart (1990) §2); same recurrence | `test_oracle_parity` (oracle has no recurrence at all) |
| 6 | `rank_stop` early-exit option | Julia only, off by default in `bsqr` | extension, not in Algorithm 1; cross-language defaults agree (no early exit) | Julia testset "Rank-stop policy"; default-behavior parity pinned by the V3 fixture suite (both languages run at defaults) |
| 7 | Short-wide fastpath (inline W-row materialization) | Julia unblocked kernel only (above the panel crossover it is reached only with `BS_PANEL_NB=0`) | pure memory-layout reordering; identical arithmetic | Julia testset "Kernel helper invariants and fastpath knobs" |
| 8 | First-minimum tie-breaking via strict `<` | all three kernels and the oracle | Alg. A.1 allows any tie rule ("any reasonable choice is fine, e.g., picking the smaller index" — exactly our first-minimum scan); consistent with Stewart's "ties by doing nothing" | dedicated exact-tie tests: `testTieBreakFirstMinimum` (MATLAB, both backends) and Julia "Exact criterion ties: first minimum wins"; also `testPivotTieStability` (MATLAB), Julia "Criterion-consistent pivot sequence" |
| 9 | `‖w_j‖²` maintained by the manuscript's eq. (A.2) recurrence instead of Stewart's direct `σ_j = ‖s_j‖` recomputation | all three kernels | algebraic expansion of `‖s_j − β_j s_k‖² + β_j²`; identical in exact arithmetic (the manuscript notes the formula is optional for the asymptotic cost) | `test_oracle_parity` (oracle computes `‖w_j‖` from scratch) |
| 10 | Panel/blocked execution (dlaqps-style): within a width-`nb` panel the trailing block of `A` and the prefix rows of `W` receive no per-step writes — Householder applications accumulate in `V`/`Fa` and W-prefix rank-1s in `Ω`/`Bsel` (panel rows of `W` eagerly in `B`), flushed as rank-`nb` `gemm`s at panel end; criterion dots come from the deferred representation, `d = W₀ᵀω − Bselᵀ(Ωᵀω) + Bᵀq` | Julia + MEX panel kernels (the default path; see preamble) | regrouped sums only: every selection quantity, every reflector, and the final `R`/`W` equal the unblocked kernel's in exact arithmetic (derivation in `docs/P3_BLOCKED_BSQR.md` §§2–3, 6) | fixture suite in forced-panel mode: Julia `test_parity_fixtures.jl` ("forced panel" mode) and MATLAB `testMexPanelKernelMatchesFixtures` |
| 11 | Batched norm-recompute refresh: a guard trip mid-panel flags the column (once per panel) and the exact refresh happens at panel flush, up to `nb−1` steps later, instead of immediately | Julia + MEX panel kernels | the one panel change that alters values selection sees in floating point: only the flagged column's criterion denominator is affected, its gap to the refreshed value stays frozen at the downdate-noise floor, and a flagged column (`s ≤ tol·s_ref`) cannot win a selection inside the window — safety analysis in `docs/P3_BLOCKED_BSQR.md` §5. Trace consequence: per-step `nrecomp` counts one flag per column per panel, not per-trip | forced-panel fixture parity (incl. the recompute-heavy zoo members); the mfile↔MEX `nrecomp` trace-equality test pins `BS_PANEL_NB=0` for this reason |
| 12 | Zero-tail-norm handling: a column with running `s_j = 0` gets criterion `c_j = Inf` (ineligible) and its norm state is frozen at 0 (downdates and recompute checks skipped); `β_j` is taken as 0 when `ρ_i = 0`; if every remaining column is zero the kernels continue with identity reflectors (`tau = 0`) | all kernels (the `Inf` guard is shared by the oracle) | extension: Algorithm A.1 requires `k ≤ rank(A)` and never divides by zero; past that point the oracle errors (`oracle_bsqr:RankDeficient`) and kernel behavior is an extension covered by contract tests, not by oracle parity | `testOracleRejectsExactRankDeficiency` (documents the oracle boundary); Julia "Rank-stop policy" (continuation past deficiency with `rank_stop = false`) |

## Validation assets

| Asset | Status |
|---|---|
| V1 oracle (`matlab/tests/oracle_bsqr.m`) | done |
| Oracle/backend parity + invariants + Osinsky canaries (`matlab/tests/test_oracle_parity.m`) | done (mfile + MEX) |
| V0 paper → writeup audit | done (see above) |
| V3 cross-language fixtures (`parity/`, generated by `matlab/tests/generate_parity_fixtures.m`; consumed by `matlab/tests/test_parity_fixtures.m` and `julia/test/test_parity_fixtures.jl`) | done |
| V4 per-step trace instrumentation: `trace` option (crit + nrecomp per step) in mfile and MEX, `recomp_history`/`frob_inv_trace` hooks in the Julia kernel; oracle crit traces in `parity/` compared by all three implementations; mfile↔MEX `nrecomp` equality pinned (in unblocked mode, `BS_PANEL_NB=0` — deviation #11) | done |
| V5 Householder unification: Julia `_householder!` is a direct `dlarfg` ccall, sharing LAPACK's construction (and rescaling safeguards) with the MEX | done |
| V2 follow-up: `wnorm2`/W cancellation stress test (V0 open item) | done — no guard needed; see resolution above |
| P3 panel kernel parity: fixture suite run in default and forced-panel dispatch modes in both languages (deviations #10–#11) | done |

Notes on test design choices:

- Exact pivot-sequence equality across implementations is demanded only on inputs screened for
  criterion near-ties (`testZooIsTieFree` asserts every step's runner-up gap exceeds `1e-6`
  relative). With fixed seeds this is deterministic; a failing member gets a new seed, never a
  looser tolerance.
- Forward agreement of `Q`, `rinv_r12`, and `inv(R11)` is condition-sensitive even with
  identical pivots, so graded-spectrum zoo members check pivot sequence + `R` agreement only;
  their `Q`/solve quality is covered by backward-style contract tests instead.
