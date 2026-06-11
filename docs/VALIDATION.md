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

Paper-level audit (Bischof 1990 / Stewart 1990 equations → writeup equations): **pending (V0)**.

## Deviation table

Sanctioned deviations from the literal Algorithm 1. "Pinned by" names the test that fails if the
equivalence claim is wrong.

| # | Deviation | Where | Why mathematically equivalent | Pinned by |
|---|---|---|---|---|
| 1 | Norm downdate `s_j -= α_j²` with recompute guard `s_j ≤ s_ref_j·tol` instead of from-scratch norms | all three kernels | recomputation returns the same mathematical quantity; guard is the writeup's Remark "Column norm downdating" | `test_oracle_parity` (oracle uses exact norms; pivot sequences must still match); Julia testset "Norm recompute tolerance knob" |
| 2 | Householder sign convention `ρ = -sign(x₁)‖x‖` (LAPACK), identity reflector when the tail is zero; writeup states `ρ = ‖x‖ ≥ 0` | all three kernels and the oracle | row sign flips cancel in `(D·R11)⁻¹(D·R12)`, so every `w_j` and the criterion are unchanged; `Q` absorbs `D` (proof sketch in `oracle_bsqr.m` header) | direct (sign-normalization-free) Q/R comparison in `test_oracle_parity` |
| 3 | Clamp `wnorm2` and `s` at 0 after downdates | all three kernels | guards `-eps`-level negatives of nonnegative quantities | indirectly by parity; dedicated near-zero unit test pending |
| 4 | `s(i) = β_i²` refresh of the pivot column post-reflection | all three kernels | exact identity after the reflection | `test_oracle_parity` |
| 5 | BLAS-2 organization of the W update (`gemv` for dots, `ger` for the rank-1) | Julia, MEX | writeup Remark "Practical considerations"; same recurrence, regrouped | `test_oracle_parity` (oracle has no recurrence at all) |
| 6 | `rank_stop` early-exit option | Julia only, off by default in `bsqr` | extension, not in Algorithm 1; cross-language defaults agree (no early exit) | Julia testset "Rank-stop policy"; cross-language default parity pending (V3 fixtures) |
| 7 | Short-wide fastpath (inline W-row materialization) | Julia only | pure memory-layout reordering; identical arithmetic | Julia testset "Kernel helper invariants and fastpath knobs" |
| 8 | First-minimum tie-breaking via strict `<` | all three kernels and the oracle | writeup's argmin leaves ties unspecified; fixed convention | `testPivotTieStability` (MATLAB); Julia "Criterion-consistent pivot sequence" |

## Validation assets

| Asset | Status |
|---|---|
| V1 oracle (`matlab/tests/oracle_bsqr.m`) | done |
| Oracle/backend parity + invariants + Osinsky canaries (`matlab/tests/test_oracle_parity.m`) | done (mfile + MEX) |
| V0 paper → writeup audit | pending |
| V3 cross-language fixtures (Julia ↔ MATLAB, incl. Julia ↔ oracle) | pending |
| V4 per-step trace parity instrumentation | pending |
| V5 Julia `larfg!` unification (or extreme-scale equivalence tests) | pending |

Notes on test design choices:

- Exact pivot-sequence equality across implementations is demanded only on inputs screened for
  criterion near-ties (`testZooIsTieFree` asserts every step's runner-up gap exceeds `1e-6`
  relative). With fixed seeds this is deterministic; a failing member gets a new seed, never a
  looser tolerance.
- Forward agreement of `Q`, `rinv_r12`, and `inv(R11)` is condition-sensitive even with
  identical pivots, so graded-spectrum zoo members check pivot sequence + `R` agreement only;
  their `Q`/solve quality is covered by backward-style contract tests instead.
