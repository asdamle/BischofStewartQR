# Publication-Readiness Plan

Working plan for the final pre-submission pass over the repository. Five
workstreams: correctness verification, bug hunting beyond the current tests,
documentation validation, comment quality, and a performance-opportunity
search. Ordered so that later phases describe/polish behavior the earlier
phases have verified; each phase has explicit acceptance criteria so progress
is checkable across sessions.

**Definition of done.** Both language test suites green with MEX built;
parity fixtures regenerated and exactly matching across Julia / MATLAB m-file
/ MEX; every command and factual claim in the docs executed or checked
against code; kernel comments state invariants rather than narration;
performance findings either applied (with perf-gate evidence) or written down
as consciously deferred; publication artifacts regenerated from the final
tree (the planned benchmark rerun) and gates passing.

---

## Phase 0 — Baseline (half a session) — DONE 2026-07-05

Establish a known-good starting point before changing anything.

- [x] `julia --project=julia -e 'using Pkg; Pkg.test()'` green.
- [x] `matlab -batch "addpath('matlab/tests'); run_tests"` green (MEX built).
- [x] `matlab -batch "addpath('matlab_rand'); addpath('matlab_rand/tests'); run_rand_tests"` green (MEX built).
- [x] Parity drift guard passes (fixtures in `parity/` current w.r.t. `parity_zoo.m`).
- [x] Record baseline timings: `run_publication_smoke_benchmark` (both languages)
      stashed for later perf-gate comparison. Stashed (local, gitignored) at
      `julia/benchmark/results/tmp_phase0_baseline/` and
      `matlab/benchmark/results/tmp_phase0_baseline/` (timings CSV +
      metadata + environment note each).
- [x] `git status` clean; toolchain recorded in the stash's
      `phase0_env.txt`: commit `c480c58`, macOS 15.7.7 arm64, Julia 1.11.6
      (Accelerate for benchmarks), MATLAB R2025b.

## Phase 1 — Correctness verification (1–2 sessions)

The existing harness is strong (oracle transcription, cross-language parity
fixtures, cancellation stress, bound checks). This phase re-validates it and
extends coverage where the harness is thin.

**Re-validate the existing harness:**
- [x] Refresh the `docs/VALIDATION.md` deviation audit: re-diff both kernels
      (unblocked + panel/MEX in-block) against `matlab/tests/oracle_bsqr.m`
      and the writeup (`notes/bischof_stewart_pivoting.tex`); confirm every
      documented deviation is still the complete list. (2026-07-05: no
      undocumented kernel/oracle mismatch; added deviations #10–#12 — panel
      deferred/regrouped updates, batched norm refresh, zero-tail-norm
      extension — and fixed stale scope/pins in #5–#8.)
- [x] Regenerate parity fixtures and confirm byte-exact pivot sequences
      across all three implementations on the full zoo. (2026-07-05:
      regenerated fixtures are byte-identical to the committed ones.)
- [x] Tie coverage. Revised: the zoo *deliberately screens out* near-ties
      (tie outcomes are not BLAS-portable, so cross-implementation tie parity
      would be an unfair contract). Added per-implementation tie tests
      instead — bitwise-duplicate columns at step 1 and mid-factorization,
      first minimum must win — in the Julia suite and both MATLAB backends.

**Extend coverage (new committed tests):**
- [x] Shape/edge grid, both languages: empty (0×n, m×0, 0×0), single
      row/column, tall, k = 0 / 1 / min−1. All passed as-is (Julia "Edge
      shapes" testset; MATLAB `testEdgeShapes`).
- [x] Early-stop + `rank_stop` accessor interactions (Julia "Accessors on a
      rank-stopped factorization" testset): all five accessors verified on
      `ksteps < k`, including `rinv_r12` and exact `reconstruct`.
- [x] Julia strided inputs — **found a real bug**: a row-strided view
      (`stride(A,1) > 1`) type-checks as `StridedMatrix` but silently
      produced a garbage factorization (residual ~1e2, orthogonality ~1e6);
      the kernel's BLAS calls assume LAPACK layout. Fixed by rejecting
      non-unit first stride in `_validate_bsqr_common_args` with a clear
      error; offset views (unit stride) verified to match dense runs
      ("Strided input contract" testset).
- [x] MATLAB input-type edges — **found a real bug + an API inconsistency**:
      (1) sparse inputs pass both MEXes' double/real/2-D validation, but
      `mxGetPr` on sparse yields nonzero storage only — reading it as dense
      is a buffer overread (masked in casual testing because fully-dense
      sparse matrices have `nzmax = m*n`). Both MEXes now reject sparse
      explicitly, as do both dispatchers. (2) single/integer inputs worked on
      the m-file backends but errored on the MEXes, making behavior depend on
      which backend was built; dispatchers now normalize to double before
      dispatch. Tests: `testInputTypeContract` in both suites.
- [x] `bsqr_mex` persistent-workspace sequence test
      (`testMexWorkspaceReuseSequence`): 13 interleaved grow/shrink/reshape/
      early-stop calls in one process, each checked for exact pivot parity
      against a fresh m-file reference plus the tight residual identity.
- [x] Extreme scaling — probed all three implementations with finite inputs
      whose squared norms overflow (`‖a‖ ~ 1e160`) / underflow (`~1e-170`):
      the factorization stays exact everywhere (residual ~2e-16; Householder
      generation is scale-safe); only the pivot criterion, and hence the
      selection guarantee, degrades. Decision: document the guarantee's
      domain (column norms roughly in `[1e-150, 1e150]`) in all three help
      texts, and codify factorization-stays-exact as tests
      (`testExtremeColumnScaling`, Julia "Extreme column scaling" testset).
- [x] Tolerance audit (2026-07-05, user-requested): swept every
      tolerance-bearing assertion in all three suites for vacuous bounds.
      One further offender fixed (`testOptionalRinvR12` used the vacuous
      `rel_resid(A(:,p), Q*R) < 1.0` on an early-stopped factorization —
      replaced with the tight identities). Everything else is eps-scaled or
      a deliberate, documented conditioning allowance.
- [x] Randomized variant, statistical pass:
      `matlab_rand/tests/stress_bsqr_rand_bounds.m` (opt-in, not
      auto-discovered). 2026-07-05 run: 100 seeds × {batched, single-select}
      × {uniform, normweighted} × {running_mean, worstcase_allowance} ×
      {gaussian, spiked_leverage, needle} = 2400 runs — per-step
      `f2 ≤ Fhat` held in every run; quality ratio median 0.263, max 0.994;
      median 16 samples tested per selection.
- [x] Doc-claims tests: the documented identities are now asserted in
      committed tests — Julia "Q accessor"/"Edge shapes"/rank-stop testsets;
      MATLAB `testOutputContracts` + `testEdgeShapes` + a strengthened
      `testEarlyStoppingByK` (its old identity check used tolerance 1.0
      because `A(:,p) ≈ Q*R` is the wrong identity under early stop; it now
      asserts the tight selected-block and R12-projection identities).

## Phase 2 — Bug hunt beyond the tests (1–2 sessions)

Techniques, in order of expected yield:

- [ ] **Differential fuzzing.** A driver that generates random matrices
      (shapes, scalings, structured families, near-ties) and compares pivot
      sequence + `R` + `tau` exactly between MATLAB m-file and MEX, and
      against Julia via CSV interchange; any divergence is a bug by the
      repo's own lockstep contract. Run a large batch (10⁴–10⁵ cases)
      overnight; keep the driver in the repo as an opt-in stress tool.
- [ ] **Sanitizer builds of both MEX kernels.** Build with
      `-fsanitize=address,undefined` (mex CXXFLAGS variant) and run the full
      MATLAB suites + the fuzz driver under them. This is the only realistic
      memory-safety check for the ~2200 lines of C++ (macOS: ASan works,
      valgrind does not).
- [ ] **Coverage-guided gap analysis.** Julia `Pkg.test(coverage=true)` and
      MATLAB profiler coverage over the suites; list uncovered branches in
      kernel files (candidates: recompute safeguards, rank-stop inside the
      panel kernel, Fenwick degenerate paths `total <= 0`, the batched
      feasibility net, swap-pop pool maintenance). Write a test per uncovered
      branch or justify why it is unreachable.
- [ ] **Targeted adversarial cases** (from reading the kernels this session):
      panel-boundary interactions (`k % BS_PANEL_NB ≠ 0`, `BS_PANEL_MIN_KN`
      crossover exactly at `k·n`), `block_size > n` and `block_size = 1` in
      the randomized MEX, all-zero-weight columns under norm-weighted
      sampling (Fenwick `total ≤ 0` fallback), the batched fallback firing on
      the *first* block (`since_last ≥ rem` at `nsel = 0` requires `b ≥ n`),
      `worstcase_allowance` at the final step (`k−i = 1` division), and
      `slack` large enough that the bound is intentionally void.
- [ ] **Optional deep review.** A multi-agent adversarial review of the two
      MEX kernels and the two pivot-criterion state machines is the highest-
      value target for `/code-review ultra` (user-triggered) or an explicit
      workflow run; the kernels are subtle enough that independent skeptical
      readers are worth the tokens. Decide whether to spend it.

Any bug found: fix with a regression test, then re-run Phase 0's baseline
checks (parity is the canary for kernel-behavior changes — both languages
must change together).

## Phase 3 — Documentation validation (1 session)

Inventory: `CLAUDE.md`, root `startup.m`/`startup.jl`, `matlab/README.md`,
`matlab/mex/README.md`, `matlab_rand/README.md`, `julia/README.md`,
`julia/docs/TESTS_AND_BENCHMARKS.md`, the seven `docs/*.md`, and inline help
(`bsqr.m`, `bsqr_rand.m`, Julia docstrings — already refreshed this cycle).

- [ ] **Execute every documented command.** Each shell/matlab/julia snippet
      in every README and CLAUDE.md gets run (smoke-scale where the real run
      is expensive). Fix or delete anything that errors or drifts.
- [ ] **Check every factual claim against code**: defaults (block sizes,
      `norm_recomp_tol`, panel width, env-var names `BS_PUB_*` /
      `BS_MATLAB_PUB_*`), output shapes, complexity claims, and the MEX
      promotion rule in `matlab/mex/README.md` (evaluate it against the final
      benchmark numbers — either promote or keep the rule and note status).
- [ ] **Flag stale measured numbers.** Tables in `docs/RANDOMIZED_BSQR_PLAN.md`
      and READMEs carrying pre-rerun timings are already annotated; after the
      final rerun, replace them with the regenerated numbers.
- [ ] **Cross-document consistency**: CLAUDE.md vs READMEs vs `docs/` plans
      (one description per fact, others link); mark completed plan docs
      (`VALIDATION_AND_PERF_PLAN.md`, `PERF_P0_FINDINGS.md`, this file when
      done) as historical records in a line at the top.
- [ ] **Manuscript-facing docs**: `docs/RANDOMIZED_BSQR_ALGORITHM.md` §6 was
      audited and compile-checked this cycle; re-verify against the final
      kernel behavior after any Phase 1/2 fixes, and confirm theorem/equation
      references into `notes/` still match the writeup.

## Phase 4 — Comment quality pass (1 session)

File-by-file over the ~4000 lines of kernel/private code plus benchmark
drivers (`julia/src/*.jl`, `matlab/private/*.m`, `matlab/mex/src/bsqr_mex.cpp`,
`matlab_rand/private/*.m`, `matlab_rand/mex/src/bsqr_rand_mex.cpp`).

Criteria (delete or rewrite anything failing them):
- [ ] A comment states an invariant, constraint, or non-obvious *why* — not
      what the next line does, not development history ("landed
      optimization", section numbers of superseded plans) unless it earns its
      place as a pointer a maintainer needs.
- [ ] Each nontrivial kernel block states its loop invariant once (`W` rows
      maintained, `s`/`s_ref` meaning, compact-WY `T` recurrence, in-block
      contiguity of active columns) — most already do; make it uniform.
- [ ] Terminology is consistent across languages and with the writeup
      (criterion, downdate, recompute safeguard, head/tail, compact-WY) —
      the two implementations are meant to be read side by side.
- [ ] No stale facts (the `[16,128]` clamp comment found this cycle is the
      cautionary example): every constant named in a comment greps back to
      code that matches.
- [ ] `grep -rn "TODO\|FIXME\|XXX"` stays at zero (currently zero).

## Phase 5 — Performance-opportunity search (1–2 sessions; search, then decide)

Goal is to *find and evaluate*, not to reflexively apply — any kernel change
re-triggers Phase 1's parity work. For each finding: measure, write down the
expected gain, and either apply (with perf-gate + parity evidence) or record
as consciously deferred.

- [ ] **Profile before guessing.** Julia: `@profile` + allocation tracking on
      publication-shaped cases (square and short-wide at the benchmark
      sizes); confirm the kernel loop is allocation-free. MEX: sample with
      Instruments on the same cases.
- [ ] **Known knobs to sanity-check rather than re-derive**: `BS_PANEL_NB`
      (default 8) and `BS_PANEL_MIN_KN` crossover on the current toolchain;
      BLAS thread pinning behavior at benchmark sizes; Accelerate vs default
      BLAS deltas already handled by `setup_accelerate.jl`.
- [ ] **Candidates spotted while reading the kernels this cycle** (each needs
      a measurement before touching):
      - deterministic MEX: whether the one-output form skipping Q is
        exploited everywhere it could be internally;
      - randomized MEX: `dger`-based in-block downdates are BLAS-2 by design —
        check whether a two-column blocking of the reflector apply pays at
        `b = k` sizes; Fenwick vs plain partial sums at small `n`;
      - `check_finite` is an `O(mn)` pre-scan both sides pay by default —
        confirm benchmark paths disable it (they do) and that the default is
        documented as safety-first;
      - Julia `R(F)` returns `triu(view(...))` (allocates) and `perm`/
        `rinv_r12` copy — fine for the API, but confirm no benchmark-timed
        path calls them inside loops.
- [ ] **Regime map vs `dgeqp3`**: one table of BSQR/baseline ratio over the
      full (m, n, k) publication grid from the rerun — identifies any regime
      where a claim in the paper would be weak, which is the real
      "performance improvement" that matters for publication.
- [ ] Apply/defer decisions recorded here; perf gate
      (`check_publication_perf_gate.*`) run against the Phase 0 baseline for
      anything applied.

## Phase 6 — Final regeneration and freeze (already planned)

- [ ] Rerun both languages' publication benchmarks on the final tree
      (the user's planned pre-submission rerun); commit artifacts.
- [ ] Perf gates pass against the recorded baseline.
- [ ] Replace the stale-flagged measured tables (Phase 3) with rerun numbers.
- [ ] Final end-to-end: fresh clone → `startup` (MATLAB) / `startup.jl`
      (Julia) → tests green → smoke benchmarks run, following only the
      READMEs. Anything that requires knowledge not in the docs is a Phase 3
      failure; fix and repeat.

---

## Sequencing and effort

Correctness before docs (docs must describe verified behavior); comments
after bug-hunt fixes land (no polishing code that is about to change); perf
search last among code phases (its changes are the riskiest and gated by the
re-validated harness); regeneration at the very end. Rough total: 6–9
working sessions. Phases 1–2 are the highest-value and least parallelizable;
Phases 3–4 are mechanical and can interleave with waiting on long fuzz /
benchmark runs.

Standing rules while executing: kernel-behavior changes require the matching
change in the other language plus fixture regeneration (CLAUDE.md contract);
`matlab_rand` stays decoupled from the deterministic kernels; every doc claim
that can be executed becomes a test.
