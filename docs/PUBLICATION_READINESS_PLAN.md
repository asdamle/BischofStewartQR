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
      and the writeup (then `notes/bischof_stewart_pivoting.tex`, now
      Appendix A of `notes/GKSevolved_draft.tex`); confirm every
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
      *Partially covered 2026-07-05 by merging the stale `audit/*` branches:*
      `matlab/tests/test_mfile_mex_fuzz.m` (seeded mfile-vs-MEX differential
      fuzz, now a standing test), a randomized unblocked-vs-panel Julia
      parity testset, and six new parity-zoo stress corners.
      **Cross-language half done 2026-07-05**: opt-in driver
      `matlab/tests/fuzz_cross_language_gen.m` (randomized zoo-format
      fixtures with the near-tie screen; five families incl. graded and
      near-tied norms) + `julia/test/fuzz_cross_language_check.jl` (both
      kernel paths per case). Batch of 2000 seeds (1866 kept, 134 near-tie
      skips) × 2 kernel paths: **0 failures**. Along the way found and fixed
      a latent writer bug in `generate_parity_fixtures.m` (MATLAB `fprintf`
      with an empty data array still emits the format's literal characters,
      so 1-column matrices gained a stray trailing comma — unreachable with
      the current zoo, fatal for any future k = 1 member). One calibration:
      criterion-trace parity is condition-amplified, so for randomly-graded
      fuzz cases the checked contract is exact pivots + tight R (the
      13 initial "failures" were final-step criterion drift ~3e-3 on
      numerical-rank tails with pivots exact and R at ~5e-16).
- [x] **Sanitizer builds of both MEX kernels.** ASan is unusable inside
      MATLAB on macOS — dyld's platform policy refuses to load the sanitizer
      runtime into a hardened-runtime process ("Sanitizer load violates
      platform policy"), and static ASan for dylibs does not exist. Used the
      strongest embeddable instrumentation instead: trap-mode UBSan
      (`-fsanitize=undefined -fsanitize-trap=undefined`, no runtime needed)
      plus libc++ hardened mode (`_LIBCPP_HARDENING_MODE_DEBUG`, real vector
      bounds checks). Build recipe:
      `build_bsqr_mex('CXXFLAGS=$CXXFLAGS -fsanitize=undefined
      -fsanitize-trap=undefined -fno-omit-frame-pointer
      -D_LIBCPP_HARDENING_MODE=_LIBCPP_HARDENING_MODE_DEBUG')` (same for
      `build_bsqr_rand_mex`). **Found and fixed one real UB**: the
      deterministic MEX's column-norm init loop formed a reference into an
      empty vector for m = 0 inputs (benign in production only because
      `dnrm2` with len 0 never dereferences). After the fix, both full
      suites plus a 720-run randomized stress sweep run clean under
      instrumentation.
- [x] **Coverage-guided gap analysis** (2026-07-05). Julia
      `Pkg.test(coverage=true)`: interface/workspace 100%; kernel 280/294,
      panel 144/151 executable lines. MATLAB `CodeCoveragePlugin` (Cobertura)
      over both suites: kernels 93–100% per file. Gaps classified:
      *(a) instrumentation hooks* (`kernel_stats` timing, `recomp_history`,
      `pivot_history` in the panel) — benchmark-only plumbing, justified;
      *(b) unreachable-by-design* — `bsqr_rand_threshold`'s `otherwise` (mode
      validated upstream), the single-select rounding-tie fallback
      (`accept_id == 0` requires the scanned global minimizer to fail the
      threshold, impossible in exact arithmetic), one `return_r12` arm only
      reachable by bypassing the dispatcher; *(c) real test gaps — all now
      closed with committed tests*: long-tail (>256) BLAS recompute branch
      (Julia), rank-stop inside the panel kernel (Julia), exact rank
      deficiency past the numerical rank mid-kernel (`beta = 0` →
      `invdiag = 0`/`tau = 0` guards, both MATLAB backends), zero-tail
      Householder in the rand kernel, single-select uniform visiting order
      and `pick = 'first'` (the knob test ran them batched, where `pick` is
      ignored), the single-select rank guard, `R12` at `k = n`, m-file
      `k = 0` rinv branch, and every parser rejection identifier in both
      dispatchers.
- [x] **Targeted adversarial cases** — probed on both backends (some under
      the instrumented MEXes) and codified as `testAdversarialCorners` in the
      randomized suite: `block_size = 1` (both kernel paths) and `10·n`,
      uniform sampling on `needle` with block 1 (forces the fallback path),
      `worstcase_allowance` at `k = 2` (final-step `k−i = 1` division),
      `slack = 1e3` (bound void, exactness preserved), and fully-underflowed
      sampling weights (squared norms leave the documented domain → the rank
      guard fires with the documented error, on both backends). All behaved
      correctly. Panel-boundary interactions are covered by the merged
      randomized unblocked-vs-panel testset plus the fuzz checker's forced-
      panel pass.
- [x] **Optional deep review.** Decided 2026-07-05: skipped by the author.
      The kernels have since been covered by the oracle deviation audit, the
      sanitizer pass, 1866-case cross-language fuzz, and coverage-gap
      closure, which together substitute for much of what an adversarial
      review panel would probe.

Any bug found: fix with a regression test, then re-run Phase 0's baseline
checks (parity is the canary for kernel-behavior changes — both languages
must change together).

## Phase 3 — Documentation validation (1 session) — DONE 2026-07-05

Inventory: `CLAUDE.md`, root `startup.m`/`startup.jl`, `matlab/README.md`,
`matlab/mex/README.md`, `matlab_rand/README.md`, `julia/README.md`,
`julia/docs/TESTS_AND_BENCHMARKS.md`, the seven `docs/*.md`, and inline help
(`bsqr.m`, `bsqr_rand.m`, Julia docstrings — already refreshed this cycle).

- [x] **Execute every documented command.** (2026-07-05: every snippet in
      CLAUDE.md and all READMEs/startups ran green — suites, builds,
      `startup`/`startup rebuild`/`startup.jl`, quick-starts (incl. the
      5-output `trace` and `return_r12` forms), smoke benchmarks, plotters
      (pointed at smoke CSVs so committed artifacts stayed untouched), perf
      gates, the Python timing compare, fixture regen (byte-identical), and
      all matlab_rand comparison suites at reduced-but-complete scale — every
      documented CSV/figure name produced; `rejection_rpqr` run against a
      local `ext_comparisons/` download (URL in the README verified live).
      Full publication runs validated via their smoke variants. Fixes: output
      names in the root-README/startup quick-starts (`[p, Q, R11]`,
      `[Q, R, E]`) and `startup.jl`'s `reconstruct(F, A)` hint.)
- [x] **Check every factual claim against code** (2026-07-05: defaults
      (`BS_PANEL_NB=8`, `BS_PANEL_MIN_KN=24576`, `norm_recomp_tol=sqrt(eps)`,
      block-size clamp `[16,64]`), both env-var lists complete and correct,
      output shapes/identities, stats/trace fields, figure/CSV names, and
      complexity statements all verified. Drifts fixed: CLAUDE.md's Julia
      test-file list, and the MEX promotion rule restated to match reality —
      `auto` already prefers the built MEX; the ≥20% short-wide criterion is
      kept and will be evaluated at the Phase 6 rerun.)
- [x] **Flag stale measured numbers.** (2026-07-05: the
      `docs/RANDOMIZED_BSQR_PLAN.md` speedup table carries the pre-lazy-Q
      annotation; `matlab_rand/README.md` has no measured tables. Replacement
      with rerun numbers stays a Phase 6 item.)
- [x] **Cross-document consistency** (2026-07-05: single-source sweep across
      CLAUDE.md / READMEs / startups clean after the fixes above;
      `VALIDATION_AND_PERF_PLAN.md` and `PERF_P0_FINDINGS.md` (all items
      done) marked as historical records at the top.)
- [x] **Manuscript-facing docs** (2026-07-05: `RANDOMIZED_BSQR_ALGORITHM.md`
      formulas re-verified against both randomized kernels (no kernel-visible
      change since the audit); `docs/VALIDATION.md`'s deviation list, test
      pins, and P3 section references all verified current. Same-day
      follow-up: `notes/` now holds the manuscript draft
      (`GKSevolved_draft.tex` / `manuscript_draft.pdf`), superseding the
      standalone note — every repo reference re-keyed to the manuscript's
      rendered numbering (Alg. A.1 and its lines, Thm. 3.1, Cor. 3.2,
      Thm. 5.1, eqs. (2.3)/(2.4)/(A.2), Remark 7), verified against the
      compiled PDF.)

## Phase 4 — Comment quality pass (1 session)

File-by-file over the ~4000 lines of kernel/private code plus benchmark
drivers (`julia/src/*.jl`, `matlab/private/*.m`, `matlab/mex/src/bsqr_mex.cpp`,
`matlab_rand/private/*.m`, `matlab_rand/mex/src/bsqr_rand_mex.cpp`).

Done 2026-07-05 via a four-way audit (parallel read-only reviews of the MEX
kernels, Julia sources, MATLAB m-files, and benchmark drivers, each finding
verified against code before editing):
- [x] Invariants-not-narration: dev-history labels removed ("P3 prototype"
      on the now-default panel kernel, "(P3)" dispatch tag); the three
      m-file kernel helpers (`householder_column`, `apply_householder_left`,
      `build_q_from_factors`) gained contract comments (packed-form
      conventions, tau = 0 identity encoding, backward accumulation).
- [x] Loop invariants: nontrivial blocks verified stated (panel deferred
      updates, in-block contiguity, Fenwick ops, wnorm2/`s`/`s_ref` meaning).
- [x] Terminology: "pivot score" → "pivot criterion" and "LAPACK-style
      guard" → "Businger-Golub recompute safeguard" in `bsqr_mex.cpp`; the
      rand MEX header now names the criterion increment and matches the
      m-file's "for every column at every step" qualifier.
- [x] No stale facts: every constant named in a comment verified against
      code (panel nb = 8, crossover 24576, block clamp [16,64], sqrt(eps),
      256 recompute threshold, P1.1 section reference); `oracle_bsqr.m`'s
      eight "Alg. 1" labels re-keyed to the manuscript's "Alg. A.1".
      Rejected findings recorded for the record: the timed-product comment
      in `run_rand_benchmarks.m` is correct (the 4-output call nearby is
      the untimed stats pass), and "norm downdates" (plural) matches
      repo-wide usage.
- [x] `grep -rn "TODO\|FIXME\|XXX"` is zero.

## Phase 5 — Performance-opportunity search (1–2 sessions; search, then decide)

Goal is to *find and evaluate*, not to reflexively apply — any kernel change
re-triggers Phase 1's parity work. For each finding: measure, write down the
expected gain, and either apply (with perf-gate + parity evidence) or record
as consciously deferred.

Done 2026-07-06. Findings and decisions:

- [x] **timeit everywhere (user-directed).** All MATLAB benchmark timings now
      go through `timeit`, so recorded variability is across random problem
      instances (the seed/trial columns), not timer noise. Applied:
      `run_publication_benchmarks.m` (`bench_function` was warmup +
      tic/toc-samples; now one `timeit` per instance — the dead `tmin_s`
      column dropped from the CSV, the now-meaningless `warmup`/`samples`
      knobs and `BS_MATLAB_PUB_{WARMUP,SAMPLES}` env vars removed, schema
      bumped to 2026-07-06.matlab.v5) and `run_approx_comparison.m` (single-
      shot tic/toc per trial → `timeit` with the selector re-seeded inside
      the thunk so every repetition times the identical instance). The four
      other drivers already used `timeit`. Verified: smoke run green, perf
      gate vs the Phase 0 baseline passes (0 violations — methodology change
      is within noise), full suite green.
- [x] **Profile/allocation check.** The preallocated `bsqr!` path is
      allocation-free after warm-up on the unblocked kernel; the panel
      kernel allocates ~`3·n·nb` doubles of panel-local scratch per call
      (~124KB at 64×512). Timing impact is negligible (µs vs ms per
      factorization) → fix deferred; the `bsqr!` docstring now states the
      caveat and the `BS_PANEL_NB=0` escape hatch. Caching the panel
      buffers inside `BSWorkspace` is the future fix if strict allocation-
      freedom ever matters.
- [x] **Knob sanity.** `BS_PANEL_NB` swept {0,4,8,16} on Accelerate over
      128×256 / 192×192 / 256×1024 / 64×640: nb=8 fastest on every shape
      (e.g. 256×1024: 16.7ms unblocked → 9.5ms at nb=8); default confirmed,
      crossover default consistent with the file's recorded boundary.
- [x] **Candidates evaluated.** `check_finite`/`check` confirmed disabled on
      every timed path in both languages; the accessor copies (`R(F)`,
      `perm`, `rinv_r12`) confirmed absent from timed loops; deterministic
      MEX one-output Q-skip confirmed in place. Randomized-MEX two-column
      blocking and Fenwick-vs-partial-sums at small n: consciously deferred
      — research-grade kernel changes with parity costs, out of scope
      pre-submission.
- [x] **Regime map vs `dgeqp3`** (preliminary, from the committed
      publication CSVs; regenerate at Phase 6): the deterministic BSQR with
      the full `[Q,R,p]` product trails `dgeqp3` everywhere — geomean ratio
      1.2–1.5× on square and moderate aspects, worst ~2–2.2× at the widest
      short-wide shapes (512×4096/5120), consistent in both languages. For
      the paper: never claim raw-speed wins for the deterministic kernel;
      its value is criterion quality + the free `R11⁻¹R12` (and the `rinv`
      product's comparison vs `qr`+`trsm` is the fair speed story); raw
      speedups belong to the randomized variant. This also re-confirms the
      MEX promotion rule's status quo.
- [x] Apply/defer decisions recorded above; perf gate run against the
      Phase 0 baseline for the applied timing change (0 violations).

## Phase R — External review response (2026-07-10; gates Phase 6)

An external review of the repository found the core algorithms clean (no
undocumented mathematical deviations; all validation metrics measure what
their labels say; zero bound violations across 9,600 conditioning rows) and
two substantive experiment-aggregation problems plus six documentation
drifts. All findings were independently verified before this plan was
written. Author guidance: documentation inconsistencies resolve toward the
current code (the figures show what they should show).

**R1 — Deduplicate the MATLAB publication baseline (code fix). DONE 2026-07-10.** Verified:
`bench_pair` is called twice per case (plain + rinv) and both calls time
`qr_pivoted`, so the committed CSV has 300 baseline rows vs 150 per BSQR
variant, and the `pair_relative_times` join (keys `family…seed`) pairs each
BSQR row with both duplicates — double-weighted plain relative-time artifacts
(up to 22.3% group sensitivity per the review). Julia is clean (one `dgeqp3`
row per case via include-flags).
- [x] Restructure `run_publication_benchmarks.m` to mirror Julia: the plain
      pass times `bsqr_full` + `qr_pivoted`; the rinv pass times `bsqr_rinv`
      + `qr_pivoted_trsm` only. Downstream joins and seed-grouped plots then
      pair 1:1 with no aggregation change needed.
- [x] Add a cheap invariant to the runner: exactly one row per
      (case, method) before the CSV is written — prevents recurrence.
- [x] Verify by smoke run (multiplicities 1:1:1:1) and a perf-gate pass
      (note: the Phase 0 baseline CSV carries the duplicate rows; the gate
      joins per method on `tmed_s`, so the comparison still works — expect
      duplicated join rows against the old baseline, harmless).
- [ ] The committed CSV/plots/tables stay as-is until the Phase 6 rerun
      regenerates them; until then, treat MATLAB *plain* relative-time
      artifacts as double-weighted (quality metrics unaffected).

**R2 — Count all randomized sampling work (code fix, both backends in
lockstep). DONE 2026-07-10.** Verified: on the batched path a sampled block that yields no
selection is evaluated but never recorded (`samples_tested`/`rounds` are
written only by selection-yielding blocks), so `total_tested`,
`blocks_sampled`, and the plotted `tested_per_k` undercount — demonstrated
at ~8× on uniform/needle (batched reports 104 vs single-select's 824 for the
same problem). The `fallback` flag records only forced above-threshold
acceptance, narrower than the help text implies.
- [x] m-file batched path: carry pending counters for failed blocks'
      candidate counts and reflector applies; attribute
      `pending + current block` to the next selection (the exhaustive-scan
      path also adds its pending). `total_tested` and `blocks_sampled`
      become true totals of candidate evaluations and block applies.
- [x] MEX batched path: identical change (`since_last` already tracks the
      candidate count; add the rounds counter), keeping the two backends'
      instrumentation semantics identical.
- [x] Documentation: `bsqr_rand.m` stats help and the README stats section
      updated to the now-true semantics; `fallback` documented precisely as
      "forced acceptance above the threshold on an exhaustive pass"
      (exhaustive scans are visible through `rounds`).
- [x] Tests: regression case where no block ever fails (block ≥ n —
      counts must equal today's) and an undercount case (uniform/needle,
      small block) asserting batched `total_tested` is now commensurate with
      the single-select count rather than ~8× below it; full rand suite +
      a bounds-stress spot check.
- [x] Affected figures (`tested/m` panels of blocksize + sampling) and any
      `tested_per_k` numbers regenerate at Phase 6; selection, timing,
      conditioning, and accuracy outputs are untouched by construction.

**R3 — Documentation drift, resolved toward the code (all verified). DONE 2026-07-10:**
- [x] `plot_approx_cond_comparison.m` header: row 4 is the
      orthogonal-projection error (the axis label is already correct);
      remove the "noisy ID" row-4 description.
- [x] `matlab_rand/README.md` + `docs/RANDOMIZED_BSQR_PLAN.md`: the
      projection-coefficient / noisy-ID curves are recorded in the CSV but
      not plotted — fix the "thin dashed `T_proj` overlaid" claims.
- [x] `plot_rand_experiments.m` + `plot_rpqr_comparison.m` comments (4
      sites): the spectral row plots `‖R₁₁⁻¹‖₂ / bound ≤ 1` (lower =
      better), not `σ_min/bound ≥ 1`.
- [x] `bench_common.jl` + `julia/docs/TESTS_AND_BENCHMARKS.md`: describe
      `tci_low/high` as empirical 2.5%/97.5% timing quantiles, not
      confidence intervals (CSV column names unchanged).
- [x] `docs/PUBLICATION_FIGURES_PLAN.md`: replace the removed
      `table_quality.csv`/p95 spec with the implemented
      `quality_summary.md` (median/max).
- [x] `slack > 1` wording in `docs/RANDOMIZED_BSQR_ALGORITHM.md` and
      `docs/RANDOMIZED_BSQR_PLAN.md`: it does **not** weaken the bound
      "proportionally" (the recursion does not telescope that way); align
      with the API help — the Osinsky guarantee no longer applies.

Order: R3 first (pure documentation, no verification burden), then R1
(runner restructure + smoke + invariant), then R2 (kernel instrumentation +
tests, both backends together). Full suites in both languages green at the
end of each code step; Phase 6's rerun then regenerates every affected
artifact on the corrected pipeline.

## Phase 6 — Final regeneration and freeze (already planned)

Driver: `./phase6_rerun.sh` at the repo root runs the whole sweep
sequentially (MEX builds → both publication runs → smoke perf gates vs the
Phase 0 baselines → all matlab_rand suites → final test suites), with
per-step logs and resumable markers in `.phase6/`, self-wrapped in
`caffeinate` so screen lock does not stop it. Launch it detached so it also
survives the terminal/session: `nohup ./phase6_rerun.sh >/dev/null 2>&1 &`,
then `./phase6_rerun.sh --status` / `tail -f .phase6/console.log`. Keep the
machine on AC with the lid open, and idle (these are timing runs).

- [x] Rerun both languages' publication benchmarks on the final tree
      (2026-07-13, via `./phase6_rerun.sh` on an idle machine: all 11 steps
      green; MATLAB CSV at exact 1:1:1:1 method multiplicity under the R1
      invariant); artifacts committed. Also rerun
      `run_largen_scaling` (author-approved; deferred 2026-07-10; the default
      sweep now extends one doubling past the regime crossing, to
      n = 2,048,000, so the figure ends on re-linearized curves). Note on the
      needle tail: the "odd" last point (n = 1024k) is *reproducible*, not
      noise — a fresh tail rerun showed every method stepping superlinearly at
      the 512k→1024k doubling (`bsqr_nw` 3.2×, `bsqr_unif` 3.0×, `rpqr` 2.3×),
      consistent with the ~½GB working set crossing a memory-bandwidth regime
      on this machine; it reads oddest on `bsqr_nw` only because that curve is
      lowest. Worth one caption sentence rather than a fix.
- [x] Perf gates pass against the recorded Phase 0 baselines (both languages, 0 violations).
- [x] Replace the stale-flagged measured tables (Phase 3) with rerun numbers (plan doc speedup table: 53×/110×/148× vs det, 19×/25×/52× vs dgeqp3 at k=64 gaussian; stress families 117×/134× at the top size).
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
