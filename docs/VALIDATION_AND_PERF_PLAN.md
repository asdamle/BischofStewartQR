# BSQR Validation and Performance Plan

Two workstreams:

1. **Validation** — establish, with executable evidence, that the Julia kernel, the MATLAB
   m-file kernel, and the MEX kernel all implement the Bischof–Stewart algorithm exactly as
   specified by the references (Bischof 1990, Stewart 1990, and Algorithm 1 of
   `notes/bischof_stewart_pivoting.tex`), modulo a short, documented list of mathematically
   equivalent stability modifications.
2. **Performance** — a structured search for optimizations that close the gap to the
   Accelerate-backed built-in pivoted QR (`dgeqp3`), under the constraint that the algorithm
   never changes and Julia/MATLAB always produce the same output for the same input.

Current gap (geomean relative time, BSQR/baseline, from the committed publication runs):

| | square | short-wide | rinv variant |
|---|---|---|---|
| Julia | 1.53–1.63 | 1.71–1.86 | 1.45–1.67 |
| MATLAB (MEX) | 1.65–1.89 | 1.90–2.19 | 1.64–1.78 |

## Definitions (what "correct" and "same output" mean)

- **Algorithmic correctness**: at every step the implementation selects
  `argmin_j (1 + ||w_j||^2) / ||a_j^(i)||^2` with first-minimum tie-breaking, where
  `w_j = R11^{-1} R12(:,j)`, and maintains `W`, `||w_j||^2`, and the running column norms by
  the recurrences of §2 of the writeup. Any reorganization that computes the same quantities
  in exact arithmetic (different rounding is acceptable) is a *mathematically equivalent*
  modification and is allowed, but must be classified (see change taxonomy below).
- **Cross-implementation parity**: for the same input, Julia, m-file, and MEX must produce
  (a) the identical pivot sequence, and (b) `R`, `tau`-equivalent `Q`, and `rinv_r12` equal to
  within `c · eps · ||A||`-scaled tolerances. Bitwise equality of floats is not required
  (different BLAS runtimes), but pivot sequences must match exactly on the curated parity
  suite. Near-ties in the criterion are the one place FP rounding can legitimately flip a
  pivot; the parity suite therefore filters/flags cases where the runner-up criterion gap is
  within FP noise, and such cases are validated by invariants instead of by sequence equality.

### Sanctioned deviations from the literal Algorithm 1 (to be confirmed and documented)

| # | Deviation | Where | Why equivalent |
|---|---|---|---|
| 1 | Norm downdate `s_j -= alpha_j^2` with LAPACK-style recompute guard `s_j <= s_ref_j * tol` | all three kernels | Remark "Column norm downdating" in the writeup; recomputation returns the exact same mathematical quantity |
| 2 | Householder sign convention (`beta = -copysign(hypot, alpha)`, LAPACK `dlarfg`) so `rho_i` may be negative vs. the writeup's `rho_i = ||·|| ≥ 0` | all three | criterion uses `rho^2`; `Q` absorbs the sign; `R` diagonal sign is a convention |
| 3 | Clamping `wnorm2` and `s` at zero after downdates | all three | guards against `-eps`-level negatives of a nonnegative quantity |
| 4 | `s(i) = beta_i^2` refresh of the pivot column after reflection | all three | exact identity `||a_i^(i)||^2 = beta_i^2` post-reflection |
| 5 | BLAS-2 organization of the W update (gemv for dots, ger for the rank-1) | Julia, MEX | Remark "Practical considerations" in the writeup |
| 6 | `rank_stop` early exit (Julia-only option, off by default in `bsqr`) | Julia | extension, not in Algorithm 1; defaults align across languages |
| 7 | Short-wide fastpath (inline W-row materialization) | Julia | pure memory-layout reordering, identical arithmetic |

Each row gets: the exact code location, the equation in the writeup it corresponds to, and the
test that pins it. Deliverable: `docs/VALIDATION.md` (the audit table, kept current).

---

## Part I — Validation plan

### V0. Literature audit (references → writeup → code)

The writeup is our own exposition; the papers are ground truth. Read Bischof (1990) and
Stewart (1990) from `notes/` and confirm:

- the criterion `(1 + ||w_j||^2)/rho_j^2` is Stewart's column-selection adaptation of
  Bischof's incremental condition estimation, with no missing scaling/normalization;
- the `w`-update recurrence (writeup eq. for `w_ell^(i)`) and the `||w||^2` recurrence match
  the papers' formulas;
- the first-step reduction to Businger–Golub (`c_j = 1/s_j`) is consistent;
- note anything the papers specify that the writeup dropped (tie handling, stopping rules,
  scaling safeguards) and decide whether it matters.

Output: a short section in `docs/VALIDATION.md` mapping paper equations → writeup equations →
code symbols (`ws.W`, `ws.wnorm2`, `ws.s`, `beta_vec`, `dots`).

### V1. Executable oracle — **done** (`matlab/tests/oracle_bsqr.m`)

A deliberately naive, literal transcription of Algorithm 1, written in MATLAB for auditability
(test-only, never benchmarked):

- recomputes `w_j = R11^{-1} r_j` from scratch each step by triangular solve (no incremental W);
- recomputes every column norm from scratch each step (no downdating);
- explicit Householder matrices `H = I - 2vv'/(v'v)`; `O(nk^3)`-ish cost is fine — it exists to
  be unambiguous, not fast;
- comments map each block to the numbered lines of Algorithm 1; deviations (LAPACK sign
  convention, first-minimum ties) documented in the header with the equivalence argument.

The oracle is the single point of truth that "code implements the math". Julia validates
against it through the V3 fixture harness rather than reimplementing it.

### V2. Property/invariant tests (per implementation)

Run on a matrix zoo: gaussian, ill-conditioned (graded singular values), orthonormal-rows
(`M M^T = I`, the GKS use case), Kahan matrices, duplicate columns, zero columns, rank-deficient,
extreme column scaling (`1e±150`), `m>n`, `m<n`, `k < min(m,n)`, `k = 0`, `n = 1`, `m = 1`.

1. **Pivot-sequence equality vs. oracle** on all zoo members small enough (`<= ~200×400`),
   with the runner-up-gap filter for near-ties (if `|c_best − c_second| <= tol_c`, accept either
   and assert invariants instead).
2. **W invariant**: at completion (and, in a debug mode, after each step) check
   `W(1:i, i+1:n) ≈ R11^{-1} R12` computed from the materialized factors by triangular solve.
3. **wnorm2 invariant**: `wnorm2[j] ≈ ||W(1:i, j)||^2`.
4. **Running-norm invariant**: `s[j] ≈ ||A^(i)(i+1:m, j)||^2` (tail norms of the reflected
   matrix), exercised at both `norm_recomp_tol = sqrt(eps)` and `0` (downdate-only) and `1`
   (recompute always) — all three must produce the same pivot sequence on the non-adversarial
   zoo (they are the same mathematical quantity).
5. **frob trace** (Julia): `sum(frob_inv_trace) ≈ ||R11^{-1}||_F^2` computed explicitly
   (Remark "Tracking the full inverse norm").
6. **Osinsky-bound canary**: for orthonormal-rows inputs, assert
   `||V_J^{-1}||_F <= sqrt(k(n−k+1))` and `||V_J^{-1}||_2 <= sqrt(1 + k(n−k))`. A correct
   implementation can never violate these; a subtly wrong criterion eventually will.
7. **Reconstruction/orthogonality**: `||A[:,p] − Q R|| <= c·eps·||A||`, `||Q^T Q − I|| <= c·eps`
   (largely covered by existing tests; extend to the zoo).

MATLAB covers 1 directly against the in-language oracle (`matlab/tests/test_oracle_parity.m`,
done for mfile + MEX, including 2, 5, 6, 7); Julia covers 1 by consuming fixtures (below) rather
than reimplementing the oracle.

### V3. Cross-implementation parity harness

- A Julia script generates the parity fixture set: input matrices plus expected outputs
  (pivot sequence, `R`, `rinv_r12`, `ksteps`) written to a `parity/` fixtures directory as
  CSV (sizes kept modest, e.g. ≤ 512×1024, so this stays fast).
- A MATLAB test (`tests/test_parity_fixtures.m`) runs both `mfile` and `mex` backends on each
  fixture and asserts: identical pivot sequence, factor agreement to scaled tolerance.
- The existing `testMexParityWhenAvailable` (mfile vs MEX, random inputs) stays; this adds the
  Julia↔MATLAB leg, which currently does not exist.
- Fixture matrices are screened at generation time for criterion-gap robustness (runner-up gap
  at every step > 1e3·eps·scale), so exact pivot equality is a fair demand across BLAS runtimes.

### V4. Per-step trace parity (debug instrumentation)

Julia already supports `pivot_history` and kernel stats. Add an option-gated equivalent to the
MEX (and m-file) kernels that records per step: chosen pivot, criterion value, `rho_i`, and
which columns triggered norm recomputes. Validation builds diff the traces across all three
implementations on the fixture suite — this localizes any divergence to the exact step and
quantity rather than discovering it post-hoc in `R`. Keep the instrumentation compiled out /
branch-predicted away in the hot path (it must not contaminate benchmarks).

### V5. Householder construction equivalence

Julia hand-rolls `_householder!` (hypot-based); MEX calls LAPACK `dlarfg`. These agree in the
common case but `dlarfg` has additional rescaling safeguards near under/overflow. Either:

- (a) switch Julia to `LAPACK.larfg!` so reflector construction is identical by construction, or
- (b) keep the hand-rolled version and add equivalence tests across extreme scales
  (`||x|| ~ 1e±300`, subnormal tails).

(a) is preferred: it removes a whole class of divergence and is performance-neutral.

### Exit criteria for Part I

- `docs/VALIDATION.md` audit table complete; every deviation classified and pinned by a test.
- Oracle parity, invariants, and Osinsky canaries green in both languages on the zoo.
- Cross-language fixture parity green (exact pivot sequences, toleranced factors).
- Trace parity green on the fixture suite.

---

## Part II — Performance optimization search

### The constraint that shapes everything: the flop floor

BSQR does strictly more arithmetic than Businger–Golub/`dgeqp3`:

- **Square (m = n = k)**: Householder work ≈ (4/3)n³; W maintenance (gemv+ger per step)
  ≈ (2/3)n³ → ≈ 2n³ total vs. (4/3)n³ for `dgeqp3` → **flop ratio ≈ 1.5**. The observed
  square-regime ratio (1.53–1.63 in Julia) is already *at* the flop floor — meaning the extra
  time is the extra math, not implementation slop.
- **Short-wide (k×n, k ≪ n)**: Householder ≈ 2nk²; W maintenance ≈ 2nk² → **flop ratio ≈ 2**.
  Observed 1.71–1.86, slightly *better* than the naive ratio (Q-materialization and pivoting
  overheads dilute it).
- **rinv variant**: baseline pays an extra `trsm`; BSQR gets `R11^{-1}R12` for free from W.
  This is the regime where true parity is most plausible (observed 1.45–1.67).

Consequence: **parity at 1.0 is not reachable for the plain variant by micro-optimization
alone.** The realistic levers are (a) making BSQR's flops run at higher throughput than
`dgeqp3`'s (blocking → BLAS-3), (b) cutting memory passes (the kernel is BLAS-2/memory-bound),
and (c) the timed-path costs outside the kernel. The plan therefore starts by *measuring the
achievable floor* so "as close as possible" becomes a number.

### P0. Measurement foundation (do this first, ~no code risk)

1. **Per-phase breakdown** across the publication grid using the existing `BSKernelStats`
   (pivot scan / householder / reflector-apply / W-update / norm-downdate). Profile the MEX
   externally (Instruments) or port the stats behind a compile flag. Output: a table of where
   the 1.5–1.9× actually lives, per regime.
2. **Unblocked-LAPACK control**: time LAPACK's *unblocked* CPQR (`dgeqpf`, or `dlaqp2` via
   `ccall`) on the same grid. This splits the gap into "blocking gap" (dgeqp3 vs dgeqpf) and
   "BSQR-extra-work gap" (BSQR vs dgeqpf). If BSQR ≈ dgeqpf, micro-optimization is exhausted
   and only P3 (blocking) can move the number. This single experiment directs all later effort.
3. **Flop/bandwidth model**: formalize the counts above per (m, n, k), including the timed
   Q-materialization and output marshaling, and publish the per-cell theoretical floor next to
   the measured ratios in the benchmark summary.
4. **Safeguard cost**: measure `norm_recomp_tol ∈ {0, sqrt(eps)}` deltas to bound what the
   stability guard costs (expected: negligible; verify).
5. **MEX overhead audit**: time `bsqr_mex` against the Julia kernel on identical inputs.
   The MATLAB short-wide gap (2.19) is notably worse than Julia's (1.83); quantify how much is
   marshaling (input copy, finite check, `W` zero-fill, mxArray output copies) vs. kernel.

### P1. Timed-path hygiene (low risk, kernel math untouched)

1. **`dorgqr` instead of `dormqr`-on-identity** for Q materialization — in Julia
   `_explicit_q` / `_materialize_bsqr_qrp`, and MEX `build_q`. The baseline's `qr` uses
   `dorgqr`, which exploits the identity structure; this is fair-comparison hygiene and a
   plausible few-percent win on square cases.
2. **MEX: fuse the finite check with the initial column-norm pass** (one read of A instead of
   two), and **drop the eager `W.assign(k·n, 0)`** if the kernel provably writes row i before
   any read (Julia already documents this contract in `_require_workspace`).
3. **Julia: fuse the next step's pivot scan into the update loop** — the criterion argmin for
   step i+1 can be accumulated during step i's wnorm2/s downdate pass over the same data
   (recomputed columns re-checked afterwards). Same arithmetic, one fewer O(n−i) pass per step.
4. **Avoid the strided beta-row re-read**: `alpha_j` (row i of the reflected trailing block) is
   computable from the larf work vector (`alpha_new = alpha_old − tau·work_j`) instead of
   re-reading the strided row. Verify bitwise-or-equivalent before adopting.

### P2. Unblocked kernel restructuring (medium risk, classify each change)

1. **Rank-2 fused updates**: per step there are two gemv+ger pairs (trailing A, W block).
   Explore expressing each as a single small-`k` `dgemm` (or grouping the two gers); Accelerate
   may schedule these better. Measure; keep only if it wins.
2. **W layout experiments**: W is `kmax×n` column-major; the per-step gemv/ger work on an
   `(i−1)×(n−i)` submatrix with lda = kmax. Try transposed or tiled layouts for cache behavior.
3. **Division-free pivot comparison** (`(1+wn_a)·s_b < (1+wn_b)·s_a`): removes n−i divisions
   per step but changes comparison rounding → this is an *algorithm-visible* change (tie flips)
   and must be adopted in all three implementations simultaneously or rejected.

### P3. Panel/blocked BSQR (the big bet — only path past the BLAS-2 ceiling)

`dgeqp3` wins by deferring trailing updates: within a panel it fully updates only the pivot
column and the current row (via the accumulated block reflectors), then applies the whole panel
to the trailing matrix as one `dgemm` (`dlaqps`). The same restructuring looks feasible for
BSQR because pivot selection needs only `s` (downdatable from row i) and `wnorm2`:

1. **Derive on paper first** (appendix to the writeup in `notes/`): a panel variant maintaining
   compact-WY reflectors, deferred W rank-1 updates (accumulate the `w*·beta^T` outer products
   over the panel), and per-step `dots = W_rem^T w*` computed from the deferred representation.
   Prove exact-arithmetic equivalence to Algorithm 1 step by step.
2. **Prototype in Julia**, behind an env flag; validate pivot-sequence parity vs. the unblocked
   kernel on the V-suite; benchmark on the publication grid.
3. **Decision gate**: adopt only on a clear win (suggested: ≥20% median improvement on the
   regimes that matter, echoing the MEX promotion rule) with zero parity regressions.
4. **If adopted, port the identical blocked algorithm to the MEX.** Blocking changes rounding,
   so per the parity policy it is an algorithm-level change: Julia and MEX must move together.
   The m-file backend stays as the simple unblocked reference (its job is correctness), with
   the documented caveat that near-tie pivots may differ from the blocked kernels; the parity
   fixture suite (tie-free by construction) must still pass across all three.

### P4. Variant- and platform-specific

1. **rinv variant**: with W free, push this comparison toward parity — it is the structurally
   winnable case. Ensure the baseline `trsm` leg is itself optimally implemented (fairness cuts
   both ways) and that no avoidable copies sit on the BSQR rinv extraction path.
2. **Threads**: BLAS-2 kernels barely use them; revisit the `BS_PUB_THREADS` grid only after P3
   (blocking is what unlocks multithreaded gemm).
3. **MATLAB marshaling**: based on the P0.5 audit — write outputs directly into `plhs` buffers,
   skip the pivot-matrix materialization in benchmarks (vector form is already the benchmark
   contract), and consider `mxCreateUninitNumericMatrix` for outputs.

### Methodology: how every candidate lands

1. Implement behind a flag (env var or option) where feasible.
2. Classify the change:
   - **Schedule/layout** (identical per-step arithmetic): may land in one implementation alone.
   - **Algorithm-visible** (different rounding/regrouping, e.g. blocking, division-free
     compare): must land in Julia + MEX together, parity suite green, `docs/VALIDATION.md`
     updated.
3. Run the full Part-I validation suite (oracle, invariants, fixtures, traces).
4. Run the smoke benchmark, then the perf gate (`check_publication_perf_gate.*`) against the
   committed baseline CSV.
5. On a win: regenerate publication artifacts (both languages if both changed), commit
   artifacts per repo convention.

### Risks and open questions

- **Near-tie pivot flips** are the fundamental tension between "optimizations may diverge" and
  "same output for same input". The fixture-suite gap filter is the mitigation; the policy
  (exact sequence parity on tie-free inputs, invariant-based acceptance otherwise) should be
  agreed before P2/P3 work starts.
- **The square regime may already be at its floor** — if P0.2 confirms BSQR ≈ unblocked LAPACK
  CPQR, the honest publication story is the flop-ratio floor plus the blocked-kernel result,
  not further micro-optimization.
- **Blocked-BSQR derivation may not close** (the `dots` maintenance under deferred W updates is
  the part most likely to resist a clean compact form). Timebox the paper derivation before
  committing to implementation.

### Suggested sequencing

1. V0–V2 (oracle + invariants) — establishes the safety net everything else relies on.
2. V3–V5 (parity harness, traces, larfg unification) in parallel with P0 (measurement).
3. P1 (hygiene wins) — cheap, lands while P3's derivation is being worked out on paper.
4. P2 selectively, guided by the P0 phase breakdown.
5. P3 derivation → prototype → gate → port.
6. P4 cleanup; final artifact regeneration; update the writeup/paper with the floor analysis.
