# Tests and Benchmarks

This document describes how to run, interpret, and extend the validation and performance workflows for `BSPivotQR`.

## Prerequisites

From repository root:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

Additional packages:

```bash
julia --project=. -e 'using Pkg; Pkg.add("AppleAccelerate")'
python3 -m pip install matplotlib
```

`AppleAccelerate` is optional (macOS backend preference). `matplotlib` is required for plotting scripts and for tier validation runners because they invoke plotting steps.

## Test Suite

Run all tests:

```bash
julia --project=. test/runtests.jl
```

The suite is defined in `test/test_bsqr.jl` and currently covers:

1. `Bischof-Stewart CPQR correctness`
   - Residual and orthogonality checks for `m<n`, `m=n`, and `m>n`.
   - Permutation validity.
   - Reconstruction consistency.
2. `Early stop and edge cases`
   - `k=0`, partial `k`, and full `k`.
   - Rank-deficient/degenerate behavior.
   - Deterministic tie handling.
3. `Criterion-consistent pivot sequence`
   - Confirms selected pivots match direct per-step `argmin` of the Bischof-Stewart criterion on small problems.
4. `Comparison with Julia qr(ColumnNorm())`
   - Compares numerical quality (residual/orthogonality) against LAPACK `DGEQP3` via Julia `qr`.
5. `Optional R11^{-1}R12 return`
   - Verifies returned `R11^{-1}R12` equals direct triangular solve from final `R`.
6. `Preallocated API and BLAS thread control`
   - Verifies low-allocation kernel API consistency and BLAS thread restoration.
7. `Rank-stop policy`
   - Verifies default rank-stop behavior on rank-deficient inputs and forced full-step behavior with `rank_stop=false`.

## Benchmark Suite

Primary benchmark script:

- `benchmark/bench_cpqr.jl`
- systematic sweep script:
  - `benchmark/bench_cpqr_sweep.jl`
- regime scaling script:
  - `benchmark/bench_cpqr_regimes.jl`
- profiling script:
  - `benchmark/profile_bsqr_breakdown.jl`

Run benchmark:

```bash
BS_USE_ACCELERATE=1 BS_REQUIRE_ACCELERATE=0 julia --project=. benchmark/bench_cpqr.jl
```

Quick smoke run:

```bash
BS_QUICK=1 BS_USE_ACCELERATE=1 BS_REQUIRE_ACCELERATE=0 julia --project=. benchmark/bench_cpqr.jl
```

### Fairness policy (important)

`bench_cpqr.jl` enforces equal footing for `bsqr!` and `qr(..., ColumnNorm())`:

1. Both benchmarked on fresh `copy(A)` per trial.
2. `bsqr!` benchmark path disables bsqr-only extras:
   - no reusable workspace,
   - no inverse-Frobenius trace tracking,
   - no `R11^{-1}R12` return.
3. `bsqr!` benchmark path always runs full-step (`k=min(m,n)`) with `rank_stop=false`.
4. BLAS thread setting is global and shared by both methods.

### Canonical benchmark policy

This is the single source of truth for benchmark behavior:

1. Benchmarks always compare full-step factorizations (`k=min(m,n)`).
2. Benchmark runs force `rank_stop=false` for BSQR so work is matched to `qr(..., ColumnNorm())`.
3. Both methods benchmark on fresh `copy(A)` allocations (no BS-only workspace reuse in baseline comparisons).
4. Optional BSQR-only outputs/features are disabled in timing paths:
   - `track_inverse_frob=false`
   - `return_rinv_r12=false`
5. Orthonormal-row regime runs are restricted to short-wide growth (`fixed_m_vary_n` only).

## Benchmark Environment Variables

1. `BS_QUICK`
   - `1`: small size set + fewer samples.
   - default `0`: full size set.
2. `BS_USE_ACCELERATE`
   - `1`: attempt `using AppleAccelerate` on macOS.
   - default `1`.
3. `BS_REQUIRE_ACCELERATE`
   - `1`: fail if Accelerate backend is not active.
   - default `0`.
4. `BS_BLAS_THREADS`
   - set BLAS thread count for benchmark process (applies to both methods).
   - Example: `BS_BLAS_THREADS=8`.
5. `BS_SWEEP_SAMPLES`
   - number of timing samples per method/case for `bench_cpqr_sweep.jl`.
   - defaults: `12` in quick mode, `24` otherwise.
6. `BS_SWEEP_WARMUP`
   - warmup evaluations per method/case for `bench_cpqr_sweep.jl`.
   - defaults: `1` in quick mode, `2` otherwise.
7. `BS_REGIME_QUICK`
   - `1`: smaller regime grids for faster smoke runs.
   - default `0`: larger regime grids.
8. `BS_REGIME_FIXED_NS`, `BS_REGIME_MS`
   - fixed-`n` and varying-`m` grids for `bench_cpqr_regimes.jl`.
9. `BS_REGIME_FIXED_MS`, `BS_REGIME_NS`
   - fixed-`m` and varying-`n` grids for `bench_cpqr_regimes.jl`.
10. `BS_REGIME_SAMPLES`, `BS_REGIME_WARMUP`
   - timing sample and warmup counts for `bench_cpqr_regimes.jl`.
11. `BS_NORM_RECOMP_TOL`
   - controls bsqr partial-norm recomputation trigger threshold used in all benchmark scripts.
   - default: `sqrt(eps(Float64))` (LAPACK-style conservative setting).
   - lower values trigger fewer exact norm recomputations.
12. `BS_PROFILE_QUICK`, `BS_PROFILE_REPS`, `BS_PROFILE_SHORT_M`, `BS_PROFILE_SHORT_N`, `BS_PROFILE_SQUARE_N`, `BS_PROFILE_FAMILIES`, `BS_PROFILE_TOLS`
   - controls profile case sizes/families/repetitions/tolerance settings for `profile_bsqr_breakdown.jl`.
13. `BS_BASELINE_TAG`
   - optional tag suffix for `benchmark/freeze_baseline.jl`; output folder is `benchmark/results/baseline_<tag>`.
14. `BS_BASELINE_DIR`
   - baseline path used by validation runners for automatic baseline-vs-candidate comparison.
15. `BS_SLOWDOWN_THRESH`
   - slowdown guardrail threshold for `benchmark/compare_results.jl` (default `0.02`, i.e., 2%).
16. `BS_ASSESS_RUN_STAGE1`, `BS_ASSESS_RUN_STAGE2`
   - controls triage/deep phases for `benchmark/assess_status.jl` (default `1` for both).
17. `BS_ASSESS_DEEP_QUICK`
   - `1`: reduced deep-assessment grid for faster turnaround.
   - default `1`; set `0` for expanded matrix/regime grid.
18. `BS_ASSESS_WARMUP`, `BS_ASSESS_SAMPLES`, `BS_ASSESS_PROFILE_REPS`
   - warmup/sample/profile repetition controls for Stage-2 deep assessment.
19. `BS_ASSESS_TOLS`
   - comma-separated tolerance sweep for deep assessment.
   - default: `sqrt(eps),1e-10,1e-12,0.0`.
20. `BS_ASSESS_FAMILIES`
   - comma-separated family subset for deep assessment (`gaussian`, `ill_conditioned`, `orthonormal_rows`).
21. `BS_REGIME_FAMILIES`
   - comma-separated family subset for `bench_cpqr_regimes.jl` (`gaussian`, `ill_conditioned`, `orthonormal_rows`).
22. `BS_SHORT_WIDE_FASTPATH`
   - controls BSQR short-wide kernel fast path (`1` enabled, `0` disabled).
   - default: `1`.

Example reproducible run:

```bash
BS_BLAS_THREADS=8 BS_USE_ACCELERATE=1 BS_REQUIRE_ACCELERATE=1 julia --project=. benchmark/bench_cpqr.jl
```

Systematic sweep run:

```bash
BS_BLAS_THREADS=8 BS_USE_ACCELERATE=1 BS_REQUIRE_ACCELERATE=1 julia --project=. benchmark/bench_cpqr_sweep.jl
```

Systematic sweep quick mode:

```bash
BS_SWEEP_QUICK=1 julia --project=. benchmark/bench_cpqr_sweep.jl
```

Default sweep grids:

1. quick mode `n` grid: `64,128,256,512,1024`
2. full mode `n` grid: `64,128,256,384,512,768,1024,1536,2048`
3. default families include `gaussian` and `ill_conditioned`.

Systematic sweep custom grids:

```bash
BS_SWEEP_NS=64,128,256,512 BS_SWEEP_ASPECTS=0.25,0.5,1,2,4 BS_SWEEP_FAMILIES=gaussian,ill_conditioned julia --project=. benchmark/bench_cpqr_sweep.jl
```

Regime benchmark quick mode:

```bash
BS_REGIME_QUICK=1 julia --project=. benchmark/bench_cpqr_regimes.jl
```

Regime defaults include `orthonormal_rows` to cover short-wide orthonormal-row behavior. For `orthonormal_rows`, only fixed-`m`/increasing-`n` cases are run.

Regime benchmark custom grids:

```bash
BS_REGIME_FIXED_NS=128,256 BS_REGIME_MS=32,64,128,256,512,1024 BS_REGIME_FIXED_MS=32,64 BS_REGIME_NS=64,128,256,512,1024,2048 julia --project=. benchmark/bench_cpqr_regimes.jl
```

Kernel breakdown profiling:

```bash
BS_PROFILE_QUICK=1 julia --project=. benchmark/profile_bsqr_breakdown.jl
```

Validation runners:

```bash
julia --project=. benchmark/run_validation_tier1.jl
julia --project=. benchmark/run_validation_tier2.jl
```

Note: tier validation runners execute plotting scripts; ensure `python3 -m pip install matplotlib` has been run.

Baseline snapshot + comparison:

```bash
julia --project=. benchmark/freeze_baseline.jl
julia --project=. benchmark/compare_results.jl benchmark/results/baseline_YYYYMMDD_HHMMSS
```

Comprehensive status + deep assessment:

```bash
julia --project=. benchmark/assess_status.jl
julia --project=. benchmark/analyze_assessment.jl
```

Skip triage reruns and execute only Stage-2 deep assessment:

```bash
BS_ASSESS_RUN_STAGE1=0 julia --project=. benchmark/assess_status.jl
```

Run expanded deep-assessment matrix (larger grids/samples):

```bash
BS_ASSESS_DEEP_QUICK=0 julia --project=. benchmark/assess_status.jl
```

## Benchmark Outputs

Outputs are written under:

- `benchmark/results/`

Files:

1. `timings.csv`
   - columns:
     - `family,m,n,k,method,tmin_s,tmed_s,alloc_bytes,residual,orthogonality`
2. `summary.md`
   - Markdown summary table with speedup versus `dgeqp3`.
3. `sweep_timings.csv`
   - systematic sweep benchmark data with `aspect` column.
   - columns:
     - `family,aspect,m,n,k,method,tmin_s,tmed_s,tci_low_s,tci_high_s,alloc_bytes,residual,orthogonality`
   - `tci_low_s` and `tci_high_s` are empirical 95% timing interval bounds from benchmark samples.
4. `sweep_summary.md`
   - Markdown summary for systematic sweeps including 95% timing intervals.
5. `regime_timings.csv`
   - regime-scaling benchmark data (fixed `n` vary `m`, fixed `m` vary `n`).
   - columns:
     - `family,regime,fixed_value,var_value,m,n,k,method,tmin_s,tmed_s,tci_low_s,tci_high_s,alloc_bytes,residual,orthogonality`
6. `regime_summary.md`
   - Markdown summary table for regime scaling runs.
7. `profile_breakdown.csv`
   - per-case bsqr kernel phase timing breakdown (`pivot`, `householder`, `apply`, `W-update`, `downdate`) and recompute counts.
8. `profile_breakdown_summary.md`
   - markdown summary of phase percentages and quality metrics.
9. `validation_report.md`
   - baseline-vs-candidate pass/fail summary including 2% slowdown guardrail checks.
10. `guardrail_failures.csv`
   - machine-readable list of guardrail or quality failures.
11. `status_triage_YYYYMMDD_HHMMSS.md`
   - stage-1 triage snapshot with environment fingerprint, speedup distribution summary, sanity checks, and stress-check results.
12. `profile_hotspots_YYYYMMDD_HHMMSS.csv`
   - stage-1 profile export with phase percentages per case (`pivot`, `householder`, `apply`, `W-update`, `downdate`).
13. `assessment_raw.csv`
   - stage-2 additive raw dataset for CI-aware speedups and phase metrics.
14. `assessment_summary.md`
   - significance-aware aggregate summary and quality-risk/hotspot correlation analysis.
15. `optimization_candidates.md`
   - ranked K1-K5 optimization backlog with go/no-go guidance.
16. `perf_assessment_conclusion.md`
   - concise executive conclusion: status, bottlenecks, top opportunities, and stability no-go areas.

Baseline freeze outputs:

1. `baseline_<tag>/metadata.txt`
   - captures Julia version, BLAS backend, BLAS thread count, CPU info, and benchmark env vars.
2. `baseline_<tag>/...`
   - snapshot copies of benchmark CSV/Markdown outputs and generated plot folders.

Interpretation:

1. `speedup > 1` means `bsqr` is faster than `dgeqp3`.
2. `residual` and `orthogonality` should stay near machine precision scale.
3. Compare both timing and quality; treat large speedups with degraded accuracy as unacceptable.

## Plotting Results

Plot script:

- `benchmark/plot_results.jl`
- sweep plot script:
  - `benchmark/plot_sweep_results.jl`
- regime plot script:
  - `benchmark/plot_regime_results.jl`

Run:

```bash
julia --project=. benchmark/plot_results.jl
```

Optional path override:

```bash
julia --project=. benchmark/plot_results.jl benchmark/results/timings.csv benchmark/results/plots
```

Sweep plotting:

```bash
julia --project=. benchmark/plot_sweep_results.jl
```

Sweep plotting with explicit paths:

```bash
julia --project=. benchmark/plot_sweep_results.jl benchmark/results/sweep_timings.csv benchmark/results/sweep_plots
```

Regime plotting:

```bash
julia --project=. benchmark/plot_regime_results.jl
```

Regime plotting with explicit paths:

```bash
julia --project=. benchmark/plot_regime_results.jl benchmark/results/regime_timings.csv benchmark/results/regime_plots
```

Sweep timing-line plots include shaded 95% timing intervals when CI columns are present. Legacy sweep CSV files without CI columns are still supported (no shaded interval in that case).
Sweep timing-line plots are rendered on log-log axes and include an `n^3` reference trend line.
Regime timing plots are log-log, include shaded 95% timing intervals, and include a linear reference trend to highlight fixed-dimension linear scaling.

Output plots (PNG) are written to:

- `benchmark/results/plots`
- `benchmark/results/sweep_plots`
- `benchmark/results/regime_plots`

Notes:

1. Plot scripts always use Python `matplotlib` (no Julia `Plots` backend path).
2. Ensure `python3` and `matplotlib` are available in the runtime environment.

## Correctness Regression Protocol

Tier 1 (every simplification/change):

1. `julia --project=. benchmark/run_validation_tier1.jl`
2. Required outputs:
   - `benchmark/results/timings.csv`
   - `benchmark/results/sweep_timings.csv`
   - `benchmark/results/regime_timings.csv`
   - `benchmark/results/profile_breakdown.csv`
3. If `BS_BASELINE_DIR` is set, require:
   - `benchmark/results/validation_report.md` marked `PASS`
   - no entries in `benchmark/results/guardrail_failures.csv`

Tier 2 (before release/merge):

1. `julia --project=. benchmark/run_validation_tier2.jl`
2. Require all tests pass and all benchmark families complete.
3. If comparing against a frozen baseline, reject if any key-case median slowdown exceeds 2%.

## Extending Tests and Benchmarks

When adding coverage:

1. Add correctness assertions before timing changes.
2. Keep benchmark fairness constraints unchanged.
3. For new matrix families/sizes, update:
   - `benchmark/matrix_generators.jl`
   - `benchmark/bench_cpqr.jl`
4. If adding a bsqr-only optimization mode (e.g., preallocation), benchmark it in a separate method label and do not use it as the baseline `bsqr_full` comparison against `dgeqp3`.
