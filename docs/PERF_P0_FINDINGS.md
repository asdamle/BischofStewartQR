# P0 Measurement Findings

Phase P0 of `docs/VALIDATION_AND_PERF_PLAN.md` (Part II). Reproduce with
`julia --project=julia julia/benchmark/perf_p0_measure.jl` (gaussian family, BLAS pinned to
1 thread, Apple Accelerate via LBT, median of 10). Measured 2026-06-11 on Apple Silicon
(macOS 24.6.0).

## P0.2 — Unblocked-LAPACK control (the decisive experiment)

`dgeqpf` is LAPACK's unblocked pivoted QR; `dgeqp3` its blocked replacement. All times are
relative to `dgeqp3` on the fair materialized-`Q,R,p` path; "flop floor" is the predicted
ratio assuming equal flop throughput: `(HH + W + Qmat)/(HH + Qmat)`.

| case | dgeqp3 (s) | bsqr | bsqr_tol0 | dgeqpf | flop floor |
|---|---:|---:|---:|---:|---:|
| square 256x256 | 0.002248 | 1.546 | 1.500 | 1.117 | 1.248 |
| square 384x384 | 0.004691 | 1.886 | 1.847 | 1.259 | 1.249 |
| square 512x512 | 0.009341 | 1.985 | 1.994 | 1.282 | 1.249 |
| short_wide 64x640 | 0.0006206 | 1.813 | 1.792 | 0.920 | 1.902 |
| short_wide 128x1024 | 0.002736 | 1.887 | 1.862 | 0.982 | 1.879 |
| short_wide 256x1024 | 0.006791 | 2.034 | 2.055 | 1.110 | 1.768 |

Conclusions:

1. **Short-wide (the GKS focus regime): BSQR is already at its flop floor.** Measured 1.887
   vs. predicted 1.879 at 128×1024. There is nothing left for micro-optimization to recover;
   only reducing the cost of the W work (BLAS-3 reorganization) or accepting the floor moves
   this number.
2. **Short-wide: blocking buys the *baseline* nothing either** — `dgeqpf` ≈ `dgeqp3`
   (0.92–1.11). The panel width is capped by the small `m`, so `dgeqp3`'s gemm advantage
   vanishes. This regime is a fair fight between unblocked kernels, and BSQR's deficit is
   exactly its extra flops.
3. **Square: two separable gaps.** The blocking gap (`dgeqpf`/`dgeqp3`) is 1.12–1.28 and grows
   with size; the BSQR-extra-work gap (`bsqr`/`dgeqpf`) is 1.38–1.55 against a kernel-flop
   prediction of ~1.25, leaving ~10–25% unblocked implementation headroom (W-update memory
   traffic and the Q-materialization asymmetry below).
4. **Caveat:** the timed path materializes Q differently on the two sides — BSQR uses
   `ormqr`-on-identity (full m×m), the baseline `orgqr` (thin). For `m ≤ n` the shapes agree
   but `orgqr` does ~1/3 fewer flops, so part of the measured BSQR deficit is materialization
   asymmetry, not kernel. Fixing this is P1.1 and is also fairness hygiene.

## P0.1 — Kernel phase breakdown (instrumented single run, kernel only)

| case | kernel (s) | pivot scan | householder | reflector apply | W update | norm downdate | other | recomputes |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| square 256x256 | 0.00272 | 0.8% | 1.1% | 57.2% | 38.2% | 0.0% | 2.7% | 0 |
| square 384x384 | 0.00712 | 0.7% | 0.8% | 56.9% | 39.8% | 0.0% | 1.8% | 0 |
| square 512x512 | 0.0152 | 0.5% | 0.6% | 53.8% | 43.2% | 0.0% | 1.8% | 0 |
| short_wide 64x640 | 0.00108 | 2.2% | 0.4% | 44.3% | 51.2% | 0.1% | 1.8% | 577 |
| short_wide 128x1024 | 0.00496 | 1.6% | 0.2% | 44.1% | 52.9% | 0.0% | 1.3% | 899 |
| short_wide 256x1024 | 0.0128 | 1.1% | 0.2% | 44.7% | 52.9% | 0.0% | 1.1% | 769 |

The kernel is two BLAS-2 walls: the trailing reflector application (44–57%) and the W update
(38–53%, dominant in short-wide). Pivot scan, norm downdates, and everything else are ≤3%
combined — the P1.3/P1.4 micro-optimizations (fused scans, beta-row reuse) target at most a
few percent and are **not worth doing**.

## P0.4 — Safeguard cost

`norm_recomp_tol = 0` vs. the default `sqrt(eps)` differs by ≤3% (within noise), even on
short-wide cases firing 577–899 recomputes. The stability safeguard is free; keep the default.

## P0.5 — MEX overhead audit (timeit probes, default MATLAB threading)

| case | qr econ (s) | mex QRp | mex R-only | mex QRp nocheck | Q+marshal share | finite-check share |
|---|---:|---:|---:|---:|---:|---:|
| square 384x384 | 0.003926 | 0.01418 | 0.01359 | 0.01524 | 4.2% | noise |
| short_wide 128x1024 | 0.005098 | 0.00502 | 0.004949 | 0.004979 | 1.4% | 0.8% |

Robust takeaways: output marshaling, Q formation, and the finite check are all small (≤4%) —
the MEX gap is in the kernel, not the mxArray plumbing. The absolute mex-vs-`qr` ratios from
this quick probe disagree with the publication numbers in both directions (MATLAB's default
multithreading and `timeit` methodology differ from the controlled harness); treat the *shares*
as the finding and re-derive ratios only from the publication runners.

## P1.1 result (orgqr materialization, landed)

Post-fix rerun: square 512 improved 1.985 → 1.929, 384 1.886 → 1.865, 256 within noise;
short-wide unchanged (Q materialization is a small share there). The remaining ratios are the
kernel story; see `docs/P3_BLOCKED_BSQR.md` for the attack on the two BLAS-2 walls.

## Revised priorities for the optimization search

1. **P1.1 (do now): materialize Q via `orgqr` on both BSQR paths** (Julia
   `_materialize_bsqr_qrp`/`_explicit_q`, MEX `build_q`). Closes the measured materialization
   asymmetry; estimated ~5% of the square timed path, and it is fairness hygiene.
2. **P3 (the main event): panel/blocked BSQR.** It attacks both BLAS-2 walls at once — the
   trailing reflector application (the square-regime blocking gap that `dgeqp3` already
   exploits) and the W update (the short-wide wall, 52%). Short-wide can only improve through
   the W term; square can plausibly reach ~1.25 (its flop floor) or better if blocked W work
   runs at gemm throughput.
3. **P2 selectively:** W-layout experiments may recover part of the 10–25% unblocked headroom
   on square (W's `lda = kmax` strided access). Skip P1.3/P1.4 and the division-free
   comparison (≤3% targets, and the latter is an algorithm-visible change).
4. **Skip:** safeguard tuning (free), MEX marshaling work (≤4%), finite-check fusion (noise).
