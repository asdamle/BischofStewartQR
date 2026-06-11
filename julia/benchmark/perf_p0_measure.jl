#!/usr/bin/env julia
#
# P0 measurement foundation (docs/VALIDATION_AND_PERF_PLAN.md, Part II).
#
# For each grid case this script reports, against the dgeqp3 baseline:
#   * bsqr_full        - the production kernel, fair materialized Q,R,p
#   * bsqr_tol0        - same with norm_recomp_tol = 0 (safeguard cost, P0.4)
#   * dgeqpf           - LAPACK's UNBLOCKED pivoted QR (the control, P0.2):
#                        the dgeqp3/dgeqpf gap is the "blocking gap"; the
#                        bsqr/dgeqpf gap is the "BSQR-extra-work gap"
#   * flop-floor model - predicted timed-path ratio assuming equal flop
#                        throughput (P0.3): (HH + W + Qmat)/(HH + Qmat)
#   * per-phase kernel breakdown via BSKernelStats (P0.1)
#
# Defaults are single-threaded BLAS on the gaussian family. Knobs:
#   BS_P0_SAMPLES (default 10), BS_P0_WARMUP (default 2)

using LinearAlgebra
using LinearAlgebra: BLAS, BlasInt
using Printf
using Random

using BSPivotQR

include("bench_common.jl")
using .BenchCommon
using .BenchCommon: make_matrix, bench_trial_basic, run_bsqr_fair, run_qr_fair

const CASES = [
    (regime = "square", m = 256, n = 256),
    (regime = "square", m = 384, n = 384),
    (regime = "square", m = 512, n = 512),
    (regime = "short_wide", m = 64, n = 640),
    (regime = "short_wide", m = 128, n = 1024),
    (regime = "short_wide", m = 256, n = 1024),
]

function run_dgeqpf_fair(A::Matrix{Float64})
    m, n = size(A)
    k = min(m, n)
    Awork = copy(A)
    jpvt = zeros(BlasInt, n)
    tau = zeros(k)
    work = zeros(3 * n)
    info = Ref{BlasInt}(0)
    ccall(
        (BLAS.@blasfunc(dgeqpf_), BLAS.libblastrampoline),
        Cvoid,
        (Ref{BlasInt}, Ref{BlasInt}, Ptr{Float64}, Ref{BlasInt},
         Ptr{BlasInt}, Ptr{Float64}, Ptr{Float64}, Ptr{BlasInt}),
        m, n, Awork, m, jpvt, tau, work, info,
    )
    info[] == 0 || error("dgeqpf failed with info=$(info[])")
    R = triu(Awork[1:k, :])
    Q = LAPACK.orgqr!(Awork[:, 1:k], tau)
    return (Q = Q, R = R, p = Int.(jpvt))
end

# Flop model (P0.3). Householder reduction and the W maintenance are exact
# sums of the per-step BLAS-2 costs; Q materialization (orgqr, m-by-k from k
# reflectors) is paid identically by every method and dilutes the ratio.
hh_flops(m, n, k) = sum(4.0 * (m - i) * (n - i) for i in 1:k)
w_flops(n, k) = sum(4.0 * (i - 1) * (n - i) for i in 1:k)
qmat_flops(m, k) = 2.0 * m * k^2 - (2.0 / 3.0) * k^3

function flop_floor(m, n, k)
    hh = hh_flops(m, n, k)
    w = w_flops(n, k)
    qm = qmat_flops(m, k)
    return (hh + w + qm) / (hh + qm)
end

function kernel_phase_stats(A::Matrix{Float64}, k::Int)
    m, n = size(A)
    run() = begin
        Awork = copy(A)
        ws = BSPivotQR.BSWorkspace(m, n, k)
        tau = zeros(k)
        jpvt = collect(1:n)
        stats = BSPivotQR.BSKernelStats()
        t0 = time_ns()
        BSPivotQR._bsqr_kernel!(Awork, tau, jpvt, ws, k;
            rank_stop = false, kernel_stats = stats)
        wall = time_ns() - t0
        (stats, wall)
    end
    run()  # warmup
    stats, wall = run()
    total = Float64(wall)
    phases = (
        pivot = stats.pivot_select_ns / total,
        hh = stats.householder_ns / total,
        reflect = stats.apply_reflector_ns / total,
        wupd = stats.w_update_ns / total,
        downdate = stats.norm_downdate_ns / total,
    )
    other = 1.0 - sum(values(phases))
    return phases, other, stats.recompute_count, total * 1e-9
end

function main()
    check_backend()
    BLAS.set_num_threads(1)
    println("BLAS threads: ", BLAS.get_num_threads())

    samples = parse_env_int("BS_P0_SAMPLES", 10)
    warmup = parse_env_int("BS_P0_WARMUP", 2)
    tol = DEFAULT_NORM_RECOMP_TOL
    rng = MersenneTwister(20260310)

    println("\n## Timings (relative to dgeqp3; median of $samples, gaussian, 1 thread)\n")
    println("| case | dgeqp3 (s) | bsqr | bsqr_tol0 | dgeqpf | flop floor |")
    println("|---|---:|---:|---:|---:|---:|")
    phase_rows = String[]
    for c in CASES
        A = make_matrix(:gaussian, c.m, c.n, rng)
        k = min(c.m, c.n)

        _, t_qp3, _ = bench_trial_basic(() -> run_qr_fair(A); warmup = warmup, samples = samples)
        _, t_bs, _ = bench_trial_basic(() -> run_bsqr_fair(A, k, tol); warmup = warmup, samples = samples)
        _, t_bs0, _ = bench_trial_basic(() -> run_bsqr_fair(A, k, 0.0); warmup = warmup, samples = samples)
        _, t_qpf, _ = bench_trial_basic(() -> run_dgeqpf_fair(A); warmup = warmup, samples = samples)

        label = @sprintf("%s %dx%d", c.regime, c.m, c.n)
        @printf("| %s | %.4g | %.3f | %.3f | %.3f | %.3f |\n",
            label, t_qp3, t_bs / t_qp3, t_bs0 / t_qp3, t_qpf / t_qp3,
            flop_floor(c.m, c.n, k))

        phases, other, nrecomp, kt = kernel_phase_stats(A, k)
        push!(phase_rows, @sprintf(
            "| %s | %.3g | %.1f%% | %.1f%% | %.1f%% | %.1f%% | %.1f%% | %.1f%% | %d |",
            label, kt, 100 * phases.pivot, 100 * phases.hh, 100 * phases.reflect,
            100 * phases.wupd, 100 * phases.downdate, 100 * other, nrecomp))
    end

    println("\n## BSQR kernel phase breakdown (single instrumented run, kernel only)\n")
    println("| case | kernel (s) | pivot scan | householder | reflector apply | W update | norm downdate | other | recomputes |")
    println("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    foreach(println, phase_rows)
end

main()
