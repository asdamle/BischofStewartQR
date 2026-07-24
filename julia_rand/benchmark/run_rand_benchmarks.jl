# Fair randomized-vs-deterministic BSQR timing (no R12) -- the Julia port of
# matlab_rand/benchmark/run_rand_benchmarks.m. Validates the implementation:
# it is not part of the publication pipeline.
#
#   julia --project=julia_rand/benchmark julia_rand/benchmark/run_rand_benchmarks.jl
#
# Fairness/measurement discipline (all methods, identical conditions):
#   * check/check_finite = false everywhere -- the O(mn) finiteness scan is
#     identical overhead for all and part of no algorithm.
#   * BenchmarkTools medians; only the factorization call is inside the timed
#     thunk (matrix generation is outside).
#   * Every side materializes its full product set, matching the publication
#     convention:
#       - randomized: [p, Q, R11] -- the "R12 not needed" product; Q is the
#         economy factor via orgqr (O(mk^2), n-independent). No R12.
#         Norm-weighted sampling, running_mean threshold (the defaults).
#       - deterministic BSQR baseline: bsqr(M; k) with Q, R, perm all
#         materialized (R12 arrives as an unavoidable byproduct of the
#         deterministic kernel's O(nk^2) scan -- the work the randomized
#         variant skips).
#       - built-in baseline: qr(M, ColumnNorm()) (LAPACK dgeqp3, a
#         vendor-tuned classical algorithm) with Q, R, p materialized.
#   * A second randomized row adds return_r12 = true -- how much speedup
#     survives when R12 is required.
#
# Configuration (env vars):
#   BS_RAND_SIZES    "32x2000,64x4000,64x8000,128x8000,128x16000" (default)
#   BS_RAND_FAMILIES "gaussian,needle" (default; see rand_test_matrix.jl)
#   BS_RAND_SEED     base seed (default 1)
#   BS_RAND_OUTDIR   output dir (default julia_rand/benchmark/results, git-ignored)
#
# Writes <outdir>/rand_timings.csv and prints the summary table.

using BenchmarkTools
using Dates
using DelimitedFiles
using LinearAlgebra
using Printf
using Random
using Statistics

import BSPivotQR
using BSRandPivotQR

if Sys.isapple()
    include(joinpath(@__DIR__, "..", "..", "julia", "benchmark", "setup_accelerate.jl"))
end

include(joinpath(@__DIR__, "rand_test_matrix.jl"))

function _parse_sizes(s::AbstractString)
    return [
        (parse(Int, split(t, 'x')[1]), parse(Int, split(t, 'x')[2]))
        for t in split(s, ',') if !isempty(strip(t))
    ]
end

const SIZES = _parse_sizes(
    get(ENV, "BS_RAND_SIZES", "32x2000,64x4000,64x8000,128x8000,128x16000"))
const FAMILIES = Symbol.(split(get(ENV, "BS_RAND_FAMILIES", "gaussian,needle"), ','))
const SEED = parse(Int, get(ENV, "BS_RAND_SEED", "1"))
const OUTDIR = get(ENV, "BS_RAND_OUTDIR", joinpath(@__DIR__, "results"))

_frobinv(R::AbstractMatrix) = sqrt(sum(abs2, 1.0 ./ svdvals(R)))

_median_s(trial) = median(trial).time / 1e9

function main()
    mkpath(OUTDIR)
    println("BLAS threads = ", BLAS.get_num_threads(),
        "  julia threads = ", Threads.nthreads())
    @printf("%-10s %-12s %-18s %10s %10s %10s %8s %8s %7s %8s\n",
        "family", "size", "mode", "t_rand(ms)", "t_det(ms)", "t_qr(ms)",
        "spd_det", "spd_qr", "cond", "tested/k")

    header = ["family" "k" "n" "mode" "t_rand_s" "t_det_s" "t_builtin_s" #=
        =# "speedup_det" "speedup_builtin" "frobinv_rand" "frobinv_det" #=
        =# "frobinv_builtin" "osinsky" "tested_per_k"]
    rows = Vector{Any}[]

    for fam in FAMILIES
        for (ci, (k, n)) in enumerate(SIZES)
            M = rand_test_matrix(fam, k, n, SEED + ci)
            s = SEED + ci
            bs = k   # batched default block

            # Deterministic BSQR: Q, R, perm materialized.
            t_det = _median_s(@benchmark begin
                F = BSPivotQR.bsqr($M; k = $k, check = false)
                (BSPivotQR.Q(F), BSPivotQR.R(F), BSPivotQR.perm(F))
            end)
            Fdet = BSPivotQR.bsqr(M; k = k, check = false)
            frobinv_det = _frobinv(BSPivotQR.R(Fdet)[1:k, 1:k])

            # Built-in column-pivoted QR (LAPACK dgeqp3): Q, R, p materialized.
            t_qr = _median_s(@benchmark begin
                F = qr($M, ColumnNorm())
                (Matrix(F.Q), F.R, F.p)
            end)
            Fqr = qr(M, ColumnNorm())
            frobinv_qr = _frobinv(Fqr.R[1:k, 1:k])

            osinsky = sqrt(k * (n - k + 1))

            # Randomized: [p, Q, R11], defaults (batched, norm-weighted,
            # running_mean), block = k.
            t_rand = _median_s(@benchmark begin
                F = bsqr_rand($M; k = $k, block_size = $bs, seed = $s)
                (F.p, BSRandPivotQR.Q(F), F.R11)
            end)
            Frand = bsqr_rand(M; k = k, block_size = bs, seed = s)
            push!(rows, _report_row(fam, k, n, "running_mean", t_rand, t_det, t_qr,
                Frand.stats.frob_inv, frobinv_det, frobinv_qr, osinsky,
                Frand.stats.total_tested / k))

            # R12 desired: selection plus the final Q'-apply to the leftovers.
            t_r12 = _median_s(@benchmark begin
                F = bsqr_rand($M; k = $k, block_size = $bs, seed = $s,
                    return_r12 = true)
                (F.p, BSRandPivotQR.Q(F), F.R11, F.R12)
            end)
            F12 = bsqr_rand(M; k = k, block_size = bs, seed = s, return_r12 = true)
            push!(rows, _report_row(fam, k, n, "running_mean+R12", t_r12, t_det, t_qr,
                F12.stats.frob_inv, frobinv_det, frobinv_qr, osinsky,
                F12.stats.total_tested / k))
        end
    end

    csv = joinpath(OUTDIR, "rand_timings.csv")
    open(csv, "w") do io
        writedlm(io, header, ',')
        writedlm(io, permutedims(reduce(hcat, [Any[r...] for r in rows])), ',')
    end
    open(joinpath(OUTDIR, "metadata.txt"), "w") do io
        println(io, "generated: ", Dates.now())
        println(io, "julia: ", VERSION)
        println(io, "blas: ", sprint(show, BLAS.get_config()))
        println(io, "blas threads: ", BLAS.get_num_threads())
        println(io, "sizes: ", SIZES)
        println(io, "families: ", FAMILIES)
        println(io, "seed: ", SEED)
    end
    println("\nWrote ", csv)
end

function _report_row(fam, k, n, mode, t_rand, t_det, t_qr,
                     fr_rand, fr_det, fr_qr, osinsky, tested)
    spd_det = t_det / t_rand
    spd_qr = t_qr / t_rand
    cond = fr_rand / fr_det
    @printf("%-10s %-12s %-18s %10.3f %10.3f %10.3f %8.2f %8.2f %7.2f %8.1f\n",
        fam, "$(k)x$(n)", mode, t_rand * 1e3, t_det * 1e3, t_qr * 1e3,
        spd_det, spd_qr, cond, tested)
    return Any[fam, k, n, mode, t_rand, t_det, t_qr, spd_det, spd_qr,
        fr_rand, fr_det, fr_qr, osinsky, tested]
end

main()
