#!/usr/bin/env julia

using Dates
using LinearAlgebra

include("bench_common.jl")
using .BenchCommon

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULTS_DIR = joinpath(@__DIR__, "results")

const SNAPSHOT_FILES = [
    "timings.csv",
    "summary.md",
    "sweep_timings.csv",
    "sweep_summary.md",
    "regime_timings.csv",
    "regime_summary.md",
    "profile_breakdown.csv",
    "profile_breakdown_summary.md",
]

const SNAPSHOT_DIRS = ["plots", "sweep_plots", "regime_plots"]

function _run_step(cmd::Cmd, envvars::Dict{String,String})
    println(">> ", cmd)
    run(addenv(Cmd(cmd; dir = REPO_ROOT), envvars...))
end

function _write_metadata(path::String, envvars::Dict{String,String})
    cpu_model = try
        cpu = first(Sys.cpu_info())
        getproperty(cpu, :model)
    catch
        "unknown"
    end

    open(path, "w") do io
        println(io, "timestamp = ", Dates.now())
        println(io, "julia_version = ", VERSION)
        println(io, "cpu_model = ", cpu_model)
        println(io, "cpu_threads = ", Sys.CPU_THREADS)
        println(io, "blas_threads = ", BLAS.get_num_threads())
        println(io, "blas_config = ", backend_string())
        println(io, "")
        println(io, "[env]")
        for k in sort!(collect(keys(envvars)))
            println(io, "$k = ", envvars[k])
        end
    end
end

function _copy_outputs!(dest::String)
    mkpath(dest)
    for f in SNAPSHOT_FILES
        src = joinpath(RESULTS_DIR, f)
        isfile(src) || error("Expected result file missing: $src")
        cp(src, joinpath(dest, f); force = true)
    end
    for d in SNAPSHOT_DIRS
        src = joinpath(RESULTS_DIR, d)
        isdir(src) || continue
        cp(src, joinpath(dest, d); force = true)
    end
end

function run_freeze()
    check_backend()
    configure_blas_threads()

    tag = strip(get(ENV, "BS_BASELINE_TAG", ""))
    if isempty(tag)
        tag = Dates.format(now(), "yyyymmdd_HHMMSS")
    end
    outdir = joinpath(RESULTS_DIR, "baseline_" * tag)
    mkpath(outdir)

    envvars = Dict(
        "BS_BLAS_THREADS" => get(ENV, "BS_BLAS_THREADS", string(BLAS.get_num_threads())),
        "BS_USE_ACCELERATE" => get(ENV, "BS_USE_ACCELERATE", "1"),
        "BS_REQUIRE_ACCELERATE" => get(ENV, "BS_REQUIRE_ACCELERATE", Sys.isapple() ? "1" : "0"),
        "BS_NORM_RECOMP_TOL" => get(ENV, "BS_NORM_RECOMP_TOL", string(DEFAULT_NORM_RECOMP_TOL)),
        "BS_QUICK" => "1",
        "BS_SWEEP_QUICK" => "1",
        "BS_REGIME_QUICK" => "1",
        "BS_PROFILE_QUICK" => "1",
    )

    _run_step(`julia --project=. test/runtests.jl`, envvars)
    _run_step(`julia --project=. benchmark/bench_cpqr.jl`, envvars)
    _run_step(`julia --project=. benchmark/bench_cpqr_sweep.jl`, envvars)
    _run_step(`julia --project=. benchmark/bench_cpqr_regimes.jl`, envvars)
    _run_step(`julia --project=. benchmark/profile_bsqr_breakdown.jl`, envvars)
    _run_step(`julia --project=. benchmark/plot_results.jl`, envvars)
    _run_step(`julia --project=. benchmark/plot_sweep_results.jl`, envvars)
    _run_step(`julia --project=. benchmark/plot_regime_results.jl`, envvars)

    _copy_outputs!(outdir)
    _write_metadata(joinpath(outdir, "metadata.txt"), envvars)
    println("Baseline snapshot written to: $outdir")
end

if abspath(PROGRAM_FILE) == @__FILE__
    run_freeze()
end
