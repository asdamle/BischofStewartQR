#!/usr/bin/env julia

using LinearAlgebra

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))

function _run_step(cmd::Cmd, envvars::Dict{String,String})
    println(">> ", cmd)
    run(addenv(Cmd(cmd; dir = REPO_ROOT), envvars...))
end

function run_tier1(; baseline_dir::Union{Nothing,String} = nothing)
    envvars = Dict(
        "BS_BLAS_THREADS" => get(ENV, "BS_BLAS_THREADS", string(BLAS.get_num_threads())),
        "BS_USE_ACCELERATE" => get(ENV, "BS_USE_ACCELERATE", "1"),
        "BS_REQUIRE_ACCELERATE" => get(ENV, "BS_REQUIRE_ACCELERATE", "0"),
        "BS_NORM_RECOMP_TOL" => get(ENV, "BS_NORM_RECOMP_TOL", string(sqrt(eps(Float64)))),
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

    if baseline_dir !== nothing && !isempty(strip(baseline_dir))
        _run_step(`julia --project=. benchmark/compare_results.jl $baseline_dir`, envvars)
    end
end

function main()
    baseline_dir = if length(ARGS) >= 1
        ARGS[1]
    else
        get(ENV, "BS_BASELINE_DIR", "")
    end
    run_tier1(; baseline_dir = baseline_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
