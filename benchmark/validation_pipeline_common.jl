module ValidationPipelineCommon

using LinearAlgebra

export REPO_ROOT, build_validation_env, run_validation_pipeline, resolve_baseline_dir

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const VALIDATION_STEPS = (
    `julia --project=. test/runtests.jl`,
    `julia --project=. benchmark/bench_cpqr.jl`,
    `julia --project=. benchmark/bench_cpqr_sweep.jl`,
    `julia --project=. benchmark/bench_cpqr_regimes.jl`,
    `julia --project=. benchmark/profile_bsqr_breakdown.jl`,
    `julia --project=. benchmark/plot_results.jl`,
    `julia --project=. benchmark/plot_sweep_results.jl`,
    `julia --project=. benchmark/plot_regime_results.jl`,
)

function _run_step(cmd::Cmd, envvars::Dict{String,String})
    println(">> ", cmd)
    run(addenv(Cmd(cmd; dir = REPO_ROOT), envvars...))
end

function build_validation_env(; quick::Bool, require_accelerate_default::String)
    envvars = Dict(
        "BS_BLAS_THREADS" => get(ENV, "BS_BLAS_THREADS", string(BLAS.get_num_threads())),
        "BS_USE_ACCELERATE" => get(ENV, "BS_USE_ACCELERATE", "1"),
        "BS_REQUIRE_ACCELERATE" => get(ENV, "BS_REQUIRE_ACCELERATE", require_accelerate_default),
        "BS_NORM_RECOMP_TOL" => get(ENV, "BS_NORM_RECOMP_TOL", string(sqrt(eps(Float64)))),
    )
    if quick
        envvars["BS_QUICK"] = "1"
        envvars["BS_SWEEP_QUICK"] = "1"
        envvars["BS_REGIME_QUICK"] = "1"
        envvars["BS_PROFILE_QUICK"] = "1"
    end
    return envvars
end

function run_validation_pipeline(
    envvars::Dict{String,String};
    baseline_dir::Union{Nothing,String} = nothing,
)
    for cmd in VALIDATION_STEPS
        _run_step(cmd, envvars)
    end
    if baseline_dir !== nothing && !isempty(strip(baseline_dir))
        _run_step(`julia --project=. benchmark/compare_results.jl $baseline_dir`, envvars)
    end
    return nothing
end

function resolve_baseline_dir(args::Vector{String})
    return length(args) >= 1 ? args[1] : get(ENV, "BS_BASELINE_DIR", "")
end

end
