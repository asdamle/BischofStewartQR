#!/usr/bin/env julia

include("validation_pipeline_common.jl")
using .ValidationPipelineCommon

function run_tier1(; baseline_dir::Union{Nothing,String} = nothing)
    envvars = build_validation_env(; quick = true, require_accelerate_default = "0")
    run_validation_pipeline(envvars; baseline_dir = baseline_dir)
end

function main()
    baseline_dir = resolve_baseline_dir(ARGS)
    run_tier1(; baseline_dir = baseline_dir)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
