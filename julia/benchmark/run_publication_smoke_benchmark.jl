#!/usr/bin/env julia

function _save_env(keys)
    old = Dict{String,Union{Nothing,String}}()
    for k in keys
        old[k] = haskey(ENV, k) ? ENV[k] : nothing
    end
    return old
end

function _restore_env!(old)
    for (k, v) in old
        if v === nothing
            haskey(ENV, k) && delete!(ENV, k)
        else
            ENV[k] = v
        end
    end
    return nothing
end

env_keys = [
    "BS_PUB_OUTDIR",
    "BS_PUB_THREADS",
    "BS_PUB_SEEDS",
    "BS_PUB_FAMILIES",
    "BS_PUB_SQUARE_MS",
    "BS_PUB_SHORT_MS",
    "BS_PUB_SHORT_ASPECTS",
    "BS_PUB_WARMUP",
    "BS_PUB_SAMPLES",
]

old = _save_env(env_keys)
try
    ENV["BS_PUB_OUTDIR"] = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "results", "publication_smoke")
    ENV["BS_PUB_THREADS"] = "1"
    ENV["BS_PUB_SEEDS"] = "20260310"
    ENV["BS_PUB_FAMILIES"] = "gaussian"
    ENV["BS_PUB_SQUARE_MS"] = "64"
    ENV["BS_PUB_SHORT_MS"] = "32"
    ENV["BS_PUB_SHORT_ASPECTS"] = "2.0"
    ENV["BS_PUB_WARMUP"] = "1"
    ENV["BS_PUB_SAMPLES"] = "3"

    include(joinpath(@__DIR__, "bench_cpqr_publication.jl"))
finally
    _restore_env!(old)
end
