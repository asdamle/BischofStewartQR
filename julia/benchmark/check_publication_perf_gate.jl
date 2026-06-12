#!/usr/bin/env julia

using DelimitedFiles
using Printf
using Statistics

const DEFAULT_INPUT = joinpath(@__DIR__, "results", "publication", "publication_timings.csv")
const METHODS = Set(["bsqr_full", "bsqr_rinv"])

_as_string(x) = String(x)
_as_int(x) = x isa Int ? x : parse(Int, _as_string(x))
_as_float(x) = x isa AbstractFloat ? Float64(x) : parse(Float64, _as_string(x))

function _load(csv_path::String)
    isfile(csv_path) || error("CSV not found: $csv_path")
    data, raw_header = readdlm(csv_path, ',', header = true)
    header = _as_string.(vec(raw_header))
    idx = Dict{String,Int}(h => i for (i, h) in enumerate(header))
    required = ["family", "regime", "m", "n", "aspect", "seed", "method", "tmed_s"]
    for c in required
        haskey(idx, c) || error("Missing required column '$c' in $csv_path")
    end

    out = Dict{Tuple{String,String,Int,Int,Float64,Int,String},Float64}()
    for i in 1:size(data, 1)
        row = view(data, i, :)
        method = _as_string(row[idx["method"]])
        method in METHODS || continue
        key = (
            _as_string(row[idx["family"]]),
            _as_string(row[idx["regime"]]),
            _as_int(row[idx["m"]]),
            _as_int(row[idx["n"]]),
            _as_float(row[idx["aspect"]]),
            _as_int(row[idx["seed"]]),
            method,
        )
        out[key] = _as_float(row[idx["tmed_s"]])
    end
    return out
end

function main()
    baseline_csv = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_INPUT
    candidate_csv = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_INPUT
    max_slowdown = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.05
    # Cells whose baseline median is below this floor cannot carry a slowdown
    # signal across runs (timer/scheduler noise dominates sub-ms cases); they
    # are reported in the aggregates but exempt from the gate.
    min_tmed = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 0.0
    max_slowdown >= 0 || error("max_slowdown must be >= 0")
    min_tmed >= 0 || error("min_tmed must be >= 0")

    base = _load(baseline_csv)
    cand = _load(candidate_csv)

    common_keys = sort!(collect(intersect(Set(Base.keys(base)), Set(Base.keys(cand)))))
    isempty(common_keys) && error("No overlapping rows between baseline and candidate.")

    by_method = Dict("bsqr_full" => Float64[], "bsqr_rinv" => Float64[])
    violations = Tuple{Tuple{String,String,Int,Int,Float64,Int,String},Float64}[]
    ngated = 0
    for k in common_keys
        tb = base[k]
        tc = cand[k]
        tb > 0 || error("Non-positive baseline tmed for key=$k")
        slowdown = tc / tb - 1.0
        push!(by_method[k[7]], slowdown)
        tb >= min_tmed || continue
        ngated += 1
        if isfinite(slowdown) && slowdown > max_slowdown
            push!(violations, (k, slowdown))
        end
    end

    println("Performance gate compared $(length(common_keys)) rows, $ngated gated (threshold $(round(100 * max_slowdown, digits=2))% slowdown, min_tmed $(min_tmed)s).")
    println("| method | median_slowdown | max_slowdown |")
    println("|---|---:|---:|")
    for m in ["bsqr_full", "bsqr_rinv"]
        vals = by_method[m]
        isempty(vals) && continue
        @printf("| %s | %.6f | %.6f |\n", m, median(vals), maximum(vals))
    end

    if !isempty(violations)
        println("\nGate failed with $(length(violations)) violating rows.")
        for (k, s) in violations
            @printf("%s slowdown=%.6f\n", repr(k), s)
        end
        error("Performance gate failed: slowdown exceeded threshold.")
    end

    println("Performance gate passed.")
end

main()
