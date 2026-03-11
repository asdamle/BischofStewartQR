#!/usr/bin/env julia

using Dates
using LinearAlgebra
using Printf
using Random
using Statistics

using BSPivotQR

include("bench_common.jl")
using .BenchCommon

const DEFAULT_PUB_OUTDIR = joinpath(@__DIR__, "results", "publication")
const ALLOWED_METHODS = Set([BSQR_METHOD_LABEL, DGEQP3_METHOD_LABEL])

_parse_families() = parse_symbol_list(
    "BS_PUB_FAMILIES",
    [:gaussian, :ill_conditioned, :orthonormal_rows],
    [:gaussian, :ill_conditioned, :orthonormal_rows],
)

function _parse_threads()
    vals = parse_int_list("BS_PUB_THREADS", [1, 4])
    isempty(vals) && error("BS_PUB_THREADS must contain at least one value")
    all(>=(1), vals) || error("BS_PUB_THREADS entries must be >= 1")
    return vals
end

function _parse_seeds()
    vals = parse_int_list("BS_PUB_SEEDS", [20260310, 20260311])
    isempty(vals) && error("BS_PUB_SEEDS must contain at least one value")
    return vals
end

function _parse_square_ms()
    vals = parse_int_list("BS_PUB_SQUARE_MS", [64, 128, 256, 384, 512])
    isempty(vals) && error("BS_PUB_SQUARE_MS must contain at least one value")
    all(>=(1), vals) || error("BS_PUB_SQUARE_MS entries must be >= 1")
    return vals
end

function _parse_short_ms()
    vals = parse_int_list("BS_PUB_SHORT_MS", [32, 64, 128, 256, 512])
    isempty(vals) && error("BS_PUB_SHORT_MS must contain at least one value")
    all(>=(1), vals) || error("BS_PUB_SHORT_MS entries must be >= 1")
    return vals
end

function _parse_short_aspects()
    vals = parse_float_list("BS_PUB_SHORT_ASPECTS", [2.0, 4.0, 8.0, 10.0])
    isempty(vals) && error("BS_PUB_SHORT_ASPECTS must contain at least one value")
    all(a -> a >= 1.0, vals) || error("BS_PUB_SHORT_ASPECTS entries must be >= 1.0")
    return vals
end

function _case_grid(square_ms::Vector{Int}, short_ms::Vector{Int}, short_aspects::Vector{Float64})
    cases = NamedTuple[]
    for m in square_ms
        push!(cases, (regime = "square", m = m, n = m, aspect = 1.0))
    end
    for m in short_ms
        for aspect in short_aspects
            n = round(Int, m * aspect)
            n >= m || error("Short-wide case generated n < m (m=$m, aspect=$aspect)")
            push!(cases, (regime = "short_wide", m = m, n = n, aspect = aspect))
        end
    end
    return cases
end

function _write_metadata(
    metadata_path::String,
    run_id::String,
    threads::Vector{Int},
    seeds::Vector{Int},
    families::Vector{Symbol},
    square_ms::Vector{Int},
    short_ms::Vector{Int},
    short_aspects::Vector{Float64},
    warmup::Int,
    samples::Int,
)
    cpu_model = try
        cpu = first(Sys.cpu_info())
        getproperty(cpu, :model)
    catch
        "unknown"
    end

    open(metadata_path, "w") do io
        println(io, "run_id = ", run_id)
        println(io, "timestamp = ", Dates.now())
        println(io, "julia_version = ", VERSION)
        println(io, "cpu_model = ", cpu_model)
        println(io, "cpu_threads = ", Sys.CPU_THREADS)
        println(io, "blas_threads_at_end = ", BLAS.get_num_threads())
        println(io, "blas_config = ", backend_string())
        println(io, "")
        println(io, "[publication_knobs]")
        println(io, "BS_PUB_THREADS = ", join(threads, ","))
        println(io, "BS_PUB_SEEDS = ", join(seeds, ","))
        println(io, "BS_PUB_FAMILIES = ", join(string.(families), ","))
        println(io, "BS_PUB_SQUARE_MS = ", join(square_ms, ","))
        println(io, "BS_PUB_SHORT_MS = ", join(short_ms, ","))
        println(io, "BS_PUB_SHORT_ASPECTS = ", join(string.(short_aspects), ","))
        println(io, "BS_PUB_WARMUP = ", warmup)
        println(io, "BS_PUB_SAMPLES = ", samples)
        println(io, "BS_NORM_RECOMP_TOL = ", get(ENV, "BS_NORM_RECOMP_TOL", string(DEFAULT_NORM_RECOMP_TOL)))
        println(io, "")
        println(io, "[environment]")
        for key in sort!(collect(keys(ENV)))
            startswith(key, "BS_") || continue
            println(io, key, " = ", ENV[key])
        end
    end
end

function _speedup_rows(rows::Vector{NamedTuple})
    idx = Dict{Tuple{Symbol,String,Int,Int,Int,Int,String},Float64}()
    for r in rows
        key = (r.family, r.regime, r.m, r.n, r.seed, r.blas_threads, r.method)
        idx[key] = r.tmed
    end

    speed = NamedTuple[]
    for r in rows
        r.method == BSQR_METHOD_LABEL || continue
        key_dg = (r.family, r.regime, r.m, r.n, r.seed, r.blas_threads, DGEQP3_METHOD_LABEL)
        haskey(idx, key_dg) || continue
        dg_t = idx[key_dg]
        sp = r.tmed == 0.0 ? NaN : (dg_t / r.tmed)
        push!(speed, (
            family = r.family,
            regime = r.regime,
            m = r.m,
            n = r.n,
            aspect = r.aspect,
            seed = r.seed,
            blas_threads = r.blas_threads,
            speedup = sp,
            bsqr_tmed_s = r.tmed,
            dgeqp3_tmed_s = dg_t,
        ))
    end
    return speed
end

function _geomean(vals::Vector{Float64})
    good = filter(v -> isfinite(v) && v > 0.0, vals)
    isempty(good) && return NaN
    return exp(mean(log.(good)))
end

function _write_summary(md_path::String, rows::Vector{NamedTuple}, run_id::String, threads::Vector{Int}, seeds::Vector{Int})
    speeds = _speedup_rows(rows)

    grouped = Dict{Tuple{Symbol,String,Int},Vector{Float64}}()
    for s in speeds
        key = (s.family, s.regime, s.blas_threads)
        get!(grouped, key, Float64[])
        push!(grouped[key], s.speedup)
    end

    open(md_path, "w") do md
        println(md, "# Publication Benchmark Summary")
        println(md, "")
        println(md, "- Run ID: `$run_id`")
        println(md, "- Generated: $(Dates.now())")
        println(md, "- BLAS: ", backend_string())
        println(md, "- Threads: ", join(threads, ", "))
        println(md, "- Seeds: ", join(seeds, ", "))
        println(md, "")
        println(md, "| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr) |")
        println(md, "|---|---|---:|---:|")
        for key in sort!(collect(keys(grouped)), by = x -> (string(x[1]), x[2], x[3]))
            g = _geomean(grouped[key])
            println(md, "| $(key[1]) | $(key[2]) | $(key[3]) | $(round(g, sigdigits=5)) |")
        end
    end
end

function _validate_pairs(rows::Vector{NamedTuple})
    keyset = Dict{Tuple{Symbol,String,Int,Int,Int,Int},Set{String}}()
    for r in rows
        key = (r.family, r.regime, r.m, r.n, r.seed, r.blas_threads)
        get!(keyset, key, Set{String}())
        push!(keyset[key], r.method)
    end
    for (k, methods) in keyset
        methods == ALLOWED_METHODS || error("Missing benchmark pair for key=$k; methods present=$(collect(methods))")
    end
end

function run_publication_benchmarks()
    check_backend()
    norm_recomp_tol = parse_env_float("BS_NORM_RECOMP_TOL", DEFAULT_NORM_RECOMP_TOL)

    outdir = strip(get(ENV, "BS_PUB_OUTDIR", DEFAULT_PUB_OUTDIR))
    isempty(outdir) && error("BS_PUB_OUTDIR cannot be empty")
    mkpath(outdir)

    threads = _parse_threads()
    seeds = _parse_seeds()
    families = _parse_families()
    square_ms = _parse_square_ms()
    short_ms = _parse_short_ms()
    short_aspects = _parse_short_aspects()
    warmup = parse_env_int("BS_PUB_WARMUP", 1; minval = 0)
    samples = parse_env_int("BS_PUB_SAMPLES", 30; minval = 1)

    run_id = Dates.format(now(), "yyyymmdd_HHMMSS")
    timestamp = string(Dates.now())
    csv_path = joinpath(outdir, "publication_timings.csv")
    summary_path = joinpath(outdir, "publication_summary.md")
    metadata_path = joinpath(outdir, "metadata.txt")

    io = open(csv_path, "w")
    println(
        io,
        "run_id,timestamp,family,regime,m,n,aspect,seed,blas_threads,method,tmin_s,tmed_s,tci_low_s,tci_high_s,alloc_bytes,residual,orthogonality",
    )

    rows = NamedTuple[]
    cases = _case_grid(square_ms, short_ms, short_aspects)
    old_threads = BLAS.get_num_threads()

    try
        for nth in threads
            BLAS.set_num_threads(nth)
            println("Configured BLAS threads: ", BLAS.get_num_threads())

            for seed in seeds
                rng = MersenneTwister(seed)
                for family in families
                    for c in cases
                        if family === :orthonormal_rows && c.m > c.n
                            continue
                        end

                        A = make_matrix(family, c.m, c.n, rng)
                        kfull = min(c.m, c.n)
                        rowset = bench_pair_ci(A, kfull, norm_recomp_tol; warmup = warmup, samples = samples)

                        for row in rowset
                            row.method in ALLOWED_METHODS || continue
                            outrow = (
                                run_id = run_id,
                                timestamp = timestamp,
                                family = family,
                                regime = c.regime,
                                m = c.m,
                                n = c.n,
                                aspect = c.aspect,
                                seed = seed,
                                blas_threads = nth,
                                method = row.method,
                                tmin = row.tmin,
                                tmed = row.tmed,
                                tci_low = row.tci_low,
                                tci_high = row.tci_high,
                                alloc = row.alloc,
                                resid = row.resid,
                                orth = row.orth,
                            )
                            push!(rows, outrow)
                            println(
                                io,
                                "$(outrow.run_id),$(outrow.timestamp),$(outrow.family),$(outrow.regime),$(outrow.m),$(outrow.n),$(@sprintf("%.6f", outrow.aspect)),$(outrow.seed),$(outrow.blas_threads),$(outrow.method),$(outrow.tmin),$(outrow.tmed),$(outrow.tci_low),$(outrow.tci_high),$(outrow.alloc),$(outrow.resid),$(outrow.orth)",
                            )
                        end
                    end
                end
            end
        end
    finally
        close(io)
        BLAS.set_num_threads(old_threads)
    end

    _validate_pairs(rows)
    _write_summary(summary_path, rows, run_id, threads, seeds)
    _write_metadata(
        metadata_path,
        run_id,
        threads,
        seeds,
        families,
        square_ms,
        short_ms,
        short_aspects,
        warmup,
        samples,
    )

    println("Wrote publication benchmark CSV to: $csv_path")
    println("Wrote publication summary to: $summary_path")
    println("Wrote publication metadata to: $metadata_path")
end

run_publication_benchmarks()
