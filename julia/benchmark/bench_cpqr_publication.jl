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
const PUBLICATION_SCHEMA_VERSION = "2026-03-12.v1"
const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const PLAIN_METHODS = Set([BSQR_METHOD_LABEL, DGEQP3_METHOD_LABEL])
const RINV_METHODS = Set([BSQR_RINV_METHOD_LABEL, DGEQP3_TRSM_METHOD_LABEL])
const SQUARE_METHODS = Set(vcat(collect(PLAIN_METHODS), collect(RINV_METHODS)))
const SHORT_WIDE_METHODS = Set(vcat(collect(PLAIN_METHODS), collect(RINV_METHODS)))
const ALLOWED_METHODS = Set(vcat(collect(SQUARE_METHODS), collect(SHORT_WIDE_METHODS)))

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

@inline function _family_code(family::Symbol)
    family === :gaussian && return UInt(1)
    family === :ill_conditioned && return UInt(2)
    family === :orthonormal_rows && return UInt(3)
    error("Unknown family: $family")
end

@inline function _case_seed(seed::Int, family::Symbol, m::Int, n::Int, aspect::Float64)
    # Deterministic per-case seed that does not depend on loop traversal order.
    aspect_key = UInt(round(Int, aspect * 1_000_000))
    x = reinterpret(UInt, seed)
    x ⊻= (_family_code(family) << 48)
    x ⊻= (UInt(m) << 24)
    x ⊻= UInt(n)
    x ⊻= (aspect_key << 8)
    return Int(mod(x, UInt(typemax(Int) - 1))) + 1
end

function _git_string(args...)
    try
        cmd = Cmd(["git", "-C", REPO_ROOT, String.(args)...])
        return strip(read(cmd, String))
    catch
        return "unknown"
    end
end

function _git_dirty_state()
    try
        wd_clean = success(`git -C $REPO_ROOT diff --quiet --ignore-submodules=all`)
        ix_clean = success(`git -C $REPO_ROOT diff --cached --quiet --ignore-submodules=all`)
        return !(wd_clean && ix_clean)
    catch
        return "unknown"
    end
end

function _expected_row_count(
    threads::Vector{Int},
    seeds::Vector{Int},
    families::Vector{Symbol},
    cases::Vector{<:NamedTuple},
)
    methods_per_case = length(ALLOWED_METHODS)
    return length(threads) * length(seeds) * length(families) * length(cases) * methods_per_case
end

function _validate_ci_and_quality(rows::Vector{NamedTuple})
    ci_enforce = strip(get(ENV, "BS_PUB_CI_ENFORCE", "0")) == "1"
    ci_warn_frac = parse_env_float("BS_PUB_CI_WARN_FRAC", 0.5; minval = 0.0, maxval = 100.0)
    ci_fail_frac = parse_env_float("BS_PUB_CI_FAIL_FRAC", 10.0; minval = 0.0, maxval = 1000.0)
    ci_min_tmed = parse_env_float("BS_PUB_CI_MIN_TMED", 1.0e-4; minval = 0.0, maxval = 1.0e6)
    resid_factor = parse_env_float("BS_PUB_RESID_FACTOR", 2.5e3; minval = 1.0, maxval = 1.0e12)
    orth_factor = parse_env_float("BS_PUB_ORTH_FACTOR", 2.5e3; minval = 1.0, maxval = 1.0e12)

    warn_keys = Tuple{Symbol,String,Int,Int,Int,Int,String,Float64}[]
    for r in rows
        tmin = r.tmin
        tmed = r.tmed
        tci_low = r.tci_low
        tci_high = r.tci_high
        alloc = r.alloc
        all(isfinite, (tmin, tmed, tci_low, tci_high, r.resid, r.orth)) ||
            error("Non-finite benchmark row encountered: $r")
        (tmin >= 0.0 && tmed > 0.0 && tci_low >= 0.0 && tci_high >= 0.0) ||
            error("Negative timing encountered: $r")
        (tci_low <= tmed <= tci_high) ||
            error("Invalid CI bounds (must satisfy tci_low <= tmed <= tci_high): $r")
        alloc >= 0 || error("Negative allocation count encountered: $r")

        ci_frac = (tci_high - tci_low) / max(tmed, eps(Float64))
        if tmed >= ci_min_tmed && ci_frac > ci_fail_frac
            if ci_enforce
                error("CI spread too wide for publication stability (frac=$(round(ci_frac, sigdigits=4)) > $ci_fail_frac): $r")
            else
                push!(warn_keys, (r.family, r.regime, r.m, r.n, r.seed, r.blas_threads, r.method, ci_frac))
            end
        end
        if tmed >= ci_min_tmed && ci_frac > ci_warn_frac
            push!(warn_keys, (r.family, r.regime, r.m, r.n, r.seed, r.blas_threads, r.method, ci_frac))
        end

        tol_resid = resid_factor * eps(Float64) * max(r.m, r.n)
        tol_orth = orth_factor * eps(Float64) * max(r.m, 1)
        r.resid <= tol_resid || error("Residual exceeds tolerance ($tol_resid): $r")
        r.orth <= tol_orth || error("Orthogonality exceeds tolerance ($tol_orth): $r")
    end

    if !isempty(warn_keys)
        uniq = unique(warn_keys)
        worst = maximum(x -> x[8], uniq)
        println(
            "Warning: high CI spread in $(length(uniq)) benchmark groups " *
            "(warn>$(ci_warn_frac), fail>$(ci_fail_frac), enforce=$(ci_enforce), worst=$(round(worst, sigdigits=4))).",
        )
    end
    return nothing
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
    expected_rows::Int,
    observed_rows::Int,
)
    cpu_model = try
        cpu = first(Sys.cpu_info())
        getproperty(cpu, :model)
    catch
        "unknown"
    end
    git_sha = _git_string("rev-parse", "HEAD")
    git_branch = _git_string("rev-parse", "--abbrev-ref", "HEAD")
    git_dirty = _git_dirty_state()

    open(metadata_path, "w") do io
        println(io, "schema_version = ", PUBLICATION_SCHEMA_VERSION)
        println(io, "run_id = ", run_id)
        println(io, "timestamp = ", Dates.now())
        println(io, "julia_version = ", VERSION)
        println(io, "git_sha = ", git_sha)
        println(io, "git_branch = ", git_branch)
        println(io, "git_dirty = ", git_dirty)
        println(io, "repo_root = ", REPO_ROOT)
        println(io, "cpu_model = ", cpu_model)
        println(io, "cpu_threads = ", Sys.CPU_THREADS)
        println(io, "blas_threads_at_end = ", BLAS.get_num_threads())
        println(io, "blas_config = ", backend_string())
        println(io, "expected_rows = ", expected_rows)
        println(io, "observed_rows = ", observed_rows)
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
        println(io, "BS_PUB_CI_WARN_FRAC = ", get(ENV, "BS_PUB_CI_WARN_FRAC", "0.5"))
        println(io, "BS_PUB_CI_FAIL_FRAC = ", get(ENV, "BS_PUB_CI_FAIL_FRAC", "10.0"))
        println(io, "BS_PUB_CI_MIN_TMED = ", get(ENV, "BS_PUB_CI_MIN_TMED", "1.0e-4"))
        println(io, "BS_PUB_CI_ENFORCE = ", get(ENV, "BS_PUB_CI_ENFORCE", "0"))
        println(io, "BS_PUB_RESID_FACTOR = ", get(ENV, "BS_PUB_RESID_FACTOR", "2500.0"))
        println(io, "BS_PUB_ORTH_FACTOR = ", get(ENV, "BS_PUB_ORTH_FACTOR", "2500.0"))
        println(io, "")
        println(io, "[environment]")
        for key in sort!(collect(keys(ENV)))
            startswith(key, "BS_") || continue
            println(io, key, " = ", ENV[key])
        end
    end
end

function _speedup_rows(
    rows::Vector{NamedTuple},
    bsqr_method::String,
    baseline_method::String;
    regime::Union{Nothing,String} = nothing,
)
    idx = Dict{Tuple{Symbol,String,Int,Int,Int,Int,String},Float64}()
    for r in rows
        key = (r.family, r.regime, r.m, r.n, r.seed, r.blas_threads, r.method)
        idx[key] = r.tmed
    end

    speed = NamedTuple[]
    for r in rows
        r.method == bsqr_method || continue
        regime !== nothing && r.regime != regime && continue
        key_baseline = (r.family, r.regime, r.m, r.n, r.seed, r.blas_threads, baseline_method)
        haskey(idx, key_baseline) || continue
        baseline_t = idx[key_baseline]
        sp = r.tmed == 0.0 ? NaN : (baseline_t / r.tmed)
        push!(speed, (
            family = r.family,
            regime = r.regime,
            m = r.m,
            n = r.n,
            aspect = r.aspect,
            seed = r.seed,
            blas_threads = r.blas_threads,
            speedup = sp,
            bsqr_method = bsqr_method,
            baseline_method = baseline_method,
            bsqr_tmed_s = r.tmed,
            baseline_tmed_s = baseline_t,
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
    plain_speeds = _speedup_rows(rows, BSQR_METHOD_LABEL, DGEQP3_METHOD_LABEL)
    rinv_speeds = _speedup_rows(rows, BSQR_RINV_METHOD_LABEL, DGEQP3_TRSM_METHOD_LABEL)

    grouped_plain = Dict{Tuple{Symbol,String,Int},Vector{Float64}}()
    for s in plain_speeds
        key = (s.family, s.regime, s.blas_threads)
        get!(grouped_plain, key, Float64[])
        push!(grouped_plain[key], s.speedup)
    end
    grouped_rinv = Dict{Tuple{Symbol,String,Int},Vector{Float64}}()
    for s in rinv_speeds
        key = (s.family, s.regime, s.blas_threads)
        get!(grouped_rinv, key, Float64[])
        push!(grouped_rinv[key], s.speedup)
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
        println(md, "| family | regime | blas_threads | geomean speedup (dgeqp3/bsqr_full) |")
        println(md, "|---|---|---:|---:|")
        for key in sort!(collect(keys(grouped_plain)), by = x -> (string(x[1]), x[2], x[3]))
            g = _geomean(grouped_plain[key])
            println(md, "| $(key[1]) | $(key[2]) | $(key[3]) | $(round(g, sigdigits=5)) |")
        end
        println(md, "")
        println(md, "| family | regime | blas_threads | geomean speedup (dgeqp3_trsm/bsqr_rinv) |")
        println(md, "|---|---|---:|---:|")
        for key in sort!(collect(keys(grouped_rinv)), by = x -> (string(x[1]), x[2], x[3]))
            g = _geomean(grouped_rinv[key])
            println(md, "| $(key[1]) | $(key[2]) | $(key[3]) | $(round(g, sigdigits=5)) |")
        end
    end
end

function _validate_pairs(rows::Vector{NamedTuple})
    keyset = Dict{Tuple{Symbol,String,Int,Int,Int,Int},Set{String}}()
    for r in rows
        r.method in ALLOWED_METHODS || error("Unexpected method in benchmark rows: $(r.method)")
        key = (r.family, r.regime, r.m, r.n, r.seed, r.blas_threads)
        get!(keyset, key, Set{String}())
        push!(keyset[key], r.method)
    end
    for (k, methods) in keyset
        expected = k[2] == "short_wide" ? SHORT_WIDE_METHODS : SQUARE_METHODS
        methods == expected || error("Missing benchmark rows for key=$k; expected=$(collect(expected)), present=$(collect(methods))")
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
    expected_rows = _expected_row_count(threads, seeds, families, cases)
    old_threads = BLAS.get_num_threads()

    try
        for nth in threads
            BLAS.set_num_threads(nth)
            println("Configured BLAS threads: ", BLAS.get_num_threads())

            for seed in seeds
                for family in families
                    for c in cases
                        if family === :orthonormal_rows && c.m > c.n
                            continue
                        end

                        case_rng = MersenneTwister(_case_seed(seed, family, c.m, c.n, c.aspect))
                        A = make_matrix(family, c.m, c.n, case_rng)
                        kfull = min(c.m, c.n)
                        plain_rows = bench_pair_ci(
                            A,
                            kfull,
                            norm_recomp_tol;
                            warmup = warmup,
                            samples = samples,
                            bsqr_return_rinv_r12 = false,
                            include_dgeqp3_trsm = false,
                        )
                        rinv_rows_raw = bench_pair_ci(
                            A,
                            kfull,
                            norm_recomp_tol;
                            warmup = warmup,
                            samples = samples,
                            bsqr_return_rinv_r12 = true,
                            include_dgeqp3_trsm = true,
                        )
                        rinv_rows = NamedTuple[]
                        for row in rinv_rows_raw
                            if row.method == BSQR_METHOD_LABEL
                                push!(rinv_rows, (
                                    method = BSQR_RINV_METHOD_LABEL,
                                    tmin = row.tmin,
                                    tmed = row.tmed,
                                    tci_low = row.tci_low,
                                    tci_high = row.tci_high,
                                    alloc = row.alloc,
                                    resid = row.resid,
                                    orth = row.orth,
                                ))
                            elseif row.method == DGEQP3_TRSM_METHOD_LABEL
                                push!(rinv_rows, row)
                            end
                        end
                        rowset = vcat(plain_rows, rinv_rows)

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

    length(rows) == expected_rows ||
        error("Row-count mismatch: expected $expected_rows rows, observed $(length(rows))")
    _validate_pairs(rows)
    _validate_ci_and_quality(rows)
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
        expected_rows,
        length(rows),
    )

    println("Wrote publication benchmark CSV to: $csv_path")
    println("Wrote publication summary to: $summary_path")
    println("Wrote publication metadata to: $metadata_path")
end

run_publication_benchmarks()
