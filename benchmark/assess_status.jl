#!/usr/bin/env julia

using Dates
using LinearAlgebra
using Random
using Statistics

using BSPivotQR

include("bench_common.jl")
using .BenchCommon
include("validation_pipeline_common.jl")
using .ValidationPipelineCommon

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULTS_DIR = joinpath(@__DIR__, "results")
const PROFILE_BREAKDOWN_CSV = joinpath(RESULTS_DIR, "profile_breakdown.csv")
const TIMINGS_CSV = joinpath(RESULTS_DIR, "timings.csv")
const SWEEP_CSV = joinpath(RESULTS_DIR, "sweep_timings.csv")
const REGIME_CSV = joinpath(RESULTS_DIR, "regime_timings.csv")
const GUARDRAIL_CSV = joinpath(RESULTS_DIR, "guardrail_failures.csv")
const VALIDATION_MD = joinpath(RESULTS_DIR, "validation_report.md")

const ASSESSMENT_RAW_CSV = joinpath(RESULTS_DIR, "assessment_raw.csv")
const STATUS_TRIAGE_PREFIX = "status_triage_"
const HOTSPOT_PREFIX = "profile_hotspots_"

const SUMMARY_PAIRS = (
    (
        source = "timings",
        csv_path = TIMINGS_CSV,
        keycols = ["family", "m", "n", "k"],
        md_path = joinpath(RESULTS_DIR, "summary.md"),
    ),
    (
        source = "sweep",
        csv_path = SWEEP_CSV,
        keycols = ["family", "aspect", "m", "n", "k"],
        md_path = joinpath(RESULTS_DIR, "sweep_summary.md"),
    ),
    (
        source = "regime",
        csv_path = REGIME_CSV,
        keycols = ["family", "regime", "fixed_value", "var_value", "m", "n", "k"],
        md_path = joinpath(RESULTS_DIR, "regime_summary.md"),
    ),
)

@inline _f64(x::AbstractString) = parse(Float64, x)
@inline _i64(x::AbstractString) = parse(Int, x)

function _display_path(path::String)
    return relpath(abspath(path), REPO_ROOT)
end

function _timestamp_tag()
    return Dates.format(now(), "yyyymmdd_HHMMSS")
end

function _cpu_model()
    try
        return getproperty(first(Sys.cpu_info()), :model)
    catch
        return "unknown"
    end
end

function _load_csv_rows(path::String)
    isfile(path) || error("Missing CSV: $path")
    lines = readlines(path)
    isempty(lines) && error("Empty CSV: $path")
    cols = split(lines[1], ',')
    rows = Dict{String,String}[]
    for line in lines[2:end]
        isempty(strip(line)) && continue
        vals = split(line, ',')
        length(vals) == length(cols) || error("Malformed CSV row in $path: $line")
        row = Dict{String,String}()
        for (c, v) in zip(cols, vals)
            row[c] = v
        end
        push!(rows, row)
    end
    return rows
end

_row_key(r::Dict{String,String}, keycols::Vector{String}) = tuple((r[c] for c in keycols)...)

function _pair_rows(rows::Vector{Dict{String,String}}, keycols::Vector{String}; source::String)
    bykey = Dict{Tuple,Dict{String,Dict{String,String}}}()
    for r in rows
        key = _row_key(r, keycols)
        slot = get!(bykey, key, Dict{String,Dict{String,String}}())
        slot[r["method"]] = r
    end

    pairs = NamedTuple[]
    for (key, methods) in bykey
        haskey(methods, BSQR_METHOD_LABEL) || continue
        haskey(methods, DGEQP3_METHOD_LABEL) || continue
        bs = methods[BSQR_METHOD_LABEL]
        dg = methods[DGEQP3_METHOD_LABEL]
        bs_t = _f64(bs["tmed_s"])
        dg_t = _f64(dg["tmed_s"])
        bs_t > 0.0 || continue
        speed = dg_t / bs_t
        push!(
            pairs,
            (
                source = source,
                key = key,
                case_label = _case_label(source, bs),
                bs_tmed = bs_t,
                dg_tmed = dg_t,
                speedup = speed,
                slowdown = bs_t / max(dg_t, eps(Float64)),
                bs_alloc = _f64(bs["alloc_bytes"]),
                dg_alloc = _f64(dg["alloc_bytes"]),
                bs_residual = _f64(bs["residual"]),
                bs_orth = _f64(bs["orthogonality"]),
                dg_residual = _f64(dg["residual"]),
                dg_orth = _f64(dg["orthogonality"]),
                bs_tci_low = haskey(bs, "tci_low_s") ? _f64(bs["tci_low_s"]) : bs_t,
                bs_tci_high = haskey(bs, "tci_high_s") ? _f64(bs["tci_high_s"]) : bs_t,
                dg_tci_low = haskey(dg, "tci_low_s") ? _f64(dg["tci_low_s"]) : dg_t,
                dg_tci_high = haskey(dg, "tci_high_s") ? _f64(dg["tci_high_s"]) : dg_t,
            ),
        )
    end
    return pairs
end

function _case_label(source::String, r::Dict{String,String})
    if source == "timings"
        return "$(r["family"]) m=$(r["m"]) n=$(r["n"]) k=$(r["k"])"
    elseif source == "sweep"
        return "$(r["family"]) aspect=$(r["aspect"]) m=$(r["m"]) n=$(r["n"]) k=$(r["k"])"
    else
        return "$(r["family"]) regime=$(r["regime"]) fixed=$(r["fixed_value"]) var=$(r["var_value"]) m=$(r["m"]) n=$(r["n"]) k=$(r["k"])"
    end
end

function _summary_consistency(csv_rows::Int, md_path::String)
    isfile(md_path) || return (ok = false, reason = "missing markdown summary")
    lines = readlines(md_path)
    md_rows = count(l -> occursin("|", l) && (occursin("bsqr_full", l) || occursin("dgeqp3", l)), lines)
    return (ok = md_rows == csv_rows, reason = "csv_rows=$csv_rows md_rows=$md_rows")
end

function _quality_sanity(rows::Vector{Dict{String,String}})
    failures = NamedTuple[]
    for r in rows
        m = _i64(r["m"])
        n = _i64(r["n"])
        resid = _f64(r["residual"])
        orth = _f64(r["orthogonality"])
        tol = max(2.5e3 * eps(Float64) * max(m, n), 1e-12)
        if !isfinite(resid) || !isfinite(orth) || resid > tol || orth > tol
            push!(failures, (family = r["family"], m = m, n = n, method = r["method"], residual = resid, orth = orth, tol = tol))
        end
    end
    return failures
end

function _ci_width_warnings(rows::Vector{Dict{String,String}})
    warnings = NamedTuple[]
    for r in rows
        haskey(r, "tci_low_s") || continue
        haskey(r, "tci_high_s") || continue
        tmed = _f64(r["tmed_s"])
        width = _f64(r["tci_high_s"]) - _f64(r["tci_low_s"])
        rel = tmed > 0.0 ? width / tmed : Inf
        if rel > 1.0
            push!(
                warnings,
                (
                    family = r["family"],
                    method = r["method"],
                    m = _i64(r["m"]),
                    n = _i64(r["n"]),
                    rel_width = rel,
                ),
            )
        end
    end
    return warnings
end

function _write_hotspots(ts::String)
    rows = _load_csv_rows(PROFILE_BREAKDOWN_CSV)
    outpath = joinpath(RESULTS_DIR, HOTSPOT_PREFIX * ts * ".csv")
    open(outpath, "w") do io
        println(io, "family,case,tol,total_med_s,pivot_pct,householder_pct,apply_pct,w_update_pct,downdate_pct,recompute_count,residual,orthogonality")
        for r in rows
            total = _f64(r["total_med_s"])
            t = max(total, eps(Float64))
            pivot_pct = 100.0 * _f64(r["pivot_s"]) / t
            hh_pct = 100.0 * _f64(r["householder_s"]) / t
            apply_pct = 100.0 * _f64(r["apply_reflector_s"]) / t
            wup_pct = 100.0 * _f64(r["w_update_s"]) / t
            down_pct = 100.0 * _f64(r["norm_downdate_s"]) / t
            println(
                io,
                "$(r["family"]),$(r["case"]),$(r["tol"]),$(r["total_med_s"]),$pivot_pct,$hh_pct,$apply_pct,$wup_pct,$down_pct,$(r["recompute_count"]),$(r["residual"]),$(r["orthogonality"])",
            )
        end
    end
    return outpath
end

function _run_assessment_stress_checks()
    rng = MersenneTwister(20260225)
    checks = NamedTuple[]

    function quality_ok(A::Matrix{Float64}; desc::String)
        F = bsqr(A; check = true)
        rb, qb = residual_bs(A, F)
        tol = max(5.0e3 * eps(Float64) * max(size(A)...), 5e-12)
        ok = all(isfinite, F.factors) && all(isfinite, F.tau) && rb <= tol && qb <= tol
        push!(checks, (name = desc, ok = ok, residual = rb, orth = qb, tol = tol, ksteps = F.ksteps))
    end

    m, n = 64, 48
    r = min(m, n)
    U = Matrix(qr(randn(rng, m, r)).Q)
    V = Matrix(qr(randn(rng, n, r)).Q)
    for kappa in (1.0e6, 1.0e8, 1.0e10, 1.0e12)
        s = exp.(range(0.0, stop = -log(kappa), length = r))
        quality_ok(U * Diagonal(s) * V'; desc = "ill_conditioned_kappa=$(kappa)")
    end

    A_rank = randn(rng, 48, 32)
    A_rank[:, 5] .= A_rank[:, 4] .+ 1e-13 * randn(rng, 48)
    A_rank[:, 7] .= A_rank[:, 2]
    quality_ok(A_rank; desc = "near_rank_deficient_duplicate_columns")

    A_tie = randn(rng, 40, 20)
    A_tie[:, 3] .= A_tie[:, 1] + 1e-14 * randn(rng, 40)
    A_tie[:, 4] .= A_tie[:, 2] + 1e-14 * randn(rng, 40)
    quality_ok(A_tie; desc = "near_tie_columns")

    A_scale = randn(rng, 60, 24)
    scales = exp10.(range(-12, 12, length = size(A_scale, 2)))
    A_scale .*= reshape(scales, 1, :)
    quality_ok(A_scale; desc = "adversarial_column_scaling")

    passed = all(c.ok for c in checks)
    return passed, checks
end

function _write_status_triage(
    ts::String,
    pairs_by_source::Dict{String,Vector{NamedTuple}},
    consistency::Dict{String,NamedTuple},
    quality_failures::Dict{String,Vector{NamedTuple}},
    ci_warnings::Dict{String,Vector{NamedTuple}},
    stress_checks::Vector{NamedTuple},
    stress_passed::Bool,
    hotspots_path::String,
    baseline_dir::Union{Nothing,String},
)
    status_path = joinpath(RESULTS_DIR, STATUS_TRIAGE_PREFIX * ts * ".md")
    all_pairs = vcat(values(pairs_by_source)...)
    sorted_best = sort(all_pairs, by = p -> p.speedup, rev = true)
    sorted_worst = sort(all_pairs, by = p -> p.speedup)

    open(status_path, "w") do io
        println(io, "# Status Triage Snapshot")
        println(io, "")
        println(io, "- Generated: $(Dates.now())")
        println(io, "- Timestamp tag: `$ts`")
        println(io, "- Julia: `$(VERSION)`")
        println(io, "- CPU: `$(_cpu_model())`")
        println(io, "- CPU threads: `$(Sys.CPU_THREADS)`")
        println(io, "- BLAS threads: `$(BLAS.get_num_threads())`")
        println(io, "- BLAS config: `$(backend_string())`")
        println(io, "- Baseline dir: `$(baseline_dir === nothing ? "none" : _display_path(baseline_dir))`")
        println(io, "- Validation report: `$(_display_path(VALIDATION_MD))`")
        println(io, "- Hotspot CSV: `$(_display_path(hotspots_path))`")
        println(io, "")
        println(io, "## Pass/Fail Status")
        println(io, "- Tests + quick validation pipeline: **PASS**")
        println(io, "- Guardrail failures file present: `$(isfile(GUARDRAIL_CSV))`")
        println(io, "- Stress checks: **$(stress_passed ? "PASS" : "FAIL")**")
        println(io, "")
        println(io, "## Speedup Summary (`dgeqp3/bsqr`)")
        println(io, "")
        println(io, "| source | cases | faster | slower | median speedup | mean speedup |")
        println(io, "|---|---:|---:|---:|---:|---:|")
        for (source, pairs) in sort(collect(pairs_by_source), by = first)
            speeds = [p.speedup for p in pairs]
            faster = count(>(1.0), speeds)
            slower = count(<(1.0), speeds)
            med = isempty(speeds) ? NaN : median(speeds)
            avg = isempty(speeds) ? NaN : mean(speeds)
            println(io, "| $source | $(length(speeds)) | $faster | $slower | $(round(med, sigdigits=5)) | $(round(avg, sigdigits=5)) |")
        end
        println(io, "")
        println(io, "## Top 10 Best Cases")
        for p in first(sorted_best, min(10, length(sorted_best)))
            println(io, "- `$(p.source)` $(p.case_label): speedup=$(round(p.speedup, sigdigits=5))")
        end
        println(io, "")
        println(io, "## Top 10 Worst Cases")
        for p in first(sorted_worst, min(10, length(sorted_worst)))
            println(io, "- `$(p.source)` $(p.case_label): speedup=$(round(p.speedup, sigdigits=5))")
        end
        println(io, "")
        println(io, "## Sanity Checks")
        for source in sort(collect(keys(pairs_by_source)))
            c = consistency[source]
            qbad = length(quality_failures[source])
            cbad = length(ci_warnings[source])
            println(io, "- `$source` summary consistency: **$(c.ok ? "PASS" : "FAIL")** (`$(c.reason)`)")
            println(io, "- `$source` quality anomalies: `$qbad`")
            println(io, "- `$source` CI-width warnings (rel width > 1): `$cbad`")
        end
        println(io, "")
        println(io, "## Stress Check Details")
        for c in stress_checks
            println(io, "- $(c.name): $(c.ok ? "PASS" : "FAIL"), residual=$(round(c.residual, sigdigits=5)), orth=$(round(c.orth, sigdigits=5)), tol=$(round(c.tol, sigdigits=5)), ksteps=$(c.ksteps)")
        end
    end

    return status_path
end

function _profile_case_percentages(
    A::Matrix{Float64},
    tol::Float64;
    reps::Int,
)
    m, n = size(A)
    k = min(m, n)
    totals = Float64[]
    pivots = Float64[]
    house = Float64[]
    apply = Float64[]
    wup = Float64[]
    down = Float64[]
    recomp = Float64[]

    for _ in 1:reps
        Awork = copy(A)
        ws = BSPivotQR.BSWorkspace(m, n, k)
        tau = zeros(Float64, k)
        jpvt = collect(1:n)
        stats = BSPivotQR.BSKernelStats()
        t0 = time_ns()
        BSPivotQR._bsqr_kernel!(
            Awork,
            tau,
            jpvt,
            ws,
            k;
            rank_stop = false,
            norm_recomp_tol = tol,
            kernel_stats = stats,
        )
        total_s = (time_ns() - t0) * 1e-9
        push!(totals, total_s)
        push!(pivots, stats.pivot_select_ns * 1e-9)
        push!(house, stats.householder_ns * 1e-9)
        push!(apply, stats.apply_reflector_ns * 1e-9)
        push!(wup, stats.w_update_ns * 1e-9)
        push!(down, stats.norm_downdate_ns * 1e-9)
        push!(recomp, stats.recompute_count)
    end

    total_med = median(totals)
    t = max(total_med, eps(Float64))
    return (
        total_med_s = total_med,
        pivot_pct = 100.0 * median(pivots) / t,
        householder_pct = 100.0 * median(house) / t,
        apply_pct = 100.0 * median(apply) / t,
        w_update_pct = 100.0 * median(wup) / t,
        downdate_pct = 100.0 * median(down) / t,
        recompute_count = Int(round(median(recomp))),
    )
end

function _parse_tolerance_list()
    dflt = [sqrt(eps(Float64)), 1.0e-10, 1.0e-12, 0.0]
    return parse_float_list("BS_ASSESS_TOLS", dflt)
end

function _deep_case_grid(quick::Bool)
    cases = NamedTuple[]
    seen = Set{Tuple{String,Int,Int}}()
    function push_case(regime::String, m::Int, n::Int)
        key = (regime, m, n)
        key in seen && return
        push!(seen, key)
        push!(cases, (regime = regime, m = m, n = n))
    end

    square_ns = quick ? [64, 128] : [64, 128, 256, 512, 1024]
    short_ns = quick ? [512, 1024] : [256, 512, 1024, 2048, 4096]
    tall_ms = quick ? [512, 1024] : [256, 512, 1024, 2048, 4096]
    fixed_ms = quick ? [32, 64] : [32, 64, 128]
    fixed_ns = quick ? [128, 256] : [128, 256, 512]
    vary_ns = quick ? [128, 512, 2048] : [64, 128, 256, 512, 1024, 2048, 4096, 8192]
    vary_ms = quick ? [32, 128, 512] : [32, 64, 128, 256, 512, 1024, 2048, 3072]

    for n in square_ns
        push_case("square", n, n)
    end
    for n in short_ns
        push_case("short_wide", 64, n)
    end
    for m in tall_ms
        push_case("tall_skinny", m, 64)
    end
    for mfix in fixed_ms, n in vary_ns
        push_case("fixed_m_vary_n", mfix, n)
    end
    for nfix in fixed_ns, m in vary_ms
        push_case("fixed_n_vary_m", m, nfix)
    end

    return sort(cases, by = c -> (c.regime, c.m, c.n))
end

function _write_assessment_raw(ts::String)
    quick = get(ENV, "BS_ASSESS_DEEP_QUICK", "1") == "1"
    warmup = parse_env_int("BS_ASSESS_WARMUP", quick ? 1 : 2; minval = 0)
    samples = parse_env_int("BS_ASSESS_SAMPLES", quick ? 10 : 24)
    profile_reps = parse_env_int("BS_ASSESS_PROFILE_REPS", quick ? 2 : 4)
    tols = _parse_tolerance_list()
    families = parse_symbol_list(
        "BS_ASSESS_FAMILIES",
        [:gaussian, :ill_conditioned, :orthonormal_rows],
        [:gaussian, :ill_conditioned, :orthonormal_rows],
    )
    cases = _deep_case_grid(quick)

    rng = MersenneTwister(20260225)
    open(ASSESSMENT_RAW_CSV, "w") do io
        println(io, "timestamp,family,regime,m,n,k,tol,method,tmin_s,tmed_s,tci_low_s,tci_high_s,alloc_bytes,residual,orthogonality,pivot_pct,householder_pct,apply_pct,w_update_pct,downdate_pct,recompute_count")
        for family in families
            for c in cases
                m, n = c.m, c.n
                if family === :orthonormal_rows && m > n
                    continue
                end
                A = make_matrix(family, m, n, rng)
                kfull = min(m, n)
                for tol in tols
                    p = _profile_case_percentages(A, tol; reps = profile_reps)
                    for row in bench_pair_ci(A, kfull, tol; warmup = warmup, samples = samples)
                        if row.method == BSQR_METHOD_LABEL
                            println(
                                io,
                                "$ts,$family,$(c.regime),$m,$n,$kfull,$tol,$(row.method),$(row.tmin),$(row.tmed),$(row.tci_low),$(row.tci_high),$(row.alloc),$(row.resid),$(row.orth),$(p.pivot_pct),$(p.householder_pct),$(p.apply_pct),$(p.w_update_pct),$(p.downdate_pct),$(p.recompute_count)",
                            )
                        else
                            println(
                                io,
                                "$ts,$family,$(c.regime),$m,$n,$kfull,$tol,$(row.method),$(row.tmin),$(row.tmed),$(row.tci_low),$(row.tci_high),$(row.alloc),$(row.resid),$(row.orth),NaN,NaN,NaN,NaN,NaN,NaN",
                            )
                        end
                    end
                end
            end
        end
    end
    return ASSESSMENT_RAW_CSV
end

function _run_stage1(ts::String)
    check_backend()
    configure_blas_threads()

    mkpath(RESULTS_DIR)
    baseline_dir = String(strip(get(ENV, "BS_BASELINE_DIR", joinpath(RESULTS_DIR, "baseline_codex_check"))))
    baseline_arg = isdir(baseline_dir) ? baseline_dir : nothing

    envvars = Dict(
        "BS_BLAS_THREADS" => "4",
        "BS_USE_ACCELERATE" => "1",
        "BS_REQUIRE_ACCELERATE" => "1",
        "BS_NORM_RECOMP_TOL" => string(DEFAULT_NORM_RECOMP_TOL),
        "BS_QUICK" => "1",
        "BS_SWEEP_QUICK" => "1",
        "BS_REGIME_QUICK" => "1",
        "BS_PROFILE_QUICK" => "1",
    )

    run_validation_pipeline(envvars; baseline_dir = baseline_arg)

    hotspots_path = _write_hotspots(ts)
    stress_passed, stress_checks = _run_assessment_stress_checks()

    pairs_by_source = Dict{String,Vector{NamedTuple}}()
    consistency = Dict{String,NamedTuple}()
    quality_failures = Dict{String,Vector{NamedTuple}}()
    ci_warnings = Dict{String,Vector{NamedTuple}}()
    for spec in SUMMARY_PAIRS
        rows = _load_csv_rows(spec.csv_path)
        pairs = _pair_rows(rows, spec.keycols; source = spec.source)
        pairs_by_source[spec.source] = pairs
        consistency[spec.source] = _summary_consistency(length(rows), spec.md_path)
        quality_failures[spec.source] = _quality_sanity(rows)
        ci_warnings[spec.source] = _ci_width_warnings(rows)
    end

    status_path = _write_status_triage(
        ts,
        pairs_by_source,
        consistency,
        quality_failures,
        ci_warnings,
        stress_checks,
        stress_passed,
        hotspots_path,
        baseline_arg,
    )
    return (
        status_path = status_path,
        hotspots_path = hotspots_path,
        stress_passed = stress_passed,
    )
end

function main()
    ts = _timestamp_tag()
    do_stage1 = get(ENV, "BS_ASSESS_RUN_STAGE1", "1") == "1"
    do_stage2 = get(ENV, "BS_ASSESS_RUN_STAGE2", "1") == "1"

    if do_stage1
        s1 = _run_stage1(ts)
        println("Wrote triage status: ", s1.status_path)
        println("Wrote hotspot table: ", s1.hotspots_path)
        s1.stress_passed || error("Stage-1 stress checks failed; investigate before continuing.")
    else
        println("Skipping stage 1 (BS_ASSESS_RUN_STAGE1=0)")
    end

    if do_stage2
        raw_path = _write_assessment_raw(ts)
        println("Wrote assessment raw dataset: ", raw_path)
    else
        println("Skipping stage 2 (BS_ASSESS_RUN_STAGE2=0)")
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
