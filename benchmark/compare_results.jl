#!/usr/bin/env julia

using Dates

const DEFAULT_CANDIDATE_DIR = joinpath(@__DIR__, "results")

const TIMINGS_FILE = "timings.csv"
const SWEEP_FILE = "sweep_timings.csv"
const REGIME_FILE = "regime_timings.csv"

const SLOWDOWN_THRESH = 0.02

function _load_csv(path::String)
    isfile(path) || error("Missing CSV: $path")
    lines = readlines(path)
    isempty(lines) && error("Empty CSV: $path")
    cols = split(lines[1], ',')
    rows = Vector{Dict{String,String}}()
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

_pint(r::Dict{String,String}, k::String) = parse(Int, r[k])
_pfloat(r::Dict{String,String}, k::String) = parse(Float64, r[k])

function _index_rows(rows::Vector{Dict{String,String}}, keycols::Vector{String})
    idx = Dict{Tuple,Dict{String,String}}()
    for r in rows
        key = tuple((r[c] for c in keycols)...)
        idx[key] = r
    end
    return idx
end

function _quality_failures(crows, brows, keycols; label::String)
    base_idx = _index_rows(brows, keycols)
    failures = NamedTuple[]
    compared = 0

    for c in crows
        key = tuple((c[k] for k in keycols)...)
        b = get(base_idx, key, nothing)
        b === nothing && continue
        compared += 1

        m = _pint(c, "m")
        n = _pint(c, "n")
        cres = _pfloat(c, "residual")
        corth = _pfloat(c, "orthogonality")
        bres = _pfloat(b, "residual")
        borth = _pfloat(b, "orthogonality")

        res_tol = max(10 * bres, 2.5e3 * eps(Float64) * max(m, n), 1.0e-12)
        orth_tol = max(10 * borth, 2.5e3 * eps(Float64) * max(m, n), 1.0e-12)

        bad = !isfinite(cres) || !isfinite(corth) || cres > res_tol || corth > orth_tol
        if bad
            push!(
                failures,
                (
                    source = label,
                    family = c["family"],
                    method = c["method"],
                    m = m,
                    n = n,
                    residual = cres,
                    residual_tol = res_tol,
                    orthogonality = corth,
                    orthogonality_tol = orth_tol,
                ),
            )
        end
    end

    return failures, compared
end

function _select_key_cases(srows, rrows)
    cases = Vector{Tuple{Symbol,String,Dict{String,String}}}()

    function push_largest!(label::String, rows)
        isempty(rows) && return
        best = rows[argmax(map(r -> _pint(r, "n"), rows))]
        push!(cases, (:sweep, label, best))
    end

    gauss_sw = filter(
        r -> r["method"] == "bsqr_full" && r["family"] == "gaussian" && _pint(r, "m") < _pint(r, "n"),
        srows,
    )
    ill_sw = filter(
        r -> r["method"] == "bsqr_full" && r["family"] == "ill_conditioned" && _pint(r, "m") < _pint(r, "n"),
        srows,
    )
    sq = filter(
        r -> r["method"] == "bsqr_full" && _pint(r, "m") == _pint(r, "n"),
        srows,
    )

    push_largest!("short_wide_gaussian", gauss_sw)
    push_largest!("short_wide_ill_conditioned", ill_sw)
    push_largest!("square_case", sq)

    ortho = filter(
        r -> r["method"] == "bsqr_full" && r["family"] == "orthonormal_rows" && r["regime"] == "fixed_m_vary_n",
        rrows,
    )
    fixed_ms = unique(map(r -> _pint(r, "fixed_value"), ortho))
    for mfix in sort(fixed_ms)
        group = filter(r -> _pint(r, "fixed_value") == mfix, ortho)
        isempty(group) && continue
        best = group[argmax(map(r -> _pint(r, "n"), group))]
        push!(cases, (:regime, "orthonormal_rows_fixed_m_vary_n", best))
    end

    return cases
end

function _slowdown_failures(c_srows, b_srows, c_rrows, b_rrows; thresh::Float64)
    skey = ["family", "aspect", "m", "n", "k", "method"]
    rkey = ["family", "regime", "fixed_value", "var_value", "m", "n", "k", "method"]

    cand_s_idx = _index_rows(c_srows, skey)
    base_s_idx = _index_rows(b_srows, skey)
    cand_r_idx = _index_rows(c_rrows, rkey)
    base_r_idx = _index_rows(b_rrows, rkey)

    cases = _select_key_cases(c_srows, c_rrows)
    failures = NamedTuple[]
    inspected = 0

    for (source, label, c) in cases
        if source === :regime
            key = tuple((c[k] for k in rkey)...)
            b = get(base_r_idx, key, nothing)
            c2 = get(cand_r_idx, key, nothing)
        elseif source === :sweep
            key = tuple((c[k] for k in skey)...)
            b = get(base_s_idx, key, nothing)
            c2 = get(cand_s_idx, key, nothing)
        end
        (b === nothing || c2 === nothing) && continue
        inspected += 1

        bt = _pfloat(b, "tmed_s")
        ct = _pfloat(c2, "tmed_s")
        slowdown = bt == 0.0 ? 0.0 : (ct / bt - 1.0)
        has_ci = haskey(b, "tci_low_s") && haskey(b, "tci_high_s") &&
                 haskey(c2, "tci_low_s") && haskey(c2, "tci_high_s")
        violates = if has_ci
            b_hi = _pfloat(b, "tci_high_s")
            c_lo = _pfloat(c2, "tci_low_s")
            c_lo > b_hi * (1 + thresh)
        else
            slowdown > thresh
        end
        if violates
            push!(
                failures,
                (
                    key_case = label,
                    family = c2["family"],
                    m = _pint(c2, "m"),
                    n = _pint(c2, "n"),
                    baseline_tmed = bt,
                    candidate_tmed = ct,
                    slowdown = slowdown,
                ),
            )
        end
    end

    return failures, inspected
end

function compare_results(
    baseline_dir::String;
    candidate_dir::String = DEFAULT_CANDIDATE_DIR,
    slowdown_thresh::Float64 = SLOWDOWN_THRESH,
)
    b_t = _load_csv(joinpath(baseline_dir, TIMINGS_FILE))
    b_s = _load_csv(joinpath(baseline_dir, SWEEP_FILE))
    b_r = _load_csv(joinpath(baseline_dir, REGIME_FILE))

    c_t = _load_csv(joinpath(candidate_dir, TIMINGS_FILE))
    c_s = _load_csv(joinpath(candidate_dir, SWEEP_FILE))
    c_r = _load_csv(joinpath(candidate_dir, REGIME_FILE))

    qf_t, cmp_t = _quality_failures(c_t, b_t, ["family", "m", "n", "k", "method"]; label = TIMINGS_FILE)
    qf_s, cmp_s = _quality_failures(c_s, b_s, ["family", "aspect", "m", "n", "k", "method"]; label = SWEEP_FILE)
    qf_r, cmp_r = _quality_failures(c_r, b_r, ["family", "regime", "fixed_value", "var_value", "m", "n", "k", "method"]; label = REGIME_FILE)
    quality_failures = vcat(qf_t, qf_s, qf_r)

    perf_failures, inspected = _slowdown_failures(c_s, b_s, c_r, b_r; thresh = slowdown_thresh)
    passed = isempty(quality_failures) && isempty(perf_failures)

    out_md = joinpath(candidate_dir, "validation_report.md")
    out_csv = joinpath(candidate_dir, "guardrail_failures.csv")
    open(out_md, "w") do io
        println(io, "# Validation Report")
        println(io, "")
        println(io, "- Generated: $(Dates.now())")
        println(io, "- Baseline: `$(abspath(baseline_dir))`")
        println(io, "- Candidate: `$(abspath(candidate_dir))`")
        println(io, "- Guardrail: slowdown <= $(100 * slowdown_thresh)% on key cases")
        println(io, "- Result: **", passed ? "PASS" : "FAIL", "**")
        println(io, "")
        println(io, "## Coverage")
        println(io, "- Quality rows compared: $(cmp_t + cmp_s + cmp_r) (`timings=$(cmp_t)`, `sweep=$(cmp_s)`, `regime=$(cmp_r)`)")
        println(io, "- Key performance cases inspected: $inspected")
        println(io, "")
        println(io, "## Performance Failures")
        if isempty(perf_failures)
            println(io, "- none")
        else
            for f in perf_failures
                println(
                    io,
                    "- $(f.key_case): $(f.family) ($(f.m)x$(f.n)) baseline=$(round(f.baseline_tmed, sigdigits=5))s candidate=$(round(f.candidate_tmed, sigdigits=5))s slowdown=$(round(100*f.slowdown, sigdigits=4))%",
                )
            end
        end
        println(io, "")
        println(io, "## Quality Failures")
        if isempty(quality_failures)
            println(io, "- none")
        else
            for f in quality_failures
                println(
                    io,
                    "- $(f.source): $(f.family) $(f.method) ($(f.m)x$(f.n)) residual=$(f.residual) tol=$(f.residual_tol), orth=$(f.orthogonality) tol=$(f.orthogonality_tol)",
                )
            end
        end
    end

    open(out_csv, "w") do io
        println(io, "kind,key_case,source,family,method,m,n,baseline_tmed,candidate_tmed,slowdown,residual,residual_tol,orthogonality,orthogonality_tol")
        for f in perf_failures
            println(io, "performance,$(f.key_case),,$(f.family),bsqr_full,$(f.m),$(f.n),$(f.baseline_tmed),$(f.candidate_tmed),$(f.slowdown),,,,")
        end
        for f in quality_failures
            println(io, "quality,,$(f.source),$(f.family),$(f.method),$(f.m),$(f.n),,,,$(f.residual),$(f.residual_tol),$(f.orthogonality),$(f.orthogonality_tol)")
        end
    end

    println("Wrote validation report: $out_md")
    println("Wrote guardrail CSV: $out_csv")
    return passed
end

function main()
    baseline_dir = if length(ARGS) >= 1
        ARGS[1]
    else
        strip(get(ENV, "BS_BASELINE_DIR", ""))
    end
    isempty(baseline_dir) && error("Provide baseline directory as arg1 or set BS_BASELINE_DIR")

    candidate_dir = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_CANDIDATE_DIR
    threshold = parse(Float64, get(ENV, "BS_SLOWDOWN_THRESH", string(SLOWDOWN_THRESH)))

    ok = compare_results(baseline_dir; candidate_dir = candidate_dir, slowdown_thresh = threshold)
    ok || exit(1)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
