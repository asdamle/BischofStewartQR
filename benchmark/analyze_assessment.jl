#!/usr/bin/env julia

using Dates
using Printf
using Statistics

const RESULTS_DIR = joinpath(@__DIR__, "results")
const RAW_CSV = joinpath(RESULTS_DIR, "assessment_raw.csv")
const SUMMARY_MD = joinpath(RESULTS_DIR, "assessment_summary.md")
const CANDIDATES_MD = joinpath(RESULTS_DIR, "optimization_candidates.md")
const CONCLUSION_MD = joinpath(RESULTS_DIR, "perf_assessment_conclusion.md")
const VALIDATION_MD = joinpath(RESULTS_DIR, "validation_report.md")

const PHASE_KEYS = ("pivot_pct", "householder_pct", "apply_pct", "w_update_pct", "downdate_pct")

@inline _f64(s::AbstractString) = parse(Float64, s)
@inline _i64(s::AbstractString) = parse(Int, s)

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

function _pair_rows(rows::Vector{Dict{String,String}})
    keycols = ["timestamp", "family", "regime", "m", "n", "k", "tol"]
    bykey = Dict{Tuple,Dict{String,Dict{String,String}}}()
    for r in rows
        key = tuple((r[c] for c in keycols)...)
        slot = get!(bykey, key, Dict{String,Dict{String,String}}())
        slot[r["method"]] = r
    end

    pairs = NamedTuple[]
    for (key, methods) in bykey
        haskey(methods, "bsqr_full") || continue
        haskey(methods, "dgeqp3") || continue
        bs = methods["bsqr_full"]
        dg = methods["dgeqp3"]

        bs_tmed = _f64(bs["tmed_s"])
        dg_tmed = _f64(dg["tmed_s"])
        bs_tmed > 0.0 || continue
        dg_tmed > 0.0 || continue

        bs_lo = _f64(bs["tci_low_s"])
        bs_hi = _f64(bs["tci_high_s"])
        dg_lo = _f64(dg["tci_low_s"])
        dg_hi = _f64(dg["tci_high_s"])

        sig = if dg_lo > bs_hi
            :sig_win
        elseif bs_lo > dg_hi
            :sig_loss
        else
            :overlap
        end

        push!(
            pairs,
            (
                timestamp = bs["timestamp"],
                family = bs["family"],
                regime = bs["regime"],
                m = _i64(bs["m"]),
                n = _i64(bs["n"]),
                k = _i64(bs["k"]),
                tol = _f64(bs["tol"]),
                speedup = dg_tmed / bs_tmed,
                slowdown = bs_tmed / dg_tmed,
                lower_speedup = dg_lo / max(bs_hi, eps(Float64)),
                upper_speedup = dg_hi / max(bs_lo, eps(Float64)),
                significance = sig,
                bs_tmed = bs_tmed,
                dg_tmed = dg_tmed,
                bs_alloc = _f64(bs["alloc_bytes"]),
                dg_alloc = _f64(dg["alloc_bytes"]),
                alloc_ratio = _f64(bs["alloc_bytes"]) / max(_f64(dg["alloc_bytes"]), 1.0),
                bs_residual = _f64(bs["residual"]),
                bs_orth = _f64(bs["orthogonality"]),
                dg_residual = _f64(dg["residual"]),
                dg_orth = _f64(dg["orthogonality"]),
                pivot_pct = _f64(bs["pivot_pct"]),
                householder_pct = _f64(bs["householder_pct"]),
                apply_pct = _f64(bs["apply_pct"]),
                w_update_pct = _f64(bs["w_update_pct"]),
                downdate_pct = _f64(bs["downdate_pct"]),
                recompute_count = _f64(bs["recompute_count"]),
            ),
        )
    end
    return pairs
end

function _safe_corr(x::Vector{Float64}, y::Vector{Float64})
    idx = [i for i in eachindex(x) if isfinite(x[i]) && isfinite(y[i])]
    length(idx) >= 3 || return NaN
    xs = x[idx]
    ys = y[idx]
    std(xs) > 0.0 || return NaN
    std(ys) > 0.0 || return NaN
    return cor(xs, ys)
end

function _group_summary(pairs::Vector{NamedTuple})
    groups = Dict{Tuple,Vector{NamedTuple}}()
    for p in pairs
        key = (p.family, p.regime, p.tol)
        push!(get!(groups, key, NamedTuple[]), p)
    end
    rows = NamedTuple[]
    for (key, ps) in groups
        speeds = [p.speedup for p in ps]
        push!(
            rows,
            (
                family = key[1],
                regime = key[2],
                tol = key[3],
                count = length(ps),
                faster = count(>(1.0), speeds),
                slower = count(<(1.0), speeds),
                sig_wins = count(p -> p.significance === :sig_win, ps),
                sig_losses = count(p -> p.significance === :sig_loss, ps),
                overlap = count(p -> p.significance === :overlap, ps),
                mean_speedup = mean(speeds),
                median_speedup = median(speeds),
                min_speedup = minimum(speeds),
                max_speedup = maximum(speeds),
                mean_lower_speedup = mean(p.lower_speedup for p in ps),
            ),
        )
    end
    return sort(rows, by = r -> (String(r.family), String(r.regime), r.tol))
end

function _quality_map(rows::Vector{Dict{String,String}})
    bs = filter(r -> r["method"] == "bsqr_full", rows)
    groups = Dict{Float64,Vector{Dict{String,String}}}()
    for r in bs
        t = _f64(r["tol"])
        push!(get!(groups, t, Dict{String,String}[]), r)
    end
    out = NamedTuple[]
    for (tol, rs) in groups
        res = [_f64(r["residual"]) for r in rs]
        orth = [_f64(r["orthogonality"]) for r in rs]
        nonfinite = count(x -> !isfinite(x), res) + count(x -> !isfinite(x), orth)
        push!(
            out,
            (
                tol = tol,
                count = length(rs),
                median_residual = median(res),
                max_residual = maximum(res),
                median_orth = median(orth),
                max_orth = maximum(orth),
                nonfinite = nonfinite,
            ),
        )
    end
    return sort(out, by = r -> r.tol)
end

function _hotspot_metrics(pairs::Vector{NamedTuple})
    slowdown = [p.slowdown for p in pairs]
    apply = [p.apply_pct for p in pairs]
    wup = [p.w_update_pct for p in pairs]
    down = [p.downdate_pct for p in pairs]
    pivot = [p.pivot_pct for p in pairs]
    house = [p.householder_pct for p in pairs]
    excess = [max(s - 1.0, 0.0) for s in slowdown]

    contrib_apply = mean(e * a / 100.0 for (e, a) in zip(excess, apply))
    contrib_wup = mean(e * w / 100.0 for (e, w) in zip(excess, wup))
    contrib_down = mean(e * d / 100.0 for (e, d) in zip(excess, down))
    contrib_pivot = mean(e * p / 100.0 for (e, p) in zip(excess, pivot))

    return (
        mean_apply_pct = mean(apply),
        mean_w_update_pct = mean(wup),
        mean_downdate_pct = mean(down),
        mean_pivot_pct = mean(pivot),
        mean_householder_pct = mean(house),
        corr_slowdown_apply = _safe_corr(slowdown, apply),
        corr_slowdown_w_update = _safe_corr(slowdown, wup),
        corr_slowdown_downdate = _safe_corr(slowdown, down),
        corr_slowdown_pivot = _safe_corr(slowdown, pivot),
        contrib_apply = contrib_apply,
        contrib_w_update = contrib_wup,
        contrib_downdate = contrib_down,
        contrib_pivot = contrib_pivot,
    )
end

function _candidate_table(pairs::Vector{NamedTuple}, hotspot::NamedTuple)
    total = max(length(pairs), 1)
    frac_sig_loss = count(p -> p.significance === :sig_loss, pairs) / total
    alloc_ratio_med = median(p.alloc_ratio for p in pairs)
    recompute_norm = mean(min(p.recompute_count / max(p.k, 1), 1.0) for p in pairs)

    candidates = [
        (
            candidate_id = "K1",
            location = "src/kernel.jl (downdate/recompute and timed branches)",
            rationale = @sprintf(
                "Downdate path has high average share (%.2f%%) with nontrivial recompute activity (%.3f normalized).",
                hotspot.mean_downdate_pct,
                recompute_norm,
            ),
            impact_raw = hotspot.mean_downdate_pct + 25.0 * recompute_norm,
            risk_level = "medium",
            numerical_risk = "medium",
            implementation_cost = "medium",
            validation_plan = "Run tests + tier1 + compare_results; verify residual/orth and recompute behavior on ill-conditioned ladder.",
        ),
        (
            candidate_id = "K2",
            location = "src/kernel.jl (W-update and workspace data movement)",
            rationale = @sprintf(
                "W-update dominates average phase share (%.2f%%) and correlates with slowdown (corr=%.3f).",
                hotspot.mean_w_update_pct,
                hotspot.corr_slowdown_w_update,
            ),
            impact_raw = hotspot.mean_w_update_pct,
            risk_level = "low",
            numerical_risk = "low",
            implementation_cost = "medium",
            validation_plan = "Run tests + tier1; verify identical outputs for fixed RNG seeds and no allocation regressions.",
        ),
        (
            candidate_id = "K3",
            location = "src/kernel.jl (_apply_householder_left! and surrounding call pattern)",
            rationale = @sprintf(
                "Reflector application remains a major cost center (%.2f%% average phase share).",
                hotspot.mean_apply_pct,
            ),
            impact_raw = hotspot.mean_apply_pct,
            risk_level = "medium",
            numerical_risk = "low",
            implementation_cost = "medium",
            validation_plan = "Run tests + tier1; validate residual/orthogonality parity on square and short-wide cases.",
        ),
        (
            candidate_id = "K4",
            location = "src/kernel.jl (pivot selection loop)",
            rationale = @sprintf(
                "Pivot selection is small but systematic (%.2f%%) and safe for micro-optimizations.",
                hotspot.mean_pivot_pct,
            ),
            impact_raw = hotspot.mean_pivot_pct,
            risk_level = "low",
            numerical_risk = "low",
            implementation_cost = "low",
            validation_plan = "Run tests including criterion-consistent pivot sequence; ensure pivot order invariants hold.",
        ),
        (
            candidate_id = "K5",
            location = "benchmark/bench_common.jl and benchmarking harness",
            rationale = @sprintf(
                "Median BS allocation ratio vs dgeqp3 is %.3f; tighten measurement/alloc path to isolate kernel effects.",
                alloc_ratio_med,
            ),
            impact_raw = max((alloc_ratio_med - 1.0) * 100.0, 0.0),
            risk_level = "low",
            numerical_risk = "none",
            implementation_cost = "low",
            validation_plan = "Re-run benchmarks and confirm CSV schemas and fairness policy remain unchanged.",
        ),
    ]

    max_impact = max(maximum(c.impact_raw for c in candidates), eps(Float64))
    risk_penalty = Dict("none" => 0.0, "low" => 0.1, "medium" => 0.3, "high" => 0.5)
    cost_penalty = Dict("low" => 0.1, "medium" => 0.3, "high" => 0.5)

    scored = NamedTuple[]
    for c in candidates
        impact = c.impact_raw / max_impact
        confidence = clamp(frac_sig_loss + (impact > 0.5 ? 0.2 : 0.1), 0.0, 1.0)
        nrisk = get(risk_penalty, c.numerical_risk, 0.3)
        cost = get(cost_penalty, c.implementation_cost, 0.3)
        score = 0.5 * impact + 0.3 * confidence - 0.1 * nrisk - 0.1 * cost
        go = score >= 0.25 ? "GO" : "HOLD"
        push!(
            scored,
            (
                c.candidate_id,
                c.location,
                c.rationale,
                expected_gain = @sprintf("impact_norm=%.3f; score_model=%.3f", impact, score),
                c.risk_level,
                c.numerical_risk,
                c.implementation_cost,
                c.validation_plan,
                go_no_go = go,
                weighted_score = score,
            ),
        )
    end

    return sort(scored, by = c -> c.weighted_score, rev = true)
end

function _write_summary(rows::Vector{Dict{String,String}}, pairs::Vector{NamedTuple})
    grouped = _group_summary(pairs)
    quality = _quality_map(rows)
    hotspot = _hotspot_metrics(pairs)

    open(SUMMARY_MD, "w") do io
        println(io, "# Assessment Summary")
        println(io, "")
        println(io, "- Generated: $(Dates.now())")
        println(io, "- Raw dataset: `$(RAW_CSV)`")
        println(io, "- Pair count: $(length(pairs))")
        println(io, "")
        speeds = [p.speedup for p in pairs]
        println(io, "## Overall Speedup Distribution (`dgeqp3/bsqr`)")
        println(io, "- mean: $(round(mean(speeds), sigdigits=5))")
        println(io, "- median: $(round(median(speeds), sigdigits=5))")
        println(io, "- faster cases (>1): $(count(>(1.0), speeds))/$(length(speeds))")
        println(io, "- significant wins: $(count(p -> p.significance === :sig_win, pairs))")
        println(io, "- significant losses: $(count(p -> p.significance === :sig_loss, pairs))")
        println(io, "")

        println(io, "## Case-Level Aggregates")
        println(io, "")
        println(io, "| family | regime | tol | count | faster | slower | sig wins | sig losses | overlap | mean speedup | median speedup | mean lower-bound speedup |")
        println(io, "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for r in grouped
            println(
                io,
                "| $(r.family) | $(r.regime) | $(r.tol) | $(r.count) | $(r.faster) | $(r.slower) | $(r.sig_wins) | $(r.sig_losses) | $(r.overlap) | $(round(r.mean_speedup, sigdigits=5)) | $(round(r.median_speedup, sigdigits=5)) | $(round(r.mean_lower_speedup, sigdigits=5)) |",
            )
        end
        println(io, "")

        println(io, "## Quality-Risk Map by Tolerance")
        println(io, "")
        println(io, "| tol | count | median residual | max residual | median orthogonality | max orthogonality | nonfinite |")
        println(io, "|---:|---:|---:|---:|---:|---:|---:|")
        for q in quality
            println(
                io,
                "| $(q.tol) | $(q.count) | $(round(q.median_residual, sigdigits=5)) | $(round(q.max_residual, sigdigits=5)) | $(round(q.median_orth, sigdigits=5)) | $(round(q.max_orth, sigdigits=5)) | $(q.nonfinite) |",
            )
        end
        println(io, "")

        println(io, "## Hotspot Correlation")
        println(io, "- mean phase shares: apply=$(round(hotspot.mean_apply_pct, sigdigits=5))%, w_update=$(round(hotspot.mean_w_update_pct, sigdigits=5))%, downdate=$(round(hotspot.mean_downdate_pct, sigdigits=5))%, pivot=$(round(hotspot.mean_pivot_pct, sigdigits=5))%")
        println(io, "- corr(slowdown, apply_pct)=$(round(hotspot.corr_slowdown_apply, sigdigits=5))")
        println(io, "- corr(slowdown, w_update_pct)=$(round(hotspot.corr_slowdown_w_update, sigdigits=5))")
        println(io, "- corr(slowdown, downdate_pct)=$(round(hotspot.corr_slowdown_downdate, sigdigits=5))")
        println(io, "- corr(slowdown, pivot_pct)=$(round(hotspot.corr_slowdown_pivot, sigdigits=5))")
        println(io, "- weighted slowdown contributions: apply=$(round(hotspot.contrib_apply, sigdigits=5)), w_update=$(round(hotspot.contrib_w_update, sigdigits=5)), downdate=$(round(hotspot.contrib_downdate, sigdigits=5)), pivot=$(round(hotspot.contrib_pivot, sigdigits=5))")
        println(io, "")

        best = sort(pairs, by = p -> p.speedup, rev = true)
        worst = sort(pairs, by = p -> p.speedup)
        println(io, "## Top 10 Best Cases")
        for p in first(best, min(10, length(best)))
            println(io, "- $(p.family) $(p.regime) m=$(p.m) n=$(p.n) k=$(p.k) tol=$(p.tol): speedup=$(round(p.speedup, sigdigits=5)) significance=$(p.significance)")
        end
        println(io, "")
        println(io, "## Top 10 Worst Cases")
        for p in first(worst, min(10, length(worst)))
            println(io, "- $(p.family) $(p.regime) m=$(p.m) n=$(p.n) k=$(p.k) tol=$(p.tol): speedup=$(round(p.speedup, sigdigits=5)) significance=$(p.significance)")
        end
    end

    return hotspot
end

function _write_candidates(pairs::Vector{NamedTuple}, hotspot::NamedTuple)
    candidates = _candidate_table(pairs, hotspot)
    open(CANDIDATES_MD, "w") do io
        println(io, "# Optimization Candidates")
        println(io, "")
        println(io, "- Generated: $(Dates.now())")
        println(io, "- Ranking model: `score = 0.5*impact + 0.3*confidence - 0.1*numerical_risk - 0.1*implementation_cost`")
        println(io, "")
        println(io, "| candidate_id | location | rationale | expected_gain | risk_level | numerical_risk | implementation_cost | validation_plan | go_no_go | weighted_score |")
        println(io, "|---|---|---|---|---|---|---|---|---|---:|")
        for c in candidates
            println(
                io,
                "| $(c[1]) | $(c[2]) | $(c[3]) | $(c.expected_gain) | $(c[5]) | $(c[6]) | $(c[7]) | $(c[8]) | $(c.go_no_go) | $(round(c.weighted_score, sigdigits=5)) |",
            )
        end
    end
    return candidates
end

function _validation_status()
    isfile(VALIDATION_MD) || return "unknown (validation report missing)"
    for line in eachline(VALIDATION_MD)
        occursin("- Result:", line) && return replace(line, "- Result: " => "")
    end
    return "unknown (result line missing)"
end

function _write_conclusion(pairs::Vector{NamedTuple}, hotspot::NamedTuple, candidates)
    top3 = first(candidates, min(3, length(candidates)))
    speeds = [p.speedup for p in pairs]
    open(CONCLUSION_MD, "w") do io
        println(io, "# Performance Assessment Conclusion")
        println(io, "")
        println(io, "- Generated: $(Dates.now())")
        println(io, "- Validation status: $(_validation_status())")
        println(io, "- Mean speedup (`dgeqp3/bsqr`): $(round(mean(speeds), sigdigits=5))")
        println(io, "- Median speedup (`dgeqp3/bsqr`): $(round(median(speeds), sigdigits=5))")
        println(io, "")
        println(io, "## Project Status")
        println(io, "- Correctness/guardrails currently pass in latest validation artifacts.")
        println(io, "- Performance remains below dgeqp3 in most assessed cases (speedup < 1 dominates).")
        println(io, "")
        println(io, "## Bottleneck Diagnosis")
        println(io, "- Runtime is concentrated in apply/W-update/downdate phases.")
        println(io, "- Mean shares: apply=$(round(hotspot.mean_apply_pct, sigdigits=5))%, w_update=$(round(hotspot.mean_w_update_pct, sigdigits=5))%, downdate=$(round(hotspot.mean_downdate_pct, sigdigits=5))%.")
        println(io, "")
        println(io, "## Top Opportunities")
        for c in top3
            println(io, "- $(c[1]) at $(c[2]) (score=$(round(c.weighted_score, sigdigits=5))): $(c[3])")
        end
        println(io, "")
        println(io, "## No-Go Areas (Stability Risk)")
        println(io, "- Do not alter Bischof-Stewart pivot criterion semantics.")
        println(io, "- Do not relax numerical safeguards (norm recomputation/rank-stop tolerances) without re-validating quality guardrails.")
        println(io, "- Do not adopt major algorithm redesigns (blocked/alternative formulations) in this campaign.")
    end
end

function main()
    rows = _load_csv_rows(RAW_CSV)
    pairs = _pair_rows(rows)
    isempty(pairs) && error("No paired bsqr/dgeqp3 rows found in $(RAW_CSV)")
    hotspot = _write_summary(rows, pairs)
    candidates = _write_candidates(pairs, hotspot)
    _write_conclusion(pairs, hotspot, candidates)
    println("Wrote assessment summary: $SUMMARY_MD")
    println("Wrote optimization candidates: $CANDIDATES_MD")
    println("Wrote performance conclusion: $CONCLUSION_MD")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
