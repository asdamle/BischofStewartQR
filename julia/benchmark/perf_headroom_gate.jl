#!/usr/bin/env julia

using LinearAlgebra
using Printf
using Random
using Statistics

using BSPivotQR

include("bench_common.jl")
using .BenchCommon: make_matrix

const DEFAULT_VARIANTS = ["baseline", "fastpath_off", "aspect6"]
const DEFAULT_FAMILIES = [:gaussian, :ill_conditioned, :orthonormal_rows]
const DEFAULT_CASES = [
    (regime = "square", m = 256, n = 256, aspect = 1.0),
    (regime = "square", m = 512, n = 512, aspect = 1.0),
    (regime = "short_wide", m = 128, n = 1024, aspect = 8.0),
    (regime = "short_wide", m = 64, n = 640, aspect = 10.0),
]

function _parse_symbol_list(envkey::String, default::Vector{Symbol})
    s = strip(get(ENV, envkey, ""))
    isempty(s) && return default
    out = Symbol[]
    for tok in split(lowercase(s), ',')
        t = Symbol(strip(tok))
        isempty(String(t)) && continue
        push!(out, t)
    end
    isempty(out) && return default
    return unique(out)
end

function _parse_variant_list(envkey::String, default::Vector{String})
    s = strip(get(ENV, envkey, ""))
    isempty(s) && return default
    out = String[]
    for tok in split(lowercase(s), ',')
        t = strip(tok)
        isempty(t) && continue
        push!(out, t)
    end
    isempty(out) && return default
    return unique(out)
end

function _save_env(keys::Vector{String})
    old = Dict{String,Union{Nothing,String}}()
    for key in keys
        old[key] = haskey(ENV, key) ? ENV[key] : nothing
    end
    return old
end

function _restore_env!(old::Dict{String,Union{Nothing,String}})
    for (key, val) in old
        if val === nothing
            haskey(ENV, key) && delete!(ENV, key)
        else
            ENV[key] = val
        end
    end
    return nothing
end

function _configure_variant!(variant::String)
    keys = [
        "BS_SHORT_WIDE_FASTPATH",
        "BS_SHORT_WIDE_FASTPATH_ASPECT",
        "BS_SHORT_WIDE_FASTPATH_MMAX",
        "BS_SHORT_WIDE_FASTPATH_NMIN",
    ]
    old = _save_env(keys)

    if variant == "baseline"
        ENV["BS_SHORT_WIDE_FASTPATH"] = "1"
        haskey(ENV, "BS_SHORT_WIDE_FASTPATH_ASPECT") && delete!(ENV, "BS_SHORT_WIDE_FASTPATH_ASPECT")
        haskey(ENV, "BS_SHORT_WIDE_FASTPATH_MMAX") && delete!(ENV, "BS_SHORT_WIDE_FASTPATH_MMAX")
        haskey(ENV, "BS_SHORT_WIDE_FASTPATH_NMIN") && delete!(ENV, "BS_SHORT_WIDE_FASTPATH_NMIN")
    elseif variant == "fastpath_off"
        ENV["BS_SHORT_WIDE_FASTPATH"] = "0"
        haskey(ENV, "BS_SHORT_WIDE_FASTPATH_ASPECT") && delete!(ENV, "BS_SHORT_WIDE_FASTPATH_ASPECT")
        haskey(ENV, "BS_SHORT_WIDE_FASTPATH_MMAX") && delete!(ENV, "BS_SHORT_WIDE_FASTPATH_MMAX")
        haskey(ENV, "BS_SHORT_WIDE_FASTPATH_NMIN") && delete!(ENV, "BS_SHORT_WIDE_FASTPATH_NMIN")
    elseif variant == "aspect6"
        ENV["BS_SHORT_WIDE_FASTPATH"] = "1"
        ENV["BS_SHORT_WIDE_FASTPATH_ASPECT"] = "6"
        ENV["BS_SHORT_WIDE_FASTPATH_MMAX"] = "256"
        ENV["BS_SHORT_WIDE_FASTPATH_NMIN"] = "512"
    else
        error("Unsupported variant: $variant")
    end
    return old
end

function _factor_quality(A::Matrix{Float64}, F::BSQRPivoted)
    Q = BSPivotQR._explicit_q(F)
    T = BSPivotQR._packed_to_qt(F)
    Aperm = A[:, perm(F)]
    resid = norm(Aperm - Q * T) / max(norm(A), eps(Float64))
    orth = norm(I - Q' * Q)
    return resid, orth
end

function _run_kernel_case(
    A::Matrix{Float64};
    warmup::Int,
    samples::Int,
)
    m, n = size(A)
    k = min(m, n)

    for _ in 1:warmup
        Aw = copy(A)
        tau = zeros(Float64, k)
        jpvt = collect(1:n)
        ws = BSPivotQR.BSWorkspace(m, n, k)
        BSPivotQR._bsqr_kernel!(Aw, tau, jpvt, ws, k; rank_stop = false)
    end

    tmeds = Float64[]
    stage_shares = NamedTuple[]
    recomputes = Int[]
    for _ in 1:samples
        Aw = copy(A)
        tau = zeros(Float64, k)
        jpvt = collect(1:n)
        ws = BSPivotQR.BSWorkspace(m, n, k)
        stats = BSPivotQR.BSKernelStats()
        rc = Ref(0)

        t = @elapsed BSPivotQR._bsqr_kernel!(
            Aw,
            tau,
            jpvt,
            ws,
            k;
            rank_stop = false,
            kernel_stats = stats,
            norm_recomp_count = rc,
        )
        push!(tmeds, t)
        push!(recomputes, rc[])

        total_ns = stats.pivot_select_ns + stats.householder_ns + stats.apply_reflector_ns + stats.w_update_ns + stats.norm_downdate_ns
        denom = max(total_ns, 1)
        push!(stage_shares, (
            pivot = stats.pivot_select_ns / denom,
            householder = stats.householder_ns / denom,
            apply = stats.apply_reflector_ns / denom,
            w_update = stats.w_update_ns / denom,
            norm_refresh = stats.norm_downdate_ns / denom,
        ))
    end

    F = bsqr(copy(A); check = false, rank_stop = false)
    resid, orth = _factor_quality(A, F)

    return (
        tmed = median(tmeds),
        recompute_median = median(recomputes),
        pivot_share = median(map(x -> x.pivot, stage_shares)),
        householder_share = median(map(x -> x.householder, stage_shares)),
        apply_share = median(map(x -> x.apply, stage_shares)),
        w_update_share = median(map(x -> x.w_update, stage_shares)),
        norm_refresh_share = median(map(x -> x.norm_refresh, stage_shares)),
        residual = resid,
        orthogonality = orth,
    )
end

function _geomean(vals::Vector{Float64})
    good = filter(v -> isfinite(v) && v > 0.0, vals)
    isempty(good) && return NaN
    return exp(mean(log.(good)))
end

function run_headroom_gate()
    BLAS.set_num_threads(parse(Int, get(ENV, "BS_HEADROOM_BLAS_THREADS", "1")))
    warmup = parse(Int, get(ENV, "BS_HEADROOM_WARMUP", "1"))
    samples = parse(Int, get(ENV, "BS_HEADROOM_SAMPLES", "8"))
    seed = parse(Int, get(ENV, "BS_HEADROOM_SEED", "20260312"))
    gate = parse(Float64, get(ENV, "BS_HEADROOM_GATE", "1.15"))
    enforce = strip(get(ENV, "BS_HEADROOM_ENFORCE", "0")) == "1"

    variants = _parse_variant_list("BS_HEADROOM_VARIANTS", DEFAULT_VARIANTS)
    families = _parse_symbol_list("BS_HEADROOM_FAMILIES", DEFAULT_FAMILIES)
    cases = DEFAULT_CASES

    println("# Perf Headroom Gate")
    println("BLAS threads: $(BLAS.get_num_threads()) | warmup=$warmup samples=$samples seed=$seed gate=$gate")
    println("variants: ", join(variants, ", "))
    println("families: ", join(string.(families), ", "))
    println("")

    results = Dict{Tuple{String,Symbol,Int,Int},NamedTuple}()
    for family in families
        for c in cases
            rng = MersenneTwister(seed + 10_000 * c.m + c.n)
            A = make_matrix(family, c.m, c.n, rng)
            for variant in variants
                old_env = _configure_variant!(variant)
                try
                    results[(variant, family, c.m, c.n)] = _run_kernel_case(
                        A;
                        warmup = warmup,
                        samples = samples,
                    )
                finally
                    _restore_env!(old_env)
                end
            end
        end
    end

    baseline = "baseline"
    haskey(results, (baseline, first(families), first(cases).m, first(cases).n)) ||
        error("baseline variant must be included")

    speedups = Dict{String,Vector{Float64}}(v => Float64[] for v in variants)
    for family in families
        for c in cases
            tbase = results[(baseline, family, c.m, c.n)].tmed
            for variant in variants
                tvar = results[(variant, family, c.m, c.n)].tmed
                push!(speedups[variant], tbase / tvar)
            end
        end
    end

    println("| variant | geomean speedup vs baseline | median residual | median orthogonality |")
    println("|---|---:|---:|---:|")
    for variant in variants
        vals = speedups[variant]
        g = _geomean(vals)
        residuals = [results[(variant, family, c.m, c.n)].residual for family in families for c in cases]
        orths = [results[(variant, family, c.m, c.n)].orthogonality for family in families for c in cases]
        @printf("| %s | %.5f | %.3e | %.3e |\n", variant, g, median(residuals), median(orths))
    end

    println("")
    println("Stage-share medians (baseline):")
    println("| family | m | n | pivot | householder | apply | w_update | norm_refresh | recompute_median |")
    println("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for family in families
        for c in cases
            r = results[(baseline, family, c.m, c.n)]
            @printf(
                "| %s | %d | %d | %.3f | %.3f | %.3f | %.3f | %.3f | %.1f |\n",
                string(family),
                c.m,
                c.n,
                r.pivot_share,
                r.householder_share,
                r.apply_share,
                r.w_update_share,
                r.norm_refresh_share,
                r.recompute_median,
            )
        end
    end

    best_variant = baseline
    best_speedup = 1.0
    for variant in variants
        g = _geomean(speedups[variant])
        if isfinite(g) && g > best_speedup
            best_speedup = g
            best_variant = variant
        end
    end

    println("")
    @printf("Best variant: %s (geomean speedup %.5f)\n", best_variant, best_speedup)
    if best_speedup >= gate
        println("Gate PASSED: substantial improvement found.")
    else
        println("Gate NOT met: no substantial improvement found; prioritize readability + benchmark quality upgrades.")
        enforce && error("BS_HEADROOM_ENFORCE=1 and no variant reached the speedup gate")
    end
end

run_headroom_gate()
