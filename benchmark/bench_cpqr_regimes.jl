using Dates
using Printf
using Random

using BSPivotQR

include("bench_common.jl")
using .BenchCommon

_parse_families() = parse_symbol_list(
    "BS_REGIME_FAMILIES",
    [:gaussian, :ill_conditioned, :orthonormal_rows],
    [:gaussian, :ill_conditioned, :orthonormal_rows],
)

function _bench_case!(rows, io, rng, family::Symbol, regime::String, fixed_value::Int, var_value::Int, m::Int, n::Int, warmup::Int, samples::Int, norm_recomp_tol::Float64)
    if family === :orthonormal_rows && m > n
        return
    end

    A = make_matrix(family, m, n, rng)
    kfull = min(m, n)

    f_bs() = run_bsqr_fair(A, kfull, norm_recomp_tol)
    tmin, tmed, tci_low, tci_high, alloc = bench_trial_ci(f_bs; warmup = warmup, samples = samples)
    Fbs = f_bs()
    rb, qb = residual_bs(A, Fbs)
    push!(rows, (
        family = family,
        regime = regime,
        fixed_value = fixed_value,
        var_value = var_value,
        m = m,
        n = n,
        k = kfull,
        method = "bsqr_full",
        tmin = tmin,
        tmed = tmed,
        tci_low = tci_low,
        tci_high = tci_high,
        alloc = alloc,
        resid = rb,
        orth = qb,
    ))
    println(io, "$(family),$(regime),$(fixed_value),$(var_value),$(m),$(n),$(kfull),bsqr_full,$(tmin),$(tmed),$(tci_low),$(tci_high),$(alloc),$(rb),$(qb)")

    f_qr() = run_qr_fair(A)
    tmin, tmed, tci_low, tci_high, alloc = bench_trial_ci(f_qr; warmup = warmup, samples = samples)
    Fq = f_qr()
    rq, qq = residual_qr(A, Fq)
    push!(rows, (
        family = family,
        regime = regime,
        fixed_value = fixed_value,
        var_value = var_value,
        m = m,
        n = n,
        k = kfull,
        method = "dgeqp3",
        tmin = tmin,
        tmed = tmed,
        tci_low = tci_low,
        tci_high = tci_high,
        alloc = alloc,
        resid = rq,
        orth = qq,
    ))
    println(io, "$(family),$(regime),$(fixed_value),$(var_value),$(m),$(n),$(kfull),dgeqp3,$(tmin),$(tmed),$(tci_low),$(tci_high),$(alloc),$(rq),$(qq)")
end

function run_regime_sweep()
    check_backend()
    configure_blas_threads()
    norm_recomp_tol = parse_env_float("BS_NORM_RECOMP_TOL", DEFAULT_NORM_RECOMP_TOL)

    quick = get(ENV, "BS_REGIME_QUICK", "0") == "1"
    families = _parse_families()

    default_m_sweep = quick ? [32, 64, 128, 256, 512, 1024] : [32, 64, 96, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096]
    default_n_sweep = quick ? [64, 128, 256, 512, 1024, 2048] : [64, 128, 192, 256, 384, 512, 768, 1024, 1536, 2048, 3072, 4096, 6144, 8192]
    default_fixed_ns = quick ? [128, 256] : [128, 256, 512]
    default_fixed_ms = quick ? [32, 64] : [32, 64, 128]

    m_sweep = parse_int_list("BS_REGIME_MS", default_m_sweep)
    n_sweep = parse_int_list("BS_REGIME_NS", default_n_sweep)
    fixed_ns = parse_int_list("BS_REGIME_FIXED_NS", default_fixed_ns)
    fixed_ms = parse_int_list("BS_REGIME_FIXED_MS", default_fixed_ms)

    warmup = parse_env_int("BS_REGIME_WARMUP", quick ? 1 : 2; minval = 0)
    samples = parse_env_int("BS_REGIME_SAMPLES", quick ? 12 : 24)

    mkpath(joinpath(@__DIR__, "results"))
    csv_path = joinpath(@__DIR__, "results", "regime_timings.csv")
    md_path = joinpath(@__DIR__, "results", "regime_summary.md")

    io = open(csv_path, "w")
    println(io, "family,regime,fixed_value,var_value,m,n,k,method,tmin_s,tmed_s,tci_low_s,tci_high_s,alloc_bytes,residual,orthogonality")

    rng = MersenneTwister(20260224)
    rows = Vector{NamedTuple}()

    for family in families
        if family !== :orthonormal_rows
            for nfix in fixed_ns
                for m in m_sweep
                    _bench_case!(rows, io, rng, family, "fixed_n_vary_m", nfix, m, m, nfix, warmup, samples, norm_recomp_tol)
                end
            end
        end

        for mfix in fixed_ms
            for n in n_sweep
                _bench_case!(rows, io, rng, family, "fixed_m_vary_n", mfix, n, mfix, n, warmup, samples, norm_recomp_tol)
            end
        end
    end
    close(io)

    open(md_path, "w") do md
        println(md, "# CPQR Regime Sweep Summary")
        println(md, "")
        println(md, "Generated: $(Dates.now())")
        println(md, "")
        println(md, "BLAS: ", backend_string())
        println(md, "")
        println(md, "| family | regime | fixed | varying | m | n | method | median (s) | 95% CI (s) | speedup vs dgeqp3 | residual | orthogonality |")
        println(md, "|---|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|")

        groups = unique((r.family, r.regime, r.fixed_value, r.var_value) for r in rows)
        for (family, regime, fixed_value, var_value) in sort(collect(groups), by = x -> (string(x[1]), x[2], x[3], x[4]))
            keyrows = filter(r -> r.family == family && r.regime == regime && r.fixed_value == fixed_value && r.var_value == var_value, rows)
            drows = filter(r -> r.method == "dgeqp3", keyrows)
            isempty(drows) && continue
            d = only(drows)
            for r in sort(keyrows, by = x -> x.method)
                speed = d.tmed / r.tmed
                ci = "[$(round(r.tci_low, sigdigits=5)), $(round(r.tci_high, sigdigits=5))]"
                println(md, "| $(r.family) | $(r.regime) | $(r.fixed_value) | $(r.var_value) | $(r.m) | $(r.n) | $(r.method) | $(round(r.tmed, sigdigits=5)) | $(ci) | $(round(speed, sigdigits=5)) | $(round(r.resid, sigdigits=5)) | $(round(r.orth, sigdigits=5)) |")
            end
        end
    end

    println("Wrote regime sweep results to: $csv_path")
    println("Wrote regime sweep summary to: $md_path")
end

run_regime_sweep()
