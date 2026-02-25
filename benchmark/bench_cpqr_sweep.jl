using Dates
using Printf
using Random

using BSPivotQR

include("bench_common.jl")
using .BenchCommon

_parse_families() = parse_symbol_list(
    "BS_SWEEP_FAMILIES",
    [:gaussian, :ill_conditioned],
    [:gaussian, :ill_conditioned, :orthonormal_rows],
)

function run_sweep()
    check_backend()
    configure_blas_threads()
    norm_recomp_tol = parse_env_float("BS_NORM_RECOMP_TOL", DEFAULT_NORM_RECOMP_TOL)

    quick = get(ENV, "BS_SWEEP_QUICK", "0") == "1"

    default_ns = quick ? [64, 128, 256, 512, 1024] : [64, 128, 256, 384, 512, 768, 1024, 1536, 2048]
    default_aspects = quick ? [0.5, 1.0, 2.0] : [0.25, 0.5, 1.0, 2.0, 4.0]

    ns = parse_int_list("BS_SWEEP_NS", default_ns)
    aspects = parse_float_list("BS_SWEEP_ASPECTS", default_aspects)
    families = _parse_families()

    warmup = parse_env_int("BS_SWEEP_WARMUP", quick ? 1 : 2; minval = 0)
    samples = parse_env_int("BS_SWEEP_SAMPLES", quick ? 12 : 24)

    mkpath(joinpath(@__DIR__, "results"))
    csv_path = joinpath(@__DIR__, "results", "sweep_timings.csv")
    md_path = joinpath(@__DIR__, "results", "sweep_summary.md")

    io = open(csv_path, "w")
    println(io, "family,aspect,m,n,k,method,tmin_s,tmed_s,tci_low_s,tci_high_s,alloc_bytes,residual,orthogonality")

    rng = MersenneTwister(20260223)
    rows = Vector{NamedTuple}()

    for family in families
        for aspect in aspects
            for n in ns
                m = max(16, round(Int, aspect * n))
                if family === :orthonormal_rows && m > n
                    continue
                end

                A = make_matrix(family, m, n, rng)
                kfull = min(m, n)
                aspect_s = @sprintf("%.3f", aspect)

                bs_row, dg_row = bench_pair_ci(A, kfull, norm_recomp_tol; warmup = warmup, samples = samples)
                for row in (bs_row, dg_row)
                    push!(rows, (
                        family = family,
                        aspect = aspect,
                        m = m,
                        n = n,
                        k = kfull,
                        method = row.method,
                        tmin = row.tmin,
                        tmed = row.tmed,
                        tci_low = row.tci_low,
                        tci_high = row.tci_high,
                        alloc = row.alloc,
                        resid = row.resid,
                        orth = row.orth,
                    ))
                    println(io, "$(family),$aspect_s,$m,$n,$kfull,$(row.method),$(row.tmin),$(row.tmed),$(row.tci_low),$(row.tci_high),$(row.alloc),$(row.resid),$(row.orth)")
                end
            end
        end
    end
    close(io)

    open(md_path, "w") do md
        println(md, "# CPQR Systematic Sweep Summary")
        println(md, "")
        println(md, "Generated: $(Dates.now())")
        println(md, "")
        println(md, "BLAS: ", backend_string())
        println(md, "")
        println(md, "| family | aspect m/n | m | n | k | method | median (s) | 95% CI (s) | min (s) | speedup vs dgeqp3 | residual | orthogonality |")
        println(md, "|---|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|---:|")

        for group in grouped_rows_with_baseline(rows, r -> (r.family, r.aspect, r.m, r.n); sortby = x -> (string(x[1]), x[2], x[4], x[3]))
            d = group.baseline
            for r in group.rows
                speed = d.tmed / r.tmed
                ci = "[$(round(r.tci_low, sigdigits=5)), $(round(r.tci_high, sigdigits=5))]"
                println(md, "| $(r.family) | $(@sprintf("%.3f", r.aspect)) | $(r.m) | $(r.n) | $(r.k) | $(r.method) | $(round(r.tmed, sigdigits=5)) | $(ci) | $(round(r.tmin, sigdigits=5)) | $(round(speed, sigdigits=5)) | $(round(r.resid, sigdigits=5)) | $(round(r.orth, sigdigits=5)) |")
            end
        end
    end

    println("Wrote sweep results to: $csv_path")
    println("Wrote sweep summary to: $md_path")
end

run_sweep()
