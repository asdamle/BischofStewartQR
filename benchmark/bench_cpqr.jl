using Dates
using Random

using BSPivotQR

include("bench_common.jl")
using .BenchCommon

function run_benchmarks()
    check_backend()
    configure_blas_threads()
    norm_recomp_tol = parse_env_float("BS_NORM_RECOMP_TOL", DEFAULT_NORM_RECOMP_TOL)
    quick = get(ENV, "BS_QUICK", "0") == "1"

    mkpath(joinpath(@__DIR__, "results"))
    csv_path = joinpath(@__DIR__, "results", "timings.csv")
    md_path = joinpath(@__DIR__, "results", "summary.md")

    io = open(csv_path, "w")
    println(io, "family,m,n,k,method,tmin_s,tmed_s,alloc_bytes,residual,orthogonality")

    rng = MersenneTwister(20260223)
    sizes = quick ?
        [(64, 64), (128, 64), (64, 128)] :
        [
            (64, 64), (128, 64), (64, 128),
            (512, 512), (2000, 256), (256, 2000),
            (2048, 2048), (8000, 512),
        ]
    orth_sizes = quick ?
        [(64, 128), (64, 256), (64, 512)] :
        [(64, 128), (64, 256), (64, 512), (64, 1024), (64, 2048), (64, 4096)]
    families = (:gaussian, :ill_conditioned, :orthonormal_rows)
    warmup = quick ? 1 : 2
    samples = quick ? 3 : 8

    rows = Vector{NamedTuple}()

    for family in families
        family_sizes = family === :orthonormal_rows ? orth_sizes : sizes
        for (m, n) in family_sizes
            if family === :orthonormal_rows && m > n
                continue
            end
            A = make_matrix(family, m, n, rng)
            kfull = min(size(A)...)

            for row in bench_pair_basic(A, kfull, norm_recomp_tol; warmup = warmup, samples = samples)
                push!(rows, (
                    family = family,
                    m = size(A, 1),
                    n = size(A, 2),
                    k = kfull,
                    method = row.method,
                    tmin = row.tmin,
                    tmed = row.tmed,
                    alloc = row.alloc,
                    resid = row.resid,
                    orth = row.orth,
                ))
                println(io, "$(family),$(size(A,1)),$(size(A,2)),$kfull,$(row.method),$(row.tmin),$(row.tmed),$(row.alloc),$(row.resid),$(row.orth)")
            end
        end
    end
    close(io)

    open(md_path, "w") do md
        println(md, "# CPQR Benchmark Summary")
        println(md, "")
        println(md, "Generated: $(Dates.now())")
        println(md, "")
        println(md, "BLAS: ", backend_string())
        println(md, "")
        println(md, "| family | m | n | k | method | median (s) | min (s) | speedup vs dgeqp3 | residual | orthogonality |")
        println(md, "|---|---:|---:|---:|---|---:|---:|---:|---:|---:|")

        for group in grouped_rows_with_baseline(rows, r -> (r.family, r.m, r.n); sortby = x -> (string(x[1]), x[2], x[3]))
            d = group.baseline
            for r in group.rows
                speed = d.tmed / r.tmed
                println(md, "| $(r.family) | $(r.m) | $(r.n) | $(r.k) | $(r.method) | $(round(r.tmed, sigdigits=4)) | $(round(r.tmin, sigdigits=4)) | $(round(speed, sigdigits=4)) | $(round(r.resid, sigdigits=4)) | $(round(r.orth, sigdigits=4)) |")
            end
        end
    end

    println("Wrote benchmark results to: $csv_path")
    println("Wrote benchmark summary to: $md_path")
end

run_benchmarks()
