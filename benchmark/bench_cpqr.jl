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

            f_bs_full() = run_bsqr_fair(A, kfull, norm_recomp_tol)
            tmin, tmed, alloc = bench_trial_basic(f_bs_full; warmup = warmup, samples = samples)
            Fbs = f_bs_full()
            rb, qb = residual_bs(A, Fbs)
            push!(rows, (family = family, m = size(A, 1), n = size(A, 2), k = kfull, method = "bsqr_full", tmin = tmin, tmed = tmed, alloc = alloc, resid = rb, orth = qb))
            println(io, "$(family),$(size(A,1)),$(size(A,2)),$kfull,bsqr_full,$tmin,$tmed,$alloc,$rb,$qb")

            f_qr() = run_qr_fair(A)
            tmin, tmed, alloc = bench_trial_basic(f_qr; warmup = warmup, samples = samples)
            Fq = f_qr()
            rq, qq = residual_qr(A, Fq)
            push!(rows, (family = family, m = size(A, 1), n = size(A, 2), k = kfull, method = "dgeqp3", tmin = tmin, tmed = tmed, alloc = alloc, resid = rq, orth = qq))
            println(io, "$(family),$(size(A,1)),$(size(A,2)),$kfull,dgeqp3,$tmin,$tmed,$alloc,$rq,$qq")
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

        groups = unique((r.family, r.m, r.n) for r in rows)
        for (family, m, n) in sort(collect(groups), by = x -> (string(x[1]), x[2], x[3]))
            keyrows = filter(r -> r.family == family && r.m == m && r.n == n, rows)
            dgeqp3_rows = filter(r -> r.method == "dgeqp3", keyrows)
            isempty(dgeqp3_rows) && continue
            d = only(dgeqp3_rows)
            for r in sort(keyrows, by = x -> x.method)
                speed = d.tmed / r.tmed
                println(md, "| $(r.family) | $(r.m) | $(r.n) | $(r.k) | $(r.method) | $(round(r.tmed, sigdigits=4)) | $(round(r.tmin, sigdigits=4)) | $(round(speed, sigdigits=4)) | $(round(r.resid, sigdigits=4)) | $(round(r.orth, sigdigits=4)) |")
            end
        end
    end

    println("Wrote benchmark results to: $csv_path")
    println("Wrote benchmark summary to: $md_path")
end

run_benchmarks()
