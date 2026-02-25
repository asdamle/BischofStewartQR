using Dates
using LinearAlgebra
using Random
using Statistics

using BSPivotQR

include("bench_common.jl")
using .BenchCommon

_parse_families() = parse_symbol_list(
    "BS_PROFILE_FAMILIES",
    [:gaussian, :ill_conditioned, :orthonormal_rows],
    [:gaussian, :ill_conditioned, :orthonormal_rows],
)
function _parse_env_float_list(envkey::String, default::Vector{Float64})
    vals = parse_float_list(envkey, default)
    for v in vals
        (0.0 <= v <= 1.0) || error("$envkey entries must satisfy 0 <= value <= 1")
    end
    return vals
end

function _profile_once(A::Matrix{Float64}, tol::Float64)
    m, n = size(A)
    k = min(m, n)
    ws = BSPivotQR.BSWorkspace(m, n, k)
    tau = zeros(Float64, k)
    jpvt = collect(1:n)
    stats = BSPivotQR.BSKernelStats()
    Awork = copy(A)

    t0 = time_ns()
    ksteps = BSPivotQR._bsqr_kernel!(
        Awork,
        tau,
        jpvt,
        ws,
        k;
        rank_stop = false,
        norm_recomp_tol = tol,
        kernel_stats = stats,
    )
    total_ns = time_ns() - t0

    F = BSPivotQR.BSQRPivoted(Awork, tau, jpvt, ksteps, nothing, nothing)
    Q = BSPivotQR._explicit_q(F)
    T = BSPivotQR._packed_to_qt(F)
    resid = norm(A[:, BSPivotQR.perm(F)] - Q * T) / max(norm(A), eps(Float64))
    orth = norm(I - Q' * Q)
    return total_ns, stats, resid, orth, ksteps
end

function _profile_case(A::Matrix{Float64}, tol::Float64, reps::Int)
    _profile_once(A, tol) # warmup

    totals = Float64[]
    pivot_ns = Int[]
    hh_ns = Int[]
    apply_ns = Int[]
    wup_ns = Int[]
    down_ns = Int[]
    rec_counts = Int[]
    residuals = Float64[]
    orths = Float64[]
    ksteps_last = 0

    for _ in 1:reps
        total_ns, stats, resid, orth, ksteps = _profile_once(A, tol)
        push!(totals, total_ns * 1e-9)
        push!(pivot_ns, stats.pivot_select_ns)
        push!(hh_ns, stats.householder_ns)
        push!(apply_ns, stats.apply_reflector_ns)
        push!(wup_ns, stats.w_update_ns)
        push!(down_ns, stats.norm_downdate_ns)
        push!(rec_counts, stats.recompute_count)
        push!(residuals, resid)
        push!(orths, orth)
        ksteps_last = ksteps
    end

    return (
        ksteps = ksteps_last,
        total_med_s = median(totals),
        total_min_s = minimum(totals),
        pivot_s = median(pivot_ns) * 1e-9,
        householder_s = median(hh_ns) * 1e-9,
        apply_s = median(apply_ns) * 1e-9,
        w_update_s = median(wup_ns) * 1e-9,
        norm_downdate_s = median(down_ns) * 1e-9,
        recompute_count = Int(round(median(rec_counts))),
        residual = median(residuals),
        orthogonality = median(orths),
    )
end

function run_profile()
    check_backend()
    configure_blas_threads()

    quick = get(ENV, "BS_PROFILE_QUICK", "1") == "1"
    reps = parse_env_int("BS_PROFILE_REPS", quick ? 3 : 5)

    short_m = parse_env_int("BS_PROFILE_SHORT_M", 64)
    short_n = parse_env_int("BS_PROFILE_SHORT_N", quick ? 2048 : 8192)
    square_n = parse_env_int("BS_PROFILE_SQUARE_N", quick ? 256 : 1024)
    tols = _parse_env_float_list("BS_PROFILE_TOLS", [sqrt(eps(Float64)), 0.0])
    families = _parse_families()

    cases = [
        ("short_wide", short_m, short_n),
        ("square", square_n, square_n),
    ]

    mkpath(joinpath(@__DIR__, "results"))
    csv_path = joinpath(@__DIR__, "results", "profile_breakdown.csv")
    md_path = joinpath(@__DIR__, "results", "profile_breakdown_summary.md")

    io = open(csv_path, "w")
    println(io, "family,case,m,n,tol,ksteps,total_med_s,total_min_s,pivot_s,householder_s,apply_reflector_s,w_update_s,norm_downdate_s,recompute_count,residual,orthogonality")

    rng = MersenneTwister(20260224)
    rows = NamedTuple[]
    for family in families
        for (case_name, m, n) in cases
            if family === :orthonormal_rows && m > n
                continue
            end
            if family === :orthonormal_rows && case_name != "short_wide"
                continue
            end
            A = make_matrix(family, m, n, rng)
            for tol in tols
                metrics = _profile_case(A, tol, reps)
                row = (
                    family = family,
                    case = case_name,
                    m = m,
                    n = n,
                    tol = tol,
                    ksteps = metrics.ksteps,
                    total_med_s = metrics.total_med_s,
                    total_min_s = metrics.total_min_s,
                    pivot_s = metrics.pivot_s,
                    householder_s = metrics.householder_s,
                    apply_s = metrics.apply_s,
                    w_update_s = metrics.w_update_s,
                    norm_downdate_s = metrics.norm_downdate_s,
                    recompute_count = metrics.recompute_count,
                    residual = metrics.residual,
                    orthogonality = metrics.orthogonality,
                )
                push!(rows, row)
                println(
                    io,
                    "$(row.family),$(row.case),$(row.m),$(row.n),$(row.tol),$(row.ksteps),$(row.total_med_s),$(row.total_min_s),$(row.pivot_s),$(row.householder_s),$(row.apply_s),$(row.w_update_s),$(row.norm_downdate_s),$(row.recompute_count),$(row.residual),$(row.orthogonality)",
                )
            end
        end
    end
    close(io)

    open(md_path, "w") do md
        println(md, "# BSQR Kernel Breakdown Profile")
        println(md, "")
        println(md, "Generated: $(Dates.now())")
        println(md, "")
        println(md, "BLAS: ", backend_string())
        println(md, "")
        println(md, "| family | case | m | n | tol | total median (s) | pivot % | householder % | apply % | W-update % | downdate % | recomputes | residual | orthogonality |")
        println(md, "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|")
        for r in sort(rows, by = x -> (string(x.family), x.case, x.m, x.n, x.tol))
            t = max(r.total_med_s, eps(Float64))
            p_pivot = 100 * r.pivot_s / t
            p_hh = 100 * r.householder_s / t
            p_apply = 100 * r.apply_s / t
            p_w = 100 * r.w_update_s / t
            p_down = 100 * r.norm_downdate_s / t
            println(
                md,
                "| $(r.family) | $(r.case) | $(r.m) | $(r.n) | $(r.tol) | $(round(r.total_med_s, sigdigits=5)) | $(round(p_pivot, sigdigits=4)) | $(round(p_hh, sigdigits=4)) | $(round(p_apply, sigdigits=4)) | $(round(p_w, sigdigits=4)) | $(round(p_down, sigdigits=4)) | $(r.recompute_count) | $(round(r.residual, sigdigits=4)) | $(round(r.orthogonality, sigdigits=4)) |",
            )
        end
    end

    println("Wrote profile CSV to: $csv_path")
    println("Wrote profile summary to: $md_path")
end

run_profile()
