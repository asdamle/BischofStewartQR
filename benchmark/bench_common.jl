module BenchCommon

if Sys.isapple() && get(ENV, "BS_USE_ACCELERATE", "1") == "1"
    try
        @eval using AppleAccelerate
    catch err
        @warn "AppleAccelerate not available; continuing with current BLAS backend" err
    end
end

using BenchmarkTools
using LinearAlgebra
using Random
using Statistics

using BSPivotQR

include("matrix_generators.jl")
using .MatrixGenerators

const DEFAULT_NORM_RECOMP_TOL = sqrt(eps(Float64))
const BSQR_METHOD_LABEL = "bsqr_full"
const BSQR_LAZY_BLAS_METHOD_LABEL = "bsqr_lazy_blas"
const DGEQP3_METHOD_LABEL = "dgeqp3"
const DGEQP3_TRSM_METHOD_LABEL = "dgeqp3_trsm"

function backend_string()
    return sprint(show, BLAS.get_config())
end

function accelerate_active()
    return occursin("Accelerate", backend_string())
end

function check_backend()
    cfg = backend_string()
    println("BLAS config: $cfg")
    if Sys.isapple() && get(ENV, "BS_REQUIRE_ACCELERATE", "0") == "1"
        accelerate_active() || error("BS_REQUIRE_ACCELERATE=1 but Accelerate backend is not active")
    end
end

function configure_blas_threads()
    s = strip(get(ENV, "BS_BLAS_THREADS", ""))
    isempty(s) && return nothing
    nth = parse(Int, s)
    nth >= 1 || error("BS_BLAS_THREADS must be >= 1")
    BLAS.set_num_threads(nth)
    println("Configured BLAS threads: ", BLAS.get_num_threads())
    return nth
end

function parse_env_int(envkey::String, default::Int; minval::Int = 1)
    s = strip(get(ENV, envkey, ""))
    isempty(s) && return default
    v = parse(Int, s)
    v >= minval || error("$envkey must be >= $minval")
    return v
end

function parse_env_float(
    envkey::String,
    default::Float64;
    minval::Float64 = 0.0,
    maxval::Float64 = 1.0,
)
    s = strip(get(ENV, envkey, ""))
    isempty(s) && return default
    v = parse(Float64, s)
    (minval <= v <= maxval) || error("$envkey must satisfy $minval <= value <= $maxval")
    return v
end

function parse_int_list(envkey::String, default::Vector{Int})
    s = strip(get(ENV, envkey, ""))
    isempty(s) && return default
    out = Int[]
    for tok in split(s, ',')
        t = strip(tok)
        isempty(t) && continue
        push!(out, parse(Int, t))
    end
    isempty(out) && return default
    return sort(unique(out))
end

function parse_float_list(envkey::String, default::Vector{Float64})
    s = strip(get(ENV, envkey, ""))
    isempty(s) && return default
    out = Float64[]
    for tok in split(s, ',')
        t = strip(tok)
        isempty(t) && continue
        push!(out, parse(Float64, t))
    end
    isempty(out) && return default
    return sort(unique(out))
end

function parse_symbol_list(
    envkey::String,
    default::Vector{Symbol},
    allowed::Vector{Symbol},
)
    s = lowercase(strip(get(ENV, envkey, join(string.(default), ","))))
    vals = Symbol[]
    for tok in split(s, ',')
        t = Symbol(strip(tok))
        t in allowed || error("Unsupported value in $envkey: $(String(t))")
        push!(vals, t)
    end
    isempty(vals) && return default
    return unique(vals)
end

function lazy_benchmark_enabled()
    return get(ENV, "BS_BENCH_INCLUDE_LAZY", "0") == "1"
end

function lazy_batch_kwargs(n::Int)
    batch_size_s = strip(get(ENV, "BS_LAZY_BATCH_SIZE", ""))
    batch_size = isempty(batch_size_s) ? nothing : parse(Int, batch_size_s)
    batch_fraction = parse_env_float("BS_LAZY_BATCH_FRACTION", 0.125; minval = 0.0, maxval = 1.0)
    batch_min = parse_env_int("BS_LAZY_BATCH_MIN", min(8, max(n, 1)))
    batch_max_s = strip(get(ENV, "BS_LAZY_BATCH_MAX", ""))
    batch_max = isempty(batch_max_s) ? nothing : parse(Int, batch_max_s)
    return (
        lazy_batch_size = batch_size,
        lazy_batch_fraction = batch_fraction,
        lazy_batch_min = batch_min,
        lazy_batch_max = batch_max,
    )
end

function bench_trial_basic(f; warmup::Int = 2, samples::Int = 8)
    for _ in 1:warmup
        f()
    end
    trial = @benchmark $f() evals = 1 samples = samples
    tmin = minimum(trial).time * 1e-9
    tmed = median(trial).time * 1e-9
    alloc = minimum(trial).memory
    return tmin, tmed, alloc
end

function bench_trial_ci(f; warmup::Int = 2, samples::Int = 24)
    for _ in 1:warmup
        f()
    end
    trial = @benchmark $f() evals = 1 samples = samples
    tmin = minimum(trial).time * 1e-9
    tmed = median(trial).time * 1e-9
    times = sort!(Float64.(trial.times) .* 1e-9)
    tci_low = quantile(times, 0.025)
    tci_high = quantile(times, 0.975)
    alloc = minimum(trial).memory
    return tmin, tmed, tci_low, tci_high, alloc
end

function _bench_and_quality_basic(fbench::F, quality::Q; warmup::Int, samples::Int) where {F<:Function,Q<:Function}
    tmin, tmed, alloc = bench_trial_basic(fbench; warmup = warmup, samples = samples)
    Ffact = fbench()
    resid, orth = quality(Ffact)
    return (tmin = tmin, tmed = tmed, alloc = alloc, resid = resid, orth = orth)
end

function _bench_and_quality_ci(fbench::F, quality::Q; warmup::Int, samples::Int) where {F<:Function,Q<:Function}
    tmin, tmed, tci_low, tci_high, alloc = bench_trial_ci(fbench; warmup = warmup, samples = samples)
    Ffact = fbench()
    resid, orth = quality(Ffact)
    return (
        tmin = tmin,
        tmed = tmed,
        tci_low = tci_low,
        tci_high = tci_high,
        alloc = alloc,
        resid = resid,
        orth = orth,
    )
end

function _bench_pair(
    A::Matrix{Float64},
    k::Int,
    norm_recomp_tol::Float64;
    warmup::Int,
    samples::Int,
    bench_quality::B,
    bsqr_return_rinv_r12::Bool = false,
    include_dgeqp3_trsm::Bool = false,
) where {B<:Function}
    rows = NamedTuple[]

    f_bs() = run_bsqr_fair(A, k, norm_recomp_tol; return_rinv_r12 = bsqr_return_rinv_r12)
    bs = bench_quality(f_bs, F -> residual_bs(A, F); warmup = warmup, samples = samples)
    push!(rows, (; method = BSQR_METHOD_LABEL, bs...))

    if lazy_benchmark_enabled()
        f_lazy() = run_bsqr_lazy_blas_fair(A, k, norm_recomp_tol)
        lazy = bench_quality(f_lazy, F -> residual_bs(A, F); warmup = warmup, samples = samples)
        push!(rows, (; method = BSQR_LAZY_BLAS_METHOD_LABEL, lazy...))
    end

    f_qr() = run_qr_fair(A)
    dg = bench_quality(f_qr, F -> residual_qr(A, F); warmup = warmup, samples = samples)
    push!(rows, (; method = DGEQP3_METHOD_LABEL, dg...))

    if include_dgeqp3_trsm
        f_qr_trsm() = run_qr_trsm_fair(A)
        dg_trsm = bench_quality(f_qr_trsm, F -> residual_qr(A, F); warmup = warmup, samples = samples)
        push!(rows, (; method = DGEQP3_TRSM_METHOD_LABEL, dg_trsm...))
    end

    return rows
end

function residual_bs(A::Matrix{Float64}, F::BSQRPivoted)
    Q = BSPivotQR._explicit_q(F)
    T = BSPivotQR._packed_to_qt(F)
    Aperm = A[:, perm(F)]
    resid = norm(Aperm - Q * T) / max(norm(A), eps(Float64))
    orth = norm(I - Q' * Q)
    return resid, orth
end

function residual_qr(A::Matrix{Float64}, F)
    p = Vector(F.p)
    Q = Matrix(F.Q)
    R = Matrix(F.R)
    resid = norm(A[:, p] - Q * R) / max(norm(A), eps(Float64))
    orth = norm(I - Q' * Q)
    return resid, orth
end

function run_bsqr_fair(
    A::Matrix{Float64},
    k::Int,
    norm_recomp_tol::Float64;
    return_rinv_r12::Bool = false,
)
    return bsqr!(
        copy(A);
        k = k,
        check = false,
        track_inverse_frob = false,
        return_rinv_r12 = return_rinv_r12,
        rank_stop = false,
        norm_recomp_tol = norm_recomp_tol,
        workspace = nothing,
        blas_threads = nothing,
    )
end

function run_bsqr_lazy_blas_fair(A::Matrix{Float64}, k::Int, norm_recomp_tol::Float64)
    _ = norm_recomp_tol
    kwargs = lazy_batch_kwargs(size(A, 2))
    return BSPivotQR._bsqr_lazy_blas(
        copy(A);
        k = k,
        check = false,
        rank_stop = false,
        blas_threads = nothing,
        kwargs...,
    )
end

run_qr_fair(A::Matrix{Float64}) = qr(copy(A), ColumnNorm())

function run_qr_trsm_fair(A::Matrix{Float64})
    F = run_qr_fair(A)
    m, n = size(A)
    k = min(m, n)
    if n > k
        Rf = Matrix(F.R)
        R11 = Matrix(view(Rf, 1:k, 1:k))
        R12 = Matrix(view(Rf, 1:k, (k + 1):n))
        BLAS.trsm!('L', 'U', 'N', 'N', 1.0, R11, R12)
    end
    return F
end

function bench_pair_basic(
    A::Matrix{Float64},
    k::Int,
    norm_recomp_tol::Float64;
    warmup::Int,
    samples::Int,
    bsqr_return_rinv_r12::Bool = false,
    include_dgeqp3_trsm::Bool = false,
)
    return _bench_pair(
        A,
        k,
        norm_recomp_tol;
        warmup = warmup,
        samples = samples,
        bench_quality = _bench_and_quality_basic,
        bsqr_return_rinv_r12 = bsqr_return_rinv_r12,
        include_dgeqp3_trsm = include_dgeqp3_trsm,
    )
end

function bench_pair_ci(
    A::Matrix{Float64},
    k::Int,
    norm_recomp_tol::Float64;
    warmup::Int,
    samples::Int,
    bsqr_return_rinv_r12::Bool = false,
    include_dgeqp3_trsm::Bool = false,
)
    return _bench_pair(
        A,
        k,
        norm_recomp_tol;
        warmup = warmup,
        samples = samples,
        bench_quality = _bench_and_quality_ci,
        bsqr_return_rinv_r12 = bsqr_return_rinv_r12,
        include_dgeqp3_trsm = include_dgeqp3_trsm,
    )
end

function grouped_rows_with_baseline(rows, keyf::F; sortby::S = identity) where {F<:Function,S<:Function}
    groups = sort(collect(unique(keyf(r) for r in rows)); by = sortby)
    out = NamedTuple[]
    for key in groups
        keyrows = filter(r -> keyf(r) == key, rows)
        drows = filter(r -> r.method == DGEQP3_METHOD_LABEL, keyrows)
        isempty(drows) && continue
        push!(out, (key = key, rows = sort(keyrows, by = x -> x.method), baseline = only(drows)))
    end
    return out
end

function make_matrix(family::Symbol, m::Int, n::Int, rng::AbstractRNG)
    if family === :gaussian
        return Matrix{Float64}(gaussian_matrix(m, n, rng))
    elseif family === :ill_conditioned
        return Matrix{Float64}(ill_conditioned_matrix(m, n, rng))
    elseif family === :orthonormal_rows
        m <= n || throw(ArgumentError("orthonormal_rows requires m <= n"))
        return Matrix{Float64}(orthonormal_row_matrix(m, n, rng))
    else
        error("Unknown family: $family")
    end
end

export DEFAULT_NORM_RECOMP_TOL
export BSQR_METHOD_LABEL, BSQR_LAZY_BLAS_METHOD_LABEL, DGEQP3_METHOD_LABEL, DGEQP3_TRSM_METHOD_LABEL
export accelerate_active, backend_string, bench_pair_basic, bench_pair_ci
export bench_trial_basic, bench_trial_ci, grouped_rows_with_baseline
export check_backend, configure_blas_threads, make_matrix
export lazy_batch_kwargs, lazy_benchmark_enabled
export parse_env_float, parse_env_int, parse_float_list, parse_int_list, parse_symbol_list
export residual_bs, residual_qr, run_bsqr_fair, run_bsqr_lazy_blas_fair, run_qr_fair, run_qr_trsm_fair

end
