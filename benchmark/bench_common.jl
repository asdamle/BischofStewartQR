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

function run_bsqr_fair(A::Matrix{Float64}, k::Int, norm_recomp_tol::Float64)
    return bsqr!(
        copy(A);
        k = k,
        check = false,
        track_inverse_frob = false,
        return_rinv_r12 = false,
        rank_stop = false,
        norm_recomp_tol = norm_recomp_tol,
        workspace = nothing,
        blas_threads = nothing,
    )
end

run_qr_fair(A::Matrix{Float64}) = qr(copy(A), ColumnNorm())

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
export accelerate_active, backend_string, bench_trial_basic, bench_trial_ci
export check_backend, configure_blas_threads, make_matrix
export parse_env_float, parse_env_int, parse_float_list, parse_int_list, parse_symbol_list
export residual_bs, residual_qr, run_bsqr_fair, run_qr_fair

end
