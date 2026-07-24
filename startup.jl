# startup.jl -- set up the BSQR Julia packages for interactive use.
#
# From the repository root:
#     julia> include("startup.jl")
#
# It activates and instantiates the BSPivotQR environment under julia/, makes
# the randomized-variant package under julia_rand/ loadable alongside it (via
# LOAD_PATH), and brings both into scope so the default calls just work:
#
#     F = bsqr(A)                 # deterministic Bischof-Stewart pivoted QR
#     R(F), perm(F)               # R factor and column permutation
#     reconstruct(F, A)           # rebuild A from the factorization (up to round-off)
#     G = bsqr_rand(M)            # randomized column selection (BSQRRandPivoted)
#
# The equivalent without this script:  `julia --project=julia` then
# `using BSPivotQR` (deterministic), or `julia --project=julia_rand` then
# `using BSRandPivotQR` (randomized).

import Pkg
Pkg.activate(joinpath(@__DIR__, "julia"))
Pkg.instantiate()

using BSPivotQR

# The randomized package is its own decoupled project; putting its root on
# LOAD_PATH makes it loadable without deactivating the BSPivotQR environment.
let randdir = joinpath(@__DIR__, "julia_rand")
    randdir in LOAD_PATH || push!(LOAD_PATH, randdir)
end

using BSRandPivotQR

println("""

BSPivotQR ready (deterministic Bischof-Stewart pivoted QR).
  F = bsqr(A)                      BSQRPivoted factorization (see `?bsqr`)
  R(F)  perm(F)  reconstruct(F,A)  R factor, column permutation, rebuilt A
  rinv_r12(F)                      R11^{-1} R12 (with return_rinv_r12=true)
BSRandPivotQR ready (randomized column selection, julia_rand/).
  G = bsqr_rand(M)                 BSQRRandPivoted (see `?bsqr_rand`)
  selected(G)  G.R11  G.stats      subset, triangular factor, instrumentation
Note: `Q` and `perm` are exported by both packages -- qualify on conflict
(e.g. BSRandPivotQR.Q(G)). The MATLAB twin lives in matlab_rand/ (startup.m).
""")
