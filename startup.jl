# startup.jl -- set up the BSQR Julia package for interactive use.
#
# From the repository root:
#     julia> include("startup.jl")
#
# It activates and instantiates the BSPivotQR environment under julia/ and
# brings the package into scope so the default call just works:
#
#     F = bsqr(A)                 # deterministic Bischof-Stewart pivoted QR
#     R(F), perm(F)               # R factor and column permutation
#     reconstruct(F, A)           # rebuild A from the factorization (up to round-off)
#
# The equivalent without this script:  `julia --project=julia` then
# `using BSPivotQR`.  The randomized variant (bsqr_rand) is MATLAB-only --
# see matlab_rand/ and the repository-root startup.m.

import Pkg
Pkg.activate(joinpath(@__DIR__, "julia"))
Pkg.instantiate()

using BSPivotQR

println("""

BSPivotQR ready (deterministic Bischof-Stewart pivoted QR).
  F = bsqr(A)                      BSQRPivoted factorization (see `?bsqr`)
  R(F)  perm(F)  reconstruct(F,A)  R factor, column permutation, rebuilt A
  rinv_r12(F)                      R11^{-1} R12 (with return_rinv_r12=true)
The randomized variant (bsqr_rand) is MATLAB-only; see startup.m.
""")
