module BSRandPivotQR

using LinearAlgebra
using Random

import LinearAlgebra: Factorization

export bsqr_rand, BSQRRandPivoted, BSQRRandStats, RankDeficientError,
    Q, R11, perm, selected, r12

include("helpers.jl")
include("kernel_rand.jl")
include("interface.jl")

end
