module BSPivotQR

using LinearAlgebra

import LinearAlgebra: Factorization

export bsqr, bsqr!, BSQRPivoted, R, perm, rinv_r12, reconstruct

include("workspace.jl")
include("kernel.jl")
include("kernel_panel.jl")
include("interface.jl")

end
