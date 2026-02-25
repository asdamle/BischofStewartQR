module MatrixGenerators

using LinearAlgebra
using Random

export gaussian_matrix, ill_conditioned_matrix, orthonormal_row_matrix

gaussian_matrix(m::Int, n::Int, rng::AbstractRNG) = randn(rng, m, n)

function ill_conditioned_matrix(m::Int, n::Int, rng::AbstractRNG; kappa::Float64 = 1.0e10)
    r = min(m, n)
    U = Matrix(qr(randn(rng, m, r)).Q)
    V = Matrix(qr(randn(rng, n, r)).Q)
    s = exp.(range(0.0, stop = -log(kappa), length = r))
    return U * Diagonal(s) * V'
end

function orthonormal_row_matrix(m::Int, n::Int, rng::AbstractRNG)
    m <= n || throw(ArgumentError("orthonormal_row_matrix requires m <= n"))
    Q = Matrix(qr(randn(rng, n, m)).Q)
    return Matrix(Q[:, 1:m]')
end

end
