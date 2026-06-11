function R(F::BSQRPivoted)
    k = F.ksteps
    if k == 0
        return zeros(Float64, 0, size(F.factors, 2))
    end
    return triu(view(F.factors, 1:k, :))
end

perm(F::BSQRPivoted) = copy(F.jpvt)
rinv_r12(F::BSQRPivoted) = F.rinv_r12 === nothing ? nothing : copy(F.rinv_r12)

function _packed_to_qt(F::BSQRPivoted)
    T = copy(F.factors)
    m = size(T, 1)
    # Factors store Householder vectors below the diagonal; zero those entries
    # to recover the explicit upper-trapezoidal factor T used in A[:,p] = Q*T.
    for i in 1:F.ksteps
        if i < m
            fill!(view(T, (i + 1):m, i), 0.0)
        end
    end
    return T
end

function _explicit_q(F::BSQRPivoted)
    m = size(F.factors, 1)
    k = F.ksteps
    k == 0 && return Matrix{Float64}(I, m, m)

    # orgqr generates Q directly from the stored reflectors (~1/3 fewer
    # flops than applying them to an identity via ormqr) and matches how
    # the dgeqp3 baseline materializes its Q (P1.1 in docs/PERF_P0_FINDINGS.md).
    Q = Matrix{Float64}(undef, m, m)
    copyto!(view(Q, :, 1:k), view(F.factors, :, 1:k))
    return LAPACK.orgqr!(Q, view(F.tau, 1:k), k)
end

function reconstruct(F::BSQRPivoted, Aorig::AbstractMatrix)
    size(Aorig) == size(F.factors) ||
        throw(DimensionMismatch("Aorig size does not match factorization"))

    Aperm = _explicit_q(F) * _packed_to_qt(F)
    Arec = similar(Aperm)
    @views Arec[:, F.jpvt] .= Aperm
    return Arec
end
