"""
    R(F::BSQRPivoted) -> Matrix{Float64}

The `k×n` upper-trapezoidal factor of the *permuted* columns
(`k = F.ksteps`): `R(F)[1:k, 1:k]` is `R11`, `R(F)[:, k+1:n]` is `R12`.
For a full factorization (`k == min(m, n)`, the default),
`A[:, perm(F)] ≈ Q(F) * R(F)`. With early stop (`k < min(m, n)`) that holds
only for the selected block, `A[:, perm(F)[1:k]] ≈ Q(F) * R(F)[:, 1:k]`,
while `R(F)[:, k+1:n]` is the unselected columns' projection
`Q(F)' * A[:, perm(F)[k+1:n]]`.
"""
function R(F::BSQRPivoted)
    k = F.ksteps
    if k == 0
        return zeros(Float64, 0, size(F.factors, 2))
    end
    return triu(view(F.factors, 1:k, :))
end

"""
    Q(F::BSQRPivoted) -> Matrix{Float64}

The economy orthogonal factor: `m×k` with orthonormal columns
(`k = F.ksteps`), satisfying `A[:, perm(F)[1:k]] ≈ Q(F) * R(F)[:, 1:k]`.
`F` stores Q only implicitly (packed reflectors in `factors` + `tau`);
each call materializes it with LAPACK `orgqr`, one `O(m k^2)` pass.
"""
function Q(F::BSQRPivoted)
    m = size(F.factors, 1)
    k = F.ksteps
    k == 0 && return Matrix{Float64}(undef, m, 0)
    Qm = Matrix{Float64}(undef, m, k)
    copyto!(Qm, view(F.factors, :, 1:k))
    return LAPACK.orgqr!(Qm, view(F.tau, 1:k), k)
end

"""
    perm(F::BSQRPivoted) -> Vector{Int}

The length-`n` column permutation (a copy of `F.jpvt`): `perm(F)[1:F.ksteps]`
are the selected columns in pivot order, the rest are the unselected columns.
"""
perm(F::BSQRPivoted) = copy(F.jpvt)

"""
    rinv_r12(F::BSQRPivoted) -> Union{Nothing, Matrix{Float64}}

The `k×(n-k)` matrix `R11⁻¹R12` (a copy), extracted from the kernel workspace
at factorization time with no extra triangular solve — or `nothing` unless
the factorization was computed with `return_rinv_r12 = true`.
"""
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

"""
    reconstruct(F::BSQRPivoted, Aorig) -> Matrix{Float64}

Rebuild the factored matrix: forms the full `m×m` Q times the full reduced
matrix (all `n` columns, including any not-yet-eliminated trailing block) and
un-permutes the columns, so the result approximates the original `A` up to
roundoff for any `ksteps`. `Aorig` is used only to check the size. Intended
for validation and tests.
"""
function reconstruct(F::BSQRPivoted, Aorig::AbstractMatrix)
    size(Aorig) == size(F.factors) ||
        throw(DimensionMismatch("Aorig size does not match factorization"))

    Aperm = _explicit_q(F) * _packed_to_qt(F)
    Arec = similar(Aperm)
    @views Arec[:, F.jpvt] .= Aperm
    return Arec
end
