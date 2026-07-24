# Statistical validation sweep for the bound guarantee -- opt-in (not
# included by runtests.jl; run it directly). Port of
# matlab_rand/tests/stress_bsqr_rand_bounds.m: sweeps seeds x kernel paths x
# sampling schemes x threshold modes x families and asserts, for every run on
# orthonormal-row input, that
#   * the factorization is exact (Q'A(:,p(1:k)) = R11, Q orthonormal), and
#   * the per-step bound f2_i <= Fhat_i holds at every step (hence the final
#     Osinsky bound ||R11^{-1}||_F <= sqrt(k(n-k+1))).
# Prints a per-configuration table of realized quality (frob_inv/osinsky
# ratio, samples tested per selection).
#
#   julia --project=julia_rand julia_rand/test/stress_bsqr_rand_bounds.jl [nseeds]

using LinearAlgebra
using Printf
using Random
using Statistics
using BSRandPivotQR

include(joinpath(@__DIR__, "..", "benchmark", "rand_test_matrix.jl"))

nseeds = isempty(ARGS) ? 100 : parse(Int, ARGS[1])
k, n = 16, 200
families = (:gaussian, :spiked_leverage, :needle)
bound2 = k * (n - k + 1)
bound_tol = 1e-8 * bound2

nrun = 0
@printf("%-18s %-8s %-14s %-20s %10s %10s %10s\n",
    "family", "batched", "sampling", "threshold", "ratio_med", "ratio_max", "tested/k")
for fam in families
    M = rand_test_matrix(fam, k, n, 4242)   # fixed matrix; only the RNG varies
    for batched in (true, false)
        for sampl in (:uniform, :normweighted)
            for mode in (:running_mean, :worstcase_allowance)
                ratios = Float64[]
                tested = Float64[]
                for s in 1:nseeds
                    global nrun += 1
                    F = bsqr_rand(M; k = k, seed = s, batched = batched,
                        sampling = sampl, threshold_mode = mode)
                    sel = selected(F)
                    Qm = Q(F)
                    scale = max(1.0, norm(M[:, sel]))
                    err = norm(Qm' * M[:, sel] - F.R11) / scale
                    orth = norm(Qm' * Qm - I)
                    err < 1e-10 || error("factorization inexact ($fam, seed=$s): $err")
                    orth < 1e-10 || error("Q not orthonormal ($fam, seed=$s): $orth")
                    all(F.stats.f2 .<= F.stats.Fhat .+ bound_tol) ||
                        error("per-step bound violated ($fam, batched=$batched, " *
                              "sampling=$sampl, mode=$mode, seed=$s)")
                    push!(ratios, F.stats.frob_inv / F.stats.osinsky_bound)
                    push!(tested, F.stats.total_tested / k)
                end
                @printf("%-18s %-8s %-14s %-20s %10.4f %10.4f %10.1f\n",
                    fam, batched, sampl, mode,
                    median(ratios), maximum(ratios), median(tested))
            end
        end
    end
end
println("\nAll $nrun runs: exact factorization, per-step bound intact.")
