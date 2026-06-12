# V3 cross-language parity (docs/VALIDATION_AND_PERF_PLAN.md): validate the
# Julia kernel against the committed fixtures under <repo_root>/parity.
#
# The fixtures are outputs of the MATLAB V1 oracle (matlab/tests/oracle_bsqr.m,
# a literal transcription of Algorithm 1 with no incremental recurrences),
# written by matlab/tests/generate_parity_fixtures.m with %.17g formatting so
# both languages read bit-identical inputs. Inputs are screened at generation
# time for criterion near-ties, which is what makes exact pivot-sequence
# equality a fair demand across BLAS runtimes. MATLAB consumes the same files
# in matlab/tests/test_parity_fixtures.m.

using DelimitedFiles
using LinearAlgebra: norm

const PARITY_DIR = normpath(joinpath(@__DIR__, "..", "..", "parity"))

read_fixture(name::AbstractString) =
    readdlm(joinpath(PARITY_DIR, name), ',', Float64)

rel_err(X, Y) = norm(X - Y) / max(norm(Y), eps(Float64))

# The fixture zoo sits below the panel kernel's k*n crossover, so the default
# dispatch runs the unblocked kernel. Each fixture is therefore checked twice:
# once at defaults and once with the panel kernel forced, keeping both
# execution paths pinned to the oracle.
const PANEL_MODES = [
    ("default", Dict{String,String}()),
    ("forced panel", Dict("BS_PANEL_NB" => "8", "BS_PANEL_MIN_KN" => "0")),
]

@testset "Cross-language parity fixtures (vs MATLAB oracle, $mode_name)" for (mode_name, mode_env) in PANEL_MODES
    @test isdir(PARITY_DIR)

    manifest = readlines(joinpath(PARITY_DIR, "manifest.csv"))
    @test first(manifest) == "name,m,n,k,rinv_check,gap_min,rtol_R,rtol_rinv,rtol_crit"

    withenv(pairs(mode_env)...) do
    for line in manifest[2:end]
        fields = split(line, ',')
        name = String(fields[1])
        m, n, k = parse.(Int, fields[2:4])
        rinv_check = parse(Int, fields[5]) == 1
        rtol_R = parse(Float64, fields[7])
        rtol_rinv = parse(Float64, fields[8])
        rtol_crit = parse(Float64, fields[9])

        @testset "$name" begin
            A = read_fixture("$(name)_A.csv")
            @test size(A) == (m, n)
            p_exp = Int.(vec(read_fixture("$(name)_p.csv")))
            R_exp = read_fixture("$(name)_R.csv")
            crit_exp = vec(read_fixture("$(name)_crit.csv"))

            F = bsqr(A; k = k, return_rinv_r12 = rinv_check,
                track_inverse_frob = true)

            @test F.ksteps == k
            @test perm(F) == p_exp
            @test rel_err(R(F), R_exp) < rtol_R

            # V4 trace parity: the per-step selected-pivot criterion values
            # (frob_inv_trace) against the oracle's from-scratch trace.
            crit = F.frob_inv_trace
            @test length(crit) == k
            @test maximum(abs.(crit .- crit_exp) ./ crit_exp) < rtol_crit

            if rinv_check
                rinv_exp = read_fixture("$(name)_rinv.csv")
                @test rel_err(rinv_r12(F), rinv_exp) < rtol_rinv
            end
        end
    end
    end
end
