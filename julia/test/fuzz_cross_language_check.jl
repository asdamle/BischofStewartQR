# Cross-language differential fuzz checker (opt-in; Phase 2 of
# docs/PUBLICATION_READINESS_PLAN.md). Not included by runtests.jl.
#
# Consumes a directory of randomized fixtures written by
# matlab/tests/fuzz_cross_language_gen.m (same format as parity/) and checks
# the Julia kernel against the MATLAB oracle's expected outputs on every
# case, on both execution paths (default dispatch and forced panel kernel).
# Any pivot-sequence divergence on these tie-screened inputs is a bug under
# the repo's lockstep contract.
#
#   BS_FUZZ_DIR=/path/to/fuzzdir julia --project=julia julia/test/fuzz_cross_language_check.jl

using DelimitedFiles
using LinearAlgebra: norm
using BSPivotQR

const FUZZ_DIR = get(ENV, "BS_FUZZ_DIR", "")
isempty(FUZZ_DIR) && error("set BS_FUZZ_DIR to the fuzz_cross_language_gen output directory")

read_fixture(name) = readdlm(joinpath(FUZZ_DIR, name), ',', Float64)
rel_err(X, Y) = norm(X - Y) / max(norm(Y), eps(Float64))

const PANEL_MODES = [
    ("default", Dict{String,String}()),
    ("forced panel", Dict("BS_PANEL_NB" => "8", "BS_PANEL_MIN_KN" => "0")),
]

manifest = readlines(joinpath(FUZZ_DIR, "manifest.csv"))
@assert first(manifest) == "name,m,n,k,rinv_check,gap_min,rtol_R,rtol_rinv,rtol_crit"

ncases = length(manifest) - 1
failures = String[]
for (mode_name, mode_env) in PANEL_MODES
    withenv(pairs(mode_env)...) do
        for line in manifest[2:end]
            fields = split(line, ',')
            name = String(fields[1])
            m, n, k = parse.(Int, fields[2:4])
            rinv_check = parse(Int, fields[5]) == 1
            rtol_R = parse(Float64, fields[7])
            rtol_rinv = parse(Float64, fields[8])
            rtol_crit = parse(Float64, fields[9])

            A = read_fixture("$(name)_A.csv")
            p_exp = Int.(vec(read_fixture("$(name)_p.csv")))
            R_exp = read_fixture("$(name)_R.csv")
            crit_exp = vec(read_fixture("$(name)_crit.csv"))

            F = bsqr(A; k = k, return_rinv_r12 = rinv_check, track_inverse_frob = true)

            crit_drift = maximum(abs.(F.frob_inv_trace .- crit_exp) ./
                                 max.(abs.(crit_exp), eps()))
            ok = F.ksteps == k &&
                 perm(F) == p_exp &&
                 rel_err(R(F), R_exp) < rtol_R &&
                 crit_drift < rtol_crit
            if ok && rinv_check
                rinv_exp = read_fixture("$(name)_rinv.csv")
                ok = rel_err(rinv_r12(F), rinv_exp) < rtol_rinv
            end
            if !ok
                push!(failures, "$name [$mode_name]: pivots match=$(perm(F) == p_exp), " *
                                "relR=$(rel_err(R(F), R_exp)), crit_drift=$crit_drift")
            end
        end
    end
end

println("fuzz check: $ncases cases x $(length(PANEL_MODES)) kernel paths, $(length(failures)) failures")
if !isempty(failures)
    foreach(f -> println("  FAIL ", f), failures[1:min(end, 20)])
    error("$(length(failures)) cross-language fuzz failures")
end
