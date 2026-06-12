#!/usr/bin/env julia
#
# Validates the publication benchmark CSV and invokes the Python figure/table
# generator (plot_publication.py) for both comparison modes. The figure spec
# lives in docs/PUBLICATION_FIGURES_PLAN.md; the MATLAB plotter conforms to
# the same spec.

using DelimitedFiles

const DEFAULT_INPUT = joinpath(@__DIR__, "results", "publication", "publication_timings.csv")
const DEFAULT_OUTDIR = joinpath(@__DIR__, "results", "publication", "plots")
const DEFAULT_TABLEDIR = joinpath(@__DIR__, "results", "publication", "tables")
const PLOTTER = joinpath(@__DIR__, "plot_publication.py")
const REQUIRED_COLUMNS = [
    "run_id",
    "timestamp",
    "family",
    "regime",
    "m",
    "n",
    "aspect",
    "seed",
    "blas_threads",
    "method",
    "tmin_s",
    "tmed_s",
    "tci_low_s",
    "tci_high_s",
    "alloc_bytes",
    "residual",
    "orthogonality",
]
const EXPECTED_METHODS = Set(["bsqr_full", "bsqr_rinv", "dgeqp3", "dgeqp3_trsm"])

_as_string(x) = String(x)
_as_int(x) = x isa Int ? x : parse(Int, _as_string(x))
_as_float(x) = x isa AbstractFloat ? Float64(x) : parse(Float64, _as_string(x))

function _validate_publication_csv(input::AbstractString)
    isfile(input) || error("Input CSV does not exist: $input")
    data, raw_header = readdlm(input, ',', header = true)
    header = _as_string.(vec(raw_header))
    isempty(header) && error("CSV header is empty: $input")

    hset = Set(header)
    for col in REQUIRED_COLUMNS
        col in hset || error("Missing required CSV column '$col' in $input")
    end
    nrows = size(data, 1)
    nrows > 0 || error("No benchmark rows found in $input")

    colidx = Dict{String,Int}(name => findfirst(==(name), header) for name in REQUIRED_COLUMNS)
    key_methods = Dict{Tuple{String,String,Int,Int,Int,Int},Set{String}}()
    methods_seen = Set{String}()
    for i in 1:nrows
        row = view(data, i, :)
        family = _as_string(row[colidx["family"]])
        regime = _as_string(row[colidx["regime"]])
        m = _as_int(row[colidx["m"]])
        n = _as_int(row[colidx["n"]])
        seed = _as_int(row[colidx["seed"]])
        nth = _as_int(row[colidx["blas_threads"]])
        method = _as_string(row[colidx["method"]])
        tmin = _as_float(row[colidx["tmin_s"]])
        tmed = _as_float(row[colidx["tmed_s"]])
        tci_low = _as_float(row[colidx["tci_low_s"]])
        tci_high = _as_float(row[colidx["tci_high_s"]])
        resid = _as_float(row[colidx["residual"]])
        orth = _as_float(row[colidx["orthogonality"]])

        all(isfinite, (tmin, tmed, tci_low, tci_high, resid, orth)) ||
            error("Non-finite numeric entry in row $i")
        tmin >= 0.0 || error("Negative tmin_s in row $i")
        tmed > 0.0 || error("Non-positive tmed_s in row $i")
        (tci_low <= tmed <= tci_high) ||
            error("Invalid CI ordering in row $i (tci_low <= tmed <= tci_high required)")

        push!(methods_seen, method)
        key = (family, regime, m, n, seed, nth)
        get!(key_methods, key, Set{String}())
        push!(key_methods[key], method)
    end

    methods_seen == EXPECTED_METHODS ||
        error("CSV methods mismatch: expected $(collect(EXPECTED_METHODS)), got $(collect(methods_seen))")
    for (key, got) in key_methods
        got == EXPECTED_METHODS ||
            error("Missing method rows for key=$key: expected $(collect(EXPECTED_METHODS)), got $(collect(got))")
    end
    return nrows
end

function main()
    input = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_INPUT
    outdir = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_OUTDIR
    tabledir = length(ARGS) >= 3 ? ARGS[3] : DEFAULT_TABLEDIR
    nrows = _validate_publication_csv(input)
    println("Validated publication CSV ($nrows rows): $input")
    mkpath(outdir)
    mkpath(tabledir)
    comparisons = [
        ("plain", "bsqr_full", "dgeqp3"),
        ("rinv", "bsqr_rinv", "dgeqp3_trsm"),
    ]
    for (mode, bsqr_method, baseline_method) in comparisons
        comp_outdir = joinpath(outdir, mode)
        comp_tabledir = joinpath(tabledir, mode)
        mkpath(comp_outdir)
        mkpath(comp_tabledir)
        run(`python3 $PLOTTER $input $comp_outdir $comp_tabledir $bsqr_method $baseline_method $mode`)
    end
end

main()
