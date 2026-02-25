#!/usr/bin/env julia

# Intentionally use Python/matplotlib for stable behavior across environments.

const DEFAULT_INPUT = joinpath(@__DIR__, "results", "timings.csv")
const DEFAULT_OUTDIR = joinpath(@__DIR__, "results", "plots")

function _python_plot(input::AbstractString, outdir::AbstractString)
    py = raw"""
import csv
import math
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

input_path = sys.argv[1]
outdir = sys.argv[2]
os.makedirs(outdir, exist_ok=True)

rows = []
with open(input_path, newline="") as f:
    rdr = csv.DictReader(f)
    for row in rdr:
        rows.append({
            "family": row["family"],
            "m": int(row["m"]),
            "n": int(row["n"]),
            "method": row["method"],
            "tmed": float(row["tmed_s"]),
            "residual": float(row["residual"]),
            "orthogonality": float(row["orthogonality"]),
        })

if not rows:
    raise RuntimeError(f"No rows in {input_path}")

families = sorted({r["family"] for r in rows})
for fam in families:
    keys = sorted({(r["m"], r["n"]) for r in rows if r["family"] == fam})
    labels = [f"{m}x{n}" for (m, n) in keys]
    x = list(range(len(keys)))

    def series(method, field):
        out = []
        for (m, n) in keys:
            vals = [r[field] for r in rows if r["family"] == fam and r["m"] == m and r["n"] == n and r["method"] == method]
            out.append(vals[0] if vals else float("nan"))
        return out

    t_dg = series("dgeqp3", "tmed")
    t_bs = series("bsqr_full", "tmed")
    r_dg = series("dgeqp3", "residual")
    r_bs = series("bsqr_full", "residual")
    o_dg = series("dgeqp3", "orthogonality")
    o_bs = series("bsqr_full", "orthogonality")

    width = 0.35
    plt.figure(figsize=(10, 4))
    plt.bar([i - 0.5 * width for i in x], t_dg, width=width, label="dgeqp3")
    plt.bar([i + 0.5 * width for i in x], t_bs, width=width, label="bsqr_full")
    plt.title(f"{fam}: timing")
    plt.xlabel("matrix size")
    plt.ylabel("median time (s)")
    plt.xticks(x, labels, rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{fam}_timing.png"), dpi=170)
    plt.close()

    plt.figure(figsize=(10, 4))
    plt.semilogy(x, r_dg, marker="o", label="dgeqp3 residual")
    plt.semilogy(x, r_bs, marker="D", label="bsqr_full residual")
    plt.semilogy(x, o_dg, marker="^", label="dgeqp3 orthogonality")
    plt.semilogy(x, o_bs, marker="v", label="bsqr_full orthogonality")
    plt.title(f"{fam}: quality")
    plt.xlabel("matrix size")
    plt.ylabel("error")
    plt.xticks(x, labels, rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{fam}_quality.png"), dpi=170)
    plt.close()

    speed = [a / b if (not math.isnan(a) and not math.isnan(b) and b != 0.0) else float("nan") for a, b in zip(t_dg, t_bs)]
    plt.figure(figsize=(10, 4))
    plt.bar(x, speed, width=0.45, label="dgeqp3 / bsqr_full")
    plt.axhline(1.0, color="black", linestyle="--", label="parity")
    plt.title(f"{fam}: speedup (>1 means BS faster)")
    plt.xlabel("matrix size")
    plt.ylabel("speedup")
    plt.xticks(x, labels, rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{fam}_speedup.png"), dpi=170)
    plt.close()

print(f"Wrote {len(families)} family plot sets to {outdir}")
"""

    run(`python3 -c $py $input $outdir`)
end

function main()
    input = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_INPUT
    outdir = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_OUTDIR
    mkpath(outdir)
    _python_plot(input, outdir)
end

main()
