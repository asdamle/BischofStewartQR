#!/usr/bin/env julia

const DEFAULT_INPUT = joinpath(@__DIR__, "results", "regime_timings.csv")
const DEFAULT_OUTDIR = joinpath(@__DIR__, "results", "regime_plots")

function _python_fallback(input::AbstractString, outdir::AbstractString)
    py = raw"""
import csv
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

input_path = sys.argv[1]
outdir = sys.argv[2]
os.makedirs(outdir, exist_ok=True)

rows = []
with open(input_path, newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        rows.append({
            "family": row["family"],
            "regime": row["regime"],
            "fixed_value": int(row["fixed_value"]),
            "var_value": int(row["var_value"]),
            "m": int(row["m"]),
            "n": int(row["n"]),
            "method": row["method"],
            "tmed": float(row["tmed_s"]),
            "tci_low": float(row["tci_low_s"]) if "tci_low_s" in row and row["tci_low_s"] else float(row["tmed_s"]),
            "tci_high": float(row["tci_high_s"]) if "tci_high_s" in row and row["tci_high_s"] else float(row["tmed_s"]),
        })

if not rows:
    raise RuntimeError(f"No rows found in {input_path}")

families = sorted({r["family"] for r in rows})
regimes = ("fixed_n_vary_m", "fixed_m_vary_n")
methods = ("bsqr_full", "dgeqp3")
linestyles = {"bsqr_full": "-", "dgeqp3": "--"}
markers = {"bsqr_full": "o", "dgeqp3": "D"}
alphas = {"bsqr_full": 0.18, "dgeqp3": 0.14}

def find_row(fam, regime, fixed, var, method):
    for r in rows:
        if (
            r["family"] == fam
            and r["regime"] == regime
            and r["fixed_value"] == fixed
            and r["var_value"] == var
            and r["method"] == method
        ):
            return r
    return None

for fam in families:
    for regime in regimes:
        rsel = [r for r in rows if r["family"] == fam and r["regime"] == regime]
        if not rsel:
            continue

        fixed_values = sorted({r["fixed_value"] for r in rsel})
        cmap = plt.get_cmap("tab10")

        plt.figure(figsize=(10, 4.8))
        ref_x = None
        ref_y = None

        for i, fixed in enumerate(fixed_values):
            var_values = sorted({r["var_value"] for r in rsel if r["fixed_value"] == fixed})
            color = cmap(i % 10)
            for method in methods:
                ys = []
                yl = []
                yh = []
                for v in var_values:
                    rr = find_row(fam, regime, fixed, v, method)
                    if rr is None:
                        ys.append(np.nan)
                        yl.append(np.nan)
                        yh.append(np.nan)
                    else:
                        ys.append(rr["tmed"])
                        yl.append(rr["tci_low"])
                        yh.append(rr["tci_high"])
                if regime == "fixed_n_vary_m":
                    fixed_label = f"n={fixed}"
                else:
                    fixed_label = f"m={fixed}"
                plt.plot(
                    var_values,
                    ys,
                    marker=markers[method],
                    linestyle=linestyles[method],
                    color=color,
                    label=f"{method} ({fixed_label})",
                )
                plt.fill_between(var_values, yl, yh, color=color, alpha=alphas[method])

                if method == "dgeqp3" and ref_x is None:
                    yarr = np.asarray(ys, dtype=float)
                    finite = np.isfinite(yarr)
                    if finite.any():
                        j0 = int(np.where(finite)[0][0])
                        x0 = float(var_values[j0])
                        y0 = float(yarr[j0])
                        xarr = np.asarray(var_values, dtype=float)
                        ref_x = xarr
                        ref_y = y0 * (xarr / x0)

        if ref_x is not None:
            plt.plot(ref_x, ref_y, color="black", linestyle=":", linewidth=1.3, label="linear reference")

        plt.xscale("log", base=2)
        plt.yscale("log")
        if regime == "fixed_n_vary_m":
            plt.title(f"{fam}: fixed n, increasing m (log-log)")
            plt.xlabel("m")
            fname = f"{fam}_fixed_n_vary_m_loglog.png"
        else:
            plt.title(f"{fam}: fixed m, increasing n (short-wide focus, log-log)")
            plt.xlabel("n")
            fname = f"{fam}_fixed_m_vary_n_loglog.png"
        plt.ylabel("median time (s)")
        plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8)
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, fname), dpi=170)
        plt.close()

print(f"Wrote regime plots to {outdir}")
"""

    cmd = `python3 -c $py $input $outdir`
    run(cmd)
end

function main()
    input = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_INPUT
    outdir = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_OUTDIR
    mkpath(outdir)
    _python_fallback(input, outdir)
end

main()
