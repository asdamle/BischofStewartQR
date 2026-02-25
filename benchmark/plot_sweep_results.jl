#!/usr/bin/env julia

# Intentionally use Python/matplotlib for stable behavior across environments.

const DEFAULT_INPUT = joinpath(@__DIR__, "results", "sweep_timings.csv")
const DEFAULT_OUTDIR = joinpath(@__DIR__, "results", "sweep_plots")

function _python_plot(input::AbstractString, outdir::AbstractString)
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
    rdr = csv.DictReader(f)
    for row in rdr:
        rows.append({
            "family": row["family"],
            "aspect": float(row["aspect"]),
            "n": int(row["n"]),
            "method": row["method"],
            "tmed": float(row["tmed_s"]),
            "tci_low": float(row["tci_low_s"]) if "tci_low_s" in row and row["tci_low_s"] else float(row["tmed_s"]),
            "tci_high": float(row["tci_high_s"]) if "tci_high_s" in row and row["tci_high_s"] else float(row["tmed_s"]),
        })

if not rows:
    raise RuntimeError(f"No rows in {input_path}")

def find_row(fam, aspect, n, method):
    for r in rows:
        if (
            r["family"] == fam
            and abs(r["aspect"] - aspect) < 1e-12
            and r["n"] == n
            and r["method"] == method
        ):
            return r
    return None

families = sorted({r["family"] for r in rows})
for fam in families:
    aspects = sorted({r["aspect"] for r in rows if r["family"] == fam})
    ns = sorted({r["n"] for r in rows if r["family"] == fam})

    if not aspects or not ns:
        continue

    zfull = np.full((len(aspects), len(ns)), np.nan)
    dge_ref = {}
    for ia, a in enumerate(aspects):
        dge_series = []
        for inx, n in enumerate(ns):
            d = find_row(fam, a, n, "dgeqp3")
            b = find_row(fam, a, n, "bsqr_full")
            if d is not None and b is not None and b["tmed"] != 0.0:
                zfull[ia, inx] = d["tmed"] / b["tmed"]
            dge_series.append(np.nan if d is None else d["tmed"])
        dge_ref[a] = dge_series

    plt.figure(figsize=(8, 4.5))
    im = plt.imshow(zfull, aspect="auto", origin="lower")
    plt.title(f"{fam}: speedup dgeqp3 / bsqr_full")
    plt.xlabel("n")
    plt.ylabel("aspect m/n")
    plt.xticks(range(len(ns)), [str(v) for v in ns], rotation=45, ha="right")
    plt.yticks(range(len(aspects)), [f"{a:.3f}" for a in aspects])
    cb = plt.colorbar(im)
    cb.set_label("speedup")
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{fam}_speedup_full_heatmap.png"), dpi=170)
    plt.close()

    plt.figure(figsize=(9, 4.5))
    for a in aspects:
        y_bs, lo_bs, hi_bs = [], [], []
        y_dg, lo_dg, hi_dg = [], [], []
        for n in ns:
            b = find_row(fam, a, n, "bsqr_full")
            d = find_row(fam, a, n, "dgeqp3")
            if b is None:
                y_bs.append(np.nan)
                lo_bs.append(np.nan)
                hi_bs.append(np.nan)
            else:
                y_bs.append(b["tmed"])
                lo_bs.append(b["tci_low"])
                hi_bs.append(b["tci_high"])
            if d is None:
                y_dg.append(np.nan)
                lo_dg.append(np.nan)
                hi_dg.append(np.nan)
            else:
                y_dg.append(d["tmed"])
                lo_dg.append(d["tci_low"])
                hi_dg.append(d["tci_high"])
        plt.plot(ns, y_bs, marker="o", label=f"bsqr_full a={a:.3f}")
        plt.fill_between(ns, lo_bs, hi_bs, alpha=0.18)
        plt.plot(ns, y_dg, marker="D", linestyle="--", label=f"dgeqp3 a={a:.3f}")
        plt.fill_between(ns, lo_dg, hi_dg, alpha=0.14)

    ref_x = np.asarray(ns, dtype=float)
    ref_y = None
    for a in aspects:
        y = np.asarray(dge_ref[a], dtype=float)
        finite = np.isfinite(y)
        if finite.any():
            i0 = int(np.where(finite)[0][0])
            n0 = ref_x[i0]
            y0 = y[i0]
            ref_y = y0 * (ref_x / n0) ** 3
            break
    if ref_y is not None:
        plt.plot(ref_x, ref_y, color="black", linestyle=":", linewidth=1.3, label="n^3 reference")

    plt.xscale("log", base=2)
    plt.yscale("log")
    plt.title(f"{fam}: median time vs n by aspect (log-log)")
    plt.xlabel("n")
    plt.ylabel("median time (s)")
    plt.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), fontsize=8)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f"{fam}_timing_lines.png"), dpi=170)
    plt.close()

print(f"Wrote {len(families)} family sweep-plot sets to {outdir}")
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
