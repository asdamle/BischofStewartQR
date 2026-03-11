#!/usr/bin/env julia

const DEFAULT_INPUT = joinpath(@__DIR__, "results", "publication", "publication_timings.csv")
const DEFAULT_OUTDIR = joinpath(@__DIR__, "results", "publication", "plots")
const DEFAULT_TABLEDIR = joinpath(@__DIR__, "results", "publication", "tables")

function _python_plot(input::AbstractString, outdir::AbstractString, tabledir::AbstractString)
    py = raw"""
import csv
import math
import os
import sys
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np

input_path = sys.argv[1]
outdir = sys.argv[2]
tabledir = sys.argv[3]
os.makedirs(outdir, exist_ok=True)
os.makedirs(tabledir, exist_ok=True)

rows = []
with open(input_path, newline="") as f:
    r = csv.DictReader(f)
    for row in r:
        rows.append({
            "run_id": row["run_id"],
            "timestamp": row["timestamp"],
            "family": row["family"],
            "regime": row["regime"],
            "m": int(row["m"]),
            "n": int(row["n"]),
            "aspect": float(row["aspect"]),
            "seed": int(row["seed"]),
            "blas_threads": int(row["blas_threads"]),
            "method": row["method"],
            "tmin": float(row["tmin_s"]),
            "tmed": float(row["tmed_s"]),
            "tci_low": float(row["tci_low_s"]),
            "tci_high": float(row["tci_high_s"]),
            "alloc_bytes": int(row["alloc_bytes"]),
            "residual": float(row["residual"]),
            "orthogonality": float(row["orthogonality"]),
        })

if not rows:
    raise RuntimeError(f"No rows found in {input_path}")

families = sorted({r["family"] for r in rows})
threads = sorted({r["blas_threads"] for r in rows})
methods = ("bsqr_full", "dgeqp3")
regimes = ("square", "short_wide")

if set(r["method"] for r in rows) - set(methods):
    raise RuntimeError("Unexpected methods in publication CSV; expected only bsqr_full and dgeqp3")

plt.rcParams.update({
    "font.size": 8,
    "axes.titlesize": 8,
    "axes.labelsize": 8,
    "legend.fontsize": 7,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
})

def finite(vals):
    return [v for v in vals if math.isfinite(v)]

def geomean(vals):
    good = [v for v in vals if math.isfinite(v) and v > 0.0]
    if not good:
        return float("nan")
    return math.exp(sum(math.log(v) for v in good) / len(good))

def grouped_rows(rows_in, keys):
    out = defaultdict(list)
    for r in rows_in:
        out[tuple(r[k] for k in keys)].append(r)
    return out

def median_min_max(vals):
    v = finite(vals)
    if not v:
        return float("nan"), float("nan"), float("nan")
    arr = np.asarray(v, dtype=float)
    return float(np.median(arr)), float(np.min(arr)), float(np.max(arr))

def paired_speedup_rows(rows_in):
    idx = {}
    for r in rows_in:
        key = (
            r["family"], r["regime"], r["m"], r["n"], r["aspect"],
            r["seed"], r["blas_threads"], r["method"]
        )
        idx[key] = r

    out = []
    for r in rows_in:
        if r["method"] != "bsqr_full":
            continue
        key_dg = (
            r["family"], r["regime"], r["m"], r["n"], r["aspect"],
            r["seed"], r["blas_threads"], "dgeqp3"
        )
        d = idx.get(key_dg)
        if d is None or r["tmed"] == 0.0:
            continue
        out.append({
            "family": r["family"],
            "regime": r["regime"],
            "m": r["m"],
            "n": r["n"],
            "aspect": r["aspect"],
            "seed": r["seed"],
            "blas_threads": r["blas_threads"],
            "speedup": d["tmed"] / r["tmed"],
            "bsqr_tmed": r["tmed"],
            "dgeqp3_tmed": d["tmed"],
        })
    return out

speed_rows = paired_speedup_rows(rows)

def save_fig(fig, stem):
    fig.savefig(os.path.join(outdir, f"{stem}.png"), dpi=300)
    fig.savefig(os.path.join(outdir, f"{stem}.pdf"))
    plt.close(fig)

def mk_axes_grid(nr, nc, width_scale=2.7, height_scale=2.1):
    fig, axes = plt.subplots(
        nr,
        nc,
        figsize=(max(3.6, width_scale * nc), max(4.0, height_scale * nr)),
        squeeze=False,
        constrained_layout=True,
    )
    return fig, axes

# Figure 1: square runtime log-log (per family, faceted by threads)
fig1, axes1 = mk_axes_grid(len(families), len(threads), width_scale=2.6, height_scale=2.0)
for i, fam in enumerate(families):
    for j, th in enumerate(threads):
        ax = axes1[i][j]
        sub = [
            r for r in rows
            if r["family"] == fam and r["regime"] == "square" and r["blas_threads"] == th
        ]
        ms = sorted({r["m"] for r in sub})
        for method, marker, ls in [("bsqr_full", "o", "-"), ("dgeqp3", "D", "--")]:
            ymed = []
            ylow = []
            yhigh = []
            for m in ms:
                vals = [r["tmed"] for r in sub if r["m"] == m and r["method"] == method]
                med, lo, hi = median_min_max(vals)
                ymed.append(med)
                ylow.append(lo)
                yhigh.append(hi)
            ax.plot(ms, ymed, marker=marker, linestyle=ls, label=method)
            ax.fill_between(ms, ylow, yhigh, alpha=0.18)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_title(f"{fam}, BLAS={th}")
        if i == len(families) - 1:
            ax.set_xlabel("m=n")
        if j == 0:
            ax.set_ylabel("median time (s)")
        ax.grid(True, alpha=0.25)
axes1[0][0].legend(loc="best")
save_fig(fig1, "figure1_square_runtime")

# Figure 2: short-wide runtime vs n (fixed m curves, per family, faceted by threads)
fig2, axes2 = mk_axes_grid(len(families), len(threads), width_scale=2.9, height_scale=2.0)
for i, fam in enumerate(families):
    for j, th in enumerate(threads):
        ax = axes2[i][j]
        sub = [
            r for r in rows
            if r["family"] == fam and r["regime"] == "short_wide" and r["blas_threads"] == th
        ]
        ms = sorted({r["m"] for r in sub})
        cmap = plt.get_cmap("tab10")
        for midx, m in enumerate(ms):
            color = cmap(midx % 10)
            for method, marker, ls in [("bsqr_full", "o", "-"), ("dgeqp3", "D", "--")]:
                rn = [r for r in sub if r["m"] == m and r["method"] == method]
                ns = sorted({r["n"] for r in rn})
                ymed = []
                ylow = []
                yhigh = []
                for n in ns:
                    vals = [r["tmed"] for r in rn if r["n"] == n]
                    med, lo, hi = median_min_max(vals)
                    ymed.append(med)
                    ylow.append(lo)
                    yhigh.append(hi)
                label = f"m={m} {method}"
                ax.plot(ns, ymed, marker=marker, linestyle=ls, color=color, label=label)
                ax.fill_between(ns, ylow, yhigh, color=color, alpha=0.12)
        ax.set_xscale("log", base=2)
        ax.set_yscale("log")
        ax.set_title(f"{fam}, BLAS={th}")
        if i == len(families) - 1:
            ax.set_xlabel("n")
        if j == 0:
            ax.set_ylabel("median time (s)")
        ax.grid(True, alpha=0.25)
handles, labels = axes2[0][0].get_legend_handles_labels()
if handles:
    fig2.legend(handles, labels, loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.02))
save_fig(fig2, "figure2_shortwide_runtime")

# Figure 3: short-wide speedup heatmaps over (m, n/m) per family/thread
fig3, axes3 = mk_axes_grid(len(families), len(threads), width_scale=2.7, height_scale=2.0)
all_speed = [s["speedup"] for s in speed_rows if s["regime"] == "short_wide" and math.isfinite(s["speedup"])]
vmin = np.nanmin(all_speed) if all_speed else 0.5
vmax = np.nanmax(all_speed) if all_speed else 1.5
if vmin == vmax:
    vmin = min(vmin, 0.95)
    vmax = max(vmax, 1.05)
if vmax <= 1.0:
    vmax = 1.01
if vmin >= 1.0:
    vmin = 0.99
norm = mcolors.TwoSlopeNorm(vmin=vmin, vcenter=1.0, vmax=vmax)
for i, fam in enumerate(families):
    for j, th in enumerate(threads):
        ax = axes3[i][j]
        sub = [
            s for s in speed_rows
            if s["family"] == fam and s["blas_threads"] == th and s["regime"] == "short_wide"
        ]
        ms = sorted({s["m"] for s in sub})
        aspects = sorted({s["aspect"] for s in sub})
        Z = np.full((len(ms), len(aspects)), np.nan)
        for im, m in enumerate(ms):
            for ia, aspect in enumerate(aspects):
                vals = [s["speedup"] for s in sub if s["m"] == m and abs(s["aspect"] - aspect) < 1e-12]
                Z[im, ia] = geomean(vals)
        imh = ax.imshow(Z, origin="lower", aspect="auto", cmap="RdYlGn", norm=norm)
        ax.set_title(f"{fam}, BLAS={th}")
        ax.set_xticks(range(len(aspects)))
        ax.set_xticklabels([f"{a:.1f}" for a in aspects], rotation=45, ha="right")
        ax.set_yticks(range(len(ms)))
        ax.set_yticklabels([str(m) for m in ms])
        if i == len(families) - 1:
            ax.set_xlabel("aspect (n/m)")
        if j == 0:
            ax.set_ylabel("m")
cbar = fig3.colorbar(imh, ax=axes3.ravel().tolist(), shrink=0.8)
cbar.set_label("speedup dgeqp3/bsqr")
save_fig(fig3, "figure3_shortwide_speedup_heatmap")

# Figure 4: quality plot (residual and orthogonality, log-scale) across regimes
def regime_method_stats(metric):
    labels = []
    meds = []
    lows = []
    highs = []
    for regime in regimes:
        for method in methods:
            vals = [r[metric] for r in rows if r["regime"] == regime and r["method"] == method]
            if vals:
                arr = np.asarray(vals, dtype=float)
                labels.append(f"{regime}\\n{method}")
                meds.append(float(np.median(arr)))
                lows.append(float(np.percentile(arr, 10)))
                highs.append(float(np.percentile(arr, 90)))
    return labels, meds, lows, highs

fig4, axes4 = plt.subplots(1, 2, figsize=(7.0, 3.0), squeeze=False, constrained_layout=True)
for ax, metric, title in [
    (axes4[0][0], "residual", "Residual"),
    (axes4[0][1], "orthogonality", "Orthogonality"),
]:
    labels, meds, lows, highs = regime_method_stats(metric)
    x = np.arange(len(labels))
    yerr = np.vstack([np.asarray(meds) - np.asarray(lows), np.asarray(highs) - np.asarray(meds)])
    ax.bar(x, meds, width=0.7, color="#4C78A8")
    ax.errorbar(x, meds, yerr=yerr, fmt="none", ecolor="black", capsize=3, linewidth=0.8)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.25)
save_fig(fig4, "figure4_quality")

# Figure 5: aggregate speedup bars with per-seed intervals
seed_group = defaultdict(list)
for s in speed_rows:
    key = (s["family"], s["regime"], s["blas_threads"], s["seed"])
    seed_group[key].append(s["speedup"])

combo_seed_vals = defaultdict(list)
for (fam, regime, th, seed), vals in seed_group.items():
    combo_seed_vals[(fam, regime, th)].append(geomean(vals))

labels = []
centers = []
lows = []
highs = []
for key in sorted(combo_seed_vals.keys()):
    vals = finite(combo_seed_vals[key])
    if not vals:
        continue
    center = geomean(vals)
    lo = min(vals)
    hi = max(vals)
    labels.append(f"{key[0]}\\n{key[1]}\\nBLAS={key[2]}")
    centers.append(center)
    lows.append(lo)
    highs.append(hi)

fig5, ax5 = plt.subplots(1, 1, figsize=(7.4, 3.4), constrained_layout=True)
x = np.arange(len(labels))
yerr = np.vstack([np.asarray(centers) - np.asarray(lows), np.asarray(highs) - np.asarray(centers)])
ax5.bar(x, centers, width=0.7, color="#72B7B2")
ax5.errorbar(x, centers, yerr=yerr, fmt="none", ecolor="black", capsize=3, linewidth=0.8)
ax5.axhline(1.0, color="black", linestyle="--", linewidth=1.0)
ax5.set_xticks(x)
ax5.set_xticklabels(labels)
ax5.set_ylabel("geomean speedup (dgeqp3/bsqr)")
ax5.set_title("Aggregate speedup by family/regime/thread")
ax5.grid(True, axis="y", alpha=0.25)
save_fig(fig5, "figure5_aggregate_speedup")

# Caption-ready tables
def write_csv(path, columns, rows_out):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        for row in rows_out:
            w.writerow(row)

square_table = []
sq_groups = grouped_rows([s for s in speed_rows if s["regime"] == "square"], ["family", "blas_threads", "m", "n"])
for key in sorted(sq_groups.keys()):
    vals = sq_groups[key]
    speedups = [v["speedup"] for v in vals]
    bs = [v["bsqr_tmed"] for v in vals]
    dg = [v["dgeqp3_tmed"] for v in vals]
    square_table.append({
        "family": key[0],
        "blas_threads": key[1],
        "m": key[2],
        "n": key[3],
        "speedup_geomean": geomean(speedups),
        "speedup_seed_min": min(finite(speedups)) if finite(speedups) else float("nan"),
        "speedup_seed_max": max(finite(speedups)) if finite(speedups) else float("nan"),
        "bsqr_tmed_geomean_s": geomean(bs),
        "dgeqp3_tmed_geomean_s": geomean(dg),
    })

short_table = []
sw_groups = grouped_rows([s for s in speed_rows if s["regime"] == "short_wide"], ["family", "blas_threads", "m", "n", "aspect"])
for key in sorted(sw_groups.keys()):
    vals = sw_groups[key]
    speedups = [v["speedup"] for v in vals]
    bs = [v["bsqr_tmed"] for v in vals]
    dg = [v["dgeqp3_tmed"] for v in vals]
    short_table.append({
        "family": key[0],
        "blas_threads": key[1],
        "m": key[2],
        "n": key[3],
        "aspect": key[4],
        "speedup_geomean": geomean(speedups),
        "speedup_seed_min": min(finite(speedups)) if finite(speedups) else float("nan"),
        "speedup_seed_max": max(finite(speedups)) if finite(speedups) else float("nan"),
        "bsqr_tmed_geomean_s": geomean(bs),
        "dgeqp3_tmed_geomean_s": geomean(dg),
    })

quality_table = []
q_groups = grouped_rows(rows, ["family", "blas_threads", "regime", "method"])
for key in sorted(q_groups.keys()):
    vals = q_groups[key]
    residuals = np.asarray([v["residual"] for v in vals], dtype=float)
    orths = np.asarray([v["orthogonality"] for v in vals], dtype=float)
    quality_table.append({
        "family": key[0],
        "blas_threads": key[1],
        "regime": key[2],
        "method": key[3],
        "residual_median": float(np.median(residuals)),
        "residual_p95": float(np.percentile(residuals, 95)),
        "orthogonality_median": float(np.median(orths)),
        "orthogonality_p95": float(np.percentile(orths, 95)),
    })

write_csv(
    os.path.join(tabledir, "table_square_speedup.csv"),
    [
        "family", "blas_threads", "m", "n",
        "speedup_geomean", "speedup_seed_min", "speedup_seed_max",
        "bsqr_tmed_geomean_s", "dgeqp3_tmed_geomean_s",
    ],
    square_table,
)
write_csv(
    os.path.join(tabledir, "table_shortwide_speedup.csv"),
    [
        "family", "blas_threads", "m", "n", "aspect",
        "speedup_geomean", "speedup_seed_min", "speedup_seed_max",
        "bsqr_tmed_geomean_s", "dgeqp3_tmed_geomean_s",
    ],
    short_table,
)
write_csv(
    os.path.join(tabledir, "table_quality.csv"),
    [
        "family", "blas_threads", "regime", "method",
        "residual_median", "residual_p95", "orthogonality_median", "orthogonality_p95",
    ],
    quality_table,
)

print(f"Wrote publication figures to {outdir}")
print(f"Wrote publication tables to {tabledir}")
"""

    run(`python3 -c $py $input $outdir $tabledir`)
end

function main()
    input = length(ARGS) >= 1 ? ARGS[1] : DEFAULT_INPUT
    outdir = length(ARGS) >= 2 ? ARGS[2] : DEFAULT_OUTDIR
    tabledir = length(ARGS) >= 3 ? ARGS[3] : DEFAULT_TABLEDIR
    mkpath(outdir)
    mkpath(tabledir)
    _python_plot(input, outdir, tabledir)
end

main()
