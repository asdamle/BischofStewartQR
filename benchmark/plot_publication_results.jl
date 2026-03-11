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
from matplotlib.lines import Line2D
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
    "font.size": 10,
    "axes.titlesize": 10,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
    "legend.title_fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
})

layout = os.getenv("BS_PUB_LAYOUT", "single_column").strip().lower()
if layout not in ("single_column", "double_column"):
    raise RuntimeError("BS_PUB_LAYOUT must be 'single_column' or 'double_column'")
single_col_width_in = float(os.getenv("BS_PUB_SINGLE_COL_WIDTH_IN", "3.35"))
double_col_width_in = float(os.getenv("BS_PUB_DOUBLE_COL_WIDTH_IN", "6.9"))
fig_width_in = single_col_width_in if layout == "single_column" else double_col_width_in

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

def _auto_height(kind, n_panels=1, n_items=0):
    if layout == "single_column":
        if kind == "line_panels":
            return max(2.8, 1.55 * n_panels + 0.85)
        if kind == "heatmap_panels":
            return max(2.8, 1.4 * n_panels + 0.95)
        if kind == "quality":
            return 5.0
        if kind == "aggregate":
            return max(3.0, 0.34 * n_items + 1.35)
    else:
        if kind == "line_panels":
            return max(3.8, 1.15 * n_panels + 0.8)
        if kind == "heatmap_panels":
            return max(3.8, 1.05 * n_panels + 1.0)
        if kind == "quality":
            return 3.8
        if kind == "aggregate":
            return max(3.2, 0.28 * n_items + 1.3)
    return max(3.0, 1.5 * n_panels + 1.0)

def mk_axes_grid(nr, nc, kind):
    fig, axes = plt.subplots(
        nr,
        nc,
        figsize=(fig_width_in, _auto_height(kind, nr * nc)),
        squeeze=False,
        constrained_layout=True,
    )
    return fig, axes

# Figure 1: square runtime log-log (per family, faceted by threads)
panel_keys = [(fam, th) for fam in families for th in threads]
fig1, axes1 = mk_axes_grid(len(panel_keys), 1, "line_panels")
for idx, (fam, th) in enumerate(panel_keys):
    ax = axes1[idx][0]
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
    ax.set_title(f"{fam} | BLAS={th}")
    if idx == len(panel_keys) - 1:
        ax.set_xlabel("m=n")
    ax.set_ylabel("median time (s)")
    ax.grid(True, alpha=0.25)
axes1[0][0].legend(loc="best", frameon=True)
save_fig(fig1, "figure1_square_runtime")

# Figure 2: short-wide runtime vs n (fixed m curves, per family, faceted by threads)
fig2, axes2 = mk_axes_grid(len(panel_keys), 1, "line_panels")
global_short_ms = sorted({r["m"] for r in rows if r["regime"] == "short_wide"})
cmap = plt.get_cmap("tab10")
for idx, (fam, th) in enumerate(panel_keys):
    ax = axes2[idx][0]
    sub = [
        r for r in rows
        if r["family"] == fam and r["regime"] == "short_wide" and r["blas_threads"] == th
    ]
    ms = sorted({r["m"] for r in sub})
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
            ax.plot(ns, ymed, marker=marker, linestyle=ls, color=color)
            ax.fill_between(ns, ylow, yhigh, color=color, alpha=0.12)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_title(f"{fam} | BLAS={th}")
    if idx == len(panel_keys) - 1:
        ax.set_xlabel("n")
    ax.set_ylabel("median time (s)")
    ax.grid(True, alpha=0.25)
method_handles = [
    Line2D([0], [0], color="black", marker="o", linestyle="-", label="BSQR"),
    Line2D([0], [0], color="black", marker="D", linestyle="--", label="DGEQP3"),
]
m_handles = [
    Line2D([0], [0], color=cmap(i % 10), linestyle="-", marker=None, label=str(m))
    for i, m in enumerate(global_short_ms)
]
leg_method = fig2.legend(
    handles=method_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.018),
    ncol=2,
    frameon=True,
    title="Method",
)
leg_m = fig2.legend(
    handles=m_handles,
    loc="upper center",
    bbox_to_anchor=(0.5, 0.993),
    ncol=max(1, min(3, len(m_handles))),
    frameon=True,
    title="m (rows)",
)
fig2.add_artist(leg_method)
save_fig(fig2, "figure2_shortwide_runtime")

# Figure 3: short-wide speedup heatmaps over (m, n/m) per family/thread
fig3, axes3 = mk_axes_grid(len(panel_keys), 1, "heatmap_panels")
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
imh = None
for idx, (fam, th) in enumerate(panel_keys):
    ax = axes3[idx][0]
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
    ax.set_title(f"{fam} | BLAS={th}")
    ax.set_xticks(range(len(aspects)))
    ax.set_xticklabels([f"{a:.1f}" for a in aspects], rotation=45, ha="right")
    ax.set_yticks(range(len(ms)))
    ax.set_yticklabels([str(m) for m in ms])
    if idx == len(panel_keys) - 1:
        ax.set_xlabel("aspect (n/m)")
    ax.set_ylabel("m")
if imh is not None:
    cbar = fig3.colorbar(imh, ax=axes3.ravel().tolist(), shrink=0.8, pad=0.01)
    cbar.set_label("speedup dgeqp3/bsqr")
save_fig(fig3, "figure3_shortwide_speedup_heatmap")

# Figure 4: quality plot (residual and orthogonality, log-scale) across regimes
def grouped_quality(metric, regime, method):
    vals = [r[metric] for r in rows if r["regime"] == regime and r["method"] == method]
    if not vals:
        return float("nan"), float("nan"), float("nan")
    arr = np.asarray(vals, dtype=float)
    return float(np.median(arr)), float(np.percentile(arr, 10)), float(np.percentile(arr, 90))

if layout == "single_column":
    fig4, axes4 = plt.subplots(2, 1, figsize=(fig_width_in, _auto_height("quality")), squeeze=False, constrained_layout=True)
    quality_axes = [axes4[0][0], axes4[1][0]]
else:
    fig4, axes4 = plt.subplots(1, 2, figsize=(fig_width_in, _auto_height("quality")), squeeze=False, constrained_layout=True)
    quality_axes = [axes4[0][0], axes4[0][1]]
regime_labels = ["square", "short_wide"]
regime_display = ["square", "short-wide"]
x = np.arange(len(regime_labels))
barw = 0.34
for ax, metric, title in zip(quality_axes, ["residual", "orthogonality"], ["Residual", "Orthogonality"]):
    y_bs = []
    y_bs_lo = []
    y_bs_hi = []
    y_dg = []
    y_dg_lo = []
    y_dg_hi = []
    for regime in regime_labels:
        med, lo, hi = grouped_quality(metric, regime, "bsqr_full")
        y_bs.append(med)
        y_bs_lo.append(lo)
        y_bs_hi.append(hi)
        med, lo, hi = grouped_quality(metric, regime, "dgeqp3")
        y_dg.append(med)
        y_dg_lo.append(lo)
        y_dg_hi.append(hi)
    yerr_bs = np.vstack([np.asarray(y_bs) - np.asarray(y_bs_lo), np.asarray(y_bs_hi) - np.asarray(y_bs)])
    yerr_dg = np.vstack([np.asarray(y_dg) - np.asarray(y_dg_lo), np.asarray(y_dg_hi) - np.asarray(y_dg)])
    ax.bar(x - barw / 2, y_bs, width=barw, color="#4C78A8", label="bsqr_full")
    ax.bar(x + barw / 2, y_dg, width=barw, color="#F58518", label="dgeqp3")
    ax.errorbar(x - barw / 2, y_bs, yerr=yerr_bs, fmt="none", ecolor="black", capsize=3, linewidth=0.8)
    ax.errorbar(x + barw / 2, y_dg, yerr=yerr_dg, fmt="none", ecolor="black", capsize=3, linewidth=0.8)
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels(regime_display)
    ax.set_title(title)
    ax.grid(True, axis="y", alpha=0.25)
quality_axes[0].legend(loc="upper left")
save_fig(fig4, "figure4_quality")

# Figure 5: aggregate speedup bars with per-seed intervals
seed_group = defaultdict(list)
for s in speed_rows:
    key = (s["family"], s["regime"], s["blas_threads"], s["seed"])
    seed_group[key].append(s["speedup"])

combo_seed_vals = defaultdict(list)
for (fam, regime, th, seed), vals in seed_group.items():
    combo_seed_vals[(fam, regime, th)].append(geomean(vals))

def short_family_name(name):
    if name == "ill_conditioned":
        return "ill-cond"
    if name == "orthonormal_rows":
        return "orth-rows"
    return name

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
    regime_short = "sq" if key[1] == "square" else "sw"
    labels.append(f"{short_family_name(key[0])} | {regime_short} | t{key[2]}")
    centers.append(center)
    lows.append(lo)
    highs.append(hi)

fig5, ax5 = plt.subplots(
    1,
    1,
    figsize=(fig_width_in, _auto_height("aggregate", n_items=len(labels))),
    constrained_layout=True,
)
y = np.arange(len(labels))
xerr = np.vstack([np.asarray(centers) - np.asarray(lows), np.asarray(highs) - np.asarray(centers)])
ax5.barh(y, centers, color="#72B7B2")
ax5.errorbar(centers, y, xerr=xerr, fmt="none", ecolor="black", capsize=3, linewidth=0.8)
ax5.axvline(1.0, color="black", linestyle="--", linewidth=1.0)
ax5.set_yticks(y)
ax5.set_yticklabels(labels)
ax5.set_xlabel("geomean speedup (dgeqp3/bsqr)")
ax5.set_title("Aggregate speedup by family/regime/thread")
ax5.grid(True, axis="x", alpha=0.25)
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
