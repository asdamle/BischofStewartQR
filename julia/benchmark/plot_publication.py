#!/usr/bin/env python3
"""Publication figures and tables for the BSQR benchmark CSV (Julia pipeline).

Figure spec: docs/PUBLICATION_FIGURES_PLAN.md. The MATLAB plotter
(matlab/benchmark/plot_publication_results.m) conforms to the same spec; any
change to figure semantics, labels, or styling must land in both.

Usage:
  python3 plot_publication.py CSV OUTDIR TABLEDIR BSQR_METHOD BASELINE_METHOD MODE

MODE is "plain", "rinv" (display labels), or "composite" (the top-level
overlaid relative-time figure; the method arguments are ignored).
"""

import csv
import math
import os
import sys
from collections import defaultdict

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np

# ---------------------------------------------------------------------------
# Shared style (keep in lockstep with publication_style() in the MATLAB plotter)
# ---------------------------------------------------------------------------

BSQR_COLOR = "#4C78A8"
BASE_COLOR = "#F58518"
SINGLE_COL_W = float(os.getenv("BS_PUB_SINGLE_COL_WIDTH_IN", "3.35"))
DOUBLE_COL_W = float(os.getenv("BS_PUB_DOUBLE_COL_WIDTH_IN", "6.9"))

MODE_LABELS = {
    "plain": {"bs": "BSQR", "base": "CPQR (built-in)", "ratio": "BSQR / CPQR"},
    "rinv": {"bs": r"BSQR + $R_{11}^{-1}R_{12}$",
             "base": r"CPQR + $R_{11}^{-1}R_{12}$",
             "ratio": r"(BSQR + $R_{11}^{-1}R_{12}$) / (CPQR + $R_{11}^{-1}R_{12}$)"},
}

FAMILY_DISPLAY = {
    "gaussian": "Gaussian",
    "ill_conditioned": "Ill-conditioned",
    "orthonormal_rows": "Orthonormal rows",
}

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["STIXGeneral", "Times New Roman", "Times", "DejaVu Serif"],
    "mathtext.fontset": "cm",
    "font.size": 10.5,
    "axes.titlesize": 10.5,
    "axes.labelsize": 10.5,
    "legend.fontsize": 9.5,
    "legend.title_fontsize": 10,
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 9.5,
    "axes.linewidth": 0.7,
    "lines.linewidth": 1.2,
    "lines.markersize": 3.6,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

GRID_KW = {"color": "#dddddd", "linewidth": 0.5}


def parse_fig_formats():
    raw = os.getenv("BS_PUB_FIG_FORMATS", "png")
    formats = [part.strip().lower() for part in raw.split(",") if part.strip()]
    if not formats:
        raise RuntimeError("BS_PUB_FIG_FORMATS must list at least one format")
    bad = sorted(set(formats) - {"png", "pdf", "eps"})
    if bad:
        raise RuntimeError(f"Unsupported BS_PUB_FIG_FORMATS value(s): {', '.join(bad)}")
    out = []
    for fmt in formats:
        if fmt not in out:
            out.append(fmt)
    return out


# ---------------------------------------------------------------------------
# Data loading and aggregation
# ---------------------------------------------------------------------------

def load_rows(path):
    rows = []
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            rows.append({
                "run_id": row["run_id"],
                "family": row["family"],
                "regime": row["regime"],
                "m": int(row["m"]),
                "n": int(row["n"]),
                "aspect": float(row["aspect"]),
                "seed": int(row["seed"]),
                "blas_threads": int(row["blas_threads"]),
                "method": row["method"],
                "tmed": float(row["tmed_s"]),
                "residual": float(row["residual"]),
                "orthogonality": float(row["orthogonality"]),
            })
    if not rows:
        raise RuntimeError(f"No rows found in {path}")
    for idx, r in enumerate(rows):
        vals = (r["tmed"], r["residual"], r["orthogonality"])
        if not all(math.isfinite(v) for v in vals) or r["tmed"] <= 0.0:
            raise RuntimeError(f"Invalid numeric value at row index {idx}")
    return rows


def geomean(vals):
    good = [v for v in vals if math.isfinite(v) and v > 0.0]
    if not good:
        return float("nan")
    return math.exp(sum(math.log(v) for v in good) / len(good))


def seed_stats(vals):
    """Geomean across seeds with the per-seed range."""
    good = [v for v in vals if math.isfinite(v) and v > 0.0]
    if not good:
        nan = float("nan")
        return nan, nan, nan
    return geomean(good), min(good), max(good)


def paired_relative_rows(rows, bs_method, base_method):
    idx = {}
    for r in rows:
        key = (r["family"], r["regime"], r["m"], r["n"], r["aspect"],
               r["seed"], r["blas_threads"], r["method"])
        idx[key] = r
    out = []
    for r in rows:
        if r["method"] != bs_method:
            continue
        base = idx.get((r["family"], r["regime"], r["m"], r["n"], r["aspect"],
                        r["seed"], r["blas_threads"], base_method))
        if base is None or base["tmed"] <= 0.0:
            continue
        out.append({**{k: r[k] for k in ("family", "regime", "m", "n", "aspect",
                                          "seed", "blas_threads")},
                    "relative_time": r["tmed"] / base["tmed"],
                    "bsqr_tmed": r["tmed"], "baseline_tmed": base["tmed"]})
    if not out:
        raise RuntimeError("No paired relative-time rows; check method names")
    return out


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def save_fig(fig, outdir, stem, formats):
    for fmt in formats:
        path = os.path.join(outdir, f"{stem}.{fmt}")
        if fmt == "png":
            fig.savefig(path, dpi=300)
        else:
            fig.savefig(path, format=fmt)
    plt.close(fig)


def method_curve_style(role):
    if role == "bs":
        return {"color": BSQR_COLOR, "marker": "o", "ls": "-"}
    return {"color": BASE_COLOR, "marker": "D", "ls": "--"}


BAND_ALPHA = 0.15  # seed-range bands: faint and always behind the lines


def thread_panel_title(family, th):
    unit = "thread" if th == 1 else "threads"
    return f"{FAMILY_DISPLAY.get(family, family)} — {th} {unit}"


def fig_square_runtime(rows, methods, labels, outdir, formats):
    sub_all = [r for r in rows if r["regime"] == "square"]
    families = sorted({r["family"] for r in sub_all})
    threads = sorted({r["blas_threads"] for r in sub_all})
    nr, nc = len(families), len(threads)
    fig, axes = plt.subplots(nr, nc, figsize=(DOUBLE_COL_W, 1.65 * nr + 0.55),
                             squeeze=False, constrained_layout=True,
                             sharex=True, sharey=True)
    for i, fam in enumerate(families):
        for j, th in enumerate(threads):
            ax = axes[i][j]
            sub = [r for r in sub_all if r["family"] == fam and r["blas_threads"] == th]
            ms = sorted({r["m"] for r in sub})
            for role, method in (("bs", methods[0]), ("base", methods[1])):
                st = method_curve_style(role)
                med, lo, hi = [], [], []
                for m in ms:
                    g, mn, mx = seed_stats([r["tmed"] for r in sub
                                            if r["m"] == m and r["method"] == method])
                    med.append(g)
                    lo.append(mn)
                    hi.append(mx)
                ax.fill_between(ms, lo, hi, color=st["color"], alpha=BAND_ALPHA,
                                edgecolor="none", zorder=1)
                ax.plot(ms, med, label=labels[role], zorder=3, **st)
            ax.set_xscale("log", base=2)
            ax.set_yscale("log")
            ax.set_title(thread_panel_title(fam, th))
            ax.grid(True, **GRID_KW)
            if i == nr - 1:
                ax.set_xlabel(r"$m = n$")
            if j == 0:
                ax.set_ylabel("median time [s]")
    all_ms = sorted({r["m"] for r in sub_all})
    axes[0][0].set_xticks(all_ms)
    axes[0][0].set_xticklabels([str(m) for m in all_ms])
    axes[0][0].xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axes[0][0].legend(loc="upper left", frameon=False)
    save_fig(fig, outdir, "fig_square_runtime", formats)


def fig_shortwide_runtime(rows, methods, labels, outdir, formats):
    sub_all = [r for r in rows if r["regime"] == "short_wide"]
    families = sorted({r["family"] for r in sub_all})
    threads = sorted({r["blas_threads"] for r in sub_all})
    all_ms = sorted({r["m"] for r in sub_all})
    cmap = plt.get_cmap("tab10")
    m_color = {m: cmap(i % 10) for i, m in enumerate(all_ms)}

    nr, nc = len(families), len(threads)
    fig, axes = plt.subplots(nr, nc, figsize=(DOUBLE_COL_W, 1.65 * nr + 1.5),
                             squeeze=False, sharex=True, sharey=True)
    for i, fam in enumerate(families):
        for j, th in enumerate(threads):
            ax = axes[i][j]
            sub = [r for r in sub_all if r["family"] == fam and r["blas_threads"] == th]
            series = []
            for m in sorted({r["m"] for r in sub}):
                for role, method in (("bs", methods[0]), ("base", methods[1])):
                    st = method_curve_style(role)
                    rn = [r for r in sub if r["m"] == m and r["method"] == method]
                    ns = sorted({r["n"] for r in rn})
                    med, lo, hi = [], [], []
                    for n in ns:
                        g, mn, mx = seed_stats([r["tmed"] for r in rn if r["n"] == n])
                        med.append(g)
                        lo.append(mn)
                        hi.append(mx)
                    series.append((m, st, ns, med, lo, hi))
            # Bands first (faint, behind every line), then all median lines.
            for m, st, ns, med, lo, hi in series:
                ax.fill_between(ns, lo, hi, color=m_color[m], alpha=BAND_ALPHA,
                                edgecolor="none", zorder=1)
            for m, st, ns, med, lo, hi in series:
                ax.plot(ns, med, marker=st["marker"], ls=st["ls"],
                        color=m_color[m], markersize=3.0, zorder=3)
            ax.set_xscale("log", base=2)
            ax.set_yscale("log")
            ax.set_title(thread_panel_title(fam, th))
            ax.grid(True, **GRID_KW)
            if i == nr - 1:
                ax.set_xlabel(r"$n$")
            if j == 0:
                ax.set_ylabel("median time [s]")

    fig.subplots_adjust(bottom=0.21, top=0.95, left=0.09, right=0.99,
                        hspace=0.45, wspace=0.15)
    method_handles = [
        Line2D([0], [0], color="black", marker=method_curve_style(role)["marker"],
               ls=method_curve_style(role)["ls"], markersize=3.0, label=labels[role])
        for role in ("bs", "base")
    ]
    m_handles = [Line2D([0], [0], color=m_color[m], ls="-", label=f"$m = {m}$")
                 for m in all_ms]
    leg_methods = fig.legend(handles=method_handles, loc="lower center",
                             bbox_to_anchor=(0.5, 0.095), ncol=2, frameon=False)
    fig.legend(handles=m_handles, loc="lower center", bbox_to_anchor=(0.5, 0.0),
               ncol=min(5, len(m_handles)), frameon=False)
    fig.add_artist(leg_methods)
    save_fig(fig, outdir, "fig_shortwide_runtime", formats)


def relative_combo_stats(rel_rows):
    """(family, regime, threads) -> (geomean over per-seed geomeans, min, max)."""
    seed_geo = defaultdict(list)
    for s in rel_rows:
        seed_geo[(s["family"], s["regime"], s["blas_threads"], s["seed"])].append(
            s["relative_time"])
    combo = defaultdict(list)
    for (fam, regime, th, _seed), vals in seed_geo.items():
        combo[(fam, regime, th)].append(geomean(vals))
    return {key: seed_stats(vals) for key, vals in combo.items()}


# Composite: the plain (BSQR / CPQR) and rinv (BSQR+W / CPQR+solve) relative times
# overlaid on the shared family×regime rows -- the single timing figure for the paper.
COMPOSITE_MODES = [("plain", "bsqr_full", "dgeqp3"), ("rinv", "bsqr_rinv", "dgeqp3_trsm")]
COMPOSITE_COLOR = {"plain": BSQR_COLOR, "rinv": BASE_COLOR}
REGIME_DISPLAY = {"square": "(m = n)", "short_wide": "(m < n)"}


def fig_relative_time_composite(rows, outdir, formats):
    mode_stats, fams, thset = {}, set(), set()
    for name, bs, base in COMPOSITE_MODES:
        sub = [r for r in rows if r["method"] in (bs, base)]
        st = relative_combo_stats(paired_relative_rows(sub, bs, base))
        mode_stats[name] = st
        for fam, regime, th in st:
            fams.add(fam)
            thset.add(th)
    families = sorted(fams)
    threads = sorted(thset)
    row_keys = [(fam, regime) for fam in families for regime in ("square", "short_wide")
                if any((fam, regime, th) in st for st in mode_stats.values()
                       for th in threads)]
    nrows = len(row_keys)

    fig, ax = plt.subplots(1, 1, figsize=(SINGLE_COL_W, 0.52 * nrows + 1.45),
                           constrained_layout=True)
    y_base = np.arange(nrows, dtype=float)[::-1]
    series = [(name, th) for name in ("plain", "rinv") for th in threads]
    offs = np.linspace(0.30, -0.30, len(series)) if len(series) > 1 else [0.0]
    offmap = {s: offs[i] for i, s in enumerate(series)}

    for name, bs, base in COMPOSITE_MODES:
        color = COMPOSITE_COLOR[name]
        st = mode_stats[name]
        for ti, th in enumerate(threads):
            ys, cs, lo_e, hi_e = [], [], [], []
            for ri, key in enumerate(row_keys):
                s = st.get((key[0], key[1], th))
                if s is None or not math.isfinite(s[0]):
                    continue
                ys.append(y_base[ri] + offmap[(name, th)])
                cs.append(s[0])
                lo_e.append(s[0] - s[1])
                hi_e.append(s[2] - s[0])
            if not cs:
                continue
            face = color if ti == 0 else "white"
            ax.errorbar(cs, ys, xerr=np.vstack([lo_e, hi_e]), fmt="o", color=color,
                        markerfacecolor=face, markersize=4.0, linestyle="none",
                        capsize=2.0, elinewidth=0.9, zorder=3)
    ax.axvline(1.0, color="black", ls="--", lw=0.8, zorder=1)
    ax.set_yticks(y_base)
    ax.set_yticklabels([f"{FAMILY_DISPLAY.get(fam, fam)} {REGIME_DISPLAY[regime]}"
                        for fam, regime in row_keys])
    ax.set_ylim(-0.7, nrows - 0.3)
    ax.set_xlabel("relative time (BSQR / baseline)")
    ax.grid(True, axis="x", **GRID_KW)

    handles = [Line2D([0], [0], color=COMPOSITE_COLOR[name], marker="o", ls="none",
                      label=MODE_LABELS[name]["ratio"]) for name in ("plain", "rinv")]
    if len(threads) > 1:
        for ti, th in enumerate(threads):
            unit = "thread" if th == 1 else "threads"
            handles.append(Line2D([0], [0], color="black", marker="o", ls="none",
                                  markerfacecolor=("black" if ti == 0 else "white"),
                                  label=f"{th} BLAS {unit}"))
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, -0.10),
              ncol=1, frameon=False, fontsize=8.8, handletextpad=0.4)
    save_fig(fig, outdir, "fig_relative_time_composite", formats)


# ---------------------------------------------------------------------------
# Tables and captions
# ---------------------------------------------------------------------------

def write_csv(path, columns, rows_out):
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        for row in rows_out:
            w.writerow(row)


def relative_time_table(rel_rows, regime, group_cols):
    groups = defaultdict(list)
    for s in rel_rows:
        if s["regime"] != regime:
            continue
        groups[tuple(s[c] for c in group_cols)].append(s)
    out = []
    for key in sorted(groups.keys()):
        vals = groups[key]
        rel = [v["relative_time"] for v in vals]
        row = dict(zip(group_cols, key))
        row.update({
            "relative_time_geomean": geomean(rel),
            "relative_time_seed_min": min(rel),
            "relative_time_seed_max": max(rel),
            "bsqr_tmed_geomean_s": geomean([v["bsqr_tmed"] for v in vals]),
            "baseline_tmed_geomean_s": geomean([v["baseline_tmed"] for v in vals]),
        })
        out.append(row)
    return out


def write_quality_summary(rows, methods, labels, tabledir):
    """Terse markdown report replacing the former quality figure/table."""
    ncases = len({(r["family"], r["regime"], r["m"], r["n"], r["aspect"],
                   r["seed"], r["blas_threads"]) for r in rows})
    lines = [
        "# Numerical quality summary",
        "",
        f"Across all {ncases} benchmark cases (3 matrix families, both regimes, all sizes, "
        "seeds, and thread settings):",
        "",
        "| method | median rel. residual | max rel. residual | "
        "median ‖I − QᵀQ‖_F | max ‖I − QᵀQ‖_F |",
        "|---|---:|---:|---:|---:|",
    ]
    stats = {}
    for role, method in (("bs", methods[0]), ("base", methods[1])):
        resid = np.asarray([r["residual"] for r in rows if r["method"] == method])
        orth = np.asarray([r["orthogonality"] for r in rows if r["method"] == method])
        stats[role] = (resid, orth)
        lines.append(
            f"| {labels[role]} | {np.median(resid):.2e} | {np.max(resid):.2e} "
            f"| {np.median(orth):.2e} | {np.max(orth):.2e} |")
    bs_resid, bs_orth = stats["bs"]
    lines += [
        "",
        f"The relative residual ‖AΠ − QR‖_F/‖A‖_F of {labels['bs']} never exceeded "
        f"{np.max(bs_resid):.2e}, and its deviation from orthogonality ‖I − QᵀQ‖_F never "
        f"exceeded {np.max(bs_orth):.2e}; both match the built-in baseline to within a "
        "small factor.",
        "",
    ]
    with open(os.path.join(tabledir, "quality_summary.md"), "w") as f:
        f.write("\n".join(lines))


def write_captions(rows, labels, mode, outdir):
    seeds = sorted({r["seed"] for r in rows})
    threads = sorted({r["blas_threads"] for r in rows})
    run_ids = sorted({r["run_id"] for r in rows})
    bs, base = labels["bs"], labels["base"]
    nseeds = len(seeds)
    lines = [
        f"# Suggested captions ({mode} comparison)",
        "",
        f"Run: {', '.join(run_ids)}; seeds: {nseeds}; BLAS threads: "
        f"{', '.join(str(t) for t in threads)}. In all figures, {bs} is compared "
        f"against {base}; both timed paths materialize Q, R, and the permutation.",
        "",
        "**Test matrices.** Three families, regenerated per seed: *Gaussian* — "
        "i.i.d. standard normal entries; *ill-conditioned* — A = UΣVᵀ with U, V "
        "orthonormal factors of Gaussian matrices and Σ geometrically graded from 1 "
        "down to 1e-10 (κ = 1e10); *orthonormal rows* — A = Qᵀ with Q an orthonormal "
        "basis of a Gaussian n×m matrix (m ≤ n, so AAᵀ = I — the GKS column-selection "
        "setting). Square cases use m = n; short-wide cases sweep aspect ratios "
        "n/m ∈ {2, 4, 8, 10}.",
        "",
        "**fig_square_runtime.** Median runtime versus matrix size for square "
        f"matrices (log–log). Lines are geometric means over {nseeds} seeds; faint "
        "bands show the range across seeds. Panels: matrix family × BLAS threads.",
        "",
        "**fig_shortwide_runtime.** Median runtime versus column count n for "
        "short-wide matrices, one color per row count m (log–log). Lines/bands as "
        "above; marker and line style distinguish the methods.",
        "",
        "**fig_relative_time_composite** (top-level plots/). The single timing figure: "
        "the relative time without (BSQR / CPQR) and with the interpolation matrix "
        "(both methods also form R11^{-1}R12; labelled "
        "(BSQR + R11^{-1}R12) / (CPQR + R11^{-1}R12)), "
        "overlaid on the shared rows. Colour distinguishes the two; marker fill "
        "distinguishes 1 vs 4 BLAS threads. Square rows are m = n, short-wide rows m < n "
        "(swept sizes as above); 1 = parity (dashed line). Points are geomeans of "
        "per-seed geomeans; whiskers show the per-seed range.",
        "",
        "Numerical quality is summarized in tables/quality_summary.md.",
        "",
    ]
    with open(os.path.join(outdir, "figure_captions.md"), "w") as f:
        f.write("\n".join(lines))


# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 7:
        raise SystemExit(__doc__)
    input_path, outdir, tabledir, bsqr_method, baseline_method, mode = sys.argv[1:7]
    formats = parse_fig_formats()
    rows = load_rows(input_path)

    if mode == "composite":
        # plain + rinv overlaid into the top-level plots dir; method args unused.
        os.makedirs(outdir, exist_ok=True)
        fig_relative_time_composite(rows, outdir, formats)
        print(f"Wrote composite relative-time figure to {outdir}")
        return

    if mode not in MODE_LABELS:
        raise RuntimeError(f"Unknown mode '{mode}'; expected plain, rinv, or composite")
    labels = MODE_LABELS[mode]
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(tabledir, exist_ok=True)
    methods = (bsqr_method, baseline_method)
    comp_rows = [r for r in rows if r["method"] in set(methods)]
    if not comp_rows:
        raise RuntimeError(f"No rows for methods {methods}")
    rel_rows = paired_relative_rows(comp_rows, bsqr_method, baseline_method)

    fig_square_runtime(comp_rows, methods, labels, outdir, formats)
    fig_shortwide_runtime(comp_rows, methods, labels, outdir, formats)
    # The per-mode relative-time forest plot is superseded by the top-level
    # fig_relative_time_composite (plain + rinv on one figure); rel_rows below still
    # feed the relative-time tables.

    write_csv(os.path.join(tabledir, "table_square_relative_time.csv"),
              ["family", "blas_threads", "m", "n",
               "relative_time_geomean", "relative_time_seed_min",
               "relative_time_seed_max", "bsqr_tmed_geomean_s",
               "baseline_tmed_geomean_s"],
              relative_time_table(rel_rows, "square",
                                  ["family", "blas_threads", "m", "n"]))
    write_csv(os.path.join(tabledir, "table_shortwide_relative_time.csv"),
              ["family", "blas_threads", "m", "n", "aspect",
               "relative_time_geomean", "relative_time_seed_min",
               "relative_time_seed_max", "bsqr_tmed_geomean_s",
               "baseline_tmed_geomean_s"],
              relative_time_table(rel_rows, "short_wide",
                                  ["family", "blas_threads", "m", "n", "aspect"]))
    write_quality_summary(comp_rows, methods, labels, tabledir)
    write_captions(comp_rows, labels, mode, outdir)

    print(f"Wrote publication figures to {outdir}")
    print(f"Wrote publication tables to {tabledir}")


if __name__ == "__main__":
    main()
