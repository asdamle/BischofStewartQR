#!/usr/bin/env python3
"""
Compare MATLAB and Julia BSQR timing rows from publication CSV files.

This script matches rows on (family, regime, m, n, seed) for one method,
then reports geomean timings for representative cases.
"""

import argparse
import csv
import math
from collections import defaultdict


def geomean(values):
    vals = [v for v in values if math.isfinite(v) and v > 0.0]
    if not vals:
        return float("nan")
    return math.exp(sum(math.log(v) for v in vals) / len(vals))


def parse_cases(raw):
    if not raw:
        return None
    out = []
    for tok in raw.split(","):
        tok = tok.strip()
        if not tok:
            continue
        regime, dims = tok.split(":", 1)
        m_s, n_s = dims.lower().split("x", 1)
        out.append((regime.strip(), int(m_s), int(n_s)))
    if not out:
        return None
    return out


def read_csv_rows(path):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError(f"No rows in {path}")
    return rows


def filter_rows(rows, method, blas_threads=None):
    out = []
    for r in rows:
        if r.get("method", "") != method:
            continue
        if blas_threads is not None:
            bt = r.get("blas_threads")
            if bt is None or int(bt) != blas_threads:
                continue
        out.append(r)
    return out


def aggregate_by_case(rows):
    case_map = defaultdict(list)
    for r in rows:
        key = (r["regime"], int(r["m"]), int(r["n"]))
        case_map[key].append(float(r["tmed_s"]))
    return {k: geomean(v) for (k, v) in case_map.items()}


def keyed_rows(rows):
    out = {}
    for r in rows:
        key = (r["family"], r["regime"], int(r["m"]), int(r["n"]), int(r["seed"]))
        out[key] = float(r["tmed_s"])
    return out


def main():
    ap = argparse.ArgumentParser(description="Compare MATLAB and Julia BSQR timing CSVs")
    ap.add_argument("--matlab-csv", required=True)
    ap.add_argument("--julia-csv", required=True)
    ap.add_argument("--method", default="bsqr_full", choices=["bsqr_full", "bsqr_rinv"])
    ap.add_argument("--julia-threads", type=int, default=1)
    ap.add_argument(
        "--cases",
        default="",
        help="Comma-separated regime:mxn list, e.g. square:256x256,short_wide:128x1024",
    )
    args = ap.parse_args()

    matlab_rows = read_csv_rows(args.matlab_csv)
    julia_rows = read_csv_rows(args.julia_csv)

    matlab = filter_rows(matlab_rows, method=args.method)
    julia = filter_rows(
        julia_rows,
        method=args.method,
        blas_threads=args.julia_threads,
    )

    if not matlab:
        raise RuntimeError("No MATLAB rows matched filters")
    if not julia:
        raise RuntimeError("No Julia rows matched filters")

    # Shared-key ratio statistics.
    mk = keyed_rows(matlab)
    jk = keyed_rows(julia)
    shared_keys = sorted(set(mk.keys()) & set(jk.keys()))
    if not shared_keys:
        raise RuntimeError("No shared (family, regime, m, n, seed) keys matched filters")
    ratios = [mk[k] / jk[k] for k in shared_keys if jk[k] > 0.0]
    ratio_geomean = geomean(ratios)
    ratio_median = sorted(ratios)[len(ratios) // 2]

    m_case = aggregate_by_case(matlab)
    j_case = aggregate_by_case(julia)
    cases = parse_cases(args.cases)
    if cases is None:
        shared_case_keys = sorted(set(m_case.keys()) & set(j_case.keys()), key=lambda x: (x[0], x[1], x[2]))
        cases = shared_case_keys[:6]

    print(f"# MATLAB vs Julia ({args.method}, julia_threads={args.julia_threads})")
    print("")
    print("| regime | m x n | MATLAB tmed (ms) | Julia tmed (ms) | MATLAB/Julia |")
    print("|---|---:|---:|---:|---:|")
    for regime, m, n in cases:
        km = (regime, m, n)
        mt = m_case.get(km, float("nan"))
        jt = j_case.get(km, float("nan"))
        rr = (mt / jt) if (jt and math.isfinite(jt) and jt > 0.0) else float("nan")
        mt_ms = mt * 1e3 if math.isfinite(mt) else float("nan")
        jt_ms = jt * 1e3 if math.isfinite(jt) else float("nan")
        print(f"| {regime} | {m}x{n} | {mt_ms:.3f} | {jt_ms:.3f} | {rr:.3f}x |")

    print("")
    print(f"Shared matched keys: {len(shared_keys)}")
    print(f"Overall MATLAB/Julia median ratio: {ratio_median:.3f}x")
    print(f"Overall MATLAB/Julia geomean ratio: {ratio_geomean:.3f}x")


if __name__ == "__main__":
    main()
