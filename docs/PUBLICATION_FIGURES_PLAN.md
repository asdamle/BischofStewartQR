# Publication Benchmark Figures: Specification

Executed 2026-06-12 (this file previously held the improvement plan; it is now the figure
spec both plotters conform to). Implementations: `julia/benchmark/plot_publication.py`
(invoked by the validating wrapper `julia/benchmark/plot_publication_results.jl`) and
`matlab/benchmark/plot_publication_results.m` (native MATLAB, no Python dependency). Any
change to figure semantics, labels, or styling must land in both. Everything here is
downstream of the timing CSVs — restyling never requires re-benchmarking.

## Figure inventory

Per comparison mode (`plain`: `bsqr_full` vs the built-in pivoted QR; `rinv`: `bsqr_rinv` vs
built-in + `trsm`), four figures plus a captions file are written to `plots/{plain,rinv}/`:

| file | content | size |
|---|---|---|
| `fig_square_runtime` | median runtime vs `m = n`, log–log; panels family × threads (Julia) / family (MATLAB); x-ticks at the measured sizes | double col (6.9 in) |
| `fig_shortwide_runtime` | median runtime vs `n`, log–log, one color per `m` (tab10/lines); power-of-two x-ticks; combined method + m legend below | double col |
| `fig_relative_time` | forest plot of geomean relative time per family × regime, parity line at 1; Julia: filled/open markers for 1/4 BLAS threads. Rows labelled `Family (m = n)` / `Family (m < n)` | single col (3.35 in) |
| `figure_captions.md` | suggested caption text with test-matrix and metric definitions plus run metadata | — |

One additional **composite** figure is written once to the top-level `plots/`
(`fig_relative_time_composite`): the `plain` (`BSQR / CPQR`) and `rinv`
(`(BSQR + R₁₁⁻¹R₁₂) / (CPQR + R₁₁⁻¹R₁₂)` — both methods also form the interpolation
matrix `R₁₁⁻¹R₁₂`) relative times overlaid on the shared `Family (m = n)` / `Family (m < n)`
rows, **colour = ratio** (plain `#4C78A8`, rinv `#F58518`). Julia additionally encodes
threads by marker fill (filled = 1, open = 4 BLAS threads). This is the single timing
figure intended for the paper; the per-mode `fig_relative_time` figures are retained.

Numerical quality is reported as numbers, not a figure: `tables/{mode}/quality_summary.md`
gives median and maximum relative residual and `‖I − QᵀQ‖_F` per method over all cases, with
a one-sentence "never exceeded" summary. The former quality figure/table and the short-wide
heatmap are intentionally absent (the heatmap's per-(m, aspect) resolution lives in
`table_shortwide_relative_time.csv`).

## Test-matrix families (identical construction in both languages)

1. **Gaussian** — i.i.d. standard normal entries (`randn(m, n)`).
2. **Ill-conditioned** — `A = U Σ Vᵀ` with `U` (m×r) and `V` (n×r), `r = min(m, n)`,
   orthonormal factors of Gaussian matrices, and `Σ` geometrically graded:
   `σ_i = exp(linspace(0, −log κ, r))` with `κ = 1e10`.
3. **Orthonormal rows** — `A = Qᵀ` where `Q` is an orthonormal basis (n×m, `m ≤ n`) of a
   Gaussian matrix, so `A Aᵀ = I_m` — the GKS column-selection setting the Bischof–Stewart
   criterion is designed for.

Square cases use `m = n ∈ {64, 128, 256, 384, 512}`; short-wide cases sweep
`m ∈ {32, 64, 128, 256, 512}` × aspect `n/m ∈ {2, 4, 8, 10}` (orthonormal-rows skips
`m > n`). Each (family, size, seed) matrix is regenerated deterministically from the seed.

## Statistics (identical in both languages)

- Runtime curves: points are **geomeans across seeds** of the per-case median time; the
  seed-range bands are drawn **first (behind every line) at alpha 0.15** so they read as
  faint context and never occlude data. (EPS does not support transparency and flattens the
  bands opaque, but the draw order still keeps them behind the lines; prefer PDF.)
- Relative time: per-seed geomean over cases, then geomean across seeds for the point and the
  per-seed range for the whiskers.
- Quality (in `quality_summary.md`): median and maximum per method over all cases. Metrics
  are Frobenius-norm quantities in both languages: `‖AΠ − QR‖_F / ‖A‖_F` and `‖I − QᵀQ‖_F`.

## Style

- Method display names: `BSQR` vs `CPQR (built-in)` (plain); `BSQR (W returned)` vs
  `CPQR + solve` (rinv). Captions name the exact routine per language; figures stay
  language-neutral. Ratio axis label: `relative time (BSQR / CPQR)` (resp. `BSQR+W /
  CPQR+solve`); the parity line is explained in the caption.
- Colors/markers: BSQR `#4C78A8`, solid line, filled circle; baseline `#F58518`, dashed,
  diamond. In `fig_shortwide_runtime` the color encodes `m` and the marker/line style the
  method. No in-figure titles beyond panel headers (`Gaussian — 1 thread` / `Gaussian`).
- Typography: serif text with TeX-rendered math — matplotlib `mathtext.fontset=cm` with a
  serif text font and `pdf.fonttype=ps.fonttype=42` (embedded, no Type-3); MATLAB uses the
  `latex` interpreter for all labels/ticks/legends (Computer Modern). Base size 8.5 pt,
  legends/ticks 7.5 pt, axis line width 0.7 pt.
- Output: `exportgraphics` (MATLAB) / `savefig` 300 dpi PNG + vector PDF/EPS. Formats via
  `BS_PUB_FIG_FORMATS` / `BS_MATLAB_PUB_FIG_FORMATS` (default `png`; committed artifacts are
  generated with `png,pdf,eps`). Known cosmetic issue: MATLAB's EPS export of
  latex-interpreted text warns `Font cmss8 is not supported` and substitutes; the PDFs are
  clean — prefer PDF for the paper's vector graphics.

## Tables (identical schemas; MATLAB omits `blas_threads`, which its CSV does not have)

- `table_square_relative_time.csv`: `family[, blas_threads], m, n, relative_time_geomean,
  relative_time_seed_min, relative_time_seed_max, bsqr_tmed_geomean_s,
  baseline_tmed_geomean_s`.
- `table_shortwide_relative_time.csv`: same plus `aspect`.
- `table_quality.csv`: `family[, blas_threads], regime, method, residual_median,
  residual_p95, orthogonality_median, orthogonality_p95`.

## Decisions adopted

1. Baseline display name `CPQR (built-in)` (language-neutral) rather than routine names.
2. Semantic filenames (`fig_*`) rather than renumbered `figureN_*`.
3. MATLAB plotter stays native; cohesion enforced by this spec rather than a shared renderer.
4. Intervals are per-seed ranges (honest with 2 seeds; definable identically in both
   languages — the MATLAB CSV has no CI columns).
