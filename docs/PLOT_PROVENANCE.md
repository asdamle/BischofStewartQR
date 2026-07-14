# Plot Provenance

Where every figure in the repository comes from: for each plot (or family of
plots), the script that produces the underlying data, the data file it writes,
and the script that renders the figure. Rerun the data script, then the plot
script, to regenerate any figure from scratch; rerunning only the plot script
re-renders from the existing data file. See CLAUDE.md for the exact
command-line invocations, and `docs/MANUSCRIPT_FIGURES.md` for which figures
appear in the manuscript under which figure number.

## Julia publication pipeline

Figures land in `julia/benchmark/results/publication/plots/` (committed);
`plain/` is the standard factorization, `rinv/` the `R11^{-1}R12`-returning
variant. All scripts live in `julia/benchmark/` and run in that environment.

| Figure | Data produced by | Data file | Plot produced by |
|---|---|---|---|
| `plain/fig_square_runtime` | `run_publication_benchmarks.jl` | `publication_timings.csv` | `plot_publication_results.jl` |
| `plain/fig_shortwide_runtime` | `run_publication_benchmarks.jl` | `publication_timings.csv` | `plot_publication_results.jl` |
| `rinv/fig_square_runtime` | `run_publication_benchmarks.jl` | `publication_timings.csv` | `plot_publication_results.jl` |
| `rinv/fig_shortwide_runtime` | `run_publication_benchmarks.jl` | `publication_timings.csv` | `plot_publication_results.jl` |
| `fig_relative_time_composite` | `run_publication_benchmarks.jl` | `publication_timings.csv` | `plot_publication_results.jl` |

`plot_publication_results.jl` also writes the derived tables
(`tables/{plain,rinv}/table_{square,shortwide}_relative_time.csv`,
`quality_summary.md`) and per-figure caption files from the same CSV.
`run_publication_smoke_benchmark.jl` is a fast, gitignored-output variant of
the data runner for iteration; it feeds the same plot script.

## MATLAB publication pipeline

Same structure as the Julia pipeline; scripts live in `matlab/benchmark/`,
figures in `matlab/benchmark/results/publication/plots/` (committed).

| Figure | Data produced by | Data file | Plot produced by |
|---|---|---|---|
| `plain/fig_square_runtime` | `run_publication_benchmarks.m` | `publication_timings.csv` | `plot_publication_results.m` |
| `plain/fig_shortwide_runtime` | `run_publication_benchmarks.m` | `publication_timings.csv` | `plot_publication_results.m` |
| `rinv/fig_square_runtime` | `run_publication_benchmarks.m` | `publication_timings.csv` | `plot_publication_results.m` |
| `rinv/fig_shortwide_runtime` | `run_publication_benchmarks.m` | `publication_timings.csv` | `plot_publication_results.m` |
| `fig_relative_time_composite` | `run_publication_benchmarks.m` | `publication_timings.csv` | `plot_publication_results.m` |

Derived tables and captions: as in the Julia pipeline, written by
`plot_publication_results.m`. `run_publication_smoke_benchmark.m` is the fast
iteration variant.

## Randomized BSQR experiments (`matlab_rand/benchmark/`)

Figures land in `matlab_rand/benchmark/plots/` (committed); data CSVs in the
gitignored `results/`. The `_k<K>` families exist for K = 64, 128, 256 — one
figure per K, all from the K-tagged CSV of the same runner (pass `'k', K` to
both scripts).

| Figure (family) | Data produced by | Data file | Plot produced by |
|---|---|---|---|
| `fig_scaling_time_k<K>` | `run_rand_experiments.m` | `exp_scaling_k<K>.csv` | `plot_rand_experiments.m` |
| `fig_scaling_speedup_k<K>` | `run_rand_experiments.m` | `exp_scaling_k<K>.csv` | `plot_rand_experiments.m` |
| `fig_scaling_quality_k<K>` | `run_rand_experiments.m` | `exp_scaling_k<K>.csv` | `plot_rand_experiments.m` |
| `fig_blocksize_k<K>` | `run_rand_experiments.m` | `exp_blocksize_k<K>.csv` | `plot_rand_experiments.m` |
| `fig_sampling_k<K>` | `run_rand_experiments.m` | `exp_sampling_k<K>.csv` | `plot_rand_experiments.m` |
| `fig_rpqr_time_k<K>` | `run_rpqr_comparison.m` | `exp_rpqr_k<K>.csv` | `plot_rpqr_comparison.m` |
| `fig_rpqr_speedup_k<K>` | `run_rpqr_comparison.m` | `exp_rpqr_k<K>.csv` | `plot_rpqr_comparison.m` |
| `fig_rpqr_quality_k<K>` | `run_rpqr_comparison.m` | `exp_rpqr_k<K>.csv` | `plot_rpqr_comparison.m` |
| `fig_largen_scaling` | `run_largen_scaling.m` | `exp_largen.csv` | `plot_largen_scaling.m` |
| `fig_approx_quality` | `run_approx_comparison.m` | `exp_approx.csv` | `plot_approx_comparison.m` |
| `fig_approx_synth_quality` | `run_approx_synth_comparison.m` (delegates to `run_approx_comparison.m`) | `exp_approx_synth.csv` | `plot_approx_comparison.m` (`'tag','_synth'`) |
| `fig_approx_cond_quality` | `run_approx_cond_comparison.m` | `exp_approx_cond.csv` | `plot_approx_cond_comparison.m` |
| `fig_approx_cond_quality_spec` | `run_approx_cond_comparison.m` | `exp_approx_cond.csv` | `plot_approx_cond_comparison.m` (`'norm','2'`) |

Notes:
- The `rpqr`, `largen`, and `approx*` runners require the third-party
  `rejection_rpqr` code under `ext_comparisons/` (not redistributed; see
  `matlab_rand/README.md` for the download step).
- `run_rand_benchmarks.m` is a data/report-only companion (writes
  `rand_timings.csv`); it produces no figure.
