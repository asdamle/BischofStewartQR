# Manuscript Figure Map

Translation table between manuscript figure numbers and the figure files in
this repository, so figure refinements can be requested by number. All
manuscript figures except Figure 8 come from the randomized comparison suite
and live in `matlab_rand/benchmark/plots/`; Figure 8 comes from the Julia
publication pipeline. See `docs/PLOT_PROVENANCE.md` for the data and plot
scripts that regenerate each one.

| Manuscript figure | Figure file | Location |
|---|---|---|
| Figure 1 | `fig_blocksize_k256` | `matlab_rand/benchmark/plots/` |
| Figure 2 | `fig_sampling_k256` | `matlab_rand/benchmark/plots/` |
| Figure 3 | `fig_scaling_speedup_k256` | `matlab_rand/benchmark/plots/` |
| Figure 4 | `fig_scaling_quality_k256` | `matlab_rand/benchmark/plots/` |
| Figure 5 | `fig_rpqr_quality_k256` | `matlab_rand/benchmark/plots/` |
| Figure 6 | `fig_largen_scaling` | `matlab_rand/benchmark/plots/` |
| Figure 7 | `fig_approx_cond_quality_spec` | `matlab_rand/benchmark/plots/` |
| Figure 8 | `fig_relative_time_composite` (Julia) | `julia/benchmark/results/publication/plots/` |

Notes:

- `fig_rpqr_speedup_k256` (the former Figure 6) was dropped from the
  manuscript: its message -- a flat ~2x speedup of randBSQR over
  rejection_rpqr across n and families -- is visible as the line gap in
  Figure 6 (`fig_largen_scaling`) and is cited in the text instead. The
  figure is still generated and committed alongside the others.
- The `_k256` figures are the k = 256 members of families also rendered at
  k = 64 and k = 128; a refinement to a manuscript figure should usually be
  propagated to the other-k variants of the same family.
- Figure 7 is the spectral-norm variant of `fig_approx_cond_quality`
  (rendered with `plot_approx_cond_comparison('norm','2')`); the Frobenius
  variant shares its layout, so refinements usually apply to both.
- Figure 8 is the Julia relative-time composite; the MATLAB pipeline renders
  the analogous `fig_relative_time_composite` in
  `matlab/benchmark/results/publication/plots/`, and layout changes should
  generally be mirrored there.
