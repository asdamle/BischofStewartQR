function run_approx_synth_comparison(varargin)
%RUN_APPROX_SYNTH_COMPARISON Synthetic-spectrum companion to the approximation
%   comparison. Runs run_approx_comparison on matrices with a prescribed
%   "interesting" spectrum and the rejection_rpqr leverage families as right
%   singular vectors (gaussian / spiked_leverage / needle; see approx_synth_matrix),
%   so one can read off how much the R11-conditioning differences on those families
%   translate into low-rank approximation error.
%
%   Writes results/exp_approx_synth.csv; plot with
%   plot_approx_comparison('tag','_synth'). Accepts the same name/value options as
%   run_approx_comparison ('ks', 'trials', 'rel_floor', 'sizes', 'outdir', ...).

run_approx_comparison('families', {'gaussian', 'spiked_leverage', 'needle'}, ...
    'tag', '_synth', varargin{:});
end
