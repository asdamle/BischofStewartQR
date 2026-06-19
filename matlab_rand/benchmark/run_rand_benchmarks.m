function results = run_rand_benchmarks(varargin)
%RUN_RAND_BENCHMARKS Compare randomized vs deterministic BSQR selection.
%
%   Measures, on short-wide orthonormal-row matrices (the GKS regime where
%   randomization should help most):
%     * t_rand  - bsqr_rand returning only [p, reflectors, R11] (no R12)
%     * t_det   - deterministic bsqr returning R (selection + full factor;
%                 it must maintain R11^{-1}R12 for every column, which is the
%                 O(n k^2) cost the randomized variant avoids)
%   and reports speedup = t_det / t_rand alongside the two quality metrics:
%     * ratio        = ||R11^{-1}||_F / sqrt(k(n-k+1))   (Osinsky bound = 1)
%     * tested_per_k = total candidate columns evaluated / k
%
% Name-value options:
%   'sizes'   - cell array of [k n] pairs (default short-wide grid)
%   'seed'    - base RNG seed (default 1)
%   'block_size' - candidates per round (default 16)
%   'outdir'  - directory for the CSV (default matlab_rand/benchmark/results)

ip = inputParser;
addParameter(ip, 'sizes', default_sizes());
addParameter(ip, 'seed', 1);
addParameter(ip, 'block_size', 16);
addParameter(ip, 'outdir', '');
parse(ip, varargin{:});
opt = ip.Results;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab_rand'));
addpath(fullfile(repo_root, 'matlab'));
if isempty(opt.outdir)
    opt.outdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'results');
end
if ~isfolder(opt.outdir)
    mkdir(opt.outdir);
end

if ~bsqr_rand_mex_available()
    error('run_rand_benchmarks:NoMex', ...
        'bsqr_rand_mex not built. Run matlab_rand/build_bsqr_rand_mex.m first.');
end
det_backend = 'mfile';
if exist('bsqr_mex_available', 'file') && bsqr_mex_available()
    det_backend = 'mex';
end

modes = {'running_mean', 'worstcase_allowance'};
rows = {};
fprintf('%-14s %-20s %10s %10s %9s %8s %9s\n', ...
    'size', 'mode', 't_rand(ms)', 't_det(ms)', 'speedup', 'ratio', 'tested/k');

for ci = 1:numel(opt.sizes)
    k = opt.sizes{ci}(1);
    n = opt.sizes{ci}(2);
    rng(opt.seed + ci);
    M = orth(randn(n, k))';

    % deterministic selection cost: factor only (R, k-by-n), no Q.
    t_det = timeit(@() bsqr(M, 'k', k, 'backend', det_backend), 1);

    for mi = 1:numel(modes)
        mode = modes{mi};
        f = @() bsqr_rand(M, 'k', k, 'backend', 'mex', 'block_size', opt.block_size, ...
            'threshold_mode', mode, 'seed', opt.seed + ci);
        t_rand = timeit(f, 3);  % time only [p, reflectors, R11]
        [~, ~, ~, st] = bsqr_rand(M, 'k', k, 'backend', 'mex', 'block_size', opt.block_size, ...
            'threshold_mode', mode, 'seed', opt.seed + ci);
        ratio = st.frob_inv / st.osinsky_bound;
        tested_per_k = st.total_tested / k;
        speedup = t_det / t_rand;

        fprintf('%-14s %-20s %10.3f %10.3f %9.2f %8.3f %9.1f\n', ...
            sprintf('%dx%d', k, n), mode, t_rand * 1e3, t_det * 1e3, ...
            speedup, ratio, tested_per_k);
        rows(end+1, :) = {k, n, mode, t_rand, t_det, speedup, ratio, tested_per_k}; %#ok<AGROW>
    end
end

results = cell2table(rows, 'VariableNames', ...
    {'k', 'n', 'mode', 't_rand_s', 't_det_s', 'speedup', 'ratio', 'tested_per_k'});
csv = fullfile(opt.outdir, 'rand_timings.csv');
writetable(results, csv);
fprintf('\nWrote %s\n', csv);
end

function sizes = default_sizes()
sizes = {[32, 2000], [64, 4000], [64, 8000], [128, 8000], [32, 16000], [128, 16000]};
end
