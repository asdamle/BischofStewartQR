function results = run_rand_benchmarks(varargin)
%RUN_RAND_BENCHMARKS Fair randomized-vs-deterministic BSQR timing (no R12).
%
%   Fairness/measurement discipline (both methods, identical conditions):
%     * The compiled MEX kernels are called DIRECTLY (bsqr_mex, bsqr_rand_mex),
%       so the m-file dispatcher and its inputParser are never timed.
%     * 'check_finite', false on both -- the O(m*n) finiteness scan is identical
%       overhead for both and is not part of either algorithm.
%     * timeit() handles warm-up and returns a robust median; only the kernel
%       call is inside the timed thunk (matrix generation is outside).
%     * Deterministic baseline: bsqr_mex(M,'k',k) with ONE output (R only).
%       That is the cheapest deterministic call that still performs the column
%       selection; it forms R12 as an unavoidable byproduct (the O(n*k^2) work
%       the randomized variant skips) but does NOT materialize Q.
%     * Randomized: bsqr_rand_mex(...) with THREE outputs [p, reflectors, R11]
%       -- the "R12 not needed" product. No Q, no R12.
%
% Name-value options: 'sizes' (cell of [k n]), 'seed', 'block_size',
%   'family' (see rand_test_matrix), 'outdir'.

ip = inputParser;
addParameter(ip, 'sizes', {[32, 2000], [64, 4000], [64, 8000], [128, 8000], [128, 16000]});
addParameter(ip, 'seed', 1);
addParameter(ip, 'block_size', 16);
addParameter(ip, 'family', 'gaussian');
addParameter(ip, 'outdir', '');
parse(ip, varargin{:});
opt = ip.Results;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab_rand'));
addpath(fullfile(repo_root, 'matlab_rand', 'mex'));
addpath(fullfile(repo_root, 'matlab_rand', 'benchmark'));
addpath(fullfile(repo_root, 'matlab'));
addpath(fullfile(repo_root, 'matlab', 'mex'));
if isempty(opt.outdir)
    opt.outdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'results');
end
if ~isfolder(opt.outdir); mkdir(opt.outdir); end

assert(exist('bsqr_rand_mex', 'file') == 3, ...
    'bsqr_rand_mex not built. Run matlab_rand/build_bsqr_rand_mex.m first.');
assert(exist('bsqr_mex', 'file') == 3, ...
    'bsqr_mex (deterministic) not built. Run matlab/build_bsqr_mex.m first.');

modes = {'running_mean', 'worstcase_allowance'};
rows = {};
fprintf('family=%s  block_size=%d  BLAS threads=%d\n', opt.family, opt.block_size, maxNumCompThreads);
fprintf('%-12s %-20s %11s %11s %8s %8s %9s\n', ...
    'size', 'mode', 't_rand(ms)', 't_det(ms)', 'speedup', 'cond', 'tested/k');

for ci = 1:numel(opt.sizes)
    k = opt.sizes{ci}(1); n = opt.sizes{ci}(2);
    M = rand_test_matrix(opt.family, k, n, opt.seed + ci);

    t_det = timeit(@() bsqr_mex(M, 'k', k, 'check_finite', false), 1);
    Rdet = bsqr_mex(M, 'k', k, 'check_finite', false);
    frobinv_det = norm(inv(triu(Rdet(1:k, 1:k))), 'fro');

    for mi = 1:numel(modes)
        mode = modes{mi};
        t_rand = timeit(@() bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
            'block_size', opt.block_size, 'threshold_mode', mode, ...
            'seed', opt.seed + ci), 3);
        [~, ~, ~, st] = bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
            'block_size', opt.block_size, 'threshold_mode', mode, 'seed', opt.seed + ci);
        cond_ratio = st.frob_inv / frobinv_det;     % >1 means worse conditioned than deterministic
        speedup = t_det / t_rand;
        fprintf('%-12s %-20s %11.3f %11.3f %8.2f %8.2f %9.1f\n', ...
            sprintf('%dx%d', k, n), mode, t_rand*1e3, t_det*1e3, speedup, cond_ratio, st.total_tested/k);
        rows(end+1, :) = {opt.family, k, n, mode, t_rand, t_det, speedup, ...
            st.frob_inv, frobinv_det, st.osinsky_bound, st.total_tested/k}; %#ok<AGROW>
    end

    % R12 desired: running_mean selection plus the final Q'-apply to the
    % leftover columns (the O(n k^2) work the deterministic baseline already
    % includes in R). Shows how much speedup survives when R12 is required.
    t_r12 = timeit(@() bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
        'block_size', opt.block_size, 'threshold_mode', 'running_mean', ...
        'return_r12', true, 'seed', opt.seed + ci), 5);
    [~, ~, ~, st12] = bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
        'block_size', opt.block_size, 'threshold_mode', 'running_mean', ...
        'return_r12', true, 'seed', opt.seed + ci);
    fprintf('%-12s %-20s %11.3f %11.3f %8.2f %8.2f %9.1f\n', ...
        sprintf('%dx%d', k, n), 'running_mean+R12', t_r12*1e3, t_det*1e3, ...
        t_det/t_r12, st12.frob_inv/frobinv_det, st12.total_tested/k);
    rows(end+1, :) = {opt.family, k, n, 'running_mean+R12', t_r12, t_det, t_det/t_r12, ...
        st12.frob_inv, frobinv_det, st12.osinsky_bound, st12.total_tested/k}; %#ok<AGROW>
end

results = cell2table(rows, 'VariableNames', {'family', 'k', 'n', 'mode', ...
    't_rand_s', 't_det_s', 'speedup', 'frobinv_rand', 'frobinv_det', ...
    'osinsky', 'tested_per_k'});
csv = fullfile(opt.outdir, 'rand_timings.csv');
writetable(results, csv);
fprintf('\nWrote %s\n', csv);
end
