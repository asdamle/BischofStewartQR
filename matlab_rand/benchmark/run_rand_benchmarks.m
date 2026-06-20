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
%       -- the "R12 not needed" product. No Q, no R12. Norm-weighted (column-norm)
%       sampling, the standard choice for cross-method comparisons.
%     * Built-in baseline: qr(M,'econ','vector') -- LAPACK dgeqp3, a different,
%       vendor-tuned classical algorithm. spd_qr = t_qr / t_rand shows whether
%       the randomized method beats the library call.
%
% Name-value options: 'sizes' (cell of [k n]), 'seed', 'block_size',
%   'family' (see rand_test_matrix), 'outdir'.

ip = inputParser;
addParameter(ip, 'sizes', {[32, 2000], [64, 4000], [64, 8000], [128, 8000], [128, 16000]});
addParameter(ip, 'seed', 1);
addParameter(ip, 'block_size', []);   % [] = auto (rand_default_block(k) per size)
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

rows = {};
if isempty(opt.block_size); bsdesc = 'auto(k)'; else; bsdesc = num2str(opt.block_size); end
fprintf('family=%s  block_size=%s  BLAS threads=%d\n', opt.family, bsdesc, maxNumCompThreads);
fprintf('%-12s %-18s %10s %10s %10s %8s %8s %7s %8s\n', ...
    'size', 'mode', 't_rand(ms)', 't_det(ms)', 't_qr(ms)', 'spd_det', 'spd_qr', 'cond', 'tested/k');

for ci = 1:numel(opt.sizes)
    k = opt.sizes{ci}(1); n = opt.sizes{ci}(2);
    M = rand_test_matrix(opt.family, k, n, opt.seed + ci);
    if isempty(opt.block_size); bs = k; else; bs = opt.block_size; end   % batched default block

    t_det = timeit(@() bsqr_mex(M, 'k', k, 'check_finite', false), 1);
    Rdet = bsqr_mex(M, 'k', k, 'check_finite', false);
    frobinv_det = norm(1 ./ svd(triu(Rdet(1:k, 1:k))));   % ||R11^{-1}||_F via singular values

    % Built-in column-pivoted QR (LAPACK dgeqp3): vendor-tuned classical baseline.
    t_qr = timeit(@() qr(M, 'econ', 'vector'), 3);
    [~, Rqr, ~] = qr(M, 'econ', 'vector');
    frobinv_qr = norm(1 ./ svd(triu(Rqr(1:k, 1:k))));   % ||R11^{-1}||_F via singular values

    osinsky = sqrt(k * (n - k + 1));
    % running_mean is the default; worstcase_allowance is a documented option (see
    % bsqr_rand help) but kept out of the headline comparison for clarity.
    t_rand = timeit(@() bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
        'block_size', bs, 'sampling', 'normweighted', 'threshold_mode', 'running_mean', 'seed', opt.seed + ci), 3);
    [~, ~, ~, st] = bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
        'block_size', bs, 'sampling', 'normweighted', 'threshold_mode', 'running_mean', 'seed', opt.seed + ci);
    rows(end+1, :) = report_row(opt.family, k, n, 'running_mean', t_rand, t_det, t_qr, ...
        st.frob_inv, frobinv_det, frobinv_qr, osinsky, st.total_tested/k); %#ok<AGROW>

    % R12 desired: running_mean selection plus the final Q'-apply to the leftover
    % columns -- shows how much speedup survives when R12 is required.
    t_r12 = timeit(@() bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
        'block_size', bs, 'sampling', 'normweighted', 'threshold_mode', 'running_mean', ...
        'return_r12', true, 'seed', opt.seed + ci), 5);
    [~, ~, ~, st12] = bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
        'block_size', bs, 'sampling', 'normweighted', 'threshold_mode', 'running_mean', ...
        'return_r12', true, 'seed', opt.seed + ci);
    rows(end+1, :) = report_row(opt.family, k, n, 'running_mean+R12', t_r12, t_det, t_qr, ...
        st12.frob_inv, frobinv_det, frobinv_qr, osinsky, st12.total_tested/k); %#ok<AGROW>
end

results = cell2table(rows, 'VariableNames', {'family', 'k', 'n', 'mode', ...
    't_rand_s', 't_det_s', 't_builtin_s', 'speedup_det', 'speedup_builtin', ...
    'frobinv_rand', 'frobinv_det', 'frobinv_builtin', 'osinsky', 'tested_per_k'});
csv = fullfile(opt.outdir, 'rand_timings.csv');
writetable(results, csv);
fprintf('\nWrote %s\n', csv);
end

function row = report_row(family, k, n, mode, t_rand, t_det, t_qr, ...
                         fr_rand, fr_det, fr_qr, osinsky, tested)
% Print one table line and return its CSV row. spd_det / spd_qr are the
% randomized speedups over deterministic BSQR and built-in pivoted QR; cond is
% ||R11^{-1}||_F(rand) / ||R11^{-1}||_F(det).
spd_det = t_det / t_rand;
spd_qr = t_qr / t_rand;
cond = fr_rand / fr_det;
fprintf('%-12s %-18s %10.3f %10.3f %10.3f %8.2f %8.2f %7.2f %8.1f\n', ...
    sprintf('%dx%d', k, n), mode, t_rand*1e3, t_det*1e3, t_qr*1e3, spd_det, spd_qr, cond, tested);
row = {family, k, n, mode, t_rand, t_det, t_qr, spd_det, spd_qr, ...
       fr_rand, fr_det, fr_qr, osinsky, tested};
end
