function run_rand_experiments(varargin)
%RUN_RAND_EXPERIMENTS Generate the data behind the publication plots.
%
%   Writes three tidy CSVs under matlab_rand/benchmark/results/, tagged by k
%   (e.g. exp_scaling_k64.csv, exp_scaling_k128.csv, ...):
%     exp_scaling_k<K>.csv   - time / conditioning / samples vs n, per family,
%                              deterministic vs randomized (both modes + R12).
%     exp_blocksize_k<K>.csv - effect of block_size (k-relative sweep) on time /
%                              samples / cond, batched + norm-weighted.
%     exp_sampling_k<K>.csv  - uniform vs norm-weighted sampling across families.
%
%   All randomized timing uses the MEX kernel called directly with
%   check_finite=false (see run_rand_benchmarks for the fairness rationale).
%   Each configuration is repeated over several seeds so the plotter can draw
%   confidence bands. Plot with plot_rand_experiments('k', K).
%
% Options: 'k' (rows to select, default 64), 'trials' (seeds, default 20 --
%   the publication count; timeit already stabilizes the timing, so 20 suffices),
%   'which' (cellstr subset of {'scaling','blocksize','sampling'}, default all),
%   'outdir'.

ip = inputParser;
addParameter(ip, 'k', 64);
addParameter(ip, 'trials', 20);
addParameter(ip, 'which', {'scaling', 'blocksize', 'sampling'});
addParameter(ip, 'outdir', '');
parse(ip, varargin{:});
opt = ip.Results;
if ischar(opt.which) || isstring(opt.which); opt.which = cellstr(opt.which); end
opt.tag = sprintf('_k%d', opt.k);

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
assert(exist('bsqr_rand_mex', 'file') == 3 && exist('bsqr_mex', 'file') == 3, ...
    'Both MEX kernels must be built first.');

fprintf('BLAS threads = %d, trials = %d\n', maxNumCompThreads, opt.trials);

if any(strcmp(opt.which, 'scaling'));   exp_scaling(opt);   end
if any(strcmp(opt.which, 'blocksize')); exp_blocksize(opt); end
if any(strcmp(opt.which, 'sampling'));  exp_sampling(opt);  end
end

% ===========================================================================
function exp_scaling(opt)
families = {'gaussian', 'spiked_leverage', 'needle'};
k = opt.k;
db = k;   % batched default block size (= k; see bsqr_rand_parse_options)
ns = [1000, 2000, 4000, 8000, 16000, 32000, 64000];
rows = {};
for fi = 1:numel(families)
    fam = families{fi};
    for n = ns
        for s = 1:opt.trials
            M = rand_test_matrix(fam, k, n, 1000 * fi + s);
            d = measure_det(M, k);
            rows(end+1, :) = {fam, k, n, s, 'det', 'na', 'na', 0, d.time, d.frobinv, d.osinsky, NaN, d.sigma_min}; %#ok<AGROW>
            b = measure_builtin(M, k);
            rows(end+1, :) = {fam, k, n, s, 'builtin', 'na', 'na', 0, b.time, b.frobinv, b.osinsky, NaN, b.sigma_min}; %#ok<AGROW>
            % Randomized BSQR. Norm-weighted (column-norm) sampling is the standard
            % choice for cross-method comparisons -- robust across leverage profiles
            % and matching the sampling used by the rejection_rpqr comparison.
            % running_mean is the default threshold; worstcase_allowance is a
            % documented option kept out of the plots for a cleaner narrative.
            r = measure_rand(M, k, 'running_mean', 'normweighted', db, 1000 * fi + s);
            rows(end+1, :) = {fam, k, n, s, 'rand', 'running_mean', 'normweighted', db, ...
                r.time, r.frobinv, r.osinsky, r.tested_per_k, r.sigma_min}; %#ok<AGROW>
            % R12 desired: same selection plus the final Q'-apply to the leftover
            % columns. Quantifies the advantage lost when R12 cannot be skipped.
            r = measure_rand(M, k, 'running_mean', 'normweighted', db, 1000 * fi + s, true);
            rows(end+1, :) = {fam, k, n, s, 'rand_r12', 'running_mean', 'normweighted', db, ...
                r.time, r.frobinv, r.osinsky, r.tested_per_k, r.sigma_min}; %#ok<AGROW>
        end
        fprintf('  scaling %-16s n=%-6d done\n', fam, n);
    end
end
write_rows(rows, fullfile(opt.outdir, ['exp_scaling' opt.tag '.csv']));
end

% ===========================================================================
function exp_blocksize(opt)
families = {'gaussian', 'spiked_leverage', 'needle'};
k = opt.k; n = 8000;
% k-relative sweep that brackets the batched default block = k (with finer steps
% at 3/4 k and 3/2 k around the knee, block=1 as the no-batching reference, and up
% to 4k for the conditioning-recovery tail). All stay < n.
blocks = unique([1, round(k .* [1/8 1/4 1/2 3/4 1 1.5 2 4])]);
samplings = {'normweighted'};   % batched default; uniform-vs-norm lives in exp_sampling
rows = {};
for fi = 1:numel(families)
    fam = families{fi};
    for s = 1:opt.trials
        M = rand_test_matrix(fam, k, n, 2000 * fi + s);
        b = measure_builtin(M, k);   % vendor-qr reference (constant across block size)
        rows(end+1, :) = {fam, k, n, s, 'builtin', 'na', 'na', 0, ...
            b.time, b.frobinv, b.osinsky, NaN, b.sigma_min}; %#ok<AGROW>
        for bi = 1:numel(blocks)
            for si = 1:numel(samplings)
                r = measure_rand(M, k, 'running_mean', samplings{si}, blocks(bi), 2000 * fi + s);
                rows(end+1, :) = {fam, k, n, s, 'rand', 'running_mean', samplings{si}, blocks(bi), ...
                    r.time, r.frobinv, r.osinsky, r.tested_per_k, r.sigma_min}; %#ok<AGROW>
            end
        end
    end
    fprintf('  blocksize %-16s done\n', fam);
end
write_rows(rows, fullfile(opt.outdir, ['exp_blocksize' opt.tag '.csv']));
end

% ===========================================================================
function exp_sampling(opt)
% chebyshev is the applied family (polynomial least-squares / optimal design,
% edge-concentrated leverage) -- included here to motivate the synthetic
% leverage profiles; kept out of the scaling/rpqr suites where it just blends.
families = {'gaussian', 'graded_leverage', 'spiked_leverage', 'coherent', 'needle', 'chebyshev'};
k = opt.k; n = 32000; block = k;   % batched default block size
samplings = {'uniform', 'normweighted'};
rows = {};
for fi = 1:numel(families)
    fam = families{fi};
    for s = 1:opt.trials
        M = rand_test_matrix(fam, k, n, 3000 * fi + s);
        b = measure_builtin(M, k);   % vendor-qr reference (constant across sampling)
        rows(end+1, :) = {fam, k, n, s, 'builtin', 'na', 'na', 0, ...
            b.time, b.frobinv, b.osinsky, NaN, b.sigma_min}; %#ok<AGROW>
        for si = 1:numel(samplings)
            r = measure_rand(M, k, 'running_mean', samplings{si}, block, 3000 * fi + s);
            rows(end+1, :) = {fam, k, n, s, 'rand', 'running_mean', samplings{si}, block, ...
                r.time, r.frobinv, r.osinsky, r.tested_per_k, r.sigma_min}; %#ok<AGROW>
        end
    end
    fprintf('  sampling %-16s done\n', fam);
end
write_rows(rows, fullfile(opt.outdir, ['exp_sampling' opt.tag '.csv']));
end

% ===========================================================================
function out = measure_rand(M, k, mode, sampling, block, seed, return_r12)
if nargin < 7; return_r12 = false; end
if return_r12
    % R12 desired: the kernel must apply the accumulated Q' to the leftover
    % columns (one O(n k^2) pass). Time the full 5-output call.
    f = @() bsqr_rand_mex(M, 'k', k, 'check_finite', false, 'block_size', block, ...
        'threshold_mode', mode, 'sampling', sampling, 'return_r12', true, 'seed', seed);
    out.time = timeit(f, 5);
    [~, ~, R11, st] = bsqr_rand_mex(M, 'k', k, 'check_finite', false, 'block_size', block, ...
        'threshold_mode', mode, 'sampling', sampling, 'return_r12', true, 'seed', seed);
else
    f = @() bsqr_rand_mex(M, 'k', k, 'check_finite', false, 'block_size', block, ...
        'threshold_mode', mode, 'sampling', sampling, 'seed', seed);
    out.time = timeit(f, 3);
    [~, ~, R11, st] = f();
end
[out.frobinv, out.sigma_min] = svd_quality(triu(R11));
out.osinsky = st.osinsky_bound;
out.tested_per_k = st.total_tested / k;
end

function out = measure_det(M, k)
out.time = timeit(@() bsqr_mex(M, 'k', k, 'check_finite', false), 1);
R = bsqr_mex(M, 'k', k, 'check_finite', false);
[out.frobinv, out.sigma_min] = svd_quality(triu(R(1:k, 1:k)));
out.osinsky = sqrt(k * (size(M, 2) - k + 1));
end

function out = measure_builtin(M, k)
% Built-in column-pivoted QR (LAPACK dgeqp3 via MATLAB qr) -- a different,
% vendor-tuned classical algorithm (Businger-Golub max-norm pivoting). It is the
% natural "can the randomized method beat the library?" baseline. Requesting the
% pivot vector forces the full [Q,R,e] factorization (MATLAB has no pivoted-QR
% interface that skips Q), so this is the standard library call a user would make
% to get a column selection.
out.time = timeit(@() qr(M, 'econ', 'vector'), 3);
[~, R, ~] = qr(M, 'econ', 'vector');
[out.frobinv, out.sigma_min] = svd_quality(triu(R(1:k, 1:k)));
out.osinsky = sqrt(k * (size(M, 2) - k + 1));
end

function [fr, sm] = svd_quality(R11)
% Both selection-quality metrics from a single SVD of R11 (no inv):
%   fr = ||R11^{-1}||_F = sqrt(sum 1/sigma_i^2),  sm = sigma_min(R11).
s = svd(R11);
fr = norm(1 ./ s);
sm = min(s);
end

function write_rows(rows, path)
T = cell2table(rows, 'VariableNames', {'family', 'k', 'n', 'seed', 'method', ...
    'mode', 'sampling', 'block_size', 'time_s', 'frobinv', 'osinsky', 'tested_per_k', 'sigma_min'});
writetable(T, path);
fprintf('Wrote %s (%d rows)\n', path, height(T));
end
