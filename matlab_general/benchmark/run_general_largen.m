function run_general_largen(varargin)
%RUN_GENERAL_LARGEN Large-n head-to-head: built-in pivoted QR vs general randBSQR.
%
%   Two methods only, one family, one m, over a wide n sweep -- probing
%   whether the wrapper (unpivoted BLAS-3 qr(A','econ') + randomized
%   selection on Qt') becomes PRACTICALLY faster than built-in column-pivoted
%   qr(A,'econ','vector') (LAPACK dgeqp3, whose per-step pivot/norm-update
%   dependency is BLAS-2-bound) as n grows at fixed m:
%     qr_builtin        - [Q,R,p] = qr(A,'econ','vector')
%     general_bsqr_rand - [Q,R,p] = bsqr_general(A,'selector','bsqr_rand',...)
%   Both timed calls materialize Q, R, and the index-vector permutation.
%   The wrapper's phase split (t_qr = qr of A', t_select = selector) is
%   recorded from one instrumented call per point.
%
%   Writes results/exp_general_largen.csv (checkpointed after every n);
%   plot with plot_general_largen.
%
% Options: 'm' (default 256), 'ns' (default [4000 8000 16000 32000 64000
%   128000]; A at the top size is ~260 MB), 'family' (default 'graded_cols',
%   the neutral representative -- see general_test_matrix), 'trials' (seeds,
%   default 5), 'outdir'.

ip = inputParser;
addParameter(ip, 'm', 256);
addParameter(ip, 'ns', [4000, 8000, 16000, 32000, 64000, 128000]);
addParameter(ip, 'family', 'graded_cols');
addParameter(ip, 'trials', 5);
addParameter(ip, 'outdir', '');
parse(ip, varargin{:});
opt = ip.Results;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab'));
addpath(fullfile(repo_root, 'matlab', 'mex'));
addpath(fullfile(repo_root, 'matlab_rand'));
addpath(fullfile(repo_root, 'matlab_rand', 'mex'));
addpath(fullfile(repo_root, 'matlab_rand', 'benchmark'));
addpath(fullfile(repo_root, 'matlab_general'));
addpath(fullfile(repo_root, 'matlab_general', 'benchmark'));
if isempty(opt.outdir)
    opt.outdir = fullfile(repo_root, 'matlab_general', 'benchmark', 'results');
end
if ~isfolder(opt.outdir); mkdir(opt.outdir); end

m = opt.m;
fprintf('large-n general comparison: family=%s, m=%d, trials=%d, BLAS threads=%d\n', ...
    opt.family, m, opt.trials, maxNumCompThreads);
rows = {};
for n = opt.ns
    bound = 1 / sqrt(m * (n - m + 1));
    for s = 1:opt.trials
        seed = 9000 + s;
        [A, info] = general_test_matrix(opt.family, m, n, seed);
        smin_A = info.sigma_min;

        % --- built-in column-pivoted QR on A ---
        t = timeit(@() qr(A, 'econ', 'vector'), 3);
        [~, ~, p] = qr(A, 'econ', 'vector');
        [smin_sel, frobinv] = quality(A, p(1:m));
        rows(end+1, :) = {opt.family, m, n, m, s, 'qr_builtin', t, NaN, NaN, ...
            smin_sel, smin_A, smin_sel / smin_A, bound, frobinv}; %#ok<AGROW>

        % --- the wrapper with the randomized selector ---
        rargs = {'selector', 'bsqr_rand', 'seed', seed, 'pivot_format', 'vector'};
        t = timeit(@() bsqr_general(A, rargs{:}), 3);
        [~, ~, p, st] = bsqr_general(A, rargs{:});
        [smin_sel, frobinv] = quality(A, p(1:m));
        rows(end+1, :) = {opt.family, m, n, m, s, 'general_bsqr_rand', t, ...
            st.t_qr_At, st.t_select, smin_sel, smin_A, smin_sel / smin_A, ...
            bound, frobinv}; %#ok<AGROW>
    end
    write_csv(rows, opt);   % checkpoint
    fprintf('  n=%-7d done\n', n);
end

csv = write_csv(rows, opt);
fprintf('Wrote %s (%d rows)\n', csv, size(rows, 1));
end

% ===========================================================================
function csv = write_csv(rows, opt)
T = cell2table(rows, 'VariableNames', {'family', 'm', 'n', 'k', 'seed', ...
    'method', 'time_s', 't_qr', 't_select', 'sigma_min_sel', 'sigma_min_A', ...
    'ratio', 'bound', 'frobinv'});
csv = fullfile(opt.outdir, 'exp_general_largen.csv');
writetable(T, csv);
end

function [sm, fr] = quality(A, sel)
% Same SVD-based helper as run_general_comparison (no inv()).
sel = sel(:).';
if numel(unique(sel)) < numel(sel)
    sm = 0; fr = Inf;
    return;
end
sv = svd(A(:, sel));
sm = min(sv);
fr = norm(1 ./ sv);
end
