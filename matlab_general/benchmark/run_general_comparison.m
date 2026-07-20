function run_general_comparison(varargin)
%RUN_GENERAL_COMPARISON bsqr_general vs direct selection on general A.
%
%   For each family/size/seed, times and scores five methods that each select
%   k = m columns of a general short-wide A. Every timed call materializes
%   the factors plus the permutation as an index VECTOR (the repo's baseline
%   convention -- never the matrix pivot format):
%     qr_builtin        - [Q,R,p] = qr(A,'econ','vector')  (column-pivoted)
%     bsqr_direct       - bsqr(A,...) on A itself (outside its theory)
%     bsqr_rand_direct  - bsqr_rand(A,...) on A itself (outside its theory)
%     general_bsqr      - bsqr_general(A): qr(A','econ') then bsqr on Qt'
%     general_bsqr_rand - bsqr_general(A,'selector','bsqr_rand')
%
%   Quality of the selected k columns S: sigma_min(A(:,S)), its ratio to
%   sigma_min(A) (the general_* methods guarantee ratio >= 1/sqrt(m(n-m+1))
%   at k = m), and frobinv = ||A(:,S)^{-1}||_F (= norm(1./svd), no inv()).
%   For the wrapper methods one instrumented call also records the phase
%   split t_qr (qr of A') / t_select (selector on Qt').
%
%   Sizes: an n-sweep at fixed m plus an m-sweep at fixed n. Writes
%   results/exp_general.csv; plot with plot_general_comparison.
%
% Options: 'm_fixed' (default 64), 'ns' (default [500 1000 2000 4000 8000]),
%   'n_fixed' (default 4000), 'ms' (default [32 128 256]; the m-sweep runs at
%   n_fixed, and m_fixed x n_fixed already comes from the n-sweep),
%   'families' (default all five in general_test_matrix), 'trials' (seeds,
%   default 10), 'k' (default [] -> m; an explicit k tags the CSV and voids
%   the bound column), 'outdir', 'exclude' (cell of {method, family} pairs
%   recorded as NaN rows instead of being run; default skips
%   bsqr_rand_direct on spectrum_needle -- bsqr_rand's MEX has a KNOWN
%   OPEN BUG on that family's concentrated-norm general input: it can return
%   pivot index n+1 and corrupt the heap, eventually crashing MATLAB
%   (observed 2026-07-20; pass 'exclude',{} to run it anyway).
%
%   The CSV is checkpointed after every family/size block, so a crash loses
%   at most one block.

ip = inputParser;
addParameter(ip, 'm_fixed', 64);
addParameter(ip, 'ns', [500, 1000, 2000, 4000, 8000]);
addParameter(ip, 'n_fixed', 4000);
addParameter(ip, 'ms', [32, 128, 256]);
addParameter(ip, 'families', {'graded_cols', 'loud_collinear', ...
    'spectrum_needle', 'spectrum_coherent', 'kahan_decoy'});
addParameter(ip, 'trials', 10);
addParameter(ip, 'k', []);
addParameter(ip, 'outdir', '');
addParameter(ip, 'exclude', {{'bsqr_rand_direct', 'spectrum_needle'}});
parse(ip, varargin{:});
opt = ip.Results;
if isempty(opt.k); opt.tag = ''; else; opt.tag = sprintf('_k%d', opt.k); end

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab'));
addpath(fullfile(repo_root, 'matlab', 'mex'));
addpath(fullfile(repo_root, 'matlab_rand'));
addpath(fullfile(repo_root, 'matlab_rand', 'mex'));
addpath(fullfile(repo_root, 'matlab_rand', 'benchmark'));   % rand_test_matrix
addpath(fullfile(repo_root, 'matlab_general'));
addpath(fullfile(repo_root, 'matlab_general', 'benchmark'));
if isempty(opt.outdir)
    opt.outdir = fullfile(repo_root, 'matlab_general', 'benchmark', 'results');
end
if ~isfolder(opt.outdir); mkdir(opt.outdir); end

% Size pairs: n-sweep at m_fixed plus m-sweep at n_fixed (deduplicated,
% short-wide only).
pairs = [repmat(opt.m_fixed, numel(opt.ns), 1), opt.ns(:)];
pairs = [pairs; opt.ms(:), repmat(opt.n_fixed, numel(opt.ms), 1)];
pairs = unique(pairs, 'rows', 'stable');
pairs = pairs(pairs(:, 1) <= pairs(:, 2), :);

fprintf('general-A comparison: %d size pairs, trials=%d, BLAS threads=%d\n', ...
    size(pairs, 1), opt.trials, maxNumCompThreads);
for i = 1:numel(opt.exclude)
    fprintf('  excluding %s on %s (NaN rows; see header note)\n', ...
        opt.exclude{i}{1}, opt.exclude{i}{2});
end
rows = {};
for fi = 1:numel(opt.families)
    fam = opt.families{fi};
    for pi_ = 1:size(pairs, 1)
        m = pairs(pi_, 1); n = pairs(pi_, 2);
        if isempty(opt.k); k = m; else; k = opt.k; end
        for s = 1:opt.trials
            seed = 7000 * fi + s;
            [A, info] = general_test_matrix(fam, m, n, seed);
            smin_A = info.sigma_min;
            if k == m; bound = 1 / sqrt(m * (n - m + 1)); else; bound = NaN; end

            % --- built-in column-pivoted QR on A ---
            rows = run_method(rows, fam, m, n, k, s, 'qr_builtin', ...
                @() qr(A, 'econ', 'vector'), 3, ...
                @() third_output(@() qr(A, 'econ', 'vector')), ...
                [], A, smin_A, bound);

            % --- bsqr directly on A ---
            rows = run_method(rows, fam, m, n, k, s, 'bsqr_direct', ...
                @() bsqr(A, 'k', k, 'pivot_format', 'vector'), 3, ...
                @() third_output(@() bsqr(A, 'k', k, 'pivot_format', 'vector')), ...
                [], A, smin_A, bound);

            % --- bsqr_rand directly on A (outside its theory) ---
            if is_excluded(opt.exclude, 'bsqr_rand_direct', fam)
                rows(end+1, :) = {fam, m, n, k, s, 'bsqr_rand_direct', ...
                    NaN, NaN, NaN, NaN, smin_A, NaN, bound, NaN}; %#ok<AGROW>
            else
                rows = run_method(rows, fam, m, n, k, s, 'bsqr_rand_direct', ...
                    @() bsqr_rand(A, 'k', k, 'seed', seed), 3, ...
                    @() first_output(@() bsqr_rand(A, 'k', k, 'seed', seed)), ...
                    [], A, smin_A, bound);
            end

            % --- the wrapper, both selectors ---
            wargs = {'k', k, 'pivot_format', 'vector'};
            rows = run_method(rows, fam, m, n, k, s, 'general_bsqr', ...
                @() bsqr_general(A, wargs{:}), 3, [], wargs, A, smin_A, bound);
            rargs = {'k', k, 'selector', 'bsqr_rand', 'seed', seed, ...
                'pivot_format', 'vector'};
            rows = run_method(rows, fam, m, n, k, s, 'general_bsqr_rand', ...
                @() bsqr_general(A, rargs{:}), 3, [], rargs, A, smin_A, bound);
        end
        write_csv(rows, opt);   % checkpoint: a crash loses at most one block
        fprintf('  %-18s m=%-4d n=%-6d done\n', fam, m, n);
    end
end

csv = write_csv(rows, opt);
fprintf('Wrote %s (%d rows)\n', csv, size(rows, 1));
end

% ===========================================================================
function csv = write_csv(rows, opt)
T = cell2table(rows, 'VariableNames', {'family', 'm', 'n', 'k', 'seed', ...
    'method', 'time_s', 't_qr', 't_select', 'sigma_min_sel', 'sigma_min_A', ...
    'ratio', 'bound', 'frobinv'});
csv = fullfile(opt.outdir, ['exp_general' opt.tag '.csv']);
writetable(T, csv);
end

function tf = is_excluded(exclude, name, fam)
tf = false;
for i = 1:numel(exclude)
    if strcmp(exclude{i}{1}, name) && strcmp(exclude{i}{2}, fam)
        tf = true;
        return;
    end
end
end

% ===========================================================================
function rows = run_method(rows, fam, m, n, k, s, name, thunk, nout, pivots, ...
    wrapper_args, A, smin_A, bound)
% Time the materializing call, then one untimed call for the pivots (and, for
% the wrapper methods, the phase split from STATS). A failing method is
% recorded as a NaN row rather than aborting the sweep.
try
    t = timeit(thunk, nout);
    t_qr = NaN; t_sel = NaN;
    if isempty(pivots)   % wrapper method: pivots and phase stats in one call
        [~, ~, p, st] = bsqr_general(A, wrapper_args{:});
        t_qr = st.t_qr_At; t_sel = st.t_select;
    else
        p = pivots();
    end
    if any(p(1:k) < 1) || any(p(1:k) > size(A, 2))
        % Observed from bsqr_rand's MEX on concentrated-norm general A after a
        % long call history (out-of-range pivot indices) -- record the failure
        % honestly instead of crashing the sweep.
        error('run_general_comparison:InvalidPivots', ...
            'selected pivot indices out of range (max %d, n = %d)', ...
            max(p(1:k)), size(A, 2));
    end
    [smin_sel, frobinv] = quality(A, p(1:k));
catch ME
    warning('run_general_comparison:MethodFailed', ...
        '%s failed on %s m=%d n=%d seed=%d: %s', name, fam, m, n, s, ME.message);
    [t, t_qr, t_sel, smin_sel, frobinv] = deal(NaN);
end
rows(end+1, :) = {fam, m, n, k, s, name, t, t_qr, t_sel, ...
    smin_sel, smin_A, smin_sel / smin_A, bound, frobinv};
end

function [sm, fr] = quality(A, sel)
% Quality of the selected columns from one SVD (no inv()):
%   sm = sigma_min(A(:,sel)),  fr = ||R11^{-1}||_F = norm(1./svd(A(:,sel))).
% Inf/0 if the selection is rank-deficient or duplicated (reported honestly).
sel = sel(:).';
if numel(unique(sel)) < numel(sel)
    sm = 0; fr = Inf;
    return;
end
sv = svd(A(:, sel));
sm = min(sv);
fr = norm(1 ./ sv);
end

function p = third_output(fn)
[~, ~, p] = fn();
end

function p = first_output(fn)
p = fn();
end
