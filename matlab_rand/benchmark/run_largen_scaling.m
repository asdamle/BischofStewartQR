function run_largen_scaling(varargin)
%RUN_LARGEN_SCALING Fixed-k, large-n runtime scaling of three column samplers.
%
%   Compares, at fixed k and growing n, the wall-clock time of:
%     bsqr_nw        - randomized BSQR, norm-weighted sampling (the default)
%     bsqr_unif      - randomized BSQR, uniform sampling
%     rejection_rpqr - Adaptive Randomized Pivoting (squared-column-norm proposal)
%   on two leverage regimes chosen to separate the samplers:
%     gaussian - near-uniform leverage; uniform sampling is fine.
%     needle   - ~k high-leverage columns hidden among many near-null ones;
%                uniform sampling keeps missing them and degenerates to the
%                global-min fallback, while norm-weighted sampling finds them.
%   The message is the time-vs-n scaling: norm-weighted BSQR and rejection_rpqr
%   (both norm-weighted, so both carry the common O(mn) norm work) scale alike,
%   while uniform BSQR is fastest on gaussian but blows up on needle.
%
%   Measurement matches run_rpqr_comparison for apples-to-apples timing: the
%   compiled kernel is called directly (bsqr_rand_mex), 'check_finite' is false
%   (the O(mn) finiteness scan is not part of either algorithm), and timeit
%   stabilizes each measurement. Writes results/exp_largen.csv; plot with
%   plot_largen_scaling.
%
% Options:
%   'k'        (default 64) columns to select.
%   'trials'   (default 20) seeds per (family, n); the plot's band is their range.
%   'ns'       (default 1000*2.^(0:11) = 1000 .. 2048000) the n sweep. The top
%              size allocates a k-by-n matrix (~1 GB at n = 2e6; generation
%              peaks near ~3 GB during orth) -- sized so the sweep ends one
%              doubling PAST the ~0.5 GB memory-bandwidth regime crossing
%              rather than exactly on it, letting the curves re-linearize in
%              the new regime. Do not extend further on 16 GB machines: the
%              next doubling's generation peak (~6 GB) risks swap-polluted
%              timings.
%   'families' (default {'gaussian','needle'}) any rand_test_matrix families.
%   'outdir'   (default matlab_rand/benchmark/results).

ip = inputParser;
addParameter(ip, 'k', 64);
addParameter(ip, 'trials', 20);
addParameter(ip, 'ns', round(1000 * 2 .^ (0:11)));
addParameter(ip, 'families', {'gaussian', 'needle'});
addParameter(ip, 'outdir', '');
parse(ip, varargin{:});
opt = ip.Results;
if ischar(opt.families) || isstring(opt.families); opt.families = cellstr(opt.families); end

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab_rand'));
addpath(fullfile(repo_root, 'matlab_rand', 'mex'));
addpath(fullfile(repo_root, 'matlab_rand', 'benchmark'));
arp = fullfile(repo_root, 'ext_comparisons', 'Adaptive-Randomized-Pivoting-main');
assert(isfolder(arp), ['rejection_rpqr not found. Download Adaptive-Randomized-Pivoting ', ...
    'into ext_comparisons/ (see matlab_rand/README.md).']);
addpath(fullfile(arp, 'code')); addpath(fullfile(arp, 'utils'));
assert(exist('rejection_rpqr', 'file') > 0, 'rejection_rpqr is not on the path.');

if isempty(opt.outdir)
    opt.outdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'results');
end
if ~isfolder(opt.outdir); mkdir(opt.outdir); end

k = opt.k;
fprintf('largen scaling: k=%d, trials=%d, ns=[%s], BLAS threads=%d\n', ...
    k, opt.trials, strjoin(string(opt.ns), ' '), maxNumCompThreads);

rows = cell(0, 6);
for fi = 1:numel(opt.families)
    fam = opt.families{fi};
    for n = opt.ns
        for s = 1:opt.trials
            seed = 9000 * fi + s;
            M = rand_test_matrix(fam, k, n, seed);

            % Direct MEX, check_finite=false, timeit -- matches run_rpqr_comparison.
            tnw = timeit(@() bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
                'sampling', 'normweighted', 'seed', seed), 3);
            tun = timeit(@() bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
                'sampling', 'uniform', 'seed', seed), 3);
            trp = timeit(@() rejection_rpqr(M, k), 3);

            rows(end+1, :) = {fam, k, n, s, 'bsqr_nw',        tnw}; %#ok<AGROW>
            rows(end+1, :) = {fam, k, n, s, 'bsqr_unif',      tun}; %#ok<AGROW>
            rows(end+1, :) = {fam, k, n, s, 'rejection_rpqr', trp}; %#ok<AGROW>
        end
        fprintf('  %-10s n=%-8d done (nw/unif/rpqr last seed: %.4g / %.4g / %.4g s)\n', ...
            fam, n, tnw, tun, trp);
    end
end

T = cell2table(rows, 'VariableNames', {'family', 'k', 'n', 'seed', 'method', 'time_s'});
csv = fullfile(opt.outdir, 'exp_largen.csv');
writetable(T, csv);
fprintf('Wrote %s\n', csv);
end
