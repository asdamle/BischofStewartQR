function run_rpqr_comparison(varargin)
%RUN_RPQR_COMPARISON Randomized BSQR vs rejection_rpqr (Adaptive Randomized Pivoting).
%
%   A standalone comparison (separate CSVs/plots from the main experiment suite)
%   between bsqr_rand and the `rejection_rpqr` selector from
%       https://github.com/eepperly/Adaptive-Randomized-Pivoting
%   which must be downloaded into ext_comparisons/ (that folder is git-ignored;
%   only the comparison harness and the resulting plots are tracked).
%
%   Same three metrics as the main suite: runtime, speedup, and selection quality
%   ||R11^{-1}||_F. Writes results/exp_rpqr_k<K>.csv; plot with
%   plot_rpqr_comparison. Caveats worth stating alongside the plots:
%     * rejection_rpqr is the published .m implementation (with a small
%       rejection_helper MEX); bsqr_rand is our MEX -- so timing reflects both the
%       algorithm and the implementation language.
%     * The two methods optimize different objectives: BSQR directly minimizes the
%       per-step growth of ||R11^{-1}||_F, while ARP/rejection_rpqr targets a
%       volume/DPP criterion. The ||R11^{-1}||_F metric therefore favors BSQR by
%       construction -- it measures BSQR's objective.
%
% Options: 'k' (default 64), 'trials' (seeds, default 5), 'ns' (vector),
%   'families' (cellstr; see rand_test_matrix), 'outdir'.

ip = inputParser;
addParameter(ip, 'k', 64);
addParameter(ip, 'trials', 5);
addParameter(ip, 'ns', [1000, 2000, 4000, 8000, 16000, 32000]);
addParameter(ip, 'families', {'gaussian', 'spiked_leverage', 'needle'});
addParameter(ip, 'outdir', '');
parse(ip, varargin{:});
opt = ip.Results;
opt.tag = sprintf('_k%d', opt.k);

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab_rand'));
addpath(fullfile(repo_root, 'matlab_rand', 'mex'));
addpath(fullfile(repo_root, 'matlab_rand', 'benchmark'));
arp = fullfile(repo_root, 'ext_comparisons', 'Adaptive-Randomized-Pivoting-main');
assert(isfolder(arp), ['rejection_rpqr not found. Download Adaptive-Randomized-Pivoting ', ...
    'into ext_comparisons/ (see matlab_rand/README.md).']);
addpath(fullfile(arp, 'code'));
addpath(fullfile(arp, 'utils'));
assert(exist('rejection_rpqr', 'file') > 0, 'rejection_rpqr is not on the path.');
assert(exist('bsqr_rand_mex', 'file') == 3, 'bsqr_rand_mex is not built.');
if isempty(opt.outdir)
    opt.outdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'results');
end
if ~isfolder(opt.outdir); mkdir(opt.outdir); end

k = opt.k;
fprintf('rejection_rpqr comparison: k=%d, trials=%d, BLAS threads=%d\n', k, opt.trials, maxNumCompThreads);
rows = {};
for fi = 1:numel(opt.families)
    fam = opt.families{fi};
    for n = opt.ns
        for s = 1:opt.trials
            seed = 7000 * fi + s;
            M = rand_test_matrix(fam, k, n, seed);
            os = sqrt(k * (n - k + 1));

            % --- randomized BSQR (MEX, [p,reflectors,R11], default block) ---
            % Norm-weighted sampling, matching rejection_rpqr's squared-column-norm
            % sampling -- the apples-to-apples choice (and the documented one for
            % concentrated leverage).
            tb = timeit(@() bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
                'sampling', 'normweighted', 'seed', seed), 3);
            pb = bsqr_rand_mex(M, 'k', k, 'check_finite', false, ...
                'sampling', 'normweighted', 'seed', seed);
            [fr_b, sm_b] = quality(M, pb(1:k));
            rows(end+1, :) = {fam, k, n, s, 'bsqr', tb, fr_b, sm_b, os}; %#ok<AGROW>

            % --- rejection_rpqr (default proposal size l = k) ---
            tr = timeit(@() rejection_rpqr(M, k), 3);
            rng(seed);
            idx = rejection_rpqr(M, k);
            [fr_r, sm_r] = quality(M, idx(1:k));
            rows(end+1, :) = {fam, k, n, s, 'rejection_rpqr', tr, fr_r, sm_r, os}; %#ok<AGROW>
        end
        fprintf('  %-16s n=%-6d done\n', fam, n);
    end
end

T = cell2table(rows, 'VariableNames', {'family', 'k', 'n', 'seed', 'method', ...
    'time_s', 'frobinv', 'sigma_min', 'osinsky'});
csv = fullfile(opt.outdir, ['exp_rpqr' opt.tag '.csv']);
writetable(T, csv);
fprintf('Wrote %s (%d rows)\n', csv, height(T));
end

function [fr, sm] = quality(M, sel)
% Selection-quality metrics of the chosen k columns from one SVD (no inv):
%   fr = ||R11^{-1}||_F,  sm = sigma_min(R11). For orthonormal-row input (m=k) the
% selected block is square and its singular values equal those of R11 (Frobenius
% and singular values are invariant under the orthogonal Q), so this is
% backend-agnostic and identical for both methods. fr = Inf / sm = 0 if the
% selection is rank-deficient (reported honestly).
sel = sel(:).';
if numel(unique(sel)) < numel(sel)
    fr = Inf; sm = 0;
    return;
end
s = svd(M(:, sel));
fr = norm(1 ./ s);
sm = min(s);
end
