function run_approx_comparison(varargin)
%RUN_APPROX_COMPARISON Low-rank approximation quality: randomized BSQR vs
%   rejection_rpqr on application-driven matrices.
%
%   A standalone, application-facing companion to run_rpqr_comparison. Instead of
%   scoring the selection by ||R11^{-1}||_F (BSQR's own objective, which favors
%   it), this scores it by genuine low-rank *approximation* error in both the
%   Frobenius and spectral (2-) norms -- the metric a practitioner cares about.
%
%   The pipeline is the standard CSSP-via-leading-singular-vectors construction,
%   exactly how ARP's own arp.m uses rejection_rpqr:
%     1. Build a full application matrix A (fixed per family; see
%        approx_test_matrix / approx_real_matrix).
%     2. Compute leading right singular vectors with svds (Lanczos):
%        [~,~,V] = svds(A, k); W = V' is k-by-n with orthonormal rows. Both
%        selectors receive the *same* W, so the comparison isolates the
%        column-selection step, not the (shared, accurate) subspace estimate.
%     3. Each method picks k of the n columns of A; the chosen subset is scored by
%        the orthogonal-projection error ||A - P_S A|| (interpolation-free, hence
%        identical scoring for both methods).
%     4. Reference: the optimal rank-k error from the *exact* singular spectrum
%        (one dense svd per family -- matrices are moderate). The Frobenius optimum
%        is the direct tail sum sqrt(sum_{i>k} sigma_i^2) and the spectral optimum
%        is sigma_{k+1}; computing the tail directly (not as ||A||_F^2 - top-k)
%        avoids catastrophic cancellation once the spectrum decays near eps, so the
%        reference is a faithful lower bound even where the error is tiny.
%
%   The matrix is fixed per family (matching the fixed-matrix design of the ARP
%   accuracy study); trials vary only the selectors' RNG, so the seed bands show
%   selection variability on one problem. svds is computed once per (family,k) and
%   reused across trials.
%
%   Alongside the accuracy metrics it also records the interpolation-coefficient
%   magnitude maxT = max|R11^{-1}R12| (in the leading-k frame W=V_k'): the
%   basis-dependent quantity an ID/CUR pays for, paired here with the projection
%   error so one figure shows "accuracy ~identical, coefficients differ".
%
%   Writes results/exp_approx<tag>.csv with one row per (family, k, seed, method):
%     family, m, n, k, seed, method, frob_err, spec_err, opt_frob, opt_spec,
%     normA_fro, normA_2, time_s, maxT
%   Plot with plot_approx_comparison.
%
% Options:
%   'ks'       (default [4 8 16 32 64 128 256]) target ranks to sweep.
%   'rel_floor'(default 1e-11) per-family cap: a k is swept only while the optimal
%              rank-k *relative* Frobenius error stays above this floor. Past it the
%              spectrum has decayed near eps and every method sits on the
%              machine-precision floor with nothing to distinguish. Dropped k's are
%              reported (no silent truncation).
%   'trials'   (default 20) selector-RNG repeats per (family,k).
%   'families' (default {'gmm_kernel','integral_skeleton','snapshots'}) any mix of
%              application families (approx_test_matrix), synthetic-spectrum
%              families {'gaussian','spiked_leverage','needle'} (approx_synth_matrix
%              -- see run_approx_synth_comparison), and real names (loaded from
%              ext_comparisons/data/; absent ones are skipped with a printed note).
%   'sizes'    (struct, default struct()) per-family size overrides forwarded to the
%              builder (e.g. struct('gmm_kernel',struct('N',1500)) or
%              struct('needle',struct('r',400))).
%   'tag'      (default '') suffix on the output file: exp_approx<tag>.csv. Use a
%              tag (e.g. '_synth') to keep companion runs from clobbering each other.
%   'outdir'   (default matlab_rand/benchmark/results).

ip = inputParser;
addParameter(ip, 'ks', [4 8 16 32 64 128 256]);
addParameter(ip, 'rel_floor', 1e-11);
addParameter(ip, 'trials', 20);
addParameter(ip, 'families', {'gmm_kernel', 'integral_skeleton', 'snapshots'});
addParameter(ip, 'sizes', struct());
addParameter(ip, 'tag', '');
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
addpath(fullfile(arp, 'code'));
addpath(fullfile(arp, 'utils'));
assert(exist('rejection_rpqr', 'file') > 0, 'rejection_rpqr is not on the path.');
assert(exist('bsqr_rand_mex', 'file') == 3, 'bsqr_rand_mex is not built.');
if isempty(opt.outdir)
    opt.outdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'results');
end
if ~isfolder(opt.outdir); mkdir(opt.outdir); end

app   = {'gmm_kernel', 'integral_skeleton', 'snapshots'};   % approx_test_matrix
synth = {'gaussian', 'spiked_leverage', 'needle'};          % approx_synth_matrix
fprintf('approx comparison%s: ks=[%s], trials=%d, BLAS threads=%d\n', ...
    opt.tag, strjoin(string(opt.ks), ' '), opt.trials, maxNumCompThreads);

rows = {};
for fi = 1:numel(opt.families)
    fam = opt.families{fi};

    % --- build the (fixed) matrix; skip absent real families ---
    szfam = struct();
    if isfield(opt.sizes, fam); szfam = opt.sizes.(fam); end
    if any(strcmpi(fam, app))
        [A, in] = approx_test_matrix(fam, 20240601 + fi, szfam);
    elseif any(strcmpi(fam, synth))
        [A, in] = approx_synth_matrix(fam, 20240601 + fi, szfam);
    else
        try
            [A, in] = approx_real_matrix(fam);
        catch err
            fprintf('  [skip] real family "%s": %s\n', fam, err.message);
            continue;
        end
    end
    [m, n] = size(A);
    normA_fro = norm(A, 'fro');
    svals = svd(A);                             % exact spectrum for the optimal
    normA_2 = svals(1);                         % rank-k reference (no cancellation)
    fprintf('  %-18s %dx%d (%s)\n', fam, m, n, in.desc);

    % Sweep only meaningful k: k < rank cap min(m,n) (need a column to drop and a
    % sigma_{k+1}), and the optimal rank-k relative Frobenius error still above
    % rel_floor (else the spectrum has decayed to ~eps and there is nothing to
    % distinguish). Capping is per family because decay rates differ; report the
    % drop so the truncation is explicit.
    cand = opt.ks(opt.ks >= 1 & opt.ks < min(m, n));
    relopt = arrayfun(@(kk) sqrt(sum(svals(kk+1:end).^2)) / normA_fro, cand);
    ks = cand(relopt >= opt.rel_floor);
    dropped = cand(relopt < opt.rel_floor);
    if ~isempty(dropped)
        fprintf('    capped: k=[%s] dropped (optimal rel. Frob < %.0e)\n', ...
            strjoin(string(dropped), ' '), opt.rel_floor);
    end
    if isempty(ks)
        fprintf('    no k above rel_floor=%.0e; skipping family\n', opt.rel_floor);
        continue;
    end
    for k = ks
        % leading right singular vectors via Lanczos (svds) -- shared by both
        % selectors, so the comparison isolates column selection, not this
        % (accurate) subspace estimate.
        [~, ~, V] = svds(A, k);
        W = V.';                                % k-by-n, orthonormal rows
        opt_frob = sqrt(sum(svals(k+1:end).^2));   % exact optimal rank-k errors
        opt_spec = svals(k+1);

        for s = 1:opt.trials
            % randomized BSQR -- public defaults (batched, block=k, norm-weighted)
            t0 = tic;
            p = bsqr_rand(W, 'k', k, 'seed', s);
            tb = toc(t0);
            [fb, sb] = approx_err(A, p(1:k));
            mtb = interp_coeff_max(W, p(1:k));     % max|R11^{-1}R12| (basis quality)
            rows(end+1, :) = {fam, m, n, k, s, 'bsqr_rand', fb, sb, ...
                opt_frob, opt_spec, normA_fro, normA_2, tb, mtb}; %#ok<AGROW>

            % rejection_rpqr -- seeded for reproducibility
            rng(s);
            t0 = tic;
            idx = rejection_rpqr(W, k);
            tr = toc(t0);
            [fr, sr] = approx_err(A, idx(1:k));
            mtr = interp_coeff_max(W, idx(1:k));
            rows(end+1, :) = {fam, m, n, k, s, 'rejection_rpqr', fr, sr, ...
                opt_frob, opt_spec, normA_fro, normA_2, tr, mtr}; %#ok<AGROW>
        end
        fprintf('    k=%-4d done\n', k);
    end
end

if isempty(rows)
    warning('No families produced data (all skipped?).'); return;
end
T = cell2table(rows, 'VariableNames', {'family', 'm', 'n', 'k', 'seed', 'method', ...
    'frob_err', 'spec_err', 'opt_frob', 'opt_spec', 'normA_fro', 'normA_2', 'time_s', 'maxT'});
csv = fullfile(opt.outdir, ['exp_approx' opt.tag '.csv']);
writetable(T, csv);
fprintf('Wrote %s (%d rows)\n', csv, height(T));
end

% ---------------------------------------------------------------------------
function [fr, sp] = approx_err(A, sel)
% Orthogonal-projection error of approximating A by its selected columns:
% with Q_S = orth(A(:,sel)), the optimal approximation in that column span is
% P_S A = Q_S (Q_S' A), and the residual norms are method-agnostic. Uses no inv;
% orth handles a rank-deficient selection (duplicate/near-dependent columns)
% honestly by projecting onto the realized span. svds(E,1) gives the 2-norm
% cheaply (E is dense but the leading singular value needs only a few matvecs).
sel = sel(:).';
Qs = orth(A(:, sel));
E  = A - Qs * (Qs' * A);
fr = norm(E, 'fro');
sp = svds(E, 1);
if ~isfinite(sp); sp = norm(E); end            % rare svds non-convergence fallback
end

% ---------------------------------------------------------------------------
function mt = interp_coeff_max(W, sel)
% Interpolation-coefficient magnitude max|R11^{-1}R12| in the leading-k frame
% W = V_k' (k-by-n, orthonormal rows). W(:,sel) is k-by-k (square), so with its QR
% the coefficients are T = R11\R12, R12 = Q1'*W(:,rest) -- the same expression used
% by run_approx_cond_comparison. No inv(): a triangular solve. This tracks
% ||R11^{-1}|| (bounded for BSQR, uncontrolled for rejection) and is what an ID/CUR
% built on these columns would carry.
sel = sel(:).';
rest = setdiff(1:size(W, 2), sel);
[Q1, R11] = qr(W(:, sel), 0);
T = R11 \ (Q1' * W(:, rest));
mt = max(abs(T), [], 'all');
end
