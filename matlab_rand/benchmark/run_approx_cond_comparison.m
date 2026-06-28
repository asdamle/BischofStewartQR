function run_approx_cond_comparison(varargin)
%RUN_APPROX_COND_COMPARISON R11-conditioning companion: where the BSQR guarantee
%   shows up. randomized BSQR vs rejection_rpqr.
%
%   The approximation-error companions (run_approx_comparison / _synth) found both
%   methods near-optimal: the orthogonal-projection error depends only on the
%   *span* of the selected columns, so ||R11^{-1}|| is invisible to it. This
%   experiment measures the quantities that *do* depend on the selected basis -- the
%   ones a CUR / interpolative-decomposition pipeline actually pays for:
%
%     frobinv    = ||R11^{-1}||_F            (BSQR guarantees <= sqrt(k(n-k+1)))
%     maxT       = max |R11^{-1} R12|        (interpolation-coefficient magnitude)
%     proj_frob  = ||A - P_S A||_F / ||A||_F (orthogonal-projection error -- the
%                                             best fit in span(A(:,S)); ~identical
%                                             per method, the conditioning-blind
%                                             LOWER BOUND for any reconstruction)
%     id_err     = ||A - A(:,S)[I,T]||_F/||A||_F  (rank-k ID with the oblique
%                                             leading-k coefficients T: even
%                                             noiseless its error exceeds proj_frob
%                                             by a factor growing with ||R11^{-1}||)
%     noisy_id_err = same ID from *noisy* selected columns: the noise (sized to
%                                             noise_rel*||A(:,S)||) is amplified by
%                                             ||T|| ~ ||R11^{-1}|| on top of id_err
%     svdopt_frob = sqrt(sum sigma_{k+1:}^2)/||A||_F  (best-possible rank-k error,
%                                             SVD truncation -- the absolute lower
%                                             bound below proj_frob; from the exact
%                                             prescribed spectrum, no cancellation)
%
%   As a coefficient-choice control we also record the *standard projection* ID
%   coefficients T_proj = A(:,S)^+ A(:,rest) (k-step pivoted QR on A with the V_k^T
%   permutation), versus the leading-k-frame T above: maxTproj = max|T_proj|, and
%   noisy_id_err_proj / noisy_id_spec_proj (same noise, T_proj coefficients). With
%   T_proj the noiseless reconstruction equals the projection error exactly, so the
%   comparison shows whether switching coefficients changes max|.| and the noisy gap.
%
%   Every conditioning/error metric is recorded in BOTH norms: the Frobenius columns
%   above and their spectral (2-norm) counterparts specinv = ||R11^{-1}||_2,
%   proj_spec, id_spec, noisy_id_spec (each error / ||A||_2), with osinsky2 =
%   sqrt(1+k(n-k)) the spectral bound. The coefficient magnitudes maxT / maxTproj are
%   max-norm quantities in both. plot_approx_cond_comparison('norm','2') draws the
%   spectral version; the default draws the Frobenius one.
%
%   Here R11/R12 are taken in the leading-k right-singular-vector frame: W = V_k'
%   (k-by-n, orthonormal rows, from svds), W(:,S) is k-by-k, R12 = Q1' W(:,rest),
%   T = R11 \ R12 (no inv -- backslash solve; frobinv via norm(1./svd(R11)), spectral
%   residual norms via svds(.,1)).
%
%   Families are A = U diag(s) V' (approx_synth_matrix) with a fixed spectrum (a
%   sharp cliff near the middle of the swept k, so the rank-k ID-error rows have a
%   clear knee) and V from leverage families: 'gaussian' (benign control),
%   'spiked_leverage' (a few dominant high-leverage columns), 'collinear_cluster' (a
%   planted high-leverage near-collinear cluster). bsqr_rand runs with its public
%   defaults. Writes results/exp_approx_cond.csv; plot with
%   plot_approx_cond_comparison (linear k axis by default).
%
% Options: 'ks' (default 1:80, every k; the spectrum knee is at k~50),
%   'trials' (20), 'families'
%   (default {'gaussian','spiked_leverage','collinear_cluster'}), 'noise_rel' (default
%   1e-2, relative noise on the measured selected columns for noisy_id_err -- large
%   enough that the noise term competes with id_err; sweep it to vary the regime),
%   'sizes' (per-family struct forwarded to approx_synth_matrix), 'outdir'.

ip = inputParser;
addParameter(ip, 'ks', 1:80);          % every k up to 80 (plotted on a linear k axis)
addParameter(ip, 'trials', 20);
addParameter(ip, 'families', {'gaussian', 'spiked_leverage', 'collinear_cluster'});
addParameter(ip, 'noise_rel', 1e-2);
addParameter(ip, 'sizes', struct());
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

fprintf('approx cond comparison: ks=[%s], trials=%d, noise_rel=%.0e, BLAS threads=%d\n', ...
    strjoin(string(opt.ks), ' '), opt.trials, opt.noise_rel, maxNumCompThreads);

rows = {};
for fi = 1:numel(opt.families)
    fam = opt.families{fi};
    szfam = struct();
    if isfield(opt.sizes, fam); szfam = opt.sizes.(fam); end
    [A, in] = approx_synth_matrix(fam, 20240601 + fi, szfam);
    [m, n] = size(A);
    normA_fro = norm(A, 'fro');
    fprintf('  %-18s %dx%d rank %d\n', fam, m, n, in.r);

    ks = opt.ks(opt.ks >= 1 & opt.ks < min(m, n) & opt.ks < in.r);
    if isempty(ks); continue; end
    % One Lanczos SVD per family (the leading max(ks) right vectors); slice V(:,1:k)
    % per k. The leading-k subspace is the same as svds(A,k)'s, and every recorded
    % metric is invariant to right-vector column signs, so this is exact -- and far
    % cheaper than recomputing svds(A,k) for every k.
    [~, Sfull, Vfull] = svds(A, max(ks));
    normA_2 = Sfull(1, 1);                        % ||A||_2 = sigma_max(A) (k-independent)
    svals = in.svals(:);                          % exact prescribed spectrum (tail sums, no cancellation)
    for k = ks
        W = Vfull(:, 1:k).';                      % k-by-n, orthonormal rows = V_k'
        osinsky  = sqrt(k * (n - k + 1));        % Frobenius bound on ||R11^{-1}||_F
        osinsky2 = sqrt(1 + k * (n - k));        % spectral bound on ||R11^{-1}||_2
        % best-possible rank-k error (SVD truncation), relative -- the absolute lower
        % bound below any column-selection projection / ID error.
        svdopt_frob = sqrt(sum(svals(k+1:end) .^ 2)) / normA_fro;
        svdopt_spec = svals(k+1) / normA_2;
        for s = 1:opt.trials
            % randomized BSQR -- public defaults (batched, block=k, norm-weighted)
            p = bsqr_rand(W, 'k', k, 'seed', s);
            mb = cond_metrics(A, W, p(1:k), normA_fro, normA_2, opt.noise_rel, s);
            rows(end+1, :) = {fam, m, n, in.r, k, s, 'bsqr_rand', mb.frobinv, osinsky, ...
                mb.specinv, osinsky2, mb.maxT, mb.normT, mb.maxTproj, mb.proj_frob, mb.proj_spec, ...
                mb.id_err, mb.id_spec, mb.noisy_id_err, mb.noisy_id_err_proj, ...
                mb.noisy_id_spec, mb.noisy_id_spec_proj, opt.noise_rel, ...
                svdopt_frob, svdopt_spec}; %#ok<AGROW>

            % rejection_rpqr -- seeded for reproducibility
            rng(s);
            idx = rejection_rpqr(W, k);
            mr = cond_metrics(A, W, idx(1:k), normA_fro, normA_2, opt.noise_rel, s);
            rows(end+1, :) = {fam, m, n, in.r, k, s, 'rejection_rpqr', mr.frobinv, osinsky, ...
                mr.specinv, osinsky2, mr.maxT, mr.normT, mr.maxTproj, mr.proj_frob, mr.proj_spec, ...
                mr.id_err, mr.id_spec, mr.noisy_id_err, mr.noisy_id_err_proj, ...
                mr.noisy_id_spec, mr.noisy_id_spec_proj, opt.noise_rel, ...
                svdopt_frob, svdopt_spec}; %#ok<AGROW>
        end
        fprintf('    k=%-4d done\n', k);
    end
end

if isempty(rows); warning('No families produced data.'); return; end
T = cell2table(rows, 'VariableNames', {'family', 'm', 'n', 'r', 'k', 'seed', 'method', ...
    'frobinv', 'osinsky', 'specinv', 'osinsky2', 'maxT', 'normT', 'maxTproj', ...
    'proj_frob', 'proj_spec', 'id_err', 'id_spec', ...
    'noisy_id_err', 'noisy_id_err_proj', 'noisy_id_spec', 'noisy_id_spec_proj', 'noise_rel', ...
    'svdopt_frob', 'svdopt_spec'});
csv = fullfile(opt.outdir, 'exp_approx_cond.csv');
writetable(T, csv);
fprintf('Wrote %s (%d rows)\n', csv, height(T));
end

% ---------------------------------------------------------------------------
function out = cond_metrics(A, W, S, normA_fro, normA_2, noise_rel, seed)
% Basis-dependent quality of the selection S (k columns), in the leading-k
% right-singular-vector frame W = V_k'. Every error is reported in both the
% Frobenius (_frob/_err) and spectral (_spec) norms so the figure can be drawn in
% either norm; the coefficient magnitude maxT stays a max-norm quantity. No inv():
% conditioning via svd, interpolation coefficients via a triangular solve, and the
% residual 2-norms via svds(.,1).
S = S(:).';
n = size(W, 2); k = numel(S);
rest = setdiff(1:n, S);

[Q1, R11] = qr(W(:, S), 0);              % W(:,S) is k-by-k (square) => Q1 orthogonal
R12 = Q1' * W(:, rest);                  % so [R11 R12] = Q1' W(:,[S rest])
sv = svd(R11);
out.frobinv = norm(1 ./ sv);             % ||R11^{-1}||_F = ||W(:,S)^{-1}||_F
out.specinv = max(1 ./ sv);              % ||R11^{-1}||_2 = 1/sigma_min(R11)
T = R11 \ R12;                           % = W(:,S)^{-1} W(:,rest), interp. coeffs
out.maxT = max(abs(T), [], 'all');       % max-norm coefficient magnitude (both plots)
out.normT = norm(T, 'fro');

% orthogonal-projection error + the *standard* (least-squares / projection) ID
% coefficients T_proj = A(:,S)^+ A(:,rest), both from one QR of A(:,S). The
% projection error is the conditioning-blind LOWER BOUND; with T_proj the noiseless
% reconstruction equals it exactly (the oblique penalty of the V_k-frame T is gone),
% so only the noisy_*_proj errors are recorded separately. No inv(): triangular solve.
[Qa, Ra] = qr(A(:, S), 0);              % A(:,S) = Qa*Ra (Qa m-by-k, Ra k-by-k)
QtA = Qa' * A;                          % k-by-n; QtA(:,rest) = R12 of A
Eproj = A - Qa * QtA;
out.proj_frob = norm(Eproj, 'fro') / normA_fro;
out.proj_spec = snorm2(Eproj) / normA_2;
Tproj = Ra \ QtA(:, rest);             % = A(:,S)^+ A(:,rest), the projection coeffs
out.maxTproj = max(abs(Tproj), [], 'all');

% rank-k interpolative decomposition A ~ A(:,S) [I, T]. The coefficients T are the
% *oblique* leading-k-frame coefficients (not the least-squares ones), so even with
% no noise this reconstruction's error exceeds the projection error by a factor that
% grows with ||R11^{-1}|| -- the conditioning penalty. id_* isolates that effect.
Ahat = zeros(size(A));
Ahat(:, S) = A(:, S);
Ahat(:, rest) = A(:, S) * T;
Eid = A - Ahat;
out.id_err  = norm(Eid, 'fro') / normA_fro;
out.id_spec = snorm2(Eid) / normA_2;

% same ID rebuilt from *noisy* measured selected columns: Ahat(:,rest) = C*T with
% C = A(:,S)+E. The noise (sized to noise_rel*||A(:,S)||) propagates as E*T, adding a
% term that scales with ||T|| ~ ||R11^{-1}|| on top of id_*, so an ill-conditioned
% R11 amplifies measurement noise too.
rng(seed + 777);
E = randn(size(A, 1), k);
E = E * (noise_rel * norm(A(:, S), 'fro') / norm(E, 'fro'));
C = A(:, S) + E;
Ahat(:, S) = C;
Ahat(:, rest) = C * T;                  % V_k-frame coefficients
En = A - Ahat;
out.noisy_id_err  = norm(En, 'fro') / normA_fro;
out.noisy_id_spec = snorm2(En) / normA_2;
Ahat(:, rest) = C * Tproj;              % projection coefficients, same noise draw
Enp = A - Ahat;
out.noisy_id_err_proj  = norm(Enp, 'fro') / normA_fro;
out.noisy_id_spec_proj = snorm2(Enp) / normA_2;
end

% ---------------------------------------------------------------------------
function s = snorm2(E)
% Spectral norm (largest singular value) of a dense residual, cheaply via svds;
% fall back to the full norm on the rare non-convergence.
s = svds(E, 1);
if ~isfinite(s); s = norm(E); end
end
