function tests = test_mfile_mex_fuzz
%TEST_MFILE_MEX_FUZZ Seeded randomized fuzz of the two MATLAB backends.
%   bsqr dispatches to the pure-MATLAB kernel (private/bsqr_mfile.m) and the
%   C++ MEX (mex/src/bsqr_mex.cpp). They implement the same unblocked
%   Bischof-Stewart algorithm -- same pivot criterion (1+||w_j||^2)/||a_j||^2
%   with strict first-minimum tie-break, same sqrt(eps) norm-recompute guard
%   -- and must agree exactly. This test hammers that contract over many
%   seeded random matrices spanning shapes (tall/square/short-wide) and
%   conditionings (Gaussian, clustered/decaying spectra, exact
%   rank-deficiency, and deliberately TIED column norms that exercise the
%   strict-< tie-break), asserting:
%     a. identical pivot sequences (integer equality) over the numerical-rank
%        block and R11 agreement to a tight bound; rinv_r12 (== R11^{-1}R12,
%        the derived W state) is compared cross-backend on the full-rank
%        wide cases where it is well-defined, under a cond(R11)-aware bound;
%     b. reconstruction residual and Q-orthogonality at machine precision
%        for each backend independently.
%   A couple of cancellation-prone cases (near-dependent columns forcing
%   norm downdates through the recompute threshold) complement
%   test_cancellation_stress.m.
%
%   The MEX defaults to a panel/blocked kernel above a k*n crossover; the
%   mfile has no panel path. To compare like for like we pin the MEX to its
%   unblocked kernel (BS_PANEL_NB=0) for the duration -- the meaningful
%   "same algorithm, two implementations" check. Under that pinning, both
%   backends evaluate the identical criterion recurrence in the identical
%   order, so exact pivot equality is a fair demand across engineered
%   near-ties (equal *nonzero* column norms): a disagreement there is a real
%   backend bug, not tie noise.
%
%   Exact equality is only demanded over the leading NUMERICAL-RANK steps.
%   Past the numerical rank every remaining working column is pure roundoff
%   (running norm collapsed to ~ (eps*||A||)^2, criterion ~ 1e30+), and
%   which roundoff-noise column "wins" the argmin legitimately depends on
%   BLAS summation order -- the interpreted mfile and the compiled MEX may
%   differ there with neither being wrong. This mirrors the rest of the
%   suite, where oracle_bsqr rejects exact rank-deficiency outright and the
%   parity zoo is screened for near-ties. On the roundoff block we still
%   require each backend independently to stay finite, reconstruct A(:,p),
%   and yield an orthogonal Q.
%
%   If the MEX is not built, every case still runs its mfile self-checks
%   (so the file is never a silent no-op); only the cross-backend comparison
%   is skipped, and each test method records that skip once via assumeTrue,
%   so run_tests reports Incomplete rather than a spurious hard failure.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab'));
testCase.TestData.repo_root = repo_root;
end

function setup(testCase)
% Pin the MEX to its unblocked kernel for every case so it runs the same
% algorithm as the mfile (which is unblocked only). Restored per-case.
old_nb = getenv('BS_PANEL_NB');
testCase.TestData.restore_nb = onCleanup(@() setenv('BS_PANEL_NB', old_nb));
setenv('BS_PANEL_NB', '0');
end

% ----------------------------------------------------------------------
% Main fuzz sweep: shapes x conditionings x seeds.
% ----------------------------------------------------------------------

function testShapeConditioningFuzz(testCase)
% Kept small so k*n stays well under the panel crossover and the run is
% quick, but broad: 3 shapes x 4 families x several seeds each.
shapes = { ...
    'tall',      [60, 24]; ...
    'square',    [40, 40]; ...
    'shortwide', [24, 72]};
families = {'gaussian', 'clustered', 'decaying', 'rankdef'};
seeds = [11, 101, 4242, 99991];

for si = 1:size(shapes, 1)
    shape_name = shapes{si, 1};
    sz = shapes{si, 2};
    m = sz(1); n = sz(2);
    for fi = 1:numel(families)
        fam = families{fi};
        for seed = seeds
            A = make_matrix(fam, m, n, seed);
            label = sprintf('%s/%s/seed%d [%dx%d]', shape_name, fam, seed, m, n);
            compare_backends(testCase, A, min(m, n), label);
        end
    end
end
note_if_mex_absent(testCase);
end

% ----------------------------------------------------------------------
% Tied column norms: exercises the strict-< first-minimum tie-break. When
% many columns share an identical criterion value, both backends must break
% the tie to the same (earliest) column, step after step.
% ----------------------------------------------------------------------

function testTiedColumnNorms(testCase)
% (1) Exactly equal columns: an all-ones block plus a duplicated Gaussian
% block. Every criterion is a hard tie at step 1 and stays tie-prone as
% reduction proceeds.
rng(20260630, 'twister');
G = randn(30, 6);
A1 = [ones(30, 6), G, G, 2 * G];
compare_backends(testCase, A1, min(size(A1)), 'tied/equal_and_duplicate_cols');

% (2) Blockwise-equal column norms with orthogonal directions: identical
% ||a_j|| across a block but different content, so early ties resolve into
% a nontrivial but still tie-studded sequence.
rng(777, 'twister');
[U, ~] = qr(randn(48, 48));
blk = 8;
A2 = zeros(48, 48);
for b = 1:6
    cols = (b - 1) * blk + (1:blk);
    A2(:, cols) = U(:, cols) * b;   % each block has constant column norm b
end
% Permute columns with a fixed shuffle so equal-norm columns are
% interleaved, stressing first-minimum selection order.
rng(778, 'twister');
A2 = A2(:, randperm(48));
compare_backends(testCase, A2, min(size(A2)), 'tied/blockwise_equal_norms');

% (3) Integer matrix with repeated rows/cols: exact ties with zero
% floating error in the initial norms.
A3 = double([1 1 2 2 3 3; 2 2 1 1 3 3; 3 3 3 3 1 1; 1 1 1 1 1 1]);
A3 = [A3; A3];   % 8x6, guaranteed exact tie structure
compare_backends(testCase, A3, min(size(A3)), 'tied/integer_repeats');
note_if_mex_absent(testCase);
end

% ----------------------------------------------------------------------
% Cancellation-prone cases beyond test_cancellation_stress.m: columns that
% become near-dependent so the running norm s(j) decays through the
% sqrt(eps) recompute threshold, triggering exact-tail refreshes that both
% backends must handle identically.
% ----------------------------------------------------------------------

function testCancellationDowndates(testCase)
% (1) A near-rank-1 backbone with tiny orthogonal perturbations: after the
% dominant direction is pivoted out, trailing running norms collapse toward
% the recompute floor.
rng(20260629, 'twister');
m = 40; n = 30;
u = randn(m, 1); u = u / norm(u);
v = randn(n, 1);
delta = 1e-9;
A1 = u * v' + delta * randn(m, n);
compare_backends(testCase, A1, min(m, n), 'cancel/near_rank1_backbone');

% (2) Graded spectrum whose columns are then near-duplicated: the update to
% wnorm2/s subtracts near-equal quantities, maximizing catastrophic
% cancellation in the recurrence right at the guard threshold.
rng(20260628, 'twister');
r = 20; m = 44;
[Uf, ~] = qr(randn(m, m));
Ur = Uf(:, 1:r);
[Vr, ~] = qr(randn(r, r));
B = Ur * diag(logspace(0, -9, r)) * Vr';
tail = Uf(:, r+1:end) * randn(m - r, r);
tail = tail ./ vecnorm(tail);
A2 = [B, B + 1e-7 * tail];   % 44x40 near-duplicate pair per column
compare_backends(testCase, A2, min(size(A2)), 'cancel/graded_near_duplicate');

% (3) Progressive collinearity: each column is the previous plus a shrinking
% orthogonal kick, so consecutive columns are increasingly parallel.
rng(20260627, 'twister');
m = 36; n = 28;
Q0 = orth(randn(m, n));
A3 = zeros(m, n);
A3(:, 1) = Q0(:, 1);
for j = 2:n
    A3(:, j) = A3(:, j - 1) + (0.5 ^ j) * Q0(:, j);
end
compare_backends(testCase, A3, min(m, n), 'cancel/progressive_collinear');
note_if_mex_absent(testCase);
end

% ----------------------------------------------------------------------
% Early-stopping sweep: k < min(m,n) changes the trailing-block bookkeeping;
% the backends must still agree on the truncated pivots and rinv_r12.
% ----------------------------------------------------------------------

function testEarlyStopFuzz(testCase)
cases = { ...
    [50, 40], 12; ...
    [40, 40], 20; ...
    [24, 80], 10; ...
    [60, 30], 30};
for ci = 1:size(cases, 1)
    sz = cases{ci, 1};
    k = cases{ci, 2};
    m = sz(1); n = sz(2);
    for seed = [5, 55, 555]
        A = make_matrix('gaussian', m, n, seed);
        label = sprintf('earlystop/k=%d/seed%d [%dx%d]', k, seed, m, n);
        compare_backends(testCase, A, k, label);
    end
end
note_if_mex_absent(testCase);
end

% ----------------------------------------------------------------------
% Core comparison: run both backends, assert exact pivot equality, tight
% agreement on R and rinv_r12, and machine-precision residual/orthogonality
% for each backend independently.
% ----------------------------------------------------------------------

function compare_backends(testCase, A, k, label)
[m, n] = size(A);

% mfile backend always runs, for EVERY case: the mfile self-checks are never
% skipped, so the file is never a silent no-op even without the MEX.
[Qm, Rm, pm, Wm] = bsqr(A, 'backend', 'mfile', 'k', k, ...
    'pivot_format', 'vector', 'return_rinv_r12', true);
check_single_backend(testCase, A, Qm, Rm, pm, k, ['mfile ', label]);

if ~bsqr_mex_available()
    % Skip only the cross-backend comparison, and do it with `return` rather
    % than assumeFail so the remaining cases in the caller's loop still run
    % their mfile self-checks. The method-level assumeFail (checked once in
    % each test method) records the skip so run_tests reports Incomplete, not
    % a spurious failure.
    return;
end

[Qx, Rx, px, Wx] = bsqr(A, 'backend', 'mex', 'k', k, ...
    'pivot_format', 'vector', 'return_rinv_r12', true);
check_single_backend(testCase, A, Qx, Rx, px, k, ['mex ', label]);

% Both backends always agree exactly on shape.
verifyEqual(testCase, size(Rx), size(Rm), ...
    sprintf('%s: R size mex vs mfile', label));
verifyEqual(testCase, size(Wx), size(Wm), ...
    sprintf('%s: rinv_r12 size mex vs mfile', label));

% Exact equality is a fair demand only over the leading numerical-rank
% steps; past that the working columns are roundoff and the argmin is
% BLAS-order dependent (see file header). kc == k for full-rank inputs.
kc = min(k, numerical_rank(A));

% (a) Exact pivot equality on the meaningful block, integer for integer.
verifyEqual(testCase, px(1:kc), pm(1:kc), ...
    sprintf('%s: pivot sequence (first %d steps) mex vs mfile', label, kc));

% (a) R agrees to a tight tolerance on the meaningful leading triangle.
% Both are the same arithmetic in the same order there, so the residual
% scales with matrix magnitude, not conditioning; a snug bound is warranted.
tolRR = 1e-11;
if kc > 0
    % Rows/cols 1:kc carry the leading triangle; the trailing columns move
    % under the shared permutation, already pinned by the exact pivot check.
    verifyLessThan(testCase, rel_err(Rx(1:kc, 1:kc), Rm(1:kc, 1:kc)), tolRR, ...
        sprintf('%s: R11 factor mex vs mfile', label));
end

% rinv_r12 (== R11^{-1}R12 for the k-split) is only well-defined when the
% whole k-block is within numerical rank (R11 nonsingular). Compare the full
% k-by-(n-k) block then; skip it when the factorization runs into the
% roundoff tail (kc < k), where R11 is singular and the solve is noise.
%
% Unlike R (a direct kernel output, agreeing at ~1e-15), rinv_r12 is a
% *derived* triangular solve: its backend-to-backend difference is the
% forward error of that solve and is amplified by cond(R11). Both backends
% compute the same solve but sum in different orders (interpreted mfile vs.
% compiled MEX), so on clustered/graded spectra the leading digits still
% match but the trailing ones need not. The condition-amplification bound is
% eps*cond(R11); measured disagreement runs ~0.05x of it, so 10*eps*cond
% gives generous headroom while staying snug. A hard ceiling of 1e-6 keeps
% the bound from ever going vacuous on the steepest spectra (decaying,
% cond ~ 2e9) -- well above the ~2e-8 actually observed there, but tight
% enough to still catch a gross W/rinv-recurrence divergence. This is the
% same condition sensitivity that makes parity_zoo disable rinv_check on
% graded members outright rather than compare it.
if kc == k && k > 0 && k < n
    condR11 = cond(Rm(1:k, 1:k));
    tol_rinv = min(1e-6, max(tolRR, 10 * eps('double') * condR11));
    verifyLessThan(testCase, rel_err(Wx, Wm), tol_rinv, ...
        sprintf('%s: rinv_r12 mex vs mfile (cond(R11)=%.2e)', label, condR11));
end
end

function note_if_mex_absent(testCase)
% Called once at the end of each test method. When the MEX is not built the
% mfile self-checks above have already run for every case (so the file is
% never a silent no-op); this records the missing cross-backend comparison
% as an assumption failure -> the method reports Incomplete (treated as
% success by assertSuccess), never a spurious hard failure.
assumeTrue(testCase, bsqr_mex_available(), ...
    'bsqr_mex not built: mfile self-checks ran, but the mfile-vs-mex comparison was skipped (run build_bsqr_mex to exercise it).');
end

function r = numerical_rank(A)
% Numerical rank via the standard relative singular-value threshold. This
% is the count of columns that carry content above roundoff; beyond it the
% Bischof-Stewart working norms have collapsed and pivot selection among the
% residual roundoff is legitimately implementation-defined.
s = svd(A);
if isempty(s) || s(1) == 0
    r = 0;
    return;
end
tol = max(size(A)) * eps(s(1));
r = sum(s > tol);
end

function check_single_backend(testCase, A, Q, R, p, k, label)
[m, n] = size(A);
verifySize(testCase, Q, [m, k], sprintf('%s: Q size', label));
verifySize(testCase, R, [k, n], sprintf('%s: R size', label));
verifyEqual(testCase, sort(p(:)'), 1:n, sprintf('%s: p is a permutation', label));

% (b) Reconstruction residual at machine precision: A(:,p) == Q*R for a
% full-rank factorization (k == min(m,n)); the criterion never selects a
% numerically-zero pivot before then, so this holds even on our
% rank-deficient / near-dependent inputs.
if k == min(m, n)
    resid = rel_resid(A(:, p), Q * R);
    verifyLessThan(testCase, resid, recon_tol([m, n]), ...
        sprintf('%s: reconstruction residual', label));
end

% (b) Q-orthogonality at machine precision.
verifyLessThan(testCase, orth_err(Q), recon_tol([m, n]), ...
    sprintf('%s: Q orthogonality', label));
end

% ----------------------------------------------------------------------
% Matrix builders (self-contained; not coupled to parity_zoo / matlab_rand).
% ----------------------------------------------------------------------

function A = make_matrix(family, m, n, seed)
rng(seed, 'twister');
switch family
    case 'gaussian'
        A = randn(m, n);
    case 'clustered'
        % Singular values in a few tight clusters (near-multiplicities).
        r = min(m, n);
        sv = build_clustered_spectrum(r);
        A = random_matrix_with_spectrum(m, n, sv);
    case 'decaying'
        % Smoothly decaying (geometric) spectrum down to ~1e-9.
        r = min(m, n);
        sv = logspace(0, -9, r)';
        A = random_matrix_with_spectrum(m, n, sv);
    case 'rankdef'
        % Exact rank deficiency: rank ~ floor(min(m,n)/2), remaining
        % singular values exactly zero. Production kernels must pivot the
        % nonzero block first and stay finite; reconstruction still holds
        % at full k. Backends are compared exactly only over that nonzero
        % block (see compare_backends / numerical_rank).
        r = min(m, n);
        rk = max(1, floor(r / 2));
        sv = zeros(r, 1);
        sv(1:rk) = logspace(0, -2, rk)';
        A = random_matrix_with_spectrum(m, n, sv);
    otherwise
        error('test_mfile_mex_fuzz:UnknownFamily', 'family=%s', family);
end
end

function sv = build_clustered_spectrum(r)
% Three clusters at 1e0, 1e-3, 1e-6 with tiny jitter -> near-multiplicities.
centers = [1e0, 1e-3, 1e-6];
sv = zeros(r, 1);
for idx = 1:r
    c = centers(mod(idx - 1, numel(centers)) + 1);
    sv(idx) = c * (1 + 1e-10 * randn());
end
sv = sort(sv, 'descend');
end

function A = random_matrix_with_spectrum(m, n, sv)
% Build A = U * diag(sv) * V' with Haar-random U, V (economy).
r = min(m, n);
sv = sv(1:r);
[U, ~] = qr(randn(m, r), 0);
[V, ~] = qr(randn(n, r), 0);
A = U * (sv .* V');
end

% ----------------------------------------------------------------------
% Metrics.
% ----------------------------------------------------------------------

function r = rel_err(X, Y)
r = norm(X - Y, 'fro') / max(norm(Y, 'fro'), eps('double'));
end

function r = rel_resid(X, Y)
r = norm(X - Y, 'fro') / max(norm(X, 'fro'), eps('double'));
end

function o = orth_err(Q)
if isempty(Q)
    o = 0;
else
    o = norm(eye(size(Q, 2)) - Q' * Q, 'fro');
end
end

function t = recon_tol(sz)
% Machine-precision residual/orthogonality bound, scaled by problem size
% (matches the scaling used across the existing MATLAB suite).
t = 8e2 * eps('double') * max(sz);
end
