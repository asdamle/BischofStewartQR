function tests = test_bsqr_rand_bounds
%TEST_BSQR_RAND_BOUNDS The headline guarantee: ||R11^{-1}||_F stays bounded.
%   These assertions hold for orthonormal-row input (the GKS setting).
tests = functiontests(localfunctions);
end

function M = orthonormal_rows(k, n, seed)
rng(seed);
M = orth(randn(n, k))';
end

function verify_bound(testCase, stats, k, n)
% Per-step running bound and the final Osinsky bound.
tol = 1e-8 * max(1, stats.Fhat(end));
verifyLessThanOrEqual(testCase, stats.f2(:), stats.Fhat(:) + tol, ...
    'running f2 must stay under the per-step worst-case Fhat');
osinsky2 = k * (n - k + 1);
verifyLessThanOrEqual(testCase, stats.f2(end), osinsky2 * (1 + 1e-8), ...
    'final ||R11^{-1}||_F^2 must respect Osinsky''s bound');
end

function testRunningMeanBound(testCase)
cases = {[8, 50], [16, 120], [32, 300], [24, 400]};
for c = 1:numel(cases)
    k = cases{c}(1); n = cases{c}(2);
    M = orthonormal_rows(k, n, 200 + c);
    [~, ~, ~, stats] = bsqr_rand(M, 'backend', 'mfile', 'seed', c, ...
        'threshold_mode', 'running_mean');
    verify_bound(testCase, stats, k, n);
end
end

function testWorstcaseAllowanceBound(testCase)
cases = {[8, 50], [16, 120], [32, 300], [24, 400]};
for c = 1:numel(cases)
    k = cases{c}(1); n = cases{c}(2);
    M = orthonormal_rows(k, n, 300 + c);
    [~, ~, ~, stats] = bsqr_rand(M, 'backend', 'mfile', 'seed', c, ...
        'threshold_mode', 'worstcase_allowance');
    verify_bound(testCase, stats, k, n);
end
end

function testF2MatchesR11Inverse(testCase)
% The implicitly tracked f2 must equal the explicit ||inv(R11)||_F^2.
k = 20; n = 140;
M = orthonormal_rows(k, n, 55);
[~, ~, R11, stats] = bsqr_rand(M, 'backend', 'mfile', 'seed', 6);
explicit = sum(1 ./ svd(R11).^2);   % ||R11^{-1}||_F^2 via singular values
verifyLessThan(testCase, abs(stats.f2(end) - explicit) / explicit, 1e-8);
end

function testSamplesReported(testCase)
k = 16; n = 200;
M = orthonormal_rows(k, n, 71);
[~, ~, ~, stats] = bsqr_rand(M, 'backend', 'mfile', 'seed', 8);
verifyGreaterThanOrEqual(testCase, stats.total_tested, k);
verifyGreaterThanOrEqual(testCase, stats.samples_tested(:), ones(k, 1));
end

function testComparableToDeterministic(testCase)
% Conditioning should be in the same ballpark as the deterministic kernel
% (not identical pivots). Compare ||R11^{-1}||_F.
k = 20; n = 160;
M = orthonormal_rows(k, n, 88);
[~, ~, ~, stats] = bsqr_rand(M, 'backend', 'mfile', 'seed', 2);
rand_frobinv = stats.frob_inv;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
det_dir = fullfile(repo_root, 'matlab');
if exist(fullfile(det_dir, 'bsqr.m'), 'file')
    addpath(det_dir);
    [~, ~, pvec, rinv] = bsqr(M, 'k', k, 'backend', 'mfile', ...
        'pivot_format', 'vector', 'return_rinv_r12', true); %#ok<ASGLU>
    Rdet = bsqr(M, 'k', k, 'backend', 'mfile');
    det_frobinv = norm(1 ./ svd(Rdet(1:k, 1:k)));   % ||R11^{-1}||_F via singular values
    % Randomized should be within a small constant factor of deterministic.
    verifyLessThan(testCase, rand_frobinv, 3 * det_frobinv + 1e-8);
end
% Always: respects the absolute Osinsky guarantee.
verifyLessThanOrEqual(testCase, rand_frobinv, sqrt(k * (n - k + 1)) * (1 + 1e-8));
end
