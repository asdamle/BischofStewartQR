function tests = test_bsqr_rand
%TEST_BSQR_RAND Correctness tests for the randomized BSQR variant.
tests = functiontests(localfunctions);
end

% ----- helpers -------------------------------------------------------------

function M = orthonormal_rows(k, n, seed)
% k-by-n matrix with orthonormal rows (M*M' = I_k), the GKS setting.
rng(seed);
M = orth(randn(n, k))';
end

function check_factorization(testCase, M, p, reflectors, R11, k)
[m, n] = size(M);
verifyEqual(testCase, sort(p), 1:n, 'p must be a permutation of 1:n');
verifySize(testCase, R11, [k, k]);
verifyEqual(testCase, R11, triu(R11), 'AbsTol', 0, 'R11 must be upper triangular');

if k == 0
    return;
end
Qk = bsqr_rand_formQ(reflectors, k);          % m-by-k
sel = p(1:k);
scale = max(1, norm(M(:, sel), 'fro'));
% Q' * A(:,sel) = R11
verifyLessThan(testCase, norm(Qk' * M(:, sel) - R11, 'fro') / scale, 1e-10);
% A(:,sel) = Q(:,1:k) * R11
verifyLessThan(testCase, norm(M(:, sel) - Qk * R11, 'fro') / scale, 1e-10);
% columns of Qk orthonormal
verifyLessThan(testCase, norm(Qk' * Qk - eye(k), 'fro'), 1e-10);
end

% ----- tests ---------------------------------------------------------------

function testReconstructionOrthonormal(testCase)
cases = {[8, 40], [16, 64], [32, 200], [10, 10]};
for c = 1:numel(cases)
    k = cases{c}(1); n = cases{c}(2);
    M = orthonormal_rows(k, n, 100 + c);
    [p, reflectors, R11] = bsqr_rand(M, 'backend', 'mfile', 'seed', c);
    check_factorization(testCase, M, p, reflectors, R11, k);
end
end

function testGeneralMatrix(testCase)
% Theory assumes orthonormal rows, but the factorization must still be exact.
rng(7);
M = randn(20, 90);
k = 12;
[p, reflectors, R11] = bsqr_rand(M, 'k', k, 'backend', 'mfile', 'seed', 3);
check_factorization(testCase, M, p, reflectors, R11, k);
end

function testReturnR12(testCase)
k = 12; n = 60;
M = orthonormal_rows(k, n, 11);
[p, reflectors, R11, ~, R12] = bsqr_rand(M, 'backend', 'mfile', ...
    'seed', 5, 'return_r12', true);
check_factorization(testCase, M, p, reflectors, R11, k);
Qk = bsqr_rand_formQ(reflectors, k);
verifySize(testCase, R12, [k, n - k]);
verifyLessThan(testCase, norm(Qk' * M(:, p(k+1:n)) - R12, 'fro'), 1e-10);
end

function testR12RequiresOptIn(testCase)
M = orthonormal_rows(6, 30, 2);
verifyError(testCase, @() call_with_five_outputs(M), 'bsqr_rand:R12NotRequested');
end

function out = call_with_five_outputs(M)
[~, ~, ~, ~, out] = bsqr_rand(M, 'backend', 'mfile');
end

function testEdgeKZero(testCase)
M = orthonormal_rows(5, 20, 1);
[p, reflectors, R11, stats] = bsqr_rand(M, 'k', 0, 'backend', 'mfile');
verifyEqual(testCase, sort(p), 1:20);
verifySize(testCase, R11, [0, 0]);
verifySize(testCase, reflectors.V, [5, 0]);
verifyEqual(testCase, stats.total_tested, 0);
end

function testEdgeKFull(testCase)
% Square orthonormal matrix: selecting all columns.
rng(21);
M = orth(randn(15, 15));
k = 15;
[p, reflectors, R11] = bsqr_rand(M, 'k', k, 'backend', 'mfile', 'seed', 1);
check_factorization(testCase, M, p, reflectors, R11, k);
end

function testDeterminismWithSeed(testCase)
M = orthonormal_rows(16, 80, 9);
[p1, ~, R11a, s1] = bsqr_rand(M, 'backend', 'mfile', 'seed', 42);
[p2, ~, R11b, s2] = bsqr_rand(M, 'backend', 'mfile', 'seed', 42);
verifyEqual(testCase, p1, p2);
verifyEqual(testCase, R11a, R11b);
verifyEqual(testCase, s1.samples_tested, s2.samples_tested);
end

function testThresholdAndSamplingKnobs(testCase)
M = orthonormal_rows(20, 120, 13);
opts = { ...
    {'threshold_mode', 'running_mean', 'sampling', 'uniform', 'pick', 'best_in_block'}, ...
    {'threshold_mode', 'worstcase_allowance', 'sampling', 'uniform', 'pick', 'first'}, ...
    {'threshold_mode', 'running_mean', 'sampling', 'normweighted', 'pick', 'best_in_block'}, ...
    {'threshold_mode', 'worstcase_allowance', 'sampling', 'normweighted', 'pick', 'first'}};
for c = 1:numel(opts)
    [p, reflectors, R11] = bsqr_rand(M, 'backend', 'mfile', 'seed', c, opts{c}{:});
    check_factorization(testCase, M, p, reflectors, R11, 20);
end
end

function testMexAgreesOnInvariants(testCase)
% MEX uses its own RNG, so pivots need not match the m-file; both must
% nonetheless produce a valid factorization and respect the bound. Exercise
% both sampling schemes (uniform / Fenwick-based norm-weighted) and both
% threshold modes.
if ~bsqr_rand_mex_available()
    return;  % nothing to check until the MEX is built
end
k = 24; n = 150;
M = orthonormal_rows(k, n, 17);
for sampling = ["uniform", "normweighted"]
    for mode = ["running_mean", "worstcase_allowance"]
        for bs = [1, 7, 32]
            [p, reflectors, R11, stats] = bsqr_rand(M, 'backend', 'mex', 'seed', 4, ...
                'sampling', char(sampling), 'threshold_mode', char(mode), 'block_size', bs);
            check_factorization(testCase, M, p, reflectors, R11, k);
            tol = 1e-9 * max(1, stats.Fhat(end));
            verifyLessThanOrEqual(testCase, stats.f2(:), stats.Fhat(:) + tol);
        end
    end
end
end

function testMexWeightedConcentratedLeverage(testCase)
% On concentrated-leverage input the Fenwick weighted sampler must still
% produce an exact factorization and respect the bound.
if ~bsqr_rand_mex_available()
    return;
end
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab_rand', 'benchmark'));
k = 16; n = 400;
M = rand_test_matrix('needle', k, n, 23);
[p, reflectors, R11, stats] = bsqr_rand(M, 'backend', 'mex', 'seed', 2, ...
    'sampling', 'normweighted');
check_factorization(testCase, M, p, reflectors, R11, k);
verifyLessThanOrEqual(testCase, stats.f2(end), k * (n - k + 1) * (1 + 1e-8));
end

function testMexReturnR12(testCase)
% The MEX R12 path (compact-WY apply of Q' to the leftover columns) must match
% Q(:,1:k)' * A(:,P(k+1:n)).
if ~bsqr_rand_mex_available()
    return;
end
k = 14; n = 80;
M = orthonormal_rows(k, n, 19);
[p, reflectors, R11, ~, R12] = bsqr_rand(M, 'backend', 'mex', 'seed', 7, 'return_r12', true);
check_factorization(testCase, M, p, reflectors, R11, k);
Qk = bsqr_rand_formQ(reflectors, k);
verifySize(testCase, R12, [k, n - k]);
verifyLessThan(testCase, norm(Qk' * M(:, p(k+1:n)) - R12, 'fro') / max(1, norm(R12, 'fro')), 1e-10);
end

function testMexPickFirst(testCase)
% pick='first' on the MEX must still factor exactly and respect the bound.
if ~bsqr_rand_mex_available()
    return;
end
k = 16; n = 100;
M = orthonormal_rows(k, n, 23);
[p, reflectors, R11, stats] = bsqr_rand(M, 'backend', 'mex', 'seed', 3, 'pick', 'first');
check_factorization(testCase, M, p, reflectors, R11, k);
tol = 1e-9 * max(1, stats.Fhat(end));
verifyLessThanOrEqual(testCase, stats.f2(:), stats.Fhat(:) + tol);
end

function testRankDeficientErrors(testCase)
% k beyond the (exact) rank: only k-1 nonzero columns, the rest literally zero,
% so at step k every remaining residual is exactly zero and the rank guard must
% fire rather than propagate Inf into f2/R11. (Near-rank-deficiency where the
% residual is ~eps rather than exactly 0 is a documented precondition, not a
% guarded case.)
rng(31);
k = 6; n = 40;
M = [randn(k, k - 1), zeros(k, n - (k - 1))];   % exact rank k-1 < k
verifyError(testCase, @() bsqr_rand(M, 'k', k, 'backend', 'mfile'), 'bsqr_rand:RankDeficient');
if bsqr_rand_mex_available()
    verifyError(testCase, @() bsqr_rand(M, 'k', k, 'backend', 'mex'), 'bsqr_rand:RankDeficient');
end
end
