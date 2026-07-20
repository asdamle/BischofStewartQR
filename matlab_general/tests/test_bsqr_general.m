function tests = test_bsqr_general
%TEST_BSQR_GENERAL Contract tests for the general-A wrapper.
tests = functiontests(localfunctions);
end

% ----- helpers -------------------------------------------------------------

function A = general_matrix(m, n, seed)
rng(seed);
A = randn(m, n);
end

function check_factorization(testCase, A, Q, R, p, k)
[m, n] = size(A);
verifyEqual(testCase, sort(p), 1:n, 'p must be a permutation of 1:n');
verifySize(testCase, Q, [m, k]);
verifySize(testCase, R, [k, n]);
verifyEqual(testCase, R(:, 1:k), triu(R(:, 1:k)), 'AbsTol', 0, 'R11 must be upper triangular');
verifyLessThan(testCase, norm(Q' * Q - eye(k), 'fro'), 1e-10);
scale = max(1, norm(A, 'fro'));
verifyLessThan(testCase, norm(A(:, p(1:k)) - Q * R(:, 1:k), 'fro') / scale, 1e-10);
% Trailing block R12 = Q'*A(:,p(k+1:n)) (bsqr's early-stop semantics; for
% k = m this coincides with the full reconstruction A(:,p) = Q*R).
verifyLessThan(testCase, norm(R(:, k+1:n) - Q' * A(:, p(k+1:n)), 'fro') / scale, 1e-10);
end

% ----- tests ---------------------------------------------------------------

function testReconstructionMatrixPivot(testCase)
shapes = {[8, 40], [16, 64], [12, 12]};
for c = 1:numel(shapes)
    m = shapes{c}(1); n = shapes{c}(2);
    A = general_matrix(m, n, 100 + c);
    [Q, R, E] = bsqr_general(A);
    verifySize(testCase, E, [n, n]);
    verifyEqual(testCase, E' * E, eye(n), 'AbsTol', 0, 'E must be a permutation matrix');
    p = (1:n) * E;   % E = I(:,p), so (1:n)*E recovers the index row
    check_factorization(testCase, A, Q, R, p, m);
    verifyLessThan(testCase, norm(A * E - Q * R, 'fro') / max(1, norm(A, 'fro')), 1e-10);
end
end

function testReconstructionVectorPivot(testCase)
A = general_matrix(10, 50, 7);
[Q, R, p] = bsqr_general(A, 'pivot_format', 'vector');
check_factorization(testCase, A, Q, R, p, 10);
end

function testOutputArities(testCase)
A = general_matrix(9, 30, 11);
[Q, R, p] = bsqr_general(A, 'pivot_format', 'vector');
R1 = bsqr_general(A);
verifyEqual(testCase, R1, R);
[Q2, R2] = bsqr_general(A);
verifyEqual(testCase, Q2, Q);
verifyEqual(testCase, R2, R);
check_factorization(testCase, A, Q, R, p, 9);
end

function testStatsOutput(testCase)
A = general_matrix(12, 80, 3);
[Q, R, p, stats] = bsqr_general(A, 'pivot_format', 'vector');
check_factorization(testCase, A, Q, R, p, 12);
verifyEqual(testCase, stats.selector, 'bsqr');
verifyEqual(testCase, stats.k, 12);
for f = {'t_qr_At', 't_select', 't_assemble', 't_total'}
    verifyGreaterThanOrEqual(testCase, stats.(f{1}), 0);
end
verifyGreaterThan(testCase, stats.rt_diag_ratio, 0);
verifyLessThanOrEqual(testCase, stats.rt_diag_ratio, 1);
verifyEmpty(testCase, stats.selector_stats);

[~, ~, ~, rstats] = bsqr_general(A, 'selector', 'bsqr_rand', 'seed', 5, ...
    'pivot_format', 'vector');
verifyEqual(testCase, rstats.selector, 'bsqr_rand');
verifyTrue(testCase, isstruct(rstats.selector_stats));
verifyTrue(testCase, isfield(rstats.selector_stats, 'frob_inv'));
end

function testRandSelector(testCase)
A = general_matrix(14, 90, 21);
[Q, R, p] = bsqr_general(A, 'selector', 'bsqr_rand', 'seed', 4, 'pivot_format', 'vector');
check_factorization(testCase, A, Q, R, p, 14);
[~, ~, p2] = bsqr_general(A, 'selector', 'bsqr_rand', 'seed', 4, 'pivot_format', 'vector');
verifyEqual(testCase, p, p2, 'same seed must reproduce the same pivots');
end

function testSelectorEquivalence(testCase)
% The wrapper's pivots must equal the selector's pivots on the orthonormal
% row basis G = Qt' from the same qr(A','econ') call -- bitwise-identical
% inputs, so the pivot sequences match exactly.
A = general_matrix(10, 60, 17);
[Qt, ~] = qr(A', 'econ');
G = Qt';
[~, ~, p] = bsqr_general(A, 'pivot_format', 'vector');
[~, ~, pref] = bsqr(G, 'k', 10, 'pivot_format', 'vector', 'check_finite', false);
verifyEqual(testCase, p, pref);
[~, ~, pr] = bsqr_general(A, 'selector', 'bsqr_rand', 'seed', 9, 'pivot_format', 'vector');
prref = bsqr_rand(G, 'k', 10, 'return_r12', true, 'check_finite', false, 'seed', 9);
verifyEqual(testCase, pr, prref);
end

function testEarlyStop(testCase)
A = general_matrix(10, 40, 5);
k = 5;
[Q, R, p] = bsqr_general(A, 'k', k, 'pivot_format', 'vector');
check_factorization(testCase, A, Q, R, p, k);
[Qr, Rr, pr] = bsqr_general(A, 'k', k, 'selector', 'bsqr_rand', 'seed', 2, ...
    'pivot_format', 'vector');
check_factorization(testCase, A, Qr, Rr, pr, k);
end

function testEdgeKZero(testCase)
A = general_matrix(6, 20, 1);
[Q, R, p, stats] = bsqr_general(A, 'k', 0, 'pivot_format', 'vector');
verifySize(testCase, Q, [6, 0]);
verifySize(testCase, R, [0, 20]);
verifyEqual(testCase, sort(p), 1:20);
verifyEqual(testCase, stats.k, 0);
end

function testBadlyScaledColumns(testCase)
% Column norms over 12 decades: the reduction must still reconstruct.
A = general_matrix(10, 60, 13) .* logspace(0, -12, 60);
[Q, R, p] = bsqr_general(A, 'pivot_format', 'vector');
check_factorization(testCase, A, Q, R, p, 10);
end

function testSquareInput(testCase)
A = general_matrix(12, 12, 19);
[Q, R, p] = bsqr_general(A, 'pivot_format', 'vector');
check_factorization(testCase, A, Q, R, p, 12);
end

function testBackendForwarding(testCase)
% 'backend' is not a wrapper option; it must forward to the selector. The
% deterministic pivot sequence is backend-invariant (parity-tested in
% matlab/tests), so mfile and auto must agree exactly here too.
A = general_matrix(10, 45, 23);
[Q, R, p] = bsqr_general(A, 'backend', 'mfile', 'pivot_format', 'vector');
check_factorization(testCase, A, Q, R, p, 10);
[~, ~, p2] = bsqr_general(A, 'pivot_format', 'vector');
verifyEqual(testCase, p, p2);
end

function testInputTypeContract(testCase)
A = general_matrix(8, 30, 29);
verifyError(testCase, @() bsqr_general(sparse(A)), 'bsqr_general:InvalidInput');
[Q, R, p] = bsqr_general(single(A), 'pivot_format', 'vector');
check_factorization(testCase, double(single(A)), Q, R, p, 8);
end

function testIllConditionedWarning(testCase)
rng(37);
v = randn(1, 30);
A = [v; 2 * v];          % exact rank 1
verifyWarning(testCase, @() bsqr_general(A), 'bsqr_general:IllConditioned');
% A mildly ill-conditioned A only warns when rank_tol is raised to catch it.
w = randn(1, 30);
B = [v; v + 1e-6 * w];
verifyWarningFree(testCase, @() bsqr_general(B));
verifyWarning(testCase, @() bsqr_general(B, 'rank_tol', 1e-3), ...
    'bsqr_general:IllConditioned');
end

function testArgumentValidation(testCase)
A = general_matrix(8, 30, 31);
verifyError(testCase, @() bsqr_general(A'), 'bsqr_general:NotShortWide');
verifyError(testCase, @() bsqr_general(A + 1i), 'bsqr_general:InvalidInput');
verifyError(testCase, @() bsqr_general(A, 'selector', 'bogus'), 'bsqr_general:InvalidSelector');
verifyError(testCase, @() bsqr_general(A, 'k', 1.5), 'bsqr_general:InvalidK');
verifyError(testCase, @() bsqr_general(A, 'k', 99), 'bsqr_general:InvalidK');
verifyError(testCase, @() bsqr_general(A, 'pivot_format', 'bogus'), 'bsqr_general:InvalidPivotFormat');
verifyError(testCase, @() bsqr_general(A, 'return_r12', true), 'bsqr_general:UnsupportedOption');
verifyError(testCase, @() bsqr_general(A, 'return_rinv_r12', true), 'bsqr_general:UnsupportedOption');
verifyError(testCase, @() bsqr_general(A, 'trace', true), 'bsqr_general:UnsupportedOption');
verifyError(testCase, @() five_outputs(A), 'bsqr_general:TooManyOutputs');
B = A; B(3, 7) = NaN;
verifyError(testCase, @() bsqr_general(B), 'bsqr_general:NonFiniteInput');
end

function testForwardedOptionRouting(testCase)
% Unknown forwards are rejected by the selector's own parser, and
% selector-specific options must reach the right selector: 'sampling' is a
% bsqr_rand option, so the bsqr selector rejects it while bsqr_rand accepts.
A = general_matrix(8, 30, 33);
verifyTrue(testCase, errors_out(@() bsqr_general(A, 'bogus_option', 1)));
verifyTrue(testCase, errors_out(@() bsqr_general(A, 'selector', 'bsqr_rand', 'bogus_option', 1)));
verifyTrue(testCase, errors_out(@() bsqr_general(A, 'sampling', 'uniform')));
[Q, R, p] = bsqr_general(A, 'selector', 'bsqr_rand', 'sampling', 'uniform', ...
    'seed', 1, 'pivot_format', 'vector');
check_factorization(testCase, A, Q, R, p, 8);
end

function out = five_outputs(A)
[~, ~, ~, ~, out] = bsqr_general(A); %#ok<ASGLU>
end

function tf = errors_out(fn)
tf = false;
try
    fn();
catch
    tf = true;
end
end
