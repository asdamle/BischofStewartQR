function tests = test_oracle_parity
%TEST_ORACLE_PARITY V1/V2 validation: production kernels vs. oracle_bsqr.
%   Checks, per docs/VALIDATION_AND_PERF_PLAN.md:
%     * the oracle is self-consistent (orthogonal Q, exact reconstruction);
%     * mfile and mex backends reproduce the oracle's pivot sequence
%       exactly, and Q/R/rinv_r12 to rounding-level tolerances, on a
%       matrix zoo screened for criterion near-ties;
%     * the running criterion sum reproduces ||R11^{-1}||_F^2;
%     * Osinsky bound canaries hold on orthonormal-rows inputs.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab'));
testCase.TestData.zoo = parity_zoo();
end

function testOracleSelfConsistency(testCase)
zoo = testCase.TestData.zoo;
for idx = 1:numel(zoo)
    z = zoo{idx};
    out = oracle_bsqr(z.A, z.k);

    verifyLessThan(testCase, orth_err(out.Q), scaled_tol(size(z.A)), z.name);
    % Q*T reconstructs A(:,p) for any k: T carries both the trapezoid R
    % and the unreduced tail block of an early-stopped factorization.
    verifyLessThan(testCase, rel_resid(z.A(:, out.p), out.Q * out.T), ...
        scaled_tol(size(z.A)), z.name);
    verifyEqual(testCase, sort(out.p), 1:size(z.A, 2), z.name);
end
end

function testZooIsTieFree(testCase)
% Exact pivot equality across implementations is only a fair demand when
% no step has a floating-point near-tie in the criterion. The zoo seeds
% are fixed, so this is a deterministic property; if it ever fails, the
% offending zoo member needs a new seed, not a looser tolerance.
zoo = testCase.TestData.zoo;
for idx = 1:numel(zoo)
    z = zoo{idx};
    out = oracle_bsqr(z.A, z.k);
    verifyGreaterThan(testCase, min(out.crit_gap), 1e-6, ...
        sprintf('%s: near-tie pivot step present', z.name));
end
end

function testMfileMatchesOracle(testCase)
verify_backend_against_oracle(testCase, 'mfile');
end

function testMexMatchesOracle(testCase)
if ~bsqr_mex_available()
    return;
end
verify_backend_against_oracle(testCase, 'mex');
end

function testCriterionTraceMatchesInverseFrobNorm(testCase)
% Writeup Remark "Tracking the full inverse norm": the running sum of the
% selected criterion values equals ||R11^{-1}||_F^2 after every step.
zoo = testCase.TestData.zoo;
for idx = 1:numel(zoo)
    z = zoo{idx};
    if ~z.frob_check
        continue;
    end
    out = oracle_bsqr(z.A, z.k);
    R11 = out.R(1:z.k, 1:z.k);
    frob2 = norm(inv(R11), 'fro')^2;
    err = abs(sum(out.crit_best) - frob2) / frob2;
    verifyLessThan(testCase, err, 1e-10, z.name);
end
end

function testOsinskyBoundCanaries(testCase)
% For M with orthonormal rows, a correct Bischof-Stewart selection always
% satisfies the Osinsky bounds; a subtly wrong criterion eventually won't.
zoo = testCase.TestData.zoo;
backends = available_backends();
for idx = 1:numel(zoo)
    z = zoo{idx};
    if ~z.orthonormal_rows
        continue;
    end
    [kk, n] = size(z.A);
    for b = 1:numel(backends)
        [~, R] = bsqr(z.A, 'backend', backends{b});
        Rinv = inv(R(1:kk, 1:kk));
        verifyLessThanOrEqual(testCase, norm(Rinv, 'fro'), ...
            sqrt(kk * (n - kk + 1)) * (1 + 1e-12), ...
            sprintf('%s/%s: Frobenius Osinsky bound', z.name, backends{b}));
        verifyLessThanOrEqual(testCase, norm(Rinv, 2), ...
            sqrt(1 + kk * (n - kk)) * (1 + 1e-12), ...
            sprintf('%s/%s: spectral Osinsky bound', z.name, backends{b}));
    end
end
end

function testOracleRejectsExactRankDeficiency(testCase)
verifyError(testCase, @() oracle_bsqr(zeros(5, 8)), 'oracle_bsqr:RankDeficient');

A = [randn_seeded(6, 3, 20260601), zeros(6, 5)];
verifyError(testCase, @() oracle_bsqr(A), 'oracle_bsqr:RankDeficient');
end

% ----------------------------------------------------------------------
% Backend-vs-oracle comparison
% ----------------------------------------------------------------------

function verify_backend_against_oracle(testCase, backend)
zoo = testCase.TestData.zoo;
for idx = 1:numel(zoo)
    z = zoo{idx};
    ref = oracle_bsqr(z.A, z.k);

    [Q, R, p, rinv] = bsqr(z.A, 'backend', backend, 'k', z.k, ...
        'pivot_format', 'vector', 'return_rinv_r12', true);

    label = sprintf('%s/%s', z.name, backend);
    verifyEqual(testCase, p(1:z.k), ref.p(1:z.k), ...
        sprintf('%s: pivot sequence', label));

    rtol = 1e-10;
    verifyLessThan(testCase, rel_err(R, ref.R), rtol, ...
        sprintf('%s: R factor', label));
    if z.q_check
        verifyLessThan(testCase, rel_err(Q, ref.Q(:, 1:z.k)), rtol, ...
            sprintf('%s: Q factor', label));
    end

    if z.rinv_check && z.k < size(z.A, 2)
        verifyLessThan(testCase, rel_err(rinv, ref.rinv_r12), 2e-8, ...
            sprintf('%s: rinv_r12', label));
    end
end
end

% ----------------------------------------------------------------------
% Matrix zoo (fixed seeds; near-tie freedom asserted by testZooIsTieFree)
% ----------------------------------------------------------------------

function zoo = parity_zoo()
zoo = {};
zoo{end+1} = zoo_member('gaussian_square_a', gaussian(40, 40, 20260401), 40);
zoo{end+1} = zoo_member('gaussian_square_b', gaussian(48, 48, 20260402), 48);
zoo{end+1} = zoo_member('gaussian_tall', gaussian(60, 24, 20260403), 24);
zoo{end+1} = zoo_member('gaussian_shortwide', gaussian(24, 96, 20260404), 24);
zoo{end+1} = zoo_member('gaussian_early_stop', gaussian(40, 40, 20260405), 16);

% Forward agreement of Q, rinv_r12, and inv(R11) is condition-sensitive
% even with identical pivots, so for graded-spectrum members the valid
% cross-implementation contract is pivot equality + R agreement only.
z = zoo_member('illcond_square', graded(48, 48, 1e-10, 20260406), 48);
z.rinv_check = false;
z.frob_check = false;
z.q_check = false;
zoo{end+1} = z;

z = zoo_member('illcond_shortwide', graded(24, 96, 1e-8, 20260407), 24);
z.rinv_check = false;
z.frob_check = false;
z.q_check = false;
zoo{end+1} = z;

z = zoo_member('orthrows_24x96', orthonormal_rows(24, 96, 20260408), 24);
z.orthonormal_rows = true;
zoo{end+1} = z;

z = zoo_member('orthrows_32x64', orthonormal_rows(32, 64, 20260409), 32);
z.orthonormal_rows = true;
zoo{end+1} = z;

zoo{end+1} = zoo_member('scaled_cols', scaled_columns(32, 48, 20260410), 32);
end

function z = zoo_member(name, A, k)
z = struct('name', name, 'A', A, 'k', k, 'rinv_check', true, ...
    'frob_check', true, 'q_check', true, 'orthonormal_rows', false);
end

function A = gaussian(m, n, seed)
A = randn_seeded(m, n, seed);
end

function A = graded(m, n, sigma_min, seed)
rng(seed, 'twister');
r = min(m, n);
[U, ~] = qr(randn(m, r), 0);
[V, ~] = qr(randn(n, r), 0);
A = U * diag(logspace(0, log10(sigma_min), r)) * V';
end

function A = orthonormal_rows(k, n, seed)
rng(seed, 'twister');
[Qf, ~] = qr(randn(n, k), 0);
A = Qf';
end

function A = scaled_columns(m, n, seed)
rng(seed, 'twister');
A = randn(m, n) .* (10 .^ linspace(-6, 6, n));
end

function A = randn_seeded(m, n, seed)
rng(seed, 'twister');
A = randn(m, n);
end

% ----------------------------------------------------------------------
% Helpers
% ----------------------------------------------------------------------

function backends = available_backends()
backends = {'mfile'};
if bsqr_mex_available()
    backends{end+1} = 'mex';
end
end

function r = rel_err(X, Y)
r = norm(X - Y, 'fro') / max(norm(Y, 'fro'), eps('double'));
end

function r = rel_resid(X, Y)
r = norm(X - Y, 'fro') / max(norm(X, 'fro'), eps('double'));
end

function o = orth_err(Q)
o = norm(eye(size(Q, 2)) - Q' * Q, 'fro');
end

function t = scaled_tol(sz)
t = 8e2 * eps('double') * max(sz);
end
