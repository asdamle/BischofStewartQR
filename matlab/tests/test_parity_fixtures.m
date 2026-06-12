function tests = test_parity_fixtures
%TEST_PARITY_FIXTURES V3 validation: committed fixtures vs. MATLAB backends.
%   The fixtures under <repo_root>/parity are oracle outputs written by
%   generate_parity_fixtures. This test asserts that
%     * the committed inputs are bit-identical to the current parity_zoo
%       (catches "zoo changed but fixtures not regenerated" drift);
%     * both MATLAB backends reproduce the expected pivot sequence exactly
%       and the factors to tolerance.
%   The Julia consumer of the same files is julia/test/test_parity_fixtures.jl.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab'));
fixdir = fullfile(repo_root, 'parity');
assertTrue(testCase, isfolder(fixdir), ...
    'parity/ fixtures missing; run matlab/tests/generate_parity_fixtures.m');
testCase.TestData.fixdir = fixdir;
testCase.TestData.manifest = readtable(fullfile(fixdir, 'manifest.csv'), ...
    'TextType', 'char');
end

function testFixturesMatchCurrentZoo(testCase)
zoo = parity_zoo();
manifest = testCase.TestData.manifest;
verifyEqual(testCase, height(manifest), numel(zoo), 'zoo size vs manifest');
for idx = 1:numel(zoo)
    z = zoo{idx};
    A_file = read_fixture(testCase, [z.name, '_A.csv']);
    verifyTrue(testCase, isequal(A_file, z.A), ...
        sprintf('%s: committed input differs from parity_zoo; regenerate fixtures', z.name));
end
end

function testMfileMatchesFixtures(testCase)
verify_backend_against_fixtures(testCase, 'mfile');
end

function testMexMatchesFixtures(testCase)
if ~bsqr_mex_available()
    return;
end
verify_backend_against_fixtures(testCase, 'mex');
end

function verify_backend_against_fixtures(testCase, backend)
manifest = testCase.TestData.manifest;
for row = 1:height(manifest)
    name = manifest.name{row};
    k = manifest.k(row);
    n = manifest.n(row);

    A = read_fixture(testCase, [name, '_A.csv']);
    p_exp = read_fixture(testCase, [name, '_p.csv']);
    R_exp = read_fixture(testCase, [name, '_R.csv']);
    crit_exp = read_fixture(testCase, [name, '_crit.csv']);

    [~, R, p, rinv, trace] = bsqr(A, 'backend', backend, 'k', k, ...
        'pivot_format', 'vector', 'return_rinv_r12', true, 'trace', true);

    label = sprintf('%s/%s', name, backend);
    verifyEqual(testCase, p, p_exp, sprintf('%s: pivot sequence', label));
    verifyLessThan(testCase, rel_err(R, R_exp), manifest.rtol_R(row), ...
        sprintf('%s: R factor', label));

    crit_err = max(abs(trace.crit - crit_exp) ./ crit_exp);
    verifyLessThan(testCase, crit_err, manifest.rtol_crit(row), ...
        sprintf('%s: criterion trace', label));

    if manifest.rinv_check(row) && k < n
        rinv_exp = read_fixture(testCase, [name, '_rinv.csv']);
        verifyLessThan(testCase, rel_err(rinv, rinv_exp), manifest.rtol_rinv(row), ...
            sprintf('%s: rinv_r12', label));
    end
end
end

function testTraceParityBetweenBackends(testCase)
% Recompute events are threshold crossings on running norms, so equality
% between the two MATLAB backends pins both the guard logic and the norm
% downdate arithmetic step by step (measured exactly equal on the zoo).
% Pinned to the unblocked kernels: the panel kernel batches guard trips
% per panel (one flag per column), so its counts differ by design — and
% this keeps the unblocked MEX path covered now that panel is the default.
if ~bsqr_mex_available()
    return;
end
old_nb = getenv('BS_PANEL_NB');
restore_nb = onCleanup(@() setenv('BS_PANEL_NB', old_nb));
setenv('BS_PANEL_NB', '0');
manifest = testCase.TestData.manifest;
for row = 1:height(manifest)
    name = manifest.name{row};
    A = read_fixture(testCase, [name, '_A.csv']);
    k = manifest.k(row);
    [~, ~, ~, ~, tm] = bsqr(A, 'backend', 'mfile', 'k', k, ...
        'pivot_format', 'vector', 'trace', true);
    [~, ~, ~, ~, tx] = bsqr(A, 'backend', 'mex', 'k', k, ...
        'pivot_format', 'vector', 'trace', true);
    verifyEqual(testCase, tx.nrecomp, tm.nrecomp, ...
        sprintf('%s: nrecomp trace mfile vs mex', name));
end
end

function X = read_fixture(testCase, fname)
X = readmatrix(fullfile(testCase.TestData.fixdir, fname));
end

function r = rel_err(X, Y)
r = norm(X - Y, 'fro') / max(norm(Y, 'fro'), eps('double'));
end
