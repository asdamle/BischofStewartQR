function tests = test_bsqr
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab'));
addpath(fullfile(repo_root, 'matlab', 'benchmark'));
testCase.TestData.repo_root = repo_root;
end

function testOutputContracts(testCase)
rng(20260314, 'twister');
A = randn(14, 10);

R = bsqr(A);
verifySize(testCase, R, [10, 10]);

[Q2, R2] = bsqr(A);
verifySize(testCase, Q2, [14, 10]);
verifySize(testCase, R2, [10, 10]);

[Qv, Rv, p] = bsqr(A, 'pivot_format', 'vector');
verifySize(testCase, p, [1, 10]);
verifyEqual(testCase, sort(p), 1:10);
verifyLessThan(testCase, rel_resid(A(:, p), Qv * Rv), scaled_tol(size(A)));
verifyLessThan(testCase, orth_err(Qv), scaled_tol(size(A)));

[Qm, Rm, E] = bsqr(A, 'pivot_format', 'matrix');
verifySize(testCase, E, [10, 10]);
verifyLessThan(testCase, rel_resid(A * E, Qm * Rm), scaled_tol(size(A)));
end

function testShapeCoverage(testCase)
rng(7, 'twister');
shapes = [12 8; 10 10; 8 12; 3 20; 20 3];
for idx = 1:size(shapes, 1)
    m = shapes(idx, 1);
    n = shapes(idx, 2);
    A = randn(m, n);
    [Q, R, p] = bsqr(A, 'pivot_format', 'vector');
    k = min(m, n);
    verifySize(testCase, Q, [m, k]);
    verifySize(testCase, R, [k, n]);
    verifyEqual(testCase, sort(p), 1:n);
    verifyLessThan(testCase, rel_resid(A(:, p), Q * R), scaled_tol([m, n]));
end
end

function testEarlyStoppingByK(testCase)
rng(17, 'twister');
A = randn(15, 11);

[Q0, R0, p0] = bsqr(A, 'k', 0, 'pivot_format', 'vector');
verifySize(testCase, Q0, [15, 0]);
verifySize(testCase, R0, [0, 11]);
verifyEqual(testCase, p0, 1:11);

k = 4;
[Qk, Rk, pk] = bsqr(A, 'k', k, 'pivot_format', 'vector');
verifySize(testCase, Qk, [15, k]);
verifySize(testCase, Rk, [k, 11]);
verifyEqual(testCase, sort(pk), 1:11);
% Early stop: A(:,p) = Q*R holds only for the selected block; the trailing
% block of R is the unselected columns' projection onto span(Q). These are
% the tight identities (documented in bsqr.m) -- not rel_resid(A(:,p), Q*R),
% which is O(1) for k < min(m,n).
verifyLessThan(testCase, rel_resid(A(:, pk(1:k)), Qk * Rk(:, 1:k)), scaled_tol(size(A)));
verifyLessThan(testCase, ...
    norm(Qk' * A(:, pk(k+1:end)) - Rk(:, k+1:end), 'fro') / max(norm(A, 'fro'), eps), ...
    scaled_tol(size(A)));
verifyLessThan(testCase, orth_err(Qk), scaled_tol(size(A)));
end

function testOptionalRinvR12(testCase)
rng(2026, 'twister');
A = randn(18, 12);
k = 8;

[Q, R, p, rinv] = bsqr(A, 'k', k, 'pivot_format', 'vector', 'return_rinv_r12', true);
verifySize(testCase, Q, [18, k]);
verifySize(testCase, R, [k, 12]);
verifySize(testCase, rinv, [k, 12 - k]);

expected = R(1:k, 1:k) \ R(1:k, k+1:end);
err = norm(rinv - expected, 'fro') / max(norm(expected, 'fro'), eps('double'));
verifyLessThan(testCase, err, 2e-8);
% k < min(m,n): the tight identities are selected-block reconstruction and
% the R12 projection (A(:,p) = Q*R does not hold under early stop).
verifyLessThan(testCase, rel_resid(A(:, p(1:k)), Q * R(:, 1:k)), scaled_tol(size(A)));
verifyLessThan(testCase, ...
    norm(Q' * A(:, p(k+1:end)) - R(:, k+1:end), 'fro') / max(norm(A, 'fro'), eps('double')), ...
    scaled_tol(size(A)));

[~, ~, ~, rinv0] = bsqr(A, 'k', 0, 'return_rinv_r12', true, 'pivot_format', 'vector');
verifySize(testCase, rinv0, [0, size(A, 2)]);
[~, ~, ~, rinv0m] = bsqr(A, 'k', 0, 'return_rinv_r12', true, ...
    'pivot_format', 'vector', 'backend', 'mfile');   % m-file k = 0 branch too
verifySize(testCase, rinv0m, [0, size(A, 2)]);

B = randn(10, 10);
[~, ~, ~, rinv_sq] = bsqr(B, 'k', 10, 'return_rinv_r12', true, 'pivot_format', 'vector');
verifySize(testCase, rinv_sq, [10, 0]);
end

function testPivotTieStability(testCase)
C = ones(6, 6);
[~, ~, p] = bsqr(C, 'pivot_format', 'vector');
verifyEqual(testCase, p, 1:6);
end

function testTieBreakFirstMinimum(testCase)
% Bitwise-identical duplicate columns give bitwise-equal criteria within one
% backend, so the strict `<` tie-break must pick the earlier candidate --
% including mid-factorization. Checked per backend; cross-backend tie
% agreement is deliberately not a contract (the parity zoo screens ties out
% because tie outcomes are not BLAS-portable).
rng(20260709, 'twister');
c = randn(10, 1); c = c * 5 / norm(c);
b = randn(10, 1); b = b * 9 / norm(b);
A1 = [c, c, randn(10, 6)];          % tie at step 1 between cols 1,2
A2 = [b, c, c, randn(10, 5)];       % tie at step 2 between cols 2,3
backends = {'mfile'};
if bsqr_mex_available(); backends{end+1} = 'mex'; end
for bi = 1:numel(backends)
    [~, ~, p1] = bsqr(A1, 'k', 3, 'backend', backends{bi}, 'pivot_format', 'vector');
    verifyEqual(testCase, p1(1), 1);
    [~, ~, p2] = bsqr(A2, 'k', 2, 'backend', backends{bi}, 'pivot_format', 'vector');
    verifyEqual(testCase, p2(1), 1);
    verifyEqual(testCase, p2(2), 2);
end
end

function testInputTypeContract(testCase)
% Both backends must accept the same inputs. Sparse is rejected outright (a
% sparse mxArray passes the MEX's double/real/2-D checks but mxGetPr yields
% only the nonzero storage -- reading it as dense is a buffer overread);
% non-double numerics are normalized to double by the dispatcher so behavior
% cannot depend on which backend happens to be available.
rng(20260710, 'twister');
A = randn(8, 12);
[~, ~, pref_single] = bsqr(double(single(A)), 'backend', 'mfile', 'pivot_format', 'vector');
Ai = int32(round(A * 10));
[~, ~, pref_int] = bsqr(double(Ai), 'backend', 'mfile', 'pivot_format', 'vector');
backends = {'mfile'};
if bsqr_mex_available(); backends{end+1} = 'mex'; end
for bi = 1:numel(backends)
    b = backends{bi};
    verifyError(testCase, @() bsqr(sparse(A), 'backend', b), 'bsqr:InvalidInput');
    [~, ~, ps] = bsqr(single(A), 'backend', b, 'pivot_format', 'vector');
    verifyEqual(testCase, ps, pref_single);
    [~, ~, pi_] = bsqr(Ai, 'backend', b, 'pivot_format', 'vector');
    verifyEqual(testCase, pi_, pref_int);
end
if bsqr_mex_available()
    verifyError(testCase, @() bsqr_mex(sparse(A)), 'bsqr:InvalidInput');
end
end

function testExtremeColumnScaling(testCase)
% The pivot criterion tracks SQUARED column norms, which over/underflow for
% finite inputs with ||a_j|| beyond ~1e154 / ~1e-154. Outside that domain the
% selection guarantee is void (documented in bsqr.m), but the factorization
% itself must remain exact -- Householder generation is scale-safe. Codifies
% that behavior on both backends.
rng(20260712, 'twister');
A = randn(8, 12);
A(:, 3) = A(:, 3) * 1e160;    % squared norm overflows to Inf
A(:, 7) = A(:, 7) * 1e-170;   % squared norm underflows to 0
backends = {'mfile'};
if bsqr_mex_available(); backends{end+1} = 'mex'; end
for bi = 1:numel(backends)
    [Q, R, p] = bsqr(A, 'backend', backends{bi}, 'pivot_format', 'vector');
    verifyEqual(testCase, sort(p), 1:12);
    verifyLessThan(testCase, rel_resid(A(:, p), Q * R), scaled_tol(size(A)));
    verifyLessThan(testCase, orth_err(Q), scaled_tol(size(A)));
end
end

function testMexWorkspaceReuseSequence(testCase)
% bsqr_mex keeps a persistent workspace across calls; results must be
% independent of call history. Interleave grow/shrink/reshape/early-stop
% calls in one process and check every result against a fresh m-file
% reference (exact pivot parity is the repo contract on tie-free input).
if ~bsqr_mex_available()
    return;
end
rng(20260713, 'twister');
cases = {randn(30, 20), randn(6, 6), randn(40, 10), randn(10, 40), randn(25, 25)};
ks = [20, 6, 10, 10, 12];
pref = cell(1, numel(cases));
for i = 1:numel(cases)
    [~, ~, pref{i}] = bsqr(cases{i}, 'k', ks(i), 'backend', 'mfile', 'pivot_format', 'vector');
end
order = [1 2 3 4 5 5 4 3 2 1 3 1 5];
for idx = order
    [Q, R, p] = bsqr(cases{idx}, 'k', ks(idx), 'backend', 'mex', 'pivot_format', 'vector');
    verifyEqual(testCase, p, pref{idx});
    k = ks(idx);
    verifyLessThan(testCase, rel_resid(cases{idx}(:, p(1:k)), Q * R(:, 1:k)), ...
        scaled_tol(size(cases{idx})));
end
end

function testExactRankDeficiencyFullK(testCase)
% k past the exact rank without rank stopping (the MATLAB API has no
% rank_stop): mid-kernel the pivot column is exactly zero, so beta = 0 and
% the invdiag = 0 / tau = 0 guards must hold (no division by zero, reflector
% apply skipped). Reconstruction must remain exact -- the trailing R rows
% are simply zero. Both backends.
rng(20260714, 'twister');
r = 4;
A = randn(12, r) * randn(r, 8);   % exact rank 4, k = min(m,n) = 8 > 4
backends = {'mfile'};
if bsqr_mex_available(); backends{end+1} = 'mex'; end
for bi = 1:numel(backends)
    [Q, R, p] = bsqr(A, 'backend', backends{bi}, 'pivot_format', 'vector');
    verifyTrue(testCase, all(isfinite(Q(:))) && all(isfinite(R(:))));
    verifyLessThan(testCase, rel_resid(A(:, p), Q * R), scaled_tol(size(A)));
    verifyLessThan(testCase, orth_err(Q), scaled_tol(size(A)));
end
end

function testArgumentValidation(testCase)
% Each parser rejection branch, by identifier.
A = randn(6, 8);
verifyError(testCase, @() bsqr(A + 1i), 'bsqr:InvalidInput');
verifyError(testCase, @() bsqr(randn(3, 3, 3)), 'bsqr:InvalidInput');
verifyError(testCase, @() bsqr(A, 'k', 1.5), 'bsqr:InvalidK');
verifyError(testCase, @() bsqr(A, 'k', 99), 'bsqr:InvalidK');
verifyError(testCase, @() bsqr(A, 'pivot_format', 'bogus'), 'bsqr:InvalidPivotFormat');
verifyError(testCase, @() bsqr(A, 'backend', 'bogus'), 'bsqr:InvalidBackend');
verifyError(testCase, @() bsqr(A, 'norm_recomp_tol', 2), 'bsqr:InvalidNormRecompTol');
verifyError(testCase, @() bsqr([1 Inf; 2 3]), 'bsqr:NonFiniteInput');
end

function testEdgeShapes(testCase)
% Degenerate and tall shapes through the public API on both backends.
rng(20260711, 'twister');
shapes = [0 5; 5 0; 1 8; 8 1; 50 5];
backends = {'mfile'};
if bsqr_mex_available(); backends{end+1} = 'mex'; end
for idx = 1:size(shapes, 1)
    m = shapes(idx, 1); n = shapes(idx, 2); k = min(m, n);
    A = randn(m, n);
    for bi = 1:numel(backends)
        [Q, R, p] = bsqr(A, 'backend', backends{bi}, 'pivot_format', 'vector');
        verifySize(testCase, Q, [m, k]);
        verifySize(testCase, R, [k, n]);
        verifyEqual(testCase, sort(p), 1:n);
        if k > 0
            verifyLessThan(testCase, rel_resid(A(:, p), Q * R), scaled_tol([m, n]));
        end
    end
end
end

function testBackendOptions(testCase)
A = randn(8, 5);
[Q, R, p] = bsqr(A, 'backend', 'auto', 'pivot_format', 'vector');
verifySize(testCase, Q, [8, 5]);
verifySize(testCase, R, [5, 5]);
verifyEqual(testCase, sort(p), 1:5);

if ~bsqr_mex_available()
    verifyError(testCase, @() bsqr(A, 'backend', 'mex'), 'bsqr:MexUnavailable');
end
end

function testMexParityWhenAvailable(testCase)
if ~bsqr_mex_available()
    return;
end

rng(123, 'twister');
A = randn(20, 14);
k = 9;
[Qm, Rm, pm, Wm] = bsqr(A, 'backend', 'mfile', 'k', k, 'pivot_format', 'vector', 'return_rinv_r12', true);
[Qx, Rx, px, Wx] = bsqr(A, 'backend', 'mex', 'k', k, 'pivot_format', 'vector', 'return_rinv_r12', true);

verifyEqual(testCase, px, pm);
verifyLessThan(testCase, norm(Rx - Rm, 'fro') / max(norm(Rm, 'fro'), eps('double')), 1e-11);
verifyLessThan(testCase, norm(Qx - Qm, 'fro') / max(norm(Qm, 'fro'), eps('double')), 1e-11);
verifyLessThan(testCase, norm(Wx - Wm, 'fro') / max(norm(Wm, 'fro'), eps('double')), 1e-11);
end

function testBenchmarkSmoke(testCase)
if ~bsqr_mex_available()
    return;
end

old_formats = getenv('BS_MATLAB_PUB_FIG_FORMATS');
clean_formats = onCleanup(@() setenv('BS_MATLAB_PUB_FIG_FORMATS', old_formats));
setenv('BS_MATLAB_PUB_FIG_FORMATS', 'png,pdf,eps');

outdir = fullfile(tempdir, ['bsqr_matlab_smoke_', char(string(randi(1e9)))]);
cfg = struct();
cfg.outdir = outdir;
cfg.seeds = 20260310;
cfg.families = {'gaussian'};
cfg.square_ms = 32;
cfg.short_ms = 16;
cfg.short_aspects = 2;
cfg.warmup = 1;
cfg.samples = 3;
cfg.norm_recomp_tol = sqrt(eps('double'));

run_publication_smoke_benchmark(cfg);

verifyTrue(testCase, isfile(fullfile(outdir, 'publication_timings.csv')));
verifyTrue(testCase, isfile(fullfile(outdir, 'publication_summary.md')));
verifyTrue(testCase, isfile(fullfile(outdir, 'metadata.txt')));
verifyTrue(testCase, isfile(fullfile(outdir, 'plots', 'plain', 'fig_square_runtime.png')));
verifyTrue(testCase, isfile(fullfile(outdir, 'plots', 'plain', 'figure_captions.md')));
verifyTrue(testCase, isfile(fullfile(outdir, 'tables', 'plain', 'table_square_relative_time.csv')));
verifyTrue(testCase, isfile(fullfile(outdir, 'plots', 'rinv', 'fig_shortwide_runtime.pdf')));
verifyTrue(testCase, isfile(fullfile(outdir, 'plots', 'fig_relative_time_composite.pdf')));
verifyTrue(testCase, isfile(fullfile(outdir, 'plots', 'fig_relative_time_composite.eps')));

tbl = readtable(fullfile(outdir, 'tables', 'plain', 'table_square_relative_time.csv'));
verifyTrue(testCase, ismember('relative_time_geomean', tbl.Properties.VariableNames));
end

function testBenchmarkRejectsBackendField(testCase)
cfg = struct();
cfg.outdir = fullfile(tempdir, ['bsqr_matlab_reject_backend_', char(string(randi(1e9)))]);
cfg.seeds = 20260310;
cfg.families = {'gaussian'};
cfg.square_ms = 8;
cfg.short_ms = 4;
cfg.short_aspects = 2;
cfg.warmup = 0;
cfg.samples = 1;
cfg.norm_recomp_tol = sqrt(eps('double'));
cfg.bsqr_backend = 'mex';

verifyError(testCase, @() run_publication_benchmarks(cfg), 'run_publication_benchmarks:UnsupportedConfigField');
end

function testBenchmarkOutdirSeparationGuard(testCase)
repo_root = testCase.TestData.repo_root;
julia_out = fullfile(repo_root, 'benchmark', 'results', 'publication');
julia_out_new = fullfile(repo_root, 'julia', 'benchmark', 'results', 'publication');

cfg = struct();
cfg.outdir = julia_out;
cfg.seeds = 20260310;
cfg.families = {'gaussian'};
cfg.square_ms = 8;
cfg.short_ms = 4;
cfg.short_aspects = 2;
cfg.warmup = 0;
cfg.samples = 1;
cfg.norm_recomp_tol = sqrt(eps('double'));

verifyError(testCase, @() run_publication_benchmarks(cfg), 'run_publication_benchmarks:SharedOutdirBlocked');
verifyError(testCase, @() plot_publication_results(fullfile(tempdir, 'dummy.csv'), julia_out, julia_out), ...
    'plot_publication_results:SharedOutdirBlocked');

cfg.outdir = julia_out_new;
verifyError(testCase, @() run_publication_benchmarks(cfg), 'run_publication_benchmarks:SharedOutdirBlocked');
verifyError(testCase, @() plot_publication_results(fullfile(tempdir, 'dummy.csv'), julia_out_new, julia_out_new), ...
    'plot_publication_results:SharedOutdirBlocked');
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

function t = scaled_tol(sz)
t = 8e2 * eps('double') * max(sz);
end
