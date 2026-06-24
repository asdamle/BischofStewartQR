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
verifyLessThan(testCase, rel_resid(A(:, pk), Qk * Rk), 1.0);
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
verifyLessThan(testCase, rel_resid(A(:, p), Q * R), 1.0);

[~, ~, ~, rinv0] = bsqr(A, 'k', 0, 'return_rinv_r12', true, 'pivot_format', 'vector');
verifySize(testCase, rinv0, [0, size(A, 2)]);

B = randn(10, 10);
[~, ~, ~, rinv_sq] = bsqr(B, 'k', 10, 'return_rinv_r12', true, 'pivot_format', 'vector');
verifySize(testCase, rinv_sq, [10, 0]);
end

function testPivotTieStability(testCase)
C = ones(6, 6);
[~, ~, p] = bsqr(C, 'pivot_format', 'vector');
verifyEqual(testCase, p, 1:6);
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
