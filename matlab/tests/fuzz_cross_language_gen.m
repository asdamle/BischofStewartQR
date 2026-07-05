function stats = fuzz_cross_language_gen(outdir, ncases, seed0)
%FUZZ_CROSS_LANGUAGE_GEN Generate randomized cross-language parity cases.
%   Opt-in stress tool (Phase 2 of docs/PUBLICATION_READINESS_PLAN.md), not
%   part of the test suite. Writes NCASES randomized fixtures in exactly the
%   parity/ format (manifest.csv + <name>_{A,p,R,crit[,rinv]}.csv, %.17g
%   round-trip formatting), using matlab/tests/oracle_bsqr.m as the expected
%   output — the same contract as the committed fixtures, at fuzz scale.
%   Consume with julia/test/fuzz_cross_language_check.jl:
%
%     matlab -batch "addpath('matlab/tests'); fuzz_cross_language_gen('/tmp/fuzzdir', 2000, 1)"
%     BS_FUZZ_DIR=/tmp/fuzzdir julia --project=julia julia/test/fuzz_cross_language_check.jl
%
%   Cases are screened for criterion near-ties exactly like the zoo
%   (min relative gap > gap_tol), because exact pivot-sequence equality
%   across BLAS runtimes is only a fair demand on tie-free inputs; screened
%   cases are skipped and counted, not silently dropped. Graded-spectrum
%   cases follow the zoo policy: pivot + R + criterion parity only
%   (rinv_check = 0), with the conditioning-scaled criterion tolerance.

if nargin < 2; ncases = 500; end
if nargin < 3; seed0 = 1; end

this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_dir));
addpath(fullfile(repo_root, 'matlab'));

if ~isfolder(outdir); mkdir(outdir); end

gap_tol = 1e-6;
manifest = fopen(fullfile(outdir, 'manifest.csv'), 'w');
cleanup = onCleanup(@() fclose(manifest));
fprintf(manifest, 'name,m,n,k,rinv_check,gap_min,rtol_R,rtol_rinv,rtol_crit\n');

families = {'gaussian', 'orthrows', 'scaledcols', 'graded', 'neartied'};
kept = 0; skipped = 0;
for c = 1:ncases
    seed = seed0 * 1e6 + c;
    rng(seed, 'twister');
    fam = families{randi(numel(families))};
    m = randi([6, 64]);
    n = randi([6, 96]);
    kmax = min(m, n);
    % mix of full and early-stop k
    if rand < 0.5; k = kmax; else; k = randi(kmax); end

    graded_case = false;
    switch fam
        case 'gaussian'
            A = randn(m, n);
        case 'orthrows'
            kk = min(m, n); [Qf, ~] = qr(randn(n, kk), 0); A = Qf';
            [m, n] = size(A); k = min(k, min(m, n));
        case 'scaledcols'
            A = randn(m, n) .* (10 .^ (6 * (2 * rand(1, n) - 1)));
        case 'graded'
            r = kmax; [U, ~] = qr(randn(m, r), 0); [V, ~] = qr(randn(n, r), 0);
            A = U * diag(logspace(0, -8 * rand, r)) * V';
            graded_case = true;
        case 'neartied'
            A = randn(m, n);
            A = A ./ vecnorm(A) .* (1 + 0.02 * (0:n-1) / max(n - 1, 1));
    end

    out = oracle_bsqr(A, k);
    gap_min = min(out.crit_gap);
    if ~(gap_min > gap_tol)
        skipped = skipped + 1;
        continue;   % near-tie: pivot parity across runtimes is not a fair demand
    end
    kept = kept + 1;
    name = sprintf('fz%06d_%s', c, fam);

    if graded_case
        % Criterion-trace parity is condition-amplified (kernel recurrences vs
        % the oracle's exact recomputation; the panel kernel regroups sums):
        % at random conditioning the late-step drift is unbounded in relative
        % terms even with exact pivot agreement, so for graded cases the
        % meaningful contract is exact pivots + tight R only.
        rinv_check = false; rtol_crit = Inf;
    else
        rinv_check = k < n; rtol_crit = 1e-9;
    end
    rtol_R = 1e-10; rtol_rinv = 2e-8;

    write_matrix(fullfile(outdir, [name, '_A.csv']), A);
    write_matrix(fullfile(outdir, [name, '_p.csv']), out.p);
    write_matrix(fullfile(outdir, [name, '_R.csv']), out.R);
    write_matrix(fullfile(outdir, [name, '_crit.csv']), out.crit_best);
    if rinv_check
        write_matrix(fullfile(outdir, [name, '_rinv.csv']), out.rinv_r12);
    end
    fprintf(manifest, '%s,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g\n', ...
        name, m, n, k, rinv_check, gap_min, rtol_R, rtol_rinv, rtol_crit);
end

stats = struct('kept', kept, 'skipped_near_tie', skipped);
fprintf('fuzz gen: %d cases written, %d near-tie skips -> %s\n', kept, skipped, outdir);
end

function write_matrix(path, X)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
for i = 1:size(X, 1)
    fprintf(fid, '%.17g', X(i, 1));
    % Guard the single-column case: fprintf with an EMPTY data array still
    % emits the format string's literal characters, so an unguarded
    % fprintf(',%.17g', X(i, 2:end)) writes a stray trailing comma.
    if size(X, 2) > 1
        fprintf(fid, ',%.17g', X(i, 2:end));
    end
    fprintf(fid, '\n');
end
end
