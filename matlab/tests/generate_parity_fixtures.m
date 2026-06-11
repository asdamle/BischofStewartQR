function generate_parity_fixtures(outdir)
%GENERATE_PARITY_FIXTURES Write cross-language parity fixtures (plan V3).
%   GENERATE_PARITY_FIXTURES() writes to <repo_root>/parity.
%   GENERATE_PARITY_FIXTURES(OUTDIR) writes to OUTDIR.
%
%   For every parity_zoo member, runs the V1 oracle (oracle_bsqr) and
%   writes the input matrix together with the oracle's expected outputs:
%
%     manifest.csv      name,m,n,k,rinv_check,gap_min,rtol_R,rtol_rinv,rtol_crit
%     <name>_A.csv      input matrix
%     <name>_p.csv      expected full permutation vector (1-indexed)
%     <name>_R.csv      expected k-by-n upper trapezoidal factor
%     <name>_rinv.csv   expected R11 \ R12 (only when rinv_check and k < n)
%     <name>_crit.csv   expected per-step pivot criterion values (V4 trace)
%
%   All values are printed with %.17g so doubles round-trip exactly: the
%   Julia and MATLAB consumers see bit-identical inputs. The acceptance
%   tolerances and the near-tie screen live on the zoo members
%   (parity_zoo.m) and flow to both consumers through the manifest;
%   generation re-asserts the screen so exact pivot-sequence equality is a
%   fair demand across BLAS runtimes.
%
%   Fixtures are committed to git; regenerate only when the zoo changes,
%   and rerun both consumers afterwards:
%     matlab/tests/test_parity_fixtures.m
%     julia/test/runtests.jl

this_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_dir));
if nargin < 1
    outdir = fullfile(repo_root, 'parity');
end
addpath(fullfile(repo_root, 'matlab'));

if ~isfolder(outdir)
    mkdir(outdir);
end

zoo = parity_zoo();

manifest = fopen(fullfile(outdir, 'manifest.csv'), 'w');
cleanup = onCleanup(@() fclose(manifest));
fprintf(manifest, 'name,m,n,k,rinv_check,gap_min,rtol_R,rtol_rinv,rtol_crit\n');

for idx = 1:numel(zoo)
    z = zoo{idx};
    [m, n] = size(z.A);
    out = oracle_bsqr(z.A, z.k);

    gap_min = min(out.crit_gap);
    assert(gap_min > z.gap_tol, 'generate_parity_fixtures:NearTie', ...
        '%s has a near-tie pivot step (gap %.3g); pick a new seed.', ...
        z.name, gap_min);

    write_matrix(fullfile(outdir, [z.name, '_A.csv']), z.A);
    write_matrix(fullfile(outdir, [z.name, '_p.csv']), out.p);
    write_matrix(fullfile(outdir, [z.name, '_R.csv']), out.R);
    write_matrix(fullfile(outdir, [z.name, '_crit.csv']), out.crit_best);

    has_rinv = z.rinv_check && z.k < n;
    if has_rinv
        write_matrix(fullfile(outdir, [z.name, '_rinv.csv']), out.rinv_r12);
    end

    fprintf(manifest, '%s,%d,%d,%d,%d,%.17g,%.17g,%.17g,%.17g\n', ...
        z.name, m, n, z.k, has_rinv, gap_min, z.rtol_R, z.rtol_rinv, z.rtol_crit);
    fprintf('wrote %-22s (%dx%d, k=%d, gap_min=%.3g)\n', z.name, m, n, z.k, gap_min);
end

fprintf('Parity fixtures written to: %s\n', outdir);
end

function write_matrix(path, X)
fid = fopen(path, 'w');
cleanup = onCleanup(@() fclose(fid));
for i = 1:size(X, 1)
    fprintf(fid, '%.17g', X(i, 1));
    fprintf(fid, ',%.17g', X(i, 2:end));
    fprintf(fid, '\n');
end
end
