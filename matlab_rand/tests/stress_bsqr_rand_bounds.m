function stats = stress_bsqr_rand_bounds(varargin)
%STRESS_BSQR_RAND_BOUNDS Statistical validation sweep for the bound guarantee.
%   Opt-in stress pass (not auto-discovered by runtests; run it directly).
%   Sweeps seeds x kernel paths x sampling schemes x threshold modes x
%   families and asserts, for every run on orthonormal-row input:
%     * the factorization is exact (Q'A(:,p(1:k)) = R11, Q orthonormal), and
%     * the per-step bound f2_i <= Fhat_i holds at every step (hence the
%       final Osinsky bound ||R11^{-1}||_F <= sqrt(k(n-k+1))).
%   Returns a table of realized quality (frob_inv / osinsky ratio, samples
%   tested per selection) for the paper's empirical claims.
%
%   stress_bsqr_rand_bounds('seeds', 100, 'backend', 'mex') % default
%   stress_bsqr_rand_bounds('seeds', 20, 'backend', 'mfile')

ip = inputParser;
addParameter(ip, 'seeds', 100);
addParameter(ip, 'backend', 'mex');
addParameter(ip, 'k', 16);
addParameter(ip, 'n', 200);
parse(ip, varargin{:});
opt = ip.Results;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab_rand'));
addpath(fullfile(repo_root, 'matlab_rand', 'benchmark'));   % rand_test_matrix

k = opt.k; n = opt.n;
families = {'gaussian', 'spiked_leverage', 'needle'};
rows = cell(0, 8);
nrun = 0;
for fi = 1:numel(families)
    fam = families{fi};
    M = rand_test_matrix(fam, k, n, 4242);        % fixed matrix; RNG varies
    bound_tol = 1e-8 * k * (n - k + 1);
    for batched = [true, false]
        for sampling = {'uniform', 'normweighted'}
            for mode = {'running_mean', 'worstcase_allowance'}
                for s = 1:opt.seeds
                    nrun = nrun + 1;
                    [p, Q, R11, st] = bsqr_rand(M, 'k', k, 'backend', opt.backend, ...
                        'seed', s, 'batched', batched, 'sampling', sampling{1}, ...
                        'threshold_mode', mode{1});
                    % exactness
                    sel = p(1:k);
                    assert(norm(Q' * M(:, sel) - R11, 'fro') < 1e-10 * max(1, norm(M(:, sel), 'fro')), ...
                        'factorization drift: %s/%d/%s/%s seed %d', fam, batched, sampling{1}, mode{1}, s);
                    assert(norm(Q' * Q - eye(k), 'fro') < 1e-10, 'Q orthonormality: seed %d', s);
                    % per-step bound (the guarantee, step by step)
                    assert(all(st.f2 <= st.Fhat + bound_tol), ...
                        'per-step bound violated: %s batched=%d %s/%s seed %d', ...
                        fam, batched, sampling{1}, mode{1}, s);
                    rows(end+1, :) = {fam, batched, sampling{1}, mode{1}, s, ...
                        st.frob_inv / st.osinsky_bound, st.total_tested / k, st.blocks_sampled}; %#ok<AGROW>
                end
            end
        end
    end
    fprintf('%-16s done (%d runs so far, all bounds held)\n', fam, nrun);
end

stats = cell2table(rows, 'VariableNames', ...
    {'family', 'batched', 'sampling', 'mode', 'seed', 'quality_ratio', 'tested_per_k', 'blocks'});
fprintf(['\n%d runs: per-step bound held in every run.\n', ...
    'quality ratio (frob_inv/osinsky): median %.4g, max %.4g\n', ...
    'samples tested per selection:     median %.3g, max %.3g\n'], ...
    nrun, median(stats.quality_ratio), max(stats.quality_ratio), ...
    median(stats.tested_per_k), max(stats.tested_per_k));
end
