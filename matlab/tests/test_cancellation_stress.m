function tests = test_cancellation_stress
%TEST_CANCELLATION_STRESS V2 follow-up to the V0 audit (docs/VALIDATION.md).
%   Stewart's footnote 2 prescribes recomputing an S column ab initio if
%   cancellation occurs in its update; our kernels carry no such guard for
%   W or the wnorm2 recurrence. This test pins the evidence that none is
%   needed: matrices engineered for maximal one-step wnorm2 cancellation
%   (near-duplicate columns whose unit tails are orthogonal to the base
%   range, on top of a graded spectrum driving ||R11^{-1}||_F to ~1e10)
%   leave pivot quality and W-derived outputs at reference accuracy.
%
%   The structural reason: a catastrophic one-step collapse of a large
%   wnorm2 requires the *pivot's* ||w*|| to be large, but the criterion
%   min (1+||w||^2)/rho^2 never selects such a column, so the subtracted
%   term is always moderate relative to machine precision. Measured drift:
%   criterion <= ~1e-7, selection quality ratio == 1, rinv at eps level.
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab'));
end

function testStressQualityAndOutputs(testCase)
for sigma_min = [1e-8, 1e-10]
    A = stress_matrix(sigma_min, 20260502);
    k = min(size(A));

    out = oracle_bsqr(A, k);
    for backend = available_backends()
        [~, R, ~, rinv, tr] = bsqr(A, 'backend', backend{1}, 'k', k, ...
            'pivot_format', 'vector', 'return_rinv_r12', true, 'trace', true);
        label = sprintf('sigma=%g/%s', sigma_min, backend{1});

        % Selected-pivot criterion values stay within the near-tie noise
        % band of the oracle's from-scratch values (measured <= ~1e-7).
        crit_drift = max(abs(tr.crit - out.crit_best) ./ out.crit_best);
        verifyLessThan(testCase, crit_drift, 1e-4, ...
            sprintf('%s: criterion drift', label));

        % Selection quality: pivot flips among engineered near-ties are
        % acceptable, a worse achieved ||R11^{-1}||_F is not.
        f_impl = norm(1 ./ svd(R(1:k, 1:k)));        % ||R11^{-1}||_F via singular values
        Ro = out.R;
        f_oracle = norm(1 ./ svd(Ro(1:k, 1:k)));
        verifyLessThan(testCase, f_impl / f_oracle, 1 + 1e-6, ...
            sprintf('%s: selection quality', label));

        % W-derived output stays consistent with the factor it ships with.
        rinv_true = R(1:k, 1:k) \ R(1:k, k+1:end);
        rinv_err = norm(rinv - rinv_true, 'fro') / max(norm(rinv_true, 'fro'), eps);
        verifyLessThan(testCase, rinv_err, 1e-12, ...
            sprintf('%s: rinv_r12 accuracy', label));
    end
end
end

function A = stress_matrix(sigma_min, seed)
% Base block B with graded spectrum; near-duplicates c*B + T where the
% unit-norm tails T are orthogonal to range(B), so the duplicates' reduced
% top rows are exactly proportional to their base columns (maximal wnorm2
% cancellation when the base pivot is reduced) while their tail norms stay
% O(1), keeping them credible pivot candidates.
rng(seed, 'twister');
m = 32; r = 24;
[Ufull, ~] = qr(randn(m, m));
U = Ufull(:, 1:r);
N = Ufull(:, r+1:end);
[V, ~] = qr(randn(r, r));
B = U * diag(logspace(0, log10(sigma_min), r)) * V';
T = N * randn(m - r, r);
T = T ./ vecnorm(T);
A = [B, (1/3) * B + T];
end

function backends = available_backends()
backends = {'mfile'};
if bsqr_mex_available()
    backends{end+1} = 'mex';
end
end
