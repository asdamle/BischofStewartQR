function zoo = parity_zoo()
%PARITY_ZOO Fixed-seed matrix zoo shared by oracle parity tests and the
%   cross-language fixture generator (docs/VALIDATION_AND_PERF_PLAN.md V2/V3).
%
%   Every member must be free of criterion near-ties (asserted by
%   testZooIsTieFree and re-asserted at fixture generation time): exact
%   pivot-sequence equality across implementations and BLAS runtimes is
%   only a fair demand under that screening. A failing member gets a new
%   seed, never a looser tolerance.
%
%   Member flags:
%     rinv_check        compare rinv_r12 against the oracle (off for
%                       graded spectra, where the solve is condition-
%                       amplified even with identical pivots)
%     frob_check        check sum(crit_best) == ||R11^{-1}||_F^2
%     q_check           compare Q columns against the oracle (off for
%                       graded spectra; forward Q error scales with cond)
%     orthonormal_rows  member is a valid Osinsky-bound canary
%
%   Member contract values (single home for every consumer, including the
%   Julia fixture tests, which receive them through manifest.csv):
%     gap_tol           minimum relative criterion gap for the near-tie
%                       screen
%     rtol_R, rtol_rinv relative tolerances for factor agreement
%     rtol_crit         relative tolerance for the per-step criterion
%                       trace vs. the oracle (kernels use the wnorm2
%                       recurrence and downdated norms, the oracle exact
%                       recomputation; drift scales with conditioning —
%                       measured ~1e-12 for well-conditioned members,
%                       ~2e-7 for the graded ones)

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
z.rtol_crit = 1e-4;
zoo{end+1} = z;

z = zoo_member('illcond_shortwide', graded(24, 96, 1e-8, 20260407), 24);
z.rinv_check = false;
z.frob_check = false;
z.q_check = false;
z.rtol_crit = 1e-4;
zoo{end+1} = z;

z = zoo_member('orthrows_24x96', orthonormal_rows(24, 96, 20260408), 24);
z.orthonormal_rows = true;
zoo{end+1} = z;

z = zoo_member('orthrows_32x64', orthonormal_rows(32, 64, 20260409), 32);
z.orthonormal_rows = true;
zoo{end+1} = z;

zoo{end+1} = zoo_member('scaled_cols', scaled_columns(32, 48, 20260410), 32);

% --- Additional stress corners (2026-06-30 parity expansion) ---------------
% These exercise pivot-criterion corners not covered above. Every member is
% still screened tie-free (min crit_gap > gap_tol); a failing screen means a
% new seed, never a looser tolerance.

% Near-tied column norms: all columns share a norm within a 2% band, broken
% only by a small deterministic geometric factor. Exactly-equal norms are
% deliberately NOT used here -- the step-1 criterion is 1/||a_j||^2, so equal
% norms force an n-way tie at step 1 (incompatible with the tie-free screen);
% that exact-tie / strict-'<' first-min case is pinned separately by
% testPivotTieStability (all-ones input) and the Julia "Criterion-consistent
% pivot sequence" testset. Here the near-tied norms stress the criterion's
% denominator while keeping every step strictly resolved.
zoo{end+1} = zoo_member('near_tied_col_norms', near_tied_norms(28, 22, 20260702), 22);

% Exact rank deficiency: 12 independent columns and 12 exact-zero columns,
% column order shuffled so the zeros are interleaved. k = 12 = rank, so the
% oracle completes 12 full-rank steps (its rank-deficiency guard never trips)
% while the kernels must still skip the zero-tail columns (criterion = Inf)
% at every step. k < n also exercises the rinv_r12 path.
zoo{end+1} = zoo_member('rank_deficient_tall', rank_deficient(40, 24, 12, 20260702), 12);

% Tightly clustered singular values: spectrum within +/-1% of 1, so the
% matrix is well conditioned but every pivot choice is a close call in the
% criterion. Screened tie-free; keeps full forward-agreement checks on.
zoo{end+1} = zoo_member('clustered_svals', clustered_svals(44, 44, 20260704), 44);

% Strongly graded / decaying spectrum, tall (m > n): sigma from 1 down to
% 1e-9. Like the illcond_* members the forward solve is condition-amplified,
% so only pivot sequence + R agreement is a fair cross-implementation
% contract; the running-norm downdate also decays enough here to fire the
% sqrt(eps) recompute path.
z = zoo_member('graded_tall', graded(64, 28, 1e-9, 20260704), 28);
z.rinv_check = false;
z.frob_check = false;
z.q_check = false;
z.rtol_crit = 1e-4;
zoo{end+1} = z;

% Norm-recompute trigger: two thirds of the columns are placed almost inside
% the span of a small fixed basis (tiny orthogonal residual), so once that
% basis is pivoted in, their running tail-norm estimate decays past
% s_ref * sqrt(eps) and forces an exact norm refresh in the downdate. The
% collapse makes the factor ill conditioned, so treat it like the graded
% members (pivot + R agreement only).
z = zoo_member('norm_recompute', norm_recompute_case(36, 30, 20260711), 20);
z.rinv_check = false;
z.frob_check = false;
z.q_check = false;
z.rtol_crit = 1e-4;
zoo{end+1} = z;

% Near-rank-deficient, wide (m < n): full-rank leading block plus four
% moderate columns and a batch of 1e-9-scale near-null columns, shuffled.
% Wide shape with k < n exercises rinv_r12; the near-null columns' tail norms
% decay enough to fire the recompute path, yet R11 stays well conditioned so
% every forward-agreement check can stay on.
zoo{end+1} = zoo_member('near_deficient_wide', near_deficient_wide(18, 40, 20260706), 18);
end

function z = zoo_member(name, A, k)
z = struct('name', name, 'A', A, 'k', k, 'rinv_check', true, ...
    'frob_check', true, 'q_check', true, 'orthonormal_rows', false, ...
    'gap_tol', 1e-6, 'rtol_R', 1e-10, 'rtol_rinv', 2e-8, 'rtol_crit', 1e-9);
end

function A = gaussian(m, n, seed)
rng(seed, 'twister');
A = randn(m, n);
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

function A = near_tied_norms(m, n, seed)
% Unit-norm random columns rescaled into a narrow 2% norm band. Norms are
% near-tied but strictly separated (the scale factors are distinct and the
% order is shuffled), so the pivot criterion is well resolved at every step.
rng(seed, 'twister');
A = randn(m, n);
A = A ./ vecnorm(A);
scal = 1 + 2e-2 * linspace(0, 1, n);
scal = scal(randperm(n));
A = A .* scal;
end

function A = rank_deficient(m, n, r, seed)
% r independent Gaussian columns padded with n-r exact-zero columns, then
% column-shuffled so the zeros are interleaved rather than trailing.
rng(seed, 'twister');
A = [randn(m, r), zeros(m, n - r)];
A = A(:, randperm(n));
end

function A = clustered_svals(m, n, seed)
% Random singular subspaces with singular values in [1-1e-2, 1+1e-2].
rng(seed, 'twister');
r = min(m, n);
[U, ~] = qr(randn(m, r), 0);
[V, ~] = qr(randn(n, r), 0);
A = U * diag(1 + 1e-2 * linspace(-1, 1, r)) * V';
end

function A = norm_recompute_case(m, n, seed)
% Two thirds of the columns are pushed to within a 1e-9 relative residual of
% the span of a fixed 6-column basis, so their running tail norm collapses
% after that basis is selected and the sqrt(eps) recompute guard fires.
rng(seed, 'twister');
B = randn(m, 6);
A = randn(m, n);
Pb = B * ((B' * B) \ B');
for j = 1:n
    if mod(j, 3) ~= 0
        res = A(:, j) - Pb * A(:, j);
        A(:, j) = Pb * A(:, j) + 1e-9 * (res / norm(res)) * norm(A(:, j));
    end
end
end

function A = near_deficient_wide(m, n, seed)
% Full-rank leading block, four moderate columns, and 1e-9-scale near-null
% columns, all column-shuffled. Wide (m < n) with near-null tails.
rng(seed, 'twister');
A = [randn(m, 14), randn(m, 4), 1e-9 * randn(m, n - 18)];
A = A(:, randperm(n));
end
