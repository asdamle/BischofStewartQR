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
    'frob_check', true, 'q_check', true, 'orthonormal_rows', false, ...
    'gap_tol', 1e-6, 'rtol_R', 1e-10, 'rtol_rinv', 2e-8);
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
