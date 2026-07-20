function tests = test_bsqr_rand_pivot_range
%TEST_BSQR_RAND_PIVOT_RANGE Regression: MEX sampler must stay index-safe on
% concentrated-norm general input.
%
% The MEX's norm-weighted Fenwick sampler keeps a separately rounded running
% total (cur_total); on inputs whose squared column norms span many orders of
% magnitude (here ~[1e-15, 1e-1]), remove/restore traffic of the large weights
% leaves cur_total above the tree's own prefix total, and find() then walked
% one past the end (returned pos == n), yielding pivot index n+1 in the output
% and a one-byte out-of-bounds write (taken[n] = 1) that corrupted the heap.
% This input family is outside the documented orthonormal-row theory but must
% be memory-safe. The m-file backend samples over an explicit index list and
% is index-safe by construction, so only the MEX is exercised here.
%
% Matches matlab_general/benchmark/general_test_matrix.m 'spectrum_needle'
% (rebuilt inline from rand_test_matrix so this suite stays self-contained).
tests = functiontests(localfunctions);
end

function A = spectrum_needle(m, n, seed)
V = rand_test_matrix('needle', m, n, seed);   % orthonormal rows (seeds rng)
[U, ~] = qr(randn(m), 'econ');                % continues the seeded stream
A = U * (logspace(0, -5, m).' .* V);
end

function testConcentratedNormPivotsAreValid(testCase)
if ~bsqr_rand_mex_available()
    return;
end
sizes = {[64, 500], [64, 1000]};
for c = 1:numel(sizes)
    m = sizes{c}(1); n = sizes{c}(2);
    for mat_seed = 0:4
        A = spectrum_needle(m, n, mat_seed);
        for s = 0:4
            p = bsqr_rand(A, 'k', m, 'seed', s, 'backend', 'mex');
            verifyEqual(testCase, sort(p), 1:n, sprintf( ...
                'invalid pivot vector (m=%d n=%d mat_seed=%d seed=%d): max %d', ...
                m, n, mat_seed, s, max(p)));
        end
    end
end
end
