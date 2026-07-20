function tests = test_bsqr_general_bounds
%TEST_BSQR_GENERAL_BOUNDS sigma_min guarantee of the general-A reduction.
%   The guarantee is a k = m statement: selecting all m columns of the
%   orthonormal-row basis G = Qt' (from A' = Qt*Rt) gives
%     sigma_min(A(:,S)) >= sigma_min(A) * sigma_min(G(:,S))
%                       >= sigma_min(A) / sqrt(m*(n-m+1)).
tests = functiontests(localfunctions);
end

function setupOnce(testCase) %#ok<INUSD>
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab_general', 'benchmark'));
addpath(fullfile(repo_root, 'matlab_rand', 'benchmark'));   % rand_test_matrix
end

function testSigmaMinGuarantee(testCase)
sizes = {[8, 100], [16, 200]};
families = {'graded_cols', 'loud_collinear', 'spectrum_needle', ...
    'spectrum_coherent', 'kahan_decoy'};
for fi = 1:numel(families)
    for c = 1:numel(sizes)
        m = sizes{c}(1); n = sizes{c}(2);
        for seed = 1:3
            [A, info] = general_test_matrix(families{fi}, m, n, seed);
            bound = info.sigma_min / sqrt(m * (n - m + 1));
            label = sprintf('%s m=%d n=%d seed=%d', families{fi}, m, n, seed);
            [~, ~, p] = bsqr_general(A, 'pivot_format', 'vector');
            verifyGreaterThanOrEqual(testCase, min(svd(A(:, p(1:m)))), ...
                bound * (1 - 1e-8), ['bsqr selector: ' label]);
            [~, ~, pr] = bsqr_general(A, 'selector', 'bsqr_rand', 'seed', seed, ...
                'pivot_format', 'vector');
            verifyGreaterThanOrEqual(testCase, min(svd(A(:, pr(1:m)))), ...
                bound * (1 - 1e-8), ['bsqr_rand selector: ' label]);
        end
    end
end
end

function testIntermediateInequalityOnG(testCase)
% The selector's own guarantee on the orthonormal-row basis:
% sigma_min(G(:,S)) >= 1/||R11^{-1}||_F >= 1/sqrt(m(n-m+1)).
m = 12; n = 150;
A = general_test_matrix('loud_collinear', m, n, 4);
[Qt, ~] = qr(A', 'econ');
G = Qt';
[~, ~, p] = bsqr_general(A, 'pivot_format', 'vector');
verifyGreaterThanOrEqual(testCase, min(svd(G(:, p(1:m)))), ...
    (1 - 1e-8) / sqrt(m * (n - m + 1)));
end
