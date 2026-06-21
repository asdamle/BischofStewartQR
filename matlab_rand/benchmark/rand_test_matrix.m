function [M, info] = rand_test_matrix(family, k, n, seed)
%RAND_TEST_MATRIX Generate k-by-n matrices with orthonormal rows (M*M'=I_k)
%   that stress the randomized column-selection scheme in different ways.
%
%   The randomized sampler's cost is governed by how many remaining columns are
%   "acceptable" (criterion c_j below the running threshold) at each step. That
%   fraction is large for benign inputs and small when the useful directions are
%   concentrated in a few columns -- the leverage scores ell_j = ||M(:,j)||^2
%   (which sum to k) are exactly the column quantities the sampler sees first.
%   The families below span uniform to highly concentrated leverage.
%
%   family (char):
%     'gaussian'        - orth(randn): near-uniform leverage (benign).
%     'graded_leverage' - leverage decays smoothly across columns.
%     'spiked_leverage' - a handful of high-leverage columns, rest flat.
%     'coherent'        - columns cluster into a few groups (many redundant).
%     'needle'          - ~k useful columns hidden among n near-null ones
%                         (hardest for uniform sampling; the case where
%                         norm-weighted sampling should help most).
%     'chebyshev'       - applied structure: first k Chebyshev polynomials
%                         sampled at n equispaced points (a polynomial
%                         least-squares / interpolation design matrix).
%                         Selecting columns = choosing k nodes from n; leverage
%                         (the Christoffel function) concentrates near the ends.
%
%   info has fields .family, .k, .n, .leverage (1-by-n column squared norms).

if nargin < 4 || isempty(seed)
    seed = 0;
end
rng(seed);

switch lower(char(family))
    case 'gaussian'
        W = randn(n, k);

    case 'graded_leverage'
        % Smoothly decaying row scaling -> graded leverage across columns.
        d = logspace(0, -3, n).';            % 1 ... 1e-3
        W = (d .* randn(n, k));

    case 'spiked_leverage'
        % A few dominant columns carry most of the leverage.
        nspike = max(1, ceil(1.5 * k));
        d = ones(n, 1);
        idx = randperm(n, min(nspike, n));
        d(idx) = 100;
        W = (d .* randn(n, k));

    case 'coherent'
        % Rows of the factor cluster into r groups -> columns of M cluster,
        % so most columns are near-duplicates of a few representatives.
        r = max(2, ceil(k / 4));
        centers = randn(r, k);
        assign = randi(r, n, 1);
        W = centers(assign, :) + 0.05 * randn(n, k);

    case 'needle'
        % ng ~ k useful columns; the rest are near-null (tiny leverage).
        ng = min(n, k + ceil(0.25 * k));
        d = 1e-6 * ones(n, 1);
        idx = randperm(n, ng);
        d(idx) = 1;
        W = (d .* randn(n, k));

    case 'chebyshev'
        % Applied structure: a degree-(k-1) Chebyshev design matrix on n random
        % sample points in [-1,1] -- polynomial least-squares / optimal design,
        % where selecting k of the n columns = picking the k most informative
        % samples. Leverage (the Christoffel function) concentrates near the
        % ends. Chebyshev polynomials are bounded by 1 on [-1,1] so the factor
        % is well-scaled; orth makes the rows orthonormal (rank k).
        x = sort(2 * rand(1, n) - 1);          % seeded random design points
        B = ones(k, n);
        if k >= 2; B(2, :) = x; end
        for i = 3:k
            B(i, :) = 2 * x .* B(i - 1, :) - B(i - 2, :);   % Chebyshev recurrence
        end
        W = B.';

    otherwise
        error('rand_test_matrix:UnknownFamily', 'Unknown family "%s".', family);
end

% Orthonormal columns of the n-by-k factor -> orthonormal rows of M (k-by-n).
Q = orth(W);
if size(Q, 2) < k
    % Degenerate factor (e.g. coherent with tiny rank); pad to full row rank.
    Q = orth([Q, randn(n, k - size(Q, 2))]);
end
M = Q(:, 1:k).';

info = struct('family', char(family), 'k', k, 'n', n, ...
    'leverage', sum(M.^2, 1));
end
