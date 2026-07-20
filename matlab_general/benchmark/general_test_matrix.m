function [A, info] = general_test_matrix(family, m, n, seed)
%GENERAL_TEST_MATRIX Generate general (non-orthonormal-row) m-by-n test
%   matrices, m <= n, full row rank, that stress DIRECT column-pivoted
%   selection on A, while the row-space reduction (bsqr_general) carries its
%   sigma_min guarantee regardless of scaling.
%
%   family (char):
%     'graded_cols'       - randn with column norms graded over 8 decades
%                           (logspace column weights): norm-driven pivot
%                           criteria see the scaling, not the conditioning.
%     'loud_collinear'    - a cluster of very large-norm, nearly collinear
%                           columns hidden among modest well-spread background
%                           columns; greedy norm pivoting is drawn to the
%                           cluster first.
%     'spectrum_needle'   - A = U*diag(logspace(0,-5,m))*V with V the 'needle'
%                           orthonormal-row family from rand_test_matrix
%                           (~m useful columns among near-null ones), U a
%                           random orthogonal factor.
%     'spectrum_coherent' - the same graded spectrum with the 'coherent' V
%                           (columns cluster into a few groups).
%     'kahan_decoy'       - an m-by-m Kahan-type upper-triangular block (the
%                           classical greedy-pivoting failure: column norms
%                           decay geometrically but the block is exponentially
%                           closer to singular than its norms suggest),
%                           followed by smaller-norm well-conditioned decoy
%                           columns that a guarantee-bearing selection must
%                           mix in. Greedy norm pivoting takes all m Kahan
%                           columns and lands on sigma_min(K) << sigma_min(A).
%
%   info fields: .family, .m, .n, .sigma_min (min singular value of A, from
%   one svd), .desc.
%
%   The spectrum_* families reuse matlab_rand/benchmark/rand_test_matrix.m
%   (must be on the path); this suite is otherwise self-contained.

if nargin < 4 || isempty(seed)
    seed = 0;
end

switch lower(char(family))
    case 'graded_cols'
        rng(seed);
        A = randn(m, n) .* logspace(0, -8, n);
        desc = 'randn with column norms graded over 8 decades';

    case 'loud_collinear'
        rng(seed);
        ncl = max(2, round(0.2 * n));
        v0 = randn(m, 1); v0 = v0 / norm(v0);
        A = randn(m, n);
        A(:, 1:ncl) = 1e3 * (v0 + 1e-4 * randn(m, ncl));
        A = A(:, randperm(n));
        desc = 'large-norm near-collinear cluster over modest background';

    case {'spectrum_needle', 'spectrum_coherent'}
        vfam = strrep(lower(char(family)), 'spectrum_', '');
        V = rand_test_matrix(vfam, m, n, seed);   % m-by-n orthonormal rows (seeds rng)
        [U, ~] = qr(randn(m), 'econ');            % continues the seeded stream
        A = U * (logspace(0, -5, m).' .* V);
        desc = sprintf('graded spectrum (5 decades) with %s right factor', vfam);

    case 'kahan_decoy'
        rng(seed);
        % Kahan block: diag(c.^(0:m-1)) * (I - s*triu(ones,1)), scaled so the
        % last column norm is c^(m-1) = 0.1. Column norms decay strictly
        % (a tiny extra grading breaks rounding ties toward no exchange), so
        % greedy norm pivoting accepts the columns left to right, yet
        % sigma_min(K) is smaller than c^(m-1) by a factor exponential in m.
        c = 0.1^(1 / max(m - 1, 1));
        s = sqrt(1 - c^2);
        K = diag(c.^(0:m-1)) * (eye(m) - s * triu(ones(m), 1));
        K = K * diag((1 - 1e-10).^(0:m-1));
        % Decoys: well-conditioned random columns at norm just below K's
        % smallest column, so they are picked only after all of K.
        delta = 0.5 * c^(m - 1);
        D = (delta / sqrt(m)) * randn(m, n - m);
        A = [K, D];
        desc = 'Kahan block (greedy-pivot trap) plus well-conditioned decoys';

    otherwise
        error('general_test_matrix:UnknownFamily', 'Unknown family "%s".', family);
end

info = struct('family', char(family), 'm', m, 'n', n, ...
    'sigma_min', min(svd(A)), 'desc', desc);
end
