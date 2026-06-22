function [A, info] = approx_test_matrix(family, seed, sz)
%APPROX_TEST_MATRIX Application-driven matrices for the low-rank approximation
%   comparison (randomized BSQR vs rejection_rpqr). Unlike rand_test_matrix --
%   which returns a k-by-n matrix with orthonormal rows (already a basis) -- this
%   returns a *full* m-by-n matrix A drawn from an application. The comparison
%   harness then computes accurate leading right singular vectors V_k via svds and
%   selects k of the n columns; the metric is the resulting low-rank
%   approximation error. Each family has a rapidly *decaying* spectrum so that a
%   rank-k approximation (and hence the choice of columns) is meaningful.
%
%   family (char):
%     'gmm_kernel'        - RBF (Gaussian) kernel K_ij = exp(-||x_i-x_j||^2 /
%                           (2 ell^2)) on a Gaussian-mixture point cloud. The
%                           smooth kernel + clustered points give fast spectral
%                           decay; choosing k columns = choosing k Nystrom
%                           landmark points (the canonical kernel-method use).
%                           Symmetric PSD, m = n = N.
%     'integral_skeleton' - A_ij = 1/||x_i - y_j|| between two *separated* point
%                           clouds in a box. A smooth (non-singular, because the
%                           clouds are separated) kernel is numerically low rank;
%                           choosing columns = skeletonization / interpolative
%                           decomposition for fast direct solvers and the FMM.
%     'snapshots'         - Parametric model-reduction snapshots: each column is a
%                           field f(x; mu_j, w_j) sampled on a spatial grid x in
%                           [0,1] (here a Gaussian bump of random center/width).
%                           Smoothness in the parameters makes the snapshot
%                           matrix low rank; choosing columns = empirical
%                           interpolation (DEIM) / representative parameters.
%
%   seed (scalar)  - construction seed (fixed per matrix; the comparison varies
%                    only the *selector* RNG across trials, matching the
%                    fixed-matrix design of the ARP accuracy study).
%   sz (struct, optional) - per-family size overrides; recognized fields:
%     gmm_kernel:        N (points, default 2000), d (features, 10), c (components, 10)
%     integral_skeleton: m (sources, 2000), n (targets, 2000), d (dim, 2), gap (1.5)
%     snapshots:         m (grid points, 1000), n (parameters, 2000)
%
%   info has fields .family, .m, .n, and .desc (one-line description).

if nargin < 2 || isempty(seed); seed = 0; end
if nargin < 3 || isempty(sz); sz = struct(); end
rng(seed);

switch lower(char(family))
    case 'gmm_kernel'
        N = getdef(sz, 'N', 2000);
        d = getdef(sz, 'd', 10);
        c = getdef(sz, 'c', 10);
        mu  = 4 * randn(c, d);                 % well-separated cluster centers
        lab = randi(c, N, 1);
        X   = mu(lab, :) + randn(N, d);        % mixture point cloud
        D2  = sq_dists(X, X);                  % squared pairwise distances
        off = D2(~eye(N, 'logical'));
        ell = sqrt(median(off));               % median-distance bandwidth
        A   = exp(-D2 / (2 * ell^2));          % N-by-N RBF kernel (PSD)
        desc = sprintf('RBF kernel, %d-pt %d-component GMM in R^%d', N, c, d);

    case 'integral_skeleton'
        m   = getdef(sz, 'm', 2000);
        n   = getdef(sz, 'n', 2000);
        d   = getdef(sz, 'd', 2);
        gap = getdef(sz, 'gap', 1.5);          % cloud separation; moderate so the
                                               % spectrum decays over a useful k
                                               % range rather than collapsing at
                                               % very low rank (still well-separated:
                                               % min distance ~ gap, kernel smooth).
        X = rand(m, d);                        % sources in [0,1]^d
        Y = rand(n, d); Y(:, 1) = Y(:, 1) + gap;  % targets shifted away
        A = 1 ./ sqrt(sq_dists(X, Y));         % 1/dist kernel, m-by-n
        desc = sprintf('1/dist kernel, %dx%d separated clouds in R^%d (gap %.1f)', ...
            m, n, d, gap);

    case 'snapshots'
        m = getdef(sz, 'm', 1000);             % spatial grid resolution
        n = getdef(sz, 'n', 2000);             % parameter samples
        x  = linspace(0, 1, m).';
        mu = rand(1, n);                        % bump center in [0,1]
        w  = 0.02 + 0.20 * rand(1, n);          % bump width
        A  = exp(-(x - mu).^2 ./ (2 * w.^2));   % m-by-n, each column a snapshot
        desc = sprintf('Gaussian-bump snapshots, %d grid x %d parameters', m, n);

    otherwise
        error('approx_test_matrix:UnknownFamily', 'Unknown family "%s".', family);
end

info = struct('family', char(family), 'm', size(A, 1), 'n', size(A, 2), 'desc', desc);
end

% ---------------------------------------------------------------------------
function D2 = sq_dists(X, Y)
% Squared Euclidean distances between rows of X (m-by-d) and Y (n-by-d), without
% the Statistics Toolbox. Clamped at 0 (round-off can make the identity term
% slightly negative).
xx = sum(X.^2, 2);
yy = sum(Y.^2, 2).';
D2 = max(0, xx + yy - 2 * (X * Y.'));
end

function v = getdef(s, name, default)
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); else; v = default; end
end
