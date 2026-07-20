function varargout = bsqr_general(A, varargin)
%BSQR_GENERAL Bischof-Stewart column selection for general short-wide A.
%   R = BSQR_GENERAL(A)
%   [Q,R] = BSQR_GENERAL(A)
%   [Q,R,E_OR_P] = BSQR_GENERAL(A)
%   [Q,R,E_OR_P,STATS] = BSQR_GENERAL(A)
%
%   For A m-by-n with m <= n and full row rank, selects columns of A with the
%   Bischof-Stewart guarantee transplanted from the orthonormal-row setting:
%   first the reduced QR factorization A' = Qt*Rt is computed, then the
%   selector (bsqr or bsqr_rand) runs on G = Qt' (m-by-n, orthonormal rows),
%   where its theory applies. Since A = Rt'*G, at k = m the selected columns
%   satisfy
%
%       sigma_min(A(:,p(1:m))) >= sigma_min(A) * sigma_min(G(:,p(1:m)))
%                              >= sigma_min(A) / sqrt(m*(n-m+1)).
%
%   Early stopping (k < m) is supported mechanically but carries no such
%   bound (the guarantee is a k = m statement).
%
% Outputs (for A m-by-n and k selection steps, default k = m):
%   R      - k-by-n upper-trapezoidal factor of the PERMUTED columns; for the
%            default k = m, A(:,p) = Q*R (equivalently A*E = Q*R). With early
%            stop (k < m) that holds for the selected block,
%            A(:,p(1:k)) = Q*R(:,1:k), and the trailing block is
%            R(:,k+1:n) = Q'*A(:,p(k+1:n)) (matching bsqr's semantics).
%   Q      - m-by-k economy factor with orthonormal columns.
%   E_OR_P - the permutation, as an n-by-n permutation matrix E with
%            A*E = A(:,p) ('pivot_format','matrix', the default, matching
%            bsqr) or as the 1-by-n index row p itself
%            ('pivot_format','vector').
%   STATS  - struct with fields
%              selector       - 'bsqr' or 'bsqr_rand'
%              k              - number of selection steps
%              t_qr_At        - seconds in phase 1, qr(A','econ')
%              t_select       - seconds in phase 2, the selector on G = Qt'
%              t_assemble     - seconds in phase 3, reassembling Q and R for A
%              t_total        - seconds across all three phases
%              rt_diag_ratio  - min|diag(Rt)| / max|diag(Rt)| (rank indicator)
%              selector_stats - bsqr_rand's STATS struct ([] for bsqr / k = 0)
%
% Name-value options (owned by the wrapper):
%   'selector'     - 'bsqr' (default) or 'bsqr_rand'
%   'k'            - number of selection steps (default m). The sigma_min
%                    guarantee is a k = m statement; k < m is mechanical only.
%   'pivot_format' - 'matrix' (default) or 'vector'
%   'check_finite' - validate that A is finite (default true). The wrapper
%                    scans A once and forwards check_finite=false to the
%                    selector (G derived from finite A is finite).
%   'rank_tol'     - relative diag(Rt) threshold (default max(m,n)*eps) below
%                    which a bsqr_general:IllConditioned warning is issued:
%                    the factorization stays valid, but the selection
%                    guarantee degrades as sigma_min(A) -> 0.
%
%   Any other name-value pair is forwarded verbatim to the selected selector
%   ('backend', 'norm_recomp_tol', and for bsqr_rand 'seed', 'sampling',
%   'block_size', 'batched', 'threshold_mode', 'slack', 'pick'); unknown
%   names are rejected by the selector's own parser. The selectors'
%   output-shape options ('return_rinv_r12', 'return_r12', 'trace') are
%   blocked -- the wrapper owns its output contract.
%
%   This layer intentionally CALLS matlab/bsqr.m and matlab_rand/bsqr_rand.m
%   (both must be on the path; run startup.m from the repo root). It is pure
%   MATLAB: the heavy work happens inside the built-in qr and the selectors'
%   MEX backends.

if nargout > 4
    error('bsqr_general:TooManyOutputs', 'bsqr_general supports at most 4 outputs.');
end

[opts, fwd] = bsqr_general_parse_options(A, varargin{:});

% Normalize to double up front, matching bsqr/bsqr_rand's dispatchers.
if ~isa(A, 'double')
    A = double(A);
end

[m, n] = size(A);
k = opts.k;

t_all = tic;

% --- Phase 1: reduced QR of A^T; G = Qt' has orthonormal rows and
%     A = Rt' * G, so selecting columns of G selects columns of A.
t0 = tic;
[Qt, Rt] = qr(A', 'econ');
t_qr_At = toc(t0);

d = abs(diag(Rt));
dmax = max([d; 0]);
if dmax == 0
    rt_ratio = 0;
else
    rt_ratio = min(d) / dmax;
end
if rt_ratio <= opts.rank_tol
    warning('bsqr_general:IllConditioned', ...
        ['A appears (near-)rank-deficient: min|diag(Rt)|/max|diag(Rt)| = %.3e ', ...
         '<= rank_tol = %.3e. The factorization is still returned, but the ', ...
         'subset-selection guarantee degrades as sigma_min(A) -> 0.'], ...
        rt_ratio, opts.rank_tol);
end

% --- Phase 2: run the selector on the orthonormal-row basis G.
t0 = tic;
G = Qt';
selector_stats = [];
if k == 0
    Qg = zeros(m, 0);
    Rg = zeros(0, n);
    p = 1:n;
elseif opts.selector == "bsqr"
    [Qg, Rg, p] = bsqr(G, 'k', k, 'pivot_format', 'vector', ...
        'check_finite', false, fwd{:});
else  % "bsqr_rand"
    [p, Qg, R11g, selector_stats, R12g] = bsqr_rand(G, 'k', k, ...
        'return_r12', true, 'check_finite', false, fwd{:});
    Rg = [R11g, R12g];
end
t_select = toc(t0);

% --- Phase 3: reassemble a genuine economy QR of the permuted A. From
%     A(:,p) = Rt'*G(:,p) = (Rt'*Qg)*Rg, the small QR Rt'*Qg = Qs*Rs
%     (O(m^2 k)) gives A(:,p) = Qs*(Rs*Rg) with Rs*Rg upper trapezoidal.
t0 = tic;
B = Rt' * Qg;
[Qs, Rs] = qr(B, 'econ');
Q = Qs;
R = Rs * Rg;
if k < m
    % Rs*Rg(:,k+1:n) is an oblique coefficient block when the selection
    % stopped early; overwrite with bsqr's documented trailing-block
    % semantics (for k = m the two coincide and this is skipped).
    R(:, k+1:n) = Q' * A(:, p(k+1:n));
end
t_assemble = toc(t0);

stats = struct('selector', char(opts.selector), 'k', k, ...
    't_qr_At', t_qr_At, 't_select', t_select, 't_assemble', t_assemble, ...
    't_total', toc(t_all), 'rt_diag_ratio', rt_ratio, ...
    'selector_stats', selector_stats);

if nargout <= 1
    varargout{1} = R;
    return;
end
varargout{1} = Q;
varargout{2} = R;
if nargout >= 3
    varargout{3} = format_pivot_output(p, opts.pivot_format, n);
end
if nargout >= 4
    varargout{4} = stats;
end

end

function piv = format_pivot_output(p, pivot_format, n)
% Mirrors bsqr's pivot formatting: the n-by-n permutation matrix E with
% A*E = A(:,p) by default, or the index row itself.
if pivot_format == "vector"
    piv = p;
    return;
end
piv = eye(n);
piv = piv(:, p);
end
