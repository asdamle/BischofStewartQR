function out = oracle_bsqr(A, k)
%ORACLE_BSQR Literal-transcription oracle for Bischof-Stewart pivoted QR.
%   OUT = ORACLE_BSQR(A)    runs min(size(A)) pivot steps.
%   OUT = ORACLE_BSQR(A, K) runs K pivot steps.
%
%   This is the V1 validation oracle from docs/VALIDATION_AND_PERF_PLAN.md:
%   a deliberately naive transcription of Algorithm A.1 (Appendix A) of the
%   manuscript, notes/GKSevolved_draft.tex. Each step recomputes from
%   scratch the quantities the production kernels maintain incrementally:
%
%     w_j = R11^{-1} R12(:,j)   by triangular solve   (Alg. A.1 line 9 state)
%     s_j = ||A(i:m, j)||^2     exact tail norms      (Alg. A.1 lines 1, 8)
%
%   so it shares no recurrence (and therefore no shared recurrence bug)
%   with bsqr_mfile, bsqr_mex, or the Julia kernel. Householder steps use
%   explicit reflector matrices H = I - 2*v*v'/(v'*v). Cost is roughly
%   O(n*k^3); this is a test oracle, never a benchmark subject.
%
%   Output struct fields:
%     R          k-by-n upper trapezoidal factor
%     T          m-by-n fully reduced matrix: T(1:k,:) carries R and
%                T(k+1:m,:) is the unreduced tail block, so A(:,p) == Q*T
%                even when k < min(m, n)
%     Q          m-by-m orthogonal factor
%     p          1-by-n permutation vector; p(1:k) is the pivot sequence
%     rinv_r12   k-by-(n-k) solve R11 \ R12 from the final factor
%     crit_best  1-by-k criterion value c(j*) of the selected pivot per step
%                (running sum reproduces ||R11^{-1}||_F^2, manuscript
%                eq. (2.4))
%     crit_gap   1-by-k relative gap (c_runnerup - c_best)/c_best per step;
%                Inf when a step has a single candidate. Tests use this to
%                screen near-ties before demanding exact pivot equality.
%
%   Documented deviations from the literal Algorithm A.1 (the same
%   sanctioned deviations the production kernels make; see
%   docs/VALIDATION.md):
%
%   * Householder sign: Alg. A.1 leaves the reflector sign free
%     (+-||x||e_1). We fix the numerically stable LAPACK convention
%     rho_i = -sign(x_1)*||x||, and
%     the identity reflector (rho_i = x_1) when x has nothing below its
%     first entry. Mathematically equivalent: flipping the sign of row t
%     flips it in both R11 and R12, and (D*R11)\(D*R12) = R11\R12 for any
%     D = diag(+-1), so every w_j -- and hence the criterion -- is
%     unchanged; Q absorbs D. Using the same convention as the kernels
%     makes R/Q directly comparable, with no sign normalization in tests.
%   * Ties: strict first-minimum (Alg. A.1 allows any tie rule, e.g. the
%     smaller index; all implementations fix first-minimum order, and
%     MATLAB's min() returns the first minimizer).
%
%   The oracle errors if every remaining column has exactly zero tail norm
%   before completing k steps: Algorithm A.1 requires k <= rank(A), and
%   the production kernels' behavior past that point is an extension
%   covered by contract tests, not by this oracle.

if nargin < 2
    k = min(size(A));
end
validateattributes(A, {'double'}, {'2d', 'real', 'finite'}, mfilename, 'A');
validateattributes(k, {'numeric'}, ...
    {'scalar', 'integer', 'nonnegative', '<=', min(size(A))}, mfilename, 'k');
k = double(k);

Awork = A;
[m, n] = size(Awork);

p = 1:n;                                    % Alg. 1, line 1
Q = eye(m);
crit_best = zeros(1, k);
crit_gap = zeros(1, k);

% Alg. 1, lines 2-3 are storage initialization for the incrementally
% maintained W and s. The oracle keeps neither: both are recomputed from
% scratch inside the loop, which is the entire point of this implementation.

for i = 1:k                                 % Alg. 1, line 4
    % --- Pivot selection (Alg. 1, lines 5-8) ---
    % R11 is the leading (i-1)x(i-1) triangle of the reduced matrix. It is
    % genuinely upper triangular here because each completed step zeroed
    % its subdiagonal explicitly (see below); triu() only asserts that.
    R11 = triu(Awork(1:i-1, 1:i-1));

    c = inf(1, n);
    for j = i:n
        s_j = norm(Awork(i:m, j))^2;        % exact tail norm rho_j^2, from scratch
        if s_j > 0
            w_j = R11 \ Awork(1:i-1, j);    % Alg. A.1's w_j, from scratch
            c(j) = (1 + norm(w_j)^2) / s_j; % Alg. A.1 line 4; at i = 1 the empty
        end                                 % solve gives 1/s_j (||w_j|| = 0)
        % s_j == 0 leaves c(j) = Inf, matching the kernels' guard.
    end

    [cbest, rel_idx] = min(c(i:n));         % line 8; min() takes the FIRST minimum
    jstar = i - 1 + rel_idx;

    if ~isfinite(cbest)
        error('oracle_bsqr:RankDeficient', ...
            ['All remaining columns have zero tail norm at step %d; ', ...
             'Algorithm 1 requires full-rank progress through k steps.'], i);
    end

    crit_best(i) = cbest;
    cands = sort(c(i:n));
    if numel(cands) > 1 && isfinite(cands(2))
        crit_gap(i) = (cands(2) - cands(1)) / cands(1);
    else
        crit_gap(i) = Inf;
    end

    % --- Swap columns i and j* (Alg. 1, lines 9-10) ---
    % Only A and p exist here; the W and s swaps of line 10 are vacuous
    % because the oracle recomputes both from scratch each step.
    if jstar ~= i
        Awork(:, [i, jstar]) = Awork(:, [jstar, i]);
        p([i, jstar]) = p([jstar, i]);
    end

    % --- Householder reduction of column i (Alg. 1, lines 11-12) ---
    [H, rho] = householder_matrix(Awork(i:m, i));
    Awork(i:m, i:n) = H * Awork(i:m, i:n);
    Q(:, i:m) = Q(:, i:m) * H;              % H symmetric: Q = H_1*...*H_k

    % In exact arithmetic H*x = rho*e_1, so the subdiagonal of column i is
    % exactly zero. Assert the computed residue is rounding-level, then
    % enforce the exact-arithmetic value so R11 stays exactly triangular.
    subdiag = Awork(i+1:m, i);
    if ~isempty(subdiag)
        assert(norm(subdiag) <= 1e3 * eps * max(abs(rho), realmin), ...
            'oracle_bsqr:ReflectorResidue', ...
            'Householder application left a non-rounding subdiagonal at step %d.', i);
        Awork(i+1:m, i) = 0;
    end

    % Alg. 1, lines 13-19 (incremental W and s updates) are intentionally
    % absent: the next iteration recomputes both quantities from scratch.
end                                          % Alg. 1, line 20 returns p

out = struct();
out.p = p;
out.Q = Q;
out.T = Awork;
out.R = triu(Awork(1:k, :));
out.rinv_r12 = Awork(1:k, 1:k) \ Awork(1:k, k+1:n);
out.crit_best = crit_best;
out.crit_gap = crit_gap;
end

function [H, rho] = householder_matrix(x)
%HOUSEHOLDER_MATRIX Explicit symmetric reflector with H*x = rho*e_1.
%   Stable sign choice rho = -sign(x_1)*||x||, so v = x - rho*e_1 has no
%   cancellation in its first entry. When x has zero norm below its first
%   entry there is nothing to reflect: H = I and rho = x_1, matching
%   LAPACK dlarfg (tau = 0) and the production kernels.
len = numel(x);
if norm(x(2:end)) == 0
    H = eye(len);
    rho = x(1);
    return;
end

sgn = sign(x(1));
if sgn == 0
    sgn = 1;
end
rho = -sgn * norm(x);

v = x;
v(1) = x(1) - rho;
H = eye(len) - (2 / (v' * v)) * (v * v');
end
