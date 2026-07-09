function theta = bsqr_rand_threshold(f2, nsel, k, n, mode)
%BSQR_RAND_THRESHOLD Per-step acceptance threshold on the c-increment.
%   A candidate column at the step that grows the selected set from NSEL to
%   NSEL+1 columns is acceptable when its increment
%       c = (1 + ||w||^2) / rho^2
%   to the squared inverse Frobenius norm f2 = ||R11^{-1}||_F^2 satisfies
%   c <= THETA. F2 is the *squared* quantity (see docs/RANDOMIZED_BSQR_PLAN.md).
%
%   mode = 'running_mean'  (default; the per-step acceptance threshold)
%       theta = (f2 + n - 2*nsel) / (k - nsel)
%       This equals the rho^2-weighted mean of c over the remaining columns
%       (for orthonormal-row input), so the minimizer always qualifies and
%       the per-step bound f_{i+1} <= ((k-i+1)/(k-i)) f_i + (n-2i)/(k-i)
%       holds, giving Osinsky's per-singular-value guarantees.
%
%   mode = 'worstcase_allowance'
%       theta = Fhat_{nsel+1} - f2,  Fhat_m = m(n-m+1)/(k-m+1)
%       More permissive (spends accumulated slack between the running f2 and
%       the deterministic worst case Fhat), still provably maintains
%       ||R11^{-1}||_F^2 <= k(n-k+1). Typically accepts with fewer samples.

switch mode
    case 'running_mean'
        theta = (f2 + n - 2 * nsel) / (k - nsel);
    case 'worstcase_allowance'
        Fhat_next = (nsel + 1) * (n - nsel) / (k - nsel);
        theta = Fhat_next - f2;
    otherwise
        error('bsqr_rand:InvalidThresholdMode', 'Unknown threshold_mode "%s".', mode);
end
end
