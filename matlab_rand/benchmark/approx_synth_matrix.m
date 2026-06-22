function [A, info] = approx_synth_matrix(family, seed, sz)
%APPROX_SYNTH_MATRIX Synthetic matrices with a *prescribed* singular spectrum and
%   *prescribed* leverage structure, the companion to approx_test_matrix. The goal
%   is to see how much the R11-conditioning results from the rejection_rpqr
%   comparison (run_rpqr_comparison, families gaussian/spiked_leverage/needle)
%   translate into low-rank *approximation* error.
%
%   A = U * diag(s) * V' is built so that:
%     * the right singular vectors V are exactly the leverage families used in the
%       rpqr comparison: V' = rand_test_matrix(family, r, n, seed) (r-by-n with
%       orthonormal rows), so the column-leverage profile the selectors see is
%       literally the same construction. The leading-k right singular vectors
%       V(:,1:k) inherit that profile for every swept k.
%     * the singular values s are a fixed, deliberately "interesting" spectrum: a
%       few large values, a decay, a flatter (near-plateau) section, then more
%       decay (see synth_spectrum). The three families share this one spectrum, so
%       the optimal rank-k error is identical across them and only the leverage
%       structure differs.
%     * the left singular vectors U are arbitrary orthonormal (they do not affect
%       column selection on V' nor the unitarily-invariant projection error).
%
%   family (char): any rand_test_matrix family. 'gaussian' | 'spiked_leverage' |
%     'needle' mirror the rpqr comparison; 'coherent' | 'collinear_cluster' are the
%     R11-conditioning stress cases used by run_approx_cond_comparison. r plays the
%     role of "k" in rand_test_matrix's recipes.
%   seed (scalar): construction seed (fixed per matrix).
%   sz (struct, optional): fields m (rows, default 2000), n (cols, default 2000),
%     r (numerical rank / spectrum length, default 300; needs r < min(m,n) and
%     r > max swept k so every k has columns to drop and a sigma_{k+1}).
%
%   info has fields .family, .m, .n, .r, .desc.

if nargin < 2 || isempty(seed); seed = 0; end
if nargin < 3 || isempty(sz); sz = struct(); end
m = getdef(sz, 'm', 2000);
n = getdef(sz, 'n', 2000);
r = getdef(sz, 'r', 300);
assert(r < min(m, n), 'approx_synth_matrix: need r < min(m,n).');

s = synth_spectrum(r);                          % r-by-1 prescribed singular values
M = rand_test_matrix(family, r, n, seed);       % r-by-n, orthonormal rows = V'
rng(seed + 99991);                              % independent stream for U
U = orth(randn(m, r));                          % arbitrary orthonormal left vectors
A = U * (s .* M);                               % = U diag(s) V', rank r

info = struct('family', char(family), 'm', m, 'n', n, 'r', r, 'desc', ...
    sprintf('%s leverage, prescribed spectrum, rank %d (%dx%d)', ...
    char(family), r, m, n));
end

% ---------------------------------------------------------------------------
function s = synth_spectrum(r)
% A deliberately structured, strictly decreasing spectrum: a few large values, a
% pronounced first drop, a near-flat "plateau" section, then more decay. The drop
% lands well below the dominant values (0.9 -> 0.02, ~45x) so the top directions
% clearly dominate -- this reads as a steep initial decline in the error-vs-k
% curve. The plateau declines gently (0.02 -> 0.015) rather than being exactly flat
% so the singular values stay distinct (no degenerate subspace -> V is exactly the
% prescribed family vectors). Breakpoints scale with r.
b1 = max(3, round(0.015 * r));                  % a few large
b2 = max(b1 + 5, round(0.13 * r));              % end of the (pronounced) first drop
b3 = max(b2 + 5, round(0.42 * r));              % end of the flatter section
s = zeros(r, 1);
s(1:b1)      = linspace(1.0, 0.9, b1);                          % few large
s(b1+1:b2)   = logspace(log10(0.9),  log10(0.02),  b2 - b1);    % pronounced first drop
s(b2+1:b3)   = logspace(log10(0.02), log10(0.015), b3 - b2);    % flatter section
s(b3+1:r)    = logspace(log10(0.015), log10(1e-5), r - b3);     % more decay
end

function v = getdef(s, name, default)
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); else; v = default; end
end
