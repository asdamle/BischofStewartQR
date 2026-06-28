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

info = struct('family', char(family), 'm', m, 'n', n, 'r', r, 'svals', s, 'desc', ...
    sprintf('%s leverage, prescribed spectrum, rank %d (%dx%d)', ...
    char(family), r, m, n));
end

% ---------------------------------------------------------------------------
function s = synth_spectrum(r)
% A deliberately structured, strictly decreasing spectrum: a steep initial decay
% (the top ~5% of directions carry most of the energy, so the relative rank-k error
% drops below 1 quickly), a gentle decay, then a SHARP cliff (a few indices, ~80x
% drop) at k ~ 0.167*r, then a gentle tail. The optimal rank-k error falls steeply
% at first and plummets across the cliff, giving the rank-k ID-error rows a clear
% knee there. Strictly decreasing so the singular values stay distinct (no
% degenerate subspace -> V is exactly the prescribed family vectors). Breakpoints
% scale with r.
b1 = max(3, round(0.05 * r));                   % steep initial decay end
b2 = max(b1 + 5, round(0.167 * r));             % cliff start (~k=50 for r=300)
bc = min(r - 1, b2 + max(3, round(0.025 * r))); % cliff end (sharp: a few indices)
s = zeros(r, 1);
s(1:b1)      = logspace(log10(1.0),  log10(0.1),  b1);          % steep initial decay
s(b1+1:b2)   = logspace(log10(0.1),  log10(0.04), b2 - b1);     % gentle decay to the cliff
s(b2+1:bc)   = logspace(log10(0.04), log10(5e-4), bc - b2);     % SHARP cliff (~80x)
s(bc+1:r)    = logspace(log10(5e-4), log10(1e-6), r - bc);      % gentle tail
end

function v = getdef(s, name, default)
if isfield(s, name) && ~isempty(s.(name)); v = s.(name); else; v = default; end
end
