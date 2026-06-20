function varargout = bsqr_rand(A, varargin)
%BSQR_RAND Randomized Bischof-Stewart column selection.
%   [P, REFLECTORS, R11] = BSQR_RAND(A) selects k = min(size(A)) columns of A
%   by the randomized acceptance rule and returns:
%     P          - permutation row vector; P(1:k) are the selected columns
%     REFLECTORS - struct('V', m-by-k unit-diagonal reflector store,
%                          'tau', k-by-1, 'm', m, 'k', k); use BSQR_RAND_FORMQ
%     R11        - k-by-k upper-triangular factor of A(:,P(1:k))
%
%   [P, REFLECTORS, R11, STATS]      also returns instrumentation (see below).
%   [P, REFLECTORS, R11, STATS, R12] also returns R12 = Q(:,1:k)'*A(:,P(k+1:n)),
%       but ONLY when 'return_r12' is true (it costs an extra O(n*k^2) pass and
%       is off by default).
%
%   By default the kernel runs in-block ('batched', see below): each sampled
%   block is brought to the current frame once, then BSQR is run within it to
%   select as many columns as meet the bound before resampling. This amortizes
%   the per-block reflector apply over many selections (O(k^3) rather than the
%   single-select O(k^4)) at a modest cost in realized conditioning.
%
% Name-value options:
%   'k'              - number of columns to select (default min(m,n))
%   'batched'        - logical (default true). true: in-block BSQR (many
%                      selections per sampled block; the fast default). false:
%                      single-select (one selection per block; tighter realized
%                      conditioning, more reflector applies). Both maintain the
%                      same per-step ||R11^{-1}||_F guarantee.
%   'block_size'     - candidates per sampled block (default k when batched, else
%                      ceil(k/2) clamped to [16,64]; see RAND_DEFAULT_BLOCK).
%                      Larger blocks improve realized selection quality (the
%                      per-step minimum is taken over more candidates) at higher
%                      cost; the guaranteed bound is unchanged.
%   'threshold_mode' - 'running_mean' (default, per-singular-value control) or
%                      'worstcase_allowance' (more permissive, fewer samples,
%                      only the final Frobenius bound -- not the per-step bound)
%   'slack'          - multiplier (>=1) loosening the threshold (default 1.0;
%                      values >1 trade conditioning for fewer samples and no
%                      longer guarantee the Osinsky bound -- experimental)
%   'sampling'       - 'uniform' (default) or 'normweighted' (by starting
%                      squared column norms; adds an O(m*n) precompute)
%   'pick'           - single-select only ('batched' false): 'best_in_block'
%                      (default) accepts the minimum-criterion column in the
%                      block when it meets the threshold; 'first' accepts the
%                      first that meets it. (Batched always takes the in-block
%                      minimizer.)
%   'seed'           - RNG seed for reproducibility (default [], leaves rng)
%   'return_r12'     - logical, compute R12 as a 5th output (default false)
%   'backend'        - 'auto' (default), 'mfile', or 'mex'
%   'check_finite'   - validate inputs are finite (default true)
%
% STATS fields (per-step arrays are 1-by-k, indexed by selection step):
%   f2             - running ||R11^{-1}||_F^2 after each step
%   Fhat           - deterministic worst-case bound Fhat_i = i(n-i+1)/(k-i+1)
%   crit           - c-increment of the accepted pivot
%   threshold      - acceptance threshold used at each step
%   samples_tested - candidate columns evaluated at each step (batched: the
%                    block's size is attributed to its first selection, 0 after)
%   rounds         - sampling rounds (blocks) at each step (batched: 1 at a
%                    block's first selection, 0 after)
%   fallback       - true where the exhaustive global-min fallback fired
%   frob_inv       - final ||R11^{-1}||_F = sqrt(f2(end))
%   osinsky_bound  - sqrt(k(n-k+1)), the Frobenius guarantee target
%   total_tested   - sum of samples_tested
%   blocks_sampled - total sampled blocks / reflector applies (sum of rounds)
%
%   The per-step bound and its guarantees assume A has orthonormal rows
%   (the GKS setting, m = k). The algorithm still runs for general A, but
%   only the orthonormal-row case is covered by the theory.
%
%   See also BSQR_RAND_FORMQ, docs/RANDOMIZED_BSQR_PLAN.md.

if nargout > 5
    error('bsqr_rand:TooManyOutputs', 'bsqr_rand supports at most 5 outputs.');
end

opts = bsqr_rand_parse_options(A, varargin{:});

if nargout >= 5 && ~opts.return_r12
    error('bsqr_rand:R12NotRequested', ...
        'The fifth output (R12) requires the option return_r12=true.');
end

nout = max(nargout, 1);
out = cell(1, nout);
if should_use_mex(opts.backend)
    ensure_bsqr_rand_mex_ready();
    [out{:}] = bsqr_rand_mex(A, varargin{:});
else
    [out{:}] = bsqr_rand_mfile(A, opts);
end
varargout = out;
end

function tf = should_use_mex(backend)
switch backend
    case 'auto'
        tf = bsqr_rand_mex_available();
    case 'mfile'
        tf = false;
    case 'mex'
        tf = true;
    otherwise
        error('bsqr_rand:InvalidBackend', 'backend must be "auto", "mfile", or "mex".');
end
end

function ensure_bsqr_rand_mex_ready()
persistent mex_ready
if ~isempty(mex_ready) && mex_ready
    return;
end

thisdir = fileparts(mfilename('fullpath'));
mexdir = fullfile(thisdir, 'mex');
if isfolder(mexdir)
    addpath(mexdir);
end

if ~bsqr_rand_mex_available()
    error('bsqr_rand:MexUnavailable', ...
        'backend="mex" requested but bsqr_rand_mex is not available. Run matlab_rand/build_bsqr_rand_mex.m first.');
end
mex_ready = true;
end
