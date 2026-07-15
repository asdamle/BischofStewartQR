function opts = bsqr_rand_parse_options(A, varargin)
%BSQR_RAND_PARSE_OPTIONS Validate bsqr_rand name-value options.
%
% This mirrors the style of matlab/private/bsqr_parse_options.m but adds the
% knobs specific to the randomized variant. The deterministic implementation
% is intentionally left untouched; nothing here is shared with it.

if ~(isnumeric(A) && ismatrix(A) && isreal(A))
    error('bsqr_rand:InvalidInput', 'A must be a real numeric matrix.');
end
if issparse(A)
    error('bsqr_rand:InvalidInput', 'A must be dense (sparse inputs are not supported).');
end

[m, n] = size(A);
kmax = min(m, n);

parser = inputParser;
parser.FunctionName = 'bsqr_rand';
addParameter(parser, 'k', [], @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(parser, 'block_size', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1));
addParameter(parser, 'threshold_mode', 'running_mean', @(x) ischar(x) || isstring(x));
addParameter(parser, 'slack', 1.0, @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 1);
addParameter(parser, 'norm_recomp_tol', sqrt(eps('double')), @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(parser, 'sampling', 'normweighted', @(x) ischar(x) || isstring(x));
addParameter(parser, 'pick', 'best_in_block', @(x) ischar(x) || isstring(x));
addParameter(parser, 'batched', true, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
addParameter(parser, 'seed', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0));
addParameter(parser, 'return_r12', false, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
addParameter(parser, 'backend', 'auto', @(x) ischar(x) || isstring(x));
addParameter(parser, 'check_finite', false, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
parse(parser, varargin{:});

opts = parser.Results;

if isempty(opts.k)
    opts.k = kmax;
else
    opts.k = double(opts.k);
    if abs(opts.k - round(opts.k)) > 0
        error('bsqr_rand:InvalidK', 'k must be an integer in [0, min(m,n)].');
    end
    opts.k = round(opts.k);
end
if opts.k < 0 || opts.k > kmax
    error('bsqr_rand:InvalidK', 'k must satisfy 0 <= k <= min(size(A)).');
end

opts.batched = logical(opts.batched);   % default kernel path: in-block BSQR
if isempty(opts.block_size)
    % Auto block: batched favours block = k (amortizes the per-block reflector
    % apply over many in-block selections); single-select favours the smaller
    % ceil(k/2)-clamped block (one selection per block, less wasted sampling).
    if opts.batched
        opts.block_size = max(opts.k, 1);
    else
        opts.block_size = rand_default_block(opts.k);   % ceil(k/2) clamped to [16,64]
    end
else
    opts.block_size = max(1, round(double(opts.block_size)));
end
opts.slack = double(opts.slack);

% Mirrors the deterministic kernel's safeguard knob (matlab/private/bsqr_parse_options.m).
% Only the MEX in-block path downdates running norms and uses it; the m-file kernel
% recomputes rho^2 / ||w||^2 exactly each in-block step, so the value is a no-op there.
opts.norm_recomp_tol = double(opts.norm_recomp_tol);
if opts.norm_recomp_tol < 0 || opts.norm_recomp_tol > 1
    error('bsqr_rand:InvalidNormRecompTol', 'norm_recomp_tol must satisfy 0 <= value <= 1.');
end

opts.threshold_mode = char(lower(string(opts.threshold_mode)));
if ~ismember(opts.threshold_mode, {'running_mean', 'worstcase_allowance'})
    error('bsqr_rand:InvalidThresholdMode', ...
        'threshold_mode must be "running_mean" or "worstcase_allowance".');
end

opts.sampling = char(lower(string(opts.sampling)));
if ~ismember(opts.sampling, {'uniform', 'normweighted'})
    error('bsqr_rand:InvalidSampling', 'sampling must be "uniform" or "normweighted".');
end

opts.pick = char(lower(string(opts.pick)));
if ~ismember(opts.pick, {'best_in_block', 'first'})
    error('bsqr_rand:InvalidPick', 'pick must be "best_in_block" or "first".');
end

if ~isempty(opts.seed)
    opts.seed = double(opts.seed);
end

opts.return_r12 = logical(opts.return_r12);
opts.check_finite = logical(opts.check_finite);

opts.backend = char(lower(string(opts.backend)));
if ~ismember(opts.backend, {'auto', 'mfile', 'mex'})
    error('bsqr_rand:InvalidBackend', 'backend must be "auto", "mfile", or "mex".');
end

if opts.check_finite && ~all(isfinite(A), 'all')
    error('bsqr_rand:NonFiniteInput', 'A contains non-finite values.');
end

end
