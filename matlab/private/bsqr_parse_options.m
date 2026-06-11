function opts = bsqr_parse_options(A, varargin)
%BSQR_PARSE_OPTIONS Validate bsqr name-value options.

if ~(isnumeric(A) && ismatrix(A) && isreal(A))
    error('bsqr:InvalidInput', 'A must be a real numeric matrix.');
end

[m, n] = size(A);
kmax = min(m, n);

parser = inputParser;
parser.FunctionName = 'bsqr';
addParameter(parser, 'k', [], @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(parser, 'return_rinv_r12', false, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
addParameter(parser, 'pivot_format', 'matrix', @(x) ischar(x) || isstring(x));
addParameter(parser, 'backend', 'auto', @(x) ischar(x) || isstring(x));
addParameter(parser, 'norm_recomp_tol', sqrt(eps('double')), @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(parser, 'check_finite', true, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
addParameter(parser, 'trace', false, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
parse(parser, varargin{:});

opts = parser.Results;

if isempty(opts.k)
    opts.k = kmax;
else
    opts.k = double(opts.k);
    if abs(opts.k - round(opts.k)) > 0
        error('bsqr:InvalidK', 'k must be an integer in [0, min(m,n)].');
    end
    opts.k = round(opts.k);
end

if opts.k < 0 || opts.k > kmax
    error('bsqr:InvalidK', 'k must satisfy 0 <= k <= min(size(A)).');
end

opts.return_rinv_r12 = logical(opts.return_rinv_r12);
opts.check_finite = logical(opts.check_finite);
opts.trace = logical(opts.trace);

opts.pivot_format = lower(string(opts.pivot_format));
if opts.pivot_format ~= "matrix" && opts.pivot_format ~= "vector"
    error('bsqr:InvalidPivotFormat', 'pivot_format must be "matrix" or "vector".');
end

opts.backend = lower(string(opts.backend));
if opts.backend ~= "auto" && opts.backend ~= "mfile" && opts.backend ~= "mex"
    error('bsqr:InvalidBackend', 'backend must be "auto", "mfile", or "mex".');
end

opts.norm_recomp_tol = double(opts.norm_recomp_tol);
if opts.norm_recomp_tol < 0 || opts.norm_recomp_tol > 1
    error('bsqr:InvalidNormRecompTol', 'norm_recomp_tol must satisfy 0 <= value <= 1.');
end

if opts.check_finite && ~all(isfinite(A), 'all')
    error('bsqr:NonFiniteInput', 'A contains non-finite values.');
end

end
