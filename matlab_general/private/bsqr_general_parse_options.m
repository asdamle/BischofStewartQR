function [opts, fwd] = bsqr_general_parse_options(A, varargin)
%BSQR_GENERAL_PARSE_OPTIONS Validate bsqr_general options; split off forwards.
%   Mirrors the style of matlab/private/bsqr_parse_options.m for the options
%   the wrapper owns. Every unmatched name-value pair is returned in FWD and
%   forwarded verbatim to the selected selector, whose own strict parser
%   rejects unknown names. The selectors' output-shape options
%   ('return_rinv_r12', 'return_r12', 'trace') are blocked here: the wrapper
%   owns its output contract.

if ~(isnumeric(A) && ismatrix(A) && isreal(A))
    error('bsqr_general:InvalidInput', 'A must be a real numeric matrix.');
end
if issparse(A)
    error('bsqr_general:InvalidInput', 'A must be dense (sparse inputs are not supported).');
end

[m, n] = size(A);
if m > n
    error('bsqr_general:NotShortWide', ...
        'A must be short-wide (m <= n); for tall matrices use bsqr directly.');
end

parser = inputParser;
parser.FunctionName = 'bsqr_general';
parser.KeepUnmatched = true;
parser.PartialMatching = false;   % forwarded selector options must never be
                                  % swallowed by a prefix match on our names
addParameter(parser, 'selector', 'bsqr', @(x) ischar(x) || isstring(x));
addParameter(parser, 'k', [], @(x) isnumeric(x) && isscalar(x) && isfinite(x));
addParameter(parser, 'pivot_format', 'matrix', @(x) ischar(x) || isstring(x));
addParameter(parser, 'check_finite', true, @(x) (islogical(x) || isnumeric(x)) && isscalar(x));
addParameter(parser, 'rank_tol', max(m, n) * eps('double'), ...
    @(x) isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0);
parse(parser, varargin{:});

opts = parser.Results;

opts.selector = lower(string(opts.selector));
if opts.selector ~= "bsqr" && opts.selector ~= "bsqr_rand"
    error('bsqr_general:InvalidSelector', 'selector must be "bsqr" or "bsqr_rand".');
end

if isempty(opts.k)
    opts.k = m;
else
    opts.k = double(opts.k);
    if abs(opts.k - round(opts.k)) > 0
        error('bsqr_general:InvalidK', 'k must be an integer in [0, size(A,1)].');
    end
    opts.k = round(opts.k);
end
if opts.k < 0 || opts.k > m
    error('bsqr_general:InvalidK', 'k must satisfy 0 <= k <= size(A,1).');
end

opts.pivot_format = lower(string(opts.pivot_format));
if opts.pivot_format ~= "matrix" && opts.pivot_format ~= "vector"
    error('bsqr_general:InvalidPivotFormat', 'pivot_format must be "matrix" or "vector".');
end

opts.check_finite = logical(opts.check_finite);
opts.rank_tol = double(opts.rank_tol);

if opts.check_finite && ~all(isfinite(A), 'all')
    error('bsqr_general:NonFiniteInput', 'A contains non-finite values.');
end

% Split the unmatched pairs into the forwarded list, blocking the selectors'
% output-shape options.
blocked = ["return_rinv_r12", "return_r12", "trace"];
un = parser.Unmatched;
names = fieldnames(un);
fwd = cell(1, 2 * numel(names));
for i = 1:numel(names)
    if any(lower(string(names{i})) == blocked)
        error('bsqr_general:UnsupportedOption', ...
            'Option "%s" is not supported by bsqr_general (the wrapper owns its output contract).', ...
            names{i});
    end
    fwd{2*i - 1} = names{i};
    fwd{2*i} = un.(names{i});
end

end
