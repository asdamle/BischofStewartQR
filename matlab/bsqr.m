function varargout = bsqr(A, varargin)
%BSQR Bischof-Stewart pivoted QR in MATLAB.
%   R = BSQR(A)
%   [Q,R] = BSQR(A)
%   [Q,R,E_OR_P] = BSQR(A)
%   [Q,R,E_OR_P,R11INV_R12] = BSQR(A,'return_rinv_r12',true)
%
% Name-value options:
%   'k'                - number of BSQR steps (default min(m,n))
%   'return_rinv_r12'  - logical flag to return R11^{-1}R12 (default false)
%   'pivot_format'     - 'matrix' (default) or 'vector'
%   'backend'          - 'auto' (default), 'mfile', or 'mex'
%   'trace'            - logical flag enabling the per-step validation
%                        trace as a fifth output: a struct with fields
%                        'crit' (criterion value of the selected pivot)
%                        and 'nrecomp' (norm-recompute events per step);
%                        see docs/VALIDATION.md (V4)

if nargout > 5
    error('bsqr:TooManyOutputs', 'bsqr supports at most 5 outputs.');
end

opts = bsqr_parse_options(A, varargin{:});

if nargout > 4 && ~opts.trace
    error('bsqr:TraceNotRequested', ...
        'A fifth output requires the option trace=true.');
end

if should_use_mex(opts.backend)
    ensure_bsqr_mex_ready();
    [varargout{1:nargout}] = bsqr_mex(A, varargin{:});
    return;
end

[Q, R, p, rinv_r12, trace] = bsqr_mfile(A, opts);

if nargout == 1
    varargout{1} = R;
    return;
end

if nargout >= 2
    varargout{1} = Q;
    varargout{2} = R;
end
if nargout >= 3
    varargout{3} = format_pivot_output(p, opts.pivot_format, size(A, 2), class(R));
end
if nargout >= 4
    if opts.return_rinv_r12
        varargout{4} = rinv_r12;
    else
        varargout{4} = [];
    end
end
if nargout >= 5
    varargout{5} = trace;
end

end

function tf = should_use_mex(backend)
switch backend
    case "auto"
        tf = bsqr_mex_available();
    case "mfile"
        tf = false;
    case "mex"
        tf = true;
    otherwise
        error('bsqr:InvalidBackend', 'backend must be "auto", "mfile", or "mex".');
end
end

function ensure_bsqr_mex_ready()
persistent mex_ready
if ~isempty(mex_ready) && mex_ready
    return;
end

thisdir = fileparts(mfilename('fullpath'));
mexdir = fullfile(thisdir, 'mex');
if isfolder(mexdir)
    addpath(mexdir);
end

if ~bsqr_mex_available()
    error('bsqr:MexUnavailable', ...
        'backend="mex" requested but bsqr_mex is not available. Run matlab/build_bsqr_mex.m first.');
end
mex_ready = true;
end

function piv = format_pivot_output(p, pivot_format, n, out_class)
if pivot_format == "vector"
    piv = p;
    return;
end

piv = eye(n, out_class);
piv = piv(:, p);
end
