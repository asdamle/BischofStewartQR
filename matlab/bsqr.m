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

if nargout > 4
    error('bsqr:TooManyOutputs', 'bsqr supports at most 4 outputs.');
end

opts = bsqr_parse_options(A, varargin{:});

use_mex = false;
switch opts.backend
    case "auto"
        use_mex = bsqr_mex_available();
    case "mfile"
        use_mex = false;
    case "mex"
        use_mex = true;
end

if use_mex
    if ~bsqr_mex_available()
        error('bsqr:MexUnavailable', ...
            'backend="mex" requested but bsqr_mex is not available. Run matlab/build_bsqr_mex.m first.');
    end
    thisdir = fileparts(mfilename('fullpath'));
    mexdir = fullfile(thisdir, 'mex');
    if isfolder(mexdir)
        addpath(mexdir);
    end
    [varargout{1:nargout}] = bsqr_mex(A, varargin{:});
    return;
end

[Q, R, p, rinv_r12] = bsqr_mfile(A, opts);

if nargout == 1
    varargout{1} = R;
    return;
end

if nargout >= 2
    varargout{1} = Q;
    varargout{2} = R;
end
if nargout >= 3
    if opts.pivot_format == "vector"
        piv = p;
    else
        n = size(A, 2);
        piv = eye(n, class(R));
        piv = piv(:, p);
    end
    varargout{3} = piv;
end
if nargout >= 4
    if opts.return_rinv_r12
        varargout{4} = rinv_r12;
    else
        varargout{4} = [];
    end
end

end
