function startup(varargin)
%STARTUP  Set up the BSQR repository for interactive MATLAB use.
%   Run this once from the repository root (MATLAB also runs it automatically
%   when it is launched from this folder). It puts the deterministic and
%   randomized Bischof-Stewart pivoted QR kernels on the path and builds their
%   MEX backends if needed, so the default calls just work:
%
%       R = bsqr(A)                            % deterministic BS pivoted QR
%       [Q, R, E] = bsqr(A)                    % E: permutation matrix (see matlab/README.md)
%       [p, Q, R11] = bsqr_rand(A)             % randomized variant
%
%   Both dispatch to the compiled MEX backend by default; pass
%   'backend','mfile' for the pure-MATLAB reference implementation.
%
%   STARTUP('rebuild') forces the MEX backends to be recompiled.
%
%   See matlab/README.md and matlab_rand/README.md for the full APIs, and the
%   Julia package under julia/ for the deterministic kernel in Julia.

root = fileparts(mfilename('fullpath'));
addpath(fullfile(root, 'matlab'),      fullfile(root, 'matlab', 'mex'), ...
        fullfile(root, 'matlab_rand'), fullfile(root, 'matlab_rand', 'mex'));

force = any(strcmpi(varargin, 'rebuild'));
ok_det  = build_if_needed('bsqr',      @bsqr_mex_available,      @build_bsqr_mex,      force);
ok_rand = build_if_needed('bsqr_rand', @bsqr_rand_mex_available, @build_bsqr_rand_mex, force);
clear bsqr_mex_available bsqr_rand_mex_available   % drop cached availability so the dispatchers re-check

backend = 'MEX';
if ~(ok_det && ok_rand); backend = 'pure-MATLAB fallback for the ones that failed'; end
fprintf(['\nBSQR ready (%s).\n', ...
    '  R = bsqr(A)              [Q,R,E] = bsqr(A)            deterministic BS pivoted QR\n', ...
    '  [p,Q,R11] = bsqr_rand(A)                             randomized variant\n', ...
    'Pass ''backend'',''mfile'' for the pure-MATLAB reference; see the matlab/ and\n', ...
    'matlab_rand/ READMEs.\n'], backend);
end

function ok = build_if_needed(name, avail_fn, build_fn, force)
if ~force && avail_fn()
    fprintf('  %-9s MEX backend: ready\n', name);
    ok = true;
    return;
end
try
    fprintf('  %-9s MEX backend: %s...\n', name, ternary(force, 'rebuilding', 'building'));
    build_fn();
    ok = true;
catch ME
    ok = false;
    warning('startup:MexBuildFailed', ...
        ['Could not build %s_mex (%s).\n%s will fall back to the slower pure-MATLAB ', ...
         'backend.\nInstall a supported C++ compiler (run "mex -setup C++") and rerun ', ...
         '"startup rebuild" to enable the MEX backend.'], name, ME.message, name);
end
end

function out = ternary(cond, a, b)
if cond; out = a; else; out = b; end
end
