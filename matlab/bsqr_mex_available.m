function tf = bsqr_mex_available()
%BSQR_MEX_AVAILABLE True when compiled bsqr_mex backend is available.

thisdir = fileparts(mfilename('fullpath'));
mex_candidate = fullfile(thisdir, 'mex', ['bsqr_mex.', mexext]);
tf = isfile(mex_candidate) || ~isempty(which('bsqr_mex'));
end
