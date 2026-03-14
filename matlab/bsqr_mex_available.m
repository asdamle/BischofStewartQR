function tf = bsqr_mex_available()
%BSQR_MEX_AVAILABLE True when compiled bsqr_mex backend is available.

persistent cached_tf cached_mexext
if ~isempty(cached_tf) && isequal(cached_mexext, mexext)
    tf = cached_tf;
    return;
end

thisdir = fileparts(mfilename('fullpath'));
mex_candidate = fullfile(thisdir, 'mex', ['bsqr_mex.', mexext]);
tf = isfile(mex_candidate) || ~isempty(which('bsqr_mex'));
cached_tf = tf;
cached_mexext = mexext;
end
