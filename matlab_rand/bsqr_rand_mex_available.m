function tf = bsqr_rand_mex_available()
%BSQR_RAND_MEX_AVAILABLE True when the compiled bsqr_rand_mex backend exists.

persistent cached_tf cached_mexext
if ~isempty(cached_tf) && isequal(cached_mexext, mexext)
    tf = cached_tf;
    return;
end

thisdir = fileparts(mfilename('fullpath'));
mex_candidate = fullfile(thisdir, 'mex', ['bsqr_rand_mex.', mexext]);
tf = isfile(mex_candidate) || ~isempty(which('bsqr_rand_mex'));
cached_tf = tf;
cached_mexext = mexext;
end
