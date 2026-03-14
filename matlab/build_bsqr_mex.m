function build_bsqr_mex(varargin)
%BUILD_BSQR_MEX Build optional MEX backend for bsqr.
%
% v1 ships with m-file BSQR as default. Add C/C++ sources under
% matlab/mex/src and run this builder to produce bsqr_mex.

thisdir = fileparts(mfilename('fullpath'));
srcdir = fullfile(thisdir, 'mex', 'src');
outdir = fullfile(thisdir, 'mex');

if ~isfolder(srcdir)
    error('build_bsqr_mex:MissingSrcDir', 'Missing source directory: %s', srcdir);
end

cc = dir(fullfile(srcdir, '*.c'));
cpp = dir(fullfile(srcdir, '*.cpp'));
sources = [cc; cpp];
if isempty(sources)
    fprintf(2, ['No C/C++ sources found in %s.\n', ...
        'MEX backend is optional in v1; continuing with m-file backend.\n'], srcdir);
    return;
end

src_paths = fullfile({sources.folder}, {sources.name});
mex('-R2018a', '-outdir', outdir, '-output', 'bsqr_mex', src_paths{:}, varargin{:});
fprintf('Built MEX backend at %s\n', fullfile(outdir, ['bsqr_mex.', mexext]));
end
