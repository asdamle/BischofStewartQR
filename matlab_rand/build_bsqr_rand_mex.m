function build_bsqr_rand_mex(varargin)
%BUILD_BSQR_RAND_MEX Build the MEX backend for bsqr_rand.
%
% Sources live under matlab_rand/mex/src. This mirrors matlab/build_bsqr_mex.m
% but targets the randomized variant and writes bsqr_rand_mex.

thisdir = fileparts(mfilename('fullpath'));
srcdir = fullfile(thisdir, 'mex', 'src');
outdir = fullfile(thisdir, 'mex');

if ~isfolder(srcdir)
    error('build_bsqr_rand_mex:MissingSrcDir', 'Missing source directory: %s', srcdir);
end

cc = dir(fullfile(srcdir, '*.c'));
cpp = dir(fullfile(srcdir, '*.cpp'));
sources = [cc; cpp];
if isempty(sources)
    fprintf(2, ['No C/C++ sources found in %s.\n', ...
        'Cannot build bsqr_rand_mex without source files.\n'], srcdir);
    return;
end

src_paths = fullfile({sources.folder}, {sources.name});
mex('-R2018a', 'CXXOPTIMFLAGS=$CXXOPTIMFLAGS -O3 -DNDEBUG', ...
    '-outdir', outdir, '-output', 'bsqr_rand_mex', ...
    src_paths{:}, '-lmwblas', '-lmwlapack', varargin{:});
fprintf('Built MEX backend at %s\n', fullfile(outdir, ['bsqr_rand_mex.', mexext]));
end
