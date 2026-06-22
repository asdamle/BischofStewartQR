function [A, info] = approx_real_matrix(name, datadir)
%APPROX_REAL_MATRIX Load a real (non-synthetic) matrix for the approximation
%   comparison from the git-ignored ext_comparisons/data/ folder. Real data is
%   never committed to this repo; see matlab_rand/README.md for where to get it.
%
%   [A, info] = approx_real_matrix(NAME) loads ext_comparisons/data/<NAME>.mat.
%   The file should contain a variable named A (a 2-D double matrix); if not, the
%   first 2-D numeric variable in the file is used. Sparse matrices are densified
%   (the approximation harness needs a dense A for svds/projection on moderate
%   sizes). Errors with concise download instructions if the file is absent.
%
%   approx_real_matrix(NAME, DATADIR) overrides the default data directory.

if nargin < 2 || isempty(datadir)
    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    datadir = fullfile(repo_root, 'ext_comparisons', 'data');
end
f = fullfile(datadir, [char(name) '.mat']);

if ~isfile(f)
    error('approx_real_matrix:Missing', ['Real matrix "%s" not found at\n  %s\n', ...
        'Real data is not committed. Place a .mat holding a 2-D double variable A\n', ...
        'in ext_comparisons/data/ (git-ignored). See matlab_rand/README.md for\n', ...
        'suggested sources (e.g. a matrix from https://sparse.tamu.edu/).'], ...
        char(name), f);
end

S = load(f);
A = pick_matrix(S, char(name), f);
if issparse(A); A = full(A); end
A = double(A);
info = struct('family', char(name), 'm', size(A, 1), 'n', size(A, 2), ...
    'desc', sprintf('real matrix from %s', f));
end

% ---------------------------------------------------------------------------
function A = pick_matrix(S, name, f)
% Prefer a variable literally named A; otherwise take the first 2-D numeric one.
if isfield(S, 'A') && isnumeric(S.A) && ismatrix(S.A) && min(size(S.A)) > 0
    A = S.A; return;
end
fn = fieldnames(S);
for i = 1:numel(fn)
    v = S.(fn{i});
    if isnumeric(v) && ismatrix(v) && min(size(v)) > 1
        A = v; return;
    end
end
error('approx_real_matrix:NoMatrix', ...
    'No 2-D numeric matrix found in %s (for "%s").', f, name);
end
