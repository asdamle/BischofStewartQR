function p = bsqr_bench_canonical_path(inpath)
%BSQR_BENCH_CANONICAL_PATH Resolve canonical path for robust prefix checks.

if isempty(inpath)
    p = '';
    return;
end

if isfolder(inpath)
    p = char(java.io.File(inpath).getCanonicalPath());
    return;
end

[parent, name, ext] = fileparts(inpath);
if isempty(parent)
    parent = pwd;
end
if isfolder(parent)
    p = fullfile(char(java.io.File(parent).getCanonicalPath()), [name, ext]);
else
    p = fullfile(parent, [name, ext]);
end
end
