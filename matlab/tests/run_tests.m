function run_tests()
%RUN_TESTS Run MATLAB BSQR test suite.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab'));
addpath(fullfile(repo_root, 'matlab', 'benchmark'));
results = runtests(fullfile(repo_root, 'matlab', 'tests'));
assertSuccess(results);
end
