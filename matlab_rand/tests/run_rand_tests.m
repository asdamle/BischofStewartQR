function run_rand_tests()
%RUN_RAND_TESTS Run the randomized BSQR test suite.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab_rand'));
addpath(fullfile(repo_root, 'matlab_rand', 'benchmark'));
results = runtests(fullfile(repo_root, 'matlab_rand', 'tests'));
assertSuccess(results);
end
