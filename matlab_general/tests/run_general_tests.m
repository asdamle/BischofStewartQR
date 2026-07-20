function run_general_tests()
%RUN_GENERAL_TESTS Run the general-A wrapper test suite.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab'));
addpath(fullfile(repo_root, 'matlab', 'mex'));
addpath(fullfile(repo_root, 'matlab_rand'));
addpath(fullfile(repo_root, 'matlab_rand', 'mex'));
addpath(fullfile(repo_root, 'matlab_rand', 'benchmark'));   % rand_test_matrix
addpath(fullfile(repo_root, 'matlab_general'));
addpath(fullfile(repo_root, 'matlab_general', 'benchmark'));
results = runtests(fullfile(repo_root, 'matlab_general', 'tests'));
assertSuccess(results);
end
