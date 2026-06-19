function run_rand_experiments(varargin)
%RUN_RAND_EXPERIMENTS Generate the data behind the publication plots.
%
%   Writes three tidy CSVs under matlab_rand/benchmark/results/:
%     exp_scaling.csv   - time / conditioning / samples vs n, per family,
%                         deterministic vs randomized (both threshold modes).
%     exp_blocksize.csv - effect of block_size on time / samples / conditioning,
%                         uniform vs norm-weighted sampling.
%     exp_sampling.csv   - uniform vs norm-weighted sampling across families.
%
%   All randomized timing uses the MEX kernel called directly with
%   check_finite=false (see run_rand_benchmarks for the fairness rationale).
%   Each configuration is repeated over several seeds so the plotter can draw
%   confidence bands. Plot with plot_rand_experiments.
%
% Options: 'trials' (seeds, default 5), 'which' (cellstr subset of
%   {'scaling','blocksize','sampling'}, default all), 'outdir'.

ip = inputParser;
addParameter(ip, 'trials', 5);
addParameter(ip, 'which', {'scaling', 'blocksize', 'sampling'});
addParameter(ip, 'outdir', '');
parse(ip, varargin{:});
opt = ip.Results;
if ischar(opt.which) || isstring(opt.which); opt.which = cellstr(opt.which); end

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'matlab_rand'));
addpath(fullfile(repo_root, 'matlab_rand', 'mex'));
addpath(fullfile(repo_root, 'matlab_rand', 'benchmark'));
addpath(fullfile(repo_root, 'matlab'));
addpath(fullfile(repo_root, 'matlab', 'mex'));
if isempty(opt.outdir)
    opt.outdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'results');
end
if ~isfolder(opt.outdir); mkdir(opt.outdir); end
assert(exist('bsqr_rand_mex', 'file') == 3 && exist('bsqr_mex', 'file') == 3, ...
    'Both MEX kernels must be built first.');

fprintf('BLAS threads = %d, trials = %d\n', maxNumCompThreads, opt.trials);

if any(strcmp(opt.which, 'scaling'));   exp_scaling(opt);   end
if any(strcmp(opt.which, 'blocksize')); exp_blocksize(opt); end
if any(strcmp(opt.which, 'sampling'));  exp_sampling(opt);  end
end

% ===========================================================================
function exp_scaling(opt)
families = {'gaussian', 'spiked_leverage', 'needle'};
k = 64;
ns = [1000, 2000, 4000, 8000, 16000, 32000, 64000];
modes = {'running_mean', 'worstcase_allowance'};
rows = {};
for fi = 1:numel(families)
    fam = families{fi};
    for n = ns
        for s = 1:opt.trials
            M = rand_test_matrix(fam, k, n, 1000 * fi + s);
            d = measure_det(M, k);
            rows(end+1, :) = {fam, k, n, s, 'det', 'na', 'na', 0, d.time, d.frobinv, d.osinsky, NaN}; %#ok<AGROW>
            for mi = 1:numel(modes)
                r = measure_rand(M, k, modes{mi}, 'uniform', 16, 1000 * fi + s);
                rows(end+1, :) = {fam, k, n, s, 'rand', modes{mi}, 'uniform', 16, ...
                    r.time, r.frobinv, r.osinsky, r.tested_per_k}; %#ok<AGROW>
            end
        end
        fprintf('  scaling %-16s n=%-6d done\n', fam, n);
    end
end
write_rows(rows, fullfile(opt.outdir, 'exp_scaling.csv'));
end

% ===========================================================================
function exp_blocksize(opt)
families = {'gaussian', 'spiked_leverage', 'needle'};
k = 64; n = 8000;
blocks = [1, 2, 4, 8, 16, 32, 64, 128];
samplings = {'uniform', 'normweighted'};
rows = {};
for fi = 1:numel(families)
    fam = families{fi};
    for s = 1:opt.trials
        M = rand_test_matrix(fam, k, n, 2000 * fi + s);
        for bi = 1:numel(blocks)
            for si = 1:numel(samplings)
                r = measure_rand(M, k, 'running_mean', samplings{si}, blocks(bi), 2000 * fi + s);
                rows(end+1, :) = {fam, k, n, s, 'rand', 'running_mean', samplings{si}, blocks(bi), ...
                    r.time, r.frobinv, r.osinsky, r.tested_per_k}; %#ok<AGROW>
            end
        end
    end
    fprintf('  blocksize %-16s done\n', fam);
end
write_rows(rows, fullfile(opt.outdir, 'exp_blocksize.csv'));
end

% ===========================================================================
function exp_sampling(opt)
families = {'gaussian', 'graded_leverage', 'spiked_leverage', 'coherent', 'needle'};
k = 64; n = 32000; block = 16;
samplings = {'uniform', 'normweighted'};
rows = {};
for fi = 1:numel(families)
    fam = families{fi};
    for s = 1:opt.trials
        M = rand_test_matrix(fam, k, n, 3000 * fi + s);
        for si = 1:numel(samplings)
            r = measure_rand(M, k, 'running_mean', samplings{si}, block, 3000 * fi + s);
            rows(end+1, :) = {fam, k, n, s, 'rand', 'running_mean', samplings{si}, block, ...
                r.time, r.frobinv, r.osinsky, r.tested_per_k}; %#ok<AGROW>
        end
    end
    fprintf('  sampling %-16s done\n', fam);
end
write_rows(rows, fullfile(opt.outdir, 'exp_sampling.csv'));
end

% ===========================================================================
function out = measure_rand(M, k, mode, sampling, block, seed)
f = @() bsqr_rand_mex(M, 'k', k, 'check_finite', false, 'block_size', block, ...
    'threshold_mode', mode, 'sampling', sampling, 'seed', seed);
out.time = timeit(f, 3);
[~, ~, ~, st] = f();
out.frobinv = st.frob_inv;
out.osinsky = st.osinsky_bound;
out.tested_per_k = st.total_tested / k;
end

function out = measure_det(M, k)
out.time = timeit(@() bsqr_mex(M, 'k', k, 'check_finite', false), 1);
R = bsqr_mex(M, 'k', k, 'check_finite', false);
out.frobinv = norm(inv(triu(R(1:k, 1:k))), 'fro');
out.osinsky = sqrt(k * (size(M, 2) - k + 1));
end

function write_rows(rows, path)
T = cell2table(rows, 'VariableNames', {'family', 'k', 'n', 'seed', 'method', ...
    'mode', 'sampling', 'block_size', 'time_s', 'frobinv', 'osinsky', 'tested_per_k'});
writetable(T, path);
fprintf('Wrote %s (%d rows)\n', path, height(T));
end
