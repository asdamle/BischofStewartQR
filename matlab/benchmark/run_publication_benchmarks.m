function run_publication_benchmarks(varargin)
%RUN_PUBLICATION_BENCHMARKS MATLAB publication benchmark workflow for BSQR.
%
% Example:
%   addpath('matlab')
%   run_publication_benchmarks

addpath(fileparts(fileparts(mfilename('fullpath'))));

cfg = parse_config(varargin{:});
ensure_mex_path_ready();
if ~isfolder(cfg.outdir)
    mkdir(cfg.outdir);
end

run_id = datestr(now, 'yyyymmdd_HHMMSS');
timestamp = string(datetime('now'));

cases = make_case_grid(cfg.square_ms, cfg.short_ms, cfg.short_aspects);
rows = table();

for seed = cfg.seeds
    for f = 1:numel(cfg.families)
        family = cfg.families{f};
        for c = 1:numel(cases)
            kase = cases(c);
            if strcmp(family, 'orthonormal_rows') && kase.m > kase.n
                continue;
            end

            case_seed = compute_case_seed(seed, family, kase.m, kase.n, kase.aspect);
            rng(case_seed, 'twister');
            A = make_matrix(family, kase.m, kase.n);
            k = min(kase.m, kase.n);

            plain = bench_pair(A, k, cfg.norm_recomp_tol, cfg.warmup, cfg.samples, false);
            rinv = bench_pair(A, k, cfg.norm_recomp_tol, cfg.warmup, cfg.samples, true);

            rows = [rows; pack_rows(plain, run_id, timestamp, family, kase, seed)]; %#ok<AGROW>
            rows = [rows; pack_rows(rinv, run_id, timestamp, family, kase, seed)]; %#ok<AGROW>
        end
    end
end

csv_path = fullfile(cfg.outdir, 'publication_timings.csv');
summary_path = fullfile(cfg.outdir, 'publication_summary.md');
metadata_path = fullfile(cfg.outdir, 'metadata.txt');

writetable(rows, csv_path);
write_summary(rows, summary_path, run_id, cfg);
write_metadata(metadata_path, run_id, cfg, rows);

plots_dir = fullfile(cfg.outdir, 'plots');
tables_dir = fullfile(cfg.outdir, 'tables');
plot_publication_results(csv_path, plots_dir, tables_dir);

fprintf('Wrote publication benchmark CSV to: %s\n', csv_path);
fprintf('Wrote publication summary to: %s\n', summary_path);
fprintf('Wrote publication metadata to: %s\n', metadata_path);
end

function cfg = parse_config(varargin)
if nargin == 1 && isstruct(varargin{1})
    cfg = varargin{1};
else
    cfg = struct();
end
has_backend_field = isfield(cfg, 'bsqr_backend');
if has_backend_field
    requested_backend = lower(string(cfg.bsqr_backend));
else
    requested_backend = "mex";
end

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
def_outdir = fullfile(repo_root, 'matlab', 'benchmark', 'results', 'publication');

cfg.outdir = getfield_default(cfg, 'outdir', bsqr_bench_getenv_default('BS_MATLAB_PUB_OUTDIR', def_outdir));
cfg.allow_shared_outdir = getfield_default(cfg, 'allow_shared_outdir', ...
    bsqr_bench_parse_bool(bsqr_bench_getenv_default('BS_MATLAB_ALLOW_SHARED_OUTDIR', '0'), ...
    'run_publication_benchmarks:InvalidBool'));
cfg.seeds = getfield_default(cfg, 'seeds', parse_int_list(bsqr_bench_getenv_default('BS_MATLAB_PUB_SEEDS', '20260310,20260311')));
cfg.families = getfield_default(cfg, 'families', parse_str_list(bsqr_bench_getenv_default('BS_MATLAB_PUB_FAMILIES', 'gaussian,ill_conditioned,orthonormal_rows')));
cfg.square_ms = getfield_default(cfg, 'square_ms', parse_int_list(bsqr_bench_getenv_default('BS_MATLAB_PUB_SQUARE_MS', '64,128,256,384,512')));
cfg.short_ms = getfield_default(cfg, 'short_ms', parse_int_list(bsqr_bench_getenv_default('BS_MATLAB_PUB_SHORT_MS', '32,64,128,256,512')));
cfg.short_aspects = getfield_default(cfg, 'short_aspects', parse_float_list(bsqr_bench_getenv_default('BS_MATLAB_PUB_SHORT_ASPECTS', '2,4,8,10')));
cfg.warmup = getfield_default(cfg, 'warmup', str2double(bsqr_bench_getenv_default('BS_MATLAB_PUB_WARMUP', '1')));
cfg.samples = getfield_default(cfg, 'samples', str2double(bsqr_bench_getenv_default('BS_MATLAB_PUB_SAMPLES', '12')));
cfg.norm_recomp_tol = getfield_default(cfg, 'norm_recomp_tol', str2double(bsqr_bench_getenv_default('BS_MATLAB_NORM_RECOMP_TOL', num2str(sqrt(eps('double'))))));

if isempty(cfg.seeds)
    error('run_publication_benchmarks:InvalidConfig', 'At least one seed is required.');
end
if isempty(cfg.families)
    error('run_publication_benchmarks:InvalidConfig', 'At least one family is required.');
end
if has_backend_field
    if requested_backend ~= "mex"
        error('run_publication_benchmarks:MexOnly', ...
            'Publication benchmarks are mex-only.');
    end
    error('run_publication_benchmarks:UnsupportedConfigField', ...
        'Publication benchmarks no longer accept cfg.bsqr_backend. This runner is mex-only.');
end
cfg.bsqr_backend = "mex";

validate_matlab_outdir(cfg.outdir, repo_root, cfg.allow_shared_outdir);
end

function rows = pack_rows(method_rows, run_id, timestamp, family, kase, seed)
rows = table();
for i = 1:numel(method_rows)
    r = method_rows(i);
    newrow = table( ...
        string(run_id), string(timestamp), string(family), string(kase.regime), ...
        kase.m, kase.n, kase.aspect, seed, string(r.method), ...
        r.tmin, r.tmed, r.residual, r.orthogonality, ...
        'VariableNames', {'run_id','timestamp','family','regime','m','n','aspect','seed','method','tmin_s','tmed_s','residual','orthogonality'});
    rows = [rows; newrow]; %#ok<AGROW>
end
end

function rows = bench_pair(A, k, norm_recomp_tol, warmup, samples, include_rinv)
rows = struct('method', {}, 'tmin', {}, 'tmed', {}, 'residual', {}, 'orthogonality', {});

f_bs = @() run_bsqr_timed(A, k, norm_recomp_tol, include_rinv);
[tmin, tmed] = bench_function(f_bs, warmup, samples);
out = run_bsqr_quality(A, k, norm_recomp_tol, include_rinv);
[resid, orth] = quality_from_factor(A, out.Q, out.R, out.p);
if include_rinv
    bs_label = "bsqr_rinv";
else
    bs_label = "bsqr_full";
end
rows(end+1) = struct('method', bs_label, 'tmin', tmin, 'tmed', tmed, 'residual', resid, 'orthogonality', orth); %#ok<AGROW>

f_qr = @() run_qr_builtin_timed(A);
[tmin, tmed] = bench_function(f_qr, warmup, samples);
out = run_qr_builtin_quality(A);
[resid, orth] = quality_from_factor(A, out.Q, out.R, out.p);
rows(end+1) = struct('method', "qr_pivoted", 'tmin', tmin, 'tmed', tmed, 'residual', resid, 'orthogonality', orth); %#ok<AGROW>

if include_rinv
    f_qr_trsm = @() run_qr_trsm_timed(A);
    [tmin, tmed] = bench_function(f_qr_trsm, warmup, samples);
    out = run_qr_trsm_quality(A);
    [resid, orth] = quality_from_factor(A, out.Q, out.R, out.p);
    rows(end+1) = struct('method', "qr_pivoted_trsm", 'tmin', tmin, 'tmed', tmed, 'residual', resid, 'orthogonality', orth); %#ok<AGROW>
end
end

function [tmin, tmed] = bench_function(f, warmup, samples)
for i = 1:warmup
    f();
end

times = zeros(samples, 1);
for i = 1:samples
    tic;
    f();
    times(i) = toc;
end
tmin = min(times);
tmed = median(times);
end

function out = run_bsqr_quality(A, k, norm_recomp_tol, return_rinv)
if return_rinv
    [Q, R, p, rinv] = bsqr_mex(A, 'k', k, 'return_rinv_r12', true, ...
        'pivot_format', 'vector', 'norm_recomp_tol', norm_recomp_tol, ...
        'check_finite', false);
    out = struct('Q', Q, 'R', R, 'p', p, 'rinv', rinv);
else
    [Q, R, p] = bsqr_mex(A, 'k', k, 'return_rinv_r12', false, ...
        'pivot_format', 'vector', 'norm_recomp_tol', norm_recomp_tol, ...
        'check_finite', false);
    out = struct('Q', Q, 'R', R, 'p', p, 'rinv', []);
end
validate_economy_vector_factor(A, out.Q, out.R, out.p, 'bsqr');
end

function run_bsqr_timed(A, k, norm_recomp_tol, return_rinv)
if return_rinv
    [Q, R, p, rinv] = bsqr_mex(A, 'k', k, 'return_rinv_r12', true, ...
        'pivot_format', 'vector', 'norm_recomp_tol', norm_recomp_tol, ...
        'check_finite', false); %#ok<ASGLU>
else
    [Q, R, p] = bsqr_mex(A, 'k', k, 'return_rinv_r12', false, ...
        'pivot_format', 'vector', 'norm_recomp_tol', norm_recomp_tol, ...
        'check_finite', false); %#ok<ASGLU>
end
end

function out = run_qr_builtin_quality(A)
[Q, R, p] = qr_vector(A);
out = struct('Q', Q, 'R', R, 'p', p);
validate_economy_vector_factor(A, out.Q, out.R, out.p, 'qr_pivoted');
end

function run_qr_builtin_timed(A)
[Q, R, p] = qr_vector(A); %#ok<ASGLU>
end

function out = run_qr_trsm_quality(A)
[Q, R, p] = qr_vector(A);
[m, n] = size(A);
k = min(m, n);
if n > k && k > 0
    R11 = R(1:k, 1:k);
    R12 = R(1:k, k+1:n);
    R11 \ R12; %#ok<VUNUS>
end
out = struct('Q', Q, 'R', R, 'p', p);
validate_economy_vector_factor(A, out.Q, out.R, out.p, 'qr_pivoted_trsm');
end

function run_qr_trsm_timed(A)
[Q, R, p] = qr_vector(A); %#ok<ASGLU>
[m, n] = size(A);
k = min(m, n);
if n > k && k > 0
    R11 = R(1:k, 1:k);
    R12 = R(1:k, k+1:n);
    R11 \ R12; %#ok<VUNUS>
end
end

function [Q, R, p] = qr_vector(A)
% Benchmark contract: built-in baseline must use economy QR + vector pivots.
try
    [Q, R, p] = qr(A, 'econ', 'vector');
catch
    % Compatibility fallback for versions that accept numeric economy flag.
    [Q, R, p] = qr(A, 0, 'vector');
end
p = p(:).';
end

function validate_economy_vector_factor(A, Q, R, p, method_name)
[m, n] = size(A);
k = min(m, n);

if ~(isnumeric(p) && isvector(p) && numel(p) == n)
    error('run_publication_benchmarks:NonVectorPermutation', ...
        '%s must return permutation as a length-n vector.', method_name);
end

if ~(isequal(size(Q), [m, k]) && isequal(size(R), [k, n]))
    error('run_publication_benchmarks:NonEconomyQR', ...
        '%s must return economy-sized Q(m,k) and R(k,n).', method_name);
end
end

function [resid, orth] = quality_from_factor(A, Q, R, p)
if isempty(Q)
    resid = norm(A(:, p) - Q * R, 'fro') / max(norm(A, 'fro'), eps('double'));
    orth = 0;
    return;
end
resid = norm(A(:, p) - Q * R, 'fro') / max(norm(A, 'fro'), eps('double'));
orth = norm(eye(size(Q, 2)) - Q' * Q, 'fro');
end

function cases = make_case_grid(square_ms, short_ms, short_aspects)
cases = struct('regime', {}, 'm', {}, 'n', {}, 'aspect', {});
for i = 1:numel(square_ms)
    m = square_ms(i);
    cases(end+1) = struct('regime', 'square', 'm', m, 'n', m, 'aspect', 1.0); %#ok<AGROW>
end
for i = 1:numel(short_ms)
    m = short_ms(i);
    for j = 1:numel(short_aspects)
        aspect = short_aspects(j);
        n = round(m * aspect);
        if n < m
            n = m;
        end
        cases(end+1) = struct('regime', 'short_wide', 'm', m, 'n', n, 'aspect', aspect); %#ok<AGROW>
    end
end
end

function A = make_matrix(family, m, n)
switch family
    case 'gaussian'
        A = randn(m, n);
    case 'ill_conditioned'
        r = min(m, n);
        U = orth(randn(m, r));
        V = orth(randn(n, r));
        s = exp(linspace(0, -log(1e10), r));
        A = U * diag(s) * V';
    case 'orthonormal_rows'
        if m > n
            error('orthonormal_rows requires m <= n.');
        end
        Q = orth(randn(n, m));
        A = Q(:, 1:m)';
    otherwise
        error('Unknown family: %s', family);
end
A = double(A);
end

function seed_out = compute_case_seed(seed, family, m, n, aspect)
family_code = 0;
switch family
    case 'gaussian'
        family_code = 1;
    case 'ill_conditioned'
        family_code = 2;
    case 'orthonormal_rows'
        family_code = 3;
end
aspect_key = round(aspect * 1e6);
seed_out = mod(seed + 1000003 * family_code + 101 * m + 17 * n + aspect_key, 2^31 - 1);
if seed_out <= 0
    seed_out = 1;
end
end

function write_summary(rows, path, run_id, cfg)
plain = relative_time_table(rows, "bsqr_full", "qr_pivoted");
rinv = relative_time_table(rows, "bsqr_rinv", "qr_pivoted_trsm");

fid = fopen(path, 'w');
if fid < 0
    error('Could not open summary path: %s', path);
end
clean = onCleanup(@() fclose(fid));

fprintf(fid, '# MATLAB Publication Benchmark Summary\n\n');
fprintf(fid, '- Run ID: `%s`\n', run_id);
fprintf(fid, '- Generated: %s\n', string(datetime('now')));
fprintf(fid, '- Seeds: %s\n', join(string(cfg.seeds), ', '));
fprintf(fid, '- Families: %s\n\n', join(string(cfg.families), ', '));
fprintf(fid, 'Relative time is BSQR median time divided by baseline median time; `1.0` is parity.\n\n');

fprintf(fid, '| family | regime | geomean relative time (bsqr_full/qr) |\n');
fprintf(fid, '|---|---|---:|\n');
for i = 1:height(plain)
    fprintf(fid, '| %s | %s | %.5f |\n', plain.family(i), plain.regime(i), plain.geomean_relative_time(i));
end
fprintf(fid, '\n');

fprintf(fid, '| family | regime | geomean relative time (bsqr_rinv/qr_trsm) |\n');
fprintf(fid, '|---|---|---:|\n');
for i = 1:height(rinv)
    fprintf(fid, '| %s | %s | %.5f |\n', rinv.family(i), rinv.regime(i), rinv.geomean_relative_time(i));
end
end

function write_metadata(path, run_id, cfg, rows)
fid = fopen(path, 'w');
if fid < 0
    error('Could not open metadata path: %s', path);
end
clean = onCleanup(@() fclose(fid));

fprintf(fid, 'schema_version = 2026-03-16.matlab.v3\n');
fprintf(fid, 'run_id = %s\n', run_id);
fprintf(fid, 'timestamp = %s\n', string(datetime('now')));
fprintf(fid, 'matlab_version = %s\n', version);
fprintf(fid, 'outdir = %s\n', cfg.outdir);
fprintf(fid, 'observed_rows = %d\n', height(rows));
fprintf(fid, 'warmup = %d\n', cfg.warmup);
fprintf(fid, 'samples = %d\n', cfg.samples);
fprintf(fid, 'norm_recomp_tol = %.17g\n', cfg.norm_recomp_tol);
fprintf(fid, 'bsqr_backend = %s\n', cfg.bsqr_backend);
end

function agg = relative_time_table(rows, bs_method, baseline_method)
rt = pair_relative_times(rows, bs_method, baseline_method);
if isempty(rt)
    agg = table(strings(0,1), strings(0,1), zeros(0,1), 'VariableNames', {'family','regime','geomean_relative_time'});
    return;
end

keys = unique(rt(:, {'family', 'regime'}), 'rows');
vals = zeros(height(keys), 1);
for i = 1:height(keys)
    mask = rt.family == keys.family(i) & rt.regime == keys.regime(i);
    vals_i = rt.relative_time(mask);
    vals(i) = exp(mean(log(vals_i(vals_i > 0 & isfinite(vals_i)))));
end
agg = table(keys.family, keys.regime, vals, 'VariableNames', {'family','regime','geomean_relative_time'});
end

function rt = pair_relative_times(rows, bs_method, baseline_method)
rows_bs = rows(rows.method == bs_method, :);
rows_base = rows(rows.method == baseline_method, :);

if isempty(rows_bs) || isempty(rows_base)
    rt = table();
    return;
end

key_cols = {'family', 'regime', 'm', 'n', 'aspect', 'seed'};
join_bs = rows_bs(:, [key_cols, {'tmed_s'}]);
join_bs.Properties.VariableNames{'tmed_s'} = 'tmed_bs';
join_base = rows_base(:, [key_cols, {'tmed_s'}]);
join_base.Properties.VariableNames{'tmed_s'} = 'tmed_base';

joined = innerjoin(join_bs, join_base, 'Keys', key_cols);
joined.relative_time = joined.tmed_bs ./ joined.tmed_base;
rt = joined(:, [key_cols, {'relative_time'}]);
end

function out = parse_int_list(str)
if isempty(str)
    out = [];
    return;
end
parts = strsplit(str, ',');
out = zeros(1, numel(parts));
for i = 1:numel(parts)
    out(i) = str2double(strtrim(parts{i}));
end
out = unique(out(~isnan(out)));
end

function out = parse_float_list(str)
if isempty(str)
    out = [];
    return;
end
parts = strsplit(str, ',');
out = zeros(1, numel(parts));
for i = 1:numel(parts)
    out(i) = str2double(strtrim(parts{i}));
end
out = unique(out(~isnan(out)));
end

function out = parse_str_list(str)
if isempty(str)
    out = {};
    return;
end
parts = strsplit(str, ',');
out = cell(1, numel(parts));
for i = 1:numel(parts)
    out{i} = char(strtrim(parts{i}));
end
out = unique(out);
end

function ensure_mex_path_ready()
persistent mex_ready
if ~isempty(mex_ready) && mex_ready
    return;
end
repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
mexdir = fullfile(repo_root, 'matlab', 'mex');
if isfolder(mexdir)
    addpath(mexdir);
end
if isempty(which('bsqr_mex'))
    error('run_publication_benchmarks:MexUnavailable', ...
        'Publication benchmarks require bsqr_mex, but it is not available.');
end
mex_ready = true;
end

function v = getfield_default(s, field, default)
if isfield(s, field)
    v = s.(field);
else
    v = default;
end
end

function validate_matlab_outdir(outdir, repo_root, allow_shared)
if allow_shared
    return;
end

julia_results_roots = { ...
    fullfile(repo_root, 'benchmark', 'results'), ...
    fullfile(repo_root, 'julia', 'benchmark', 'results') ...
};
matlab_results_root = fullfile(repo_root, 'matlab', 'benchmark', 'results');

out_abs = bsqr_bench_canonical_path(outdir);
matlab_abs = bsqr_bench_canonical_path(matlab_results_root);

for i = 1:numel(julia_results_roots)
    jr = julia_results_roots{i};
    julia_abs = bsqr_bench_canonical_path(jr);
    if startsWith(out_abs, julia_abs)
        error('run_publication_benchmarks:SharedOutdirBlocked', ...
            ['MATLAB benchmarks cannot write into Julia results (%s). ', ...
             'Use matlab/benchmark/results/... or set allow_shared_outdir=true explicitly.'], ...
            jr);
    end
end

if ~startsWith(out_abs, matlab_abs)
    warning('run_publication_benchmarks:OutdirOutsideMatlabTree', ...
        'MATLAB benchmark output is outside matlab/benchmark/results: %s', outdir);
end
end
