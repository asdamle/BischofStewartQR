function summary = check_publication_perf_gate(varargin)
%CHECK_PUBLICATION_PERF_GATE Compare two publication CSVs with slowdown gate.
%
%   summary = check_publication_perf_gate()
%   summary = check_publication_perf_gate(baseline_csv, candidate_csv, max_slowdown)

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
def_csv = fullfile(repo_root, 'matlab', 'benchmark', 'results', 'publication', 'publication_timings.csv');

if nargin >= 1
    baseline_csv = varargin{1};
else
    baseline_csv = def_csv;
end
if nargin >= 2
    candidate_csv = varargin{2};
else
    candidate_csv = def_csv;
end
if nargin >= 3
    max_slowdown = varargin{3};
else
    max_slowdown = 0.05;
end

if ~(isnumeric(max_slowdown) && isscalar(max_slowdown) && isfinite(max_slowdown) && max_slowdown >= 0)
    error('check_publication_perf_gate:InvalidThreshold', ...
        'max_slowdown must be a nonnegative finite scalar.');
end

base = normalize_text_columns(readtable(baseline_csv));
cand = normalize_text_columns(readtable(candidate_csv));

methods = ["bsqr_full", "bsqr_rinv"];
base = base(ismember(base.method, methods), :);
cand = cand(ismember(cand.method, methods), :);

key_cols = {'family', 'regime', 'm', 'n', 'aspect', 'seed', 'method'};
if ismember('bench_surface', base.Properties.VariableNames) || ismember('bench_surface', cand.Properties.VariableNames)
    if ~ismember('bench_surface', base.Properties.VariableNames)
        base.bench_surface = repmat("factor_only", height(base), 1);
    end
    if ~ismember('bench_surface', cand.Properties.VariableNames)
        cand.bench_surface = repmat("factor_only", height(cand), 1);
    end
    key_cols = [key_cols, {'bench_surface'}];
end
base_t = base(:, [key_cols, {'tmed_s'}]);
base_t.Properties.VariableNames{'tmed_s'} = 'tmed_base';
cand_t = cand(:, [key_cols, {'tmed_s'}]);
cand_t.Properties.VariableNames{'tmed_s'} = 'tmed_candidate';

joined = innerjoin(base_t, cand_t, 'Keys', key_cols);
if isempty(joined)
    error('check_publication_perf_gate:NoComparableRows', ...
        'No overlapping method/case rows were found between baseline and candidate CSVs.');
end

joined.slowdown = joined.tmed_candidate ./ joined.tmed_base - 1.0;
joined = sortrows(joined, {'method', 'family', 'regime', 'm', 'n', 'aspect', 'seed'});

viol = joined(isfinite(joined.slowdown) & joined.slowdown > max_slowdown, :);
summary = summarize_slowdown(joined);
summary.max_slowdown = max_slowdown;
summary.num_compared = height(joined);
summary.num_violations = height(viol);
summary.violations = viol;

fprintf('Performance gate compared %d rows (threshold %.2f%% slowdown).\n', ...
    summary.num_compared, 100 * max_slowdown);
disp(summary.by_method);
if ~isempty(viol)
    fprintf('Found %d violating rows.\n', height(viol));
    disp(viol(:, {'method', 'family', 'regime', 'm', 'n', 'aspect', 'seed', 'slowdown'}));
    error('check_publication_perf_gate:Regression', ...
        'Performance gate failed: slowdown exceeded threshold.');
end
fprintf('Performance gate passed.\n');
end

function out = summarize_slowdown(joined)
methods = unique(joined.method);
med = zeros(numel(methods), 1);
mx = zeros(numel(methods), 1);
for i = 1:numel(methods)
    mask = joined.method == methods(i);
    vals = joined.slowdown(mask);
    med(i) = median(vals);
    mx(i) = max(vals);
end
out = struct();
out.by_method = table(methods, med, mx, ...
    'VariableNames', {'method', 'median_slowdown', 'max_slowdown'});
end

function rows = normalize_text_columns(rows)
text_cols = {'family', 'regime', 'method', 'run_id', 'timestamp', 'bench_surface'};
for i = 1:numel(text_cols)
    c = text_cols{i};
    if ~ismember(c, rows.Properties.VariableNames)
        continue;
    end
    v = rows.(c);
    if iscellstr(v)
        rows.(c) = string(v);
    elseif ischar(v)
        rows.(c) = string(cellstr(v));
    elseif iscategorical(v)
        rows.(c) = string(v);
    end
end
end
