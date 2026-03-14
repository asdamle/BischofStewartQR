function plot_publication_results(varargin)
%PLOT_PUBLICATION_RESULTS Generate plots and tables from MATLAB benchmark CSV.
%
%   plot_publication_results()
%   plot_publication_results(csv_path, plots_dir, tables_dir)

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
def_csv = fullfile(repo_root, 'matlab', 'benchmark', 'results', 'publication', 'publication_timings.csv');
def_plots = fullfile(repo_root, 'matlab', 'benchmark', 'results', 'publication', 'plots');
def_tables = fullfile(repo_root, 'matlab', 'benchmark', 'results', 'publication', 'tables');
allow_shared_outdir = bsqr_bench_parse_bool(bsqr_bench_getenv_default('BS_MATLAB_ALLOW_SHARED_OUTDIR', '0'), ...
    'plot_publication_results:InvalidBool');

if nargin >= 1
    csv_path = varargin{1};
else
    csv_path = def_csv;
end
if nargin >= 2
    plots_dir = varargin{2};
else
    plots_dir = def_plots;
end
if nargin >= 3
    tables_dir = varargin{3};
else
    tables_dir = def_tables;
end

validate_plot_targets(plots_dir, tables_dir, repo_root, allow_shared_outdir);

rows = readtable(csv_path);
rows = normalize_text_columns(rows);

if ~isfolder(plots_dir)
    mkdir(plots_dir);
end
if ~isfolder(tables_dir)
    mkdir(tables_dir);
end

make_mode_artifacts(rows, "plain", "bsqr_full", "qr_pivoted", plots_dir, tables_dir);
make_mode_artifacts(rows, "rinv", "bsqr_rinv", "qr_pivoted_trsm", plots_dir, tables_dir);
end

function make_mode_artifacts(rows, mode_name, bs_method, baseline_method, plots_dir, tables_dir)
mode_plots = fullfile(plots_dir, char(mode_name));
mode_tables = fullfile(tables_dir, char(mode_name));
if ~isfolder(mode_plots)
    mkdir(mode_plots);
end
if ~isfolder(mode_tables)
    mkdir(mode_tables);
end

sp = pair_speedups(rows, bs_method, baseline_method);
quality = rows(rows.method == bs_method | rows.method == baseline_method, :);

square = sp(sp.regime == "square", :);
short = sp(sp.regime == "short_wide", :);

write_speedup_table(square, fullfile(mode_tables, 'table_square_speedup.csv'));
write_speedup_table(short, fullfile(mode_tables, 'table_shortwide_speedup.csv'));
write_quality_table(quality, fullfile(mode_tables, 'table_quality.csv'));

plot_square_runtime(rows, bs_method, baseline_method, fullfile(mode_plots, 'figure1_square_runtime'));
plot_shortwide_runtime(rows, bs_method, baseline_method, fullfile(mode_plots, 'figure2_shortwide_runtime'));
plot_shortwide_heatmap(short, fullfile(mode_plots, 'figure3_shortwide_speedup_heatmap'));
plot_quality(quality, bs_method, baseline_method, fullfile(mode_plots, 'figure4_quality'));
plot_aggregate_speedup(sp, fullfile(mode_plots, 'figure5_aggregate_speedup'));
end

function plot_square_runtime(rows, bs_method, baseline_method, outstem)
bs = rows(rows.method == bs_method & rows.regime == "square", :);
base = rows(rows.method == baseline_method & rows.regime == "square", :);
if isempty(bs) || isempty(base)
    return;
end

[ux, bs_t] = grouped_median(bs.m, bs.tmed_s);
[~, base_t] = grouped_median(base.m, base.tmed_s);

fig = figure('Visible', 'off');
loglog(ux, bs_t, '-o', 'LineWidth', 1.5);
hold on;
loglog(ux, base_t, '-s', 'LineWidth', 1.5);
grid on;
xlabel('m = n');
ylabel('median runtime (s)');
legend(char(bs_method), char(baseline_method), 'Location', 'northwest');
title('Square Runtime');
save_figure(fig, outstem);
end

function plot_shortwide_runtime(rows, bs_method, baseline_method, outstem)
bs = rows(rows.method == bs_method & rows.regime == "short_wide", :);
base = rows(rows.method == baseline_method & rows.regime == "short_wide", :);
if isempty(bs) || isempty(base)
    return;
end

[ux, bs_t] = grouped_median(bs.aspect, bs.tmed_s);
[~, base_t] = grouped_median(base.aspect, base.tmed_s);

fig = figure('Visible', 'off');
plot(ux, bs_t, '-o', 'LineWidth', 1.5);
hold on;
plot(ux, base_t, '-s', 'LineWidth', 1.5);
grid on;
xlabel('aspect ratio n/m');
ylabel('median runtime (s)');
legend(char(bs_method), char(baseline_method), 'Location', 'northwest');
title('Short-Wide Runtime');
save_figure(fig, outstem);
end

function plot_shortwide_heatmap(short, outstem)
if isempty(short)
    return;
end

m_vals = unique(short.m);
a_vals = unique(short.aspect);
H = nan(numel(m_vals), numel(a_vals));
for i = 1:numel(m_vals)
    for j = 1:numel(a_vals)
        mask = short.m == m_vals(i) & short.aspect == a_vals(j);
        vals = short.speedup(mask);
        vals = vals(vals > 0 & isfinite(vals));
        if ~isempty(vals)
            H(i, j) = exp(mean(log(vals)));
        end
    end
end

fig = figure('Visible', 'off');
imagesc(a_vals, m_vals, H);
set(gca, 'YDir', 'normal');
colorbar;
xlabel('aspect ratio n/m');
ylabel('m');
title('Short-Wide Speedup Heatmap');
save_figure(fig, outstem);
end

function plot_quality(rows, bs_method, baseline_method, outstem)
if isempty(rows)
    return;
end

fig = figure('Visible', 'off');
mask_bs = rows.method == bs_method;
mask_base = rows.method == baseline_method;
loglog(rows.residual(mask_bs), rows.orthogonality(mask_bs), 'o', 'DisplayName', char(bs_method));
hold on;
loglog(rows.residual(mask_base), rows.orthogonality(mask_base), 's', 'DisplayName', char(baseline_method));
grid on;
xlabel('relative residual');
ylabel('orthogonality error');
legend('Location', 'best');
title('Quality');
save_figure(fig, outstem);
end

function plot_aggregate_speedup(sp, outstem)
if isempty(sp)
    return;
end
keys = unique(sp(:, {'family', 'regime'}), 'rows');
labels = strings(height(keys), 1);
vals = zeros(height(keys), 1);
for i = 1:height(keys)
    mask = sp.family == keys.family(i) & sp.regime == keys.regime(i);
    v = sp.speedup(mask);
    v = v(v > 0 & isfinite(v));
    vals(i) = exp(mean(log(v)));
    labels(i) = keys.family(i) + "/" + keys.regime(i);
end

fig = figure('Visible', 'off');
bar(vals);
xticks(1:numel(vals));
xticklabels(cellstr(labels));
xtickangle(30);
ylabel('geomean speedup (baseline / bsqr)');
title('Aggregate Speedup');
grid on;
save_figure(fig, outstem);
end

function write_speedup_table(sp, path)
if isempty(sp)
    writetable(table(), path);
    return;
end
keys = unique(sp(:, {'family', 'regime', 'm', 'n', 'aspect'}), 'rows');
vals = zeros(height(keys), 1);
for i = 1:height(keys)
    mask = sp.family == keys.family(i) & sp.regime == keys.regime(i) & ...
        sp.m == keys.m(i) & sp.n == keys.n(i) & sp.aspect == keys.aspect(i);
    v = sp.speedup(mask);
    v = v(v > 0 & isfinite(v));
    vals(i) = exp(mean(log(v)));
end
out = [keys, table(vals, 'VariableNames', {'geomean_speedup'})];
writetable(out, path);
end

function write_quality_table(rows, path)
if isempty(rows)
    writetable(table(), path);
    return;
end
keys = unique(rows(:, {'family', 'regime', 'method'}), 'rows');
resid = zeros(height(keys), 1);
orth = zeros(height(keys), 1);
for i = 1:height(keys)
    mask = rows.family == keys.family(i) & rows.regime == keys.regime(i) & rows.method == keys.method(i);
    resid(i) = median(rows.residual(mask));
    orth(i) = median(rows.orthogonality(mask));
end
out = [keys, table(resid, orth, 'VariableNames', {'median_residual', 'median_orthogonality'})];
writetable(out, path);
end

function [u, med] = grouped_median(keys, values)
u = unique(keys);
med = zeros(numel(u), 1);
for i = 1:numel(u)
    med(i) = median(values(keys == u(i)));
end
end

function save_figure(fig, outstem)
set(fig, 'Color', 'w');
print(fig, [outstem, '.png'], '-dpng', '-r180');
print(fig, [outstem, '.pdf'], '-dpdf');
close(fig);
end

function sp = pair_speedups(rows, bs_method, baseline_method)
rows_bs = rows(rows.method == bs_method, :);
rows_base = rows(rows.method == baseline_method, :);
if isempty(rows_bs) || isempty(rows_base)
    sp = table();
    return;
end

key_cols = {'family', 'regime', 'm', 'n', 'aspect', 'seed'};
join_bs = rows_bs(:, [key_cols, {'tmed_s'}]);
join_bs.Properties.VariableNames{'tmed_s'} = 'tmed_bs';
join_base = rows_base(:, [key_cols, {'tmed_s'}]);
join_base.Properties.VariableNames{'tmed_s'} = 'tmed_base';
joined = innerjoin(join_bs, join_base, 'Keys', key_cols);
joined.speedup = joined.tmed_base ./ joined.tmed_bs;
sp = joined(:, [key_cols, {'speedup'}]);
end

function rows = normalize_text_columns(rows)
text_cols = {'family', 'regime', 'method', 'run_id', 'timestamp'};
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

function validate_plot_targets(plots_dir, tables_dir, repo_root, allow_shared)
if allow_shared
    return;
end

julia_results_root = fullfile(repo_root, 'benchmark', 'results');
targets = {plots_dir, tables_dir};
for i = 1:numel(targets)
    tgt_abs = bsqr_bench_canonical_path(targets{i});
    julia_abs = bsqr_bench_canonical_path(julia_results_root);
    if startsWith(tgt_abs, julia_abs)
        error('plot_publication_results:SharedOutdirBlocked', ...
            ['MATLAB plot/table outputs cannot be written into Julia results (%s). ', ...
             'Use matlab/benchmark/results/... or set BS_MATLAB_ALLOW_SHARED_OUTDIR=1 explicitly.'], ...
            julia_results_root);
    end
end
end
