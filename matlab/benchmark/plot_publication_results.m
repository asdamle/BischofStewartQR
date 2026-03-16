function plot_publication_results(varargin)
%PLOT_PUBLICATION_RESULTS Generate publication plots and tables from benchmark CSV.
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
rows = validate_bench_surface_column(rows);
style = publication_style();

if ~isfolder(plots_dir)
    mkdir(plots_dir);
end
if ~isfolder(tables_dir)
    mkdir(tables_dir);
end

make_mode_artifacts(rows, "plain", "bsqr_full", "qr_pivoted", plots_dir, tables_dir, style);
make_mode_artifacts(rows, "rinv", "bsqr_rinv", "qr_pivoted_trsm", plots_dir, tables_dir, style);
end

function make_mode_artifacts(rows, mode_name, bs_method, baseline_method, plots_dir, tables_dir, style)
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

plot_square_runtime(rows, bs_method, baseline_method, fullfile(mode_plots, 'figure1_square_runtime'), style);
plot_shortwide_runtime(rows, bs_method, baseline_method, fullfile(mode_plots, 'figure2_shortwide_runtime'), style);
plot_shortwide_heatmap(short, baseline_method, bs_method, fullfile(mode_plots, 'figure3_shortwide_speedup_heatmap'), style);
plot_quality(quality, bs_method, baseline_method, fullfile(mode_plots, 'figure4_quality'), style);
plot_aggregate_speedup(sp, baseline_method, bs_method, fullfile(mode_plots, 'figure5_aggregate_speedup'), style);
end

function plot_square_runtime(rows, bs_method, baseline_method, outstem, style)
bs = rows(rows.method == bs_method & rows.regime == "square", :);
base = rows(rows.method == baseline_method & rows.regime == "square", :);
if isempty(bs) || isempty(base)
    return;
end

[x_bs, med_bs, lo_bs, hi_bs] = grouped_stats(bs.m, bs.tmed_s);
[x_base, med_base, lo_base, hi_base] = grouped_stats(base.m, base.tmed_s);

fig = figure('Visible', 'off', 'Color', 'w');
ax = axes(fig);
hold(ax, 'on');

h_bs = plot_with_band(ax, x_bs, med_bs, lo_bs, hi_bs, method_style(style, bs_method));
h_base = plot_with_band(ax, x_base, med_base, lo_base, hi_base, method_style(style, baseline_method));

set(ax, 'XScale', 'log', 'YScale', 'log');
apply_axes_style(ax, style);
xlabel(ax, 'm=n', 'Interpreter', 'none');
ylabel(ax, 'median time (s)', 'Interpreter', 'none');
title(ax, 'Square Runtime', 'Interpreter', 'none');
lgd = legend(ax, [h_bs, h_base], {method_label(bs_method), method_label(baseline_method)}, ...
    'Location', 'northwest', 'Interpreter', 'none', 'Box', 'on');
set(lgd, 'FontSize', style.legend_font_size);

save_figure(fig, outstem);
end

function plot_shortwide_runtime(rows, bs_method, baseline_method, outstem, style)
short = rows(rows.regime == "short_wide" & (rows.method == bs_method | rows.method == baseline_method), :);
if isempty(short)
    return;
end

m_vals = unique(short.m);
m_vals = sort(m_vals);
if isempty(m_vals)
    return;
end

cmap = lines(numel(m_vals));
fig = figure('Visible', 'off', 'Color', 'w');
ax = axes(fig);
hold(ax, 'on');

methods = [bs_method, baseline_method];
for mi = 1:numel(m_vals)
    m = m_vals(mi);
    this_color = cmap(mi, :);
    for method = methods
        block = short(short.m == m & short.method == method, :);
        if isempty(block)
            continue;
        end
        [x_n, med_n, lo_n, hi_n] = grouped_stats(block.n, block.tmed_s);
        st = method_style(style, method);
        st.color = this_color;
        plot_with_band(ax, x_n, med_n, lo_n, hi_n, st);
    end
end

set(ax, 'XScale', 'log', 'YScale', 'log');
apply_axes_style(ax, style);
xlabel(ax, 'n', 'Interpreter', 'none');
ylabel(ax, 'median time (s)', 'Interpreter', 'none');
title(ax, 'Short-Wide Runtime', 'Interpreter', 'none');

method_handles = gobjects(0, 1);
method_labels = strings(0, 1);
for method = methods
    st = method_style(style, method);
    method_handles(end+1, 1) = plot(ax, nan, nan, ...
        'LineStyle', st.linestyle, 'Marker', st.marker, ...
        'Color', [0 0 0], 'LineWidth', style.line_width, ...
        'MarkerSize', style.marker_size); %#ok<AGROW>
    method_labels(end+1, 1) = method_label(method); %#ok<AGROW>
end

m_handles = gobjects(0, 1);
m_labels = strings(0, 1);
for mi = 1:numel(m_vals)
    m_handles(end+1, 1) = plot(ax, nan, nan, '-', ...
        'Color', cmap(mi, :), 'LineWidth', style.line_width); %#ok<AGROW>
    m_labels(end+1, 1) = "m=" + string(m_vals(mi)); %#ok<AGROW>
end

legend_handles = [method_handles; m_handles];
legend_labels = [method_labels; m_labels];
lgd = legend(ax, legend_handles, cellstr(legend_labels), ...
    'Location', 'southoutside', 'NumColumns', min(4, numel(legend_labels)), ...
    'Interpreter', 'none', 'Box', 'on');
set(lgd, 'FontSize', style.legend_font_size);

save_figure(fig, outstem);
end

function plot_shortwide_heatmap(short, baseline_method, bs_method, outstem, style)
if isempty(short)
    return;
end

m_vals = sort(unique(short.m));
a_vals = sort(unique(short.aspect));
H = nan(numel(m_vals), numel(a_vals));
for i = 1:numel(m_vals)
    for j = 1:numel(a_vals)
        mask = short.m == m_vals(i) & short.aspect == a_vals(j);
        vals = short.speedup(mask);
        vals = vals(vals > 0 & isfinite(vals));
        if ~isempty(vals)
            H(i, j) = geomean(vals);
        end
    end
end

fig = figure('Visible', 'off', 'Color', 'w');
ax = axes(fig);
imagesc(ax, a_vals, m_vals, H);
set(ax, 'YDir', 'normal');
apply_axes_style(ax, style);

finite_vals = H(isfinite(H));
if isempty(finite_vals)
    caxis(ax, [0.9, 1.1]);
else
    lo = min(finite_vals);
    hi = max(finite_vals);
    if lo == hi
        lo = min(lo, 0.95);
        hi = max(hi, 1.05);
    end
    lo = min(lo, 0.99);
    hi = max(hi, 1.01);
    caxis(ax, [lo, hi]);
end
colormap(ax, red_white_green(256));
cb = colorbar(ax);
ylabel(cb, "speedup " + baseline_method + "/" + bs_method, 'Interpreter', 'none');

xlabel(ax, 'aspect (n/m)', 'Interpreter', 'none');
ylabel(ax, 'm', 'Interpreter', 'none');
title(ax, 'Short-Wide Speedup Heatmap', 'Interpreter', 'none');

save_figure(fig, outstem);
end

function plot_quality(rows, bs_method, baseline_method, outstem, style)
if isempty(rows)
    return;
end

regimes = ["square", "short_wide"];
regime_labels = {'square', 'short-wide'};
metrics = {'residual', 'orthogonality'};
titles = {'Residual', 'Orthogonality'};

fig = figure('Visible', 'off', 'Color', 'w');
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
axes_handles = gobjects(numel(metrics), 1);

for mi = 1:numel(metrics)
    metric = metrics{mi};
    ax = nexttile(tl);
    axes_handles(mi) = ax;
    hold(ax, 'on');

    y_bs = zeros(1, numel(regimes));
    y_base = zeros(1, numel(regimes));
    y_bs_lo = zeros(1, numel(regimes));
    y_bs_hi = zeros(1, numel(regimes));
    y_base_lo = zeros(1, numel(regimes));
    y_base_hi = zeros(1, numel(regimes));

    for ri = 1:numel(regimes)
        regime = regimes(ri);
        bs_vals = extract_metric(rows, metric, regime, bs_method);
        base_vals = extract_metric(rows, metric, regime, baseline_method);

        [y_bs(ri), y_bs_lo(ri), y_bs_hi(ri)] = metric_stats(bs_vals);
        [y_base(ri), y_base_lo(ri), y_base_hi(ri)] = metric_stats(base_vals);
    end

    x = 1:numel(regimes);
    w = 0.36;
    bar(ax, x - w/2, y_bs, w, 'FaceColor', style.bsqr_color, 'DisplayName', method_label(bs_method));
    bar(ax, x + w/2, y_base, w, 'FaceColor', style.baseline_color, 'DisplayName', method_label(baseline_method));

    err_bs = errorbar(ax, x - w/2, y_bs, y_bs - y_bs_lo, y_bs_hi - y_bs, ...
        'k.', 'CapSize', 5, 'LineWidth', 0.9);
    err_base = errorbar(ax, x + w/2, y_base, y_base - y_base_lo, y_base_hi - y_base, ...
        'k.', 'CapSize', 5, 'LineWidth', 0.9);
    set(err_bs, 'HandleVisibility', 'off');
    set(err_base, 'HandleVisibility', 'off');

    set(ax, 'YScale', 'log');
    xticks(ax, x);
    xticklabels(ax, regime_labels);
    apply_axes_style(ax, style);
    title(ax, titles{mi}, 'Interpreter', 'none');
    if mi == numel(metrics)
        xlabel(ax, 'regime', 'Interpreter', 'none');
    end
    ylabel(ax, metric, 'Interpreter', 'none');
    grid(ax, 'on');
end

lgd = legend(axes_handles(1), 'Location', 'northwest', 'Interpreter', 'none', 'Box', 'on');
set(lgd, 'FontSize', style.legend_font_size);

save_figure(fig, outstem);
end

function plot_aggregate_speedup(sp, baseline_method, bs_method, outstem, style)
if isempty(sp)
    return;
end

has_threads = ismember('blas_threads', sp.Properties.VariableNames);

if has_threads
    combo_keys = unique(sp(:, {'family', 'regime', 'blas_threads'}), 'rows');
else
    combo_keys = unique(sp(:, {'family', 'regime'}), 'rows');
end

labels = strings(height(combo_keys), 1);
center = nan(height(combo_keys), 1);
lo = nan(height(combo_keys), 1);
hi = nan(height(combo_keys), 1);

for i = 1:height(combo_keys)
    if has_threads
        combo_mask = sp.family == combo_keys.family(i) & sp.regime == combo_keys.regime(i) & ...
            sp.blas_threads == combo_keys.blas_threads(i);
    else
        combo_mask = sp.family == combo_keys.family(i) & sp.regime == combo_keys.regime(i);
    end

    combo_rows = sp(combo_mask, :);
    seeds = unique(combo_rows.seed);
    seed_speedups = nan(numel(seeds), 1);
    for sidx = 1:numel(seeds)
        vals = combo_rows.speedup(combo_rows.seed == seeds(sidx));
        vals = vals(vals > 0 & isfinite(vals));
        seed_speedups(sidx) = geomean(vals);
    end
    seed_speedups = seed_speedups(isfinite(seed_speedups));
    if isempty(seed_speedups)
        continue;
    end

    center(i) = geomean(seed_speedups);
    lo(i) = min(seed_speedups);
    hi(i) = max(seed_speedups);

    family_short = short_family_name(combo_keys.family(i));
    regime_short = "sq";
    if combo_keys.regime(i) == "short_wide"
        regime_short = "sw";
    end
    if has_threads
        labels(i) = family_short + " | " + regime_short + " | t" + string(combo_keys.blas_threads(i));
    else
        labels(i) = family_short + " | " + regime_short;
    end
end

valid = isfinite(center);
labels = labels(valid);
center = center(valid);
lo = lo(valid);
hi = hi(valid);

if isempty(center)
    return;
end

fig = figure('Visible', 'off', 'Color', 'w');
ax = axes(fig);
y = 1:numel(center);
barh(ax, y, center, 'FaceColor', style.aggregate_color);
hold(ax, 'on');
errorbar(ax, center, y, center - lo, hi - center, 'horizontal', ...
    'k.', 'CapSize', 4, 'LineWidth', 0.9);
xline(ax, 1.0, '--k', 'LineWidth', 1.0);
yticks(ax, y);
yticklabels(ax, cellstr(labels));
apply_axes_style(ax, style);
xlabel(ax, "geomean speedup (" + baseline_method + "/" + bs_method + ")", 'Interpreter', 'none');
title(ax, 'Aggregate Speedup by Family/Regime', 'Interpreter', 'none');

save_figure(fig, outstem);
end

function write_speedup_table(sp, path)
if isempty(sp)
    writetable(table(), path);
    return;
end

key_cols = {'family', 'regime', 'm', 'n', 'aspect'};
if ismember('blas_threads', sp.Properties.VariableNames)
    key_cols = [key_cols, {'blas_threads'}];
end

keys = unique(sp(:, key_cols), 'rows');
vals = zeros(height(keys), 1);
for i = 1:height(keys)
    mask = true(height(sp), 1);
    for j = 1:numel(key_cols)
        col = key_cols{j};
        mask = mask & sp.(col) == keys.(col)(i);
    end
    v = sp.speedup(mask);
    v = v(v > 0 & isfinite(v));
    vals(i) = geomean(v);
end
out = [keys, table(vals, 'VariableNames', {'geomean_speedup'})];
writetable(out, path);
end

function write_quality_table(rows, path)
if isempty(rows)
    writetable(table(), path);
    return;
end

key_cols = {'family', 'regime', 'method'};
if ismember('blas_threads', rows.Properties.VariableNames)
    key_cols = {'family', 'regime', 'method', 'blas_threads'};
end

keys = unique(rows(:, key_cols), 'rows');
resid = zeros(height(keys), 1);
orth = zeros(height(keys), 1);
for i = 1:height(keys)
    mask = true(height(rows), 1);
    for j = 1:numel(key_cols)
        col = key_cols{j};
        mask = mask & rows.(col) == keys.(col)(i);
    end
    resid(i) = median(rows.residual(mask));
    orth(i) = median(rows.orthogonality(mask));
end
out = [keys, table(resid, orth, 'VariableNames', {'median_residual', 'median_orthogonality'})];
writetable(out, path);
end

function [x, med, lo, hi] = grouped_stats(keys, values)
x = sort(unique(keys));
med = nan(size(x));
lo = nan(size(x));
hi = nan(size(x));
for i = 1:numel(x)
    v = values(keys == x(i));
    v = v(isfinite(v));
    if isempty(v)
        continue;
    end
    med(i) = median(v);
    lo(i) = min(v);
    hi(i) = max(v);
end
end

function vals = extract_metric(rows, metric, regime, method)
vals = rows.(metric)(rows.regime == regime & rows.method == method);
vals = vals(isfinite(vals));
end

function [m, p10, p90] = metric_stats(vals)
if isempty(vals)
    m = nan;
    p10 = nan;
    p90 = nan;
    return;
end
m = median(vals);
p10 = prctile(vals, 10);
p90 = prctile(vals, 90);
end

function h_line = plot_with_band(ax, x, med, lo, hi, st)
h_line = gobjects(1, 1);
good = isfinite(x) & isfinite(med) & isfinite(lo) & isfinite(hi) & med > 0 & lo > 0 & hi > 0;
if ~any(good)
    return;
end
x = x(good);
med = med(good);
lo = lo(good);
hi = hi(good);

h_fill = fill(ax, [x; flipud(x)], [lo; flipud(hi)], st.color, ...
    'FaceAlpha', 0.14, 'EdgeColor', 'none');
set(h_fill, 'HandleVisibility', 'off');
h_line = plot(ax, x, med, ...
    'LineStyle', st.linestyle, ...
    'Marker', st.marker, ...
    'Color', st.color, ...
    'LineWidth', st.linewidth, ...
    'MarkerSize', st.markersize);
end

function st = method_style(style, method)
st = struct('color', style.bsqr_color, 'linestyle', '-', 'marker', 'o', ...
    'linewidth', style.line_width, 'markersize', style.marker_size);
if method == "qr_pivoted" || method == "qr_pivoted_trsm"
    st.color = style.baseline_color;
    st.linestyle = '--';
    st.marker = 's';
end
if method == "bsqr_rinv" || method == "qr_pivoted_trsm"
    st.linestyle = ':';
end
end

function label = method_label(method)
switch char(method)
    case 'bsqr_full'
        label = 'BSQR';
    case 'bsqr_rinv'
        label = 'BSQR+RINV';
    case 'qr_pivoted'
        label = 'QR';
    case 'qr_pivoted_trsm'
        label = 'QR+TRSM';
    otherwise
        label = char(method);
end
end

function name = short_family_name(family)
switch char(family)
    case 'ill_conditioned'
        name = "ill-cond";
    case 'orthonormal_rows'
        name = "orth-rows";
    otherwise
        name = string(family);
end
end

function v = geomean(x)
x = x(isfinite(x) & x > 0);
if isempty(x)
    v = nan;
else
    v = exp(mean(log(x)));
end
end

function c = red_white_green(n)
if nargin < 1
    n = 256;
end
n = max(3, round(n));
half = floor(n / 2);
upper = n - half;
red = [ones(half, 1), linspace(0, 1, half)', linspace(0, 1, half)'];
green = [linspace(1, 0, upper)', ones(upper, 1), linspace(1, 0, upper)'];
c = [red; green];
c = c(1:n, :);
end

function apply_axes_style(ax, style)
set(ax, ...
    'Box', 'on', ...
    'LineWidth', style.axis_line_width, ...
    'FontSize', style.axis_font_size, ...
    'GridAlpha', 0.25, ...
    'MinorGridAlpha', 0.15, ...
    'TickLabelInterpreter', 'none');
grid(ax, 'on');
end

function style = publication_style()
style = struct();
style.bsqr_color = [76, 120, 168] / 255;
style.baseline_color = [245, 133, 24] / 255;
style.aggregate_color = [114, 183, 178] / 255;
style.axis_font_size = 10;
style.legend_font_size = 9;
style.line_width = 1.5;
style.marker_size = 5;
style.axis_line_width = 0.9;
end

function save_figure(fig, outstem)
set(fig, 'Color', 'w');
print(fig, [outstem, '.png'], '-dpng', '-r220');
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
if ismember('blas_threads', rows.Properties.VariableNames)
    key_cols = [key_cols, {'blas_threads'}];
end
if ismember('bench_surface', rows.Properties.VariableNames)
    key_cols = [key_cols, {'bench_surface'}];
end

join_bs = rows_bs(:, [key_cols, {'tmed_s'}]);
join_bs.Properties.VariableNames{'tmed_s'} = 'tmed_bs';
join_base = rows_base(:, [key_cols, {'tmed_s'}]);
join_base.Properties.VariableNames{'tmed_s'} = 'tmed_base';
joined = innerjoin(join_bs, join_base, 'Keys', key_cols);
joined.speedup = joined.tmed_base ./ joined.tmed_bs;
sp = joined(:, [key_cols, {'speedup'}]);
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

function rows = validate_bench_surface_column(rows)
if ~ismember('bench_surface', rows.Properties.VariableNames)
    return;
end
vals = unique(rows.bench_surface);
if numel(vals) ~= 1
    error('plot_publication_results:MixedBenchSurface', ...
        ['CSV contains multiple bench_surface values. ', ...
         'Plot one surface at a time to avoid mixed speedup aggregates.']);
end
end

function validate_plot_targets(plots_dir, tables_dir, repo_root, allow_shared)
if allow_shared
    return;
end

julia_results_roots = { ...
    fullfile(repo_root, 'benchmark', 'results'), ...
    fullfile(repo_root, 'julia', 'benchmark', 'results') ...
};
targets = {plots_dir, tables_dir};
for i = 1:numel(targets)
    tgt_abs = bsqr_bench_canonical_path(targets{i});
    for j = 1:numel(julia_results_roots)
        jr = julia_results_roots{j};
        julia_abs = bsqr_bench_canonical_path(jr);
        if startsWith(tgt_abs, julia_abs)
            error('plot_publication_results:SharedOutdirBlocked', ...
                ['MATLAB plot/table outputs cannot be written into Julia results (%s). ', ...
                 'Use matlab/benchmark/results/... or set BS_MATLAB_ALLOW_SHARED_OUTDIR=1 explicitly.'], ...
                jr);
        end
    end
end
end
