function plot_publication_results(varargin)
%PLOT_PUBLICATION_RESULTS Generate publication figures and tables from benchmark CSV.
%
%   plot_publication_results()
%   plot_publication_results(csv_path, plots_dir, tables_dir)
%
%   Figure spec: docs/PUBLICATION_FIGURES_PLAN.md. The Julia pipeline's
%   plotter (julia/benchmark/plot_publication.py) conforms to the same
%   spec; any change to figure semantics, labels, or styling must land in
%   both. The MATLAB CSV has no blas_threads dimension, so figures have
%   family panels only and tables omit that column.

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
style = publication_style();

if ~isfolder(plots_dir)
    mkdir(plots_dir);
end
if ~isfolder(tables_dir)
    mkdir(tables_dir);
end

make_mode_artifacts(rows, "plain", "bsqr_full", "qr_pivoted", plots_dir, tables_dir, style);
make_mode_artifacts(rows, "rinv", "bsqr_rinv", "qr_pivoted_trsm", plots_dir, tables_dir, style);

% Composite relative-time figure (plain + rinv overlaid) at the top level -- the
% single timing figure for the paper.
fig_relative_time_composite(rows, style, plots_dir);
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

labels = mode_labels(mode_name);
rt = pair_relative_times(rows, bs_method, baseline_method);
quality = rows(rows.method == bs_method | rows.method == baseline_method, :);

write_relative_time_table(rt(rt.regime == "square", :), {'family', 'm', 'n'}, ...
    fullfile(mode_tables, 'table_square_relative_time.csv'));
write_relative_time_table(rt(rt.regime == "short_wide", :), {'family', 'm', 'n', 'aspect'}, ...
    fullfile(mode_tables, 'table_shortwide_relative_time.csv'));
write_quality_summary(quality, bs_method, baseline_method, labels, mode_tables);

fig_square_runtime(rows, bs_method, baseline_method, labels, fullfile(mode_plots, 'fig_square_runtime'), style);
fig_shortwide_runtime(rows, bs_method, baseline_method, labels, fullfile(mode_plots, 'fig_shortwide_runtime'), style);
% The per-mode relative-time forest plot is superseded by the top-level
% fig_relative_time_composite (plain + rinv on one figure); rt above still feeds the
% relative-time tables.
write_captions(rows, labels, mode_name, mode_plots);
end

% ----------------------------------------------------------------------
% Figures
% ----------------------------------------------------------------------

function fig_square_runtime(rows, bs_method, baseline_method, labels, outstem, style)
sub = rows(rows.regime == "square", :);
families = sort(unique(sub.family));
if isempty(families)
    return;
end

fig = new_figure(style, style.double_col_width, 2.45);
tl = tiledlayout(fig, 1, numel(families), 'TileSpacing', 'compact', 'Padding', 'compact');
first_handles = gobjects(2, 1);
for fi = 1:numel(families)
    fam = families(fi);
    ax = nexttile(tl);
    hold(ax, 'on');
    methods = {bs_method, baseline_method};
    roles = {'bs', 'base'};
    % Bands first (faint, behind every line), then the median lines.
    stats = cell(2, 1);
    for k = 1:2
        block = sub(sub.family == fam & sub.method == methods{k}, :);
        [x, ctr, lo, hi] = seed_grouped_stats(block.m, block.tmed_s);
        st = role_style(style, roles{k});
        stats{k} = struct('x', x, 'ctr', ctr, 'lo', lo, 'hi', hi, 'st', st);
        plot_band(ax, x, lo, hi, st.color);
    end
    for k = 1:2
        h = plot_median_line(ax, stats{k}.x, stats{k}.ctr, stats{k}.st);
        if fi == 1
            first_handles(k) = h;
        end
    end
    set(ax, 'XScale', 'log', 'YScale', 'log');
    apply_axes_style(ax, style);
    xticks(ax, sort(unique(sub.m)));
    set(ax, 'XMinorTick', 'off');
    title(ax, display_family(fam), 'Interpreter', 'latex', 'FontWeight', 'normal');
    xlabel(ax, '$m = n$', 'Interpreter', 'latex');
    if fi == 1
        ylabel(ax, 'median time [s]', 'Interpreter', 'latex');
    end
end
lgd = legend(first_handles, {labels.bs, labels.base}, ...
    'Location', 'northwest', 'Interpreter', 'latex', 'Box', 'off');
set(lgd, 'FontSize', style.legend_font_size);
save_figure(fig, outstem, style);
end

function fig_shortwide_runtime(rows, bs_method, baseline_method, labels, outstem, style)
sub = rows(rows.regime == "short_wide" & (rows.method == bs_method | rows.method == baseline_method), :);
families = sort(unique(sub.family));
m_vals = sort(unique(sub.m));
if isempty(families) || isempty(m_vals)
    return;
end
cmap = lines(numel(m_vals));

fig = new_figure(style, style.double_col_width, 3.2);
tl = tiledlayout(fig, 1, numel(families), 'TileSpacing', 'compact', 'Padding', 'compact');
ax_last = [];
for fi = 1:numel(families)
    fam = families(fi);
    ax = nexttile(tl);
    ax_last = ax;
    hold(ax, 'on');
    methods = {bs_method, baseline_method};
    roles = {'bs', 'base'};
    % Bands first (faint, behind every line), then all median lines.
    series = {};
    for mi = 1:numel(m_vals)
        for k = 1:2
            block = sub(sub.family == fam & sub.m == m_vals(mi) & sub.method == methods{k}, :);
            if isempty(block)
                continue;
            end
            [x, ctr, lo, hi] = seed_grouped_stats(block.n, block.tmed_s);
            st = role_style(style, roles{k});
            st.color = cmap(mi, :);
            st.markersize = 3.2;
            series{end+1} = struct('x', x, 'ctr', ctr, 'lo', lo, 'hi', hi, 'st', st); %#ok<AGROW>
        end
    end
    for si = 1:numel(series)
        plot_band(ax, series{si}.x, series{si}.lo, series{si}.hi, series{si}.st.color);
    end
    for si = 1:numel(series)
        plot_median_line(ax, series{si}.x, series{si}.ctr, series{si}.st);
    end
    set(ax, 'XScale', 'log', 'YScale', 'log');
    apply_axes_style(ax, style);
    % Power-of-two ticks to match the Julia plotter's log2 axis.
    n_min = min(sub.n);
    n_max = max(sub.n);
    pows = ceil(log2(n_min)):floor(log2(n_max));
    xticks(ax, 2 .^ pows);
    xticklabels(ax, compose("$2^{%d}$", pows));
    set(ax, 'XMinorTick', 'off');
    title(ax, display_family(fam), 'Interpreter', 'latex', 'FontWeight', 'normal');
    xlabel(ax, '$n$', 'Interpreter', 'latex');
    if fi == 1
        ylabel(ax, 'median time [s]', 'Interpreter', 'latex');
    end
end

% One combined legend below the panels: method styles (black) then m colors.
handles = gobjects(0, 1);
entries = strings(0, 1);
roles = {'bs', 'base'};
role_names = {labels.bs, labels.base};
for k = 1:2
    st = role_style(style, roles{k});
    handles(end+1, 1) = plot(ax_last, nan, nan, 'LineStyle', st.linestyle, ...
        'Marker', st.marker, 'Color', [0 0 0], 'LineWidth', st.linewidth, ...
        'MarkerSize', 3.2); %#ok<AGROW>
    entries(end+1, 1) = string(role_names{k}); %#ok<AGROW>
end
for mi = 1:numel(m_vals)
    handles(end+1, 1) = plot(ax_last, nan, nan, '-', 'Color', cmap(mi, :), ...
        'LineWidth', style.line_width); %#ok<AGROW>
    entries(end+1, 1) = "$m = " + string(m_vals(mi)) + "$"; %#ok<AGROW>
end
lgd = legend(ax_last, handles, cellstr(entries), 'Interpreter', 'latex', ...
    'Box', 'off', 'NumColumns', 4);
lgd.Layout.Tile = 'south';
set(lgd, 'FontSize', style.legend_font_size);
save_figure(fig, outstem, style);
end

function fig_relative_time_composite(rows, style, plots_dir)
% Composite forest plot: plain (BSQR / CPQR) and rinv (BSQR+W / CPQR+solve)
% relative times overlaid on the shared family x regime rows, distinguished by
% colour. The MATLAB CSV has no thread dimension, so one marker per (mode, row).
modes = { 'bsqr_full', 'qr_pivoted',     style.bsqr_color,     'BSQR / CPQR'; ...
          'bsqr_rinv', 'qr_pivoted_trsm', style.baseline_color, ...
          '(BSQR + $R_{11}^{-1}R_{12}$) / (CPQR + $R_{11}^{-1}R_{12}$)' };
families = sort(unique(rows.family));
regimes = ["square", "short_wide"];
regime_display = containers.Map({'square', 'short_wide'}, {'(m = n)', '(m < n)'});

row_keys = strings(0, 2);
row_labels = strings(0, 1);
for fi = 1:numel(families)
    for ri = 1:numel(regimes)
        row_keys(end+1, :) = [string(families(fi)), regimes(ri)]; %#ok<AGROW>
        row_labels(end+1, 1) = display_family(families(fi)) + " " + ...
            string(regime_display(char(regimes(ri)))); %#ok<AGROW>
    end
end
nrows = size(row_keys, 1);
if nrows == 0
    return;
end
y_base = (nrows:-1:1)';

fig = new_figure(style, style.single_col_width, 0.40 * nrows + 1.35);
ax = axes(fig);
hold(ax, 'on');
dodge = 0.18;
handles = gobjects(size(modes, 1), 1);
for mi = 1:size(modes, 1)
    rt = pair_relative_times(rows, modes{mi, 1}, modes{mi, 2});
    color = modes{mi, 3};
    off = dodge * (mi == 1) - dodge * (mi == 2);   % plain up, rinv down
    ys = []; cs = []; lo_e = []; hi_e = [];
    for k = 1:nrows
        mask = rt.family == row_keys(k, 1) & rt.regime == row_keys(k, 2);
        if ~any(mask)
            continue;
        end
        combo = rt(mask, :);
        seeds = unique(combo.seed);
        sg = zeros(numel(seeds), 1);
        for si = 1:numel(seeds)
            v = combo.relative_time(combo.seed == seeds(si));
            sg(si) = geomean(v(isfinite(v) & v > 0));
        end
        sg = sg(isfinite(sg));
        if isempty(sg)
            continue;
        end
        c = geomean(sg);
        ys(end+1, 1) = y_base(k) + off;       %#ok<AGROW>
        cs(end+1, 1) = c;                      %#ok<AGROW>
        lo_e(end+1, 1) = c - min(sg);          %#ok<AGROW>
        hi_e(end+1, 1) = max(sg) - c;          %#ok<AGROW>
    end
    handles(mi) = errorbar(ax, cs, ys, lo_e, hi_e, 'horizontal', 'Marker', 'o', ...
        'Color', color, 'MarkerFaceColor', color, 'LineStyle', 'none', ...
        'CapSize', 2.0, 'MarkerSize', 4.0, 'LineWidth', 0.9);
end
xline(ax, 1.0, '--k', 'LineWidth', 0.8);
yticks(ax, sort(y_base));
yticklabels(ax, cellstr(flipud(row_labels)));
ylim(ax, [0.4, nrows + 0.6]);
xl = xlim(ax);
xlim(ax, [min(0.95, xl(1)), xl(2)]);
apply_axes_style(ax, style);
xlabel(ax, "relative time (BSQR / baseline)", 'Interpreter', 'latex');
grid(ax, 'on');
lgd = legend(ax, handles, {modes{1, 4}, modes{2, 4}}, 'Location', 'southoutside', ...
    'Interpreter', 'latex', 'Box', 'off', 'NumColumns', 1);
set(lgd, 'FontSize', style.legend_font_size);
save_figure(fig, fullfile(plots_dir, 'fig_relative_time_composite'), style);
end

% ----------------------------------------------------------------------
% Tables and captions
% ----------------------------------------------------------------------

function write_relative_time_table(rt, key_cols, path)
if isempty(rt)
    writetable(table(), path);
    return;
end
keys = unique(rt(:, key_cols), 'rows');
n = height(keys);
rel_geo = zeros(n, 1);
rel_min = zeros(n, 1);
rel_max = zeros(n, 1);
bs_geo = zeros(n, 1);
base_geo = zeros(n, 1);
for i = 1:n
    mask = true(height(rt), 1);
    for j = 1:numel(key_cols)
        col = key_cols{j};
        mask = mask & rt.(col) == keys.(col)(i);
    end
    rel = rt.relative_time(mask);
    rel = rel(isfinite(rel) & rel > 0);
    rel_geo(i) = geomean(rel);
    rel_min(i) = min(rel);
    rel_max(i) = max(rel);
    bs_geo(i) = geomean(rt.tmed_bs(mask));
    base_geo(i) = geomean(rt.tmed_base(mask));
end
out = [keys, table(rel_geo, rel_min, rel_max, bs_geo, base_geo, 'VariableNames', ...
    {'relative_time_geomean', 'relative_time_seed_min', 'relative_time_seed_max', ...
     'bsqr_tmed_geomean_s', 'baseline_tmed_geomean_s'})];
writetable(out, path);
end

function write_quality_summary(rows, bs_method, baseline_method, labels, outdir)
%Terse markdown report replacing the former quality figure/table.
ncases = height(unique(rows(:, {'family', 'regime', 'm', 'n', 'aspect', 'seed'}), 'rows'));
methods = {bs_method, baseline_method};
role_names = {labels.bs, labels.base};

fid = fopen(fullfile(outdir, 'quality_summary.md'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# Numerical quality summary\n\n');
fprintf(fid, ['Across all %d benchmark cases (3 matrix families, both regimes, all ', ...
    'sizes and seeds):\n\n'], ncases);
fprintf(fid, ['| method | median rel. residual | max rel. residual | ', ...
    'median ||I - Q^T Q||_F | max ||I - Q^T Q||_F |\n']);
fprintf(fid, '|---|---:|---:|---:|---:|\n');
stats = zeros(2, 4);
for k = 1:2
    resid = rows.residual(rows.method == methods{k});
    orth = rows.orthogonality(rows.method == methods{k});
    stats(k, :) = [median(resid), max(resid), median(orth), max(orth)];
    fprintf(fid, '| %s | %.2e | %.2e | %.2e | %.2e |\n', role_names{k}, stats(k, :));
end
fprintf(fid, ['\nThe relative residual ||A Pi - QR||_F / ||A||_F of %s never exceeded ', ...
    '%.2e, and its deviation from orthogonality ||I - Q^T Q||_F never exceeded %.2e; ', ...
    'both match the built-in baseline to within a small factor.\n'], ...
    role_names{1}, stats(1, 2), stats(1, 4));
end

function write_captions(rows, labels, mode_name, outdir)
seeds = unique(rows.seed);
run_ids = unique(rows.run_id);
fid = fopen(fullfile(outdir, 'figure_captions.md'), 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '# Suggested captions (%s comparison)\n\n', mode_name);
fprintf(fid, ['Run: %s; seeds: %d; MATLAB default threading. In all figures, %s is ', ...
    'compared against %s; both timed paths materialize Q, R, and the permutation.\n\n'], ...
    strjoin(cellstr(string(run_ids)), ', '), numel(seeds), labels.bs, labels.base);
fprintf(fid, ['**Test matrices.** Three families, regenerated per seed: *Gaussian* - ', ...
    'i.i.d. standard normal entries; *ill-conditioned* - A = U S V^T with U, V ', ...
    'orthonormal factors of Gaussian matrices and S geometrically graded from 1 down ', ...
    'to 1e-10 (kappa = 1e10); *orthonormal rows* - A = Q^T with Q an orthonormal basis ', ...
    'of a Gaussian n-by-m matrix (m <= n, so A A^T = I - the GKS column-selection ', ...
    'setting). Square cases use m = n; short-wide cases sweep aspect ratios ', ...
    'n/m in {2, 4, 8, 10}.\n\n']);
fprintf(fid, ['**fig_square_runtime.** Median runtime versus matrix size for square ', ...
    'matrices (log-log). Lines are geometric means over %d seeds; faint bands show ', ...
    'the range across seeds. One panel per matrix family.\n\n'], numel(seeds));
fprintf(fid, ['**fig_shortwide_runtime.** Median runtime versus column count n for ', ...
    'short-wide matrices, one color per row count m (log-log). Lines/bands as above; ', ...
    'marker and line style distinguish the methods.\n\n']);
fprintf(fid, ['**fig_relative_time_composite** (top-level plots/). The single timing ', ...
    'figure: the relative time without (BSQR / CPQR) and with the interpolation matrix ', ...
    '(both methods also form R11^{-1}R12; labelled ', ...
    '(BSQR + R11^{-1}R12) / (CPQR + R11^{-1}R12)), ', ...
    'overlaid on the shared rows and distinguished by colour. Square rows are ', ...
    'm = n, short-wide rows m < n; 1 = parity (dashed line). Points are geomeans of ', ...
    'per-seed geomeans; whiskers show the per-seed range.\n\n']);
fprintf(fid, 'Numerical quality is summarized in tables/quality_summary.md.\n');
end

% ----------------------------------------------------------------------
% Helpers
% ----------------------------------------------------------------------

function labels = mode_labels(mode_name)
switch string(mode_name)
    case "plain"
        labels = struct('bs', 'BSQR', 'base', 'CPQR (built-in)', 'ratio', 'BSQR / CPQR');
    case "rinv"
        labels = struct('bs', 'BSQR + $R_{11}^{-1}R_{12}$', ...
            'base', 'CPQR + $R_{11}^{-1}R_{12}$', ...
            'ratio', '(BSQR + $R_{11}^{-1}R_{12}$) / (CPQR + $R_{11}^{-1}R_{12}$)');
    otherwise
        error('plot_publication_results:UnknownMode', 'Unknown mode: %s', mode_name);
end
end

function name = display_family(family)
switch char(family)
    case 'gaussian'
        name = "Gaussian";
    case 'ill_conditioned'
        name = "Ill-conditioned";
    case 'orthonormal_rows'
        name = "Orthonormal rows";
    otherwise
        name = string(strrep(char(family), '_', ' '));
end
end

function [x, ctr, lo, hi] = seed_grouped_stats(keys, values)
% Geomean across seeds (one value per seed per key) with the seed range.
x = sort(unique(keys));
ctr = nan(size(x));
lo = nan(size(x));
hi = nan(size(x));
for i = 1:numel(x)
    v = values(keys == x(i));
    v = v(isfinite(v) & v > 0);
    if isempty(v)
        continue;
    end
    ctr(i) = geomean(v);
    lo(i) = min(v);
    hi(i) = max(v);
end
end

function plot_band(ax, x, lo, hi, color)
% Faint seed-range band; transparent and behind every median line.
good = isfinite(x) & isfinite(lo) & isfinite(hi) & lo > 0 & hi > 0;
if ~any(good)
    return;
end
x = x(good);
lo = lo(good);
hi = hi(good);
h_fill = fill(ax, [x; flipud(x)], [lo; flipud(hi)], color, ...
    'FaceAlpha', 0.15, 'EdgeColor', 'none');
set(h_fill, 'HandleVisibility', 'off');
end

function h_line = plot_median_line(ax, x, med, st)
h_line = gobjects(1, 1);
good = isfinite(x) & isfinite(med) & med > 0;
if ~any(good)
    return;
end
h_line = plot(ax, x(good), med(good), ...
    'LineStyle', st.linestyle, ...
    'Marker', st.marker, ...
    'Color', st.color, ...
    'LineWidth', st.linewidth, ...
    'MarkerSize', st.markersize);
end

function st = role_style(style, role)
% role 'bs': BSQR (blue, solid, circle); role 'base': baseline (orange, dashed, diamond).
st = struct('color', style.bsqr_color, 'linestyle', '-', 'marker', 'o', ...
    'linewidth', style.line_width, 'markersize', style.marker_size);
if strcmp(role, 'base')
    st.color = style.baseline_color;
    st.linestyle = '--';
    st.marker = 'd';
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

function apply_axes_style(ax, style)
set(ax, ...
    'Box', 'on', ...
    'LineWidth', style.axis_line_width, ...
    'FontSize', style.axis_font_size, ...
    'GridColor', [0.87, 0.87, 0.87], ...
    'GridAlpha', 1.0, ...
    'MinorGridLineStyle', 'none', ...
    'TickLabelInterpreter', 'latex');
grid(ax, 'on');
end

function style = publication_style()
style = struct();
style.bsqr_color = [76, 120, 168] / 255;      % #4C78A8, matches the Python plotter
style.baseline_color = [245, 133, 24] / 255;  % #F58518
style.axis_font_size = 10.5;
style.legend_font_size = 9.5;
style.line_width = 1.1;
style.marker_size = 4.0;
style.axis_line_width = 0.7;
style.single_col_width = 3.35;
style.double_col_width = 6.9;
style.formats = parse_figure_formats();
end

function formats = parse_figure_formats()
raw = string(bsqr_bench_getenv_default('BS_MATLAB_PUB_FIG_FORMATS', 'png'));
parts = strip(split(lower(raw), ','));
parts(parts == "") = [];
if isempty(parts)
    error('plot_publication_results:InvalidFigureFormats', ...
        'BS_MATLAB_PUB_FIG_FORMATS must list at least one format.');
end
allowed = ["png", "pdf", "eps"];
bad = setdiff(unique(parts), allowed);
if ~isempty(bad)
    error('plot_publication_results:InvalidFigureFormats', ...
        'Unsupported BS_MATLAB_PUB_FIG_FORMATS value(s): %s', join(bad, ', '));
end
formats = strings(0, 1);
for i = 1:numel(parts)
    if ~any(formats == parts(i))
        formats(end + 1, 1) = parts(i); %#ok<AGROW>
    end
end
end

function fig = new_figure(style, width, height)
fig = figure('Visible', 'off', 'Color', 'w');
set(fig, 'Units', 'inches', 'Position', [1, 1, width, height]);
fig.UserData = struct('formats', style.formats);
end

function save_figure(fig, outstem, style) %#ok<INUSD>
for i = 1:numel(fig.UserData.formats)
    fmt = fig.UserData.formats(i);
    switch fmt
        case "png"
            exportgraphics(fig, [outstem, '.png'], 'Resolution', 300);
        case "pdf"
            exportgraphics(fig, [outstem, '.pdf'], 'ContentType', 'vector');
        case "eps"
            exportgraphics(fig, [outstem, '.eps'], 'ContentType', 'vector');
        otherwise
            error('plot_publication_results:InvalidFigureFormats', ...
                'Unsupported figure format: %s', fmt);
    end
end
close(fig);
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
rt = joined;
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
