function plot_general_comparison(varargin)
%PLOT_GENERAL_COMPARISON Plots for the general-A comparison suite.
%
%   Reads results/exp_general.csv (from run_general_comparison) and writes
%   four figures to benchmark/plots/:
%     fig_general_time.png    - runtime vs n (all five methods, per family)
%     fig_general_quality.png - sigma_min(A(:,S))/sigma_min(A) vs n, with the
%                               1/sqrt(m(n-m+1)) guarantee as a dashed line
%                               (the general_* methods must sit above it)
%     fig_general_phase.png   - wrapper phase split vs n: qr(A','econ') cost,
%                               per-selector selection cost, and total
%     fig_general_msweep.png  - runtime and quality vs m at fixed n
%   Line = seed mean; shaded band = seed min/max.
%
% Options: 'k' (default [] -> untagged CSV), 'm_fixed' (default 64),
%   'n_fixed' (default 4000), 'resultsdir', 'plotdir', 'formats' (default
%   {'png','pdf'}; pass {'png'} for quick iteration).

ip = inputParser;
addParameter(ip, 'k', []);
addParameter(ip, 'm_fixed', 64);
addParameter(ip, 'n_fixed', 4000);
addParameter(ip, 'resultsdir', '');
addParameter(ip, 'plotdir', '');
addParameter(ip, 'formats', {'png', 'pdf'});
parse(ip, varargin{:});
opt = ip.Results;
if isempty(opt.k); opt.tag = ''; else; opt.tag = sprintf('_k%d', opt.k); end

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if isempty(opt.resultsdir)
    opt.resultsdir = fullfile(repo_root, 'matlab_general', 'benchmark', 'results');
end
if isempty(opt.plotdir)
    opt.plotdir = fullfile(repo_root, 'matlab_general', 'benchmark', 'plots');
end
if ~isfolder(opt.plotdir); mkdir(opt.plotdir); end

f = fullfile(opt.resultsdir, ['exp_general' opt.tag '.csv']);
if ~isfile(f); warning('missing %s', f); return; end
T = readtable(f, 'TextType', 'string');
fams = unique(T.family, 'stable');

% method -> {display name, color}
M = { ...
    'qr_builtin',        'qr (built-in)',      [0.40 0.40 0.40]; ...
    'bsqr_direct',       'BSQR on A',          [0.00 0.45 0.74]; ...
    'bsqr_rand_direct',  'randBSQR on A',      [0.30 0.75 0.93]; ...
    'general_bsqr',      'general: BSQR',      [0.85 0.33 0.10]; ...
    'general_bsqr_rand', 'general: randBSQR',  [0.49 0.18 0.56]};

Tn = T(T.m == opt.m_fixed, :);   % the n-sweep
if isempty(Tn)
    warning('no rows with m = %d; skipping the n-sweep figures', opt.m_fixed);
else
    sweep_figure(opt, M, Tn, fams, 'time');
    sweep_figure(opt, M, Tn, fams, 'quality');
    phase_figure(opt, M, Tn, fams);
end
Tm = T(T.n == opt.n_fixed, :);   % the m-sweep
if numel(unique(Tm.m)) > 1
    msweep_figure(opt, M, Tm, fams);
end
end

% ===========================================================================
function sweep_figure(opt, M, T, fams, which)
nf = numel(fams);
fig = figure('Position', [100 100 max(360 * nf, 720) 440], 'Color', 'w');
tl = tiledlayout(fig, 1, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = [];
for fi = 1:nf
    fam = fams(fi); base = T.family == fam;
    ax = nexttile(tl); hold(ax, 'on');
    set(ax, 'XScale', 'log', 'YScale', 'log');
    switch which
        case 'time'
            for mi = 1:size(M, 1)
                band_line(ax, T, base & T.method == M{mi, 1}, T.n, 'time_s', M{mi, 3}, M{mi, 2});
            end
            if fi == 1; ylabel(ax, 'time (s)'); end
        case 'quality'
            for mi = 1:size(M, 1)
                band_line(ax, T, base & T.method == M{mi, 1}, T.n, 'ratio', M{mi, 3}, M{mi, 2});
            end
            % the k = m guarantee for the general_* methods, a curve in n
            xs = unique(T.n(base));
            plot(ax, xs, 1 ./ sqrt(double(opt.m_fixed) * (xs - opt.m_fixed + 1)), ...
                'k--', 'LineWidth', 1.2, 'DisplayName', 'guarantee');
            if fi == 1; ylabel(ax, '\sigma_{min}(A(:,S)) / \sigma_{min}(A)'); end
    end
    grid(ax, 'on'); xlabel(ax, 'n'); title(ax, fam, 'Interpreter', 'none');
end
lg = legend(ax, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
switch which
    case 'time';    title(tl, sprintf('Runtime vs n   (m = k = %d)', opt.m_fixed));
    case 'quality'; title(tl, sprintf('Subset quality vs n   (m = k = %d; higher is better)', opt.m_fixed));
end
save_fig(fig, fullfile(opt.plotdir, ['fig_general_' which opt.tag]), opt.formats);
end

% ===========================================================================
function phase_figure(opt, M, T, fams)
% Where the wrapper's time goes: the shared qr(A','econ') preprocessing, the
% per-selector selection phase, and the wrappers' end-to-end timeit totals.
nf = numel(fams);
fig = figure('Position', [100 100 max(360 * nf, 720) 440], 'Color', 'w');
tl = tiledlayout(fig, 1, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
cb = M{4, 3}; cr = M{5, 3};   % general_bsqr / general_bsqr_rand colors
ax = [];
for fi = 1:nf
    fam = fams(fi); base = T.family == fam;
    ax = nexttile(tl); hold(ax, 'on');
    set(ax, 'XScale', 'log', 'YScale', 'log');
    gb = base & T.method == "general_bsqr";
    gr = base & T.method == "general_bsqr_rand";
    band_line(ax, T, gb, T.n, 't_qr', [0.4 0.4 0.4], 'qr(A^T) (shared)');
    band_line(ax, T, gb, T.n, 't_select', cb, 'select: BSQR');
    band_line(ax, T, gr, T.n, 't_select', cr, 'select: randBSQR');
    dashed_line(ax, T, gb, T.n, 'time_s', cb, 'total: BSQR');
    dashed_line(ax, T, gr, T.n, 'time_s', cr, 'total: randBSQR');
    grid(ax, 'on'); xlabel(ax, 'n'); title(ax, fam, 'Interpreter', 'none');
    if fi == 1; ylabel(ax, 'time (s)'); end
end
lg = legend(ax, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
title(tl, sprintf('bsqr\\_general phase split vs n   (m = k = %d)', opt.m_fixed));
save_fig(fig, fullfile(opt.plotdir, ['fig_general_phase' opt.tag]), opt.formats);
end

% ===========================================================================
function msweep_figure(opt, M, T, fams)
% Two rows at fixed n: runtime vs m, and quality vs m.
nf = numel(fams);
fig = figure('Position', [100 100 max(360 * nf, 720) 620], 'Color', 'w');
tl = tiledlayout(fig, 2, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = [];
for row = {'time', 'quality'}
    for fi = 1:nf
        fam = fams(fi); base = T.family == fam;
        ax = nexttile(tl); hold(ax, 'on');
        set(ax, 'XScale', 'log', 'YScale', 'log');
        for mi = 1:size(M, 1)
            mask = base & T.method == M{mi, 1};
            if strcmp(row{1}, 'time')
                band_line(ax, T, mask, T.m, 'time_s', M{mi, 3}, M{mi, 2});
            else
                band_line(ax, T, mask, T.m, 'ratio', M{mi, 3}, M{mi, 2});
            end
        end
        grid(ax, 'on');
        if strcmp(row{1}, 'time')
            title(ax, fam, 'Interpreter', 'none');
            if fi == 1; ylabel(ax, 'time (s)'); end
        else
            xs = unique(T.m(base));
            plot(ax, xs, 1 ./ sqrt(xs .* (opt.n_fixed - xs + 1)), ...
                'k--', 'LineWidth', 1.2, 'DisplayName', 'guarantee');
            xlabel(ax, 'm');
            if fi == 1; ylabel(ax, '\sigma_{min}(A(:,S)) / \sigma_{min}(A)'); end
        end
    end
end
lg = legend(ax, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
title(tl, sprintf('Runtime and subset quality vs m   (n = %d, k = m)', opt.n_fixed));
save_fig(fig, fullfile(opt.plotdir, ['fig_general_msweep' opt.tag]), opt.formats);
end

% ===========================================================================
function band_line(ax, T, mask, xcol, yname, color, name)
% Mean line + min/max band across seeds, grouped by the x column values.
if ~any(mask); return; end
xs = xcol(mask); ys = T.(yname)(mask);
x = unique(xs); ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(x)
    v = ys(xs == x(i)); v = v(isfinite(v));
    if isempty(v); ym(i) = NaN; ylo(i) = NaN; yhi(i) = NaN; continue; end
    ym(i) = mean(v); ylo(i) = min(v); yhi(i) = max(v);
end
ok = isfinite(ym);
x = x(ok); ym = ym(ok); ylo = ylo(ok); yhi = yhi(ok);
if isempty(x); return; end
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
end

function dashed_line(ax, T, mask, xcol, yname, color, name)
% Seed-mean only, dashed (used for the wrapper totals in the phase figure).
if ~any(mask); return; end
xs = xcol(mask); ys = T.(yname)(mask);
x = unique(xs); ym = zeros(size(x));
for i = 1:numel(x)
    v = ys(xs == x(i)); v = v(isfinite(v));
    if isempty(v); ym(i) = NaN; else; ym(i) = mean(v); end
end
plot(ax, x, ym, '--', 'Color', color, 'LineWidth', 1.2, 'DisplayName', name);
end

function save_fig(fig, stem, formats)
bump_fonts(fig, 4);
for i = 1:numel(formats)
    switch formats{i}
        case 'png'; exportgraphics(fig, [stem, '.png'], 'Resolution', 150);
        case 'pdf'; exportgraphics(fig, [stem, '.pdf'], 'ContentType', 'vector');
        case 'eps'; exportgraphics(fig, [stem, '.eps'], 'ContentType', 'vector');
    end
end
fprintf('Wrote %s.%s\n', stem, strjoin(formats, ','));
close(fig);
end

function bump_fonts(fig, delta)
% Enlarge all figure text by delta points (same helper as the matlab_rand
% plot scripts): ticks/labels/titles scale with the axes font, legends are
% bumped directly, the tiledlayout super-title matches the axis-label size.
labelsz = 0;
for ax = reshape(findall(fig, 'Type', 'axes'), 1, [])
    ax.FontSize = ax.FontSize + delta;
    labelsz = max(labelsz, ax.FontSize * ax.LabelFontSizeMultiplier);
end
for lg = reshape(findall(fig, 'Type', 'legend'), 1, [])
    lg.FontSize = lg.FontSize + delta;
end
for tl = reshape(findall(fig, 'Type', 'tiledlayout'), 1, [])
    tl.Title.FontSize = labelsz;
end
end
