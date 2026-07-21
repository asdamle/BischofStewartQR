function plot_general_largen(varargin)
%PLOT_GENERAL_LARGEN Plot for the large-n qr_builtin vs general randBSQR sweep.
%
%   Reads results/exp_general_largen.csv (from run_general_largen) and writes
%   fig_general_largen to benchmark/plots/: three panels vs n --
%     time      - both methods, plus the wrapper's qr(A') / selection split
%     speedup   - t_builtin / t_general (> 1 means the wrapper is faster)
%     quality   - sigma_min(A(:,S))/sigma_min(A) with the guarantee line
%   Line = seed mean; shaded band = seed min/max.
%
% Options: 'resultsdir', 'plotdir', 'formats' (default {'png','pdf'}).

ip = inputParser;
addParameter(ip, 'resultsdir', '');
addParameter(ip, 'plotdir', '');
addParameter(ip, 'formats', {'png', 'pdf'});
parse(ip, varargin{:});
opt = ip.Results;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if isempty(opt.resultsdir)
    opt.resultsdir = fullfile(repo_root, 'matlab_general', 'benchmark', 'results');
end
if isempty(opt.plotdir)
    opt.plotdir = fullfile(repo_root, 'matlab_general', 'benchmark', 'plots');
end
if ~isfolder(opt.plotdir); mkdir(opt.plotdir); end

f = fullfile(opt.resultsdir, 'exp_general_largen.csv');
if ~isfile(f); warning('missing %s', f); return; end
T = readtable(f, 'TextType', 'string');
m = T.m(1); fam = T.family(1);
cq = [0.40 0.40 0.40]; cg = [0.49 0.18 0.56];

fig = figure('Position', [100 100 1180 420], 'Color', 'w');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log', 'YScale', 'log');
qb = T.method == "qr_builtin"; gr = T.method == "general_bsqr_rand";
band_line(ax, T, qb, 'time_s', cq, 'qr (built-in)');
band_line(ax, T, gr, 'time_s', cg, 'general: randBSQR');
comp_line(ax, T, gr, 't_qr', [0.85 0.33 0.10], 'wrapper: qr(A^T)');
comp_line(ax, T, gr, 't_select', [0.00 0.45 0.74], 'wrapper: select');
grid(ax, 'on'); xlabel(ax, 'n'); ylabel(ax, 'time (s)');
lg = legend(ax, 'Location', 'northwest'); lg.FontSize = lg.FontSize - 1;

ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log');
speedup_line(ax, T, cg);
yline(ax, 1, 'k--', 'HandleVisibility', 'off');
grid(ax, 'on'); xlabel(ax, 'n'); ylabel(ax, 'speedup  (t_{builtin} / t_{general})');

ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log', 'YScale', 'log');
band_line(ax, T, qb, 'ratio', cq, 'qr (built-in)');
band_line(ax, T, gr, 'ratio', cg, 'general: randBSQR');
xs = unique(T.n);
plot(ax, xs, 1 ./ sqrt(double(m) * (xs - double(m) + 1)), 'k--', ...
    'LineWidth', 1.2, 'DisplayName', 'guarantee');
grid(ax, 'on'); xlabel(ax, 'n'); ylabel(ax, '\sigma_{min}(A(:,S)) / \sigma_{min}(A)');

title(tl, sprintf('Built-in pivoted QR vs general randBSQR   (%s, m = k = %d)', fam, m), ...
    'Interpreter', 'none');
save_fig(fig, fullfile(opt.plotdir, 'fig_general_largen'), opt.formats);
end

% ===========================================================================
function band_line(ax, T, mask, yname, color, name)
sub = T(mask, :); x = unique(sub.n); ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(x)
    v = sub.(yname)(sub.n == x(i)); v = v(isfinite(v));
    ym(i) = mean(v); ylo(i) = min(v); yhi(i) = max(v);
end
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
end

function comp_line(ax, T, mask, yname, color, name)
sub = T(mask, :); x = unique(sub.n); ym = zeros(size(x));
for i = 1:numel(x)
    v = sub.(yname)(sub.n == x(i)); v = v(isfinite(v));
    ym(i) = mean(v);
end
plot(ax, x, ym, ':', 'Color', color, 'LineWidth', 1.5, 'DisplayName', name);
end

function speedup_line(ax, T, color)
x = unique(T.n); ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(x)
    sp = [];
    for s = unique(T.seed)'
        tb = T.time_s(T.method == "qr_builtin" & T.n == x(i) & T.seed == s);
        tg = T.time_s(T.method == "general_bsqr_rand" & T.n == x(i) & T.seed == s);
        if ~isempty(tb) && ~isempty(tg); sp(end+1) = tb(1) / tg(1); end %#ok<AGROW>
    end
    ym(i) = mean(sp); ylo(i) = min(sp); yhi(i) = max(sp);
end
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4);
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
