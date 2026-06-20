function plot_rpqr_comparison(varargin)
%PLOT_RPQR_COMPARISON Plots for the bsqr_rand vs rejection_rpqr comparison.
%
%   Reads results/exp_rpqr_k<K>.csv (from run_rpqr_comparison) and writes three
%   single-metric figures to benchmark/plots/:
%     fig_rpqr_time_k<K>.png    - runtime vs n (both methods, per family)
%     fig_rpqr_speedup_k<K>.png - speedup t_rpqr / t_bsqr vs n
%     fig_rpqr_quality_k<K>.png - ||R11^{-1}||_F / Osinsky bound vs n
%   Line = seed mean; shaded band = seed min/max.
%
% Options: 'k' (default 64), 'resultsdir', 'plotdir', 'formats'.

ip = inputParser;
addParameter(ip, 'k', 64);
addParameter(ip, 'resultsdir', '');
addParameter(ip, 'plotdir', '');
addParameter(ip, 'formats', {'png'});
parse(ip, varargin{:});
opt = ip.Results;
opt.tag = sprintf('_k%d', opt.k);

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if isempty(opt.resultsdir)
    opt.resultsdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'results');
end
if isempty(opt.plotdir)
    opt.plotdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'plots');
end
if ~isfolder(opt.plotdir); mkdir(opt.plotdir); end

f = fullfile(opt.resultsdir, ['exp_rpqr' opt.tag '.csv']);
if ~isfile(f); warning('missing %s', f); return; end
T = readtable(f, 'TextType', 'string');
fams = unique(T.family, 'stable');
C = struct('bsqr', [0.00 0.45 0.74], 'rpqr', [0.85 0.33 0.10]);

rpqr_figure(opt, C, T, fams, 'time');
rpqr_figure(opt, C, T, fams, 'speedup');
rpqr_figure(opt, C, T, fams, 'quality');
end

% ===========================================================================
function rpqr_figure(opt, C, T, fams, which)
nf = numel(fams);
fig = figure('Position', [100 100 max(440 * nf, 720) 440], 'Color', 'w');
tl = tiledlayout(fig, 1, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = [];
for fi = 1:nf
    fam = fams(fi); base = T.family == fam;
    ax = nexttile(tl); hold(ax, 'on');
    switch which
        case 'time'
            set(ax, 'XScale', 'log', 'YScale', 'log');
            band_line(ax, T, base & T.method == "bsqr", 'time_s', C.bsqr, 'randomized BSQR');
            band_line(ax, T, base & T.method == "rejection_rpqr", 'time_s', C.rpqr, 'rejection\_rpqr');
            ylabel(ax, 'time (s)');
        case 'speedup'
            set(ax, 'XScale', 'log', 'YScale', 'log');
            speedup_line(ax, T, base, C.bsqr);
            yline(ax, 1, 'k--', 'HandleVisibility', 'off');
            ylabel(ax, 'speedup (t_{rejection\_rpqr} / t_{BSQR})');
        case 'quality'
            set(ax, 'XScale', 'log');
            cond_line(ax, T, base & T.method == "bsqr", C.bsqr, 'randomized BSQR');
            cond_line(ax, T, base & T.method == "rejection_rpqr", C.rpqr, 'rejection\_rpqr');
            yline(ax, 1, 'k--', 'HandleVisibility', 'off');   % Osinsky bound
            ylabel(ax, '||R_{11}^{-1}||_F / Osinsky bound');
    end
    grid(ax, 'on'); xlabel(ax, 'n'); title(ax, fam, 'Interpreter', 'none');
end
if ~strcmp(which, 'speedup')
    lg = legend(ax, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
end
heads = struct('time', 'Runtime vs n', ...
    'speedup', 'Speedup of randomized BSQR over rejection\_rpqr (t_{rpqr}/t_{BSQR})', ...
    'quality', 'Selection quality: ||R_{11}^{-1}||_F relative to the Osinsky bound (lower = better)');
title(tl, sprintf('%s   (k=%d)', heads.(which), opt.k));
save_fig(fig, fullfile(opt.plotdir, ['fig_rpqr_' which opt.tag]), opt.formats);
end

% ===========================================================================
function [x, ym, ylo, yhi] = agg(T, mask, yname)
sub = T(mask, :);
x = unique(sub.n);
ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(x)
    v = sub.(yname)(sub.n == x(i));
    v = v(isfinite(v));
    if isempty(v); ym(i) = NaN; ylo(i) = NaN; yhi(i) = NaN; continue; end
    ym(i) = mean(v); ylo(i) = min(v); yhi(i) = max(v);
end
end

function band_line(ax, T, mask, yname, color, name)
if ~any(mask); return; end
[x, ym, ylo, yhi] = agg(T, mask, yname);
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
end

function cond_line(ax, T, mask, color, name)
if ~any(mask); return; end
sub = T(mask, :); sub.ratio = sub.frobinv ./ sub.osinsky;
x = unique(sub.n); ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(x)
    v = sub.ratio(sub.n == x(i)); v = v(isfinite(v));
    ym(i) = mean(v); ylo(i) = min(v); yhi(i) = max(v);
end
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
end

function speedup_line(ax, T, base, color)
ns = unique(T.n(base));
x = ns; ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(ns)
    sp = [];
    for s = unique(T.seed(base))'
        tb = T.time_s(base & T.method == "bsqr" & T.n == ns(i) & T.seed == s);
        tr = T.time_s(base & T.method == "rejection_rpqr" & T.n == ns(i) & T.seed == s);
        if ~isempty(tb) && ~isempty(tr); sp(end+1) = tr(1) / tb(1); end %#ok<AGROW>
    end
    ym(i) = mean(sp); ylo(i) = min(sp); yhi(i) = max(sp);
end
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, 'LineWidth', 1.5, 'MarkerSize', 4);
end

function save_fig(fig, stem, formats)
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
