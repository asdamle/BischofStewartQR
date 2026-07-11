function plot_rpqr_comparison(varargin)
%PLOT_RPQR_COMPARISON Plots for the bsqr_rand vs rejection_rpqr comparison.
%
%   Reads results/exp_rpqr_k<K>.csv (from run_rpqr_comparison) and writes three
%   figures to benchmark/plots/:
%     fig_rpqr_time_k<K>.png    - runtime vs n (both methods, per family)
%     fig_rpqr_speedup_k<K>.png - speedup t_rpqr / t_bsqr vs n
%     fig_rpqr_quality_k<K>.png - two rows: ||R11^{-1}||_F / bound (<=1) and
%                                 ||R11^{-1}||_2 / bound (<=1), vs n
%   Line = seed mean; shaded band = seed min/max.
%
% Options: 'k' (default 64), 'resultsdir', 'plotdir', 'formats' (default
%   {'png','pdf'} -- vector PDF for publication alongside PNG; pass {'png'} to
%   skip the PDF for quick iteration).

ip = inputParser;
addParameter(ip, 'k', 64);
addParameter(ip, 'resultsdir', '');
addParameter(ip, 'plotdir', '');
addParameter(ip, 'formats', {'png', 'pdf'});
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
quality_figure(opt, C, T, fams);
end

% ===========================================================================
function rpqr_figure(opt, C, T, fams, which)
nf = numel(fams);
if strcmp(which, 'speedup')
    height = 300;   % rectangular (shorter than wide) panels compress the figure
else
    height = 440;
end
fig = figure('Position', [100 100 max(440 * nf, 720) height], 'Color', 'w');
tl = tiledlayout(fig, 1, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = [];
for fi = 1:nf
    fam = fams(fi); base = T.family == fam;
    ax = nexttile(tl); hold(ax, 'on');
    switch which
        case 'time'
            set(ax, 'XScale', 'log', 'YScale', 'log');
            band_line(ax, T, base & T.method == "bsqr", 'time_s', C.bsqr, 'randBSQR');
            band_line(ax, T, base & T.method == "rejection_rpqr", 'time_s', C.rpqr, 'rejection\_rpqr');
            if fi == 1; ylabel(ax, 'time (s)'); end
        case 'speedup'
            set(ax, 'XScale', 'log', 'YScale', 'log');
            speedup_line(ax, T, base, C.bsqr);
            yline(ax, 1, 'k--', 'HandleVisibility', 'off');
            if fi == 1; ylabel(ax, 'speedup (t_{rpqr} / t_{randBSQR})'); end
    end
    grid(ax, 'on'); xlabel(ax, 'n'); title(ax, fam, 'Interpreter', 'none');
end
if ~strcmp(which, 'speedup')
    lg = legend(ax, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
end
switch which
    case 'time';    title(tl, sprintf('Runtime vs n   (m=%d)', opt.k));
    case 'speedup'; title(tl, sprintf('Speedup of randBSQR over rejection\\_rpqr (m=%d)', opt.k));
end
save_fig(fig, fullfile(opt.plotdir, ['fig_rpqr_' which opt.tag]), opt.formats);
end

% ===========================================================================
function quality_figure(opt, C, T, fams)
% Two rows, both (inverse-norm / bound) <= 1 with lower better: the Frobenius
% ratio ||R11^{-1}||_F / bound and the spectral ratio ||R11^{-1}||_2 / bound.
nf = numel(fams);
fig = figure('Position', [100 100 max(440 * nf, 720) 480], 'Color', 'w');
tl = tiledlayout(fig, 2, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
methods = {'bsqr', 'randBSQR', C.bsqr; 'rejection_rpqr', 'rejection\_rpqr', C.rpqr};
ax = [];
for row = {'frob', 'smin'}
    for fi = 1:nf
        fam = fams(fi); base = T.family == fam;
        ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log');
        for mi = 1:size(methods, 1)
            ratio_line(ax, T, base & T.method == methods{mi, 1}, row{1}, methods{mi, 3}, methods{mi, 2});
        end
        yline(ax, 1, 'k--', 'HandleVisibility', 'off');   % the bound
        ylim(ax, [0 1.5]);   % consistent axis across families; rejection_rpqr bands clip here
        grid(ax, 'on');
        if strcmp(row{1}, 'frob')
            title(ax, fam, 'Interpreter', 'none');
            if fi == 1; ylabel(ax, '||R_{11}^{-1}||_F / bound'); end
        else
            xlabel(ax, 'n');
            if fi == 1; ylabel(ax, '||R_{11}^{-1}||_2 / bound'); end
        end
    end
end
lg = legend(ax, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
title(tl, sprintf('Selection quality vs bounds (m=%d)', opt.k));
save_fig(fig, fullfile(opt.plotdir, ['fig_rpqr_quality' opt.tag]), opt.formats);
end

% ===========================================================================
function ratio_line(ax, T, mask, which, color, name)
% Both rows are (inverse-norm / bound) <= 1, lower = better:
%   'frob' -> ||R11^{-1}||_F / sqrt(k(n-k+1));
%   'smin' -> ||R11^{-1}||_2 / sqrt(1+k(n-k)) = bound / sigma_min (lower better).
if ~any(mask); return; end
sub = T(mask, :);
if strcmp(which, 'frob')
    sub.r = sub.frobinv ./ sub.osinsky;
else
    sub.r = 1 ./ (sub.sigma_min .* sqrt(1 + sub.k .* (sub.n - sub.k)));
end
x = unique(sub.n); ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(x)
    v = sub.r(sub.n == x(i)); v = v(isfinite(v));
    ym(i) = mean(v); ylo(i) = min(v); yhi(i) = max(v);
end
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
end

function band_line(ax, T, mask, yname, color, name)
if ~any(mask); return; end
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
