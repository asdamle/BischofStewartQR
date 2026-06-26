function plot_largen_scaling(varargin)
%PLOT_LARGEN_SCALING Publication figure from run_largen_scaling.
%
%   Reads results/exp_largen.csv and writes fig_largen_scaling.{png,pdf}: runtime
%   vs n on log-log axes, one panel per family (gaussian | needle), with three
%   sampler curves -- randomized BSQR (norm-weighted), randomized BSQR (uniform),
%   and rejection_rpqr. Line = seed median; shaded band = seed min/max.
%
%   The intended reading: on gaussian all three scale the same way (uniform is
%   fastest, lacking the O(mn) norm precompute); on needle uniform sampling
%   degenerates and its curve climbs away from the two norm-weighted methods,
%   which stay together.
%
% Options: 'resultsdir', 'plotdir', 'formats' (default {'png','pdf'} -- vector PDF
%   for publication alongside PNG; pass {'png'} to skip the PDF for quick iteration).

ip = inputParser;
addParameter(ip, 'resultsdir', '');
addParameter(ip, 'plotdir', '');
addParameter(ip, 'formats', {'png', 'pdf'});
parse(ip, varargin{:});
opt = ip.Results;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if isempty(opt.resultsdir); opt.resultsdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'results'); end
if isempty(opt.plotdir);    opt.plotdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'plots'); end
if ~isfolder(opt.plotdir); mkdir(opt.plotdir); end

f = fullfile(opt.resultsdir, 'exp_largen.csv');
if ~isfile(f); warning('missing %s', f); return; end
T = readtable(f, 'TextType', 'string');
fams = unique(T.family, 'stable');

% method key -> (legend label, color, marker)
M = { 'bsqr_nw',        'randomized BSQR (norm-weighted)', [0.00 0.45 0.74], 'o'; ...
      'bsqr_unif',      'randomized BSQR (uniform)',       [0.00 0.60 0.30], 's'; ...
      'rejection_rpqr', 'rejection\_rpqr',                 [0.85 0.33 0.10], 'd' };

nf = numel(fams);
fig = figure('Position', [100 100 max(460 * nf, 760) 460], 'Color', 'w');
tl = tiledlayout(fig, 1, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = [];
for fi = 1:nf
    fam = fams(fi); base = T.family == fam;
    ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log', 'YScale', 'log');
    for mi = 1:size(M, 1)
        band_line(ax, T, base & T.method == M{mi, 1}, M{mi, 3}, M{mi, 4}, M{mi, 2});
    end
    grid(ax, 'on'); xlabel(ax, 'n'); ylabel(ax, 'time (s)');
    title(ax, fam, 'Interpreter', 'none');
end
lg = legend(ax, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
title(tl, sprintf('Runtime vs n at fixed k = %d', T.k(1)));
save_fig(fig, fullfile(opt.plotdir, 'fig_largen_scaling'), opt.formats);
end

% ===========================================================================
function band_line(ax, T, mask, color, marker, name)
if ~any(mask); return; end
sub = T(mask, :); x = unique(sub.n); ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(x)
    v = sub.time_s(sub.n == x(i)); v = v(isfinite(v) & v > 0);
    ym(i) = median(v); ylo(i) = min(v); yhi(i) = max(v);
end
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, ['-' marker], 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
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
