function plot_approx_cond_comparison(varargin)
%PLOT_APPROX_COND_COMPARISON Figure for the R11-conditioning companion.
%
%   Reads results/exp_approx_cond.csv (from run_approx_cond_comparison) and writes
%   a 4-row figure (one tiled column per family) to benchmark/plots/:
%     fig_approx_cond_quality.{png,pdf}
%       row 1: ||R11^{-1}||_F / bound   (BSQR guaranteed <= 1; rejection uncontrolled)
%       row 2: max |R11^{-1} R12|       (interpolation-coefficient magnitude)
%       row 3: rank-k ID error / ||A||_F (oblique-coefficient penalty, amplified
%              by ||R11^{-1}||)
%       row 4: orthogonal-projection error / ||A||_F (conditioning-blind best fit
%              in span(A(:,S)); near-identical per method)
%   The dotted black line on rows 3-4 is the best rank-k error (SVD truncation),
%   the absolute lower bound below any column selection. Lines = seed median;
%   bands = seed min/max. All y-axes log. (A noisy-ID variant and the projection
%   coefficients T_proj are recorded in the CSV but not plotted.)
%
% Options: 'norm' ('fro' default, or '2'/'spec' for the spectral version -- rows
%   1/3/4 become ||R11^{-1}||_2 and the 2-norm approximation errors; row 2 stays
%   max-norm; written as fig_approx_cond_quality_spec), 'resultsdir', 'plotdir',
%   'formats' (default {'png','pdf'}).

ip = inputParser;
addParameter(ip, 'norm', 'fro');
addParameter(ip, 'resultsdir', '');
addParameter(ip, 'plotdir', '');
addParameter(ip, 'formats', {'png', 'pdf'});
addParameter(ip, 'kscale', 'linear');  % x-axis (rank k): 'linear' or 'log'
parse(ip, varargin{:});
opt = ip.Results;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
if isempty(opt.resultsdir)
    opt.resultsdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'results');
end
if isempty(opt.plotdir)
    opt.plotdir = fullfile(repo_root, 'matlab_rand', 'benchmark', 'plots');
end
if ~isfolder(opt.plotdir); mkdir(opt.plotdir); end

f = fullfile(opt.resultsdir, 'exp_approx_cond.csv');
if ~isfile(f); warning('missing %s', f); return; end
T = readtable(f, 'TextType', 'string');
fams = unique(T.family, 'stable');
C = struct('bsqr', [0.00 0.45 0.74], 'rpqr', [0.85 0.33 0.10]);
methods = {'bsqr_rand', 'randBSQR', C.bsqr; ...
           'rejection_rpqr', 'rejection\_rpqr', C.rpqr};
% Four rows: the cause (conditioning), the mechanism (coefficients), the rank-k
% ID error (oblique V_k-frame coefficients, conditioning-amplified), and the
% orthogonal-projection error (the conditioning-blind reference). The dotted
% black SVD line on rows 3-4 is the absolute lower bound. The 'norm' option
% swaps rows 1/3/4 (and the reference) between Frobenius and spectral; row 2
% (max|R11^{-1}R12|) is max-norm in both.
is2 = any(strcmpi(opt.norm, {'2', 'spec', 'spectral'}));
if is2
    need = {'specinv', 'osinsky2', 'proj_spec', 'id_spec', 'svdopt_spec'};
    if ~all(ismember(need, T.Properties.VariableNames))
        warning('CSV lacks spectral columns (%s); re-run run_approx_cond_comparison.', ...
            strjoin(need, ', ')); return;
    end
    T.condratio = T.specinv ./ T.osinsky2;
    stem = 'fig_approx_cond_quality_spec';
    rowspec = {'condratio',  '||R_{11}^{-1}||_2 / bound',      ''; ...
               'maxT',       'max |R_{11}^{-1}R_{12}|',        ''; ...
               'id_spec',    'ID-error / ||A||_2',       'svdopt_spec'; ...
               'proj_spec',  'proj. error / ||A||_2',    'svdopt_spec'};
    normlabel = 'spectral';
else
    T.condratio = T.frobinv ./ T.osinsky;
    stem = 'fig_approx_cond_quality';
    rowspec = {'condratio',  '||R_{11}^{-1}||_F / bound',      ''; ...
               'maxT',       'max |R_{11}^{-1}R_{12}|',        ''; ...
               'id_err',     'ID-error / ||A||_F',       'svdopt_frob'; ...
               'proj_frob',  'proj. error / ||A||_F',    'svdopt_frob'};
    normlabel = 'Frobenius';
end
nr = size(rowspec, 1);
nf = numel(fams);
fig = figure('Position', [100 100 max(440 * nf, 720) 300 * nr], 'Color', 'w');
tl = tiledlayout(fig, nr, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = [];
for ri = 1:nr
    col = rowspec{ri, 1};
    svdoptcol = rowspec{ri, 3};      % SVD best-rank-k lower-bound column ('' = none)
    for fi = 1:nf
        fam = fams(fi); base = T.family == fam;
        ax = nexttile(tl); hold(ax, 'on');
        set(ax, 'XScale', opt.kscale, 'YScale', 'log', 'TickLabelInterpreter', 'none');
        if strcmp(col, 'condratio')
            yline(ax, 1, 'k--', 'HandleVisibility', 'off');     % the Osinsky bound
        end
        for mi = 1:size(methods, 1)                              % solid per-method curve
            band_line(ax, T(base & T.method == methods{mi, 1}, :), col, ...
                methods{mi, 3}, methods{mi, 2});
        end
        if ~isempty(svdoptcol)                                   % black dotted: SVD lower bound
            svdopt_line(ax, T(base, :), svdoptcol);
        end
        grid(ax, 'on');
        if ri == 1; title(ax, fam, 'Interpreter', 'none'); end
        if ri == nr; xlabel(ax, 'k'); end
        if fi == 1; ylabel(ax, rowspec{ri, 2}); end
    end
end
lg = legend(ax, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
title(tl, 'Matrix approximation quality');
save_fig(fig, fullfile(opt.plotdir, stem), opt.formats);
end

% ===========================================================================
function band_line(ax, sub, col, color, name)
if isempty(sub); return; end
[x, ym, ylo, yhi] = agg_by_k(sub, col);
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
% Lines only (no per-point markers): with the dense k grid markers read as
% clutter/jaggedness; method is distinguished by colour.
plot(ax, x, ym, '-', 'Color', color, 'LineWidth', 1.5, 'DisplayName', name);
end

function svdopt_line(ax, sub, col)
% Best-possible rank-k error (SVD truncation): the absolute lower bound, below any
% column-selection projection or ID error. Method-independent (constant per k),
% black dotted.
if isempty(sub); return; end
[x, ym] = agg_by_k(sub, col);
plot(ax, x, ym, 'k:', 'LineWidth', 1.0, 'DisplayName', 'best rank-k (SVD)');
end

function [x, ym, ylo, yhi] = agg_by_k(sub, col)
% Center = seed median; band = seed min/max (with 20 trials a 5/95 band barely
% differs). Median is robust to the occasional rejection_rpqr outlier seed; the
% min/max band still shows the full spread.
x = unique(sub.k); ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(x)
    v = sub.(col)(sub.k == x(i)); v = v(isfinite(v));
    if isempty(v); v = NaN; end
    ym(i) = median(v); ylo(i) = min(v); yhi(i) = max(v);
end
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
% Enlarge all figure text by delta points: ticks, axis labels, and panel
% titles scale with the axes font; legends are bumped directly; the
% tiledlayout super-title is set to match the axis-label size.
labelsz = 0;
for ax = reshape(findall(fig, 'Type', 'axes'), 1, [])
    ax.FontSize = ax.FontSize + delta;   % ticks, labels, and panel title
    labelsz = max(labelsz, ax.FontSize * ax.LabelFontSizeMultiplier);
end
for lg = reshape(findall(fig, 'Type', 'legend'), 1, [])
    lg.FontSize = lg.FontSize + delta;
end
for tl = reshape(findall(fig, 'Type', 'tiledlayout'), 1, [])
    tl.Title.FontSize = labelsz;
end
end
