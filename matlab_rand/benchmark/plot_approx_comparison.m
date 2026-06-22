function plot_approx_comparison(varargin)
%PLOT_APPROX_COMPARISON Figure for the low-rank approximation comparison.
%
%   Reads results/exp_approx<tag>.csv (from run_approx_comparison) and writes a
%   figure to benchmark/plots/:
%     fig_approx<tag>_quality.{png,pdf} -- relative approximation error vs rank k,
%       one tiled column per family. Two accuracy rows, plus a third row when the CSV
%       carries the maxT column (the default runner writes it):
%         row 1 = Frobenius  ||A - P_S A||_F / ||A||_F
%         row 2 = spectral   ||A - P_S A||_2 / ||A||_2
%         row 3 = max |R_{11}^{-1}R_{12}|   (interpolation-coefficient magnitude)
%       randomized BSQR and rejection_rpqr as lines (seed mean) with min/max bands;
%       the accuracy rows add the optimal rank-k error as a dashed lower bound.
%
% Options: 'tag' (default '' -- the application run; pass '_synth' for the
%   synthetic-spectrum companion), 'resultsdir', 'plotdir', 'formats' (default
%   {'png','pdf'} -- vector PDF for publication alongside PNG; pass {'png'} to skip).

ip = inputParser;
addParameter(ip, 'tag', '');
addParameter(ip, 'resultsdir', '');
addParameter(ip, 'plotdir', '');
addParameter(ip, 'formats', {'png', 'pdf'});
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

f = fullfile(opt.resultsdir, ['exp_approx' opt.tag '.csv']);
if ~isfile(f); warning('missing %s', f); return; end
T = readtable(f, 'TextType', 'string');
fams = unique(T.family, 'stable');
C = struct('bsqr', [0.00 0.45 0.74], 'rpqr', [0.85 0.33 0.10]);
methods = {'bsqr_rand', 'randomized BSQR', C.bsqr; ...
           'rejection_rpqr', 'rejection\_rpqr', C.rpqr};

% Accuracy rows always present; the interpolation-coefficient row is added when
% the CSV carries maxT (newer runs), so the figure pairs accuracy with the
% basis-quality metric on the same matrices. Older CSVs still plot (2 rows).
rowdefs = {'frob', '||A - P_SA||_F / ||A||_F'; ...
           'spec', '||A - P_SA||_2 / ||A||_2'};
has_maxT = ismember('maxT', T.Properties.VariableNames);
if has_maxT
    rowdefs(end+1, :) = {'maxT', 'max |R_{11}^{-1}R_{12}|'};
end
nr = size(rowdefs, 1);

nf = numel(fams);
fig = figure('Position', [100 100 max(440 * nf, 720) 360 * nr], 'Color', 'w');
tl = tiledlayout(fig, nr, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = [];
for ri = 1:nr
    which = rowdefs{ri, 1};
    for fi = 1:nf
        fam = fams(fi); base = T.family == fam;
        ax = nexttile(tl); hold(ax, 'on');
        set(ax, 'XScale', 'log', 'YScale', 'log', 'TickLabelInterpreter', 'none');
        if strcmp(which, 'maxT')
            for mi = 1:size(methods, 1)
                metric_line(ax, T(base & T.method == methods{mi, 1}, :), 'maxT', ...
                    methods{mi, 3}, methods{mi, 2});
            end
        else
            opt_ref_line(ax, T(base, :), which);         % optimal rank-k (lower bound)
            for mi = 1:size(methods, 1)
                relerr_line(ax, T(base & T.method == methods{mi, 1}, :), ...
                    which, methods{mi, 3}, methods{mi, 2});
            end
        end
        grid(ax, 'on');
        if ri == 1; title(ax, fam, 'Interpreter', 'none'); end
        if ri == nr; xlabel(ax, 'rank k'); end
        if fi == 1; ylabel(ax, rowdefs{ri, 2}); end
    end
end
lg = legend(ax, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
if has_maxT
    title(tl, ['Approximation accuracy vs interpolation coefficients (vs rank k): ', ...
        'Frobenius / spectral error (top two, dashed = optimal rank-k), ', ...
        'max|R_{11}^{-1}R_{12}| (bottom)']);
else
    title(tl, ['Low-rank approximation error vs rank k (lower = better): ', ...
        'Frobenius (top), spectral (bottom); dashed = optimal rank-k']);
end
save_fig(fig, fullfile(opt.plotdir, ['fig_approx' opt.tag '_quality']), opt.formats);
end

% ===========================================================================
function relerr_line(ax, sub, which, color, name)
% Relative error per norm, aggregated over seeds at each k:
%   'frob' -> frob_err / normA_fro;  'spec' -> spec_err / normA_2.
if isempty(sub); return; end
if strcmp(which, 'frob')
    sub.r = sub.frob_err ./ sub.normA_fro;
else
    sub.r = sub.spec_err ./ sub.normA_2;
end
[x, ym, ylo, yhi] = agg_by_k(sub, 'r');
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
end

function metric_line(ax, sub, col, color, name)
% Band line for a raw (already-absolute) metric column, aggregated over seeds.
if isempty(sub); return; end
[x, ym, ylo, yhi] = agg_by_k(sub, col);
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
end

function opt_ref_line(ax, sub, which)
% Optimal rank-k error (same for both methods); one value per k, no band.
if isempty(sub); return; end
if strcmp(which, 'frob')
    sub.r = sub.opt_frob ./ sub.normA_fro;
else
    sub.r = sub.opt_spec ./ sub.normA_2;
end
[x, ym] = agg_by_k(sub, 'r');
plot(ax, x, ym, 'k--', 'LineWidth', 1.0, 'DisplayName', 'optimal rank-k');
end

function [x, ym, ylo, yhi] = agg_by_k(sub, col)
x = unique(sub.k); ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(x)
    v = sub.(col)(sub.k == x(i)); v = v(isfinite(v));
    if isempty(v); v = NaN; end
    ym(i) = mean(v); ylo(i) = min(v); yhi(i) = max(v);
end
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
