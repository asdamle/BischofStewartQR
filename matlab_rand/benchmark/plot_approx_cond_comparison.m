function plot_approx_cond_comparison(varargin)
%PLOT_APPROX_COND_COMPARISON Figure for the R11-conditioning companion.
%
%   Reads results/exp_approx_cond.csv (from run_approx_cond_comparison) and writes
%   a 4-row figure (one tiled column per family) to benchmark/plots/:
%     fig_approx_cond_quality.{png,pdf}
%       row 1: ||R11^{-1}||_F / bound   (BSQR guaranteed <= 1; rejection uncontrolled)
%       row 2: max |R11^{-1} R12|       (interpolation-coefficient magnitude)
%       row 3: rank-k ID error / ||A||_F (noiseless -- oblique-coefficient penalty,
%              amplified by ||R11^{-1}|| above the dashed projection optimum)
%       row 4: noisy ID error / ||A||_F (row 3 + measurement noise propagated
%              through T ~ ||R11^{-1}||), same dashed projection lower bound
%   The dashed line (orthogonal-projection error) is the method-independent lower
%   bound for any reconstruction from the selected columns. Lines = seed mean;
%   bands = seed min/max. All y-axes log.
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
methods = {'bsqr_rand', 'randomized BSQR', C.bsqr; ...
           'rejection_rpqr', 'rejection\_rpqr', C.rpqr};
% Four rows: the cause (conditioning), the mechanism (coefficients), then the
% rank-k ID error split into its noiseless part (oblique-coefficient penalty,
% conditioning-amplified) and its noisy part (measurement noise propagated through
% T). The orthogonal-projection error is the dashed, method-independent lower bound
% on both ID rows. The 'norm' option swaps rows 1/3/4 (and the reference) between
% Frobenius and spectral; row 2 (max|R11^{-1}R12|) is max-norm in both.
is2 = any(strcmpi(opt.norm, {'2', 'spec', 'spectral'}));
if is2
    need = {'specinv', 'osinsky2', 'proj_spec', 'id_spec', 'noisy_id_spec'};
    if ~all(ismember(need, T.Properties.VariableNames))
        warning('CSV lacks spectral columns (%s); re-run run_approx_cond_comparison.', ...
            strjoin(need, ', ')); return;
    end
    T.condratio = T.specinv ./ T.osinsky2;
    stem = 'fig_approx_cond_quality_spec';
    rowspec = {'condratio',     '||R_{11}^{-1}||_2 / bound',          '',                   ''; ...
               'maxT',          'max |R_{11}^{-1}R_{12}|',            '',                   ''; ...
               'id_spec',       'rank-k ID error (2-norm) / ||A||_2', 'proj_spec',          'svdopt_spec'; ...
               'noisy_id_spec', 'noisy ID error (2-norm) / ||A||_2',  'noisy_id_spec_proj', 'svdopt_spec'};
    normlabel = 'spectral';
else
    T.condratio = T.frobinv ./ T.osinsky;
    stem = 'fig_approx_cond_quality';
    rowspec = {'condratio',    '||R_{11}^{-1}||_F / bound',          '',                  ''; ...
               'maxT',         'max |R_{11}^{-1}R_{12}|',            '',                  ''; ...
               'id_err',       'rank-k ID error / ||A||_F',          'proj_frob',         'svdopt_frob'; ...
               'noisy_id_err', 'noisy ID error / ||A||_F',           'noisy_id_err_proj', 'svdopt_frob'};
    normlabel = 'Frobenius';
end
nr = size(rowspec, 1);
nf = numel(fams);
fig = figure('Position', [100 100 max(440 * nf, 720) 300 * nr], 'Color', 'w');
tl = tiledlayout(fig, nr, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = [];
for ri = 1:nr
    col = rowspec{ri, 1};
    projcol = rowspec{ri, 3};        % T_proj overlay column ('' = none)
    svdoptcol = rowspec{ri, 4};      % SVD best-rank-k lower-bound column ('' = none)
    for fi = 1:nf
        fam = fams(fi); base = T.family == fam;
        ax = nexttile(tl); hold(ax, 'on');
        set(ax, 'XScale', opt.kscale, 'YScale', 'log', 'TickLabelInterpreter', 'none');
        if strcmp(col, 'condratio')
            yline(ax, 1, 'k--', 'HandleVisibility', 'off');     % the Osinsky bound
        end
        for mi = 1:size(methods, 1)                              % solid: V_k-frame T
            band_line(ax, T(base & T.method == methods{mi, 1}, :), col, ...
                methods{mi, 3}, methods{mi, 2});
        end
        if ~isempty(projcol)                                     % thin dashed: T_proj
            for mi = 1:size(methods, 1)
                proj_coeff_line(ax, T(base & T.method == methods{mi, 1}, :), projcol, ...
                    methods{mi, 3}, [methods{mi, 2}, ' (proj. coeffs)']);
            end
        end
        if ~isempty(svdoptcol)                                   % black dotted: SVD lower bound
            svdopt_line(ax, T(base, :), svdoptcol);
        end
        grid(ax, 'on');
        if ri == 1; title(ax, fam, 'Interpreter', 'none'); end
        if ri == nr; xlabel(ax, 'rank k'); end
        if fi == 1; ylabel(ax, rowspec{ri, 2}); end
    end
end
lg = legend(ax, 'Orientation', 'horizontal'); lg.Layout.Tile = 'south';
eps_str = '';
if ismember('noise_rel', T.Properties.VariableNames)
    eps_str = sprintf(' (noise \\epsilon=%.0e)', T.noise_rel(1));
end
title(tl, {['Conditioning, coefficient magnitude, and rank-k ID error (noiseless & noisy', ...
    eps_str, ') -- ', normlabel, ' norm'], ...
    'ID rows: solid = T (V_k-frame), dashed = T_{proj} (projection), dotted = best rank-k (SVD)'});
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

function proj_coeff_line(ax, sub, col, color, name)
% Thin dashed, mean-only line (per method) for the projection-coefficient (T_proj)
% variant of an ID error, overlaid on the solid V_k-frame (T) line so the two
% coefficient choices read off the same panel. On the noiseless ID row this is the
% orthogonal-projection error itself (P_S A = A(:,S)[I, T_proj]), i.e. that method's
% own conditioning-blind lower bound -- so the solid T line always sits above it.
if isempty(sub); return; end
[x, ym] = agg_by_k(sub, col);
plot(ax, x, ym, '--', 'Color', color, 'LineWidth', 1.0, 'Marker', 'none', ...
    'DisplayName', name);
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
% Center = seed mean; band = seed min/max, consistent with the other rand figures
% (and with 20 trials a 5/95 band would barely differ).
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
