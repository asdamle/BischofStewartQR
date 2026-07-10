function plot_rand_experiments(varargin)
%PLOT_RAND_EXPERIMENTS Publication plots from run_rand_experiments CSVs.
%
%   Reads the k-tagged exp_*_k<K>.csv from the results dir and writes k-tagged
%   figures fig_*_k<K>.png to matlab_rand/benchmark/plots/:
%     fig_scaling_time_k<K>.png    - runtime vs n (per family): deterministic
%                                    BSQR, built-in qr, randomized, randomized+R12.
%     fig_scaling_speedup_k<K>.png - randomized speedup over each baseline vs n.
%     fig_scaling_quality_k<K>.png - two rows vs n: ||R11^{-1}||_F / bound and
%                                    sigma_min(R11) / bound.
%     fig_blocksize_k<K>.png       - time, columns sampled, conditioning vs block
%                                    size (per family), batched + norm-weighted.
%     fig_sampling_k<K>.png        - uniform vs norm-weighted across families.
%   Each line is the seed mean; the shaded band spans the seed min/max. The three
%   scaling figures are one-metric-each (one former column apiece) so each has a
%   single clear legend.
%
% Options: 'k' (default 64), 'resultsdir', 'plotdir', 'formats' (default
%   {'png','pdf'} -- vector PDF for publication alongside PNG; pass {'png'} to
%   skip the PDF for quick iteration).
%   (cellstr, default {'png'}).

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

C = struct('det', [0 0 0], 'builtin', [0.55 0.55 0.55], ...
    'rm', [0.00 0.45 0.74], 'wc', [0.85 0.33 0.10], 'r12', [0.49 0.18 0.56], ...
    'uniform', [0.00 0.45 0.74], 'normweighted', [0.47 0.67 0.19]);

plot_scaling(opt, C);
plot_blocksize(opt, C);
plot_sampling(opt, C);
end

% ===========================================================================
function plot_scaling(opt, C)
% Three single-metric figures (one per former column), each with its own clear
% legend: runtime, speedup, and selection quality.
f = fullfile(opt.resultsdir, ['exp_scaling' opt.tag '.csv']);
if ~isfile(f); warning('missing %s', f); return; end
T = readtable(f, 'TextType', 'string');
fams = unique(T.family, 'stable');
scaling_figure(opt, C, T, fams, 'time');
scaling_figure(opt, C, T, fams, 'speedup');
scaling_quality_figure(opt, C, T, fams);
end

function scaling_figure(opt, C, T, fams, which)
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
            band_line(ax, T, base & T.method == "det", 'n', 'time_s', C.det, 'deterministic BSQR');
            band_line(ax, T, base & T.method == "builtin", 'n', 'time_s', C.builtin, 'built-in qr (dgeqp3)');
            band_line(ax, T, base & T.method == "rand" & T.mode == "running_mean", 'n', 'time_s', C.rm, 'randomized BSQR');
            band_line(ax, T, base & T.method == "rand_r12", 'n', 'time_s', C.r12, 'randomized + R_{12}');
            ylabel(ax, 'time (s)');
        case 'speedup'
            set(ax, 'XScale', 'log', 'YScale', 'log');
            speedup_line(ax, T, base, 'det', 'na', 'rand', 'running_mean', C.rm, 'randBSQR vs deterministic BSQR');
            speedup_line(ax, T, base, 'builtin', 'na', 'rand', 'running_mean', C.builtin, 'randBSQR vs built-in qr');
            speedup_line(ax, T, base, 'builtin', 'na', 'rand_r12', 'running_mean', C.r12, 'randBSQR + R_{12} vs built-in qr');
            yline(ax, 1, 'k--', 'HandleVisibility', 'off');
            ylabel(ax, 'speedup (t_{baseline} / t_{randBSQR})');
    end
    grid(ax, 'on'); xlabel(ax, 'n'); title(ax, fam, 'Interpreter', 'none');
end
lg = legend(ax, 'Orientation', 'horizontal');   % shared: every panel has the same lines
lg.Layout.Tile = 'south';
switch which
    case 'time';    title(tl, sprintf('Runtime vs n   (k=%d)', opt.k));
    case 'speedup'; title(tl, sprintf('Speedup over baselines (m=%d)', opt.k));
end
save_fig(fig, fullfile(opt.plotdir, ['fig_scaling_' which opt.tag]), opt.formats);
end

function scaling_quality_figure(opt, C, T, fams)
% Two rows of selection quality: Frobenius ||R11^{-1}||_F / bound (<=1, lower
% better) and spectral sigma_min(R11) / bound (>=1, higher better).
nf = numel(fams);
fig = figure('Position', [100 100 max(440 * nf, 720) 480], 'Color', 'w');
tl = tiledlayout(fig, 2, nf, 'TileSpacing', 'compact', 'Padding', 'compact');
methods = {"builtin", 'na', 'built-in qr', C.builtin; ...
           "det", 'na', 'deterministic BSQR', C.det; ...
           "rand", 'running_mean', 'randBSQR', C.rm};
ax = [];
for row = {'frob', 'smin'}
    for fi = 1:nf
        fam = fams(fi); base = T.family == fam;
        ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log');
        for mi = 1:size(methods, 1)
            mask = base & T.method == methods{mi, 1} & ...
                (methods{mi, 2} == "na" | T.mode == methods{mi, 2});
            quality_ratio_line(ax, T, mask, row{1}, methods{mi, 4}, methods{mi, 3});
        end
        ylim(ax, [0 0.6]);   % compressed axis; both ratios sit well under the bound here
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
save_fig(fig, fullfile(opt.plotdir, ['fig_scaling_quality' opt.tag]), opt.formats);
end

function quality_ratio_line(ax, T, mask, which, color, name)
% Both rows are (inverse-norm / bound) <= 1, lower = better:
%   'frob' -> ||R11^{-1}||_F / sqrt(k(n-k+1));
%   'smin' -> ||R11^{-1}||_2 / sqrt(1+k(n-k)) = 1/(sigma_min*sqrt(1+k(n-k)))
%             = bound / sigma_min, where bound = 1/sqrt(1+k(n-k)).
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

% ===========================================================================
function plot_blocksize(opt, C)
f = fullfile(opt.resultsdir, ['exp_blocksize' opt.tag '.csv']);
if ~isfile(f); warning('missing %s', f); return; end
T = readtable(f, 'TextType', 'string');
fams = unique(T.family, 'stable');
fig = figure('Position', [100 100 1200 360 * numel(fams)], 'Color', 'w');
tl = tiledlayout(fig, numel(fams), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for fi = 1:numel(fams)
    fam = fams(fi); base = T.family == fam;

    ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log', 'YScale', 'log');
    band_line(ax, T, base & T.sampling == "normweighted", 'block_size', 'time_s', C.normweighted, 'normweighted');
    grid(ax, 'on'); xlabel(ax, 'k_b'); ylabel(ax, 'time (s)');
    title(ax, sprintf('%s: runtime', fam), 'Interpreter', 'none');
    if fi == 1; legend(ax, 'Location', 'best'); end

    ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log', 'YScale', 'log');
    band_line(ax, T, base & T.sampling == "normweighted", 'block_size', 'tested_per_k', C.normweighted, 'normweighted');
    grid(ax, 'on'); xlabel(ax, 'k_b'); ylabel(ax, 'columns tested / m');
    title(ax, sprintf('%s: sampling cost', fam), 'Interpreter', 'none');

    ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log');
    cond_ratio_line2(ax, T, base & T.sampling == "normweighted", 'block_size', C.normweighted, 'normweighted');
    yline(ax, 1, 'k--', 'HandleVisibility', 'off');
    grid(ax, 'on'); xlabel(ax, 'k_b'); ylabel(ax, '||R_{11}^{-1}||_F / bound');
    title(ax, sprintf('%s: conditioning', fam), 'Interpreter', 'none');
end
title(tl, sprintf('Effect of block size (m=%d, n=%d)', T.k(1), T.n(1)));
save_fig(fig, fullfile(opt.plotdir, ['fig_blocksize' opt.tag]), opt.formats);
end

% ===========================================================================
function plot_sampling(opt, C)
f = fullfile(opt.resultsdir, ['exp_sampling' opt.tag '.csv']);
if ~isfile(f); warning('missing %s', f); return; end
T = readtable(f, 'TextType', 'string');
fams = unique(T.family, 'stable');
metrics = {'tested_per_k', 'columns tested / m'; 'frobinv', '||R_{11}^{-1}||_F'; 'time_s', 'time (s)'};
fig = figure('Position', [100 100 1300 380], 'Color', 'w');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for mi = 1:size(metrics, 1)
    ax = nexttile(tl); hold(ax, 'on');
    set(ax, 'TickLabelInterpreter', 'none');   % family names verbatim (no TeX subscripts)
    mu = zeros(numel(fams), 2); elo = zeros(numel(fams), 2); ehi = zeros(numel(fams), 2);
    for fi = 1:numel(fams)
        for si = 1:2
            samp = ["uniform", "normweighted"]; samp = samp(si);
            v = T.(metrics{mi, 1})(T.family == fams(fi) & T.sampling == samp);
            % min/max whiskers (not +/- std): honest for floored/skewed quantities
            % like tested/k (>= 1), and matches the min/max bands in the other figures.
            mu(fi, si) = mean(v); elo(fi, si) = mean(v) - min(v); ehi(fi, si) = max(v) - mean(v);
        end
    end
    cats = categorical(fams, fams);
    b = bar(ax, cats, mu);
    b(1).FaceColor = C.uniform; b(1).DisplayName = 'uniform';
    b(2).FaceColor = C.normweighted; b(2).DisplayName = 'normweighted';
    for si = 1:2
        xc = b(si).XEndPoints;
        errorbar(ax, xc, mu(:, si), elo(:, si), ehi(:, si), 'k', 'linestyle', 'none', 'HandleVisibility', 'off');
    end
    show_legend = (mi == 1);
    if strcmp(metrics{mi, 1}, 'time_s')   % overlay vendor-qr time per family
        bt = arrayfun(@(f) mean(T.time_s(T.family == f & T.method == "builtin")), fams);
        plot(ax, cats, bt, 'd', 'MarkerSize', 7, 'MarkerFaceColor', C.builtin, ...
            'MarkerEdgeColor', 'k', 'LineStyle', 'none', 'DisplayName', 'built-in qr');
        show_legend = true;
    end
    ylabel(ax, metrics{mi, 2}); grid(ax, 'on');
    title(ax, metrics{mi, 2});
    if show_legend; legend(ax, 'Location', 'best'); end
    if strcmp(metrics{mi, 1}, 'tested_per_k')
        % Log y-axis: keep the (small, ~block) norm-weighted bars visible by
        % putting the bar baseline and the lower y-limit below the data minimum.
        set(ax, 'YScale', 'log');
        lo = min(mu(mu > 0)) / 4; hi = max(mu(:)) * 1.6;
        for si = 1:2; b(si).BaseValue = lo; end
        ylim(ax, [lo, hi]);
    end
end
title(tl, sprintf('Uniform vs norm-weighted sampling (m=%d, n=%d, k_b=%d)', ...
    T.k(1), T.n(1), max(T.block_size)));
save_fig(fig, fullfile(opt.plotdir, ['fig_sampling' opt.tag]), opt.formats);
end

% ===========================================================================
% helpers
function [x, ym, ylo, yhi] = agg(T, mask, xname, yname)
sub = T(mask, :);
x = unique(sub.(xname));
ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(x)
    v = sub.(yname)(sub.(xname) == x(i));
    ym(i) = mean(v); ylo(i) = min(v); yhi(i) = max(v);
end
end

function band_line(ax, T, mask, xname, yname, color, name)
if ~any(mask); return; end
[x, ym, ylo, yhi] = agg(T, mask, xname, yname);
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
end

function speedup_line(ax, T, base, denMethod, denMode, numMethod, numMode, color, name)
% Per-seed speedup of the target (numMethod/numMode) over the baseline
% (denMethod/denMode): t_baseline / t_target.
ns = unique(T.n(base));
x = ns; ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(ns)
    sd = sort(T.seed(base & T.method == string(denMethod) & T.n == ns(i)));
    sp = [];
    for s = sd'
        tb = T.time_s(base & T.method == string(denMethod) & T.mode == string(denMode) & T.n == ns(i) & T.seed == s);
        tt = T.time_s(base & T.method == string(numMethod) & T.mode == string(numMode) & T.n == ns(i) & T.seed == s);
        if ~isempty(tb) && ~isempty(tt); sp(end+1) = tb(1) / tt(1); end %#ok<AGROW>
    end
    ym(i) = mean(sp); ylo(i) = min(sp); yhi(i) = max(sp);
end
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
end

function cond_ratio_line2(ax, T, mask, xname, color, name)
if ~any(mask); return; end
sub = T(mask, :);
sub.ratio = sub.frobinv ./ sub.osinsky;
x = unique(sub.(xname)); ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(x)
    v = sub.ratio(sub.(xname) == x(i));
    ym(i) = mean(v); ylo(i) = min(v); yhi(i) = max(v);
end
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
end

function save_fig(fig, stem, formats)
for i = 1:numel(formats)
    fmt = formats{i};
    switch fmt
        case 'png'; exportgraphics(fig, [stem, '.png'], 'Resolution', 150);
        case 'pdf'; exportgraphics(fig, [stem, '.pdf'], 'ContentType', 'vector');
        case 'eps'; exportgraphics(fig, [stem, '.eps'], 'ContentType', 'vector');
        otherwise; warning('unknown format %s', fmt);
    end
end
fprintf('Wrote %s.%s\n', stem, strjoin(formats, ','));
close(fig);
end
