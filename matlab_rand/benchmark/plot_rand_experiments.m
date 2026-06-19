function plot_rand_experiments(varargin)
%PLOT_RAND_EXPERIMENTS Publication plots from run_rand_experiments CSVs.
%
%   Reads the k-tagged exp_*_k<K>.csv from the results dir and writes k-tagged
%   figures fig_*_k<K>.png to matlab_rand/benchmark/plots/:
%     fig_scaling_k<K>.png   - time, speedup, conditioning vs n (per family),
%                              deterministic vs randomized (both modes + R12).
%     fig_blocksize_k<K>.png - time, columns sampled, conditioning vs block size
%                              (per family), uniform vs norm-weighted sampling.
%     fig_sampling_k<K>.png  - uniform vs norm-weighted across families.
%   Each line is the seed mean; the shaded band spans the seed min/max.
%
% Options: 'k' (default 64), 'resultsdir', 'plotdir', 'formats'
%   (cellstr, default {'png'}).

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

C = struct('det', [0 0 0], 'rm', [0.00 0.45 0.74], 'wc', [0.85 0.33 0.10], ...
    'r12', [0.49 0.18 0.56], ...
    'uniform', [0.00 0.45 0.74], 'normweighted', [0.47 0.67 0.19]);

plot_scaling(opt, C);
plot_blocksize(opt, C);
plot_sampling(opt, C);
end

% ===========================================================================
function plot_scaling(opt, C)
f = fullfile(opt.resultsdir, ['exp_scaling' opt.tag '.csv']);
if ~isfile(f); warning('missing %s', f); return; end
T = readtable(f, 'TextType', 'string');
fams = unique(T.family, 'stable');
fig = figure('Position', [100 100 1200 360 * numel(fams)], 'Color', 'w');
tl = tiledlayout(fig, numel(fams), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for fi = 1:numel(fams)
    fam = fams(fi);
    base = T.family == fam;

    % (1) time vs n
    ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log', 'YScale', 'log');
    band_line(ax, T, base & T.method == "det", 'n', 'time_s', C.det, 'deterministic');
    band_line(ax, T, base & T.method == "rand" & T.mode == "running_mean", 'n', 'time_s', C.rm, 'rand running\_mean');
    band_line(ax, T, base & T.method == "rand" & T.mode == "worstcase_allowance", 'n', 'time_s', C.wc, 'rand worstcase');
    band_line(ax, T, base & T.method == "rand_r12", 'n', 'time_s', C.r12, 'rand + R_{12}');
    grid(ax, 'on'); xlabel(ax, 'n'); ylabel(ax, 'time (s)');
    title(ax, sprintf('%s: runtime', fam), 'Interpreter', 'none');
    if fi == 1; legend(ax, 'Location', 'northwest'); end

    % (2) speedup vs n (per seed det/rand ratio)
    ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log', 'YScale', 'log');
    speedup_line(ax, T, base, 'rand', 'running_mean', C.rm, 'running\_mean');
    speedup_line(ax, T, base, 'rand', 'worstcase_allowance', C.wc, 'worstcase');
    speedup_line(ax, T, base, 'rand_r12', 'running_mean', C.r12, 'rand + R_{12}');
    yline(ax, 1, 'k--', 'HandleVisibility', 'off');
    grid(ax, 'on'); xlabel(ax, 'n'); ylabel(ax, 'speedup (t_{det}/t_{rand})');
    title(ax, sprintf('%s: speedup', fam), 'Interpreter', 'none');

    % (3) conditioning vs n: ||R11^{-1}||_F / Osinsky
    ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log');
    cond_ratio_line(ax, T, base & T.method == "det", C.det, 'deterministic');
    cond_ratio_line(ax, T, base & T.method == "rand" & T.mode == "running_mean", C.rm, 'rand running\_mean');
    cond_ratio_line(ax, T, base & T.method == "rand" & T.mode == "worstcase_allowance", C.wc, 'rand worstcase');
    yline(ax, 1, 'k--', 'HandleVisibility', 'off');  % Osinsky bound
    grid(ax, 'on'); xlabel(ax, 'n'); ylabel(ax, '||R_{11}^{-1}||_F / Osinsky');
    title(ax, sprintf('%s: conditioning', fam), 'Interpreter', 'none');
end
title(tl, sprintf('Randomized vs deterministic BSQR, k=%d (with and without R_{12})', opt.k));
save_fig(fig, fullfile(opt.plotdir, ['fig_scaling' opt.tag]), opt.formats);
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
    band_line(ax, T, base & T.sampling == "uniform", 'block_size', 'time_s', C.uniform, 'uniform');
    band_line(ax, T, base & T.sampling == "normweighted", 'block_size', 'time_s', C.normweighted, 'normweighted');
    grid(ax, 'on'); xlabel(ax, 'block size'); ylabel(ax, 'time (s)');
    title(ax, sprintf('%s: runtime', fam), 'Interpreter', 'none');
    if fi == 1; legend(ax, 'Location', 'northwest'); end

    ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log', 'YScale', 'log');
    band_line(ax, T, base & T.sampling == "uniform", 'block_size', 'tested_per_k', C.uniform, 'uniform');
    band_line(ax, T, base & T.sampling == "normweighted", 'block_size', 'tested_per_k', C.normweighted, 'normweighted');
    grid(ax, 'on'); xlabel(ax, 'block size'); ylabel(ax, 'columns tested / k');
    title(ax, sprintf('%s: sampling cost', fam), 'Interpreter', 'none');

    ax = nexttile(tl); hold(ax, 'on'); set(ax, 'XScale', 'log');
    cond_ratio_line2(ax, T, base & T.sampling == "uniform", 'block_size', C.uniform, 'uniform');
    cond_ratio_line2(ax, T, base & T.sampling == "normweighted", 'block_size', C.normweighted, 'normweighted');
    yline(ax, 1, 'k--', 'HandleVisibility', 'off');
    grid(ax, 'on'); xlabel(ax, 'block size'); ylabel(ax, '||R_{11}^{-1}||_F / Osinsky');
    title(ax, sprintf('%s: conditioning', fam), 'Interpreter', 'none');
end
title(tl, sprintf('Effect of block size (running\\_mean threshold, k=%d, n=%d)', T.k(1), T.n(1)));
save_fig(fig, fullfile(opt.plotdir, ['fig_blocksize' opt.tag]), opt.formats);
end

% ===========================================================================
function plot_sampling(opt, C)
f = fullfile(opt.resultsdir, ['exp_sampling' opt.tag '.csv']);
if ~isfile(f); warning('missing %s', f); return; end
T = readtable(f, 'TextType', 'string');
fams = unique(T.family, 'stable');
metrics = {'tested_per_k', 'columns tested / k'; 'frobinv', '||R_{11}^{-1}||_F'; 'time_s', 'time (s)'};
fig = figure('Position', [100 100 1300 380], 'Color', 'w');
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for mi = 1:size(metrics, 1)
    ax = nexttile(tl); hold(ax, 'on');
    mu = zeros(numel(fams), 2); er = zeros(numel(fams), 2);
    for fi = 1:numel(fams)
        for si = 1:2
            samp = ["uniform", "normweighted"]; samp = samp(si);
            v = T.(metrics{mi, 1})(T.family == fams(fi) & T.sampling == samp);
            mu(fi, si) = mean(v); er(fi, si) = std(v);
        end
    end
    b = bar(ax, categorical(fams, fams), mu); b(1).FaceColor = C.uniform; b(2).FaceColor = C.normweighted;
    for si = 1:2
        xc = b(si).XEndPoints;
        errorbar(ax, xc, mu(:, si), er(:, si), 'k', 'linestyle', 'none', 'HandleVisibility', 'off');
    end
    ylabel(ax, metrics{mi, 2}); grid(ax, 'on');
    title(ax, metrics{mi, 2});
    if mi == 1; legend(ax, {'uniform', 'normweighted'}, 'Location', 'northwest'); end
    if strcmp(metrics{mi, 1}, 'tested_per_k'); set(ax, 'YScale', 'log'); end
end
title(tl, sprintf('Uniform vs norm-weighted sampling across families (k=%d, n=%d, block=%d)', ...
    T.k(1), T.n(1), T.block_size(1)));
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

function speedup_line(ax, T, base, method, mode, color, name)
ns = unique(T.n(base));
x = ns; ym = zeros(size(x)); ylo = ym; yhi = ym;
for i = 1:numel(ns)
    sd = sort(T.seed(base & T.method == "det" & T.n == ns(i)));
    sp = [];
    for s = sd'
        td = T.time_s(base & T.method == "det" & T.n == ns(i) & T.seed == s);
        tr = T.time_s(base & T.method == string(method) & T.mode == string(mode) & T.n == ns(i) & T.seed == s);
        if ~isempty(td) && ~isempty(tr); sp(end+1) = td(1) / tr(1); end %#ok<AGROW>
    end
    ym(i) = mean(sp); ylo(i) = min(sp); yhi(i) = max(sp);
end
fill(ax, [x; flipud(x)], [ylo; flipud(yhi)], color, 'FaceAlpha', 0.15, ...
    'EdgeColor', 'none', 'HandleVisibility', 'off');
plot(ax, x, ym, '-o', 'Color', color, 'MarkerFaceColor', color, ...
    'LineWidth', 1.5, 'MarkerSize', 4, 'DisplayName', name);
end

function cond_ratio_line(ax, T, mask, color, name)
cond_ratio_line2(ax, T, mask, 'n', color, name);
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
