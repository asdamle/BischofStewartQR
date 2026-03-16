function run_publication_smoke_benchmark(varargin)
%RUN_PUBLICATION_SMOKE_BENCHMARK Fast publication benchmark sanity run.
%
%   run_publication_smoke_benchmark()
%   run_publication_smoke_benchmark(cfg_struct)

if nargin >= 1 && isstruct(varargin{1})
    cfg = varargin{1};
else
    cfg = struct();
end

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
default_outdir = fullfile(repo_root, 'matlab', 'benchmark', 'results', 'publication_mex_smoke');

cfg.outdir = getfield_default(cfg, 'outdir', default_outdir);
cfg.seeds = getfield_default(cfg, 'seeds', 20260314);
cfg.families = getfield_default(cfg, 'families', {'gaussian'});
cfg.square_ms = getfield_default(cfg, 'square_ms', 32);
cfg.short_ms = getfield_default(cfg, 'short_ms', 16);
cfg.short_aspects = getfield_default(cfg, 'short_aspects', 2);
cfg.warmup = getfield_default(cfg, 'warmup', 1);
cfg.samples = getfield_default(cfg, 'samples', 3);
cfg.norm_recomp_tol = getfield_default(cfg, 'norm_recomp_tol', sqrt(eps('double')));
cfg.bench_surface = getfield_default(cfg, 'bench_surface', "materialize_qrp");

run_publication_benchmarks(cfg);
end

function v = getfield_default(s, field, default)
if isfield(s, field)
    v = s.(field);
else
    v = default;
end
end
