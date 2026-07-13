#!/bin/zsh
# Phase 6 publication rerun driver (docs/PUBLICATION_READINESS_PLAN.md).
#
# Runs the full benchmark regeneration sequentially with per-step logs and
# resumable step markers, wrapped in caffeinate so the system does not sleep
# when the screen locks. Launch it detached so it also survives the terminal
# (and any Claude session) going away:
#
#   nohup ./phase6_rerun.sh > /dev/null 2>&1 &
#   tail -f .phase6/console.log          # watch progress
#   ./phase6_rerun.sh --status           # step checklist
#
# Notes:
#   * Keep the machine on AC power with the lid open (caffeinate keeps the
#     system awake through screen lock, but battery + closed lid still sleeps).
#   * Do not use the machine for other work while it runs -- these are timing
#     benchmarks. Total wall time is expected to be a few hours.
#   * Re-running the script skips completed steps (markers in .phase6/).
#     Delete .phase6/DONE_<step> to force a step to rerun, or the whole
#     .phase6/ directory to start over.

set -uo pipefail
cd "$(dirname "$0")"

STATE=.phase6
mkdir -p "$STATE"

if [[ "${1:-}" == "--status" ]]; then
  for s in build_mex matlab_publication julia_publication matlab_smoke_gate \
           julia_smoke_gate rand_experiments rpqr_comparison largen_scaling \
           approx_suites rand_benchmarks final_suites; do
    if [[ -f "$STATE/DONE_$s" ]]; then echo "[done] $s"; else echo "[    ] $s"; fi
  done
  exit 0
fi

# Self-wrap in caffeinate (-i idle sleep, -s system sleep on AC) so the run
# continues through screen lock. nohup at launch handles terminal hangup.
if [[ -z "${PHASE6_CAFFEINATED:-}" ]]; then
  export PHASE6_CAFFEINATED=1
  exec caffeinate -is "$0" "$@"
fi

exec > "$STATE/console.log" 2>&1

# Vector figure formats to match the committed publication artifacts.
export BS_MATLAB_PUB_FIG_FORMATS=png,pdf,eps
export BS_PUB_FIG_FORMATS=png,pdf,eps

log() { echo "[$(date '+%H:%M:%S')] $*"; }

run_step() {
  local name=$1; shift
  if [[ -f "$STATE/DONE_$name" ]]; then log "skip $name (done)"; return 0; fi
  log "start $name"
  if "$@" > "$STATE/$name.log" 2>&1; then
    touch "$STATE/DONE_$name"; log "done  $name"
  else
    log "FAIL  $name (see $STATE/$name.log)"; exit 1
  fi
}

# --- preflight -------------------------------------------------------------
command -v matlab >/dev/null || { log "matlab not on PATH"; exit 1; }
command -v julia  >/dev/null || { log "julia not on PATH"; exit 1; }
[[ -f ext_comparisons/Adaptive-Randomized-Pivoting-main/code/rejection_rpqr.m ]] \
  || { log "ext_comparisons missing (needed for rpqr/largen/approx suites); see matlab_rand/README.md"; exit 1; }
[[ -z "$(git status --porcelain)" ]] || log "WARNING: working tree not clean"
pmset -g ps | grep -q "AC Power" || log "WARNING: on battery -- plug in AC (lid-closed battery sleep ignores caffeinate)"

# --- steps -----------------------------------------------------------------
run_step build_mex matlab -batch "addpath('matlab'); build_bsqr_mex; addpath('matlab_rand'); build_bsqr_rand_mex"

run_step matlab_publication matlab -batch "addpath('matlab'); addpath('matlab/benchmark'); run_publication_benchmarks"

run_step julia_publication zsh -c "julia --project=julia/benchmark julia/benchmark/run_publication_benchmarks.jl && julia --project=julia/benchmark julia/benchmark/plot_publication_results.jl"

# Gates compare fresh smoke runs against the recorded Phase 0 smoke baselines
# (the full-run case grid does not overlap the smoke grid, so smoke-vs-smoke
# is the meaningful comparison).
run_step matlab_smoke_gate matlab -batch "addpath('matlab'); addpath('matlab/benchmark'); run_publication_smoke_benchmark; check_publication_perf_gate('matlab/benchmark/results/tmp_phase0_baseline/publication_timings.csv', 'matlab/benchmark/results/publication_mex_smoke/publication_timings.csv', 0.10)"

run_step julia_smoke_gate zsh -c "julia --project=julia/benchmark julia/benchmark/run_publication_smoke_benchmark.jl && julia --project=julia/benchmark julia/benchmark/check_publication_perf_gate.jl julia/benchmark/results/tmp_phase0_baseline/publication_timings.csv julia/benchmark/results/publication_smoke/publication_timings.csv 0.10"

run_step rand_experiments matlab -batch "addpath('matlab'); addpath('matlab/mex'); addpath('matlab_rand'); addpath('matlab_rand/mex'); addpath('matlab_rand/benchmark'); for K = [64 128 256]; run_rand_experiments('k', K); plot_rand_experiments('k', K); end"

run_step rpqr_comparison matlab -batch "addpath('matlab_rand'); addpath('matlab_rand/mex'); addpath('matlab_rand/benchmark'); for K = [64 128 256]; run_rpqr_comparison('k', K); plot_rpqr_comparison('k', K); end"

run_step largen_scaling matlab -batch "addpath('matlab_rand'); addpath('matlab_rand/mex'); addpath('matlab_rand/benchmark'); run_largen_scaling; plot_largen_scaling"

run_step approx_suites matlab -batch "addpath('matlab_rand'); addpath('matlab_rand/mex'); addpath('matlab_rand/benchmark'); run_approx_comparison; plot_approx_comparison; run_approx_synth_comparison; plot_approx_comparison('tag','_synth'); run_approx_cond_comparison; plot_approx_cond_comparison; plot_approx_cond_comparison('norm','2')"

run_step rand_benchmarks matlab -batch "addpath('matlab'); addpath('matlab/mex'); addpath('matlab_rand'); addpath('matlab_rand/mex'); addpath('matlab_rand/benchmark'); run_rand_benchmarks"

run_step final_suites zsh -c "matlab -batch \"addpath('matlab/tests'); run_tests; addpath('matlab_rand'); addpath('matlab_rand/tests'); run_rand_tests\" && julia --project=julia -e 'using Pkg; Pkg.test()'"

log "PHASE 6 SWEEP COMPLETE -- review results, then commit artifacts and update the stale-flagged doc tables."
