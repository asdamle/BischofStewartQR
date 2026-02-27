# Optimization Candidates

- Generated: 2026-02-25T12:58:28.444
- Ranking model: `score = 0.5*impact + 0.3*confidence - 0.1*numerical_risk - 0.1*implementation_cost`

| candidate_id | location | rationale | expected_gain | risk_level | numerical_risk | implementation_cost | validation_plan | go_no_go | weighted_score |
|---|---|---|---|---|---|---|---|---|---:|
| K2 | src/kernel.jl (W-update and workspace data movement) | W-update dominates average phase share (39.58%) and correlates with slowdown (corr=0.920). | impact_norm=1.000; score_model=0.760 | low | low | medium | Run tests + tier1; verify identical outputs for fixed RNG seeds and no allocation regressions. | GO | 0.76 |
| K1 | src/kernel.jl (downdate/recompute and timed branches) | Downdate path has high average share (21.89%) with nontrivial recompute activity (0.692 normalized). | impact_norm=0.990; score_model=0.735 | medium | medium | medium | Run tests + tier1 + compare_results; verify residual/orth and recompute behavior on ill-conditioned ladder. | GO | 0.73519 |
| K3 | src/kernel.jl (_apply_householder_left! and surrounding call pattern) | Reflector application remains a major cost center (33.13% average phase share). | impact_norm=0.837; score_model=0.678 | medium | low | medium | Run tests + tier1; validate residual/orthogonality parity on square and short-wide cases. | GO | 0.67846 |
| K4 | src/kernel.jl (pivot selection loop) | Pivot selection is small but systematic (1.16%) and safe for micro-optimizations. | impact_norm=0.029; score_model=0.281 | low | low | low | Run tests including criterion-consistent pivot sequence; ensure pivot order invariants hold. | GO | 0.28118 |
| K5 | benchmark/bench_common.jl and benchmarking harness | Median BS allocation ratio vs dgeqp3 is 0.723; tighten measurement/alloc path to isolate kernel effects. | impact_norm=0.000; score_model=0.276 | low | none | low | Re-run benchmarks and confirm CSV schemas and fairness policy remain unchanged. | GO | 0.2765 |
