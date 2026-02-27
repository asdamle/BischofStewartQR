# Performance Assessment Conclusion

- Generated: 2026-02-25T12:58:28.517
- Validation status: **PASS**
- Mean speedup (`dgeqp3/bsqr`): 0.57793
- Median speedup (`dgeqp3/bsqr`): 0.51975

## Project Status
- Correctness/guardrails currently pass in latest validation artifacts.
- Performance remains below dgeqp3 in most assessed cases (speedup < 1 dominates).

## Bottleneck Diagnosis
- Runtime is concentrated in apply/W-update/downdate phases.
- Mean shares: apply=33.129%, w_update=39.583%, downdate=21.894%.

## Top Opportunities
- K2 at src/kernel.jl (W-update and workspace data movement) (score=0.76): W-update dominates average phase share (39.58%) and correlates with slowdown (corr=0.920).
- K1 at src/kernel.jl (downdate/recompute and timed branches) (score=0.73519): Downdate path has high average share (21.89%) with nontrivial recompute activity (0.692 normalized).
- K3 at src/kernel.jl (_apply_householder_left! and surrounding call pattern) (score=0.67846): Reflector application remains a major cost center (33.13% average phase share).

## No-Go Areas (Stability Risk)
- Do not alter Bischof-Stewart pivot criterion semantics.
- Do not relax numerical safeguards (norm recomputation/rank-stop tolerances) without re-validating quality guardrails.
- Do not adopt major algorithm redesigns (blocked/alternative formulations) in this campaign.
