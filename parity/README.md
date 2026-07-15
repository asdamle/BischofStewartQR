# Cross-language parity fixtures

Generated artifacts — do not edit by hand. See `docs/VALIDATION_AND_PERF_PLAN.md` (V3).

Each fixture is a `parity_zoo` input matrix plus the expected outputs of the validation
oracle (`matlab/tests/oracle_bsqr.m`): full pivot vector, `R`, and (where flagged)
`R11 \ R12`. Values are printed with `%.17g` so doubles round-trip exactly; MATLAB and
Julia consumers therefore see bit-identical inputs. Inputs are screened at generation
time so no pivot step has a criterion near-tie (`gap_min` column of `manifest.csv`),
which is what makes exact pivot-sequence comparison fair across BLAS runtimes. The
acceptance tolerances (`rtol_R`, `rtol_rinv`, `rtol_crit` columns) are defined once on
the zoo members in `matlab/tests/parity_zoo.m` and flow to both consumers through the
manifest.

- Regenerate: `matlab -batch "addpath('matlab/tests'); generate_parity_fixtures"`
- MATLAB consumer: `matlab/tests/test_parity_fixtures.m`
- Julia consumer: `julia/test/test_parity_fixtures.jl`

Regenerate only when `matlab/tests/parity_zoo.m` changes; `testFixturesMatchCurrentZoo`
fails if fixtures and zoo drift apart.
