# Lazy BSQR Scoring Assessment

This note assesses whether the Bischof-Stewart pivot score

`c_j = (1 + ||w_j||^2) / s_j`

can be evaluated lazily while preserving the exact same pivot sequence as the
current eager BSQR kernel in exact arithmetic.

## Current Eager Rule

At step `i`, the kernel in `src/kernel.jl` selects the remaining column with minimal

`c_j = (1 + ||w_j||^2) / s_j`

where:

- `s_j = ||a_j^(i)||^2` is the squared trailing column norm,
- `w_j = R11^{-1} r_j`,
- ties are resolved by the first minimum because the loop updates the pivot only
  when `cj < best_c`.

The state that matters for exact pivot equivalence is therefore:

- current permutation order,
- exact `s_j`,
- exact `w_j` or exact `||w_j||^2`,
- the strict first-min tie rule.

## BSQR-Specific Lazy Certificate

Suppose a pivot with exact `||w_*||^2` has been chosen, and a remaining column
has the exact update

`w_j^+ = [w_j - beta_j w_*; beta_j]`.

Then

`||w_j^+||^2 = ||w_j||^2 - 2 beta_j d_j + beta_j^2 (||w_*||^2 + 1)`

with `d_j = w_j^T w_*`. Minimizing the right-hand side over admissible `d_j`
and `beta_j` gives the exact lower bound

`||w_j^+||^2 >= ||w_j||^2 / (1 + ||w_*||^2)`.

This is the key BSQR-specific fact. It does not come from DGEQP3 machinery; it
comes directly from the BSQR `w` update.

If a column carries:

- a lower bound `Lw_j <= ||w_j||^2`, and
- an upper bound `Us_j >= s_j`,

then after one skipped pivot it still carries the exact-safe certificate

`c_j >= (1 + Lw_j / (1 + ||w_*||^2)) / Us_j`.

The denominator bound is cheap because `s_j` only decreases under a pivot step,
so any stale exact value remains a valid upper bound.

## Exact-Safe Lazy Policy

A viable exact-safe lazy policy is a branch-and-bound refresh loop:

1. Each remaining column stores `Lw_j` and `Us_j`.
2. Columns in a tracked/refreshed subset also have an exact current score.
3. The next column to refresh is the unrefreshed column with smallest certified
   lower bound `(1 + Lw_j) / Us_j`.
4. After each exact refresh, update the current best exact score.
5. Commit a pivot only when

   `best_exact_score < min_unrefreshed_lower_bound`.

This commit rule is exact-safe:

- every unrefreshed column has true score at least its lower bound,
- strict `<` preserves the current first-min tie behavior,
- if equality is possible, the algorithm must refresh more columns rather than
  commit early.

If all remaining columns are refreshed, the method degenerates to the eager
rule, so correctness is never at risk.

## Verdict

The assessment outcome is `conditionally viable`.

There is a coherent BSQR-specific lazy policy with an exact correctness
argument, but its practical upside appears limited:

- untouched columns start with `Lw_j = 0`, so their certificate is initially
  just `1 / Us_j`,
- the numerator bound becomes informative only after a column has been refreshed
  at least once,
- if many columns have similar residual norms or near-tied scores, the lazy
  loop may need to refresh most columns and collapse back to eager work.

The most promising regime is one where denominator separation is already strong,
or where a relatively small tracked set is reused across several steps.

## Implementation Sketch

If this is promoted into production, the internal state should be:

- `s_upper[j]`: stale upper bound on current `s_j`,
- `wnorm2_lb[j]`: stale lower bound on current `||w_j||^2`,
- transient exact score cache for refreshed columns in the current selection
  round,
- a lower-bound priority structure over remaining columns.

The pivot loop would:

1. refresh columns on demand in order of smallest certified lower bound,
2. keep the best exact score seen so far,
3. commit only under the strict certificate test above,
4. update certificates after each pivot with
   `wnorm2_lb[j] /= (1 + ||w_*||^2)` and
   `s_upper[j] = min(s_upper[j], exact_s_j_if_refreshed)`.

The fallback is explicit: if certificate separation never occurs, refresh all
remaining columns and recover the eager BSQR pivot exactly.

## BLAS-Only Prototype

The current repository now includes an internal prototype in
`src/kernel_lazy_blas.jl`:

- entrypoint: `BSPivotQR._bsqr_lazy_blas(...)`,
- workspace: `BSPivotQR.BSLazyBLASWorkspace`,
- optional stats: `BSPivotQR.BSLazyBLASStats`.

This prototype is intentionally not part of the public API. It exists for
benchmarking and evaluation only.

The refresh path is BLAS-oriented rather than Julia-threaded:

- candidate columns are refreshed in batches,
- `R11^{-1}R12_batch` is formed with a triangular solve on a block of columns,
- exact `||w_j||^2` values are taken from the solved block,
- exact trailing norms are then computed for the same refreshed columns.

The exact commit rule is unchanged:

`best_exact_score < min_unrefreshed_lower_bound`.

So the prototype preserves the eager first-min pivot semantics by construction,
and falls back to full refresh when certificates are not separating enough.

The benchmark harness can opt into the prototype with:

- `BS_BENCH_INCLUDE_LAZY=1`,
- `BS_LAZY_BATCH_SIZE=<int>` for a fixed refresh batch,
- `BS_LAZY_BATCH_FRACTION=<float>` for a size proportional to remaining columns,
- `BS_LAZY_BATCH_MIN=<int>`,
- `BS_LAZY_BATCH_MAX=<int>`.

By default the benchmark suite still runs only the eager BSQR kernel and
`dgeqp3`.

## Scope Limit

This note only establishes an exact-arithmetic assessment. It does not claim a
floating-point safe lazy rule. Any production implementation would need a
separate stability analysis before enabling lazy selection in the main kernel.
