# Randomized Bischof–Stewart Column Selection

The randomized column-selection algorithm implemented in `matlab_rand/`
(`bsqr_rand`). It selects `k` well-conditioned columns with the *same*
worst-case guarantee on `‖R₁₁⁻¹‖` as deterministic Bischof–Stewart pivoted QR,
but without the `O(n k²)` all-column scan that dominates the deterministic cost
when `k ≪ n`. Sections 1–5 are prose; §6 is a LaTeX draft of the algorithm
statement (`algorithm` / `algpseudocode`).

---

## 1. Problem and guarantee

Given `A ∈ ℝ^{m×n}` and a target `k ≤ min(m,n)`, choose an index set `J` of `k`
columns. Let `A(:,J) = Q [R₁₁; 0]` be the (unpivoted) QR factorization of the
selected columns, with `R₁₁ ∈ ℝ^{k×k}` upper triangular. Selection quality is
measured by `‖R₁₁⁻¹‖`: the smaller it is, the farther the selected columns are
from linear dependence.

The canonical setting is **orthonormal-row input**: `A` is `k×n` with `A Aᵀ = I_k`
(for example, `A = Vₖᵀ`, the leading-`k` right singular vectors of some matrix to
be column-subset-selected). There, a Bischof–Stewart selection achieves Osinsky's
bounds

```
‖R₁₁⁻¹‖_F ≤ √( k (n − k + 1) ),        ‖R₁₁⁻¹‖₂ ≤ √( 1 + k (n − k) ).
```

The randomized algorithm maintains these same bounds. (The spectral bound
follows from the Frobenius one: orthonormal rows force `σ_max(R₁₁) ≤ 1`, so the
other `k − 1` singular values of `R₁₁⁻¹` each contribute at least 1 to
`‖R₁₁⁻¹‖_F²`.)

## 2. The deterministic criterion (background)

Columns are selected one at a time. After `i` selections we have the leading
triangular block `R₁₁ ∈ ℝ^{i×i}` and an accumulated orthogonal `Qᵢ` (a product
of `i` Householder reflectors). Bring a candidate column `aⱼ` into the current
frame, `x = Qᵢᵀ aⱼ`, and split it into a **head** `x_{1:i}` and a **tail**
`x_{i+1:m}`. Define

```
ρⱼ²  = ‖x_{i+1:m}‖²            (the squared residual; ±ρⱼ is the next diagonal of R₁₁)
wⱼ   = R₁₁⁻¹ x_{1:i}           (the candidate's coefficients in the selected basis)
cⱼ   = (1 + ‖wⱼ‖²) / ρⱼ².      (the Bischof–Stewart criterion)
```

The criterion is exactly the **growth of the squared inverse Frobenius norm**:
writing `f_i = ‖R₁₁⁻¹‖_F²` after `i` selections, the rank-one bordered-inverse
formula gives, when column `j` is appended,

```
f_{i+1} = f_i + cⱼ.            (exact — Eq. (growth))
```

So the criterion penalizes a candidate that is nearly in the span of the already
selected columns (`‖wⱼ‖` large) or that has a small residual (`ρⱼ` small): either
inflates `‖R₁₁⁻¹‖_F`. The deterministic kernel selects `argminⱼ cⱼ` at every step,
which requires maintaining `R₁₁⁻¹R₁₂` and running residual norms for **all**
trailing columns — the `O(n k²)` total work that dominates the short-wide
(`k ≪ n`) regime.

## 3. Randomization: sample, don't scan

The guarantee does not need the exact `argmin` — only that each step's increment
`cⱼ` keeps the running `f` under a per-step threshold. Define the
**running-mean threshold**

```
θᵢ = (f_i + n − 2i) / (k − i).
```

For orthonormal-row input this is precisely the `ρ²`-weighted mean of the criterion
over the remaining columns (using `Σ ρⱼ² = k − i` and `Σ ‖wⱼ‖² = f_i − i`). Two
consequences follow:

1. **Feasibility.** The minimum is at most the mean, so at least one remaining
   column satisfies `cⱼ ≤ θᵢ` — accepting *any* column that meets the threshold
   can never stall the algorithm.
2. **The bound.** Enforcing `cⱼ ≤ θᵢ` at every step preserves
   `f_i ≤ F̂_i = i(n−i+1)/(k−i+1)` by induction (the per-step recursion
   `f_{i+1} ≤ ((k−i+1) f_i + n − 2i)/(k−i)` maps `F̂_i` to exactly `F̂_{i+1}`),
   hence `f_k ≤ F̂_k = k(n−k+1)` — Osinsky's bound, together with its
   per-singular-value refinement (writeup Thm. 6.6, `thm:greedy` /
   Cor. 6.1, `cor:sv-hierarchy`).

In practice a large fraction of the remaining columns clears a mean threshold, so
it suffices to **sample** a small block of candidates, bring only those up to
date, and accept one that qualifies — never touching the columns we do not sample.
(The rare block with no qualifying column is handled by the feasibility net of
§4.) This removes the `O(n)` per-step scan; per-step work depends on the block
size, not on `n`.

Within a block, candidates are drawn **without replacement, weighted by their
initial squared column norms** `gⱼ = ‖aⱼ‖²` (robust across leverage profiles;
uniform sampling is available when leverage is known to be flat), using
Efraimidis–Spirakis keys `−log(uⱼ)/gⱼ`, `uⱼ ∼ Unif(0,1)`, smallest keys first.
Sampled-but-rejected columns stay in the pool and can be drawn again in later
blocks.

## 4. The algorithm (batched in-block selection — the default)

Bringing a sampled block into the current frame costs one application of the
accumulated reflectors — the dominant per-block cost. To amortize it, once a
block is reduced we run Bischof–Stewart **within the block**: repeatedly take
the in-block minimizer, append it, update `f` and `θ`, and downdate the rest of
the block by the new reflector — continuing while the in-block minimizer still
clears the updated threshold. One expensive apply thus yields several
selections; for the canonical `k×n` input this cuts the dominant cost from
`O(k⁴)` to `O(k³)`. The default block size on this path is `b = k`.

The loop, per sampled block `B` (with reduced columns `X = Qᵢᵀ A(:,B)`):

1. Let `i = |J|`, `θ = (f + n − 2i)/(k − i)`.
2. For each still-active column `x` of `X`: `ρ² = ‖x_{i+1:m}‖²`,
   `w = R₁₁⁻¹ x_{1:i}`, `c = (1 + ‖w‖²)/ρ²`.
3. Let `(c★, x★)` be the active column of smallest `c`.
4. If `c★ > θ`, the block is exhausted — break and sample a new block.
5. Otherwise append `x★`: set the new `R₁₁` column from `x★_{1:i}`, form a
   Householder reflector from `x★_{i+1:m}` (its `β` is the new diagonal), store the
   reflector, add the original index of `x★` to `J`, and set `f ← f + c★`.
6. Downdate the remaining active columns of `X` by the new reflector and repeat
   from step 1.

Stop when `|J| = k`.

**Single-select variant (`batched = false`).** Each step draws one weighted
visiting order over *all* remaining columns and scans it block by block
(default block size `⌈k/2⌉` clamped to `[16, 64]`) until a column qualifies:
`pick = 'best_in_block'` (the default) accepts the first block whose minimizer
has `c ≤ θ`, `pick = 'first'` accepts the first qualifying column. That column
is appended and the next step resamples. This costs one reflector apply per
*scanned block*, at least one block per selection — the `O(k⁴)` path — but the
realized `‖R₁₁⁻¹‖_F` lands slightly closer to the deterministic value; the
guarantee is identical.

**Threshold variants.** `running_mean` (above) is the default. The more
permissive `worstcase_allowance` mode sets `θᵢ = F̂_{i+1} − f_i`, spending the
slack between the *actual* running `f_i` and the deterministic worst case
`F̂_i`; it accepts with fewer samples and still keeps `f_k ≤ k(n−k+1)`, at the
cost of a looser realized `‖R₁₁⁻¹‖_F`. A `slack ≥ 1` multiplier loosens either
threshold (and the bound) proportionally.

**Safeguards.**
- *Feasibility net.* Once as many candidates have been sampled since the last
  selection as there are remaining columns, the next "block" is a forced pass
  over **all** remaining columns. Its minimizer clears the threshold in exact
  arithmetic (feasibility above), so it is accepted even if rounding leaves it
  marginally over.
- *Rank guard.* If in that forced pass every remaining column has `ρ² ≈ 0` (so
  `c = ∞`), the input is numerically rank-deficient for this `k`; the algorithm
  stops with an error rather than propagating `∞` into `f` and `R₁₁`.

**Outputs.** The permutation with `J` first; the economy orthogonal factor `Q`
(`m×k`), accumulated lazily from the stored Householder reflectors — one
`O(mk²)` `dorgqr` pass, skipped when only the permutation is requested; `R₁₁`;
instrumentation; and — on request, via one extra reflector apply to the
unselected columns — `R₁₂`.

## 5. Cost

Deterministic BSQR maintains the criterion for all trailing columns: `O(n i)`
work at step `i`, `O(n k²)` total. The randomized algorithm touches only the
columns it samples:

- one-time: `O(n)` pool setup, plus `O(m n)` column norms for norm-weighted
  sampling (skipped when sampling uniformly);
- per sampled block: one reflector apply, `O(m b i)`, plus one triangular solve,
  `O(b i²)` — independent of `n`;
- total, for the canonical `k×n` input on the batched path: `O(k³)` when blocks
  yield a constant fraction of their columns (the typical case; the
  instrumentation reports realized sample counts).

For `k ≪ n` everything that grows with `n` is a single linear pass, so the
speedup over the deterministic `O(n k²)` grows with `n` while `‖R₁₁⁻¹‖` stays
under the same guarantee.

---

## 6. LaTeX draft

Compile with `\usepackage{algorithm}`, `\usepackage{algpseudocode}`, and
`\usepackage{amsmath,amssymb}` (the latter for `\mathbb`). The growth identity
referenced in the main routine is (for a candidate column brought into the
current frame, `x = Qᵢᵀ aⱼ`, after `i` selections):

```latex
Appending column $j$ grows the squared inverse Frobenius norm by exactly
\begin{equation}\label{eq:growth}
  \bigl\|R_{11}^{-1}\bigr\|_F^2
  \;\longleftarrow\;
  \bigl\|R_{11}^{-1}\bigr\|_F^2 + c_j,
  \qquad
  c_j = \frac{1 + \|w_j\|^2}{\rho_j^2},
\end{equation}
where $w_j = R_{11}^{-1} x_{1:i}$ and $\rho_j = \|x_{i+1:m}\|$.
```

```latex
\begin{algorithm}[t]
\caption{Randomized Bischof--Stewart column selection (batched; default path)}
\label{alg:randbsqr}
\begin{algorithmic}[1]
\Require $A\in\mathbb{R}^{m\times n}$; subset size $k\le\min(m,n)$; block size $b$ (default $b=k$);
         weights $g_j$ (default $g_j=\|a_j\|^2$, or $g_j\equiv1$ for uniform sampling)
\Ensure index set $J$, $|J|=k$, and upper-triangular $R\in\mathbb{R}^{k\times k}$ with
        $\|R^{-1}\|_F\le\sqrt{k(n-k+1)}$ for orthonormal-row $A$
\State $J\gets(\,)$;\quad $\mathcal{R}\gets\{1,\dots,n\}$;\quad $R\gets 0$;\quad $f\gets 0$
       \Comment{$f=\|R_{1:i,1:i}^{-1}\|_F^2$, the squared inverse Frobenius norm}
\While{$|J|<k$}
  \State $B\gets\Call{SampleBlock}{\mathcal{R},\,b,\,g}$
         \Comment{$\le b$ distinct columns, without replacement}
  \State $X\gets\Call{Reduce}{A_{:,B},\,|J|}$
         \Comment{one application of the $|J|$ accumulated reflectors}
  \State mark every column of $X$ \emph{active}
  \While{$|J|<k$ \textbf{and} some column of $X$ is active}
    \State $i\gets|J|$;\qquad $\theta\gets(f+n-2i)/(k-i)$
           \Comment{$=\rho^2$-weighted mean of $c$ over $\mathcal{R}$ for orthonormal-row $A$}
    \ForAll{active columns $x$ of $X$}
      \State $\rho^2\gets\|x_{i+1:m}\|^2$;\quad
             $w\gets R_{1:i,\,1:i}^{-1}\,x_{1:i}$;\quad
             $c\gets(1+\|w\|^2)/\rho^2$
    \EndFor
    \State $(c^\star,x^\star)\gets$ active column of least $c$
    \If{$c^\star>\theta$}
      \State \textbf{break} \Comment{no in-block column meets the bound; resample}
    \EndIf
    \State $R_{1:i,\,i+1}\gets x^\star_{1:i}$;\quad
           $(\beta,v,\tau)\gets\Call{Householder}{x^\star_{i+1:m}}$;\quad
           $R_{i+1,\,i+1}\gets\beta$
    \State store reflector $(v,\tau)$;\quad append $x^\star$'s column index to $J$;\quad
           remove it from $\mathcal{R}$;\quad deactivate $x^\star$
    \State $f\gets f+c^\star$ \Comment{exact, by Eq.~\eqref{eq:growth}}
    \State downdate the remaining active columns of $X$ by $(v,\tau)$
           \Comment{in-block reflector apply}
  \EndWhile
\EndWhile
\State \Return $J$, $R$, and the stored reflectors (an implicit $Q$)
\end{algorithmic}
\end{algorithm}
```

```latex
\begin{algorithm}[t]
\caption{Subroutines}
\begin{algorithmic}[1]
\Procedure{SampleBlock}{$\mathcal{R},\,b,\,g$}
  \Comment{weighted sampling without replacement (Efraimidis--Spirakis)}
  \State for each $j\in\mathcal{R}$ draw $u_j\sim\mathrm{Unif}(0,1)$ and set
         $\mathrm{key}_j\gets-\log(u_j)/\max(g_j,\varepsilon)$
         \Comment{$\varepsilon>0$ guards zero weights}
  \State \Return the indices of the $\min(b,|\mathcal{R}|)$ smallest keys
         \Comment{uniform sampling is $g_j\equiv1$}
\EndProcedure
\Procedure{Reduce}{$M,\,i$}
  \Comment{bring columns into the current frame}
  \State \Return $Q_i^\top M$, where $Q_i=H_1H_2\cdots H_i$ is the product of the
         $i$ stored reflectors \Comment{one compact-WY block apply}
\EndProcedure
\Procedure{Householder}{$y$}
  \Comment{reflector zeroing $y$ below its first entry}
  \State choose $(v,\tau)$ with $(I-\tau vv^\top)\,y=\beta e_1$;\quad \Return $(\beta,v,\tau)$
\EndProcedure
\end{algorithmic}
\end{algorithm}
```

> **Remark (single-select variant).** Replacing the inner `\While` by a scan that
> stops at the first sampled column with `$c\le\theta$` and then resamples gives
> the `batched=false` path: each `\Call{Reduce}{}` yields at most one selection,
> instead of the several per block above. The outer structure, threshold, and
> guarantee are unchanged.

> **Remark (safeguards).** Two safeguards are omitted from the listing for
> clarity. Once as many candidates have been sampled since the last selection as
> `$|\mathcal{R}|$`, the next block is a forced pass over all of `$\mathcal{R}$`;
> for orthonormal-row `$A$` its minimizer satisfies `$c\le\theta$` in exact
> arithmetic (the minimum is at most the $\rho^2$-weighted mean, which equals
> $\theta$), so it is accepted even on a rounding tie. If that pass finds
> `$\rho^2\approx 0$` for every remaining column, `$A$` is numerically
> rank-deficient for this `$k$` and the algorithm stops with an error.
