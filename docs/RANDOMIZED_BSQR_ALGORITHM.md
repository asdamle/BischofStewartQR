# Randomized Bischof–Stewart Column Selection

A description of the randomized column-selection algorithm implemented in
`matlab_rand/` (`bsqr_rand`). It selects a well-conditioned set of `k` columns
with the *same* worst-case guarantee on `‖R₁₁⁻¹‖` as the deterministic
Bischof–Stewart pivoted QR, but at a cost that is essentially independent of the
number of candidate columns `n`. The first part is a prose description; the second
is a LaTeX draft of the same algorithm using the `algorithm` / `algpseudocode`
packages.

---

## 1. Problem and guarantee

Given `A ∈ ℝ^{m×n}` and a target `k ≤ min(m,n)`, choose an index set `J` of `k`
columns. Let `A(:,J) = Q [R₁₁; 0]` be the (unpivoted) QR factorization of the
selected columns, with `R₁₁ ∈ ℝ^{k×k}` upper triangular. The quality of the
selection is measured by how small `‖R₁₁⁻¹‖` is — equivalently, how far the
selected columns are from being linearly dependent.

The canonical setting is **orthonormal-row input**: `A` is `k×n` with `A Aᵀ = I_k`
(for example, `A = Vₖᵀ`, the leading-`k` right singular vectors of some matrix to
be column-subset-selected). There, a Bischof–Stewart selection achieves Osinsky's
bounds

```
‖R₁₁⁻¹‖_F ≤ √( k (n − k + 1) ),        ‖R₁₁⁻¹‖₂ ≤ √( 1 + k (n − k) ).
```

The randomized algorithm maintains these same bounds.

## 2. The deterministic criterion (background)

Process columns one at a time. Suppose `i` columns have been selected, with
triangular factor `R₁₁ ∈ ℝ^{i×i}` and accumulated orthogonal `Qᵢ` (a product of
`i` Householder reflectors). Bring a candidate column `aⱼ` into the current frame,
`x = Qᵢᵀ aⱼ`, and split it into a **head** `x_{1:i}` and a **tail** `x_{i+1:m}`.
Define

```
ρⱼ²  = ‖x_{i+1:m}‖²            (the squared residual; ρⱼ is the new diagonal of R₁₁)
wⱼ   = R₁₁⁻¹ x_{1:i}           (how the candidate is expressed in the current basis)
cⱼ   = (1 + ‖wⱼ‖²) / ρⱼ².      (the Bischof–Stewart criterion)
```

The criterion is exactly the **growth of the squared inverse Frobenius norm**.
Writing `f = ‖R₁₁⁻¹‖_F²`, the rank-one bordered-inverse formula gives, when column
`j` is appended,

```
f_{i+1} = f_i + cⱼ.            (exact — Eq. (growth))
```

So the criterion penalizes a candidate that is nearly in the span of the already
selected columns (`‖wⱼ‖` large) or that has a small residual (`ρⱼ` small): either
inflates `‖R₁₁⁻¹‖_F`. The deterministic kernel selects `argminⱼ cⱼ` at every step.
Computing it requires `R₁₁⁻¹R₁₂` and the running residual norms for **all** trailing
columns at every step — the `O(n k²)` work that dominates the short-wide (`k ≪ n`)
regime.

## 3. Randomization: sample, don't scan

The guarantee does not need the exact `argmin`. It only needs each step's increment
`cⱼ` to keep the running `f` under a per-step threshold. Define the
**running-mean threshold**

```
θᵢ = (f_i + n − 2i) / (k − i).
```

For orthonormal-row input this is precisely the `ρ²`-weighted mean of the criterion
over the remaining columns (using `Σ ρⱼ² = k − i` and `Σ ‖wⱼ‖² = f_i − i`). Two
consequences follow:

1. **Feasibility.** Since the minimum is at most the mean, at least one remaining
   column satisfies `cⱼ ≤ θᵢ` — so accepting *any* column that meets the threshold
   never stalls the algorithm.
2. **The bound.** Enforcing `cⱼ ≤ θᵢ` at every step telescopes (`f_{i+1} ≤
   ((k−i+1)/(k−i)) f_i + (n−2i)/(k−i)`) to `f_k ≤ k(n−k+1)`, i.e. Osinsky's bound,
   together with its per-singular-value refinement.

Because a constant fraction of columns clears the mean, we can **sample** a small
block of candidates, bring only those up to date, and accept one that meets the
threshold — never touching the columns we do not sample. This removes the `O(n)`
per-step scan; the per-step work depends on the block size, not on `n`.

Candidates are drawn **without replacement, weighted by their starting squared
column norms** `gⱼ = ‖aⱼ‖²` (robust across leverage profiles), or uniformly when
leverage is known to be flat. Weighted sampling uses Efraimidis–Spirakis keys
`−log(uⱼ)/gⱼ`, `uⱼ ∼ Unif(0,1)`, taking the smallest keys.

## 4. The algorithm (batched in-block selection — the default)

Sampling a block and bringing it into the current frame costs one application of
the accumulated reflectors — the dominant per-block cost. To amortize it, once a
block is reduced we run Bischof–Stewart **within the block**: repeatedly take the
in-block minimizer, append it, update `f` and re-derive `θ`, and downdate the rest
of the block by the new reflector — continuing while the in-block minimizer still
clears the (now updated) threshold. One expensive apply thus yields several
selections, cutting the dominant cost from `O(k⁴)` to `O(k³)`.

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

**Single-select variant (`batched = false`).** Instead of selecting repeatedly
inside a block, scan sampled blocks in weighted order until a column with `c ≤ θ` is
found (best-in-block or first-fit), accept that one column, and resample for the
next step. This makes one reflector apply per *selection* (the `O(k⁴)` path) but
yields realized `‖R₁₁⁻¹‖_F` slightly closer to the deterministic value; the
guarantee is identical.

**Threshold variants.** `running_mean` (above) is the default. A more permissive
`worstcase_allowance` mode sets `θᵢ = F̂_{i+1} − f_i` with `F̂_m = m(n−m+1)/(k−m+1)`,
spending the slack between the *actual* running `f_i` and the deterministic worst
case `F̂_i`; it accepts with fewer samples and still keeps `f_k ≤ k(n−k+1)`. A
`slack ≥ 1` multiplier loosens either threshold (and the bound) proportionally.

**Safeguards.**
- *Feasibility net.* If many columns have been sampled since the last selection
  without success (only near-rounding ties), force one pass over **all** remaining
  columns; the global minimizer is then guaranteed to clear the threshold in exact
  arithmetic.
- *Rank guard.* If every remaining column has `ρ² ≈ 0` (so `c = ∞`), the input is
  numerically rank-deficient for this `k`; the algorithm stops with an error rather
  than propagating `∞` into `f` and `R₁₁`.

**Outputs.** The permutation with `J` first, the accumulated Householder reflectors
(an implicit `Q`), `R₁₁`, instrumentation, and — on request, via one extra
`O(n k²)` apply — `R₁₂`.

## 5. Cost

| | deterministic BSQR | randomized (this algorithm) |
|---|---|---|
| per step | `O(n k)` (scan all columns) | `O(b·k)` candidates, `b = O(1)` blocks |
| total | `O(n k²)` | `O(k³)` + `O(mn)` (norm precompute) + `O(n)` |

For `k ≪ n` the randomized cost is essentially independent of `n`, so the speedup
grows linearly in `n` while `‖R₁₁⁻¹‖` stays under the same guarantee.

---

## 6. LaTeX draft

Compile with `\usepackage{algorithm}`, `\usepackage{algpseudocode}`, and
`\usepackage{amsmath,amssymb}` (the latter for `\mathbb`). The growth identity
referenced in the main routine is

```latex
\begin{equation}\label{eq:growth}
  \bigl\|R_{11}^{-1}\bigr\|_F^2 \;\text{increases by exactly}\;
  c_j = \frac{1 + \|w_j\|^2}{\rho_j^2}
  \quad\text{when column $j$ is appended.}
\end{equation}
```

```latex
\begin{algorithm}[t]
\caption{Randomized Bischof--Stewart column selection (batched; default path)}
\label{alg:randbsqr}
\begin{algorithmic}[1]
\Require $A\in\mathbb{R}^{m\times n}$; rank $k\le\min(m,n)$; block size $b$ (default $b=k$);
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
           \Comment{per-step threshold $=$ weighted mean of $c$}
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
> accepts the first sampled column with `$c\le\theta$` and then resamples gives the
> `batched=false` path: one `\Call{Reduce}{}` per selected column rather than per
> block. The outer structure, threshold, and guarantee are unchanged.
