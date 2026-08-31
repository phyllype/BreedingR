# APY and the Mendelian residual

Theory note: why APY works and what the Mendelian residual means. Nothing here is
needed to use the package; it is the reasoning behind `apy_core=`.

## Why APY works

Single-step needs G^-1, and inverting G costs n^3: with tens of thousands of genotyped
animals it stops being a matrix and becomes a wall. APY (Misztal et al., 2014) gets
past it with an observation about the data, not about computing: **the genome of a
population carries a limited amount of independent information**. Once a set of
animals covers that variation (the core), every additional animal is almost entirely a
recombination of what the core already has, plus a small remainder of its own.

Think of a library: after the reference collection is catalogued, each new book is
mostly a recombination of chapters already seen. You do not catalogue the new book at
the same depth; you note how it maps onto the references and the little it brings that
is new. That little is the Mendelian residual.

The population reading of "limited information": the genome is transmitted in a bounded
number of effectively independent chromosome segments, which grows with the effective
population size Ne and the genome length. A closed, selected line has small Ne, hence a
modest number of segments. Pocrnic and colleagues showed the optimal core size tracks
the number of eigenvalues explaining about 98% of the variance of G: **the right core
has the size of the realized genomic dimensionality of the population.** Note the unit
of that rule: it counts eigenvectors, and a core of n animals always captures less
variance than the top n eigenvectors, so the eigenvalue count is a floor for the
animal count at equal coverage, not an equivalence.

## The statistics, tied to Henderson

For an animal i outside the core:

```
u_i = g_iC . G_CC^-1 . u_C + phi_i        Var(phi_i) = m_i = g_ii - g_iC G_CC^-1 g_Ci
```

The breeding value is the conditional expectation given the core, plus an own
deviation. Compare with the pedigree recursion:

```
u_i = 1/2 u_sire + 1/2 u_dam + phi_i      Var(phi_i) ~ 1/2 sigma2_a
```

Same structure. In matrix form the recursion is `u = T phi`, so `A = T D T'` with T
triangular and D diagonal: an LDL' (square-root-free Cholesky) factorization of A in
pedigree order, written analytically by the biology instead of computed numerically.
The famous A^-1 rules exist because T^-1 is trivially sparse (1 on the diagonal, -1/2
at the parents), and that sparsity exists because of the Markov property of
inheritance: conditional on its two parents, an animal is independent of every other
ancestor.

The pedigree hands over three things that G does not:

1. **Who explains each animal is known a priori, and it is two animals.** The factor
   row has known positions (the parents' columns) and known coefficients (1/2 and 1/2,
   the expectation of meiosis, not estimates). In G there is no such pair: the
   regression vector `g_iC G_CC^-1` is dense over the whole core and has to be
   computed.
2. **How much is left over is a formula, not a computation.** The Mendelian sampling
   variance is analytic (1/2 sigma2_a, adjusted for parental inbreeding). Its genomic
   analogue m_i must be computed from the matrix itself.
3. **The sparsity is exact.** Henderson's A^-1 is the true inverse. The APY G^-1
   imposes conditional independence of the non-core residuals given the core; that is
   the approximation, and it is why APY warrants validation and A^-1 does not.

APY is the attempt to recreate in G the luck the pedigree gives for free in A: it has
to *choose* the universal parents (the core), *compute* the coefficients, *compute*
the remainder, and *validate* the result. Four steps that cost zero in A because
meiosis already answered them.

## The Mendelian residual, piece by piece

```
m_i = g_ii - g_iC . G_CC^-1 . g_Ci  =  Var(u_i | u_core)
```

`g_ii` is the animal's own diagonal of G (its total genetic variance, about 1+F in
sigma2_a units). `g_iC` are its genomic relationships with the core. The quadratic
form is the variance of the **best possible prediction** of the animal from the core.
So m_i is a subtraction of variances, total minus explained, which makes it a
conditional variance: the uncertainty about the animal that remains even knowing the
whole core perfectly. `m_i / g_ii` is the fraction of the animal that the panel cannot
explain.

The ruler, with concrete landmarks (sigma2_a units, no inbreeding):

| situation of animal i                 | m_i   |
|---------------------------------------|-------|
| clone/identical twin of a core animal | ~0    |
| both parents in the core              | ~0.5  |
| one parent in the core                | ~0.75 |
| unrelated to the core                 | ~1    |

The 0.5 landmark is what gives the quantity its name: with both parents in the core
the best prediction is the parent average, and what remains is **exactly the Mendelian
sampling** of classical theory, the lottery of meiosis. Even knowing the parents
completely, the offspring draws one of 2^n possible gamete combinations, and half of
the genetic variance sits in that draw. m_i generalizes the idea: when the explainer
is not the parents but an arbitrary panel, it measures the animal's genomic novelty
relative to that panel.

Two jobs inside APY:

1. **In the inverse.** The m_i form the diagonal M_NN, and each non-core animal enters
   weighted by 1/m_i: an animal that is almost fully explained carries almost no
   independent information. This is also where the package's declared error lives:
   m_i ~ 0 (a clone among the non-core) blows up 1/m_i, and the fit stops with
   "degenerate Mendelian residual in APY" instead of returning garbage.
2. **As a diversity meter.** In a sequential construction, `sum(log m_i) = log|G|`:
   the sum of the log-novelties is the total genomic variance captured. When new
   cohorts arrive with rapidly shrinking m_i, selection is narrowing the line, and the
   funnel has a number.

## Practical note on core choice

Core composition matters little when the size is adequate: random cores at the 98 to
99% eigenvalue size perform at the level of the exact inverse, a result established at
scale by Fragomeni et al. (2015, J. Dairy Sci. 98:4090, genomic recursions on 6.9M
Holsteins) and confirmed across many core definitions since (Genetics Selection
Evolution, 2022). The analogous result in numerical linear algebra: APY is a
Nystrom-type approximation, core choice is landmark selection, and randomized pivoting
is near-optimal for column Nystrom approximations. The effective lever is the size,
and the size question is answered by the eigenvalue spectrum of the population's own
G, not by how the members are picked.

## The general form: per-animal conditioning

APY conditions every non-core animal on one global core. The general form of the same
recursion gives each animal its OWN conditioning set: in a chosen order, animal i is
regressed on the k previous animals with the strongest relationship to it,

    b = G[c,c]^-1 G[c,i],    d_i = g_ii - G[c,i]' b,

the column of the sparse factor is U[c,i] = -b / sqrt(d_i), U[i,i] = 1 / sqrt(d_i),
and G^-1 is approximated by U U'. This is the Vecchia approximation of spatial
statistics, and two special cases place it firmly in this field:

1. **Henderson's A^-1 is Vecchia with the parents as conditioning set.** In a pedigree
   the conditional distribution of a breeding value given ALL previous animals depends
   only on sire and dam (the pedigree is Markovian), so conditioning on those two is
   not an approximation: run the recursion on the dense A of a pedigree without full
   sibs with k = 2 and the sparse A^-1 of the textbook comes out exactly. The famous
   rules (diagonal 1/d, -1/2d on parent-progeny, 1/4d between mates) are the outer
   products U U' of this recursion.
2. **APY is Vecchia with one shared neighborhood.** Every young animal conditions on
   the same core; the Mendelian residual m_i of APY is this d_i.

Between the two extremes, per-animal neighborhoods use the budget where each animal
needs it. Schafer, Katzfuss and Owhadi (2021, SIAM J. Sci. Comput.) prove that, given
the sparsity pattern, this factor minimizes the Kullback-Leibler divergence to the true
G^-1, and that enlarging the conditioning sets never increases it — which is why the
package's gate on nested neighborhoods (error shrinking as k grows, exactness at
k = n - 1) is a theorem check, not a hope.

The order matters: conditioning on the past is what the recursion means, so rows come
ancestors first — in genotype files sorted by birth date or by a renumbered pedigree,
that is the order the data already has. In the package: `vecchia_inverse()` for the
standalone factor and diagnostics, `vecchia_k =` in the fitters for the single step,
mutually exclusive with `apy_core =`.
