# References

The published work this package implements, followed by material read for methods that
are not implemented yet.

Aguilar, I., Misztal, I., Johnson, D.L., Legarra, A., Tsuruta, S. & Lawlor, T.J. (2010).
A unified approach to utilize phenotypic, full pedigree, and genomic information for
genetic evaluation of Holstein final score. *Journal of Dairy Science* 93:743-752.

Anderson, E., Bai, Z., Bischof, C., Blackford, S., Demmel, J., Dongarra, J., Du Croz,
J., Greenbaum, A., Hammarling, S., McKenney, A. & Sorensen, D. (1999). *LAPACK Users'
Guide*, 3rd ed. SIAM, Philadelphia.

Bijma, P., Muir, W.M. & Van Arendonk, J.A.M. (2007). Multilevel selection 1:
quantitative genetics of inheritance and response to selection. *Genetics* 175:277-288.

Christensen, O.F. & Lund, M.S. (2010). Genomic prediction when some animals are not
genotyped. *Genetics Selection Evolution* 42:2.

Fragomeni, B.O., Lourenco, D.A.L., Tsuruta, S., Masuda, Y., Aguilar, I., Legarra, A.,
Lawlor, T.J. & Misztal, I. (2015). Use of genomic recursions in single-step genomic best
linear unbiased predictor with a large number of genotypes. *Journal of Dairy Science*
98:4090-4094.

George, A. & Liu, J.W.H. (1989). The evolution of the minimum degree ordering algorithm.
*SIAM Review* 31:1-19.

Gilmour, A.R., Thompson, R. & Cullis, B.R. (1995). Average information REML: an
efficient algorithm for variance parameter estimation in linear mixed models.
*Biometrics* 51:1440-1450.

Henderson, C.R. (1975). Best linear unbiased estimation and prediction under a selection
model. *Biometrics* 31:423-447.

Henderson, C.R. (1976). A simple method for computing the inverse of a numerator
relationship matrix used in prediction of breeding values. *Biometrics* 32:69-83.

Kirkpatrick, M., Lofsvold, D. & Bulmer, M. (1990). Analysis of the inheritance,
selection and evolution of growth trajectories. *Genetics* 124:979-993.

Legarra, A., Christensen, O.F., Vitezica, Z.G., Aguilar, I. & Misztal, I. (2015).
Ancestral relationships using metafounders: finite ancestral populations and across
population relationships. *Genetics* 200:455-468.

Liu, Z., Goddard, M.E., Reinhardt, F. & Reents, R. (2014). A single-step genomic model
with direct estimation of marker effects. *Journal of Dairy Science* 97:5833-5850.

McQuillan, R., Leutenegger, A.-L., Abdel-Rahman, R., Franklin, C.S., Pericic, M.,
Barac-Lauc, L. et al. (2008). Runs of homozygosity in European populations. *American
Journal of Human Genetics* 83:359-372.

Meuwissen, T.H.E. & Luo, Z. (1992). Computing inbreeding coefficients in large
populations. *Genetics Selection Evolution* 24:305-313.

Misztal, I., Legarra, A. & Aguilar, I. (2014). Using recursion to compute the inverse of
the genomic relationship matrix. *Journal of Dairy Science* 97:3943-3952.

Muir, W.M. (2005). Incorporation of competitive effects in forest tree or animal
breeding programs. *Genetics* 170:1247-1259.

Muir, W.M. & Schinckel, A.P. (2002). Incorporation of competitive effects in breeding
programs to improve productivity and animal well being. *Proceedings of the 7th World
Congress on Genetics Applied to Livestock Production*, Montpellier.

National Research Council (1971). *A Guide to Environmental Research on Animals*.
National Academy of Sciences, Washington DC.

Patterson, H.D. & Thompson, R. (1971). Recovery of inter-block information when block
sizes are unequal. *Biometrika* 58:545-554.

Pocrnic, I., Lourenco, D.A.L., Masuda, Y., Legarra, A. & Misztal, I. (2016). The
dimensionality of genomic information and its effect on genomic prediction. *Genetics*
203:573-581.

Quaas, R.L. (1976). Computing the diagonal elements and inverse of a large numerator
relationship matrix. *Biometrics* 32:949-953.

Schafer, F., Katzfuss, M. & Owhadi, H. (2021). Sparse Cholesky factorization by
Kullback-Leibler minimization. *SIAM Journal on Scientific Computing* 43:A2019-A2046.

Takahashi, K., Fagan, J. & Chen, M.-S. (1973). Formation of a sparse bus impedance
matrix and its application to short circuit study. *Proceedings of the 8th PICA
Conference*, 63-69.

Vandenplas, J., Calus, M.P.L., Eding, H. & Vuik, C. (2019). A second-level diagonal
preconditioner for single-step SNPBLUP. *Genetics Selection Evolution* 51:30.

Vandenplas, J., Eding, H., Calus, M.P.L. & Vuik, C. (2018). Deflated preconditioned
conjugate gradient method for solving single-step BLUP models efficiently. *Genetics
Selection Evolution* 50:51.

Vandenplas, J., Gengler, N., Bijma, P., Misztal, I. & Legarra, A. (2022). A comprehensive
study on size and definition of the core group in the proven and young algorithm for
single-step GBLUP. *Genetics Selection Evolution* 54:34.

VanRaden, P.M. (2008). Efficient methods to compute genomic predictions. *Journal of
Dairy Science* 91:4414-4423.

Wade, K.M. & Quaas, R.L. (1993). Solutions to a system of equations involving a
first-order autoregressive process. *Journal of Dairy Science* 76:3026-3034.

Weir, B.S. & Cockerham, C.C. (1984). Estimating F-statistics for the analysis of
population structure. *Evolution* 38:1358-1370.

Willham, R.L. (1972). The role of maternal effects in animal breeding: III. Biometrical
aspects of maternal effects in animals. *Journal of Animal Science* 35:1288-1293.


## Surveyed, not implemented

- Vandenplas, J. et al. (2020) on second-level preconditioners for ssSNPBLUP,
  *Genetics Selection Evolution* (https://doi.org/10.1186/s12711-020-00543-9). — the
  upgrade path if the diagonal preconditioner ever becomes the bottleneck.
- Greedy conditional selection for Vecchia, arXiv:2307.11648 (2023); correlation-based
  selection, arXiv:2112.14591. — smarter neighborhoods than nearest-by-relationship,
  if the k needed for a target accuracy ever matters.
- Chen, Y., Epperly, E.N., Tropp, J.A. & Webber, R.J. (2024). Randomly pivoted
  Cholesky. *Communications on Pure and Applied Mathematics*
  (https://doi.org/10.1002/cpa.22234); Epperly, E.N., Tropp, J.A. & Webber, R.J.
  (2024). XTrace: making the most of every sample in stochastic trace estimation.
  *SIMAX* (arXiv:2301.07825). — randomized numerical linear algebra for kernel
  approximation and trace/log-det estimation.
- Arakawa, A. et al. (2022). Performance of the No-U-Turn sampler in multi-trait
  variance component estimation using genomic data. *Genetics Selection Evolution*
  (https://pmc.ncbi.nlm.nih.gov/articles/PMC9275044/). — the gradient-based direction
  for the Bayesian half.
- Bermann, M. et al. (2022) on reliabilities from block sparse inversion of the APY
  G^-1, and the 2024 comparison of approximation algorithms, *Journal of Animal
  Science* (https://doi.org/10.1093/jas/skab353, https://doi.org/10.1093/jas/skae195).
