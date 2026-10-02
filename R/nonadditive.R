# The non-additive covariance matrices of Mrode & Pocrnic (2023, 4th ed.), chapter 13. None of
# this is a fitter: each function builds a K that enters the model through the generic
# declared-covariance marker, kernel(id, K = ...), and the engine treats it exactly like
# A or H, the same kron(C, K^-1) penalty, the same score, the same AI.
#
# The split of labor is deliberate. The book publishes the MATRICES (D on p.227, the
# genomic G and D on p.233, G_AA on p.238), so building them in R keeps them inspectable
# and testable against the printed pages; the inversion and everything that follows lives
# in the engine, under the gates the rest of the package already answers to.

#' Dominance relationship matrix from a pedigree
#'
#' The D of Cockerham (1954), Eqn 13.1 of Mrode & Pocrnic (4th ed., p.226): between
#' animal x with parents s, d and animal y with parents f, m,
#' `d_xy = 0.25 (a_sf a_dm + a_sd a_fm)`, with a the additive relationship and the
#' diagonal equal to 1. The formula assumes a non-inbred population; with inbreeding in
#' the pedigree the off-diagonals are the classic approximation, not an exact identity.
#'
#' The matrix is DENSE and is built from the full tabular A, so this is for pedigrees of
#' moderate size. For a large pedigree use [dominance_inverse()], the sparse inverse of
#' the same D by sire x dam subclasses (Hoeschele & VanRaden, 1991, generalized to
#' inbreeding and overlapping generations), which enters the model as
#' `kernel(id, Kinv = dominance_inverse(ped))` with the same likelihood as
#' `kernel(id, K = dominance_matrix(ped))` whenever that representation exists (see
#' its page).
#'
#' @param ped data.frame with animal, sire and dam, as in [model()]; unknown parent 0 or NA
#' @return dominance relationship matrix, rows and columns named by animal in the
#'   topological order of [pedigree()], ready for `kernel(id, K = )`
#' @examples
#' ped <- data.frame(animal = c("1", "2", "3", "4"),
#'                   sire   = c("0", "0", "1", "1"),
#'                   dam    = c("0", "0", "2", "2"))
#' dominance_matrix(ped)["3", "4"]   # full sibs: 0.25
#' @export
dominance_matrix <- function(ped) {
  recusa_mgs(ped, "dominance_matrix()")
  p <- pedigree(ped)
  n <- nrow(p)
  s <- p$sire
  d <- p$dam
  A <- matrix(0, n, n)
  for (i in seq_len(n)) {
    if (i > 1L) for (j in seq_len(i - 1L)) {
      A[i, j] <- A[j, i] <- 0.5 * ((if (!is.na(s[i])) A[j, s[i]] else 0) +
                                   (if (!is.na(d[i])) A[j, d[i]] else 0))
    }
    A[i, i] <- 1 + if (!is.na(s[i]) && !is.na(d[i])) 0.5 * A[s[i], d[i]] else 0
  }
  aij <- function(a, b) if (!is.na(a) && !is.na(b)) A[a, b] else 0
  D <- diag(1, n)
  for (i in seq_len(n)) if (i > 1L) for (j in seq_len(i - 1L)) {
    D[i, j] <- D[j, i] <- 0.25 * (aij(s[i], s[j]) * aij(d[i], d[j]) +
                                  aij(s[i], d[j]) * aij(d[i], s[j]))
  }
  dimnames(D) <- list(p$id, p$id)
  D
}

#' Sparse inverse of the dominance relationship matrix, by sire x dam subclasses
#'
#' The method of Hoeschele & VanRaden (1991), generalized to inbreeding, overlapping
#' generations, unknown parents and repeated matings. The dominance deviation of an animal
#' with both parents known is split into the effect of its sire x dam subclass and a
#' deviation within it, `d_i = h_SD + delta_i`, with `Var(h) = (F / 4) s2d`,
#' `F_{SD,XY} = a_SX a_DY + a_SY a_DX`, and `Var(delta_i) = Delta_i s2d`,
#' `Delta_i = 1 - [(1 + F_S)(1 + F_D) + 4 F_i^2] / 4`; an animal without both parents has
#' no subclass and `Delta_i = 1`. Then `D = W (F / 4) W' + diag(Delta)` is the Cockerham D
#' of [dominance_matrix()] as an identity, inbreeding included, and when every
#' `Delta_i > 0` (see below) the precision of the augmented vector `[d; h]`,
#' `Q = [[Delta^-1, -Delta^-1 W], [-W' Delta^-1, 4 F^-1 + W' Delta^-1 W]]`, exists and is
#' sparse. With `kernel(id, Kinv = dominance_inverse(ped))` the component is the same `s2d`
#' as with `kernel(id, K = dominance_matrix(ped))` and the subclass levels are latent
#' effects, like the ancestors without records in A^-1, so the likelihood is the same.
#' The tests check, at fixed components, -2logL, score and average information in
#' [model()] and [model_mt()], -2logL and score in [model_ar1()], and in [gibbs()] the
#' starting point and, by Monte Carlo, the posterior mean of the effects against the BLUP.
#' The REML optimum is the same; the path to it need not be (the EM warm-up of [model()]
#' moves differently with the latent levels).
#'
#' Two exact routes build `F^-1`. `"dense"` forms F of the subclasses from the additive
#' relationships between the parents (Colleau) and inverts it, which costs `n^3 / 3` in
#' the number of subclasses: it wins with few subclasses of many offspring (litters).
#' `"sparse"` uses the recursion of the paper on pairs of animals, which is the row of
#' `(I - P) (x) (I - P)` with residual variance `d_x d_y` (`2 d_x^2` for the own pair of an
#' animal), `d_x` the Mendelian coefficient of A, over the ancestral closure of the
#' subclasses, and integrates out exactly every pair without animals that has at most one
#' child: it wins when each subclass has few offspring (dairy). The paper derives the
#' recursion for non-inbred animals and prunes pairs by a rule it does not prove; the
#' residual variances under inbreeding and the pruning rule here are this package's
#' generalization, checked against [dominance_matrix()] on inbred pedigrees (sire x
#' daughter, self-fertilization, reciprocal crosses, animals used as sire and as dam,
#' overlapping generations) and not taken from the paper.
#'
#' `"auto"` compares the symbolic Cholesky cost of the joint mixed-model equations of a
#' reference model, the additive animal effect with A^-1 plus this precision with one
#' record per animal, and takes the cheaper. The cost is the sum of squared column counts
#' of the factor, with the sparse columns weighted 8 times the dense final block, which
#' the package factors in register-blocked tiles; the weight is a fixed constant, so the
#' route, and with it the levels of the term, depend neither on the machine nor on the
#' number of threads. The reference model is a proxy: fixed effects, permanent
#' environment and repeated records are not in it. The dense candidate is costed first
#' and bounds the work spent on the sparse one, which is dropped for the dense route
#' (`decision` says which exit) when its pair closure holds more pairs than the dense
#' block has entries or more memory than the dense route needs at its peak; when its pair
#' block of Q would be assembled from more triplets than the full triangle of the dense
#' block (`entries`: both routes assemble Q in triplets, and A^-1 and the animal part,
#' common to both, stay out of the count); or when the minimum-degree ordering of its
#' equations reads more list entries than the ordering of the dense candidate plus 0.05
#' times the dense cost (`work`): an ordering of a candidate with heavy fill is far slower
#' per operation than the tiled factorization, and the fit would pay that same ordering
#' before its first iteration. Every exit is a deterministic count, never a timing. The
#' choice is skipped when the dense route does not fit, which leaves the sparse one, and
#' when `route =` is given.
#'
#' `validation/dominance_hv91_routes.R` times the choice and both routes (build plus one
#' evaluation of the mixed-model equations at fixed components) on eight pedigrees of
#' 3,230 to 60,100 animals, from litters to one offspring per pair, with three replicates
#' except two single runs that its output marks. The decision took 0.02 to 1.9 s. In five
#' of the seven pedigrees where both routes were timed, the route chosen was the cheaper one
#' for a fit (build plus ten evaluations): with litters the sparse route cost 2 times more
#' (litters of 8, four generations) and 140 times more (litters of 6, six generations, 9,280
#' animals: 4.0 s against 572 s, the sparse route a single run), with two offspring per pair
#' 1.3 times, and with one offspring per pair the dense route cost 11 and 5.6 times more
#' (3,230 and 15,040 animals; on the latter 233 s against 42 s, the dense route a single
#' run). On 25,100 animals with 300 records in distinct subclasses the two cost the same
#' (2.1 and 2.2 s) and `"auto"` took the sparse one. On 60,100 animals with 800 records
#' after five generations the ordering of the sparse candidate passed the work ceiling and
#' `"auto"` took the dense route, which cost 1.9 times the sparse one (18.1 s against
#' 9.5 s). In the litter pedigree of 12,280 animals the sparse route is not timed: run once
#' on its own, it was still building Q after 608 s, past 2 GB of memory, when it was
#' stopped (forcing `route = "sparse"` on large litters can cost that much); `"auto"` takes
#' the dense one there. It is a heuristic over a proxy model, not a guarantee: where the
#' two routes cost about the same it may take either, and the work ceiling can send a
#' cheaper sparse route to the dense one, as in the 60,100-animal pedigree.
#'
#' The dense route needs at most about `18 n_sub^2 + 8 n_parents^2` bytes at its peak
#' (`memory`), and is not offered above 4 GB. `validation/dominance_hv91_memory.R`
#' measures the rise of a fresh R process above where it stood, three runs each: 72 MB for
#' the dense route with 2,400 subclasses, against 119 MB estimated; and in two litter
#' pedigrees where `"auto"` picks the dense route, 79 against 50 MB on 14,730 animals and
#' 112 against 88 MB on 19,530, so the decision adds 24 to 29 MB and stays under the
#' estimate of the dense peak (85 and 152 MB).
#'
#' A fit holds the precision in R (the triplets, 16 bytes per entry), as the penalty of the
#' term in the engine (12 bytes per entry; the engine's input copy is released once the
#' equations are set up), and builds the factor and the selective inverse of the
#' mixed-model equations on top. One evaluation at fixed components with 2,000 subclasses
#' and 1.8 million entries in Q rose 209 MB, about 116 bytes per entry of Q (the same
#' script). The 4 GB ceiling above guards only the construction; check the memory before
#' fitting a precision with tens of millions of entries. The long steps (the closure, the
#' orderings, the factorizations) check for an interrupt between columns or panels: a time
#' limit set with `setTimeLimit()` stops them, which the tests check; Ctrl-C goes through
#' the same check and was not tested interactively.
#'
#' WHEN IT DOES NOT EXIST. `Delta_i <= 0` leaves nothing for the within-subclass deviation,
#' and the function stops, naming the animal. With self-fertilization it happens in the
#' second generation (`Delta = -0.125`), with full-sib mating in the sixth (parents with
#' F = 0.59). Two cases. When the subclass has two or more animals in Q, D itself is not
#' positive-definite on them (two full sibs have `D_ij = 1 - Delta`, an eigenvalue of
#' `Delta`), so no route fits it. When the animal is alone in its subclass the Cockerham D
#' may still be positive-definite (the selfing example is, smallest eigenvalue 0.5) and
#' only the subclass representation is missing: if that animal has no record, leave it out
#' with `animals =` (its Delta enters nothing else); otherwise use
#' `kernel(id, K = dominance_matrix(ped))`.
#'
#' WHAT IT DOES NOT FIX. The Cockerham D with unit diagonal is the classic approximation
#' under inbreeding (the full dominance covariance needs further identity components); this
#' is the inverse of THAT matrix, with the same approximation as the dense route. A sire /
#' maternal-grandsire pedigree has no sire x dam subclass and is refused; a parent code
#' with no line of its own (a metafounder or unknown-parent group label) is refused as in
#' [pedigree()], since the Cockerham D is defined on a base of unrelated, non-inbred
#' founders.
#'
#' @param ped data.frame with animal, sire and dam, as in [model()]; unknown parent 0 or NA
#' @param animals the animals that get a dominance level, by default all of the pedigree.
#'   Pass the animals with records: the marginal of D on them is exact and the precision
#'   carries only their subclasses. A fit refuses a record with an observed response of an
#'   animal left out here, which would otherwise be dropped from the analysis without
#'   notice
#' @param route `"auto"`, `"dense"` or `"sparse"`, see above; all give the same marginal
#'   on the animals
#' @param max_pairs ceiling on the pairs the sparse closure may create, a guard on memory:
#'   a pair costs a few hundred bytes while the closure is built (an estimate from the
#'   layout, about 320, not measured), so the default allows on the order of 1.5 GB. Values
#'   above `2^32 - 1`, the most pairs the closure can index, act as `2^32 - 1`
#' @return the lower triangle of `Q` in triplets, `list(i, j, x, n, id)` as [a_inverse()]
#'   returns, plus: `type`, per level, `"animal"`, `"subclass"` (the sire x dam pair of an
#'   animal in Q) or `"ancestral"` (a pair of the sparse closure kept as a latent level);
#'   `prior`, the prior variance of each level at unit `s2d` (1 for an animal, `F_cc / 4`
#'   for a pair), which [accuracy()] divides by; `logdet`, `log|K|` of the augmented
#'   covariance by the closed form, which a fit uses instead of factoring Q again while
#'   the triplets and `logdet` still match `signature` (any change to either and the fit
#'   factors Q); `start_scale`, the geometric mean of the
#'   eigenvalues of `D[animals, animals]` (of D itself without `animals =`), which the
#'   multi-trait, AR(1) and Gibbs fitters use for their starting values so that they start
#'   where `kernel(K = D[animals, animals])` starts; `route`, `"dense"`, `"sparse"`, or
#'   `"none"` when no animal of Q has both parents (then `Q = I`); `decision`, why that
#'   route: `"given"` (`route =`), `"compared"` (both costed), `"closure"`, `"memory"`,
#'   `"ordering_work"` and `"ordering_cost"` (the early exits above, all to the dense
#'   route), `"dense_too_big"` (the sparse route is the only one) or `"no_subclass"`;
#'   `n_subclasses`; `n_pairs`;
#'   `pairs`, the two animals of each pair level in pedigree order, one row per level named
#'   by its label (a level `h` of pair x, y has `Var(h) = (a_xx a_yy + a_xy^2) / 4` and
#'   `Cov(h, h') = (a_xx' a_yy' + a_xy' a_yx') / 4` at unit `s2d`);
#'   `closure` (pairs created and processed, and whether it was abandoned); `flops`, the
#'   sum of squared column counts of the Cholesky factor of the reference model under each
#'   route, and `cost`, the same with the sparse part of the factor weighted 8 times the
#'   dense final block, which is what `"auto"` compares (NA when not computed, and the
#'   sparse ones `Inf` when its ordering was abandoned); `work`, the list entries read by
#'   the ordering of the sparse candidate (`ordering`, up to where it stopped) and by that
#'   of the dense one (`dense`), and the `ceiling` of the former, the latter plus 0.05
#'   times the dense cost; `memory`, the estimated peak bytes of building Q by the dense
#'   route; `entries`, the triplets each route assembles for the pair block of Q, the full
#'   triangle for the dense one and, counted before assembly, the sparse candidate's (NA
#'   when not counted), which the memory exit compares; `seconds`, the wall time of the
#'   pair closure, of the choice of route and of building Q; `signature`, a hash of the
#'   triplets and of `logdet`; and `left_out`, the pedigree animals without a level. A
#'   subclass is labelled `"sire x dam"` (the parents of its first animal in pedigree
#'   order) and an ancestral pair `"a x b"` in pedigree order. In [ebv()] the `"animal"`
#'   levels are the dominance deviations d-hat and the pair levels the h-hat;
#'   `k_level_type` in the fit and the `type` column of [solutions()] tell them apart. An
#'   animal of the pedigree left out of Q has no d-hat; its expected dominance is the
#'   h-hat of its sire x dam pair when that pair is a level (a full sib is in Q, or the
#'   sparse route kept the pair), and is not reported otherwise.
#' @references Hoeschele, I. & VanRaden, P.M. (1991). Rapid inversion of dominance
#'   relationship matrices for noninbred populations by including sire by dam subclass
#'   effects. Journal of Dairy Science 74:557-569.
#'
#'   Cockerham, C.C. (1954). An extension of the concept of partitioning hereditary
#'   variance for analysis of covariances among relatives when epistasis is present.
#'   Genetics 39:859-882.
#' @examples
#' ped <- data.frame(animal = as.character(1:8),
#'                   sire   = c("0", "0", "0", "1", "1", "3", "3", "4"),
#'                   dam    = c("0", "0", "0", "2", "2", "2", "2", "6"))
#' Dinv <- dominance_inverse(ped)
#' table(Dinv$type)
#' Q <- matrix(0, Dinv$n, Dinv$n)
#' Q[cbind(Dinv$i, Dinv$j)] <- Dinv$x
#' Q[cbind(Dinv$j, Dinv$i)] <- Dinv$x
#' # the animal block of Q^-1 is the dense D
#' max(abs(solve(Q)[1:8, 1:8] - dominance_matrix(ped)))
#' @export
dominance_inverse <- function(ped, animals = NULL, route = c("auto", "dense", "sparse"),
                              max_pairs = 5e6) {
  if (!is.data.frame(ped)) stop("expected a pedigree data.frame")
  recusa_mgs(ped, "dominance_inverse()")
  route <- match.arg(route)
  if (!is.numeric(max_pairs) || length(max_pairs) != 1L || !is.finite(max_pairs) ||
      max_pairs < 1)
    stop("max_pairs must be a single positive number")
  cp <- colunas_pedigree(ped)
  quem <- character(0)
  if (!is.null(animals)) {
    quem <- unique(rotulo_motor(animals))
    if (!length(quem) || anyNA(quem))
      stop("animals = must be a vector of animal ids, without NA")
  }
  r <- .Call(R_dominancia_inversa, cp$id, cp$sire, cp$dam, quem,
             match(route, c("auto", "dense", "sparse")) - 1L, as.double(max_pairs))
  pid <- r$pedigree_id
  an <- pid[r$animal_row]
  # o rotulo de uma subclasse e "pai x mae" do primeiro animal dela em ordem do pedigree; o
  # de um par ancestral, os dois animais em ordem do pedigree
  px <- pid[r$pair_x]
  py <- pid[r$pair_y]
  ia <- match(an, cp$id)
  sa <- cp$sire[ia]
  da <- cp$dam[ia]
  canon <- function(a, b) ifelse(a <= b, paste(a, b, sep = "\r"), paste(b, a, sep = "\r"))
  primeiro <- ifelse(r$pair_full, match(canon(px, py), canon(sa, da)), NA_integer_)
  # sem par nenhum o paste0 devolveria um " x " solto, um nivel inventado
  rot <- if (length(px)) paste0(ifelse(is.na(primeiro), px, sa[primeiro]), " x ",
                                ifelse(is.na(primeiro), py, da[primeiro])) else character(0)
  id <- c(an, rot)
  if (anyDuplicated(id) || any(rot %in% pid))
    stop("a subclass label (\"sire x dam\") coincides with an animal id or with another ",
         "label: the levels could not be told apart in the results", call. = FALSE)
  num <- r$num
  n_q <- length(an)
  list(i = r$q$i, j = r$q$j, x = r$q$x, n = r$q$n, id = id,
       type = c(rep("animal", n_q), ifelse(r$pair_full, "subclass", "ancestral")),
       prior = stats::setNames(c(rep(1, n_q), r$pair_prior), id),
       logdet = num[1], start_scale = exp(num[3] / n_q),
       route = if (num[5] == 0) "none" else c("dense", "sparse")[num[4]],
       decision = c("given", "compared", "closure", "memory", "ordering_work",
                    "ordering_cost", "dense_too_big", "no_subclass")[num[15] + 1],
       n_subclasses = num[5], n_pairs = length(rot),
       pairs = matrix(c(px, py), ncol = 2L, dimnames = list(rot, c("first", "second"))),
       closure = list(created = num[6], processed = num[7], aborted = r$aborted),
       flops = c(dense = num[8], sparse = num[9]), cost = c(dense = num[10], sparse = num[11]),
       work = c(ordering = num[16], dense = num[21], ceiling = num[17]),
       memory = c(dense = num[18]), entries = c(dense = num[19], sparse = num[20]),
       seconds = c(closure = num[12], decision = num[13], build = num[14]),
       signature = r$signature,
       left_out = if (length(quem)) setdiff(pid, an) else character(0))
}

#' Genomic relationship matrix (VanRaden)
#'
#' The G of VanRaden (2008), first method: genotypes centered by twice the allele
#' frequency, `G = Z Z' / (2 sum p q)`, frequencies from the genotyped animals
#' themselves. This is the RAW G of the book's worked examples (Mrode & Pocrnic, 4th
#' ed., Example 13.3, p.233), no blending with A22 and no affine adjustment, which are
#' the single-step steps that `model(genotypes = )` performs internally. With few
#' animals G is singular; add a small ridge before declaring it, `G + diag(0.01, n)`,
#' as the book does in its examples.
#'
#' @param genotypes list with `ids` and `m` (animals x markers, 0/1/2 coded); NA is
#'   imputed with the marker mean, and monomorphic markers are left out
#' @return genomic relationship matrix, rows and columns named by `ids`
#' @examples
#' geno <- list(ids = c("a", "b"), m = rbind(c(0, 1, 2), c(2, 1, 0)))
#' g_matrix(geno)
#' @export
g_matrix <- function(genotypes) {
  g <- valida_genotipos(genotypes)
  if (!length(g$gid)) stop("genotypes must be a list with 'ids' and 'm'")
  p <- .Call(R_freq_genotipos, g$gm) / 2
  if (!any(is.finite(p) & p > 0 & p < 1))
    stop("every marker is monomorphic: there is no G to build")
  # o mesmo nucleo do passo unico, lendo a matriz do R por blocos de marcadores, sem copia
  G <- .Call(R_g_matrix, g$gm)
  dimnames(G) <- list(g$gid, g$gid)
  G
}

#' Genomic dominance relationship matrix (Vitezica)
#'
#' The D of Vitezica et al. (2013), the parametrization the book adopts for its GBLUP
#' dominance model (Mrode & Pocrnic, 4th ed., Eqn 13.6, p.231, and Example 13.3, p.233):
#' each genotype 0, 1, 2 is coded `-2p^2`, `2pq`, `-2q^2` and
#' `D = W W' / sum((2 p q)^2)`, frequencies from the genotyped animals themselves. Under
#' Hardy-Weinberg equilibrium this partition is orthogonal to the additive one, which is
#' why the additive and the dominance term enter the model as separate covariance
#' groups. Like [g_matrix()], the result is raw and may need a ridge before entering
#' `kernel()`.
#'
#' @param genotypes list with `ids` and `m` (animals x markers, 0/1/2 coded); an NA is
#'   imputed with the mean of the observed dominance codes of its marker, and
#'   monomorphic markers are left out
#' @return genomic dominance relationship matrix, rows and columns named by `ids`
#' @examples
#' geno <- list(ids = c("a", "b"), m = rbind(c(0, 1, 2), c(2, 1, 0)))
#' g_dominance(geno)
#' @export
g_dominance <- function(genotypes) {
  g <- valida_genotipos(genotypes)
  if (!length(g$gid)) stop("genotypes must be a list with 'ids' and 'm'")
  gn <- genotipos_numericos(g$gm)
  p <- colMeans(gn, na.rm = TRUE) / 2
  ok <- is.finite(p) & p > 0 & p < 1
  if (!any(ok)) stop("every marker is monomorphic: there is no D to build")
  m <- gn[, ok, drop = FALSE]
  q <- 1 - p[ok]
  W <- matrix(0, nrow(m), ncol(m))
  for (j in seq_len(ncol(m))) {
    cod <- c(-2 * p[ok][j]^2, 2 * p[ok][j] * q[j], -2 * q[j]^2)[m[, j] + 1]
    cod[is.na(cod)] <- mean(cod, na.rm = TRUE)
    W[, j] <- cod
  }
  D <- tcrossprod(W) / sum((2 * p[ok] * q)^2)
  dimnames(D) <- list(g$gid, g$gid)
  D
}

#' Additive-by-additive epistatic relationship matrix
#'
#' The G_AA behind the book's epistatic GBLUP (Mrode & Pocrnic, 4th ed., Eqn 13.13 and
#' Example 13.5, p.237-238): the Hadamard square of the additive genomic relationship,
#' rescaled so the diagonal averages 1, `G_AA = (G * G) / mean(diag(G * G))`. Cheap on
#' purpose, once a G exists, epistasis is one elementwise product away. The book notes
#' this is the one non-additive term that visibly reordered the animals in its example.
#'
#' @param g additive relationship matrix, usually [g_matrix()]; any square symmetric
#'   relationship works
#' @return epistatic relationship matrix with the dimnames of `g`, raw like the input,
#'   ridge it before `kernel()` if `g` was singular
#' @examples
#' geno <- list(ids = c("a", "b"), m = rbind(c(0, 1, 2), c(2, 1, 0)))
#' g_epistasis(g_matrix(geno))
#' @export
g_epistasis <- function(g) {
  if (!is.matrix(g) || !is.numeric(g) || nrow(g) != ncol(g))
    stop("expected a square numeric relationship matrix, like the one from g_matrix()")
  gaa <- g * g
  gaa / mean(diag(gaa))
}

#' Mixed and higher-order epistatic relationship matrices
#'
#' The Hadamard construction of [g_epistasis()] is not restricted to additive-by-additive.
#' Under the orthogonal partition the relationship matrix of an interaction between two
#' effect types is the elementwise product of their relationship matrices, so with `G` the
#' additive and `D` the dominance matrix, `G_AD = G * D` and `G_DD = D * D`, each rescaled
#' so its diagonal averages one. Higher orders follow the same rule: `G_AAA = G * G * G`.
#'
#' Two warnings that belong with the formula rather than after it. The orthogonality that
#' makes these terms separable holds under Hardy-Weinberg and linkage equilibrium; away
#' from it the terms correlate and the partition of variance stops being clean. And every
#' order added is a covariance component the data has to support: a fit with additive,
#' dominance and three epistatic terms asks five variances of a design that usually
#' struggles to give two, and the components collapse to the boundary rather than
#' informing anything. Add an order because a hypothesis calls for it, never because the
#' function exists.
#'
#' @param g additive relationship matrix, as from [g_matrix()]
#' @param d dominance relationship matrix, as from [g_dominance()]
#' @param order for `g_epistasis_order()`, how many additive factors to multiply: 2 is
#'   [g_epistasis()], 3 is additive-by-additive-by-additive
#' @return relationship matrix with the dimnames of the input, rescaled to unit mean
#'   diagonal; ridge it before `kernel()` if the input was singular
#' @examples
#' geno <- list(ids = c("a", "b", "c"),
#'              m = rbind(c(0, 1, 2), c(2, 1, 0), c(1, 1, 1)))
#' G <- g_matrix(geno)
#' D <- g_dominance(geno)
#' g_epistasis_ad(G, D)["a", "b"]
#' @name epistasia
NULL

escala_diagonal <- function(m, quem) {
  if (!is.matrix(m) || !is.numeric(m) || nrow(m) != ncol(m))
    stop(quem, ": expected a square numeric relationship matrix")
  md <- mean(diag(m))
  if (!is.finite(md) || md == 0)
    stop(quem, ": the product has a zero or non-finite mean diagonal, so it cannot be ",
         "rescaled to unit diagonal; check the inputs for monomorphic markers")
  m / md
}

#' @rdname epistasia
#' @export
g_epistasis_ad <- function(g, d) {
  if (!is.matrix(d) || nrow(d) != nrow(g) || ncol(d) != ncol(g))
    stop("g_epistasis_ad(): G and D must have the same dimensions")
  if (!is.null(dimnames(g)) && !is.null(dimnames(d)) &&
      !identical(rownames(g), rownames(d)))
    stop("g_epistasis_ad(): G and D are indexed by different animals; reorder one of them")
  out <- escala_diagonal(g * d, "g_epistasis_ad")
  dimnames(out) <- dimnames(g)
  out
}

#' @rdname epistasia
#' @export
g_epistasis_dd <- function(d) {
  out <- escala_diagonal(d * d, "g_epistasis_dd")
  dimnames(out) <- dimnames(d)
  out
}

#' @rdname epistasia
#' @export
g_epistasis_order <- function(g, order = 2L) {
  order <- as.integer(order)
  if (length(order) != 1L || is.na(order) || order < 2L)
    stop("g_epistasis_order(): order must be an integer >= 2; order 2 is g_epistasis()")
  out <- g
  for (k in seq_len(order - 1L)) out <- out * g
  out <- escala_diagonal(out, "g_epistasis_order")
  dimnames(out) <- dimnames(g)
  out
}


#' Genomic inbreeding from marker homozygosity
#'
#' `f = 1 - h/N`, the proportion of homozygous SNPs per animal (Mrode & Pocrnic, 4th
#' ed., Eqn 13.10-13.11, p.234-235). Fitted as a fixed covariate, `cov(f)`, its
#' regression coefficient is the inbreeding depression, the number Example 13.4
#' publishes, read back from `fit$b`. This is input data, not a model: join the vector
#' to the data by animal before the fit.
#'
#' @param genotypes list with `ids` and `m` (animals x markers, 0/1/2 coded); markers
#'   with NA are left out of that animal's denominator
#' @return named vector of genomic inbreeding coefficients, one per genotyped animal
#' @examples
#' geno <- list(ids = c("a", "b"), m = rbind(c(0, 1, 2), c(1, 1, 0)))
#' genomic_inbreeding(geno)   # a: 1 - 1/3; b: 1 - 2/3
#' @export
genomic_inbreeding <- function(genotypes) {
  g <- valida_genotipos(genotypes)
  if (!length(g$gid)) stop("genotypes must be a list with 'ids' and 'm'")
  gn <- genotipos_numericos(g$gm)
  het <- rowSums(gn == 1, na.rm = TRUE)
  obs <- rowSums(!is.na(gn))
  if (any(obs == 0))
    stop("animal(s) with no observed marker: ",
         paste(g$gid[obs == 0], collapse = ", "))
  f <- 1 - het / obs
  names(f) <- g$gid
  f
}
