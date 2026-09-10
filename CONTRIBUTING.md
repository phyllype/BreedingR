# Contributing

## Questions and suggestions

Write to felipeoliveirafreitas@usp.br. Questions about how to specify a model, about what
a number means, or about a model the package does not fit yet are all welcome, and so are
suggestions.

If something looks like a defect, an issue is the better place, because the answer stays
visible to whoever runs into the same thing next:

https://github.com/phyllype/BreedingR/issues

## Reporting a problem

What makes a report easy to act on:

- a small script that reproduces it. `simulate_breeding()` ships with the package and
  produces a pedigree with records, so most reports can be written on simulated data.
  Please do not send data you are not free to share.
- the output of `sessionInfo()` and of `br_version()`.
- what you expected and what came out. For a fit, `summary()` of the object together with
  `converged`, `theta` and the iteration count usually says enough.

Convergence deserves a word of its own. A fit reporting `converged = FALSE` is not
necessarily a defect: a component resting on its floor, a correlation at 1, or a
covariance group the data does not identify are all legitimate places to stop. Sending
the fit along with the design, how many animals, how many groups, how many records per
animal, is usually enough to tell which of the two it is.

## Sending a change

Fork, branch, open a pull request. Two things are looked at before a merge:

- `R CMD check` passes. GitHub Actions runs it on macOS, on Windows and on three versions
  of R under Ubuntu, for every push and every pull request.
- new behaviour arrives with a test under `tests/testthat/`. The bar is that the test
  fails without the change; one that would pass either way does not gate anything.

The numerics are C++ under `src/`, compiled by `R CMD INSTALL`, and the package has no
run-time dependency. A change that would add one is a design question rather than a
detail, so it is worth raising in an issue before the code gets written.

## Building from a clone

```r
install.packages(".", repos = NULL, type = "source")   # Windows needs Rtools
```
