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

## Copyright on what you contribute

You keep the copyright on what you write. By opening a pull request you also grant the
maintainer the right to license your contribution under the licence this project carries
now, GPL-3, and under any other licence the project may adopt later.

The reason is worth stating plainly, because a clause like this is easy to resent. A
project whose copyright is spread across many authors cannot change its licence without
tracking down every one of them and getting each to agree; projects have been frozen that
way, sometimes by a single contributor nobody could reach years later. This keeps that
door open. It takes nothing away from you: your name stays on the commit, and whatever
version your work was released in stays under the licence it was released under, for
everyone who received it. Relicensing is never retroactive.

If you would rather not grant that, say so in the pull request. The contribution can
still be discussed, and often the same result is reached by describing the problem well
enough that the fix is written here.

## Building from a clone

```r
install.packages(".", repos = NULL, type = "source")   # Windows needs Rtools
```
