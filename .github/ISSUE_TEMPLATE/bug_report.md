---
name: Something is wrong
about: A result that looks incorrect, an error, or a crash
labels: ''
---

**What happened, and what you expected instead**

**A script that reproduces it**

`simulate_breeding()` ships with the package and produces a pedigree with records,
so most reports can be written on simulated data. Please do not paste data you are
not free to share.

```r
library(BreedingR)

```

**The fit, if there is one**

`summary()` of the object, or at least `converged`, `theta` and the iteration count.

A fit reporting `converged = FALSE` is not necessarily a defect: a component resting
on its floor, a correlation at 1, or a covariance group the data does not identify
are all legitimate places to stop. Saying how many animals, how many groups and how
many records per animal the design has usually settles which of the two it is.

**Versions**

```r
br_version()
sessionInfo()
```
