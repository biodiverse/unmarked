# R package unmarked

<!-- badges: start -->
[![R build status](https://github.com/hmecology/unmarked/workflows/R-CMD-check/badge.svg)](https://github.com/hmecology/unmarked/actions)
[![CRAN status](https://www.r-pkg.org/badges/version/unmarked)](https://cran.r-project.org/package=unmarked)
<!-- badges: end -->

`unmarked` is an [R](https://www.r-project.org/) package for analyzing ecological data arising from several popular sampling techniques.  The sampling methods include point counts, occurrence sampling, distance sampling, removal, double observer, and many others. `unmarked` uses hierarchical models to incorporate covariates of the latent abundance (or occupancy) and imperfect detection processes.

## Installation

The latest stable version of unmarked can be downloaded from [CRAN](https://cran.r-project.org/package=unmarked):

```{r, eval=FALSE}
install.packages("unmarked")
```

The latest development version can be installed from Github:

```{r, eval=FALSE}
install.packages("remotes")
remotes::install_github("hmecology/unmarked")
```

## Support

Support is provided through the [unmarked Google group](http://groups.google.com/group/unmarked).
The package [website](https://hmecology.github.io/unmarked) has more information.
You can report bugs [here](https://github.com/hmecology/unmarked/issues), by posting to the Google group, or by emailing [the current maintainer](https://kenkellner.com).

## Get Started

See the following vignettes for an introduction to `unmarked` and some example analyses:

[Overview of `unmarked`](https://hmecology.github.io/unmarked/articles/unmarked.html)

[Dynamic occupancy models](https://hmecology.github.io/unmarked/articles/colext.html)

[Multispecies occupancy models](https://hmecology.github.io/unmarked/articles/occuMulti.html)

[Distance sampling](https://hmecology.github.io/unmarked/articles/distsamp.html)

[Contributing to `unmarked`](https://hmecology.github.io/unmarked/articles/contributing_to_unmarked.html)
