# iec

This package provides a suite of functions for calculating the Index of Ecological Condition (IEC), a biotic indicator of ecological health first described by Howe et al. (2007a,b) and modified by us for Giese et al. (2014). Calculation of an IEC involves two steps 1) modeling responses of species to a reference gradient or quantitative environmental stressor (typically completed by prior research) and 2) calculating IEC values for new sites based occurrences (e.g., presence/absence, abundance, frequency) of multiple species or taxonomic groups at the site. The method applies an iterative maximum likelihood approach for calculating both species response functions and IEC values. Functions for calculating the biotic responses to environmental stressors (BR models) are useful as stand-alone applications of environmental gradient analysis.


## Installation

The recommended way to install this package is using `install_github` from the [devtools](https://github.com/hadley/devtools) package.  Your computer will need to be able to compile from source.  Install [Rtools](http://cran.r-project.org/bin/windows/Rtools/) on Windows, or Xcode on Mac.  If you're on Linux, you're probably ready to go.

After installing Rtools or Xcode, install `devtools` and test your installation with the following:

```R
install.packages("devtools")
devtools::has_devel()
```

If `has_devel()` returns true, your system is ready to install from GitHub.

```R
devtools::install_github("ngwalton/iec")
```

## Using package iec

The main functions in package `iec` are `est_brc` and `est_iec`.  The first is used to estimate Biotic Response Curves (BRC) for the taxa of interest, and the second is used to estimate Index of Ecological Condition (IEC) scores for sites using the BRCs.

Use `help(package="iec")` to view the help index of all functions available in package `iec`.

Two vignettes will be included with package `iec` (but aren't completed).  These will be accessible from the help index or with one of the following:

```R
vignette("minimal_iec", package = "iec")
vignette("extended_iec", package = "iec")
```

## References
Gnass Giese, E.E., R.W. Howe, A.T. Wolf, N.A. Miller, and N.G. Walton. 2014. Sensitivity of breeding birds to the "human footprint" in western Great Lakes forest landscapes. In review.
 
Howe, R.W., R. R. Regal, J.M. Hanowski, G.J. Niemi, N.P. Danz, and C.R. Smith.  2007a.  An index of biotic condition based on bird assemblages in Great Lakes coastal wetlands.  Journal of Great Lakes Research 33 (Special Issue 3): 93-105. 
 
Howe, R.W., R. R. Regal, G.J. Niemi, N.P. Danz, J.M. Hanowski. 2007b.  A probability-based indicator of ecological condition. Ecological Indicators 7:793-806.
