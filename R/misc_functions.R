# Miscellaneous function to support working with iec functions.
#
# Authors: Nicholas G. Walton and Robert W. Howe
# Created 14 Oct 2014
# Last modified: 23 Sept 2015


## Nonlinear R^2 ----

#' Nonlinear R-squared
#'
#' \code{nlr2} calculates a nonlinear R-squared for the observed values.
#'
#' This function calculates the non-linear R-squared.  The equation used is:
#' \deqn{R^{2} = 1 - \frac{SS_{reg}}{SS_{tot}}}{R^2 = 1 - (SS_reg / SS_tot)}
#' (see p. 34-35
#' \href{http://www.mcb5068.wustl.edu/MCB/Lecturers/Baranski/Articles/RegressionBook.pdf}{here}).
#' Note that this nonlinear R-squared is not anything squared so it can result
#' a negative value. Negative values indicate poor model fit.
#'
#' @param observed a numeric vector of observations.
#' @param expected a numeric vector of expected values.
#' @return The nonlinear R-squared value (numeric scalar).
nlr2 <- function(observed, expected) {
  # This function returns the nonlinear R^2 for the solutions found by
  # nlminb(). The nonlinear R^2 is not really anything squared so it can
  # result a negative value. Negative values indicate very poor model fit.

  # Sum of squared deviances from the nonlinear regression line.
  ss_reg <- sum((observed - expected) ^ 2)

  # Sum of squared deviances from the mean of the observed values.
  ss_tot <- sum((observed - mean(observed)) ^ 2)

  # Return the nonlinear R^2.
  1 - (ss_reg / ss_tot)
}


## sensitivity function ----

#' Sensitivity of taxa to an environmental gradient.
#'
#' \code{get_sens} returns the sensitivity of each taxon in a community data
#' frame (\code{sp}) to an environmental gradient (\code{env_grad}).
#'
#' Sensitivity is defined here by first finding how many detections were made
#' above and below the midpoint (\code{mid}) of the gradient, taking their
#' difference, and finally dividing this by the number of detections. Taxa with
#' more detections below (\code{mid}) than above will have a negative
#' sensitivity and vice versa. Note that while function \code{\link{est_brc}}
#' requires that a reference gradient scaled to 0-10, \code{get_sens} does not
#' require this.
#'
#' @param mid numeric scalar taken as the gradient's midpoint for sensitivity
#'   (default is 5).
#' @inheritParams scale10
#' @inheritParams est_brc
#' @return A data frame containing each taxon's sensitivity (\code{Sens}) to the
#'   gradient, and the number of sites the they were detected at (\code{n}).
#' @examples
#' x <- function() sample(0:1, 10, replace = TRUE)
#' sp <- data.frame(sp1 = x(), sp2 = x(), sp3 = x())
#' grad <- runif(10, 0, 10)
#' get_sens(sp, grad)
get_sens <- function(sp, env_grad, mid = 5) {

  sens <- function(sp_vec, grad) {
    pres <- sp_vec > 0
    n <- sum(pres)
    low <- sum(pres & grad < mid)
    high <- n - low
    s <- (high - low) / n
    c(Sens = s, n = n)
  }

  res_matrix <- t(apply(sp, 2, sens, env_grad))
  data.frame(res_matrix)
}


## scale function ----

#' Scale environmental gradient to 0-10
#'
#' \code{scale10} scales an environmental gradient to the range 0-10.
#' Optionally, it can invert the scale.  This is essential for use with
#' \code{\link{est_brc}}.
#'
#' \code{scale10} is used to scale an environmental gradient for use in
#' generating Biotic Response Curves with \code{\link{est_brc}}.  The reference
#' gradient returned by \code{scale10} has minimum 0 and maximum 10. If the
#' input environmental gradient has larger values for less desirable sites, use
#' \code{invert = TRUE} to invert the gradient so that 10 represents the most
#' desirable condition.
#'
#' @param env_grad numeric vector of environmental gradient scores, one for each
#'   site.
#' @param invert logical indicating if \code{env_grad} should be inverted
#'   (default is \code{FALSE}).
#' @return A numeric vector with \code{env_grad} scaled to 0-10.
#' @examples
#' grad <- runif(10)
#' scale10(grad, TRUE)
scale10 <- function(env_grad, invert = FALSE) {
  ref_grad <- 10 * ((env_grad - min(env_grad)) / (max(env_grad) - min(env_grad)))
  if (invert) ref_grad <- 10 - ref_grad
  ref_grad
}


## Bin function ----

#' Bin observations
#'
#' \code{bin} combines species observations based on input environmental
#'  gradient. Number of bins or number of records can be specified. Resulting
#'  bins can also be examined using the \code{summary} option.
#'
#' \code{bin} is used to "bin" (i.e., combine by taking the mean) species
#' species observations for taxa used to create Biotic Response functions
#' (\code{\link{est_brc}}). Unlike other functions in package \code{iec},
#' function \code{bin} takes a data frame that contains the environmental
#' gradient as the first column.
#'
#' @param sp data frame containing species observations as columns with an
#'   environmental gradient as the first column (needs work).
#' @param n scalar value indicating the number of bins or number of
#'   records per bin (see \code{n_bins}).
#' @param n_bins logical indicating if \code{n} is the number of bins
#'   (\code{TURE}; default) or records per bin.
#' @param summary logical indicating if \code{bin} should return the binned
#'   records (\code{FALSE}; default), or a summary bins
#' @return If \code{summary = FALSE} (default), \code{bin} returns data frame
#'   containing the binned species records and environmental gradient.  Both of
#'   these are means. If presence/absence data are used (codes as 0/1), the
#'   mean is also a probability. If \code{summary = TRUE}, \code{bin} returns
#'   a table of the number of records per bin.
#' @examples
#' x <- data.frame(gradient = round(abs(rnorm(107)), 2),
#'                 sp1 = rpois(107, 0.7),
#'                 sp2 = rpois(107, 0.5))
#'
#' bin(x, n = 7, summary = TRUE)
#' x_bin <- bin(x, n = 7)
bin <- function(sp, n, n_bins = TRUE, summary = FALSE) {
  # assumes the first column in sp contains the environmental gradient
  # note that the returned species values will be probabilities if the input
  # is pres/abs

  # setting n_bins = TRUE allows setting the number of bins used
  # setting n_bins = FALSE allows setting the number of records per bin
  # order records by gradient scores
  sp <- sp[order(sp[, 1]), ]

  # if n_bins was set to FALSE, adjust n
  if (!n_bins) {
    n <- floor(nrow(sp) / n)
  }

  bin_levels <- cut(1:nrow(sp), n, labels = FALSE, include.lowest = TRUE)

  if (summary) {
    # just return the number of records per bin
    table(bin_levels)
  } else {
    out <- t(sapply(unique(bin_levels),
                    function(x) colMeans(sp[bin_levels == x, ])))
    data.frame(out)  # t returns a matrix
  }
}
