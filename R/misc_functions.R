# Miscellanious function to support working with iec functions.
#
# Authors: Nicholas G. Walton and Robert W. Howe
# Created 14 Oct 2014
# Last modified: 15 Oct 2014


## sensitivity function ----

#' Sensitivity of taxa to a gradient.
#'
#' \code{get_sens} returns the sensisitivity of each taxon in a community data
#'  frame (\code{sp}) to a reference gradient (\code{gradient}). Taxa with more
#'  detections below the midpoint (\code{mid}) will have a negative sensitivity
#'  and vice versa.
#'
#' @param sp A community data frame (sites as rows, taxa as columns).
#' @param gradient A numeric vector of reference gradient scores (0-10),
#'   one for each site.
#' @param mid A numeric scalar taken as the gradient's midpoint for sensitivity.
#' @return A data frame containing each taxon's sensitivity (\code{Sens}) to the
#'   gradient, and the number of sites the they were detected at (\code{n}).
#' @examples
#' x <- function() sample(0:1, 10, replace = TRUE)
#' sp <- data.frame(sp1 = x(), sp2 = x(), sp3 = x())
#' grad <- runif(10, 0, 10)
#' get_sens(sp, grad)
get_sens <- function(sp, gradient, mid = 5) {
  # sp is a data frame containing species observations with
  # rows as sites, and columns as taxa.
  # gradient is a vector of reference gradient values
  # at each site.
  # Returns a data frame with taxa as rows, and columns Sens (sensitivity)
  # N (numer of cites where the taxon was detected).

  sens <- function(sp_vec, grad) {
    pres <- sp_vec > 0
    n <- sum(pres)
    low <- sum(pres & grad < mid)
    high <- n - low
    s <- (high - low) / n
    x <- c(Sens = s, n = n)
    x
  }

  res_matrix <- t(apply(sp, 2, sens, gradient))
  data.frame(res_matrix)
}


## scale funciton ----

#' Scale gradient to 0-10
#'
#' \code{scale10} scales an environmental gradient to the range 0-10.
#' Optionally, it can invert the scale.  This is essential for use with
#' \link{get_brc}.
#'
#' @param gradient A numeric vector containing gradient values.
#' @param invert Logical scalar indicatingin if the gradient should be inverted
#'   (defaults to \code{FALSE}).
#' @return A numeric vector with \code{gradient} scaled to 0-10.
#' @examples
#' grad <- runif(10)
#' scale10(grad, TRUE)
scale10 <- function(gradient, invert = FALSE) {
  # Scales a numeric vector from 0 to 10.
  # If 'invert = TRUE', the scale will be reversed.
  s <- 10 * ((gradient - min(gradient)) / (max(gradient) - min(gradient)))
  if (invert) s <- 10 - s
  s
}
