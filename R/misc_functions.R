# Miscellaneous function to support working with iec functions.
#
# Authors: Nicholas G. Walton and Robert W. Howe
# Created 14 Oct 2014
# Last modified: 21 Oct 2014


## sensitivity function ----

#' Sensitivity of taxa to a gradient.
#'
#' \code{get_sens} returns the sensitivity of each taxon in a community data
#' frame (\code{sp}) to a reference gradient (\code{gradient}).
#'
#' Sensitivity is defined here by first finding how many detections were made
#' above and below the midpoint (\code{mid}) of the gradient, taking their
#' difference, and finally dividing this by the number of detections. Taxa with
#' more detections below (\code{mid}) than above will have a negative
#' sensitivity and vice versa. Note that while function \code{\link{est_brc}}
#' requires that the reference gradient be scaled to 0-10, \code{get_sens} does
#' not.
#'
#' @param mid numeric scalar taken as the gradient's midpoint for sensitivity
#'   (default is 5).
#' @inheritParams est_brc
#' @return A data frame containing each taxon's sensitivity (\code{Sens}) to the
#'   gradient, and the number of sites the they were detected at (\code{n}).
#' @examples
#' x <- function() sample(0:1, 10, replace = TRUE)
#' sp <- data.frame(sp1 = x(), sp2 = x(), sp3 = x())
#' grad <- runif(10, 0, 10)
#' get_sens(sp, grad)
get_sens <- function(sp, gradient, mid = 5) {

  sens <- function(sp_vec, grad) {
    pres <- sp_vec > 0
    n <- sum(pres)
    low <- sum(pres & grad < mid)
    high <- n - low
    s <- (high - low) / n
    c(Sens = s, n = n)
  }

  res_matrix <- t(apply(sp, 2, sens, gradient))
  data.frame(res_matrix)
}


## scale function ----

#' Scale gradient to 0-10
#'
#' \code{scale10} scales an environmental gradient to the range 0-10.
#' Optionally, it can invert the scale.  This is essential for use with
#' \code{\link{est_brc}}.
#'
#' @param invert logical indicating if the gradient should be inverted
#'   (defaults to \code{FALSE}).
#' @inheritParams est_brc
#' @return A numeric vector with \code{gradient} scaled to 0-10.
#' @examples
#' grad <- runif(10)
#' scale10(grad, TRUE)
scale10 <- function(gradient, invert = FALSE) {
  s <- 10 * ((gradient - min(gradient)) / (max(gradient) - min(gradient)))
  if (invert) s <- 10 - s
  s
}
