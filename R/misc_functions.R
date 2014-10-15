# Miscellanious function to support working with iec functions.
#
# Authors: Nicholas G. Walton and Robert W. Howe
# Created 14 Oct 2014
# Last modified: 14 Oct 2014


## sensitivity function ----

# Example:
# sp <- data.frame(a = runif(10, -1), b = runif(10, -1), c = runif(10, -1))
# grad <- runif(10, 0, 10)
# res <- get_sens(sp, grad)

get_sens <- function(obs_df, gradient, mid = 5) {
  # obs_df is a data frame containing species observations with
  # rows as sites, and columns as taxa.
  # gradient is a vector or data frame of reference gradient values
  # at each site.
  # Returns a data frame with taxa as rows, and columns Sens (sensitivity)
  # N (numer of cites where the taxon was detected).

  sens <- function(obs_vect, grad) {
    pres <- obs_vect > 0
    n <- sum(pres)
    low <- sum(pres & grad < mid)
    high <- n - low
    s <- (high - low) / n
    x <- c(Sens = s, n = n)
    x
  }

  res_matrix <- t(apply(obs_df, 2, sens, gradient))
  data.frame(res_matrix)
}


## scale funciton ----

scale10 <- function(x, invert = FALSE) {
  # Scales a numeric vector from 0 to 10.
  # If 'invert = TRUE', the scale will be reversed.
  s <- 10 * ((x - min(x)) / (max(x) - min(x)))
  if (invert) s <- 10 - s
  s
}
