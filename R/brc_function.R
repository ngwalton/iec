# This script defines the function "est_brc" which calculates Biotic Response
# functions, AKA Biotic Response Curves (BRC). The function "est_brc"
# ("iec_builder" at the time) was originally used for Erin E. G. Giese's thesis
# (2012).  Giese's thesis used the script "IEC_Builder20120328.r" which this
# current script is based on.  The purpose of this script is to calculate
# Biotic Response (BR) functions for many species or other taxa using a
# single function.
#
# Original file (used for E. Giese's thesis): "IEC_Builder20120328.r"
#   File "IEC_Builder20120328.r" was based on March Excel files but most
#   closely resembles April files.
# "est_brc" script differs from E. Giese's thesis in that observations in "sp"
# are not constrained to be probabilities.
#
# Structure of function "est_brc":
# * "Error checking" checks that "sp" is a data frame, that "ref_grad" is
#   scaled to the range 0-10, and that "sp" and "ref_grad" have the same
#   number records.
# * "Declarations" contains variables that will be used later in the
#   "est_brc".  These include the constraints placed on the optimization
#   function "nlminb".  This is where you could modify the constraints if
#   desired.
#   In general, these are the only variables you are likely to modify in this
#   function.
# * "Functions" contains function "f" which returns the lack-of-fit criterion
#   which is minimized by "nlminb".
# * "for loop over each species" contains a for loop that
#   processes each taxa in turn. "nlminb" tries to minimize the lack-of-fit
#   score returned by function "f".
#
# Authors: Nicholas G. Walton and Robert W. Howe
# Created: Mar 2012
# Last updated: 22 Oct 2014

#' Make biotic response curves (BRC).
#'
#' \code{est_brc} generates biotic response curves (BRC) for each species or
#' taxon in \code{sp} in relation to environmental reference gradient
#' \code{ref_grad}.
#'
#' The biotic response curves (BRCs) or functions returned by \code{est_brc} are
#' normal curves fit to the observations in \code{sp} using a lack-of-fit (LOF)
#' criteria. BRCs consist of a normal curve defined by mean (mu) and standard
#' deviation (sigma), which is multiplied (scaled) a height factor (H).  The
#' reference gradient (\code{ref_grad}) input to \code{est_brc} must be a be a
#' numeric vector scaled to 0-10 where 10 represents the least impacted site.
#' Use \code{\link{scale10}} to scale the reference gradient if needed. Note
#' that \code{ref_grad} must have the same order by site as \code{sp}. The
#' results of \code{est_brc} are used to give sites an Index of Ecological
#' Condition score using function \code{\link{est_iec}}.
#'
#' @param sp community data frame (sites as rows, taxa as columns, observations
#'   as values).
#' @param ref_grad numeric vector of reference gradient scores (0-10), one for
#'   each site.
#' @return A data frame defining BRCs for each species or other taxa in
#'   \code{sp}.
#' @seealso \code{\link{scale10}}
est_brc <- function(sp, ref_grad) {
  # The input "sp" is a data frame containing the observations of
  # each species or other taxa.  The row names of "sp" must
  # be site names.  All columns must contain
  # the observations of each species at each site (rows).

  # Example of "sp":
  # row.names   species(1)    species(2)    species(...)   species(n)
  # "111"           0.104         0.291             ...        0.195
  # "156"           0.266         0.391             ...        0.040
  # "116"           0.023         0.663             ...        0.199
  #  ...              ...           ...             ...          ...
  # "120"           0.693         0.737             ...        0.532

  # Input "ref_grad" in a numeric vector containing the environmental reference
  # gradient for each site. "ref_grad" must be scaled from 0 - 10 with 10 being
  # the desirable condition. "sp" and "ref_grad" must have the same order by
  # site.


  ## Error checking ----

  # Check that the input "sp" is a data frame.
  if (!is.data.frame(sp)) {
    stop("sp must be a data frame.")
  }

  # Check that input "ref_grad" is a vector.
  if (!is.vector(ref_grad)) {
    stop("ref_grad must be a vector.")
  }

  # Check that "ref_grad" is scaled correctly.
  if (min(ref_grad) != 0 | max(ref_grad) != 10) {
    stop("ref_grad is not scaled from 0 to 10.")
  }

  # Check that "ref_grad" and "sp" have the same number of records.
  if (length(ref_grad) != nrow(sp)) {
    stop("sp and ref_grad do not have the same numer of records.")
  }


  ## Declarations ----

  # Constraints on optimization - nlminb()
  mu_min <- -10  # Mean lower bound
  mu_max <-  20  # Mean upper bound
  sd_min <-   2  # SD lower bound
                 #   Results appear to be better with higher SD (was 1e-10)
  sd_max <-  10  # SD upper bound
  ht_min <-   0  # Height scaler lower bound
  ht_max <- Inf  # Height scaler upper bound

  # Generate a data frame to hold the output parameters of BR function.
  # This is an empty data frame with 6 columns - data types defined.
  brc_pars <- data.frame(character(0), rep(list(numeric(0)), 5),
                         stringsAsFactors = FALSE)
  # Add column names.
  # Consider changing to: c("taxon", "lof", "r2", "mu", "sigma", "ht")
  names(brc_pars) <- c("Taxon", "LOF", "R2", "Mean", "SD", "H")

  # Expand rows in brc_pars to same length as ncol in sp
  brc_pars[ncol(sp), ]  <- NA


  ## Functions ----

  f <- function(x) {
    # This function returns the Lack of Fit (LOF) score which will be
    # minimized with function "nlminb".
    # x is a numeric vector containing c(Mean, SD, H).
    # Note that ref_grad and observed are called from outside the function.
    # Consider passing from "nlminb", but may slow down processing.

    # LOF equation:
    # [(observed-expected)^2]/expected

    mu_f <- x[1]  # Mean
    sd_f <- x[2]  # Standard deviation
    ht_f <- x[3]  # A scaling factor for height of normal curve

    # Calculate the expected values (Exp) based on the current curve.
    expected <- dnorm(ref_grad, mu_f, sd_f) * ht_f

    # Calculate and return the LOF statistic.
    # Ensure that expected is 0.001 or greater.
    expected[expected < 0.001] <- 0.001  # Faster than using pmax.

    # Calculate the squared deviation for each observation.
    obs_ex2 <- (observed - expected) ^ 2

    sum(obs_ex2 / expected) # Return the Lack of Fit score.
  }


  get_strat <- function(mu_min, mu_max, sd_min, sd_max) {
    # returns a 48 x 3 data frame of stratified starting values for f

    # define a list to hold output
    start_val <- list()

    # generate breaks for mean, SD, and ht
    mu_br <- seq(from = mu_min, to = mu_max, length.out = 5) # mean breaks
    sd_br <- seq(from = sd_min, to = sd_max, length.out = 4) # sd breaks
    ht_br <- seq(from = 0, to = 50, length.out = 5)          # h breaks

    vals <- function(breaks) {
      # returns random numbers stratified by the bins specified in breaks.
      num_bins <- length(breaks) - 1
      runif(num_bins, min = breaks[1:num_bins], max = breaks[-1])
    }

    get_vals <- function(n, breaks) {
      unlist(lapply(1:n, function(x) vals(breaks)))
    }

    start_val$mu <- get_vals(12, mu_br)
    start_val$sd <- get_vals(16, sd_br)
    start_val$ht <- get_vals(12, ht_br)

    # order ht
    start_val$ht <- start_val$ht[order(rep(1:(length(ht_br) - 1), 12))]

    # return a data frame of stratified starting values for nlminb
    data.frame(start_val)
  }


  ## for loop over each species ----

  # This for loop iteratively processes each species in the data frame
  # "sp".
  for (species in 1:ncol(sp)) {

    # A list to contain the best solution found by "nlminb()".
    best <- NULL

    # Extract the observed probabilities for the current species.
    observed <- sp[, species]

    # generate random starting values
    br_start <- get_strat(mu_min, mu_max, sd_min, sd_max)

    # function to add R^2
    add_r2 <- function() {
      expect <- dnorm(ref_grad, result$par[1], result$par[2]) * result$par[3]
      iec::nlr2(observed, expect)
    }

    for (i in 1:nrow(br_start)) {

      # Run "nlminb"
      result <- nlminb(start = br_start[i, ], objective = f,
                       lower = c(mu_min, sd_min, ht_min),
                       upper = c(mu_max, sd_max, ht_max))

      # On the 1st iteration set "best" to "result".
      # After the 1st iteration, check if newest "result" is better than
      # the current best result ("best").
      if (i == 1) {
        best <- result

        # Add nonlinear R^2.
        best$R2 <- add_r2()
      } else if (result$objective < best$objective) {
        # If the new "result" is better than the current best, update
        # "best" to "result".
        best <- result
        best$R2 <- add_r2()
      }
    }
    # After fitting a best solution, add the solution to the data frame
    # "brc_pars".
    brc_pars[species, ] <- list(names(sp)[species],
                                      best$objective, best$R2,
                                      best$par[1], best$par[2],
                                      best$par[3])
  }


  ## Add sensitivity ----

  sens <- iec::get_sens(sp, ref_grad)

  brc_pars <- cbind(brc_pars, sens)


  # return brc_pars
  brc_pars
}
