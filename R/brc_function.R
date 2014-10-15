### This script works as described but is NOT THE FINAL PRODUCT.
### All that remains:
### I still need to get file names for original Excel files (see comments)
#
# This script defines the function "mk_brc" which calculates Biotic Response
# (BR) functions, AKA Biotic Response Curves (BRC). The function "mk_brc"
# ("iec_builder" a the time) was originally used for Erin E. G. Giese's thesis
# (2012).  Giese's thesis used the script "IEC_Builder20120328.r" which this
# current script is based on.  The purpose of this script is to calculate
# Biotic Response (BR) functions for many species or other taxa using a
# single function.
#
# Input:
# The first argument ("sp") of the function "mk_brc" excepts a
# data frame with sites as rows and taxa as columns (community data frame).
# The second argument ("gradient") is a numeric vector containing the gradient
# score for each site.  This vector must has the same order as "sp" and must
# be scaled to 0 - 10 with 10 indicated the desirable condition.
#
# Output:
# The function "mk_brc" returns a data frame with taxa in rows (1st column) and
# the calculated values for each taxa in columns.
#
# Structure of function "mk_brc":
# * Section "Declarations" contains variables that will be used later in the
#   function.  These include the constraints placed on the optimization
#   function nlminb().  This is where you could modify the constraints if
#   desired.
#   In general, these are the only variables you are likely to modify in this
#   script.
# * Section "Functions" contains two functions embedded in "mk_brc".  Function
#   "f" returns the lack-of-fit criterion which is minimized by nlminb().
#   Function "nl_r2" returns the none-linear R-square value used for evaluating
#   the BRC after analysis.
# * Section "for loop going over all species" contains a for loop that
#   processes each taxa in turn. nlminb() tries to minimize the lack-of-fit
#   score returned by function "f".
#
# Authors: Nicholas G. Walton and Robert W. Howe
# Created: Mar 2012
# Last updated: 14 Oct 2014
#
# Original file (used for E. Giese's thesis): "IEC_Builder20120328.r"
#   File "IEC_Builder20120328.r" was based on March Excel files but most
#   closely resembles April files.
# This script differs from E. Giese's thesis in that the observations in "sp"
# are not constrained to be probabilities.


mk_brc <- function(sp, gradient) {
  # The input "sp" is a data frame containing the observations of
  # each species or other taxa.  The column names of "sp" must
  # be site names.  All columns must contain
  # the probability of each species at each site (rows).

  # Example of "sp":
  # row.names   species(1)    species(2)    species(...)   species(n)
  # "111"           0.104         0.291             ...        0.195
  # "156"           0.266         0.391             ...        0.040
  # "116"           0.023         0.663             ...        0.199
  #  ...              ...           ...             ...          ...
  # "120"           0.693         0.737             ...        0.532

  # Input "gradient" in a numeric vector containing the environmental gradient
  # for each site. "gradient" must be scaled from 0 - 10 with 10 being the
  # desirable condition. "sp" and "gradient" must have the same order by site.


  ## Error checking ----

  # Check that the input "sp" is a data frame.
  # If not, print an error message and exit function "mk_brc".
  if (!is.data.frame(sp)) {
    stop("The first input variable must be a data frame.")
  }

  # Check that "gradient" is scaled correctly.
  if (min(gradient) != 0 | max(gradient) != 10) {
    stop("gradient is not scaled from 0 to 10.")
  }

  # Check that "gradient" and "sp" have the same number of rows.
  if (length(gradient) != nrow(sp)) {
    stop("sp and gradient do not have the same numer of rows.")
  }


  ## Declarations ----

  # Set the number of times the nlminb() function will be run on each species.
  # n_reps <- 100

  # Constraints on optimization - nlminb()
  mu_min <- -10  # Mean lower bound
  mu_max <-  20  # Mean upper bound
  sd_min <-   2  # SD lower bound
                 #   Results appear to be better with higher SD (was 1e-10)
  sd_max <-  10  # SD upper bound
  ht_min <-   0  # Height scaler lower bound
  ht_max <- Inf  # Height scaler upper bound

  # Set variable "gradient" to the first column in the input data frame which
  # contains the condition gradient.
  # gradient <- sp[, 1]

  # Strip the gradient from "sp" - simplifies later for-loops.
  # sp <- sp[, -1]

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
    # Note that gradient and observed are called from outside the function,
    # but there doesn't seem to be a better way to do this with "nlminb".

    # Long version of LOF score:
    # [(observed-expected)^2]/expected

    mu_f <- x[1]  # Mean
    sd_f <- x[2]  # Standard deviation
    ht_f <- x[3]  # A scaling factor for height of normal curve

    # Calculate the expected values (Exp) based on the current curve.
    expected <- dnorm(gradient, mu_f, sd_f) * ht_f

    # Calculate the squared deviation for each observation - (Obs-Exp)^2.
    obs_ex2 <- (observed - expected) ^ 2

    # Calculate and return the LOF statistic.
    # Ensure that expected is 0.001 or greater.
    expected[expected < 0.001] <- 0.001  # faster than using pmax

    sum(obs_ex2 / expected) # Return the Lack of Fit score (O-E)2/(E).
  }

  nl_r2 <- function(observed, mu, sd, ht) {
    # This function returns the nonlinear R^2 for the solutions found by
    # nlminb(). The nonlinear R^2 is not really anything squared so it can
    # result a negative value. Negative values indicate very poor model fit.
    # Definition for nonlinear R^2 used:
    # http://www.mcb5068.wustl.edu/MCB/Lecturers/Baranski/Articles/RegressionBook.pdf (p.34-35)

    # Expected values for each site.
    expected <- dnorm(gradient, mu, sd) * ht

    # Sum of squared deviances from the nonlinear regression line.
    ss_reg <- sum((observed - expected) ^ 2)

    # Sum of squared deviances from the mean of the observed values.
    ss_tot <- sum((observed - mean(observed)) ^ 2)

    # Return the nonlinear R^2.
    1 - (ss_reg / ss_tot)
  }

  # Note that this function should be removed when we switch to making this a
  # function.
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


  ## Calculate sensistivity ----

  sens <- get_sens(sp, gradient)


  ## for loop going over all species ----

  # This for loop iteratively processes each species in the data frame
  # "sp".
  for (species in 1:ncol(sp)) {

    # A list to contain the best solution found by "nlminb()".
    best <- NULL

    # Extract the observed probabilities for the current species.
    observed <- sp[, species]

    # generate random starting values
    br_start <- get_strat(mu_min, mu_max, sd_min, sd_max)

    for (i in 1:nrow(br_start)) {

      # "n_reps" and the following bounds are set in the "Declarations"
      # section of this function.

      # The following values are random staring values for the mean,
      # standard deviation, and height.  These will be used by function
      # "nlminb" as a starting point for its algorithm.
      # mu_rand <- runif(1, min = mu_min, max = mu_max)
      # sd_rand <- runif(1, min = sd_min, max = sd_max)
      # ht_rand <- runif(1, min = ht_min, max = 50)



      # Run "nlminb"
      result <- nlminb(start = br_start[i, ], objective = f,
                       lower = c(mu_min, sd_min, ht_min),
                       upper = c(mu_max, sd_max, ht_max))

      # On the 1st iteration set "best" to "result".
      # After the 1st iteration, check if newest "result" is better than
      # the current best result ("best").
      if (i == 1) {
        best <- result

        # Add R^2.
        best$R2 <- nl_r2(observed, result$par[1], result$par[2], result$par[3])
      } else if (result$objective < best$objective) {
        # If the new "result" is better than the current best, update
        # "best" to "result".
        best <- result
        best$R2 <- nl_r2(observed, result$par[1], result$par[2], result$par[3])
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

  brc_pars <- cbind(brc_pars, sens)


  # return brc_pars
  brc_pars
}
