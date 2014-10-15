### This script works as described but is NOT THE FINAL PRODUCT.
### I still need to get file names for original Excel files
### (see comments in BR function script).
### Add comment indicating that this is tab CalcLik in Excel files.
### Does it make more sense for the input observations to be structured as
### taxa in rows/sites in columns (current) or sites in rows/taxa in columns?
### Consider renaming main function "calcIEC", "estIEC", "estSiteIEC", or
### "siteIEC".
#
# This script defines the function "score_sites" which calculates Index of
# Ecological Condition (IEC) scores. The function "score_sites"
# ("iec_site_score" at the time) was originally used for Erin E. G. Giese's
# thesis (2012).  Giese's thesis used the script
# "IEC_Application_NoForLoopExp20120328.r" which this current script is based
# on.  The purpose of this script is to estimate IEC scores based on Biotic
# Response (BR) functions of species or other taxa for many sites at once.
# BR functions can be calculated using function "mk_brc" in the R script
# "BR_function.r".
#
# Input:
# The first argument ("sp") of the function "score_sites" takes a data frame
# with the following structure: (1) the column names must contain the species
# names (or other taxa as defined by the analyst), (2) all rows must
# contain the observations for each species at each site (one row per site;
# see example at the beginning of the function).  Observations for method
# "calcLSq" should only be probabilities, but observations for method
# "calcLik" (default) can be probabilities, abundance, or presence/absence.
# The second argument ("brc") to "score_sites" is a data frame produced by
# function "mk_brc".  Data frame "brc" contains the parameters describing the
# BR function for each species (or other taxonomic groups) to be used in
# estimating site IEC scores.
# The remaining arguments to "score_sites" are optional.
# "method" is used to set the criteria to use for estimating IEC.  "method"
# can be set to "calcLik" (default), or "calcLSq".  Method "calcLik" is used
# where present/absents data are available, and method "calcLSq" is used if
# probability data are available.
# "n_reps" defines the number of random starting values that will be tried by
# function "nlminb".  It defaults to 30, but can be set to any positive value.
#
# Note: As written, input data frames "sp" and "brc" must contain the same
# species.  Error checking will alert the analyst if there are differences
# between the two data frames.  The order of the species is not important as
# the function reorders them before starting the analysis.  Future versions
# of "score_sites" may be written so that only species in "brc" are used for
# scoring sites.
#
# Output:
# The function "score_sites" returns a data frame with sites in rows and the
# estimated IEC and log-likelihood for each site in columns.  It also writes
# this data frame to CSV file in the current working directory.  The naming
# format of the output files includes year, month, date, hour, minute, and
# second.
#
# Structure of function "score_sites":
# * Section "Error checking" ensures that input variables "sp" and "brc"
#   are both data frames, and that they contain the same species.  If any of
#   these criteria are not met, an error is raised and the function exits with
#   a message describing the issue.  This section also sorts both input data
#   frames based on species.
# * Section "Declarations" contains variables that will be used later in the
#   function.  In general, these are the only variables you are likely to
#   modify in this script.
# * Section "Main Function" contains function "f" which returns the likelihood
#   function which is minimized by "nlminb".  Note that Dr. Howe's MS Excel
#   based versions of this function maximize this function. Because "nlminb"
#   can only minimize functions, the function here has the sign reversed
#   during analysis, but the sign is corrected before adding it to the output
#   so the numbers generated are comparable to those generated by Solver in
#   Excel.
# * Section "for loop going over all sites" contains a for loop that processes
#   each taxa in turn.  Each taxa has function "nlminb" called on it the
#   number of times specified by input variable "n_reps".  "nlminb" tries to
#   minimize  the likelihood function returned by function "f".
# * Section "Save output" saves the resulting IEC scores estimated for each
#   site to CSV file.  Should you forget to set the output of "score_sites" to
#   a variable, the results can simply be loaded from the output csv file.
#
# Authors: Nicholas G. Walton and Robert W. Howe
#
# File used for E. Giese's thesis: "IEC_Application_NoForLoopExp20120328.r"
# File "IEC_Builder20120328.r" was based on March Excel files but most closely
# resembles April files (right?  based on BRC script comments)
#
# Created: Mar 2013(?)
# Last modified: 14 Oct 2014


score_sites <- function(sp, brc, method = "calcLik", n_reps = 30) {

  # The input "sp" is a data frame containing the observations for each
  # species or other taxa (abundance or presence/absence) at each site.
  # The column names of "sp" must be the species or other taxa.
  # All rows must contain the observations of each species at each site
  # (one column per site).  Row names are site names.

  # Example of "sp":
  # row.names  Sp(1)   Sp(2)   Sp(?)   Sp(n)
  # Site(1)       0       0       0       0
  # Site(2)       1       0       1       0
  # Site(?)       0       0       0       0
  # Site(n)       0       3       0       0

  # The input "brc" is a data frame output produced by function "mk_brc".
  # "n_reps" is the number of times that "nlminb" will be run with
  # different starting values.



  ## Error checking ----

  # Check that the input "sp" and "brc" are both data frames.
  # If not, print an error message and exit function "score_sites".
  if (!is.data.frame(sp) | !is.data.frame(brc)) {
    stop("The first and second input variables must be data frames.")
  }

  # Check that "sp" and "brc" contain the same taxa.
  # If not, print an error message and exit function "score_sites".
  if (ncol(sp) != nrow(brc)) {  # Check for same number of taxa
    stop("Species observation and BRC tables have diffents number of taxa.")
  } else if (!identical(as.character(names(sp)), as.character(brc[, 1]))) {
    stop(paste("Species observation and BRC tables do not contain",
               "the same taxa \nor are not ordered the same."))
  }

  # Check that "method" has been set correctly.
  if (!method %in% c("calcLSq", "calcLik")) {
    stop("Method must be set to 'calcLSq' or 'calcLik' (default).")
  }


  ## Declarations ----

  # Generate a data frame to hold the output IEC site score data;
  # an empty data frame with 3 columns - data type defined.
  iec_scores <- data.frame(character(0), rep(list(numeric(0)), 2),
                              stringsAsFactors = FALSE)
  names(iec_scores) <- c("Site", "IEC", "LogLik")  # Add column names.

  # Function to generate stratified random selection of starting IEC values.
  get_strat <- function(n) {
    # Generalized to allow for other numbers to be set for "n" ("n_reps").
    # "n" does not have to be a multiple of 10.
    # Note that it's ok to set "nReps" as low as 10,
    # and it can be set even lower if you're just running tests.

    # Lower limit for each stratum (upper will be llimit + 1).
    llimit <- 0:9

    # This is the number of random values that will be selected from
    # each stratum.
    pick <- n %/% length(llimit)

    # Select stratified random values.
    strat <- unlist(lapply(llimit,
                           function(x) runif(pick, min = x, max = x + 1)))

    # If n has been set to a value that is not a multiple of 10,
    # the remainder will be filled using random values in [1, 10].
    # Note that only this part will be used to assign starting IEC values
    # if n is set to less then length(llimit) (AKA 10).
    if (length(strat) < n) {
      remain <- n - length(strat)
      strat[(length(strat) + 1):n] <- runif(remain, min = 0, max = 10)
    }
    strat   # Return the vector of starting values for "nlminb".
  }

  # Set method/criteria for estimating IEC
  # Return function for calcLSq or calcLik based on "method"
  # Method must be set to "calcLSq" or "calcLik"
  if (method == "calcLSq") {
    criteria <- function(pc, observed) {
      # Least-squares method (AKA "calcLSq").
      # Note that currently the output from this will have the
      # Oposite sign of that from the Excel spread sheet.
      # This needs to be fixed, but will also require modifications to
      # "calcLik".
      numerator <- (observed - pc) ^ 2
      sum(numerator / pc)   # (Obs-Exp)^2 / Exp
    }
  } else {
    criteria <- function(pc, observed) {
      # Likelihood method (AKA "calcLik").
      # Calculate lack-of-fit expression ("lof") for each species.
      lof <- rep(NA, length(pc))  # set lof to length of pc.
      # If species is present at the current site, set to Log P(C).
      lof[observed > 0] <- log10(pc)[observed > 0]
      # If species is absent, set to Log (1-P(C)).
      # Per the Excel file, it's really Log(1.0001-P(C)) to avoid log of 0.
      lof[is.na(lof)] <- log10(1.0001 - pc)[is.na(lof)]

      # Return the negative sum of the lof.
      # I'm using the negative because the orignial function in Excel was
      # set to maximize, but "nlminb" can only minimize a function.
      # The result is the same.
      -sum(lof)
    }
  }


  ## Main Function ----
  f <- function(x) {
    # This function returns the function which
    # will be minimized with function "nlminb".
    # "x" is the current IEC score being tried by "nlminb".
    # Note that observed is called from outside the function, but there
    # doesn't seem to be a better way to do this with "nlminb".

    # Calculate P(C) for each species.
    # In the original version, P(C) was constrained to > 0 and < or = 1.
    # P(C) is the probability or value of each species at given IEC.
    # Input values "Mean", "SD", and "H" are part of data frame "brc".
    pc <- with(brc, {dnorm(x, mean = Mean, sd = SD) * H})
    pc[pc == 0] <- .001        # Set 0 probabilities to .001.

    # criteria is set in the Declarations section of this script.
    # It is set to either calcLik (default) or calcLSq.
    criteria(pc, observed)
  }


  ## for loop going over all sites ----
  # This for loop iteratively processes each site in the data frame "sp".

  for (site in 1:nrow(sp)) {
    # A list to contain the best solution found by "nlminb".
    best <- NULL

    # Extract the observed probabilities for the current site.
    observed <- sp[site, ]

    # Function "get_strat" is defined in the "Declarations" section.
    iec_start <- get_strat(n_reps)

    # Estimate IEC for current site using "n_reps" different starting values.
    for (i in 1:n_reps) {
      # "n_reps" is set in the call to function "siteScore" and
      # defaults to 30.

      # For help with function "nlminb" type "?nlminb" in R.
      # This is the optimization function.
      # The result is constrained to be between 0 and 10 inclusive.
      result <- nlminb(start = iec_start[i], objective = f,
                       lower = 0, upper = 10)

      # Check if the current random starting value resulted in an
      # improved estimate of IEC.
      if (i == 1) {
        # On the 1st iteration set "best" to "result".
        best <- result
      } else if (result$objective < best$objective) {
        # Otherwise check if newest "result" is better than the current
        # best result ("best").
        # If the new "result" is better than the current best,
        # update "best" to "result".
        best <- result
      }
    }

    # After fitting a best solution, add the solution to the
    # data frame "iec_scores".
    # best$objective is multiplied by -1 to put it on the same scale as
    # the original Excel based versions.
    iec_scores[site, ] <- list(row.names(sp)[site], best$par, -best$objective)
  }

  # Return "iec_scores"
  iec_scores
}