# Tests results from est_brc and est_iec.

# Author: Nicholas G. Walton
# Created: 22 Oct 2014
# Last updated: 25 Sept 2015


context("BRC and IEC")

# file <- system.file("R", "sysdata.rda", package = "iec")
# load(file)

data(list = c("fish_grad","fish_sp"))

br <- function() {
  grad10 <- scale10(fish_grad, TRUE)
  brc <- iec::est_brc(fish_sp, grad10)
  # need to update in light of adding Direction and Magnitude to est_brc output
  brc[, -c(1, 7, length(brc))] <- round(brc[, -c(1, 7, length(brc))], digits = 4)
  brc <- brc[, -c(7, 8)]  # temporarily ignore Direction and Magnitude
  brc
}

brc <- br()

iec <- function(brc, method){
  score <- est_iec(fish_sp, brc, method, n_reps = 10)
  score[, -1] <- round(score[, -1], digits = 4)
  score
}


test_that("est_brc returns the same thing each time", {
  expect_equal_to_reference(brc, "./refs/brc_res.rds")
})

test_that("est_iec 'pa' returns the same thing each time", {
  expect_equal_to_reference(iec(brc, "pa"), "./refs/iec_pa_res.rds")
})

test_that("est_iec 'quant' returns the same thing each time", {
  expect_equal_to_reference(iec(brc, "q"), "./refs/iec_quant_res.rds")
})
