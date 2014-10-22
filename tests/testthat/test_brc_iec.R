# Tests results from est_brc and est_iec.

# Author: Nicholas G. Walton
# Created: 22 Oct 2014
# Last updated: 22 Oct 2014


context("BRC and IEC")

file <- system.file("R", "sysdata.rda", package = "iec")
load(file)

br <- function() {
  grad10 <- scale10(fish_grad[, 1], TRUE)
  brc <- iec::est_brc(fish_sp, grad10)
  brc[, -c(1, length(brc))] <- round(brc[, -c(1, length(brc))], digits = 4)
  brc
}

brc <- br()

iec <- function(brc){
  score <- est_iec(fish_sp, brc, n_reps = 10)
  score[, -1] <- round(score[, -1], digits = 4)
  score
}


test_that("est_brc returns the same thing each time", {
  expect_equal_to_reference(brc, "./refs/brc_res.rds")
})

test_that("est_iec returns the same thing each time", {
  expect_equal_to_reference(iec(brc), "./refs/iec_res.rds")
})
