# Script to compile data for package iec.

# Created: 22 Oct 2014
# Last modified: 25 Sept 2015
# Author: Nicholas G. Walton


## export a file for use in testing ----

fish0 <- read.csv("./data-raw/fish.csv", row.names = 1, comment.char = "#")

fish_sp <- fish0[1:25, c(1, 4, 5, 9, 10)]

fish_grad <- fish_sp[, 1]

fish_sp <- fish_sp[, -1]

devtools::use_data(fish_sp, fish_grad, overwrite = TRUE)
