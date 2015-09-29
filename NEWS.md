# iec v0.1.1

* Package `iec` now passes `R CMD check`.

* Updated Gnass Giese, et. al. (2015) Ecosphere reference

* Added new function `bin` for binning of `ref_grad` (can specify number of bins or number of approx. number of records per bin)

* Relaxed error checking in `est_brc` to allow for `ref_grad` min > 0 and max < 10 to allow for binning

* `est_brc` now fails when NAs are detected in input data

* Added option to `est_iec` to score sites based on only observed species at that specific site

* Added arguments to `plot_brc` and `brc_pdf` to allow for sub titles to main and x-axis labels

* Added the following to output from `plot_brc`: Mean square error (RMSE); Pearson correlation between reference gradient and species response; Max `ylim` is plotted at no less than 1

* Removed non-linear R^2^ and Sensitivity from output of `plot_brc`

* Fixed `plot_iec_cor` to work with vector `env_grad` (was n x 1 data frame)

* Changed reference data `fish_grad` from single column data frame to vector

* Fixed typos in documentation

# iec v0.1

* First release
