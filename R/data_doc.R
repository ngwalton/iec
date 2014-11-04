#' Fish observations
#'
#' Small dataset of fish observations.
#'
#' @format \code{fish_sp} is a data frame containing records of 4 fish species
#'   (columns) from 25 sites (rows).
#' @seealso \code{\link{fish_grad}}
'fish_sp'

#' Environmental gradient for fish sites
#'
#' Environmental gradient scores for sites in \code{\link{fish_sp}}. Larger
#' values indicate less desirable sites, so the scale will need to be inverted
#' for use with \code{\link{est_brc}}.
#'
#' @format \code{fish_grad} is a single column data frame containing
#'   environmental gradient scores from 25 sites (rows).  Site corresponds to
#'   sites (rows) in \code{\link{fish_sp}}.
#' @examples
#' data(fish_grad)
#' grad <- scale10(fish_grad, invert = TRUE)
#' # grad is now scaled to 0-10 with 10 being the least degraded.
#' summary(grad)
#' @seealso \code{\link{scale10}}, \code{\link{est_brc}}
'fish_grad'
