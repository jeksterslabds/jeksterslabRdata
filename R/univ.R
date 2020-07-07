#' Generate Univariate Data
#'
#' Generates an \eqn{n \times 1}
#' univariate data vector.
#'
#' The univariate distribution and parameters used in the
#' data generating process can be specified using `rFUN` and `...`.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family univariate data functions
#' @keywords univariate
#' @inheritParams jeksterslabRutils::util_lapply
#' @param n Integer.
#'   Sample size.
#' @param rFUN Function.
#'   Data generating function to generate univariate data.
#' @param R Integer.
#'   Number of Monte Carlo replications.
#'   By default,
#'   `R = NULL` produces a single random data set.
#'   If `R` is defined
#'   (e.g., `R = 10`),
#'   produces a list of length `R`
#'   of multiple random data sets.
#' @param ... Arguments to pass to `rFUN`.
#' @return
#' If `R = NULL`,
#' returns an \eqn{n \times 1} univariate data vector
#' generated using `rFUN` and parameters passed to `...`.
#' If `R` is an integer,
#' returns a list of length `R` of
#' \eqn{n \times 1} univariate data vector
#' generated using `rFUN` and parameters passed to `...`.
#' @examples
#' n <- 5
#' R <- 5
#' #################################
#' # normal distribution
#' #################################
#' mu <- 100
#' sigma2 <- 225
#' sigma <- sqrt(sigma2)
#' x <- univ(
#'   n = n,
#'   rFUN = rnorm,
#'   mean = mu,
#'   sd = sigma
#' )
#' str(x)
#' x <- univ(
#'   n = n,
#'   rFUN = rnorm,
#'   R = R,
#'   mean = mu,
#'   sd = sigma
#' )
#' str(x)
#' #################################
#' # binomial distribution
#' #################################
#' n_trials <- 1
#' p <- 0.50
#' x <- univ(
#'   n = n,
#'   rFUN = rbinom,
#'   size = n_trials,
#'   prob = p
#' )
#' str(x)
#' x <- univ(
#'   n = n,
#'   rFUN = rbinom,
#'   R = R,
#'   size = n_trials,
#'   prob = p
#' )
#' str(x)
#' @importFrom stats rnorm
#' @importFrom jeksterslabRutils util_lapply
#' @export
univ <- function(n,
                 rFUN = rnorm,
                 R = NULL,
                 par = FALSE,
                 ncores = NULL,
                 ...) {
  foo <- function(iter,
                  rFUN,
                  n,
                  ...) {
    rFUN(
      n = n,
      ...
    )
  }
  if (is.null(R)) {
    out <- foo(
      rFUN = rFUN,
      n = n,
      ...
    )
  } else {
    args <- list(
      iter = 1:R,
      rFUN = rFUN,
      n = n,
      ...
    )
    out <- util_lapply(
      FUN = foo,
      args = args,
      par = par,
      ncores = ncores
    )
  }
  out
}
