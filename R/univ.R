#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Generate Univariate Data
#'
#' @description Generates an \eqn{n \times 1} univariate data vector
#'   or a list of \eqn{n \times 1} univariate data vectors of length `R`.
#'   The default data generating function
#'   is the normal distribution
#'   \deqn{
#'     X \sim \mathcal{N} \left( \mu, \sigma^2 \right) .
#'     %(\#eq:dist-X-norm)
#'   }
#'
#' @details The univariate distribution and parameters used
#'   in the data generating process can be specified using `rFUN` and `...`.
#'
#'   Options for explicit parallelism are provided when `R > 1`
#'   especially when `R` is large. See `par` and suceeding arguments.
#'
#' @return If `R = NULL` or `R = 1`,
#'   returns an \eqn{ n \times 1 } univariate data vector generated
#'   using `rFUN` and parameters passed to `...`.
#'   If `R` is an integer greater than 1, (e.g., `R = 10`)
#'   returns a list of length `R` of \eqn{ n \times 1 } univariate data vector
#'   generated using `rFUN` and parameters passed to `...`.
#'
#' @family univariate data functions
#' @keywords univariate
#' @inheritParams jeksterslabRpar::par_lapply
#' @importFrom stats rnorm
#' @importFrom jeksterslabRpar par_lapply
#' @param n Integer.
#'   Sample size.
#' @param rFUN Function.
#'   Data generating function to generate univariate data.
#' @param ... Arguments to pass to `rFUN`.
#' @param R Integer.
#'   Number of Monte Carlo replications. If `R` is not provided,
#'   the function produces a single random data set.
#'   If `R` is an integer greater than 1, (e.g., `R = 10`),
#'   the function produces multiple random data sets
#'   stored in each element of a list of length `R`.
#'   `par` and all succeeding arguments are only relevant when `R > 1`.
#' @examples
#' x <- univ(n = 100, rFUN = rnorm, mean = 100, sd = sqrt(225))
#' str(x)
#' xstar <- univ(n = 100, rFUN = rnorm, mean = 100, sd = sqrt(225), R = 100)
#' str(xstar, list.len = 6)
#' @export
univ <- function(n, # rFUN
                 rFUN = rnorm,
                 ...,
                 # iter
                 R = NULL,
                 # par_lapply
                 par = FALSE,
                 ncores = NULL,
                 mc = TRUE,
                 lb = FALSE,
                 cl_eval = FALSE,
                 cl_export = FALSE,
                 cl_expr,
                 cl_vars,
                 rbind = NULL) {
  # main function----------------------------------------------------------------
  foo <- function(iter,
                  n,
                  rFUN,
                  ...) {
    rFUN(
      n = n,
      ...
    )
  }
  #------------------------------------------------------------------------------
  # negotiate par and R----------------------------------------------------------
  if (is.null(R)) {
    single_rep <- TRUE
  } else {
    if (R == 1) {
      single_rep <- TRUE
    } else {
      single_rep <- FALSE
    }
  }
  #------------------------------------------------------------------------------
  if (single_rep) {
    # single rep-----------------------------------------------------------------
    return(
      foo(
        n = n,
        rFUN = rFUN,
        ...
      )
    )
    #----------------------------------------------------------------------------
  } else {
    # multiple reps--------------------------------------------------------------
    return(
      par_lapply(
        # iter
        X = 1:R,
        # foo
        FUN = foo,
        # pass to foo
        n = n,
        rFUN = rFUN,
        ...,
        # par_lapply
        par = par,
        ncores = ncores,
        mc = mc,
        lb = lb,
        cl_eval = cl_eval,
        cl_export = cl_export,
        cl_expr = cl_expr,
        cl_vars = cl_vars,
        rbind = rbind
      )
    )
    #----------------------------------------------------------------------------
  }
}
