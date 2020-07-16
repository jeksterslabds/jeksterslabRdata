#' Fit
#'
#' @description Fit a function to each element
#' of a list of data.
#'
#' @details The first argument of `fitFUN` should be `data`.
#' The output of `fitFUN` should be a vector.
#'
#' Options for explicit parallelism are provided.
#' See `par` and suceeding arguments.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family model fit functions
#' @keywords model fit
#' @inheritParams univ
#' @param Xstar List.
#' A list of length `R` with a data set of length `n`
#' in each element of the list.
#' @param fitFUN Function.
#' Fit function to be applied to each element of `Xstar`.
#' @param ... Argument to pass to `fitFUN`.
#' @return Returns a list of parameter estimates.
#' @examples
#' n <- 4
#' R <- 1000
#' # normal distribution---------------------------------------------------------
#' mu <- 100
#' sigma2 <- 225
#' sigma <- sqrt(sigma2)
#' xstar <- univ(
#'   n = n,
#'   rFUN = rnorm,
#'   R = R,
#'   mean = mu,
#'   sd = sigma
#' )
#' str(xstar, list.len = 6)
#' # fit the var function--------------------------------------------------------
#' thetahatstar <- fit(
#'   Xstar = xstar,
#'   fitFUN = var,
#'   rbind = TRUE,
#'   par = FALSE
#' )
#' str(thetahatstar)
#' # The sampling distribution of sample variances
#' # follows a chi-square distribution with df = n - 1
#' # when the population is normally distributed
#' hist(thetahatstar)
#' @export
fit <- function(Xstar,
                fitFUN,
                ...,
                par = FALSE,
                ncores = NULL,
                mc = TRUE,
                lb = FALSE,
                cl_eval = FALSE,
                cl_export = FALSE,
                cl_expr,
                cl_vars,
                rbind = TRUE) {
  # execute loop-----------------------------------------------------------------
  par_lapply(
    X = Xstar,
    FUN = fitFUN,
    ...,
    rbind = rbind,
    par = par,
    ncores = ncores,
    mc = mc,
    lb = lb,
    cl_eval = cl_eval,
    cl_export = cl_export,
    cl_expr = cl_expr,
    cl_vars = cl_vars
  )
  #---------------------------------------------------------------------------------
}

# #' Fit Central Moments
# #'
# #' Fit central moments to each element
# #' of a list of data.
# #'
# #' @author Ivan Jacob Agaloos Pesigan
# #' @inheritParams fit
# #' @importFrom stats var
# #' @export
# fit_moments <- function(Xstar) {
#  if (is.vector(Xstar)) {
#    univariate <- TRUE
#    var_name <- names(Xstar)
#  }
#  if (is.data.frame(Xstar) | is.matrix(Xstar)) {
#    var_name <- colnames(Xstar)
#    if (ncol(Xstar) == 1 | nrow(Xstar) == 1) {
#      univariate <- FALSE
#      Xstar <- as.vector(Xstar[, 1])
#    }
#  }
#  if (univariate) {
#    mean_data <- mean(Xstar)
#    var_data <- var(Xstar)
#    sd_data <- sqrt(var_data)
#    # give names
#    out <- c(
#      mean_data,
#      var_data,
#      sd_data
#    )
#  } else {
#    # give names
#    means_data <- colMeans(Xstar)
#    as.vector(cov(Xstar))
#    as.vector(cor(Xstar))
#  }
# }
