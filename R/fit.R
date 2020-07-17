#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Fit a Function to Each Element of a List of Data
#'
#' @description Fits a function to each element of a list of data.
#'
#' @details The first argument of `fitFUN` should be `data`
#'   The output of `fitFUN` should be a vector.
#'
#'   Options for explicit parallelism are provided.
#'   See `par` and suceeding arguments.
#'
#' @family model fit functions
#' @keywords model fit
#' @inheritParams univ
#' @param Xstar List.
#'   A list of length `R` with a data set of length `n` in each element of the list.
#' @param fitFUN Function.
#'   Fit function to be applied to each element of `Xstar`.
#' @param ... Argument to pass to `fitFUN`.
#' @return Returns a list of parameter estimates.
#' @examples
#' xstar <- univ(n = 100, rFUN = rnorm, mean = 100, sd = sqrt(225), R = 100)
#' thetahatstar <- fit(Xstar = xstar, fitFUN = mean, rbind = TRUE)
#' mean(thetahatstar)
#' Xstar <- mvn(
#'   n = 100,
#'   Sigma = matrix(c(1, .5, .5, 1), nrow = 2),
#'   R = 100
#' )
#' foo <- function(X) {
#'   as.vector(cov(X))
#' }
#' thetahatstar <- fit(Xstar = Xstar, fitFUN = foo, rbind = TRUE)
#' thetahatstar <- thetahatstar[, 2]
#' mean(thetahatstar)
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
