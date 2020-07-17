#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Generate Multivariate Normal Data
#'   \eqn{
#'     \mathbf{X} \sim \mathcal{N}_{k}
#'     \left( \boldsymbol{\mu}, \boldsymbol{\Sigma} \right)
#'     %(\#eq:dist-X-mvn)
#'   }
#'
#' @description Generates an \eqn{n \times k} multivariate data matrix
#'   or a list of \eqn{n \times k} multivariate data matrices of length `R`
#'   from the multivariate normal distribution
#'   \deqn{
#'     \mathbf{X} \sim \mathcal{N}_{k}
#'     \left( \boldsymbol{\mu}, \boldsymbol{\Sigma} \right) .
#'     %(\#eq:dist-X-mvn)
#'   }
#'   This function is a wrapper around [`MASS::mvrnorm()`].
#'
#' @details The multivariate normal distribution has two parameters, namely
#'   the \eqn{k \times 1} mean vector `mu` \eqn{\left( \boldsymbol{\mu} \right)}
#'   and the \eqn{k \times k} variance-covariance matrix `Sigma`
#'   \eqn{\left( \boldsymbol{\Sigma} \right)}.
#'   If `mu` is not provided,
#'   it is set to a vector of zeroes with the appropriate length.
#'
#'   Options for explicit parallelism are provided when `R > 1`
#'   especially when `R` is large. See `par` and suceeding arguments.
#'
#' @return If `R = NULL` or `R = 1`, returns an \eqn{n \times k}
#'   multivariate normal data matrix or data frame .
#'   If `R` is an integer greater than 1, (e.g., `R = 10`)
#'   returns a list of length `R` of \eqn{n \times k}
#'   multivariate normal data matrix or data frame.
#'
#' @references Venables, W. N., Ripley, B. D., & Venables, W. N. (2002).
#'   *Modern applied statistics with S*.  New York, N.Y: Springer.
#'
#'   [Wikipedia: Multivariate normal distribution](https://en.wikipedia.org/wiki/Multivariate_normal_distribution)
#'
#' @seealso [`jeksterslabRdist::mvnpdf()`], [`jeksterslabRdist::mvnll()`],
#'   and [`jeksterslabRdist::mvn2ll()`],
#'   for more information on the multivariate normal distribution.
#'
#' @family multivariate data functions
#' @keywords multivariate, normal
#' @inheritParams jeksterslabRdist::mvnpdf
#' @inheritParams univ
#' @inherit jeksterslabRdist::mvnpdf details references
#' @importFrom MASS mvrnorm
#' @param tol Numeric.
#'   Tolerance (relative to largest variance)
#'   for numerical lack of positive-definiteness in `Sigma`.
#' @param empirical Logical.
#'   If `TRUE`, `mu` and `Sigma` specify the empirical,
#'   not population. mean and covariance matrix.
#' @param df Logical.
#'   If `TRUE`, the function returns a data frame.
#'   If `FALSE`, the function returns a matrix.
#' @param varnames Character string.
#'   Optional column names with the same length as `mu`.
#' @examples
#' mu <- c(100, 100, 100)
#' Sigma <- matrix(
#'   data = c(225, 112.50, 56.25, 112.5, 225, 112.5, 56.25, 112.50, 225),
#'   ncol = 3
#' )
#' X <- mvn(n = 100, mu = mu, Sigma = Sigma)
#' Xstar <- mvn(n = 100, mu = mu, Sigma = Sigma, R = 100)
#' str(Xstar, list.len = 6)
#' @export
mvn <- function(n, # mvrnorm
                mu = NULL,
                Sigma,
                tol = 1e-6,
                empirical = FALSE,
                df = FALSE,
                varnames = NULL,
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
                cl_vars) {
  # main function----------------------------------------------------------------
  foo <- function(iter,
                  n,
                  mu,
                  Sigma,
                  tol,
                  empirical,
                  df) {
    out <- mvrnorm(
      n = n,
      mu = mu,
      Sigma = Sigma,
      tol = tol,
      empirical = empirical
    )
    if (df) {
      out <- as.data.frame(out)
    }
    out
  }
  #------------------------------------------------------------------------------
  # mu to a vector of zeros------------------------------------------------------
  if (is.null(mu)) {
    mu <- rep(
      x = 0,
      times = nrow(Sigma)
    )
    message(
      paste0(
        "mu = NULL. mu is set to a vector of zeroes of length ",
        nrow(Sigma),
        "."
      )
    )
  }
  #------------------------------------------------------------------------------
  # optional column names------------------------------------------------------
  if (!is.null(varnames)) {
    names(mu) <- varnames
  }
  #------------------------------------------------------------------------------
  # negotiate single rep of multiple rep
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
  # execute----------------------------------------------------------------------
  if (single_rep) {
    # single rep-----------------------------------------------------------------
    return(
      foo(
        n = n,
        mu = mu,
        Sigma = Sigma,
        tol = tol,
        empirical = empirical,
        df = df
      )
    )
    # ---------------------------------------------------------------------------
  } else {
    # multiple reps--------------------------------------------------------------
    return(
      par_lapply(
        # iter
        X = 1:R,
        # foo
        FUN = foo,
        # foo
        n = n,
        mu = mu,
        Sigma = Sigma,
        tol = tol,
        empirical = empirical,
        df = df,
        # par_lapply
        par = par,
        ncores = ncores,
        mc = mc,
        lb = lb,
        cl_eval = cl_eval,
        cl_export = cl_export,
        cl_expr = cl_expr,
        cl_vars = cl_vars,
        rbind = NULL # always null
      )
    )
    #----------------------------------------------------------------------------
  }
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Generate Multivariate Normal Data
#'   \eqn{
#'     \mathbf{X} \sim \mathcal{N}_{k}
#'     \left( \boldsymbol{\mu}, \boldsymbol{\Sigma} \right)
#'     %(\#eq:dist-X-mvn)
#'   }
#'   from the Reticular Action Model Matrices
#'
#' @description Generates an \eqn{n \times k} multivariate data matrix
#'   or a list of \eqn{n \times k} multivariate data matrices of length `R`
#'   from the multivariate normal distribution
#'   \deqn{
#'     \mathbf{X} \sim \mathcal{N}_{k}
#'     \left( \boldsymbol{\mu}, \boldsymbol{\Sigma} \right) .
#'     %(\#eq:dist-X-mvn)
#'   }
#'   The model-implied matrices used to generate data
#'   is derived from the Reticular Action Model (RAM) Matrices.
#'
#' @details The multivariate normal distribution has two parameters,
#'   namely the \eqn{k \times 1} mean vector
#'   `mu` \eqn{\left( \boldsymbol{\mu} \right)}
#'   and the \eqn{k \times k} variance-covariance matrix
#'   `Sigma` \eqn{\left( \boldsymbol{\Sigma} \right)}.
#'
#'   The mean vector `mu` can be supplied directly using the `mu` argument.
#'   It can also be derived using `M`.
#'   Note that the argument `mu` takes precedence over `M`.
#'   If `mu` is not provided, it is computed using `M`
#'   with the [`jeksterslabRsem::ram_mutheta()`] function.
#'   If both `mu` and `M` are not provided,
#'   `mu` is set to a vector of zeroes with the appropriate length.
#'
#'   The variance-covariance matrix `Sigma` is derived
#'   from the RAM matrices `A`, `S`, `F`, and `I`.
#'
#'   `mu` and `Sigma` are then used by the [`mvn()`] function
#'   to generate multivariate normal data.
#'
#'   Options for explicit parallelism are provided when `R > 1`
#'   especially when `R` is large. See `par` and suceeding arguments.
#'
#' @seealso
#'   [`jeksterslabRsem::ram_Sigmatheta()`], and [`jeksterslabRsem::ram_mutheta()`]
#'   for more information on the Reticular Action Model (RAM)
#'   and [`mvn()`] for multivariate normal data generation.
#'
#' @family multivariate data functions
#' @keywords multivariate normal
#' @inheritParams mvn
#' @inheritParams jeksterslabRsem::ram_mutheta
#' @inheritParams jeksterslabRsem::ram_Sigmatheta
#' @inherit mvn return
#' @inherit jeksterslabRsem::ram_Sigmatheta references
#' @importFrom jeksterslabRsem ram_mutheta
#' @importFrom jeksterslabRsem ram_Sigmatheta
#' @examples
#' mu <- c(100, 100, 100)
#' A <- matrix(
#'   data = c(0, sqrt(0.26), 0, 0, 0, sqrt(0.26), 0, 0, 0),
#'   ncol = 3
#' )
#' S <- diag(c(225, 166.5, 116.5))
#' F <- I <- diag(3)
#' X <- mvnram(n = 100, mu = mu, A = A, S = S, F = F, I = I)
#' Xstar <- mvnram(n = 100, mu = mu, A = A, S = S, F = F, I = I, R = 100)
#' str(Xstar, list.len = 6)
#' @export
mvnram <- function(n,
                   mu = NULL,
                   M = NULL,
                   A,
                   S,
                   F,
                   I,
                   tol = 1e-6,
                   empirical = FALSE,
                   df = FALSE,
                   varnames = NULL,
                   R = NULL,
                   par = FALSE,
                   ncores = NULL,
                   mc = TRUE,
                   lb = FALSE,
                   cl_eval = FALSE,
                   cl_export = FALSE,
                   cl_expr,
                   cl_vars) {
  # derive Sigma-----------------------------------------------------------------
  Sigma <- ram_Sigmatheta(
    A = A,
    S = S,
    F = F,
    I = I
  )
  #------------------------------------------------------------------------------
  # negotiate derivation of mu---------------------------------------------------
  if (is.null(mu)) {
    # mu is not provided
    if (is.null(M)) {
      message(
        paste0(
          "mu = NULL and M = NULL. mu is set to a vector of zeroes of length ",
          nrow(Sigma),
          "."
        )
      )
      # both mu and M are not provided
      # generate mu from zero
      generate_mu_from_M <- FALSE
      generate_mu_from_zero <- TRUE
    } else {
      message(
        "mu = NULL. mu is computed using M."
      )
      # mu is not provided
      # M is provided
      # generate mu from M
      generate_mu_from_M <- TRUE
      generate_mu_from_zero <- FALSE
    }
  } else {
    # mu is provided
    # even if M is provided mu takes precedence
    generate_mu_from_M <- FALSE
    generate_mu_from_zero <- FALSE
  }
  # derive mu from M-------------------------------------------------------------
  if (generate_mu_from_M) {
    mu <- ram_mutheta(
      A = A,
      F = F,
      I = I,
      M = M
    )
  }
  #------------------------------------------------------------------------------
  # generate mu from zeros-------------------------------------------------------
  if (generate_mu_from_zero) {
    mu <- rep(
      x = 0,
      times = nrow(Sigma)
    )
  }
  #------------------------------------------------------------------------------
  # execute----------------------------------------------------------------------
  mvn(
    # mvrnorm
    n = n,
    mu = mu,
    Sigma = Sigma,
    tol = tol,
    empirical = empirical,
    df = df,
    varnames = varnames,
    # iter
    R = R,
    # par_lapply
    par = par,
    ncores = ncores,
    mc = mc,
    lb = lb,
    cl_eval = cl_eval,
    cl_export = cl_export,
    cl_expr = cl_expr,
    cl_vars = cl_vars
  )
  #------------------------------------------------------------------------------
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Generate Multivariate Normal Data
#'   \eqn{
#'     \mathbf{X} \sim \mathcal{N}_{k}
#'     \left( \boldsymbol{\mu}, \boldsymbol{\Sigma} \right)
#'     %(\#eq:dist-X-mvn)
#'   }
#'   from the Reticular Action Model  Matrices
#'   \eqn{\mathrm{A}}, \eqn{\mathrm{F}}, \eqn{\mathrm{I}},
#'   and the vector of variable variances \eqn{\sigma^2}
#'
#' @details The \eqn{\mathbf{S}} matrix is derived
#'   from a vector of variances sigma2 \eqn{\left( \sigma^2 \right)}
#'   and the proceeds to generating data using the [`mvnram()`] function.
#'   **The first element in sigma2 should be the variance of an exogenous variable.**
#'
#' @family multivariate data functions
#' @keywords multivariate normal
#' @inheritParams mvnram
#' @inheritParams jeksterslabRsem::ram_S
#' @inherit mvnram description return
#' @importFrom jeksterslabRsem ram_S
#' @examples
#' mu <- c(100, 100, 100)
#' A <- matrix(
#'   data = c(0, sqrt(0.26), 0, 0, 0, sqrt(0.26), 0, 0, 0),
#'   ncol = 3
#' )
#' sigma2 <- c(225, 225, 225)
#' F <- I <- diag(3)
#' X <- mvnramsigma2(n = 100, mu = mu, A = A, sigma2 = sigma2, F = F, I = I)
#' str(X)
#' Xstar <- mvnramsigma2(n = 100, mu = mu, A = A, sigma2 = sigma2, F = F, I = I, R = 100)
#' str(Xstar, list.len = 6)
#' @export
mvnramsigma2 <- function(n,
                         mu = NULL,
                         M = NULL,
                         A,
                         sigma2,
                         F,
                         I,
                         tol = 1e-6,
                         empirical = FALSE,
                         df = FALSE,
                         varnames = NULL,
                         R = NULL,
                         par = FALSE,
                         ncores = NULL,
                         mc = TRUE,
                         lb = FALSE,
                         cl_eval = FALSE,
                         cl_export = FALSE,
                         cl_expr,
                         cl_vars) {
  S <- ram_S(
    A = A,
    sigma2 = sigma2,
    F = F,
    I = I,
    SigmaMatrix = FALSE
  )
  mvnram(
    n = n,
    A = A,
    S = S,
    F = F,
    I = I,
    M = M,
    mu = mu,
    tol = tol,
    empirical = empirical,
    df = df,
    varnames = varnames,
    R = R,
    par = par,
    ncores = ncores,
    mc = mc,
    lb = lb,
    cl_eval = cl_eval,
    cl_export = cl_export,
    cl_expr = cl_expr,
    cl_vars = cl_vars
  )
}
