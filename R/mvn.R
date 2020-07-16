#' Generate Multivariate Normal Data
#' \eqn{
#'   \mathbf{X}
#'   \sim
#'   \mathcal{N}_{k}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'   \right)
#'   %(\#eq:dist-X-mvn)
#' }
#'
#' @description Generates an \eqn{n \times k} multivariate data matrix
#' or a list of \eqn{n \times k} multivariate data matrices of length `R`
#' from the multivariate normal distribution
#' \deqn{
#'   \mathbf{X}
#'   \sim
#'   \mathcal{N}_{k}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'   \right) .
#'   %(\#eq:dist-X-mvn)
#' }
#' This function is a wrapper around
#' [`MASS::mvrnorm()`].
#'
#' @details The multivariate normal distribution has two parameters,
#' namely the
#' \eqn{k \times 1} mean vector
#' `mu` \eqn{\left( \boldsymbol{\mu} \right)}
#' and the
#' \eqn{k \times k} variance-covariance matrix
#' `Sigma` \eqn{\left( \boldsymbol{\Sigma} \right)}.
#' If `mu` is not provided,
#' it is set to a vector of zeroes with the appropriate length.
#'
#' Options for explicit parallelism are provided
#' when `R > 1` especially when `R` is large.
#' See `par` and suceeding arguments.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family multivariate data functions
#' @keywords multivariate, normal
#' @inheritParams jeksterslabRdist::mvnpdf
#' @inheritParams univ
#' @inherit jeksterslabRdist::mvnpdf details references
#' @param tol Numeric.
#' Tolerance (relative to largest variance)
#' for numerical lack of positive-definiteness in `Sigma`.
#' @param empirical Logical.
#' If `TRUE`, `mu` and `Sigma` specify the empirical
#' not population mean and covariance matrix.
#' @return
#' If `R = NULL` or `R = 1`,
#' returns an \eqn{n \times k} multivariate normal data matrix generated
#' using the mean vector `mu`
#' \eqn{\left( \boldsymbol{\mu} \right)}
#' and the variance-covariance matrix `Sigma`
#' \eqn{\left( \boldsymbol{\Sigma} \right)}
#' provided.
#' If `R` is an integer greater than 1,
#' (e.g., `R = 10`)
#' returns a list of length `R` of
#' \eqn{n \times k} multivariate normal data matrices generated
#' using the mean vector `mu`
#' and the variance-covariance matrix
#' `Sigma` provided.
#' @examples
#' Sigma <- matrix(
#'   data = c(
#'     225, 112.50, 56.25,
#'     112.5, 225, 112.5,
#'     56.25, 112.50, 225
#'   ),
#'   ncol = 3
#' )
#' mu <- c(100, 100, 100)
#' # single random data set------------------------------------------------------
#' X <- mvn(
#'   n = 100,
#'   mu = mu,
#'   Sigma = Sigma
#' )
#' str(X, list.len = 6)
#' # multiple random data sets---------------------------------------------------
#' X <- mvn(
#'   n = 100,
#'   mu = mu,
#'   Sigma = Sigma,
#'   R = 10
#' )
#' str(X, list.len = 6)
#' @seealso [`jeksterslabRdist::mvnpdf()`], [`jeksterslabRdist::mvnll()`],
#' and [`jeksterslabRdist::mvn2ll()`],
#' for more information on the multivariate normal distribution.
#' @references
#' Venables, W. N., Ripley, B. D., & Venables, W. N. (2002).
#' *Modern applied statistics with S*.
#' New York, N.Y: Springer.
#'
#' [Wikipedia: Multivariate normal distribution](https://en.wikipedia.org/wiki/Multivariate_normal_distribution)
#' @importFrom MASS mvrnorm
#' @export
mvn <- function(n, # mvrnorm
                mu = NULL,
                Sigma,
                tol = 1e-6,
                empirical = FALSE,
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
                  empirical) {
    mvrnorm(
      n = n,
      mu = mu,
      Sigma = Sigma,
      tol = tol,
      empirical = empirical
    )
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
        empirical = empirical
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

#' Generate Multivariate Normal Data
#' \eqn{
#'   \mathbf{X}
#'   \sim
#'   \mathcal{N}_{k}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'   \right)
#'   %(\#eq:dist-X-mvn)
#' }
#' from the Reticular Action Model (RAM) Matrices
#'
#' @description Generates an \eqn{n \times k} multivariate data matrix
#' or a list of \eqn{n \times k} multivariate data matrices of length `R`
#' from the multivariate normal distribution
#' \deqn{
#'   \mathbf{X}
#'   \sim
#'   \mathcal{N}_{k}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'   \right)
#'   %(\#eq:dist-X-mvn)
#' }
#' The model-implied matrices used to generate data
#' is derived from the Reticular Action Model (RAM) Matrices.
#'
#' @details The multivariate normal distribution has two parameters,
#' namely the
#' \eqn{k \times 1} mean vector
#' `mu` \eqn{\left( \boldsymbol{\mu} \right)}
#' and the
#' \eqn{k \times k} variance-covariance matrix
#' `Sigma` \eqn{\left( \boldsymbol{\Sigma} \right)}.
#'
#' The mean vector `mu`
#' \eqn{\left( \boldsymbol{\mu} \right)}
#' can be supplied directly
#' using the `mu` argument.
#' It can also be derived
#' using `M`.
#' Note that the argument `mu`
#' \eqn{\left( \boldsymbol{\mu} \right)}
#' takes precedence over `M`
#' If `mu` is not provided,
#' it is computed using `M`
#' with the [`jeksterslabRsem::ram_mutheta()`] function.
#' If both `mu` and `M`
#' are not provided,
#' `mu` is set to a vector of zeroes with the appropriate length.
#'
#' The variance-covariance matrix `Sigma`
#' \eqn{\left(\boldsymbol{\Sigma}\right)}
#' is derived
#' from the RAM matrices `A`, `S`, `F`, and `I`.
#'
#' `mu` \eqn{\left( \boldsymbol{\mu} \right)}
#' and `Sigma` \eqn{\left(\boldsymbol{\Sigma}\right)}
#' are then used by the [`mvn()`]
#' to generate multivariate normal data.
#'
#' Options for explicit parallelism are provided
#' when `R > 1` especially when `R` is large.
#' See `par` and suceeding arguments.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family multivariate data functions
#' @keywords multivariate normal
#' @inheritParams mvn
#' @inheritParams jeksterslabRsem::ram_mutheta
#' @inheritParams jeksterslabRsem::ram_Sigmatheta
#' @inherit jeksterslabRsem::ram_Sigmatheta references
#' @return
#' If `R = NULL` or `R = 1`,
#' returns an \eqn{n \times k} multivariate normal data matrix generated
#' using the provided or derived mean vector `mu`
#' \eqn{\left( \boldsymbol{\mu} \right)}
#' and the derived variance-covariance matrix `Sigma`
#' \eqn{\left( \boldsymbol{\Sigma} \right)}.
#' If `R` is an integer greater than 1,
#' (e.g., `R = 10`)
#' returns a list of length `R` of
#' \eqn{n \times k} multivariate normal data matrices generated
#' using the provided or derived mean vector `mu`
#' and the derived variance-covariance matrix `Sigma`.
#' @examples
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' S <- F <- I <- diag(3)
#' S[1, 1] <- 225
#' S[2, 2] <- 166.5
#' S[3, 3] <- 166.5
#' mu <- c(100, 100, 100)
#' # single random data set------------------------------------------------------
#' X <- mvnram(
#'   n = 100,
#'   A = A,
#'   S = S,
#'   F = F,
#'   I = I,
#'   mu = mu
#' )
#' str(X, list.len = 6)
#' # multiple random data sets---------------------------------------------------
#' X <- mvnram(
#'   n = 100,
#'   A = A,
#'   S = S,
#'   F = F,
#'   I = I,
#'   mu = mu,
#'   R = 10
#' )
#' str(X, list.len = 6)
#' @seealso
#'   [`jeksterslabRsem::ram_Sigmatheta()`], and [`jeksterslabRsem::ram_mutheta()`]
#'   for more information on the Reticular Action Model (RAM)
#'   and [`mvn()`] for multivariate normal data generation.
#' @importFrom jeksterslabRsem ram_mutheta
#' @importFrom jeksterslabRsem ram_Sigmatheta
#' @export
mvnram <- function(n,
                   A,
                   S,
                   F,
                   I,
                   M = NULL,
                   mu = NULL,
                   tol = 1e-6,
                   empirical = FALSE,
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

#' Generate Multivariate Normal Data
#' \eqn{
#'   \mathbf{X}
#'   \sim
#'   \mathcal{N}_{k}
#'   \left(
#'     \boldsymbol{\mu},
#'     \boldsymbol{\Sigma}
#'   \right)
#'   %(\#eq:dist-X-mvn)
#' }
#' from the Reticular Action Model (RAM) Matrices.
#'
#' The \eqn{\mathbf{S}} matrix
#' is derived from a vector of variances.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family multivariate data functions
#' @keywords multivariate normal
#' @inheritParams mvnram
#' @inheritParams jeksterslabRsem::ram_S
#' @inherit mvnram return
#' @examples
#' # One-factor CFA model--------------------------------------------------------
#' A <- matrix(data = 0, ncol = 6, nrow = 6)
#' for (i in 2:6) {
#'   A[i, 1] <- 0.5
#' }
#' sigma2 <- c(1, 1, 1, 1, 1, 1)
#' F <- I <- diag(nrow(A))
#' F <- diag(nrow(A) - 1)
#' F <- cbind(0, F)
#' X <- mvnramsigma2(
#'   n = 100,
#'   A = A,
#'   sigma2 = sigma2,
#'   F = F,
#'   I = I,
#'   empirical = TRUE
#' )
#' str(X, list.len = 6)
#' cov(X)
#' @importFrom jeksterslabRsem ram_S
#' @export
mvnramsigma2 <- function(n,
                         A,
                         sigma2,
                         F,
                         I,
                         M = NULL,
                         mu = NULL,
                         tol = 1e-6,
                         empirical = FALSE,
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
