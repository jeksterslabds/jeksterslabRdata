#' Generate Multivariate Normal Data
#' \eqn{\mathbf{X} \sim \mathcal{N}_{k} \left( \boldsymbol{\mu}, \boldsymbol{\Sigma}\right)}
#'
#' Generates an \eqn{n \times k}
#' multivariate normal data
#' matrix
#' \eqn{\mathbf{X}}
#' from a
#' \eqn{k \times 1} mean vector
#' \eqn{\boldsymbol{\mu}}
#' and a
#' \eqn{k \times k} variance-covariance matrix
#' \eqn{\boldsymbol{\Sigma}}.
#' This function is a wrapper around
#' [`MASS::mvrnorm()`].
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family multivariate data functions
#' @keywords multivariate normal
#' @inheritParams jeksterslabRdist::mvnpdf
#' @inheritParams univ
#' @inherit jeksterslabRdist::mvnpdf details references
#' @param ... Arguments that can be passed to [`MASS::mvrnorm()`].
#' @return
#' If `R = NULL`,
#' returns an \eqn{n \times k} multivariate normal data matrix generated
#' using the mean vector
#' `mu`
#' \eqn{\left( \boldsymbol{\mu} \right)}
#' and the variance-covariance matrix
#' `Sigma`
#' \eqn{\left( \boldsymbol{\Sigma} \right)}
#' provided.
#' If `R` is an integer,
#' returns a list of length `R` of
#' \eqn{n \times k} multivariate normal data matrices generated
#' using the mean vector
#' `mu``
#' and the variance-covariance matrix
#' `Sigma`
#' provided.
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
#' # single random data set
#' X <- mvn(
#'   n = 100,
#'   mu = mu,
#'   Sigma = Sigma
#' )
#' str(X, list.len = 6)
#' # multiple random data sets
#' X <- mvn(
#'   n = 100,
#'   mu = mu,
#'   Sigma = Sigma,
#'   R = 10
#' )
#' str(X, list.len = 6)
#' @seealso
#' [`jeksterslabRdist::mvnpdf()`],
#' [`jeksterslabRdist::mvnll()`],
#' and
#' [`jeksterslabRdist::mvn2ll()`],
#' for more information on the multivariate normal distribution.
#' @references
#' Venables, W. N., Ripley, B. D., & Venables, W. N. (2002).
#' *Modern applied statistics with S*.
#' New York, N.Y: Springer.
#'
#' [Wikipedia: Multivariate normal distribution](https://en.wikipedia.org/wiki/Multivariate_normal_distribution)
#' @importFrom MASS mvrnorm
#' @export
mvn <- function(n,
                mu = NULL,
                Sigma,
                R = NULL,
                par = FALSE,
                ncores = NULL,
                ...) {
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
  foo <- function(iter,
                  n,
                  mu,
                  Sigma,
                  ...) {
    mvrnorm(
      n = n,
      mu = mu,
      Sigma = Sigma,
      ...
    )
  }
  if (is.null(R)) {
    out <- foo(
      n = n,
      mu = mu,
      Sigma = Sigma,
      ...
    )
  } else {
    args <- list(
      iter = 1:R,
      n = n,
      mu = mu,
      Sigma = Sigma,
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

#' Generate Multivariate Normal Data from RAM Matrices
#'
#' Generates multivariate normal data from
#' the Reticular Action Model (RAM) matrices.
#'
#' The multivariate normal distribution has two parameters,
#' namely the
#' \eqn{k \times 1} mean vector
#' `mu`
#' \eqn{\left( \boldsymbol{\mu} \right)}
#' and the
#' \eqn{k \times k} variance-covariance matrix
#' `Sigma`
#' \eqn{\left( \boldsymbol{\Sigma} \right)}.
#'
#' The mean vector
#' `mu`
#' \eqn{\left( \boldsymbol{\mu} \right)}
#' can be supplied directly using the
#' `mu` argument.
#' It can also be derived using
#' `M` .
#' Note that the argument
#' `mu`
#' \eqn{\left( \boldsymbol{\mu} \right)}
#' takes precedence over
#' `M`
#' If
#' `mu`
#' is not provided,
#' it is computed using
#' `M`
#' with the
#' [`jeksterslabRsem::ram_mutheta()`]
#' function.
#' If both
#' `mu`
#' and
#' `M`
#' are not provided,
#' `mu`
#' is set to a vector of zeroes with the appropriate length.
#'
#' The variance-covariance matrix `Sigma`
#' \eqn{\left(\boldsymbol{\Sigma}\right)}
#' is derived
#' from the RAM matrices `A`, `S`, `F`, and `I`.
#'
#' `mu`
#' \eqn{\left( \boldsymbol{\mu} \right)}
#' and
#' `Sigma`
#' \eqn{\left(\boldsymbol{\Sigma}\right)}
#' are then used
#' by the
#' [`mvn()`]
#' to generate multivariate normal data.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family multivariate data functions
#' @keywords multivariate normal
#' @inheritParams mvn
#' @inheritParams jeksterslabRsem::ram_mutheta
#' @inheritParams jeksterslabRsem::ram_Sigmatheta
#' @inherit jeksterslabRsem::ram_Sigmatheta references
#' @return
#' If `R = NULL`,
#' returns an \eqn{n \times k} multivariate normal data matrix generated
#' using the provided or derived mean vector
#' `mu`
#' \eqn{\left( \boldsymbol{\mu} \right)}
#' and the derived variance-covariance matrix
#' `Sigma`
#' \eqn{\left( \boldsymbol{\Sigma} \right)}.
#' If `R` is an integer,
#' returns a list of length `R` of
#' \eqn{n \times k} multivariate normal data matrices generated
#' using the provided or derived mean vector
#' `mu`
#' and the derived variance-covariance matrix
#' `Sigma`.
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
#' # single random data set
#' X <- mvnram(
#'   n = 100,
#'   A = A,
#'   S = S,
#'   F = F,
#'   I = I,
#'   mu = mu
#' )
#' str(X, list.len = 6)
#' # multiple random data sets
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
#' [`jeksterslabRsem::ram_Sigmatheta()`],
#' and
#' [`jeksterslabRsem::ram_mutheta()`]
#' for more information on the Reticular Action Model
#' and
#' [`mvn()`]
#' for multivariate normal data generation.
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
                   R = NULL,
                   par = FALSE,
                   ncores = NULL,
                   ...) {
  # derive Sigma
  Sigma <- ram_Sigmatheta(
    A = A,
    S = S,
    F = F,
    I = I
  )
  # derive mu
  if (is.null(mu)) {
    if (is.null(M)) {
      message(
        paste0(
          "mu = NULL and M = NULL. mu is set to a vector of zeroes of length ",
          nrow(Sigma),
          "."
        )
      )
      mu <- rep(
        x = 0,
        times = nrow(Sigma)
      )
    } else {
      message(
        "mu = NULL. mu is computed using M."
      )
      mu <- ram_mutheta(
        A = A,
        F = F,
        I = I,
        M = M
      )
    }
  }
  # generate data
  mvn(
    n = n,
    mu = mu,
    Sigma,
    R = R,
    par = par,
    ncores = ncores,
    ...
  )
}
