#' ---
#' title: "Tests for the mvn(), mvnram(), and mvnramsigma2 functions"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' description: >
#'   Tests for the mvn(), mvnram(), and mvnramsigma2 functions.
#' output:
#'   rmarkdown::html_vignette:
#'     toc: true
#' vignette: >
#'   %\VignetteIndexEntry{Tests for the mvn(), mvnram(), and mvnramsigma2 functions}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
#'
#+ include = FALSE, cache = FALSE
knitr::opts_chunk$set(
  error = TRUE,
  collapse = TRUE,
  comment = "#>",
  result.width = "100%"
)
#'
#+ plot-arguments, echo = FALSE
par(pty = "s")
breaks <- 100
#'
#+ setup
library(testthat)
library(jeksterslabRsem)
library(jeksterslabRdata)
context("Test univ Normal Distribution.")
#'
#' ## Multivariate Normal Distribution
#'
#' ### Sample Size and Monte Carlo replications
#'
#+ sizes, echo = FALSE
n <- 1000
R <- 1000
breaks
Variable <- c(
  "`n`",
  "`R`"
)
Description <- c(
  "Sample size.",
  "Monte Carlo replications."
)
Notation <- c(
  "$n$",
  "$R$"
)
Values <- c(
  n,
  R
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Notation,
    Values
  ),
  caption = "Sample Size and Monte Carlo Replications"
)
#'
#' ### Population Parameters
#'
#+ theta, echo = FALSE
varnames <- c("x", "y", "z")
mu <- c(70.18, 3.06, 3.24)
names(mu) <- varnames
M <- c(70.18, -20.70243, -12.71288)
Sigma <- matrix(
  data = c(
    1.2934694, 0.4379592, 0.4661224, 0.4379592,
    1.0779592, 0.5771429, 0.4661224, 0.5771429, 1.2881633
  ),
  nrow = 3
)
colnames(Sigma) <- varnames
rownames(Sigma) <- varnames
Sigma_vector <- as.vector(Sigma)
A <- matrix(
  data = c(
    0, 0.3385926, 0.2076475, 0, 0, 0.4510391, 0, 0, 0
  ),
  nrow = 3
)
S <- diag(c(1.2934694, 0.9296694, 0.9310601))
sigma2 <- c(1.2934694, 1.0779592, 1.2881633)
F <- I <- diag(3)
Sigmaram <- ram_Sigmatheta(
  A = A,
  S = S,
  F = F,
  I = I
)
knitr::kable(
  x = mu,
  caption = "Mean Vector $\\mu$"
)
knitr::kable(
  x = Sigma,
  caption = "Variance-Covariance Matrix $\\Sigma$"
)
#'
#' ### Monte Carlo Simulation
#'
#+ monte_carlo, echo = FALSE
Xstar <- mvn(n = n, mu = mu, Sigma = Sigma, R = R)
foo <- function(X) {
  as.vector(cov(X))
}
thetahatstar <- fit(Xstar = Xstar, fitFUN = foo, rbind = TRUE)
hist(thetahatstar[, 1],
  main = expression("Sampling distribution of" ~ hat(sigma)^2 ~ (x)),
  xlab = expression(hat(sigma)^2 ~ (x)),
  freq = FALSE,
  breaks = breaks
)
qqnorm(thetahatstar[, 1])
qqline(thetahatstar[, 1])
hist(thetahatstar[, 2],
  main = expression("Sampling distribution of" ~ hat(sigma)(list(x, y))),
  xlab = expression(hat(sigma)(list(x, y))),
  freq = FALSE,
  breaks = breaks
)
qqnorm(thetahatstar[, 2])
qqline(thetahatstar[, 2])
hist(thetahatstar[, 3],
  main = expression("Sampling distribution of" ~ hat(sigma)(list(x, z))),
  xlab = expression(hat(sigma)(list(x, z))),
  freq = FALSE,
  breaks = breaks
)
qqnorm(thetahatstar[, 3])
qqline(thetahatstar[, 3])
hist(thetahatstar[, 5],
  main = expression("Sampling distribution of" ~ hat(sigma)^2 ~ (y)),
  xlab = expression(hat(sigma)^2 ~ (y)),
  freq = FALSE,
  breaks = breaks
)
qqnorm(thetahatstar[, 5])
qqline(thetahatstar[, 5])
hist(thetahatstar[, 6],
  main = expression("Sampling distribution of" ~ hat(sigma)(list(y, z))),
  xlab = expression(hat(sigma)(list(y, z))),
  freq = FALSE,
  breaks = breaks
)
qqnorm(thetahatstar[, 6])
qqline(thetahatstar[, 6])

hist(thetahatstar[, 9],
  main = expression("Sampling distribution of" ~ hat(sigma)^2 ~ (z)),
  xlab = expression(hat(sigma)^2 ~ (z)),
  freq = FALSE,
  breaks = breaks
)
qqnorm(thetahatstar[, 9])
qqline(thetahatstar[, 9])

Sigmahatthetahatstar <- colMeans(thetahatstar)
test_that("theta", {
  expect_equivalent(
    round(Sigma_vector, digits = 1),
    round(Sigmahatthetahatstar, digits = 1)
  )
})
#'
#+ Other Tests, echo = FALSE
# matrix
set.seed(42)
X_mvn_matrix <- mvn(n = n, mu = mu, Sigma = Sigma, varnames = varnames)
X_mvn_matrix_cov <- as.vector(cov(X_mvn_matrix))
set.seed(42)
X_mvnram_mu_matrix <- mvnram(n = n, mu = mu, A = A, S = S, F = F, I = I, varnames = varnames)
X_mvnram_mu_matrix_cov <- as.vector(cov(X_mvnram_mu_matrix))
set.seed(42)
X_mvnram_M_matrix <- mvnram(n = n, M = M, A = A, S = S, F = F, I = I, varnames = varnames)
X_mvnram_M_matrix_cov <- as.vector(cov(X_mvnram_M_matrix))
set.seed(42)
X_mvnramsigma2_mu_matrix <- mvnramsigma2(n = n, mu = mu, A = A, sigma2 = sigma2, F = F, I = I, varnames = varnames)
X_mvnramsigma2_mu_matrix_cov <- as.vector(cov(X_mvnramsigma2_mu_matrix))
set.seed(42)
X_mvnramsigma2_M_matrix <- mvnramsigma2(n = n, M = M, A = A, sigma2 = sigma2, F = F, I = I, varnames = varnames)
X_mvnramsigma2_M_matrix_cov <- as.vector(cov(X_mvnramsigma2_M_matrix))
# df
set.seed(42)
X_mvn_df <- mvn(n = n, mu = mu, Sigma = Sigma, empirical = TRUE, df = TRUE, varnames = varnames)
X_mvn_df_cov <- as.vector(cov(X_mvn_df))
set.seed(42)
X_mvnram_mu_df <- mvnram(n = n, mu = mu, A = A, S = S, F = F, I = I, empirical = TRUE, df = TRUE, varnames = varnames)
X_mvnram_mu_df_cov <- as.vector(cov(X_mvnram_mu_df))
set.seed(42)
X_mvnram_M_df <- mvnram(n = n, M = M, A = A, S = S, F = F, I = I, empirical = TRUE, df = TRUE, varnames = varnames)
X_mvnram_M_df_cov <- as.vector(cov(X_mvnram_M_df))
set.seed(42)
X_mvnramsigma2_mu_df <- mvnramsigma2(n = n, mu = mu, A = A, sigma2 = sigma2, F = F, I = I, empirical = TRUE, df = TRUE, varnames = varnames)
X_mvnramsigma2_mu_df_cov <- as.vector(cov(X_mvnramsigma2_mu_df))
set.seed(42)
X_mvnramsigma2_M_df <- mvnramsigma2(n = n, M = M, A = A, sigma2 = sigma2, F = F, I = I, empirical = TRUE, df = TRUE, varnames = varnames)
X_mvnramsigma2_M_df_cov <- as.vector(cov(X_mvnramsigma2_M_df))
# sum elemens of covariance matrix should be equivalent
test_that("matrix", {
  expect_equivalent(
    sum(X_mvn_matrix_cov),
    sum(X_mvnram_mu_matrix_cov),
    sum(X_mvnram_M_matrix_cov),
    sum(X_mvnramsigma2_mu_matrix_cov),
    sum(X_mvnramsigma2_M_matrix_cov)
  )
})
test_that("df", {
  expect_equivalent(
    sum(X_mvn_df_cov),
    sum(X_mvnram_mu_df_cov),
    sum(X_mvnram_M_df_cov),
    sum(X_mvnramsigma2_mu_df_cov),
    sum(X_mvnramsigma2_M_df_cov)
  )
})
