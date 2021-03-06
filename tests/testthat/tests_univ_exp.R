#' ---
#' title: "Tests for the univ() function - Exponential Distribution"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' description: >
#'   Tests for the univ() function - Exponential Distribution.
#' output:
#'   rmarkdown::html_vignette:
#'     toc: true
#' vignette: >
#'   %\VignetteIndexEntry{Tests for the univ() function - Exponential Distribution}
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
#+ setup
library(testthat)
library(MASS)
library(jeksterslabRdata)
context("Test univ Exponential Distribution.")
#'
#+ plot-arguments, echo = FALSE
par(pty = "s")
breaks <- 100
#'
#' ## Univariate Data Generation - Exponential Distribution
#'
#' ### Sample Size and Monte Carlo replications
#'
#+ sizes, echo = FALSE
n <- 1000
R <- 1000
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
lambda <- 1
mu <- 1 / lambda
sigma2 <- 1 / lambda^2
sigma <- sqrt(sigma2)
Variable <- c(
  "`lambda`",
  "`mu`",
  "`sigma2`",
  "`sigma`"
)
Description <- c(
  "Population rate.",
  "Population mean.",
  "Population variance.",
  "Population standard deviation."
)
Notation <- c(
  "$\\lambda$",
  "$\\mu$",
  "$\\sigma^2$",
  "$\\sigma$"
)
Values <- c(
  lambda,
  mu,
  sigma2,
  sigma
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Notation,
    Values
  ),
  caption = "Population Parameters"
)
#'
#' ### Sample Data
#'
#' \begin{equation}
#'   X
#'   \sim
#'   \mathrm{exp}
#'   \left(
#'     \lambda
#'   \right)
#'   %(\#eq:dist-exp-X-r)
#' \end{equation}
#'
#' The random variable $X$
#' takes on the value $x$
#' for each observation $\omega \in \Omega$.
#'
#' - $\omega$ refers to units or observations
#' - $\Omega$ refers to the collection of all units,
#'   that is, a set of possible outcomes or the sample space
#'   contained in the set $S$.
#'
#' A random variable acts as a function,
#' inasmuch as,
#' it maps each observation $\omega \in \Omega$
#' to a value $x$.
#'
#' \begin{equation}
#'   x
#'   =
#'   X
#'   \left(
#'     \omega
#'   \right)
#'   %(\#eq:dist-random-variable-1)
#' \end{equation}
#'
#' \begin{equation}
#'   \mathbf{x}
#'   =
#'   \begin{bmatrix}
#'     x_1 = X \left( \omega_1 \right) \\
#'     x_2 = X \left( \omega_2 \right) \\
#'     x_3 = X \left( \omega_3 \right) \\
#'     x_i = X \left( \omega_i \right) \\
#'     \vdots \\
#'     x_n = X \left( \omega_n \right)
#'   \end{bmatrix}, \\
#'   i
#'   =
#'   \left\{
#'     1, 2, 3, \dots, n
#'   \right\}
#'   %(\#eq:dist-random-variable-2)
#' \end{equation}
#'
#'
#+ single_sample, echo = FALSE
x <- univ(
  n = n,
  rFUN = rexp,
  rate = lambda
)
hist(
  x,
  main = expression(
    paste(
      "Histogram of ",
      x
    )
  ),
  xlab = expression(
    x
  ),
  freq = FALSE
)
qqnorm(x)
qqline(x)
plot(
  ecdf(x),
  main = "Empirical Cumulative Density Function",
  xlab = "x"
)
#'
#' ### Monte Carlo Simulation
#'
#+ monte_carlo, echo = FALSE
x_k <- univ(
  n = n,
  rFUN = rexp,
  R = R,
  rate = lambda
)
fitlambda <- function(x) {
  out <- fitdistr(x, "exponential")
  out[["estimate"]]
}
functions <- list(
  fitlambda,
  mean,
  var,
  sd
)
out <- list(
  lambdas = rep(x = NA, times = n),
  means = rep(x = NA, times = n),
  vars = rep(x = NA, times = n),
  sds = rep(x = NA, times = n)
)
sds <- vars <- means <- lambdas <- rep(
  x = NA,
  times = length(out)
)
for (i in 1:length(out)) {
  out[[i]] <- fit(
    Xstar = x_k,
    fitFUN = functions[[i]],
    par = FALSE
  )
}
for (i in 1:length(out)) {
  means[i] <- mean(out[[i]])
  vars[i] <- var(out[[i]])
  sds[i] <- sd(out[[i]])
}
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Notation,
    Values,
    means,
    vars,
    sds
  ),
  col.names = c(
    "Variable",
    "Description",
    "Notation",
    "Values",
    "Mean of estimates",
    "Variance of estimates",
    "Standard deviation of estimates"
  ),
  caption = "Monte Carlo Simulation Results"
)
#'
#' ## Sampling Distribution of Sample Lambdas
#'
#+ sampling-dist-lambda, echo = FALSE
hist(
  out[["lambdas"]],
  main = expression(
    paste(
      "Sampling distribution of estimated lambdas ",
      (bold(hat(lambda)))
    )
  ),
  xlab = expression(
    hat(mu)
  ),
  freq = FALSE,
  breaks = breaks
)
qqnorm(out[["lambdas"]])
qqline(out[["lambdas"]])
#'
#' ## Sampling Distribution of Sample Means
#'
#' \begin{equation}
#'   \hat{\mu}
#'   \sim
#'   \mathcal{N}
#'   \left(
#'     \mu,
#'     \frac{\sigma^2}{n}
#'   \right)
#'   %(\#eq:dist-sampling-dist-mean)
#' \end{equation}
#'
#' \begin{equation}
#'   \mathbb{E}
#'   \left[
#'     \hat{\mu}
#'   \right]
#'   =
#'   \mu
#'   %(\#eq:dist-sampling-dist-mean-expected-value)
#' \end{equation}
#'
#' \begin{equation}
#'   \mathrm{Var}
#'   \left(
#'     \hat{\mu}
#'   \right)
#'   =
#'   \frac{\sigma^2}{n}
#'   %(\#eq:dist-sampling-dist-mean-var)
#' \end{equation}
#'
#' \begin{equation}
#'   \mathrm{se}
#'   \left(
#'     \hat{\mu}
#'   \right)
#'   =
#'   \frac{\sigma}{\sqrt{n}}
#'   %(\#eq:dist-sampling-dist-mean-se)
#' \end{equation}
#'
#+ sampling-dist-mean, echo = FALSE
hist(
  out[["means"]],
  main = expression(
    paste(
      "Sampling distribution of estimated means ",
      (bold(hat(mu)))
    )
  ),
  xlab = expression(
    hat(mu)
  ),
  freq = FALSE,
  breaks = breaks
)
qqnorm(out[["means"]])
qqline(out[["means"]])
#'
#' ## Sampling Distribution of Sample Variances
#'
#' \begin{equation}
#'   \hat{\sigma}^2
#'   \sim
#'   \chi^2
#'   \left(
#'     k
#'     =
#'     n - 1
#'   \right),
#'   \quad
#'   \text{when}
#'   \quad
#'   X
#'   \sim
#'   \left(
#'     \mu,
#'     \sigma^2
#'   \right)
#'   %(\#eq:dist-sampling-dist-var)
#' \end{equation}
#'
#' \begin{equation}
#'   \mathbb{E}
#'   \left[
#'     \hat{\sigma}^2
#'   \right]
#'   =
#'   \sigma^2
#'   %(\#eq:dist-sampling-dist-var-expected-value)
#' \end{equation}
#'
#' \begin{equation}
#'   \mathrm{Var}
#'   \left(
#'     \hat{\sigma}^2
#'   \right)
#'   =
#'   \frac{
#'     \left(
#'       n
#'       -
#'       1
#'     \right)
#'     \hat{\sigma}^2
#'   }{
#'     \sigma^2
#'   },
#'   \quad
#'   \text{where}
#'   \quad
#'   X
#'   \sim
#'   \mathcal{N}
#'   \left(
#'     \mu,
#'     \sigma^2
#'   \right)
#'   %(\#eq:dist-sampling-dist-var-var)
#' \end{equation}
#'
#' When $X \sim \mathcal{N} \left( \mu, \sigma^2 \right)$,
#' $\mathrm{Var} \left( \hat{\sigma}^2 \right) \sim \chi^2 \left( k = n - 1 \right)$.
#' As $k$ tends to infinity,
#' the distribution of $\mathrm{Var} \left( \hat{\sigma}^2 \right)$
#' converges to normality.
#' When the $X$ is not normally distributed,
#' the variance of the sampling distribution of the sample variances
#' takes on a slightly different functional form.
#' Increase in sample size,
#' however large it may be,
#' does not help in approximating the normal distribution.
#'
#+ sampling-dist-var, echo = FALSE
hist(
  out[["vars"]],
  main = expression(
    paste(
      "Sampling distribution of estimated variances ",
      (bold(hat(sigma)^2))
    )
  ),
  xlab = expression(
    hat(sigma)^2
  ),
  freq = FALSE,
  breaks = breaks
)
qqnorm(out[["vars"]])
qqline(out[["vars"]])
#'
#+ sampling-dist-sd, echo = FALSE
hist(
  out[["sds"]],
  main = expression(
    paste(
      "Sampling distribution of estimated standard deviations ",
      (bold(hat(sigma)))
    )
  ),
  xlab = expression(
    hat(sigma)
  ),
  freq = FALSE,
  breaks = breaks
)
qqnorm(out[["sds"]])
qqline(out[["sds"]])
#'
#+ testhat01
test_that("expected values (means) of muhat, sigma2hat, sigmahat", {
  expect_equivalent(
    round(
      means,
      digits = 0
    ),
    c(
      lambda,
      mu,
      sigma2,
      sigma
    )
  )
})
