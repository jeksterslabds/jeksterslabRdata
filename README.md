jeksterslabRdata
================
Ivan Jacob Agaloos Pesigan
2020-07-17

<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![Travis build
status](https://travis-ci.com/jeksterslabds/jeksterslabRdata.svg?branch=master)](https://travis-ci.com/jeksterslabds/jeksterslabRdata)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/jeksterslabds/jeksterslabRdata?branch=master&svg=true)](https://ci.appveyor.com/project/jeksterslabds/jeksterslabRdata)
[![codecov](https://codecov.io/github/jeksterslabds/jeksterslabRdata/branch/master/graphs/badge.svg)](https://codecov.io/github/jeksterslabds/jeksterslabRdata)
<!-- badges: end -->

`jeksterslabRdata` is a collection of functions that I find useful in
studying data generation and sampling.

## Installation

You can install the released version of `jeksterslabRdata` from
[GitHub](https://github.com/jeksterslabds/jeksterslabRdata) with:

``` r
library(devtools)
install_github("jeksterslabds/jeksterslabRdata")
```

## Documentation

See [GitHub
Pages](https://jeksterslabds.github.io/jeksterslabRdata/index.html) for
package documentation.

## Main functions

``` r
library(jeksterslabRdata)
```

### `univ()`

Generates an `n x 1` univariate data vector or a list of `n x 1`
univariate data vectors of length `R`. The default data generating
function is the normal distribution.

#### Single Random Data Set

Run the function.

``` r
x <- univ(n = 100, rFUN = rnorm, mean = 100, sd = sqrt(225))
```

Explore the output.

``` r
str(x, list.len = 6)
#>  num [1:100] 86.9 112.7 99.8 104.4 110.6 ...
hist(x, main = expression(italic(N)(list(mu == 100, sigma^2 == 225))))
```

<img src="man/figures/README-unnamed-chunk-4-1.png" width="100%" />

#### Multiple Random Data Sets

Run the function.

``` r
xstar <- univ(n = 100, rFUN = rnorm, mean = 100, sd = sqrt(225), R = 100)
```

Explore the output.

``` r
str(xstar, list.len = 6)
#> List of 100
#>  $ : num [1:100] 86.8 120.9 108 86.7 86.7 ...
#>  $ : num [1:100] 84 93.6 128.6 90.2 99.1 ...
#>  $ : num [1:100] 114.6 81.7 104.1 116 115.9 ...
#>  $ : num [1:100] 87.1 78.1 92 87.2 111.5 ...
#>  $ : num [1:100] 120 106 117 106 105 ...
#>  $ : num [1:100] 95.4 107.1 93.6 117.2 80 ...
#>   [list output truncated]
```

### `mvn()`

Generates an `n x k` multivariate data matrix or a list of `n x k`
multivariate data matrices of length `R` from the multivariate normal
distribution. This function is a wrapper around `MASS::mvrnorm()`.

#### Single Random Data Set

Set `mu` and `Sigma`.

``` r
mu <- c(100, 100, 100)
Sigma <- matrix(
  data = c(225, 112.50, 56.25, 112.5, 225, 112.5, 56.25, 112.50, 225),
  ncol = 3
)
```

Run the function.

``` r
X <- mvn(n = 100, mu = mu, Sigma = Sigma)
```

Explore the output.

``` r
str(X)
#>  num [1:100, 1:3] 116.7 105 88.7 71.5 92.3 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : NULL
pairs(X)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

``` r
colMeans(X)
#> [1]  99.94866 101.27051  98.96092
cov(X)
#>          [,1]      [,2]      [,3]
#> [1,] 237.9677 113.91144  37.26470
#> [2,] 113.9114 240.14685  94.42719
#> [3,]  37.2647  94.42719 187.28192
cor(X)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.4765077 0.1765187
#> [2,] 0.4765077 1.0000000 0.4452569
#> [3,] 0.1765187 0.4452569 1.0000000
```

#### Multiple Random Data Sets

Run the function.

``` r
Xstar <- mvn(n = 100, mu = mu, Sigma = Sigma, R = 100)
```

Explore the output.

``` r
str(Xstar, list.len = 6)
#> List of 100
#>  $ : num [1:100, 1:3] 82.3 96.2 90.9 101.7 123.1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 112.3 91.4 84.2 114 84.5 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 102.6 99.1 90.7 101.6 114 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 90.7 57.8 93.9 105 115.7 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 88 101.4 103.5 89.3 75.4 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 98.7 104.2 114.7 90.6 103.6 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   [list output truncated]
```

### `mvnram()`

Generates an `n x k` multivariate data matrix or a list of `n x k`
multivariate data matrices of length `R` from the multivariate normal
distribution. The model-implied matrices used to generate data is
derived from the Reticular Action Model (RAM) Matrices.

#### Single Random Data Set

Set matrices.

``` r
mu <- c(100, 100, 100)
A <- matrix(
  data = c(0, sqrt(0.26), 0, 0, 0, sqrt(0.26), 0, 0, 0),
  ncol = 3
)
S <- diag(c(225, 166.5, 116.5))
F <- I <- diag(3)
```

Run the function.

``` r
X <- mvnram(n = 100, mu = mu, A = A, S = S, F = F, I = I)
```

Explore the output.

``` r
str(X)
#>  num [1:100, 1:3] 86.6 92.5 99 104.6 87.5 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : NULL
pairs(X)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

``` r
colMeans(X)
#> [1] 98.55041 98.70115 99.52416
cov(X)
#>            [,1]      [,2]       [,3]
#> [1,] 225.128854  90.95131   7.837521
#> [2,]  90.951308 200.01940  66.685527
#> [3,]   7.837521  66.68553 142.354447
cor(X)
#>            [,1]      [,2]       [,3]
#> [1,] 1.00000000 0.4286051 0.04378018
#> [2,] 0.42860508 1.0000000 0.39519370
#> [3,] 0.04378018 0.3951937 1.00000000
```

#### Multiple Random Data Sets

Run the function.

``` r
Xstar <- mvnram(n = 100, mu = mu, A = A, S = S, F = F, I = I, R = 100)
```

Explore the output.

``` r
str(Xstar, list.len = 6)
#> List of 100
#>  $ : num [1:100, 1:3] 95.2 95.9 91.1 95.5 126.8 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 93.1 104.1 108.1 83.4 89 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 118.9 86.5 102.5 92.8 126.4 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 101.2 88.4 79.9 83.7 88.3 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 101.4 99.5 106.7 114.3 116.1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 104.1 65.8 88.1 78.1 81.4 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   [list output truncated]
```
