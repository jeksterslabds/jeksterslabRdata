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
#>  num [1:100] 87 86.8 97 86.3 108.5 ...
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
#>  $ : num [1:100] 96.3 113.2 109 91.7 99.7 ...
#>  $ : num [1:100] 101 104 117 100 112 ...
#>  $ : num [1:100] 99.6 107.4 129.5 101.2 114.4 ...
#>  $ : num [1:100] 100.5 104.5 105.1 101.1 98.5 ...
#>  $ : num [1:100] 115.6 111.7 115.9 91.5 93.2 ...
#>  $ : num [1:100] 74 105.2 94.6 121.6 118.3 ...
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
#>  num [1:100, 1:3] 118.6 108.9 101.1 96.4 100.7 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : NULL
pairs(X)
```

<img src="man/figures/README-unnamed-chunk-9-1.png" width="100%" />

``` r
colMeans(X)
#> [1] 101.60171 102.77214  99.25424
cov(X)
#>           [,1]     [,2]      [,3]
#> [1,] 234.07533 115.2056  85.62259
#> [2,] 115.20561 235.3096 143.49818
#> [3,]  85.62259 143.4982 258.58048
cor(X)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.4908807 0.3480268
#> [2,] 0.4908807 1.0000000 0.5817397
#> [3,] 0.3480268 0.5817397 1.0000000
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
#>  $ : num [1:100, 1:3] 108.7 105.4 82.8 88.1 104 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 111.2 94.5 98.3 91.5 58.6 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 109.8 97.5 85.1 75 81.3 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 79.5 99.7 128.8 115.4 98.3 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 78 111.3 120.7 90.4 84.3 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 105 93.2 112 95.2 95.1 ...
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
#>  num [1:100, 1:3] 79.4 96.6 97 91.1 113.1 ...
#>  - attr(*, "dimnames")=List of 2
#>   ..$ : NULL
#>   ..$ : NULL
pairs(X)
```

<img src="man/figures/README-unnamed-chunk-14-1.png" width="100%" />

``` r
colMeans(X)
#> [1] 101.03277  98.97819  98.94974
cov(X)
#>           [,1]     [,2]      [,3]
#> [1,] 220.99760 146.1968  73.66028
#> [2,] 146.19678 258.3556 132.28990
#> [3,]  73.66028 132.2899 191.79115
cor(X)
#>           [,1]      [,2]      [,3]
#> [1,] 1.0000000 0.6118360 0.3577876
#> [2,] 0.6118360 1.0000000 0.5942969
#> [3,] 0.3577876 0.5942969 1.0000000
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
#>  $ : num [1:100, 1:3] 109.4 121.2 126.9 86.8 95.9 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 89 90.1 79.3 114.1 107.8 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 78 99 116 124 128 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 83.6 103.7 120.9 80 119.1 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 95.1 82.4 109.9 87.6 95.7 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>  $ : num [1:100, 1:3] 85.8 126 89 119.9 72.3 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : NULL
#>   [list output truncated]
```
