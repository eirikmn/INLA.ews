
<!-- README.md is generated from README.Rmd. Please edit that file -->

# INLA.ews

<!-- badges: start -->

[![R-CMD-check](https://github.com/eirikmn/INLA.ews/workflows/R-CMD-check/badge.svg)](https://github.com/eirikmn/INLA.ews/actions)
<!-- badges: end -->

This Repository contains the INLA.ews package for Bayesian detection of
early warning signals.

## Installation

You can install the development version of INLA.ews from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("eirikmn/INLA.ews")
```

## Example

This is a basic example which shows you how to use the package to
perform Bayesian analysis of a simulated time series which exhibit early
warning signals in the form of increasing correlation. The ‘true’
development of the autocorrelation coefficient is included in the
resulting plot.

``` r
library(INLA.ews)
n = 1000
sigma = 1
a=0.2
b=0.7/n
time = 1:n
phis = a+b*time
data=numeric(n)
data[1] = rnorm(1,mean=0,sd=sigma)
for(i in 2:n){
  data[i] = rnorm(1, mean=phis[i]*data[i-1],sd=sigma)
}
object = inla.ews(data,model="ar1", memory.true=phis)
plot(object)
```

<embed src="man/figures/README-plot-1.pdf" width="100%" type="application/pdf" />

## Attribution

This code is associated and written for an upcoming paper. Feel free to
use the code, but please cite the accompanying paper (when they are
published).

## License

The code in this repository is made available under the terms of the MIT
License. For details, see LICENSE file.
