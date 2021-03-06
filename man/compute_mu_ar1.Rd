% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{compute_mu_ar1}
\alias{compute_mu_ar1}
\title{Compute mu in C for AR(1) model}
\usage{
compute_mu_ar1(mu, n, tforcing, phis)
}
\arguments{
\item{mu}{Vector that will be overwritten with the forced response.}

\item{n}{Integer denoting the dimension of \code{mu}.}

\item{tforcing}{Vector with the shifted and scaled known forcing.}

\item{phis}{Vector describing the time-dependent lag-one correlation.}
}
\description{
Void function that gives the forced response for the AR(1) model with 
time-dependent memory in C.
}
\seealso{
\code{\link{mu.computer}}
}
\author{
Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
}
\keyword{C}
\keyword{ar1}
\keyword{forcing}
