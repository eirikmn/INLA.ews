% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgeneric.ews.fgn2.R
\name{rgeneric.ews.fgn.forcing2}
\alias{rgeneric.ews.fgn.forcing2}
\title{rgeneric model for fgn time-dependent process}
\usage{
rgeneric.ews.fgn.forcing2(
  cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),
  theta = NULL
)
}
\arguments{
\item{cmd}{Vector containing list of function names necessary for the rgeneric model.}

\item{theta}{Vector describing the hyperparameters in internal scaling.}
}
\description{
Defines the rgeneric model structure for the fractional Gaussian noise model 
with time-dependent correlation. This model includes forcing. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
This is intended for internal use only, but the documentation is included here in case someone want to change something.
}
