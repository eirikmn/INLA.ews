% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpfulfunctions.R
\name{set.options}
\alias{set.options}
\title{Import default arguments}
\usage{
set.options(opt, default.opt)
}
\arguments{
\item{opt}{List object with different specifications.}

\item{default.opt}{List of default variables corresponding to \code{opt}.}
}
\value{
Returns the \code{opt} list, but with values from \code{default.opt} inserted
in missing values.
}
\description{
Fills out missing arguments in list \code{opt} with default arguments in list
\code{default.opt}.
}
\seealso{
\code{\link{inla.ews}}
}
\author{
Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
}
\keyword{default}
\keyword{inla.ews}
