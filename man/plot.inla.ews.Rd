% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.inla.ews.R
\name{plot.inla.ews}
\alias{plot.inla.ews}
\title{Plot inla.ews model}
\usage{
\method{plot}{inla.ews}(
  x,
  plot.options = list(plot.hyper = TRUE, plot.forced = TRUE, plot.memory = TRUE,
    use.median = TRUE, memory.true = TRUE),
  postscript = FALSE,
  pdf = FALSE,
  prefix = "INLA.ews.plots/figure-",
  ...
)
}
\arguments{
\item{x}{\code{inla.ews} S3 class. Output of \code{\link{inla.ews}} function}

\item{plot.options}{list with settings for plot.}

\item{postscript}{Boolean variable indicating if postscript files should be produced instead.}

\item{pdf}{Boolean variable indicating if pdf files should be produced instead.}

\item{prefix}{The prefix for created files. Additional numbering is added.}

\item{...}{Additional arguments to \code{postscripts()}, \code{pdf()} or \code{dev.new()}.}
}
\description{
Plots results from inla.ews S3 class.
}
\examples{
\donttest{
n = 200
sigma = 1
a=0.2
b=0.7/n
time = 1:n
phis = a+b*time
data=ar1_timedep_sim(n,phis=phis)

object = inla.ews(data,model="ar1", print.progress=TRUE,
                    memory.true=phis)
plot(object)

}

}
\seealso{
\code{\link{inla.ews}}
}
\author{
Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
}
\keyword{inla.ews}
\keyword{plot}
