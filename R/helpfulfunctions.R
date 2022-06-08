#' Construct covariance matrix of time-dependent fGn model
#' 
#' This function constructs the covariance matrix of the fGn process under the 
#' assumption that the Hurst exponent changes over time. In the context of the 
#' \code{inla.ews} package this function is only used for simulating examples.
#' 
#' @param sigma Numeric of length 1. The innovation standard deviation
#' @param Hs Numeric of length corresponding to the length of the process. This gives 
#' the value of the Hurst exponent at each time point.
#' @return Returns the covariance matrix of the time dependent fGn process. This is 
#' a dense \code{matrix} object.
#' \code{object\$results}.
#' @examples 
#' n=100
#' Hs = seq(from=0.6,to=0.8,length.out=n)
#' sigma = 1
#' cov.matrix <- sigmaHmaker(sigma,Hs)
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords covariance matrix fgn
#' @export
sigmaHmaker = function(sigma,Hs){
  n=length(Hs)
  H2 = 2*Hs
  k=0:(n-1)
  sigmat = matrix(NA,n,n)
  for(i in 1:n){
    for(j in 1:n){
      t = min(i,j)
      k=abs(i-j)
      sigmat[i,j] = sigma^2/2*( abs(k-1)^H2[t]-2*abs(k)^H2[t]+abs(k+1)^H2[t] )
    }
  }
  return(sigmat)
}


