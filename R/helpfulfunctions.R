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

#' Simulate forced response
#' 
#' This function employs Monte Carlo simulations to estimate the forced response.
#' 
#' @param object S3 object of type \code{inla.ews} which includes result from 
#' \code{inla}-program.
#' @param quick Boolean indicating whether or not the function should run without
#' storing simulations (faster, and less memory). If quantiles are desired this must
#' be set to \code{FALSE}.
#' @param seed Seed used for the random number generator.
#' @param print.progress boolean indicating if progress should be printed to screen.
#' @return Returns the same S3 object of class \code{inla.ews} as included in the 
#' input arguments, but appends summary statistics \code{object\$forced}.
#' 
#' @examples 
#' n = 200
#' time=1:n
#' a = 0.6
#' b = 0.35/n
#' Hs = a+b*time
#' F0 = -3
#' sigmaf=0.1
#' sigma = 1.2
#' #Hs = rep(0.75,n)
#' library(INLA.ews)
#' sigmat = sigmaHmaker(sigma,Hs)
#' sigmachol = chol(sigmat)
#' noise = sigmachol%*%rnorm(n)
#' forcing = arima.sim(model=list(ar=c(0.9)),n)+1:n/n*10
#' z = sigmaf*(forcing+F0)
#' 
#' struct = (1:n-0.5)^(Hs-3/2)
#' muvek=numeric(n)
#' for(i in 1:n){
#'   muvek[i] = rev(struct[1:i])%*%z[1:i]
#' }
#' 
#' y=muvek+noise + 200
#' plot(y,type="l",col="grey",lwd=1.1)
#' 
#' object = INLA.ews::inla.ews(y,forcing,model="fgn",memory.true=Hs)
#' object = compute.mu(object,quick=FALSE,print.progress=TRUE)
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords forcing response mean 
#' @export
#' @importFrom stats rnorm
forcingmaker <- function(object,quick=FALSE,seed=1234,
                       print.progress=FALSE){
  time.start = Sys.time()
  model=object$.args$model
  set.seed(seed)
  if(length(object$.args$forcing)==0){
    return(object)
  }
  nsamples = ncol(object$results$simulations$H_sims)
  if(print.progress){
    cat("Starting mu Monte Carlo sampling with n = ",format(nsamples,scientific=F)," simulations..\n",sep="")
  }
  #We sample intercepts. We know this is Gaussian since we use Gaussian priors
  interceptsims = stats::rnorm(nsamples,object$results$fixed$mean,object$results$fixed$sd) 
                               
  n=nrow(object$.args$inladata)
  hypersamples = INLA::inla.hyperpar.sample(nsamples,object$inlafit)
  
  sigmaf_samples = 1/sqrt(exp(hypersamples[,4]))
  F0_samples = hypersamples[,5]
  
  if(model %in% c("fgn","lrd")){
    memorysamples = object$results$simulations$H_sims
  }else if(model %in% c("ar1","ar(1)","1")){
    memorysamples = object$results$simulations$phi_sims
  }
  
  forcing = object$.args$forcing
  if(!quick){
    meanmat = matrix(NA,n,nsamples)
  }
  xsumvec = 0
  x2sumvec = 0
  if(model %in% c("ar1","ar(1)","1")){
    for(iter in 1:nsamples){
      if(print.progress){
        if(iter %% 5000 == 0) cat("Sampling mean vector ",iter,"/",nsamples,"\n",sep="")
      }
      lambdas = memorysamples[,iter]-1
      zz = sigmaf_samples[iter]*(F0_samples[iter]+forcing)
      struct = exp(lambdas*(1:n)-0.5)
      muvek=numeric(n)
      for(i in 1:n){
        muvek[i] = rev(struct[1:i])%*%zz[1:i]
      }
      if(!quick){
        meanmat[,iter]=muvek
      }
      xsumvec = xsumvec + muvek
      x2sumvec = x2sumvec + muvek^2
      
    }
    
  }else if(model %in% c("fgn","lrd")){
    for(iter in 1:nsamples){
      if(print.progress){
        if(iter %% 5000 == 0) cat("Sampling mean vector #",iter,"/",nsamples,"\n",sep="")
      }
      zz = sigmaf_samples[iter]*(F0_samples[iter]+forcing)
      struct = (1:n-0.5)^(memorysamples[iter]-3/2)
      muvek=numeric(n)
      for(i in 1:n){
        muvek[i] = rev(struct[1:i])%*%zz[1:i]
      }
      muvek = interceptsims[iter]+muvek
      if(!quick){
        meanmat[,iter]=muvek
      }
      xsumvec = xsumvec + muvek
      x2sumvec = x2sumvec + muvek^2
    }
    
  }
  if(print.progress){
    cat("Completed sampling!\n",sep="")
  }
  
  mu.mean = as.numeric(xsumvec/nsamples)
  mu.sd=as.numeric(sqrt( 1/(nsamples-1)*( 
    x2sumvec -2*mu.mean*xsumvec + nsamples*mu.mean^2 ) ))
  
  object$forced = list(mean = mu.mean, sd = mu.sd)
  
  if(!quick){
    if(print.progress){
      cat("Computing quantiles","...",sep="")
    }
    mu.quant0.025 = numeric(n)
    mu.quant0.25 = numeric(n)
    mu.quant0.5 = numeric(n)
    mu.quant0.75 = numeric(n)
    mu.quant0.975 = numeric(n)
    for(iter in 1:n){
      dens = density(meanmat[iter,])
      mu.quant0.025[iter]=INLA::inla.qmarginal(0.025,dens)
      mu.quant0.25[iter]=INLA::inla.qmarginal(0.25,dens)
      mu.quant0.5[iter]=INLA::inla.qmarginal(0.5,dens)
      mu.quant0.75[iter]=INLA::inla.qmarginal(0.75,dens)
      mu.quant0.975[iter]=INLA::inla.qmarginal(0.975,dens)
    }
    object$forced$quant0.025 = mu.quant0.025
    object$forced$quant0.25 = mu.quant0.25
    object$forced$quant0.5 = mu.quant0.5
    object$forced$quant0.75 = mu.quant0.75
    object$forced$quant0.975 = mu.quant0.975
    if(print.progress){
      cat(" completed!\n",sep="")
    }
  }
  
  time.total = difftime(Sys.time(), time.start,units="secs")[[1]]
 
  object$cpu.used$meansim = time.total
  
# plot(object$.args$data,ylim=range(object$forced$quant0.025,object$forced$quant0.975),type="l")
# lines(object$forced$mean,col="blue")
# lines(object$forced$quant0.025,col="red")
# lines(object$forced$quant0.975,col="red")
# 
# mean(object$forced$mean)-mean(object$.args$data)

#  
return(object)
  
}

