#' rgeneric model for AR(1) time-dependent process
#'
#' Defines the rgeneric model structure for the fractional Gaussian noise model 
#' with time-dependent correlation. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
#' This is intended for internal use only, but the documentation is included here in case someone want to change something.
#'
#' @param cmd Vector containing list of function names necessary for the rgeneric model.
#' @param theta Vector describing the hyperparameters in internal scaling.
#'
#'
#' @references{
#'   Jelena Ryvkina, Fractional Brownian Motion with Variable Hurst Parameter: Definition and Properties, Journal of Theoretical Probability, Volume 28, 2015, Pages 866-891, https://http://doi.org/10.1007/s10959-013-0502-3
#' }
#' @importFrom stats dnorm
rgeneric.ews.fgn = function(
    cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{
  
  

  envir = environment(sys.call()[[1]])
  
  
  ch = function(H){
    return(sqrt( H*(2*H-1)/(beta(2-2*H, H-0.5)) ))
  }
  
  cts = function(Ht,Hs){
    return(ch(Ht)*ch(Hs))
  }
  
  beta1 = function(Ht,Hs){
    return(beta(Ht-0.5, 2-Ht-Hs))
  }
  beta2 = function(Ht,Hs){
    return(beta(Ht-Hs+1, Ht+Hs-1))
  }
  hyperg2F1 = function(a,b,c,z, trunc=300){
    kk = 0:(trunc-1)
    cp = cumprod((a+kk)*(b+kk)/(c+kk)*z/(1:trunc)  ) 
    return( 1 + sum(cp))
  }
  
  R_H = function(t,s, Ht, Hs){
    return( cts(Ht,Hs)/(Ht+Hs) * ( beta1(Ht,Hs)*beta2(Ht,Hs)*s^(Ht+Hs) +
                                     beta1(Hs,Ht)*beta2(Hs,Ht)*t^(Ht+Hs) +
                                     beta1(Hs,Ht)*t*(t-s)^(Ht+Hs-1)/(Ht+Hs-1)*(t/s)^(Ht-Hs) *
                                     (hyperg2F1(1,2*Ht,Hs+Ht,(s-t)/s) -
                                        hyperg2F1(1,Ht-Hs,Hs+Ht,(s-t)/s)
                                     ) 
    )
    )
  }
  
  
  
  
  covmatmaker = function(HH, sx=1){
    
    nn = length(HH)
    cormat = matrix(NA,nrow=n,ncol=n)
    
    for(t in 1:nn){
      for(s in 1:nn){
        Ht = HH[t]
        Hs = HH[s]
        if (s<t){
          cormat[t,s] = R_H(t, s, Ht, Hs)
        }else{
          cormat[t,s] = R_H(s, t, Hs, Ht)
        }
        
      }
    }
    covmat = sx^2*cormat
    
    return(covmat)
  }
  
  
  
  interpret.theta = function() {
    if(!is.null(envir)){
      timee=get("time",envir)
      nn=get("n",envir)
    }
    kappa = exp(theta[1])
    bmax = 0.5/(timee[nn]-timee[1])
    bmin = -bmax
    b = bmin+(bmax-bmin)/(1+exp(-theta[2]))
    
    low = min(b*timee[1],b*timee[nn])
    high = max(b*timee[1],b*timee[nn])
    amin = 0.5-low
    amax = 1-high
    a = amin + (amax-amin)/(1+exp(-theta[3]))
    
    Hs = a+b*timee
    return(list(Hs = Hs, kappa = kappa, a=a,b=b,amin=amin,amax=amax,bmin=bmin,bmax=bmax))
  }
  mu = function() {
    
    return(numeric(0))
  }
  
  graph = function()
  {
    if(!is.null(envir)){
      nn=get("n",envir)
    }else{
      nn=get("n",environment())
    }
    return (matrix(1,nn,nn))
  }
  
  
  
  Q = function()  {
    if(!is.null(envir)){
      nn=get("n",envir)
    }
    
    hyperparam = interpret.theta()
    Hs = hyperparam$Hs
    kappa = hyperparam$kappa
    sx = 1/sqrt(hyperparam$kappa)
    a = hyperparam$a
    b = hyperparam$b
    
    covmat = covmatmaker(Hs, sx=sx)
    
    return (solve(covmat))
  }
  
  log.norm.const = function(){
    return(numeric(0))
  }
  
  log.prior = function()  {
    # if(!is.null(envir)){
    #   nn=get("n",envir)
    #   NN=get("N",envir)
    #   z=get("forcing",envir)
    # }
    params = interpret.theta()
    lprior = INLA::inla.pc.dprec(params$kappa, u=1, alpha=0.01, log=TRUE) + log(params$kappa) #kappa
    lprior = lprior + dnorm(theta[2],log=TRUE) #theta_a
    lprior = lprior + dnorm(theta[3],log=TRUE) #theta_b
    
    # if(length(z)>0){
    #   lprior = lprior + INLA::inla.pc.dprec(params$kappa_f, u=1, alpha=0.01, log=TRUE) + log(params$kappa_f) #kappa_f
    #   lprior = lprior + dnorm(theta[5],sd=5,log=TRUE) #F0
    # }
    
    return (lprior)
  }
  initial = function(){
    
    # if(length(z)>0){
    #   ini = c(0.,0.,0.,0.,0.)
    # }else{
    ini = c(0.,0.,0.)
    # }
    return (ini)
  }
  quit = function()  {
    return ()
  }
  if(is.null(theta)){
    theta = initial()
    
  }
  cmd = match.arg(cmd)
  val = do.call(cmd, args = list())
  return (val)
}

