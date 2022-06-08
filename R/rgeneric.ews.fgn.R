#' rgeneric model for AR(1) time-dependent process
#'
#' Defines the rgeneric model structure for the fractional Gaussian noise model 
#' with time-dependent correlation. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
#' This is intended for internal use only, but the documentation is included here in case someone want to change something.
#'  See example below for how rgeneric is used can be used.
#'
#' @param cmd Vector containing list of function names necessary for the rgeneric model.
#' @param theta Vector describing the hyperparameters in internal scaling.
#'
#'
#' @importFrom INLA inla.pc.dprec
#' @importFrom stats dnorm
rgeneric.ews.fgn = function(
    cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{

  # require("INLA.climate2",quietly=TRUE)

  tau = exp(15)
  envir = environment(sys.call()[[1]])

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
    if(!is.null(envir)){
      nn=get("n",envir)
      
    }else{
      nn=get("n",environment())
      
    }
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
    H2 = 2*Hs
    k=0:(nn-1)
    sigmat = matrix(NA,nn,nn)
    for(i in 1:nn){
      for(j in 1:nn){
        t = min(i,j)
        k=abs(i-j)
        sigmat[i,j] = sx^2/2*( abs(k-1)^H2[t]-2*abs(k)^H2[t]+abs(k+1)^H2[t] )
      }
    }
    #sigmat = sigmamaker(nn,sx,Hs)

    return (solve(sigmat))
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



#' rgeneric model for fgn time-dependent process
#'
#' Defines the rgeneric model structure for the fractional Gaussian noise model 
#' with time-dependent correlation. This model includes forcing. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
#' This is intended for internal use only, but the documentation is included here in case someone want to change something.
#'  See example below for how rgeneric is used can be used.
#'
#' @param cmd Vector containing list of function names necessary for the rgeneric model.
#' @param theta Vector describing the hyperparameters in internal scaling.
#'
#'
#' @importFrom INLA inla.pc.dprec
#' @importFrom stats dnorm
rgeneric.ews.fgn.forcing = function(
    cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{

  

  tau = exp(15)
  envir = environment(sys.call()[[1]])

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
    kappa_f = exp(theta[4])
    F0 = theta[5]
    return(list(Hs = Hs, kappa = kappa, a=a,b=b,amin=amin,amax=amax,bmin=bmin,bmax=bmax,kappa_f=kappa_f,F0=F0))
  }
  mu = function() {
    if(!is.null(envir)){
      nn=get("n",envir)
      z = get("forcing",envir)
    }else{
      nn=get("n",environment())
      z = get("forcing",environment())
    }

    params = interpret.theta()
    #cat("F0: ",params$F0," sigmaf: ",1/sqrt(params$kappa_f),"\n")
    zz = 1/sqrt(params$kappa_f)*(params$F0+z)
    struct = (1:nn-0.5)^(params$Hs-3/2)
    muvek=numeric(nn)
    for(i in 1:nn){
      muvek[i] = rev(struct[1:i])%*%zz[1:i]
    }


    return(muvek)
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
    H2 = 2*Hs
    k=0:(nn-1)
    sigmat = matrix(NA,nn,nn)
    for(i in 1:nn){
      for(j in 1:nn){
        t = min(i,j)
        k=abs(i-j)
        sigmat[i,j] = sx^2/2*( abs(k-1)^H2[t]-2*abs(k)^H2[t]+abs(k+1)^H2[t] )
      }
    }
    #sigmat = sigmamaker(nn,sx,Hs)

    return (solve(sigmat))
  }

  log.norm.const = function(){
    return(numeric(0))
  }

  log.prior = function()  {
    params = interpret.theta()
    lprior = INLA::inla.pc.dprec(params$kappa, u=1, alpha=0.01, log=TRUE) + log(params$kappa) #kappa_x
    lprior = lprior + dnorm(theta[2],log=TRUE) #theta_a
    lprior = lprior + dnorm(theta[3],log=TRUE) #theta_b
    lprior = lprior + INLA::inla.pc.dprec(params$kappa_f, u=1, alpha=0.01, log=TRUE) + log(params$kappa_f) #kappa_f
    lprior = lprior + dnorm(theta[5],sd=5,log=TRUE) #F0
    return (lprior)
  }
  initial = function(){
    ini = c(0.,0.,0.,0.,0.)
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

