#' rgeneric model for AR(1) time-dependent process
#'
#' Defines the rgeneric model structure for the AR(1) model with time-dependent
#' correlation. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
#' This is intended for internal use only, but the documentation is included here in case someone want to change something.
#'
#' @param cmd Vector containing list of function names necessary for the rgeneric model.
#' @param theta Vector describing the hyperparameters in internal scaling.
#'
#' @importFrom stats dnorm
rgeneric.ar1 = function(
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
    kappa_eps = exp(theta[1])
    r = diff(range(timee))
    b = -1/r+2/r*1/(1+exp(-theta[2]))
    
    low = min(b*timee[1],b*timee[nn])
    high = max(b*timee[1],b*timee[nn])
    
    a = -low + (1-high+low)/(1+exp(-theta[3]))
    
    lambdas = -log(a+b*timee)
    kappa2s = kappa_eps*2*lambdas
    
    cc = 1/(nn-1)
    phis = c(exp(-lambdas*c(1,diff(timee)/cc)) ) #rescale
    
    #kappa_f = exp(theta[4])
    #F0 = theta[5]
    
    return(list(phis = phis, kappa_eps = kappa_eps, a=a,b=b,#kappa_f=kappa_f,F0=F0
                lambdas=lambdas, kappa2s=kappa2s, sigma2s=1/kappa2s))
  }
  
  mu = function() {
    return(numeric(0))
    
  }
  
  graph = function() {
    if(!is.null(envir)){
      nn=get("n",envir)
    }else{
      nn=get("n",environment())
    }
    ii = c(1,nn,2:(nn-1),1:(nn-1));jj=c(1,nn,2:(nn-1),2:nn)
    xx=c(1,1,rep(1,nn-2),rep(1,nn-1))
    G = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
    G[G != 0] = 1
    return (G)
  }
  Q = function(){
    if(!is.null(envir)){
      nn=get("n",envir)
      timee=get("time",envir)
    }
    params = interpret.theta()
    phis = params$phis
    kappa_eps = params$kappa_eps
    kappa1=kappa_eps
    kappa2s = params$kappa2s
    ii=c(1,nn,2:(nn-1),1:(nn-1))
    jj=c(1,nn,2:(nn-1),2:nn)
    #xx = kappa_eps*c(1+phis[2]^2,1,1+phis[3:nn]^2,-phis[2:nn])
    xx = c(kappa2s[1]*(1-phis[1]^2)+kappa2s[2]*phis[2]^2, kappa2s[nn],
           kappa2s[2:(nn-1)]+kappa2s[3:nn]*phis[3:nn]^2,
           -phis[2:nn]*kappa2s[2:nn])
    
    Q = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
    return (Q)
  }
  log.norm.const = function(){return(numeric(0))}
  log.prior = function(){
    params = interpret.theta()
    lprior = INLA::inla.pc.dprec(params$kappa_eps, u=1, alpha=0.01, log=TRUE) + log(params$kappa_eps) #kappa
    #lprior = lprior + dnorm(theta[2],log=TRUE) #theta_b
    b=params$b; ra = 0; rb=1
    lprior = -log(rb-ra)-theta[2]-log(1+exp(-theta[2]))
    #lprior = lprior -log(rb-ra) + log(b-ra)+log(rb-b) -log(rb-ra)
    lprior = lprior + dnorm(theta[3],sd=3,log=TRUE) #theta_a
    
    return (lprior)
  }
  initial = function(){
    return (c(0.,0.,0.)) 
  }
  
  quit = function(){return ()  }
  if(is.null(theta)){
    theta = initial()
  }
  cmd = match.arg(cmd)
  val = do.call(cmd, args = list())
  return (val)
}






#' rgeneric forcing model for AR(1) time-dependent process
#'
#' Defines the rgeneric model structure for the AR(1) model with time-dependent
#' correlation. This model includes forcing. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
#' This is intended for internal use only, but the documentation is included here in case someone want to change something.
#'
#' @param cmd Vector containing list of function names necessary for the rgeneric model.
#' @param theta Vector describing the hyperparameters in internal scaling.
#'
#'
#' @importFrom stats dnorm
rgeneric.ar1.forcing = function(
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
    kappa_eps = exp(theta[1])
    r = diff(range(timee))
    b = -1/r+2/r*1/(1+exp(-theta[2]))
    
    low = min(b*timee[1],b*timee[nn])
    high = max(b*timee[1],b*timee[nn])
    
    a = -low + (1-high+low)/(1+exp(-theta[3]))
    
    lambdas = -log(a+b*timee)
    kappa2s = kappa_eps*2*lambdas
    
    cc = 1/(nn-1)
    phis = c(exp(-lambdas*c(1,diff(timee)/cc)) ) #rescale
    
    kappa_f = exp(theta[4])
    kappa2fs = kappa_f*2*lambdas
    F0 = theta[5]
    
    return(list(phis = phis, kappa_eps = kappa_eps, a=a,b=b,#kappa_f=kappa_f,F0=F0
                lambdas=lambdas, kappa2s=kappa2s, sigma2s=1/kappa2s,
                kappa_f=kappa_f, kappa2fs=kappa2fs, F0=F0))
  }
  
  mu = function() {
    if(!is.null(envir)){
      nn=get("n",envir)
      timee=get("time",envir)
      fforcing=get("forcing",envir)
    }
    hyperparam = interpret.theta()
    #print(hyperparam)
    
    sf = 1/sqrt(hyperparam$kappa_f)
    kappa2fs = hyperparam$kappa2fs
    F0 = hyperparam$F0
    
    zz = (fforcing+F0)/sqrt(kappa2fs)
    
    innerstruct = timee
    #solution 1:
    llambdas = hyperparam$lambdas
    
    #solution 2:
    #innerstruct = 0.5+seq(0,nn-1,length.out=nn)
    pp = hyperparam$phis
    #llambdas = pp-1
    
    
    ##
    
    struktur = exp(-hyperparam$lambdas*innerstruct)
    
    muvek = numeric(nn)
    
    for(i in 1:nn){
      muvek[i] = rev(struktur[1:i])%*%zz[1:i]
    }
    
    return(muvek)
    
  }
  
  graph = function() {
    if(!is.null(envir)){
      nn=get("n",envir)
    }else{
      nn=get("n",environment())
    }
    ii = c(1,nn,2:(nn-1),1:(nn-1));jj=c(1,nn,2:(nn-1),2:nn)
    xx=c(1,1,rep(1,nn-2),rep(1,nn-1))
    G = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
    G[G != 0] = 1
    return (G)
  }
  Q = function(){
    if(!is.null(envir)){
      nn=get("n",envir)
      timee=get("time",envir)
    }
    params = interpret.theta()
    phis = params$phis
    kappa_eps = params$kappa_eps
    kappa1=kappa_eps
    kappa2s = params$kappa2s
    ii=c(1,nn,2:(nn-1),1:(nn-1))
    jj=c(1,nn,2:(nn-1),2:nn)
    #xx = kappa_eps*c(1+phis[2]^2,1,1+phis[3:nn]^2,-phis[2:nn])
    xx = c(kappa2s[1]*(1-phis[1]^2)+kappa2s[2]*phis[2]^2, kappa2s[nn],
           kappa2s[2:(nn-1)]+kappa2s[3:nn]*phis[3:nn]^2,
           -phis[2:nn]*kappa2s[2:nn])
    
    Q = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
    return (Q)
  }
  log.norm.const = function(){return(numeric(0))}
  log.prior = function(){
    params = interpret.theta()
    lprior = INLA::inla.pc.dprec(params$kappa_eps, u=1, alpha=0.01, log=TRUE) + log(params$kappa_eps) #kappa
    #lprior = lprior + dnorm(theta[2],log=TRUE) #theta_b
    b=params$b; ra = 0; rb=1
    lprior = -log(rb-ra)-theta[2]-log(1+exp(-theta[2]))
    #lprior = lprior -log(rb-ra) + log(b-ra)+log(rb-b) -log(rb-ra)
    lprior = lprior + dnorm(theta[3],sd=3,log=TRUE) #theta_a
    lprior = lprior + INLA::inla.pc.dprec(params$kappa_f, u=1, alpha=0.01, log=TRUE) + log(params$kappa_f) #kappa_f
    lprior = lprior + dnorm(theta[5],sd=3,log=TRUE) #F0
    
    return (lprior)
  }
  initial = function(){
    return (c(0.,0.,0.,0.,0.)) 
  }
  
  quit = function(){return ()  }
  if(is.null(theta)){
    theta = initial()
  }
  cmd = match.arg(cmd)
  val = do.call(cmd, args = list())
  return (val)
}

