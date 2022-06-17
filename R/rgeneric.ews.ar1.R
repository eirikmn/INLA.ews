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
rgeneric.ews.ar1 = function(
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

    phis = a+b*timee

    kappa_f = exp(theta[4])
    F0 = theta[5]

    return(list(phis = phis, kappa_eps = kappa_eps, a=a,b=b,kappa_f=kappa_f,F0=F0))
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
    phis = params$a+params$b*timee
    kappa_eps = params$kappa_eps
    kappa1=kappa_eps
    ii=c(1,nn,2:(nn-1),1:(nn-1));jj=c(1,nn,2:(nn-1),2:nn)
    xx = kappa_eps*c(1+phis[2]^2,1,1+phis[3:nn]^2,-phis[2:nn])
    Q = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
    return (Q)
  }
  log.norm.const = function(){return(numeric(0))}
  log.prior = function(){
    params = interpret.theta()
    lprior = INLA::inla.pc.dprec(params$kappa_eps, u=1, alpha=0.01, log=TRUE) + log(params$kappa_eps) #kappa
    lprior = lprior + dnorm(theta[2],log=TRUE) #theta_b
    lprior = lprior + dnorm(theta[3],log=TRUE) #theta_a

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
rgeneric.ews.ar1.forcing = function(
    cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{
  require(INLA.ews,quietly=TRUE)
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
    # a = min(-1-b*timee) + (max(1-b*timee)-min(-1-b*timee))/(1+exp(-theta[3]))
    
    phis = a+b*timee
    
    kappa_f = exp(theta[4])
    F0 = theta[5]
    
    return(list(phis = phis, kappa_eps = kappa_eps, a=a,b=b,kappa_f=kappa_f,F0=F0))
  }
  
  
  mu = function() {
    if(!is.null(envir)){
      nn=get("n",envir)
      z = get("forcing",envir)
    }else{
      nn=get("n",environment())
      z = get("forcing",environment())
    }
    params=interpret.theta()
    zz = 1/sqrt(params$kappa_f)*(params$F0+z)
    
    muvek = numeric(nn)
    compute_mu_ar1(muvek, nn, zz,  params$phis)
    
    
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
    phis = params$a+params$b*timee
    kappa_eps = params$kappa_eps
    kappa1=kappa_eps
    ii=c(1,nn,2:(nn-1),1:(nn-1));jj=c(1,nn,2:(nn-1),2:nn)
    xx = kappa_eps*c(1+phis[2]^2,1,1+phis[3:nn]^2,-phis[2:nn])
    Q = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
    return (Q)
  }
  log.norm.const = function(){return(numeric(0))}
  log.prior = function(){
    params = interpret.theta()
    lprior = INLA::inla.pc.dprec(params$kappa_eps, u=1, alpha=0.01, log=TRUE) + log(params$kappa_eps) #kappa
    lprior = lprior + dnorm(theta[2],log=TRUE) #theta_b
    lprior = lprior + dnorm(theta[3],log=TRUE) #theta_a
    lprior = lprior + INLA::inla.pc.dprec(params$kappa_f, u=1, alpha=0.01, log=TRUE) + log(params$kappa_f) #kappa
    lprior = lprior + dnorm(theta[5],log=TRUE) #F0
    
    return (lprior)
  }
  initial = function(){return (c(0.,0.,0.,0.,0.)) }
  quit = function(){return ()  }
  if(is.null(theta)){
    theta = initial()
  }
  cmd = match.arg(cmd)
  val = do.call(cmd, args = list())
  return (val)
}

