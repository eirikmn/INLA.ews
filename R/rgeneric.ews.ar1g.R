#' rgeneric model for AR(1) time-dependent process using Green's function
#'
#' Defines the rgeneric model structure for the fractional Gaussian noise model 
#' with time-dependent correlation. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
#' This is intended for internal use only, but the documentation is included here in case someone want to change something.
#'
#' @param cmd Vector containing list of function names necessary for the rgeneric model.
#' @param theta Vector describing the hyperparameters in internal scaling.
#'
#'
#' @importFrom stats dnorm
rgeneric.ews.ar1g = function(
    cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{
  
  
  # require("INLA.climate2",quietly=TRUE)
  
  tau = exp(15)
  envir = environment(sys.call()[[1]])
  
  interpret.theta = function() {
    if(!is.null(envir)){
      timme=get("time",envir)
      nn=get("n",envir)
    }
    timee = timme-timme[1]
    timee = timee/timee[nn]
    kappa = exp(theta[1])
    bmax = 0.5/(timee[nn]-timee[1])
    bmin = -bmax
    b = bmin+(bmax-bmin)/(1+exp(-theta[2]))
    
    low = min(b*timee[1],b*timee[nn])
    high = max(b*timee[1],b*timee[nn])
    amin = 0.5-low
    amax = 1-high
    a = amin + (amax-amin)/(1+exp(-theta[3]))
    
    phis = a+b*timee
    return(list(phis = phis, kappa = kappa, a=a,b=b,amin=amin,amax=amax,bmin=bmin,bmax=bmax))
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
  
  greensar1 = function(t,s,a,b,n){
    if(t-s<0){
      return(0)
    }else{
      # phi = a+b*max(t,s)/n
      # phi = a+b*(t+s)/2/n
      phi = a+b*min(t,s)/n
      # H = a+b*(t+s)/2/n
      ret = phi^(t-s)
      return(ret)
    }
  }
  
  
  Q = function()  {
    if(!is.null(envir)){
      nn=get("n",envir)
    }
    
    hyperparam = interpret.theta()
    sx = 1/sqrt(hyperparam$kappa)
    a = hyperparam$a
    b=hyperparam$b
    phis = hyperparam$phis
    lambdas = -log(phis)
    
    kappas = numeric(nn)
    
    kappas[1] = 2*lambdas[1]/sx^2
    kappas[2:nn] = 1/sx^2*1/( 1/(2*lambdas[2:nn])-phis[2:nn]^2/(2*lambdas[1:(nn-1)]) )
    
    ii = c(1:nn,2:nn) ; jj = c(1:nn,1:(nn-1))
    xx = c(kappas[1:(nn-1)]+phis[2:nn]^2*kappas[2:nn],kappas[nn],
           -phis[2:nn]*kappas[2:nn])
    Q = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
    return(Q)
    # Hs = hyperparam$phis
    # kappa = hyperparam$kappa
    # 
    # 
    # H2 = 2*Hs
    # k=0:(nn-1)
    # # sigmat = matrix(NA,nn,nn)
    # # for(i in 1:nn){
    # #   for(j in 1:nn){
    # #     #t = max(i,j)
    # #     t = (i+j)/2/nn
    # #     H2 = 2*(hyperparam$a+hyperparam$b*t)
    # #     k=abs(i-j)
    # #     sigmat[i,j] = sx^2/2*( abs(k-1)^H2-2*abs(k)^H2+abs(k+1)^H2 )
    # #   }
    # # }
    # # #sigmat = sigmamaker(nn,sx,Hs)
    # 
    # Gmat = matrix(NA,nn,nn)
    # for(i in 1:nn){
    #   for(j in 1:nn){
    #     Gmat[i,j] = greensar1(i,j,a,b,nn)
    #   }
    # }
    # covmat = sx^2*Gmat%*%t(Gmat)
    
    # return (solve(covmat))
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