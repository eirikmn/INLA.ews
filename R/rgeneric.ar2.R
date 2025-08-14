#' rgeneric model for AR(2) time-dependent process
#'
#' Defines the rgeneric model structure for the AR(2)/nested AR(1) model with time-dependent
#' correlation. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
#' This is intended for internal use only, but the documentation is included here in case someone want to change something.
#'
#' @param cmd Vector containing list of function names necessary for the rgeneric model.
#' @param theta Vector describing the hyperparameters in internal scaling.
#' @importFrom stats dnorm dgamma
rgeneric.ar2 = function(
    cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{
  
  
  envir = environment(sys.call()[[1]])
  
  
  interpret.theta = function() {
    if(!is.null(envir)){
      timee=get("time",envir)
      nn=get("n",envir)
    }
    tau = exp(15) 
    r = diff(range(timee))
    b_Y = -1/r+2/r*1/(1+exp(-theta[1]))
    low_Y = min(b_Y*timee[1],b_Y*timee[nn])
    high_Y = max(b_Y*timee[1],b_Y*timee[nn])
    a_Y = -low_Y + (1-high_Y+low_Y)/(1+exp(-theta[2]))
    lambdas_Y = -log(a_Y+b_Y*timee)
    #cc_Y = 1/(nn-1)
    cc_Y = max(timee)/(nn-1)
    #cc_Y = 1
    phis = c(exp(-lambdas_Y*c(1,diff(timee)/cc_Y)) ) 
    
    kappa_eps_V = exp(theta[3])
    kappa2_eps_V=kappa_eps_V*2*lambdas_Y
    
    
    
    
    C=1
    b_V = -1/r+2/r*1/(1+exp(-theta[4]))
    low_V = min(b_V*timee[1],b_V*timee[nn])
    high_V = max(b_V*timee[1],b_V*timee[nn])
    a_V = -low_V + (1-high_V+low_V)/(1+exp(-theta[5]))
    lambdas = -log(a_V+b_V*timee)
    rhos=a_V+b_V*timee
    #cc = 1/(nn-1)
    cc = max(timee)/(nn-1)
    #cc = 1
    rhos = c(exp(-lambdas*c(1,diff(timee)/cc)) ) 
    
    
    print(theta)
    return(list(phis = phis,a_Y=a_Y,b_Y=b_Y,
                tau=tau,
                C=C,
                rhos=rhos,a_V=a_V,b_V=b_V,
                kappa_eps_V = kappa_eps_V, 
                kappa2_eps_V=kappa2_eps_V))
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
    nn=nn*2
    
    
    
    ii = c(1,(nn/2),((nn/2)+1),nn,
           2:(nn/2-1), (nn/2+2):(nn-1),
           1:((nn/2)-1),((nn/2)+1):(nn-1),
           1:((nn/2)), 
           1:((nn/2)-1))
    
    jj = c(1,(nn/2),((nn/2)+1),nn,
           2:(nn/2-1), (nn/2+2):(nn-1),
           2:(nn/2),((nn/2)+2):nn, 
           ((nn/2)+1):(nn),  
           ((nn/2)+2):(nn))
    
    
    
    xx=c(rep(1,length(ii)))
    G = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
    G[G != 0] = 1
    return (G)
  }
  
  
  
  Q = function(){
    if(!is.null(envir)){
      nn=get("n",envir)
      timee=get("time",envir)
    }
    
    nn=nn*2
    
    params = interpret.theta()
    phis = params$phis
    kappa2s = params$kappa2_eps_V
    tau=params$tau
    C=params$C
    rhos=params$rhos
    
    
    # nn=8
    # phis = c(12, 9, 2, 4)
    # kappa2s = c(3, 2, 6, 5)
    # tau=2
    # C=3
    # rhos = c(6, 16, 9, 3)
    
    
    
    ii = c(1,(nn/2),((nn/2)+1),nn,
           2:(nn/2-1), (nn/2+2):(nn-1),
           1:((nn/2)-1),((nn/2)+1):(nn-1),
           1:((nn/2)), 
           1:((nn/2)-1))
    
    jj = c(1,(nn/2),((nn/2)+1),nn,
           2:(nn/2-1), (nn/2+2):(nn-1),
           2:(nn/2),((nn/2)+2):nn, 
           ((nn/2)+1):(nn),  
           ((nn/2)+2):(nn))
    
    
    
    print("Q")
    cat("range: ", range(kappa2s),"\n")
    cat("range: ", range(phis),"\n")
    
    
    
    xx = c(tau+tau*phis[2]^2,tau,C^2*tau+kappa2s[2]*rhos[2]^2+kappa2s[1],kappa2s[nn/2]+tau*C^2,
           tau+tau*phis[3:(nn/2)]^2, C^2*tau+kappa2s[2:(nn/2-1)]+kappa2s[3:(nn/2)]*rhos[3:(nn/2)]^2,
           -tau*phis[2:(nn/2)],-kappa2s[2:(nn/2)]*rhos[2:(nn/2)],
           rep(-C*tau,(nn/2)),
           C*tau*phis[2:(nn/2)])
    
    cat("range x: ", range(xx),"\n")
    Q = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
    return (Q)
  }
  
  
  log.norm.const = function(){return(numeric(0))}
  
  log.prior = function(){
    if(!is.null(envir)){
      if(!is.null(envir[["my.log.prior"]])){
        my.prior=get("my.log.prior",envir)
      }else{
        my.prior=NULL
      }
    }else{
      my.prior=NULL
    }
    if(!is.null(my.prior)){
      lprior = my.prior(theta)
    }else{
      params = interpret.theta()
      
      #Gaussian_priors
      
      # lprior = dnorm(theta[1], mean=0, sd=1)
      # lprior = lprior + dnorm(theta[2], mean=0, sd=1)
      # lprior = lprior + dnorm(theta[3], mean=0, sd=1)
      # lprior = lprior + dnorm(theta[4], mean=0, sd=1)
      # lprior = lprior + dnorm(theta[5], mean=0, sd=1)
      
      #Uniform_priors
      
      lprior = -theta[1] -2*log(1+exp(-theta[1]))
      lprior = lprior -theta[2] -2*log(1+exp(-theta[2]))
      lprior = lprior + dgamma(exp(theta[3]), shape=1, rate=0.1) + theta[3]
      lprior = lprior -theta[4] -2*log(1+exp(-theta[4]))
      lprior = lprior -theta[5] -2*log(1+exp(-theta[5]))
      
    }
    return (lprior)
  }
  
  
  initial = function(){
    return (numeric(5)) 
  }
  quit = function(){return ()  }
  if(is.null(theta)){
    theta = initial()
  }
  cmd = match.arg(cmd)
  val = do.call(cmd, args = list())
  return (val)
  
  
}




#' rgeneric model for AR(2) time-dependent process
#'
#' Defines the rgeneric model structure for the AR(2)/nested AR(1) model with time-dependent
#' correlation. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
#' This is intended for internal use only, but the documentation is included here in case someone want to change something.
#'
#' @param cmd Vector containing list of function names necessary for the rgeneric model.
#' @param theta Vector describing the hyperparameters in internal scaling.
#' @importFrom stats dnorm dgamma
rgeneric.ar2.forcing = function(
    cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{
  
  
  envir = environment(sys.call()[[1]])
  
  
  interpret.theta = function() {
    if(!is.null(envir)){
      timee=get("time",envir)
      nn=get("n",envir)
    }
    tau = exp(15) 
    r = diff(range(timee))
    b_Y = -1/r+2/r*1/(1+exp(-theta[1]))
    low_Y = min(b_Y*timee[1],b_Y*timee[nn])
    high_Y = max(b_Y*timee[1],b_Y*timee[nn])
    a_Y = -low_Y + (1-high_Y+low_Y)/(1+exp(-theta[2]))
    lambdas_Y = -log(a_Y+b_Y*timee)
    cc_Y = 1/(nn-1)
    phis = c(exp(-lambdas_Y*c(1,diff(timee)/cc_Y)) ) 
    
    kappa_eps_V = exp(theta[3])
    kappa2_eps_V=kappa_eps_V*2*lambdas_Y
    
    
    
    
    C=1
    b_V = -1/r+2/r*1/(1+exp(-theta[4]))
    low_V = min(b_V*timee[1],b_V*timee[nn])
    high_V = max(b_V*timee[1],b_V*timee[nn])
    a_V = -low_V + (1-high_V+low_V)/(1+exp(-theta[5]))
    lambdas = -log(a_V+b_V*timee)
    rhos=a_V+b_V*timee
    cc = 1/(nn-1)
    rhos = c(exp(-lambdas*c(1,diff(timee)/cc)) ) 
    
    
    kappa_f = exp(theta[6])
    kappa2fs = kappa_f*2*lambdas_Y
    #F0 = theta[7]
    F0 = 0
    
    print(theta)
    return(list(phis = phis,a_Y=a_Y,b_Y=b_Y,
                tau=tau,
                C=C,
                lambdas_Y=lambdas_Y,
                rhos=rhos,a_V=a_V,b_V=b_V,
                kappa_eps_V = kappa_eps_V, 
                kappa2_eps_V=kappa2_eps_V,
                kappa_f=kappa_f, kappa2fs=kappa2fs,F0=F0))
  }
  
  
  mu = function() {
    
    if(!is.null(envir)){
      nn=get("n",envir)
      timee=get("time",envir)
      fforcing=get("forcing",envir)
    }
    #require("INLA.ews",quietly=TRUE)
    hyperparam = interpret.theta()
    
    sf = 1/sqrt(hyperparam$kappa_f)
    kappa_f = hyperparam$kappa_f
    kappa2fs = hyperparam$kappa2fs
    #F0 = fforcing[1]
    F0= hyperparam$F0
    
    zz = (fforcing+F0)/sqrt(kappa_f)
    
    innerstruct = timee
    llambdas = hyperparam$lambdas_Y
    
    
    muvek = numeric(nn)
    
    
    for(k in 1:nn){
      for(s in 1:k){
        muvek[k] = muvek[k] + zz[s]*exp(-llambdas[k]*(timee[k]-timee[s]))
      }
    }
    
    
    muvek = muvek*(1/(sqrt(2*hyperparam$lambdas_Y)))
    
    return(c(muvek,numeric(nn)))
    
    
  }
  
  
  graph = function() {
    if(!is.null(envir)){
      nn=get("n",envir)
    }else{
      nn=get("n",environment())
    }
    nn=nn*2
    
    
    
    ii = c(1,(nn/2),((nn/2)+1),nn,
           2:(nn/2-1), (nn/2+2):(nn-1),
           1:((nn/2)-1),((nn/2)+1):(nn-1),
           1:((nn/2)), 
           1:((nn/2)-1))
    
    jj = c(1,(nn/2),((nn/2)+1),nn,
           2:(nn/2-1), (nn/2+2):(nn-1),
           2:(nn/2),((nn/2)+2):nn, 
           ((nn/2)+1):(nn),  
           ((nn/2)+2):(nn))
    
    
    
    xx=c(rep(1,length(ii)))
    G = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
    G[G != 0] = 1
    return (G)
  }
  
  
  
  Q = function(){
    if(!is.null(envir)){
      nn=get("n",envir)
      timee=get("time",envir)
    }
    
    nn=nn*2
    
    params = interpret.theta()
    phis = params$phis
    kappa2s = params$kappa2_eps_V
    tau=params$tau
    C=params$C
    rhos=params$rhos
    
    
    
    
    ii = c(1,(nn/2),((nn/2)+1),nn,
           2:(nn/2-1), (nn/2+2):(nn-1),
           1:((nn/2)-1),((nn/2)+1):(nn-1),
           1:((nn/2)), 
           1:((nn/2)-1))
    
    jj = c(1,(nn/2),((nn/2)+1),nn,
           2:(nn/2-1), (nn/2+2):(nn-1),
           2:(nn/2),((nn/2)+2):nn, 
           ((nn/2)+1):(nn),  
           ((nn/2)+2):(nn))
    
    
    
    print("Q")
    cat("range: ", range(kappa2s),"\n")
    cat("range: ", range(phis),"\n")
    
    
    
    xx = c(tau+tau*phis[2]^2,tau,C^2*tau+kappa2s[2]*rhos[2]^2+kappa2s[1],kappa2s[nn/2]+tau*C^2,
           tau+tau*phis[3:(nn/2)]^2, C^2*tau+kappa2s[2:(nn/2-1)]+kappa2s[3:(nn/2)]*rhos[3:(nn/2)]^2,
           -tau*phis[2:(nn/2)],-kappa2s[2:(nn/2)]*rhos[2:(nn/2)],
           rep(-C*tau,(nn/2)),
           C*tau*phis[2:(nn/2)])
    
    cat("range x: ", range(xx),"\n")
    Q = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
    return (Q)
  }
  
  
  log.norm.const = function(){return(numeric(0))}
  log.prior = function(){
    if(!is.null(envir)){
      if(!is.null(envir[["my.log.prior"]])){
        my.prior=get("my.log.prior",envir)
      }else{
        my.prior=NULL
      }
    }else{
      my.prior=NULL
    }
    if(!is.null(my.prior)){
      lprior = my.prior(theta)
    }else{
      params = interpret.theta()
      
      #Gaussian_priors
      
      if(FALSE){
        lprior = dnorm(theta[1], mean=0, sd=1)
        lprior = lprior + dnorm(theta[2], mean=0, sd=1)
        lprior = lprior + dnorm(theta[3], mean=0, sd=1)
        lprior = lprior + dnorm(theta[4], mean=0, sd=1)
        lprior = lprior + dnorm(theta[5], mean=0, sd=1)
        lprior = lprior + dnorm(theta[6], mean=0, sd=1)
      }else{
        lprior = -theta[1] -2*log(1+exp(-theta[1])) #b_phi
        lprior = lprior -theta[2] -2*log(1+exp(-theta[2]))#a_phi
        lprior = lprior + dgamma(exp(theta[3]), shape=1, rate=0.1) + theta[3] #kappa_x
        lprior = lprior -theta[4] -2*log(1+exp(-theta[4])) #b_rho
        lprior = lprior -theta[5] -2*log(1+exp(-theta[5]))#a_rho
        # lprior = lprior + dnorm(theta[6], mean=0, sd=1) #
        lprior = lprior + dgamma(exp(theta[6]), shape=1, rate=0.1) + theta[6] #kappa_f
      }
      
      
      
      #lprior = lprior + dnorm(theta[7], mean=0, sd=1)
      #Uniform_priors
      # 
      
      # 
    }
    
    return (lprior)
  }
  
  
  initial = function(){
    return (numeric(6)) 
  }
  quit = function(){return ()  }
  if(is.null(theta)){
    theta = initial()
  }
  cmd = match.arg(cmd)
  val = do.call(cmd, args = list())
  return (val)
  
  
}

