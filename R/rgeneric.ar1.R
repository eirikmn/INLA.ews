#' rgeneric model for AR(1) time-dependent process
#'
#' Defines the rgeneric model structure for the AR(1) model with time-dependent
#' correlation. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
#' This is intended for internal use only, but the documentation is included here in case someone want to change something.
#'
#' @param cmd Vector containing list of function names necessary for the rgeneric model.
#' @param theta Vector describing the hyperparameters in internal scaling.
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
    cat("a: ",a,"b: ",b,"sigma: ",1/sqrt(kappa_eps),"\n",sep="")
    cc = 1/(nn-1)
    phis = c(exp(-lambdas*c(1,diff(timee)/cc)) ) #rescale
    
    #kappa_f = exp(theta[4])
    #F0 = theta[5]
    print(theta)
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
    print("Q")
    cat("range: ", range(kappa2s),"\n")
    cat("range: ", range(phis),"\n")
    
    #xx = kappa_eps*c(1+phis[2]^2,1,1+phis[3:nn]^2,-phis[2:nn])
    xx = c(kappa2s[1]*(1-phis[1]^2)+kappa2s[2]*phis[2]^2, kappa2s[nn],
           kappa2s[2:(nn-1)]+kappa2s[3:nn]*phis[3:nn]^2,
           -phis[2:nn]*kappa2s[2:nn])
    
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
      lprior = INLA::inla.pc.dprec(params$kappa_eps, u=1, alpha=0.01, log=TRUE) + log(params$kappa_eps) #kappa
      lprior = lprior + dnorm(theta[2],sd=1,log=TRUE)
      lprior = lprior + dnorm(theta[3],sd=1,log=TRUE) #theta_a
      
    }
    
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
#' @importFrom Rcpp sourceCpp
#' @useDynLib INLA.ews, .registration = TRUE
rgeneric.ar1.forcing = function(
    cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{
  require(INLA.ews, quietly=TRUE)
  #print(c_mu_ar1)
  tau = exp(15)
  envir = parent.env(environment())
  
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
    require("INLA.ews",quietly=TRUE)
    #print("mu")
    hyperparam = interpret.theta()
    #print(hyperparam)
    
    sf = 1/sqrt(hyperparam$kappa_f)
    kappa2fs = hyperparam$kappa2fs
    F0 = hyperparam$F0
    
    zz = (fforcing+F0)/sqrt(kappa2fs)
    
    innerstruct = timee
    #solution 1:
    llambdas = hyperparam$lambdas
    
    #soluti
    
    struktur = exp(-hyperparam$lambdas*innerstruct)
    
    muvek = numeric(nn)
    
    
    
    if(diff(range(diff(timee)))<10^(-12)){
      #print("regular")
      for(i in 1:nn){
        muvek[i] = rev(struktur[1:i])%*%zz[1:i]
      }
    }else{
      if(TRUE){
        #print("irregular")
        for(k in 1:nn){
          for(s in 1:k){
            muvek[k] = muvek[k] + zz[s]*exp(-llambdas[k]*(timee[k]-timee[s]))
          }
        }
      }else{
        
        
        c_mu_ar1(muvek,zz, nn,llambdas,timee)
        #cfunk(muvek,zz, nn,llambdas,timee)
        #if(!is.loaded('Rc_mu_ar1')){
          #dyn.load(file.path(.Library,"INLA.climate/libs/Rc_Q.so"))
        #  dyn.load(file.path("INLA.ews.so"))
        #}
        #res = .C('c_mu_ar1',mu=as.matrix(muvek,ncol=1),
        #         as.double(zz),as.integer(nn),
        #         as.double(llambdas),
        #         as.double(timee))
        #print(res$mu[1:3])
        #muvek=res$mu
        
        #.Call('_INLA_ews_c_mu_ar1', PACKAGE = 'INLA.ews', 
        #      muvek, zz, nn, llambdas, timee)
        #cfunk(muvek,z=zz, n=nn,lambda=llambdas,time=timee)
        #if(!is.loaded('Rc_mu_ar1')){
          #dyn.load(file.path(.Library,"INLA.climate/libs/Rc_Q.so"))
        #  dyn.load(file.path("Rc_mu_ar1.so"))
        #  cat("dyn.load\n")
        #}
        #res = .C('Rc_mu_ar1',mu=as.matrix(means,ncol=1),as.double(fforcing),as.integer(nn),as.integer(NN),
        #         as.double(weights),as.double(llambdas),as.double(sf),
        #         as.double(hyperparam$F0))
        
        #Cres = .C('Rc_mu_ar1',mu=as.matrix(numeric(nn),ncol=1),as.double(zz),as.integer(nn),
        #         as.double(llambdas), as.double(timee)) #new
        #cat("done C\n")
        #cat("res: ", Cres$mu,"\n")
        #muvek=Cres$mu 
      }
      
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
      lprior = INLA::inla.pc.dprec(params$kappa_eps, u=1, alpha=0.01, log=TRUE) + log(params$kappa_eps) #kappa
      #lprior = lprior + dnorm(theta[2],log=TRUE) #theta_b
      b=params$b; ra = 0; rb=1
      lprior = -log(rb-ra)-theta[2]-log(1+exp(-theta[2]))
      #lprior = lprior -log(rb-ra) + log(b-ra)+log(rb-b) -log(rb-ra)
      lprior = lprior + dnorm(theta[3],sd=3,log=TRUE) #theta_a
      lprior = lprior + INLA::inla.pc.dprec(params$kappa_f, u=1, alpha=0.01, log=TRUE) + log(params$kappa_f) #kappa_f
      lprior = lprior + dnorm(theta[5],sd=3,log=TRUE) #F0
    }
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

