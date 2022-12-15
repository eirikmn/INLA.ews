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
#' @importFrom stats dnorm
rgeneric.ews.fgn = function(
    cmd = c("graph", "Q","mu", "initial", "log.norm.const", "log.prior", "quit"),
    theta = NULL)
{

  
  # require("INLA.climate2",quietly=TRUE)

  tau = exp(15)
  envir = environment(sys.call()[[1]])

  hyperg = function(a,b,c,z,trunc){
    
    kk = 0:(trunc-1)
    
    cp = cumprod((a+kk)*(b+kk)/(c+kk)*z/(1:trunc)  ) 
    
    s = 1 + sum(cp)
    
    return(s)
    
  }
  
  
  incBeta = function(z,a,b,trim){
    
    hg = 1/a*z^a*hyperg(a,1-b,1+a,z,trunc=trim) #take this with a grain of salt, found by trial-and-error
    
    return(hg)
  }
  
  gammas = function(a,b,n,t,s){
    if((b *(s + t - 2 *max(s, t)))/n==0){
      0
    }
    else{
      ( gamma(-1 + 2 *a + (b *(s + t))/n)*gamma( 1 - 2 *a - (2*b *max(s, t))/n))/gamma((b *(s + t - 2 *max(s, t)))/n)
    }
  }
  
  
  
  R2test = function(a,b,n,t,s){
    
    (1/(2*a*n+b*(s+t)))*n*( sqrt( (a + (b* max(s,t)/n))*(-1 + 2 *a + (2* b* max(s,t)/n))/beta(
      2 - 2 *(a + (b *max(s,t)/n)), -(1/2) + a + (b *max(s,t)/n)))*sqrt( (a + (b* min(s,t)/n))*(-1 + 2 *a + (2* b* min(s,t)/n))/beta(
        2 - 2 *(a + (b *min(s,t)/n)), -(1/2) + a + (b *min(s,t)/n))))*(beta(-(1/2) + a + (b *max(s, t))/n, -((
          2* (-1 + a)* n + b* max(s, t) + b *min(s, t))/n))* beta((n + b* max(s, t) - b* min(s, t))/ n, ((-1 + 2* a)* n + b* max(s, t) + b* min(s, t))/n)* min(s, t)^( 2* a + (b* (s + t))/n) + 
            beta(-(1/2) + a + (b *min(s, t))/n, -((2*(-1 + a)*n +b*max(s, t) + b*min(s, t))/n))*
            ((-incBeta(min(s, t)/max(s, t),  1-2*a-(2*b*max(s, t))/n, -1 + 2 *a + (b *(s + t))/n,100) 
              + gammas(a,b,n,t,s)
              
            )*min(s, t)^(2 *a + (b* (s + t))/n) + (n *
                                                     hyperg(2 - 2* a - (b *(s + t))/n, (n - b* max(s, t) + b *min(s, t))/n, (2 *n - b* max(s, t) + b* min(s, t))/n, min(s, t)/max(s, t),100) 
                                                   *max(s, t)^(2*a+(b*(s + t))/n)*(min(s, t)/max(s, t))^(1+(b*(-max(s, t)+min(s, t)))/n))/(n-b*max(s, t)+b*min(s, t))))
  }
  
  Rfgn = function(a,b,n,t,s){
    R2test(a,b,n,t+1,s+1) - R2test(a,b,n,t+1,s) - R2test(a,b,n,t,s+1) + R2test(a,b,n,t,s)
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

  # greensH = function(t,s,a,b,n){
  #   if(t-s<0){
  #     return(0)
  #   }else{
  #     H = a+b*max(t,s)/n
  #     # H = a+b*min(t,s)/n
  #     # H = a+b*(t+s)/2/n
  #     ret = (t-s+0.5)^(H-3/2)
  #     return(ret)
  #   }
  # }
  
  
  Q = function()  {
    if(!is.null(envir)){
      nn=get("n",envir)
    }

    hyperparam = interpret.theta()
    Hs = hyperparam$Hs
    kappa = hyperparam$kappa
    sx = 1/sqrt(hyperparam$kappa)
    a = hyperparam$a
    b=hyperparam$b
    # H2 = 2*Hs
    # k=0:(nn-1)
    # # sigmat = matrix(NA,nn,nn)
    # for(i in 1:nn){
    #   for(j in 1:nn){
    #     #t = max(i,j)
    #     t = (i+j)/2/nn
    #     H2 = 2*(hyperparam$a+hyperparam$b*t)
    #     k=abs(i-j)
    #     sigmat[i,j] = sx^2/2*( abs(k-1)^H2-2*abs(k)^H2+abs(k+1)^H2 )
    #   }
    # }
    # #sigmat = sigmamaker(nn,sx,Hs)

    # Gmat = matrix(NA,nn,nn)
    # for(i in 1:nn){
    #   for(j in 1:nn){
    #     Gmat[i,j] = greensH(i,j,a,b,nn)
    #   }
    # }
    # covmat = sx^2*Gmat%*%t(Gmat)
    
    covmat = data.frame(row.names = NULL,col.names=NULL)
    
    for (i in 1:nn) {
      for (j in 1:nn) {
        covmat[i,j] =round( Rfgn(a,b,nn,i,j),3)
      }
    }    
    
    # url = paste0("https://www.wolframcloud.com/obj/e965ae9a-187b-469b-9d5c-c9f5eaa14656?n1=",
    #              nn,"&a1=",a,"&b1=",b,"&sigma1=",sx)
    # 
    # df = scan(url)
    # 
    # covmat = matrix(df,ncol=n)
    
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



#' rgeneric model for fgn time-dependent process
#'
#' Defines the rgeneric model structure for the fractional Gaussian noise model 
#' with time-dependent correlation. This model includes forcing. This includes key functions that provide information on precision matrix, mean vector, priors, graph etc.
#' This is intended for internal use only, but the documentation is included here in case someone want to change something.
#'
#' @param cmd Vector containing list of function names necessary for the rgeneric model.
#' @param theta Vector describing the hyperparameters in internal scaling.
#'
#'
#' @importFrom stats dnorm
rgeneric.ews.fgn.forcing = function(
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
    Hs = params$Hs
    #cat("F0: ",params$F0," sigmaf: ",1/sqrt(params$kappa_f),"\n")
    zz = 1/sqrt(params$kappa_f)*(params$F0+z)
    
      # struct = (1:nn-0.5)^(Hs-3/2)
      # muvek0=numeric(nn)
      # for(i in 1:nn){
      #   muvek0[i] = rev(struct[1:i])%*%zz[1:i]
      # }
    # # greensmat = greensmaker_fgn(params$Hs)
    # muvek = greensmat%*%zz
    muvek = numeric(nn)
    for(i in 1:nn){
      for(j in 1:i){
        muvek[i] = muvek[i] + zz[j]*(i-j+0.5)^(Hs[i]-3/2) #struct[i-j+1]
      }
    }
    # 
    # if(Sys.info()[['sysname']] == "Darwin"){
    # compute_mu_fgn(muvek, nn, zz,  Hs)
    # }else if(Sys.info()[['sysname']] == "Linux"){
    #    for(i in 1:nn){
    #      for(j in 1:i){
    #        muvek[i] = muvek[i] + zz[j]*(i-j+0.5)^(Hs[i]-3/2) #struct[i-j+1]
    #      }
    #    }
    # }else if(Sys.info()[['sysname']] == "Windows"){
    #   for(i in 1:nn){
    #     for(j in 1:i){
    #       muvek[i] = muvek[i] + zz[j]*(i-j+0.5)^(Hs[i]-3/2) #struct[i-j+1]
    #     }
    #   }
    # }
     # muvek=numeric(nn)
     # #struct=numeric(nn)
     # struct = (1:nn-0.5)^(Hs-3/2)
     # for(i in 1:nn){
     #   for(j in 1:i){
     #     muvek[i] = muvek[i] + zz[j]*(i-j+0.5)^(Hs[i]-3/2) #struct[i-j+1]
     #   }
     # }
    
    
    #cat("mu range:",range(muvek))
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
    # H2 = 2*Hs
    k=0:(nn-1)
    sigmat = matrix(NA,nn,nn)
    for(i in 1:nn){
      for(j in 1:nn){
        # t = max(i,j)
        t = (i+j)/2/nn
        H2 = 2*(hyperparam$a+hyperparam$b*t)
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

