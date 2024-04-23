if(FALSE){
  set.seed(123)
  n=500
  a = 0.3
  b = 0.2
  
  ttime = sort(1:n + rnorm(n,sd=0.1))
  ttime = 1:n
  time_norm = (ttime-min(ttime))/(max(ttime)-min(ttime))
  lambdas = -log(a+b*time_norm)
  cc=1/(n-1)
  phis = exp(-lambdas*c(1,diff(time_norm)/cc))
  
  plot(phis)
  lines(exp(-lambdas),col="red",lwd=3)
  
  #forcing = arima.sim(list(ar=c(0.95)), n=n)
  
  #phis=rep(phi,n)
  F0=0
  sigma=5
  #sigma_f = 0.1
  
  
  sims=numeric(n)
  sims[1] = rnorm(1,sd=sigma/sqrt(2*lambdas[1]))
  for(i in 2:n){
    sims[i] = phis[i]*sims[i-1] + rnorm(1,sd=sigma/sqrt(lambdas[i]))
  }
  #sims = sigma * arima.sim(list(ar=c(phi)), n=n, sd=sqrt(1-phi^2))
  #muvek = mucomputer(pars=c(sigma,b,a,sigma_f,F0),forcing,time_norm,as.theta=FALSE)
  
  intercept=0
  y = sims #+ muvek + intercept
  plot(y)
  #lines(muvek,col="red",lwd=5)
  #
  
  rrr = inla.ews(data=y,formula=y~1, stepsize=0.01)
  summary(rrr)
  #
  #
  
  rgm = inla.rgeneric.define(rgeneric.ar1, n=n, time=time_norm)
  
  rrm = inla(formula=y~-1 + f(idy,model=rgm), data=data.frame(y=y,idy=1:n), 
             control.family=list(initial=12,fixed=TRUE), verbose=TRUE,
             num.threads=inla.options$num.threads,
             control.mode=list(restart=TRUE),
             control.inla=list(h=0.005)
             )
  
  summary(rrm)
  rres = from.theta(rrm$summary.hyperpar$mean,time_norm = time_norm)
  rres$a
  rres$b
  #
  
  
  cmodel <- INLA::inla.cgeneric.define(model = "inla_cgeneric_timedep",
                                       shlib = "src/cgeneric.so", 
                                       n = n, time=as.numeric(time)
  )
  rc <- inla(
    y ~ -1 + f(idx, model = cmodel),
    data = data.frame(y, idx = 1:n),
    control.family = list(hyper = list(prec = list(initial = 12, fixed = TRUE))))
  
  summary(rc)
  
  
  ####
  
  
  rrr = inla.ews(y,forcing=forcing,formula=-1,timesteps=time_norm)
  
  rrm = inla.rgeneric.define(rgeneric.ar1.2, 
                             n=n,
                             time=time
  )
  
  rr <- inla(
    y ~ -1 + f(idx, model = rrm),
    data = data.frame(y, idx = 1:n),verbose=FALSE,
    control.family = list(hyper = list(prec = list(initial = 12, fixed = TRUE))))
  summary(rr)
  
  #
  
  
  
  
  nsims=30000
  testr = rr
  hyps = inla.hyperpar.sample(nsims,testr)
  
  bsims = -1 + 2/(1+exp(-hyps[,2]))
  asims= numeric(0)
  
  phisims = matrix(NA,nrow=n,ncol=nsims)
  phi0sims = matrix(NA,nrow=n,ncol=nsims)
  
  
  for(i in 1:nsims){
    alower = -min(bsims[i],0)
    aupper = 1-max(bsims[i],0)
    asims[i] = alower + (aupper-alower)/(1+exp(-hyps[i,3]))
    phi0sims[,i] = asims[i]+bsims[i]*time_norm
    
  }
  resmatr = data.frame(a=mean(asims),b=mean(bsims))
  
  testr = rc
  hyps = inla.hyperpar.sample(nsims,testr)
  
  bsims = -1 + 2/(1+exp(-hyps[,2]))
  asims= numeric(0)
  
  phisims = matrix(NA,nrow=n,ncol=nsims)
  phi0sims = matrix(NA,nrow=n,ncol=nsims)
  
  
  for(i in 1:nsims){
    alower = -min(bsims[i],0)
    aupper = 1-max(bsims[i],0)
    asims[i] = alower + (aupper-alower)/(1+exp(-hyps[i,3]))
    phi0sims[,i] = asims[i]+bsims[i]*time_norm
    
  }
  
  resmatc = data.frame(a=mean(asims),b=mean(bsims))
  print(resmatr); print(resmatc)
  
  
  rgeneric.ar1.2 = function(
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
      phis = c(exp(-lambdas*c(0,diff(timee)/cc)) ) #rescale
      
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
      lprior = -params$kappa_eps + theta[1] - 
        0.5 * log(2.0 * pi) - 0.5 * theta[3]^2 - 
        0.5 * log(2.0 * pi) - 0.5 * theta[2]^2 
      return (lprior)
    }
    initial = function(){
      return (c(1.0,1.0,1.0)) 
    }
    
    quit = function(){return ()  }
    if(is.null(theta)){
      theta = initial()
    }
    cmd = match.arg(cmd)
    val = do.call(cmd, args = list())
    return (val)
  }
  
}
