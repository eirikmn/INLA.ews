if(FALSE){
  rgeneric.ar1.bp = function(
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
      
      
      
      bp = 1/(1+exp(-theta[4]))
      #bp=theta[4]
      cat("bp:",bp,"\n")
      cat("theta:",theta,"\n")
      #bpafter = 1-theta[4]
      nbefore = floor(bp*nn) +1
      nafter = nn-nbefore
      if(nbefore <= 1){
        timecov = timee
      }else if(nbefore >=nn){
        timecov = numeric(nn)
      }else{
        timecov = timee-timee[nbefore]
        timecov[1:nbefore]=0
      }
      cat("nbefore",nbefore,"\n")
      
      cat("nafter",nafter,"\n")
      
      cat("b:",b,"\n")
      cat("a:",a,"\n")
      cat("high",high,"low",low,"\n")
      
      # timecov = c(numeric(nbefore), timee[(nbefore+1):nn])
      # timecov[(nbefore+1):nn] = timecov[(nbefore+1):nn]-timecov[nbefore+1]
      # timecov[(nbefore+1):nn] = timecov[(nbefore+1):nn]/timecov[nn]
      # 
      #lambdas = -log(a+b*timee)
      lambdas = -log(a+b*timecov)
      kappa2s = kappa_eps*2*lambdas
      
      cc = 1/(nn-1)
      #phis = c(exp(-lambdas*c(1,diff(timee)/cc)) ) #rescale
      phis = a+b*timecov #does not account for non-constant time steps!!!
        
      cat("a:",a,"b:",b,"lambdas:",range(lambdas),"\n")
      cat("-loga",-log(a),"-log a+b*timecov",-log(a+b*timecov[nn]))
      cat("range:",range(kappa2s))
      cat("range:",range(phis))
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
      #cat("range: ", range(kappa2s),"\n")
      #cat("range: ", range(phis),"\n")
      
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
      
        params = interpret.theta()
        lprior = INLA::inla.pc.dprec(params$kappa_eps, u=1, alpha=0.01, log=TRUE) + log(params$kappa_eps) #kappa
        lprior = lprior + dnorm(theta[2],sd=1,log=TRUE)
        lprior = lprior + dnorm(theta[3],sd=1,log=TRUE) #theta_a
        lprior = lprior + dnorm(theta[4],mean=2,sd=1,log=TRUE) #bp
        #lprior = lprior + dnorm(theta[4],mean=0.5,sd=0.1, log=TRUE)
      
      return (lprior)
    }
    initial = function(){
      return (c(0.,0.,0.,1.)) 
    }
    
    quit = function(){return ()  }
    if(is.null(theta)){
      theta = initial()
    }
    cmd = match.arg(cmd)
    val = do.call(cmd, args = list())
    return (val)
  }
  
  set.seed(1)
  n = 1000
  sigma = 1
  a=0.2
  b=0.5
  bp = 0.25
  nbefore = floor(n*bp)
  nafter=n-nbefore
  time=1:n
  time_norm = seq(0,1,length.out=n)
  nbefore = floor(bp*nn) +1
  nafter = nn-nbefore
  if(nbefore <= 1){
    timecov = timee
  }else if(nbefore >=nn){
    timecov = numeric(nn)
  }else{
    timecov = timee-timee[nbefore]
    timecov[1:nbefore]=0
  }
  phis = a+b*timecov
  
  data=ar1_timedep_sim(n,sigma=sigma,phis=phis)
  plot(data)  
  lines(phis*max(data),col="red",lwd=4)
  df = data.frame(y=data, idy=1:n)
  
  rgen_model = inla.rgeneric.define(rgeneric.ar1.bp,n=n,time=time_norm)
  formula = y~-1+f(idy,model=rgen_model)
  r = inla(formula,data=df,control.family = list(initial=12,fixed=TRUE))
  summary(r)
  
  sigma_est=INLA::inla.emarginal(function(x)1/sqrt(exp(x)),r$marginals.hyperpar$`Theta1 for idy`)
  sigma_marg=INLA::inla.tmarginal(function(x)1/sqrt(exp(x)),r$marginals.hyperpar$`Theta1 for idy`)
  
  rekke = diff(range(time_norm))
  b_est=INLA::inla.emarginal(function(x) -1/rekke+2/rekke*1/(1+exp(-x)) ,r$marginals.hyperpar$`Theta2 for idy` )
  b_marg=INLA::inla.tmarginal(function(x) -1/rekke+2/rekke*1/(1+exp(-x)) ,r$marginals.hyperpar$`Theta2 for idy`)
  b_positive_prob = 1-INLA::inla.pmarginal(0,b_marg)
  b_zmarg = INLA::inla.zmarginal(b_marg,silent=TRUE)
  b_zmarg$mode = INLA::inla.mmarginal(b_marg)
  b_low = b_zmarg$quant0.025
  b_high= b_zmarg$quant0.975
  bp_marg = inla.tmarginal(function(x)1/(1+exp(-x)), r$marginals.hyperpar$`Theta4 for idy`)
  bp_est = inla.emarginal(function(x)1/(1+exp(-x)), r$marginals.hyperpar$`Theta4 for idy`)
  
  nsims = 10000
  hypersamples = INLA::inla.hyperpar.sample(nsims,r)

  b_sims = -1/rekke+2/rekke*1/(1+exp(-hypersamples[,2]))
  #b = -1/r+2/r*1/(1+exp(-theta[2]))
  a_sims = numeric(nsims)
  
  phi_sims = matrix(NA,nrow=n,ncol=nsims)
  phi0_sims = matrix(NA,nrow=n,ncol=nsims)
  lambda_sims = matrix(NA,nrow=n,ncol=nsims)
  bp_sims = 1/(1+exp(-hypersamples[,4]))
  #bp_sims=hypersamples[,4]
  for(i in 1:nsims){
    bps = bp_sims[i]
    nbefore = floor(bps*n) +1
    nafter = nn-nbefore
    if(nbefore <= 1){
      timecov = time_norm
    }else if(nbefore >=n){
      timecov = numeric(n)
    }else{
      timecov = time_norm-time_norm[nbefore]
      timecov[1:nbefore]=0
    }
    
    low=min(b_sims[i]*time_norm[1],b_sims[i]*time_norm[n])
    high=max(b_sims[i]*time_norm[1],b_sims[i]*time_norm[n])
    a_sims[i] = -low + (1-high+low)/(1+exp(-hypersamples[i,3]))
    
    lambda_sims[,i] = -log(a_sims[i]+b_sims[i]*timecov)
    #phi_sims[,i] = a_sims[i]+b_sims[i]*df$time_normalized
    cc = 1/(length(time_norm-1))
    
    phi_sims[,i] = exp(-lambda_sims[,i]*c(0,diff(timecov))/cc)
    phi0_sims[,i] = a_sims[i]+b_sims[i]*timecov
    
    
  }
  phi_sims=phi0_sims
  a_marg = density(a_sims); a_marg=data.frame(x=a_marg$x,y=a_marg$y)
  a_zmarg = INLA::inla.zmarginal(a_marg,silent=TRUE)
  a_zmarg$mode = INLA::inla.mmarginal(a_marg)
  
  phi_means = rowMeans(phi_sims)
  phi_median = numeric(n)
  phi_sd = numeric(n)
  phi_lower = numeric(n)
  phi_upper = numeric(n)
  phi_qlower = numeric(n)
  phi_qmid = numeric(n)
  phi_qupper = numeric(n)
  phi_mode = numeric(n)
  phi0_means = rowMeans(phi0_sims)
  phi0_median = numeric(n)
  phi0_sd = numeric(n)
  phi0_lower = numeric(n)
  phi0_upper = numeric(n)
  phi0_qlower = numeric(n)
  phi0_qmid = numeric(n)
  phi0_qupper = numeric(n)
  phi0_mode = numeric(n)
  for(i in 1:n){
    dens = density(phi_sims[i,])
    dens = data.frame(x=dens$x,y=dens$y)
    zmarg = INLA::inla.zmarginal(dens,silent=TRUE)
    phi_median[i]=zmarg$quant0.5; phi_sd[i] = zmarg$sd
    hpds = INLA::inla.hpdmarginal(0.95,dens)
    phi_lower[i]=hpds[1]; phi_upper[i]=hpds[2]
    phi_qlower[i]=zmarg$quant0.025
    phi_qmid[i]=zmarg$quant0.5
    phi_qupper[i]=zmarg$quant0.975
    phi_mode[i] = INLA::inla.mmarginal(dens)
    
    dens = density(phi0_sims[i,])
    dens = data.frame(x=dens$x,y=dens$y)
    zmarg = INLA::inla.zmarginal(dens,silent=TRUE)
    phi0_median[i]=zmarg$quant0.5; phi0_sd[i] = zmarg$sd
    hpds = INLA::inla.hpdmarginal(0.95,dens)
    phi0_lower[i]=hpds[1]; phi0_upper[i]=hpds[2]
    phi0_qlower[i]=zmarg$quant0.025
    phi0_qmid[i]=zmarg$quant0.5
    phi0_qupper[i]=zmarg$quant0.975
    phi0_mode[i] = INLA::inla.mmarginal(dens)
  }
  phi_means[1] = NA
  phi_median[1] = NA
  phi_lower[1] = NA
  phi_upper[1] = NA
  phi_qmid[1] = NA
  phi_qlower[1] = NA
  phi_qupper[1] = NA
  phi_mode[1] = NA
  

  
  cat("a:  ",mean(a_sims)," (",a,")\n","b:  ",b_est," (",b,")\n","bp: ",bp_est," (",bp,")\n",sep="")

  library(ggplot2)
  
  lambdas = -log(a+b*timecov)
  stds = sqrt(sigma^2/lambdas)
  
  ggd = data.frame(time=time_norm, phis=phis, phimean = phi_means, philower=phi_lower, 
                   phiupper=phi_upper, y=data, stdpos = stds, stdneg = -stds)
  ggpy = ggplot(data=ggd,aes(x=time)) + theme_bw() + xlab("Time") + ylab("Observation") +
    geom_ribbon(aes(ymax=stdpos,  ymin=stdneg),col="gray",fill="gray",alpha=0.3)+
    geom_line(aes(y=y),col="black", linewidth=0.5) +
    ggtitle("(a) Simulated data", subtitle=paste0("Memory breakpoint at t = ",bp))
  ggpy
  
  ggp = ggplot(data=ggd,aes(x=time)) + theme_bw() + xlab("Time") + ylab(expression(paste(phi,"(t)"))) +
    geom_ribbon(aes(ymin=philower,ymax=phiupper),col="red",fill="red",alpha=0.3) +
    geom_line(aes(y=phis),col="black", linewidth=0.8) +
    geom_line(aes(y=phimean),col="blue",linewidth=0.8) +
    ylim(c(0,1)) +
    ggtitle("(b) Memory evolution", subtitle=paste0("Breakpoint at t = ",bp))
  ggp
  
  library(ggpubr)
  ggboth = ggarrange(ggpy,ggp,nrow=1)
  ggsave("breakpoint-plot-8x4.eps",plot=ggboth, device=cairo_ps, width=8,
         height=4)
  
}