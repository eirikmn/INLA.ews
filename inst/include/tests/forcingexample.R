if(FALSE){
  from.theta = function(theta, time_norm, do.forcing=FALSE){
    n=length(time_norm)
    if(do.forcing){
      pars = numeric(5)
      
      kappa_f = exp(theta[4])
      kappafs = kappa_f*2*lambdas
      F0 = theta[5]
      pars[4]=1/sqrt(kappa_f)
      pars[5]=F0
    }else{
      pars = numeric(3)
    }
    
    kappa_eps = exp(theta[1])
    r = diff(range(time_norm))
    b = -1/r+2/r*1/(1+exp(-theta[2]))
    
    low = min(b*time_norm[1],b*time_norm[n])
    high = max(b*time_norm[1],b*time_norm[n])
    
    a = -low + (1-high+low)/(1+exp(-theta[3]))
    
    pars[1]=1/sqrt(kappa_eps)
    pars[2]=b
    pars[3]=a
    lambdas = -log(a+b*time_norm)
    kappas = kappa_eps*2*lambdas
    
    cc = 1/(n-1)
    phis = c(exp(-lambdas*c(1,diff(time_norm)/cc)) ) #rescale
    if(do.forcing){
      return(list(pars=pars, phis = phis, kappa_eps = kappa_eps, a=a,b=b,#kappa_f=kappa_f,F0=F0
                  lambdas=lambdas, kappas=kappas, sigmas=1/kappas,
                  kappa_f=kappa_f, kappafs=kappafs, sigmafs=1/sqrt(kappafs), 
                  F0=F0))
    }else{
      return(list(pars=pars, phis = phis, kappa_eps = kappa_eps, a=a,b=b,#kappa_f=kappa_f,F0=F0
                  lambdas=lambdas, kappas=kappas, sigmas=1/kappas))
    }
    
  }
  
  to.theta = function(pars, time_norm, do.forcing=FALSE){
    n=length(time_norm)
    if(do.forcing){
      theta = numeric(5)
      
      theta[4] = log(1/pars[4]^2)
      theta[5] = pars[5]
    }else{
      theta = numeric(3)
    }
    kappa_eps = 1/pars[1]^2
    a=pars[3]
    b=pars[2]
    
    theta[1] = log(1/pars[1]^2)
    
    r = diff(range(time_norm))
    theta[2] = log((1+pars[2])/(1-pars[2]))
    low = min(pars[2]*time_norm[1],pars[2]*time_norm[n])
    high = max(pars[2]*time_norm[1],pars[2]*time_norm[n])
    theta[3] = log((pars[3]+low)/(1-high-pars[3]))
    
    lambdas = -log(pars[3]+pars[2]*time_norm)
    kappas = kappa_eps*2*lambdas
    
    cc = 1/(n-1)
    phis = c(exp(-lambdas*c(1,diff(time_norm)/cc)) ) #rescale
    
    kappa_f = exp(theta[4])
    kappafs = kappa_f*2*lambdas
    F0 = theta[5]
    
    if(do.forcing){
      return(list(theta=theta, phis = phis, kappa_eps = kappa_eps, a=a,b=b,#kappa_f=kappa_f,F0=F0
                  lambdas=lambdas, kappas=kappas, sigmas=1/kappas,
                  kappa_f=kappa_f, kappafs=kappafs, sigmafs=1/sqrt(kappafs), 
                  F0=F0))
    }else{
      return(list(theta=theta, phis = phis, kappa_eps = kappa_eps, a=a,b=b,#kappa_f=kappa_f,F0=F0
                  lambdas=lambdas, kappas=kappas, sigmas=1/kappas))
    }
  }
  mucomputer = function(pars,forcing,time_norm,as.theta=FALSE){
    if(as.theta){
      allinfo = from.theta(pars,time_norm,do.forcing=TRUE)
      parstheta = pars
      parsttemp = allinfo$pars
    }else{
      parstemp=pars
      allinfo = to.theta(pars,time_norm,do.forcing=TRUE)
      parstheta = allinfo$theta
    }
    a=allinfo$a; b=allinfo$b
    sigmafs=allinfo$sigmafs
    F0=allinfo$F0
    zz = sigmafs*(F0+forcing)
    lambdas = -log(a+b*time_norm)
    muvek = numeric(n)
    if(diff(range(diff(time_norm)))==0){
      struktur = exp(-lambdas*time_norm)
      for(i in 1:n){
        muvek[i] = rev(struktur[1:i])%*%zz[1:i]
      }
    }else{
      for(k in 1:n){
        for(s in 1:k){
          muvek[k] = muvek[k] + zz[s]*exp(-lambdas[k]*(time_norm[k]-time_norm[s]))
        }
      }
    }
    
    return(muvek)
  }
}

if(FALSE){
  
  ##### START HERE
  set.seed(123)
  n=1000
  a = 0.3
  b = 0.2
  
  ttime = sort(1:n + rnorm(n,sd=0.1))
  #ttime = 1:n
  time_norm = (ttime-min(ttime))/(max(ttime)-min(ttime))
  lambdas = -log(a+b*time_norm)
  cc=1/(n-1)
  phis = exp(-lambdas*c(1,diff(time_norm)/cc))
  
  plot(phis)
  lines(exp(-lambdas),col="red",lwd=3)
  my.log.prior <- function(theta) {
    lprior = dgamma(exp(theta[1]), shape=1, rate=0.1) + theta[1]
    lprior = lprior -theta[2] -2*log(1+exp(-theta[2]))
    lprior = lprior -theta[3] -2*log(1+exp(-theta[3])) 
    lprior = lprior + dgamma(exp(theta[4]), shape=1, rate=0.1) + theta[4]
    lprior = lprior + dnorm(theta[5], sd=10, log=TRUE)
    return(lprior)
  }
  
  # my.log.prior <- function(theta) {return(
  #   dnorm(theta[1], sd=10, log=TRUE) +
  #     dnorm(theta[2], sd=10, log=TRUE) +
  #     dnorm(theta[3], sd=10, log=TRUE) +
  #     dnorm(theta[4], sd=10, log=TRUE) +
  #     dnorm(theta[5], sd=10, log=TRUE)
  # )
  # }
  
  forcing = arima.sim(list(ar=c(0.95)), n=n)
  
  #phis=rep(phi,n)
  F0=0
  sigma=5
  sigma_f = 0.1
  
  
  sims=numeric(n)
  sims[1] = rnorm(1,sd=sigma/sqrt(2*lambdas[1]))
  for(i in 2:n){
    sims[i] = phis[i]*sims[i-1] + rnorm(1,sd=sigma/sqrt(lambdas[i]))
  }
  #sims = sigma * arima.sim(list(ar=c(phi)), n=n, sd=sqrt(1-phi^2))
  muvek = mucomputer(pars=c(sigma,b,a,sigma_f,F0),forcing,time_norm,as.theta=FALSE)
  
  intercept=0
  y = sims + muvek + intercept
  plot(y)
  lines(muvek,col="red",lwd=5)
  
  
  #control.family=list(initial=12,fixed=TRUE)
  #inla.options = list(control.family=control.family)
  res = inla.ews(data=y,forcing=forcing,formula=y ~ -1, log.prior=my.log.prior,
                 timesteps=ttime,print.progress = TRUE, nsims=10000,
                 #stepsize=c(0.01,0.001)
                 inla.options = list(verbose=TRUE),
                 stepsize=c(0.005)
  )
  summary(res)
  #res = inla.ews(data=y,forcing=z, timesteps=time,formula=y ~ -1)
  
  #ggsave("forcedsims-19200x7200.eps",plot=ggy0, device=cairo_ps, width=19200,
  #       height=7200, units="px", dpi=1800)
  lines(res$results$summary$alltrend$mean,col="blue",lwd=5)
  #library(ggplot2)
  
  library(ggplot2)
  ggy0 = ggplot(data=data.frame(x=ttime,y=y,mu=muvek)) + theme_bw()+
    theme(text=element_text(size=18), plot.title = element_text(size=24)) + 
    xlab("Time (yr)") + ylab("Simulated data") +
    ggtitle("(a) Simulated forced data") + xlim(range(ttime)) +
    geom_line(aes(x=x,y=y),col="gray", linewidth=0.8) +
    geom_line(aes(x=x,y=mu), col="black", linewidth=1.2)
  
  print(ggy0)
  
  
  ggm = ggplot(data=as.data.frame(res$results$summary$alltrend),aes(x=ttime)) + theme_bw() + 
    xlab("Time (yr)") + ylab("Observation") +
    theme(text=element_text(size=18), plot.title = element_text(size=24)) + 
    xlim(c(ttime[1],ttime[n])) +
    #geom_point(aes(y=y), size=2,col="gray") +
    geom_line(aes(y=y),col="gray", linewidth=0.8) +
    geom_ribbon(aes(ymin=quant0.025,ymax=quant0.975),fill="red",col="red",alpha=0.3,linewidth=1.2) +
    # geom_line(data=data.frame(time=ttime,my=muvek),
    #            aes(x=time,y=my),col="white",linewidth=1.2) +
    geom_line(aes(y=mean),col="blue",linewidth=1.2) +
    ggtitle("(b) Estimated forcing response")
  print(ggm)
  
  ggp = ggplot(data=as.data.frame(res$results$summary$phi0),aes(x=ttime)) +
    theme_bw() + xlim(range(ttime))+
    theme(text=element_text(size=18), plot.title = element_text(size=24)) + 
    ylim(range(0,1,res$results$summary$phi0$q0.025[2:n],res$results$summary$phi0$q0.975[2:n])) +
    geom_ribbon(aes(ymin=q0.025,ymax=q0.975),fill="red",col="red",alpha=0.3,linewidth=0)+
    geom_line(data=data.frame(ttime=ttime[2:n],py = res$results$summary$phi$mean[2:n]),
              aes(x=ttime,y=py),col="gray") +
    geom_line(data=data.frame(time=ttime,py=a+b*time_norm),
              mapping=aes(x=time,y=py),col="black",linewidth=1.2)+
    geom_line(aes(y=mean),col="blue",linewidth=1.2) +
    geom_line(aes(y=q0.025),col="red",linewidth=1.2) +
    geom_line(aes(y=q0.975),col="red",linewidth=1.2) +
    xlab("Time (yr)") + ylab(expression(paste(phi,"(t)"))) +
    ggtitle("(c) Memory evolution")
  print(ggp)
  
  library(ggpubr)
  ggboth = ggarrange(ggarrange(ggy0,ggm,ncol=2),ggp,nrow=2, ncol=1)
  print(ggboth)
  
  ggsave("forcingfit-both-14x10.eps",plot=ggboth, device=cairo_ps, width=14,
         height=10)
  ## hvorfor blir forcing response annerledes?
  
  cat("Estimates (Credible interval) [True value]:\n\ta = ",
      round(res$results$summary$a$mean,digits=3)," (",
      round(res$results$summary$a$quant0.025,digits=3),", ",
      round(res$results$summary$a$quant0.975,digits=3),") [",a,"]\n\tb = ",
      round(res$results$summary$b$mean,digits=3)," (",
      round(res$results$summary$b$quant0.025,digits=3),", ",
      round(res$results$summary$b$quant0.975,digits=3),") [",b,"]\n\tsigma_x = ",
      round(res$results$summary$sigma$mean,digits=3)," (",
      round(res$results$summary$sigma$quant0.025,digits=3),", ",
      round(res$results$summary$sigma$quant0.975,digits=3),") [",sigma,"]\n\tsigma_f = ",
      round(res$results$summary$sigmaf$mean,digits=3)," (",
      round(res$results$summary$sigmaf$quant0.025,digits=3),", ",
      round(res$results$summary$sigmaf$quant0.975,digits=3),
      ") [",sigma_f,"]\n\tF0 = ",
      round(res$results$summary$F0$mean,digits=3)," (",
      round(res$results$summary$F0$quant0.025,digits=3),", ",
      round(res$results$summary$F0$quant0.975,digits=3),") [",F0,"]\n",sep="")
  
  
  cat("Estimates (Credible interval) [True value]:\n\t$a$ & ",a, " & ",
      round(res$results$summary$a$mean,digits=3)," & (",
      round(res$results$summary$a$quant0.025,digits=3),", ",
      round(res$results$summary$a$quant0.975,digits=3),") \\\\ \n\t$b$ &",b, "& ",
      round(res$results$summary$b$mean,digits=3)," & (",
      round(res$results$summary$b$quant0.025,digits=3),", ",
      round(res$results$summary$b$quant0.975,digits=3),") \\\\ \n\t$\\sigma$ & ",sigma," & ",
      round(res$results$summary$sigma$mean,digits=3)," & (",
      round(res$results$summary$sigma$quant0.025,digits=3),", ",
      round(res$results$summary$sigma$quant0.975,digits=3),") \\\\ \n\t$\\sigma_f$ & ",sigma_f, " & ",
      round(res$results$summary$sigmaf$mean,digits=3)," & (",
      round(res$results$summary$sigmaf$quant0.025,digits=3),", ",
      round(res$results$summary$sigmaf$quant0.975,digits=3),
      ") \\\\ \n\t$F_0$ & ",F0," & ",
      round(res$results$summary$F0$mean,digits=3)," & (",
      round(res$results$summary$F0$quant0.025,digits=3),", ",
      round(res$results$summary$F0$quant0.975,digits=3),") \\\\ \n",sep="")
  
  
  if(FALSE){
    
    set.seed(123)
    n=500
    a = 0.3
    b = 0.2
    
    ttime = sort(1:n + rnorm(n,sd=0.1))
    time_norm = (ttime-min(ttime))/(max(ttime)-min(ttime))
    lambdas = -log(a+b*time_norm)
    cc=1/(n-1)
    phis = exp(-lambdas*c(1,diff(time_norm)/cc))
    
    plot(phis)
    lines(exp(-lambdas),col="red",lwd=3)
    
    forcing = arima.sim(list(ar=c(0.95)), n=n)
    
    #phis=rep(phi,n)
    F0=0
    sigma=5
    sigma_f = 0.1
    
    
    sims=numeric(n)
    sims[1] = rnorm(1,sd=sigma/sqrt(2*lambdas[1]))
    for(i in 2:n){
      sims[i] = phis[i]*sims[i-1] + rnorm(1,sd=sigma/sqrt(lambdas[i]))
    }
    #sims = sigma * arima.sim(list(ar=c(phi)), n=n, sd=sqrt(1-phi^2))
    muvek = mucomputer(pars=c(sigma,b,a,sigma_f,F0),forcing,time_norm,as.theta=FALSE)
    
    intercept=0
    y = sims + muvek + intercept
    plot(y)
    lines(muvek,col="red",lwd=5)
    
    library(numDeriv)
    param = c(sigma, b,a,sigma_f,F0)
    tparam = to.theta(param, time_norm=time_norm,do.forcing = TRUE)
    #theta = tparam$theta[2:5]
    theta = numeric(4)
    args = list(forcing=forcing,y=y,time_norm=time_norm)
    
    minfun.grad = function(param, args = NULL){
      return (grad(minfun, param, args=args, method.args = list(r=6)))
    }
    minfun = function(param, args = NULL){
      
      
      muvek = mucomputer(c(0,param),forcing,time_norm,as.theta=TRUE)
      mse = sum((args$y-muvek)^2)
      print(param)
      print(mse)
      return(sqrt(mse))
    }
    
    ## perform optimization to find good starting values for INLA
    #param=optparams
    #args=list(y=y)
    # fit = optim(theta,
    #             fn = minfun,
    #             gr = minfun.grad,
    #             method = "BFGS",
    #             control = list(
    #               abstol = 0,
    #               maxit = 100000,
    #               reltol = 1e-11),
    #             args = args)
    
    fit = optim(par = theta, 
                fn = minfun, method = "Nelder-Mead", hessian = FALSE, args=args)
    
    ### use least squares estimates for fixed effects as initial values in inla
    
    muvekfit = mucomputer(c(0,fit$par),forcing,time_norm,as.theta=TRUE)
    #muvekfit = mucomputer(a=fit$par[2],b=fit$par[1],sigmaf=fit$par[3],F0=fit$par[4],forcing=forcing,time_norm=time_norm)
    plot(y)
    allinfo = from.theta(c(0,fit$par),time_norm,do.forcing=TRUE)
    lines(muvek)
    lines(muvekfit,col="blue",lwd=5)
    #theta_init = c(0,fit$par)
    nhyps=5
    theta_init = c(numeric(nhyps-2),tail(fit$par,2))
    
    
    res = inla.ews(data=y,forcing=forcing,formula=y ~ -1, 
                   timesteps=ttime,print.progress = TRUE, nsims=1000#,
                   #stepsize=c(0.01,0.001)
                   #stepsize=c(0.01)#,
                   #inla.options = list(control.mode=list(fixed=FALSE) )
    )
    
    
    rgen_model = inla.rgeneric.define(rgeneric.ar1.forcing, n=n,time=time_norm,forcing=forcing)
    ff = y~-1+f(idy,model=rgen_model)
    rrr = inla(ff,data=data.frame(y=y,idy=1:n),control.family=list(initial=12,fixed=TRUE),verbose=TRUE,
               control.mode=list(theta=c(0,0,0.1,0.63,0.0512),
                                 fixed=FALSE#, restart=TRUE
               )
    )
    summary(rrr)
    raa = inla(ff,data=data.frame(y=y,idy=1:n),#control.family=list(initial=12,fixed=TRUE),
               #control.mode=list(theta=theta_init)
               #num.threads=inla.options$num.threads,
               #control.mode=inla.options$control.mode,
               #control.inla=inla.options$control.inla,
               control.family=inla.options$control.family,
               #control.compute=inla.options$control.compute
    )
    summary(raa)
    aaaa = from.theta(rrr$summary.hyperpar$mean,time_norm,do.forcing=TRUE)
    resultgather()
  }
  
  
  
  
  
}
