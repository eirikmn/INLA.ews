if(FALSE){
  
  sigmas = c(0.1, 1, 10)
  #hyperhyperparams = data.frame(sigma = sigmas)
  
  str=""
  for(i in 1:length(sigmas)){
    for(j in 1:length(sigmas)){
      for(k in 1:length(sigmas)){
        sd1=sigmas[i]; sd2=sigmas[j]; sd3=sigmas[k]
        fbody = paste0("\ndnorm(theta[1], sd=",sd1,", log=TRUE) +\n",
                       " dnorm(theta[2], sd=",sd2,", log=TRUE) +\n",
                       " dnorm(theta[3], sd=",sd3,", log=TRUE)")
        str=paste0(str, '\nlp',i,j,k,' <- function(','theta',') {return(',fbody,')\n}',sep="")
      }
    }
  }
  cat(str) #run this to get function definitions (horribly lazy way to implement this, but it works..)
  
  #######
  if(TRUE){ #this should load all of them in one click
    lp111 <- function(theta) {return(
      dnorm(theta[1], sd=0.1, log=TRUE) +
        dnorm(theta[2], sd=0.1, log=TRUE) +
        dnorm(theta[3], sd=0.1, log=TRUE))
    }
    lp112 <- function(theta) {return(
      dnorm(theta[1], sd=0.1, log=TRUE) +
        dnorm(theta[2], sd=0.1, log=TRUE) +
        dnorm(theta[3], sd=1, log=TRUE))
    }
    lp113 <- function(theta) {return(
      dnorm(theta[1], sd=0.1, log=TRUE) +
        dnorm(theta[2], sd=0.1, log=TRUE) +
        dnorm(theta[3], sd=10, log=TRUE))
    }
    lp121 <- function(theta) {return(
      dnorm(theta[1], sd=0.1, log=TRUE) +
        dnorm(theta[2], sd=1, log=TRUE) +
        dnorm(theta[3], sd=0.1, log=TRUE))
    }
    lp122 <- function(theta) {return(
      dnorm(theta[1], sd=0.1, log=TRUE) +
        dnorm(theta[2], sd=1, log=TRUE) +
        dnorm(theta[3], sd=1, log=TRUE))
    }
    lp123 <- function(theta) {return(
      dnorm(theta[1], sd=0.1, log=TRUE) +
        dnorm(theta[2], sd=1, log=TRUE) +
        dnorm(theta[3], sd=10, log=TRUE))
    }
    lp131 <- function(theta) {return(
      dnorm(theta[1], sd=0.1, log=TRUE) +
        dnorm(theta[2], sd=10, log=TRUE) +
        dnorm(theta[3], sd=0.1, log=TRUE))
    }
    lp132 <- function(theta) {return(
      dnorm(theta[1], sd=0.1, log=TRUE) +
        dnorm(theta[2], sd=10, log=TRUE) +
        dnorm(theta[3], sd=1, log=TRUE))
    }
    lp133 <- function(theta) {return(
      dnorm(theta[1], sd=0.1, log=TRUE) +
        dnorm(theta[2], sd=10, log=TRUE) +
        dnorm(theta[3], sd=10, log=TRUE))
    }
    lp211 <- function(theta) {return(
      dnorm(theta[1], sd=1, log=TRUE) +
        dnorm(theta[2], sd=0.1, log=TRUE) +
        dnorm(theta[3], sd=0.1, log=TRUE))
    }
    lp212 <- function(theta) {return(
      dnorm(theta[1], sd=1, log=TRUE) +
        dnorm(theta[2], sd=0.1, log=TRUE) +
        dnorm(theta[3], sd=1, log=TRUE))
    }
    lp213 <- function(theta) {return(
      dnorm(theta[1], sd=1, log=TRUE) +
        dnorm(theta[2], sd=0.1, log=TRUE) +
        dnorm(theta[3], sd=10, log=TRUE))
    }
    lp221 <- function(theta) {return(
      dnorm(theta[1], sd=1, log=TRUE) +
        dnorm(theta[2], sd=1, log=TRUE) +
        dnorm(theta[3], sd=0.1, log=TRUE))
    }
    lp222 <- function(theta) {return(
      dnorm(theta[1], sd=1, log=TRUE) +
        dnorm(theta[2], sd=1, log=TRUE) +
        dnorm(theta[3], sd=1, log=TRUE))
    }
    lp223 <- function(theta) {return(
      dnorm(theta[1], sd=1, log=TRUE) +
        dnorm(theta[2], sd=1, log=TRUE) +
        dnorm(theta[3], sd=10, log=TRUE))
    }
    lp231 <- function(theta) {return(
      dnorm(theta[1], sd=1, log=TRUE) +
        dnorm(theta[2], sd=10, log=TRUE) +
        dnorm(theta[3], sd=0.1, log=TRUE))
    }
    lp232 <- function(theta) {return(
      dnorm(theta[1], sd=1, log=TRUE) +
        dnorm(theta[2], sd=10, log=TRUE) +
        dnorm(theta[3], sd=1, log=TRUE))
    }
    lp233 <- function(theta) {return(
      dnorm(theta[1], sd=1, log=TRUE) +
        dnorm(theta[2], sd=10, log=TRUE) +
        dnorm(theta[3], sd=10, log=TRUE))
    }
    lp311 <- function(theta) {return(
      dnorm(theta[1], sd=10, log=TRUE) +
        dnorm(theta[2], sd=0.1, log=TRUE) +
        dnorm(theta[3], sd=0.1, log=TRUE))
    }
    lp312 <- function(theta) {return(
      dnorm(theta[1], sd=10, log=TRUE) +
        dnorm(theta[2], sd=0.1, log=TRUE) +
        dnorm(theta[3], sd=1, log=TRUE))
    }
    lp313 <- function(theta) {return(
      dnorm(theta[1], sd=10, log=TRUE) +
        dnorm(theta[2], sd=0.1, log=TRUE) +
        dnorm(theta[3], sd=10, log=TRUE))
    }
    lp321 <- function(theta) {return(
      dnorm(theta[1], sd=10, log=TRUE) +
        dnorm(theta[2], sd=1, log=TRUE) +
        dnorm(theta[3], sd=0.1, log=TRUE))
    }
    lp322 <- function(theta) {return(
      dnorm(theta[1], sd=10, log=TRUE) +
        dnorm(theta[2], sd=1, log=TRUE) +
        dnorm(theta[3], sd=1, log=TRUE))
    }
    lp323 <- function(theta) {return(
      dnorm(theta[1], sd=10, log=TRUE) +
        dnorm(theta[2], sd=1, log=TRUE) +
        dnorm(theta[3], sd=10, log=TRUE))
    }
    lp331 <- function(theta) {return(
      dnorm(theta[1], sd=10, log=TRUE) +
        dnorm(theta[2], sd=10, log=TRUE) +
        dnorm(theta[3], sd=0.1, log=TRUE))
    }
    lp332 <- function(theta) {return(
      dnorm(theta[1], sd=10, log=TRUE) +
        dnorm(theta[2], sd=10, log=TRUE) +
        dnorm(theta[3], sd=1, log=TRUE))
    }
    lp333 <- function(theta) {return(
      dnorm(theta[1], sd=10, log=TRUE) +
        dnorm(theta[2], sd=10, log=TRUE) +
        dnorm(theta[3], sd=10, log=TRUE))
    }
  }
  
  
  ########## Which data to test it on? ###########
  
  n = 500
  nsims = 1000
  
  ttime = 1:n
  time_norm = (ttime-min(ttime))/(max(ttime)-min(ttime))
  
  
  lps = c(lp111, lp222, lp333) #same prior for all thetas
  
  bgrid = seq(from=-0.8, to = 0.8, by=0.1)
  #agrid = seq(from=0.1,to=0.9,by=0.1)
  bests1 = matrix(NA, nrow=length(bgrid), ncol=nsims)
  bests2 = matrix(NA, nrow=length(bgrid), ncol=nsims)
  bests3 = matrix(NA, nrow=length(bgrid), ncol=nsims)
  bposs1 = matrix(NA, nrow=length(bgrid), ncol=nsims)
  bposs2 = matrix(NA, nrow=length(bgrid), ncol=nsims)
  bposs3 = matrix(NA, nrow=length(bgrid), ncol=nsims)
  
  bestest = list(bests1=bests1, bests2=bests2, bests3=bests3)
  bposes = list(bposs1=bposs1, bposs2=bposs2, bposs3=bposs3)
  library(INLA)
  for(sigmaiter in 1:length(sigmas)){
    
    for( iter in 1:length(bgrid)){
      cat("\n\n  SIGMAITER=",sigmaiter,": ITERATION #",iter, " of ",length(bgrid)," \n\n",sep="")
      b= bgrid[iter]
      low=min(0,b)
      high=max(0,b)
      a = -low + (1-high+low)/(1+exp(-0))
      
      lambdas = -log(a+b*time_norm)
      cc=1/(n-1)
      phis = exp(-lambdas*c(1,diff(time_norm)/cc))
      sigma = 1
      
      for(siter in 1:nsims){
        cat("Simulation ",siter," / ",nsims," (iter ",iter,"/",length(bgrid),")","\n",sep="")
        sims=numeric(n)
        sims[1] = rnorm(1,sd=sigma/sqrt(2*lambdas[1]))
        for(i in 2:n){
          sims[i] = phis[i]*sims[i-1] + rnorm(1,sd=sigma/(sqrt(2*lambdas[i])))
        }
        
        
        intercept=0
        y = sims #+ muvek + intercept
        #
        if(siter==1){
          rgm = inla.rgeneric.define(rgeneric.ar1, n=n, time=time_norm, my.log.prior=lps[[sigmaiter]])
          rrm = inla(formula=y~-1 + f(idy,model=rgm), data=data.frame(y=y,idy=1:n), 
                     control.family=list(initial=12,fixed=TRUE), verbose=FALSE,
                     control.mode=list(restart=FALSE),
                     control.inla=list(h=0.005)
          )
        }else{
          rgm = inla.rgeneric.define(rgeneric.ar1, n=n, time=time_norm)
          rrm = inla(formula=y~-1 + f(idy,model=rgm), data=data.frame(y=y,idy=1:n), 
                     control.family=list(initial=12,fixed=TRUE), verbose=FALSE,
                     control.mode=list(restart=TRUE, result=rrm),
                     control.inla=list(h=0.005)
          )
        }
        
        
        rekke = diff(range(time_norm))
        b_est=INLA::inla.emarginal(function(x) -1/rekke+2/rekke*1/(1+exp(-x)) ,
                                   rrm$marginals.hyperpar$`Theta2 for idy` )
        
        bpos = 1-inla.pmarginal(0,rrm$marginals.hyperpar$`Theta2 for idy`)
        
        bests[iter, siter] = b_est
        bposs[iter, siter] = bpos
        
        }
    
    }
    ggd = data.frame(b = factor(rep(bgrid, each=nsims)), best = c(t(bests)), bpositive = c(t(bposs)))
    #write.table(ggd,paste0("lpthetab,sigmas[sigmaiter],"_n500_bfits.txt"))
    bestest[[sigmaiter]] = bests
    bposes[[sigmaiter]] = bposs
  }
  
}

sigmaiter=1
bests = bestest[[sigmaiter]]
bposs = bposes[[sigmaiter]]
ggd1 = data.frame(b = factor(rep(bgrid, each=nsims)), best = c(t(bests)), bpositive = c(t(bposs)))
library(ggplot2)
library(litteR)
gg1 = ggplot(data=ggd1,aes(x=b,y=best)) + theme_bw() + xlab("True b") + ylab("Posterior marginal mean b") +
  #geom_boxplot(aes(x=b,y=best)) +
  stat_adj_boxplot() +
  stat_adj_boxplot_outlier()+
  geom_line(data=data.frame(bx=factor(seq(from=-0.8,to=0.8,by=0.1)),bgrid=seq(from=-0.8,to=0.8,by=0.1)),aes(x=factor(bx),y=bgrid, group=1), 
            linewidth=1.2, col="blue") +
  
  ggtitle("(a) Posterior marginal mean", subtitle=paste0("n = ",500))

gg2 = ggplot(data=ggd1,aes(x=b,y=bpositive)) + theme_bw() + xlab("True b") + ylab("P(b>0)") +
  #geom_boxplot() #+
  stat_adj_boxplot() +
  stat_adj_boxplot_outlier()+
  ggtitle("(c) Probability of positive trend", subtitle=paste0("n = ",500))

ggboth = ggarrange(gg1,gg2)
ggboth
#ggsave(paste0("accuracytest-n",n,"-sigma",sigma,"-28800x9600.eps"),plot=ggboth, device=cairo_ps, width=28800,
#       height=9600, units="px", dpi=1800, limitsize=FALSE)