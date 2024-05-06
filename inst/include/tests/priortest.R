if(FALSE){
  nsimmed=10000
  n=1000
  bsimmed = runif(nsimmed,-1,1)
  #simmed2 = exp(bsimmed)/(1+exp(bsimmed))^2
  bdens = function(x) exp(-x)/(1+exp(-x))^2
  xx = seq(from=-10,to=10,length.out=n)
  plot(xx,bdens(xx),type="l")
  lines(xx,dnorm(xx,sd=1))
  #simmed2 = -1 + 2/(1+exp(-bsimmed))
  simmed2 = log((1+bsimmed)/(1-bsimmed))
  hist(simmed2,breaks=100)
  
  b = 0.25
  L = -min(b,0)
  U = 1-max(b,0)
  adens = function(x, L, U) exp(-x)/(1+exp(-x))^2
  plot(xx,adens(xx,L,U),type="l")
  xxx = seq(0,30,length.out=n)
  plot(xxx,dgamma(xxx,shape=1,rate=0.1))
  
  mylp <- function(theta) {
    #dgamma(exp(theta[1]), 1, 5e-05, log=TRUE) + theta[1] + #vague
    lprior = dgamma(exp(theta[1]), shape=1, rate=0.1) + theta[1]
    lprior = lprior -theta[2] -2*log(1+exp(theta[2]))
    lprior = lprior -theta[3] -2*log(1+exp(-theta[3]))
      #dnorm(theta[2], sd=0.1, log=TRUE) +
      #dnorm(theta[3], sd=0.1, log=TRUE))
    return(lprior)
  }
  
  simmed1 = rnorm(nsimmed, sd=10)
  simmed1 = rgamma(nsimmed, 1, 5e-05)
  simmed2 = rnorm(nsimmed, sd=1)
  simmed3 = rnorm(nsimmed, sd=1)
  
  ssimmed = exp(simmed1)
  bsimmed = -1+2/(1+exp(-simmed2))
  asimmed = numeric(nsimmed)
  aLs = numeric(nsimmed)
  aUs = numeric(nsimmed)
  for(s in 1:nsimmed){
    aL = -min(bsimmed[s],0)
    aU = 1-max(bsimmed[s],0)
    aLs[s]=aL
    aUs[s]=aU
    asimmed[s] = aL + (aU-aL)/(1+exp(-simmed[3]))
  }
  hist(ssimmed)
  hist(bsimmed,xlim=c(-1,1))
  hist(asimmed,xlim=c(0,1))
  ggd = data.frame(s=ssimmed,b=bsimmed,a=asimmed)
  ggprior = ggplot(ggd)
  ggps = ggprior + geom_histogram(aes(s))
  ggpb = ggprior + geom_histogram(aes(b))
  ggpa = ggprior + geom_histogram(aes(a))
  
  ggps
  ggpb
  ggpa
  
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
  
  n = 1000
  nsims = 1000
  
  ttime = 1:n
  time_norm = (ttime-min(ttime))/(max(ttime)-min(ttime))
  
  
  lps = c(lp111, lp222, lp333) #same prior for all thetas
  lps = c(mylp)
  
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
          rgm = inla.rgeneric.define(rgeneric.ar1, n=n, time=time_norm, my.log.prior=lps[[sigmaiter]])
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
        
        bests[siter, iter] = b_est
        bposs[siter, iter] = bpos
        
        }
    
    }
    ggd = data.frame(b = factor(rep(bgrid, each=nsims)), best = c(t(bests)), bpositive = c(t(bposs)))
    #write.table(ggd,paste0("lpthetab,sigmas[sigmaiter],"_n500_bfits.txt"))
    bestest[[sigmaiter]] = bests
    bposes[[sigmaiter]] = bposs
  }
 
  #write.table(bposs, "bposs_n1000_unif.txt")
  #write.table(bests, "bestest_n1000_unif.txt")
  
  best1 = t(read.table("bestest_n500_unif.txt"))
  bpo1 = t(read.table("bposs_n500_unif.txt"))
  
  best2 = t(read.table("bestest_n1000_unif.txt"))
  bpo2 = t(read.table("bposs_n1000_unif.txt"))
  
  #
  #write.table(ggd,paste0("lpthetab,sigmas[sigmaiter],"_n500_bfits.txt"))
  sigmaiter=1
  bests = bestest[[sigmaiter]]
  bposs = bposes[[sigmaiter]]
  ggd1 = data.frame(b = factor(rep(bgrid, each=nsims)), best = c(t(best1)), bpositive = c(t(bpo1)))
  library(ggplot2)
  library(litteR)
  library(ggpubr)
  gg1 = ggplot(data=ggd1,aes(x=b,y=best)) + theme_bw() + xlab("True b") + ylab("Posterior marginal mean b") +
    #geom_boxplot(aes(x=b,y=best)) +
    stat_adj_boxplot() +
    theme(text=element_text(size=16), plot.title = element_text(size=22)) +
    stat_adj_boxplot_outlier()+
    geom_line(data=data.frame(bx=factor(seq(from=-0.8,to=0.8,by=0.1)),bgrid=seq(from=-0.8,to=0.8,by=0.1)),aes(x=factor(bx),y=bgrid, group=1), 
              linewidth=1.2, col="blue") +
    
    ggtitle("(a) Posterior marginal mean", subtitle=paste0("n = ",500))
  
  gg2 = ggplot(data=ggd1,aes(x=b,y=bpositive)) + theme_bw() + xlab("True b") + ylab("P(b>0)") +
    #geom_boxplot() #+
    stat_adj_boxplot() +
    theme(text=element_text(size=16), plot.title = element_text(size=22)) +
    stat_adj_boxplot_outlier()+
    ggtitle("(c) Probability of positive trend", subtitle=paste0("n = ",500))
  
  ggd2 = data.frame(b = factor(rep(bgrid, each=nsims)), best = c(t(best2)), bpositive = c(t(bpo2)))
  
  gg3 = ggplot(data=ggd2,aes(x=b,y=best)) + theme_bw() + xlab("True b") + ylab("Posterior marginal mean b") +
    #geom_boxplot(aes(x=b,y=best)) +
    stat_adj_boxplot() +
    theme(text=element_text(size=16), plot.title = element_text(size=22)) +
    stat_adj_boxplot_outlier()+
    geom_line(data=data.frame(bx=factor(seq(from=-0.8,to=0.8,by=0.1)),bgrid=seq(from=-0.8,to=0.8,by=0.1)),aes(x=factor(bx),y=bgrid, group=1), 
              linewidth=1.2, col="blue") +
    
    ggtitle("(b) Posterior marginal mean", subtitle=paste0("n = ",1000))
  
  gg4 = ggplot(data=ggd2,aes(x=b,y=bpositive)) + theme_bw() + xlab("True b") + ylab("P(b>0)") +
    #geom_boxplot() #+
    stat_adj_boxplot() +
    theme(text=element_text(size=16), plot.title = element_text(size=22)) +
    stat_adj_boxplot_outlier()+
    ggtitle("(d) Probability of positive trend", subtitle=paste0("n = ",1000))
  
  ggall = ggarrange(ggarrange(gg1,gg3),ggarrange(gg2,gg4),nrow=2)
  #ggboth = ggarrange(gg3,gg4)
  #print(ggboth)
  print(ggall)
  ggsave("accuracytest-unif-sigma1.eps",plot=ggall, device=cairo_ps, width=18, 
         height=12, limitsize=FALSE)
  
  btrue = bgrid
  bmean500 = rowMeans(best1)
  bmean1000 = rowMeans(best2)
  bpos500 = rowMeans(bpo1)
  bpos1000 = rowMeans(bpo2)
  positives500 = rowSums(bpo1>0.95)
  positives1000 = rowSums(bpo2>0.95)
  negatives500 = rowSums(bpo1<0.95)
  negatives1000 = rowSums(bpo2<0.95)
  
  rstr1 = paste0("True b")
  rstr2 = paste0("Mean b (n=500)")
  rstr3 = paste0("Mean b (n=1000)")
  rstr4 = paste0("P(b>0) (n=500)")
  rstr5 = paste0("P(b>0) (n=1000)")
  DIGITS=3
  #topstr = paste0("True b")
  #str = paste0(rstr1, " & ",rstr2, " & ",rstr3, " & ",rstr4, " & ", rstr5, "\\\\ \\hline \n ")
  str = ""
  for(i in 1:17){
    str = paste0(str, round(btrue[i], digits=DIGITS), " & ", round(bmean500[i], digits=DIGITS), 
                 " & ", round(bmean1000[i], digits=DIGITS), " & ", round(bpos500[i], digits=DIGITS), 
                 " & ",round(bpos1000[i], digits=DIGITS), " & ",round(positives500[i], digits=DIGITS), 
                 " & ", round(positives1000[i], digits=DIGITS), " \\\\ \n")
  }
  cat(str)
  
  
  best0 = t(read.table("bestest_n500_unif.txt"))
  bpos0 = t(read.table("bposs_n500_unif.txt"))
  
  #best0 = t(read.table("bestest_n1000_unif.txt"))
  #bpos0 = t(read.table("bposs_n1000_unif.txt"))
  
  Prior_bestest = readRDS("Prior_bestest.RData")
  Prior_bposes = readRDS("Prior_bposes.RData")
  bgrid=seq(-0.8,0.8,by=0.1)
  nsims=1000
  
  
  bMeans = matrix(NA,17,4)
  bLower = matrix(NA,17,4)
  bUpper = matrix(NA,17,4)
  pMeans = matrix(NA,17,4)
  pLower = matrix(NA,17,4)
  pUpper = matrix(NA,17,4)
  
  bMeans[,1] = rowMeans(best0)
  pMeans[,1] = rowMeans(bpos0)
  
  blist = list()
  blist[[1]] = best0
  blist[[2]] = Prior_bestest[[1]]
  blist[[3]] = Prior_bestest[[2]]
  blist[[4]] = Prior_bestest[[3]]
  plist = list()
  plist[[1]] = bpos0
  plist[[2]] = Prior_bposes[[1]]
  plist[[3]] = Prior_bposes[[2]]
  plist[[4]] = Prior_bposes[[3]]
  library(litteR)
  ggdtest = data.frame(b = factor(rep(bgrid, each=nsims)), 
                       bpositive = c(t(plist[[4]])))
  ggplot(data=ggdtest,aes(x=b,y=bpositive)) + theme_bw() + xlab("True b") + ylab("P(b>0)") +
    #geom_boxplot() #+
    stat_adj_boxplot() +
    theme(text=element_text(size=16), plot.title = element_text(size=22)) +
    stat_adj_boxplot_outlier()
  
  
  for(i in 1:17){
    for(j in 1:4){
      bb = blist[[j]]
      pp = plist[[j]]
      bLower[i,j] = quantile(bb[i,],probs=0.025)
      pLower[i,j] = quantile(pp[i,],probs=0.025)
      bUpper[i,j] = quantile(bb[i,],probs=0.975)
      pUpper[i,j] = quantile(pp[i,],probs=0.975)
      bMeans[i,j] = mean(bb[i,])
      pMeans[i,j] = mean(pp[i,])
    }
    
  }
  bseq = seq(-0.8,0.8,by=0.1)
  ggb = data.frame(b=bseq, bMeans=bMeans-bseq, bLower=bLower-bseq, bUpper=bUpper-bseq, bseq2=numeric(length(bseq)))
  ggb = data.frame(b=bseq, bMeans=bMeans, bLower=bLower, bUpper=bUpper, bseq2=bseq)
  
  library(ggplot2)
  priorstring = c("True b","Default prior", "N(0, 0.01)", "N(0, 1)", "N(0, 100)")
  colors = c("blue","black", "red", "orange", "yellow2")
  colorvals = c("True b" = "blue","Default prior"="black", "N(0, 0.01)" = "red", 
                "N(0, 1)" = "orange", "N(0, 100)" = "yellow2")
  ggpb = ggplot(data=ggb,aes(x=bseq)) + theme_bw() + xlab("b") + 
    #ylab(expression(paste(hat(b)," - b"))) +
    ylab("Posterior marginal mean b - true value") +
    theme(text=element_text(size=16), plot.title = element_text(size=22)) +
    geom_ribbon(aes(ymin=bLower.1,ymax=bUpper.1),fill=colors[2], alpha=0.3) +
    geom_line(aes(y=.data[["bMeans.1"]], col=priorstring[2])) +
    scale_shape_manual(name="Priors", labels = priorstring) +
    scale_color_manual(name="Priors", values = colorvals) +
    ggtitle("(a) Prior comparison",subtitle="Posterior marginal mean")
  
  ggpb <- ggpb + geom_line(aes(y=bMeans.2, col=priorstring[3]))
  ggpb <- ggpb + geom_line(aes(y=bMeans.3, col=priorstring[4]))
  ggpb <- ggpb + geom_line(aes(y=bMeans.4, col=priorstring[5]))
  
  #ggpb
  ggpb <- ggpb + geom_line(aes(y=bseq2, col=priorstring[1]))
  #
  ggp = data.frame(b=bseq, pMeans=pMeans, pLower=pLower, pUpper=pUpper, bseq2=bseq)
  ggpp = ggplot(data=ggp,aes(x=bseq)) + theme_bw() + xlab("b") + 
    #ylab(expression(paste(hat(b)," - b"))) +
    ylab("P(b>0)") +
    geom_ribbon(aes(ymin=pLower.1,ymax=pUpper.1),fill=colors[2], alpha=0.3) +
    geom_line(aes(y=.data[["pMeans.1"]], col=priorstring[2])) +
    scale_shape_manual(name="Priors", labels = priorstring) +
    scale_color_manual(name="Priors", values = colorvals)+
    theme(text=element_text(size=16), plot.title = element_text(size=22)) +
    #scale_fill_manual(name="95% CI")
    ggtitle("(b) Prior comparison",subtitle="Probability of positive trend")
  
  ggpp <- ggpp + geom_line(aes(y=pMeans.2, col=priorstring[3]))#+
    #geom_ribbon(aes(ymin=pLower.2,ymax=pUpper.2),fill=colors[3], alpha=0.3) 
  
  ggpp <- ggpp + geom_line(aes(y=pMeans.3, col=priorstring[4]))
  ggpp <- ggpp + geom_line(aes(y=pMeans.4, col=priorstring[5]))
  
  ggpp
  library(ggpubr)
  ggpbp = ggarrange(ggpb,ggpp)
  ggpbp
  ggsave("accuracytest-priorcomparison.eps",plot=ggpbp, device=cairo_ps, width=18, 
         height=6.5, limitsize=FALSE)
  
  falsepos = numeric(4)
  falseneg = numeric(4)
  for(i in 1:4){
    pp = plist[[i]]
    #bgrid[10]=0.1
    falseneg[i] = sum(pp[10:17,]<0.95)
    falsepos[i] = sum(pp[1:9,]>0.95)
  }
  
  for(i in 1:16){
    rstr1 = paste0(rstr1, " & ", round(btrue[i], digits=DIGITS))
    rstr2 = paste0(rstr2, " & ", round(bmean500[i], digits=DIGITS))
    rstr3 = paste0(rstr3, " & ", round(bmean1000[i], digits=DIGITS))
    rstr4 = paste0(rstr4, " & ", round(bpos500[i], digits=DIGITS))
    rstr5 = paste0(rstr5, " & ", round(bpos1000[i], digits=DIGITS))
  }
  rstr1 = paste0(rstr1, " & ", round(btrue[i], digits=DIGITS), "\\\\")
  rstr2 = paste0(rstr2, " & ", round(bmean500[i], digits=DIGITS), "\\\\")
  rstr3 = paste0(rstr3, " & ", round(bmean1000[i], digits=DIGITS), "\\\\")
  rstr4 = paste0(rstr4, " & ", round(bpos500[i], digits=DIGITS), "\\\\")
  rstr5 = paste0(rstr5, " & ", round(bpos1000[i], digits=DIGITS), "\\\\")
  
  cat(rstr1,rstr2,rstr3,rstr4,rstr5,sep="\n")
  
}
