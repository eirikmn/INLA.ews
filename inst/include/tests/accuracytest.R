if(FALSE){
  library(INLA)
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
        lprior = dgamma(exp(theta[1]), shape=1, rate=0.1, log=TRUE) + theta[1]
        lprior = lprior -theta[2] -2*log(1+exp(-theta[2]))
        lprior = lprior -theta[3] -2*log(1+exp(-theta[3]))
        
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
  
  
  n = 1000
  nsims = 1000
  #psims = 2000
  ttime = 1:n
  time_norm = (ttime-min(ttime))/(max(ttime)-min(ttime))
  
  
  bgrid = seq(from=-0.8, to = 0.8, by=0.1)
  #agrid = seq(from=0.1,to=0.9,by=0.1)
  bests = matrix(NA, nrow=length(bgrid), ncol=nsims)
  bposs = matrix(NA, nrow=length(bgrid), ncol=nsims)
  
  for( iter in 1:length(bgrid)){
    cat("\n\n  ITERATION #",iter, " of ",length(bgrid)," \n\n",sep="")
    b= bgrid[iter]
    low=min(0,b)
    high=max(0,b)
    a = -low + (1-high+low)/(1+exp(-0))
    # a = 0.3
    # b = 0.2
    # 
    #ttime = sort(1:n + rnorm(n,sd=0.1))
    lambdas = -log(a+b*time_norm)
    cc=1/(n-1)
    phis = exp(-lambdas*c(1,diff(time_norm)/cc))
    sigma=10
    
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
        rgm = inla.rgeneric.define(rgeneric.ar1, n=n, time=time_norm)
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
    write.table(bests,paste0("tempbest-n",n,".txt"))
    write.table(bposs,paste0("tempbpos-n",n,".txt"))
    
  }
  
  #temp = read.table("tempbest-n1000.txt")
  nsims=1000
  n=500
  #bgrid = matrix(temp$b,ncol=nsims)
  bgrid=seq(-0.8,0.8,by=0.1)
  #bests=t(matrix(temp$best,ncol=nsims))
  #bposs=t(matrix(temp$bpositive,ncol=nsims))
  
  if(FALSE){
    #bests=t(matrix(temp$best,ncol=nsims))
    #bposs=t(matrix(temp$bpositive,ncol=nsims)) 
  }else{
    tempbest = read.table(paste0("tempbest-n",500,".txt"))
    tempbpos = read.table(paste0("tempbpos-n",500,".txt"))
    
    bests = as.matrix(tempbest)
    bposs = as.matrix(tempbpos)
  }
  
  
  library(ggplot2)
  library(litteR)
  ggd = data.frame(b = factor(rep(bgrid, each=nsims)), best = c(t(bests)), bpositive = c(t(bposs)))
  bmeanavg = rowMeans(bests)
  
  ggd1 = data.frame(b = factor(rep(seq(from=-0.8,to=0.8,by=0.1), each=nsims)), best = c(t(bests)), bpositive = c(t(bposs)))
  bmeanavg = rowMeans(bests)
  
  gg1 = ggplot(data=ggd1,aes(x=b,y=best)) + theme_bw() + xlab("True b") + ylab("Posterior marginal mean b") +
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
  
  
  
  if(FALSE){
    temp1best = read.table("tempbest-n500.txt") #from 0.1 to 0.8
    temp1pos = read.table("tempbpos-n500.txt") #from 0.1 to 0.8
    temp2best = read.table("tempbest-n1000.txt") #from 0.1 to 0.8
    temp2pos = read.table("tempbpos-n1000.txt") #from 0.1 to 0.8
  }else{
    temp2best = read.table(paste0("tempbest-n",1000,".txt"))
    temp2bpos = read.table(paste0("tempbpos-n",1000,".txt"))
    
    bestfull = as.matrix(temp2best)
    bposfull = as.matrix(temp2bpos)
  }
  

  
  roww = rowMeans(temp2best)
  
  
  ggd2 = data.frame(b=factor(rep(seq(from=-0.8,to=0.8,by=0.1), each=nsims)), best=c(t(bestfull)),
                    bpositive=c(t(bposfull)))
  bmeanavg2 = rowMeans(bestfull)
  
  gg3 = ggplot(data=ggd2,aes(x=b,y=best)) + theme_bw() + xlab("True b") + ylab("Posterior marginal mean b") +
    #geom_boxplot() +
    stat_adj_boxplot() +
    stat_adj_boxplot_outlier()+
    geom_line(data=data.frame(bx=factor(bgrid),bgrid=bgrid),aes(x=factor(bx),y=bgrid, group=1), 
              linewidth=1.2, col="blue") +
    ggtitle("(b) Posterior marginal mean", subtitle=paste0("n = ",1000))
  
  gg4 = ggplot(data=ggd2,aes(x=b,y=bpositive)) + theme_bw() + xlab("True b") + ylab("P(b>0)") +
    stat_adj_boxplot() +
    stat_adj_boxplot_outlier()+
    ggtitle("(d) Probability of positive trend", subtitle=paste0("n = ",1000))
  library(ggpubr)
  ggboth = ggarrange(gg1,gg3,gg2,gg4, nrow=2,ncol=2)
  ggboth
  ggsave(paste0("accuracytest-full","-28800x19200-rev.eps"),plot=ggboth, device=cairo_ps, width=28800,
         height=19200, units="px", dpi=1800, limitsize=FALSE)
  
  
  ggd1 = ggd
  gg2 = ggplot(data=ggd) + theme_bw() + xlab("True b") + ylab("Posterior marginal mean b") +
    geom_boxplot(aes(x=b,y=bpositive)) +
    theme(text=element_text(size=18), plot.title = element_text(size=24)) + 
    ggtitle("(b) Probability of positive trend", subtitle=paste0("n = ",n))
  
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
  ggsave(paste0("accuracytest-n",n,"-sigma",sigma,"-28800x9600-rev.eps"),plot=ggboth, device=cairo_ps, width=28800,
         height=9600, units="px", dpi=1800, limitsize=FALSE)
  
}




if(FALSE){
  
  set.seed(432432432)
  
  n = 1000
  nsims = 1000
  #psims = 2000
  ttime = 1:n
  time_norm = (ttime-min(ttime))/(max(ttime)-min(ttime))
  
  
  #bgrid = seq(from=-0.8, to = 0.8, by=0.1)
  bgrid = seq(from=0.8, to = 0.1, by=-0.1)
  #agrid = seq(from=0.1,to=0.9,by=0.1)
  bests2 = matrix(NA, nrow=length(bgrid), ncol=nsims)
  bposs2 = matrix(NA, nrow=length(bgrid), ncol=nsims)
  
  for( iter in 1:length(bgrid)){
    b= bgrid[iter]
    cat("\n\n  ITERATION #",iter, " of ",length(bgrid),". b = ",b," \n\n",sep="")
    low=min(0,b)
    high=max(0,b)
    a = -low + (1-high+low)/(1+exp(-0))
    # a = 0.3
    # b = 0.2
    # 
    #ttime = sort(1:n + rnorm(n,sd=0.1))
    lambdas = -log(a+b*time_norm)
    cc=1/(n-1)
    phis = exp(-lambdas*c(1,diff(time_norm)/cc))
    sigma=1
    
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
      
      rgm = inla.rgeneric.define(rgeneric.ar1, n=n, time=time_norm)
      rrm = inla(formula=y~-1 + f(idy,model=rgm), data=data.frame(y=y,idy=1:n), 
                 control.family=list(initial=12,fixed=TRUE), verbose=FALSE,
                 control.mode=list(restart=TRUE),
                 control.inla=list(h=0.005)
      )
      
      rekke = diff(range(time_norm))
      b_est=INLA::inla.emarginal(function(x) -1/rekke+2/rekke*1/(1+exp(-x)) ,
                                 rrm$marginals.hyperpar$`Theta2 for idy` )
      
      bpos = 1-inla.pmarginal(0,rrm$marginals.hyperpar$`Theta2 for idy`)
      
      bests2[iter, siter] = b_est
      bposs2[iter, siter] = bpos
      
    }
    
    write.table(bests2,paste0("tempbest-n",n,".txt"))
    write.table(bposs2,paste0("tempbpos-n",n,".txt"))
    
    
  }
  
  bestfull = as.data.frame(matrix(0,nrow=17,ncol=nsims))
  bposfull = as.data.frame(matrix(0,nrow=17,ncol=nsims))
  bestfull[1:10,] = bests[1:10,]
  bposfull[1:10,] = bposs[1:10,]
  for(i in 1:8){
    bestfull[18-i,] = bests2[i,]
    bposfull[18-i,] = bposs2[i,]
  }
  
  
  sum(bposs[1:9,]>0.95)
  sum(bposs[9:17,]<0.05)
  
  sum(bposfull[1:9,]>0.95)
  sum(bposfull[9:17,]<0.05)
  
  
  btrue = seq(from=-0.8,to=0.8,by=0.1)
  #temp1best = read.table("tempbest-n500.txt") #from 0.1 to 0.8
  #temp1pos
  bmean500 = rowMeans(temp1best)#bests)
  bpos500 = rowMeans(temp1pos)#bposs)
  bmean1000 = rowMeans(temp2best)#bestfull)
  bpos1000 = rowMeans(temp2pos)#bposfull)
  
  bpos500_num = bpos1000_num = numeric(17)
  
  # bL500 = numeric(17)
  # bposL500 = numeric(17)
  # bL1000 = numeric(17)
  # bposL1000 = numeric(17)
  # bU500 = numeric(17)
  # bposU500 = numeric(17)
  # bU1000 = numeric(17)
  # bposU1000 = numeric(17)
   for(i in 1:17){
    #denst11 = density(bests[i,]); denst1 = data.frame(x=denst11$x,y=denst11$y)
    #densp11 = density(bposs[i,]); densp1 = data.frame(x=densp11$x,y=densp11$y)
    #denst22 = density(bestfull[i,]); denst2 = data.frame(x=denst22$x,y=denst22$y)
    #densp22 = density(bposfull[i,]); densp2 = data.frame(x=densp11$x,y=densp22$y)
    #zmt1 = inla.zmarginal(denst1,silent=TRUE)
    #zmt2 = inla.zmarginal(denst2,silent=TRUE)
    #zmp1 = inla.zmarginal(densp1,silent=TRUE)
    #zmp2 = inla.zmarginal(densp2,silent=TRUE)
    #bL500 = zmt1$quant0.025
    #bU500 = zmt1$quant0.975
    #bL1000 = zmt2$quant0.025
    #bU1000 = zmt2$quant0.975
    #bposL500 = zmp1$quant0.025
    #bposU500 = zmp1$quant0.975
    #bposL1000 = zmp2$quant0.025
    #bposU1000 = zmp2$quant0.975
    bpos500_num[i] = sum(temp1pos[i,]>0.95)
    bpos1000_num[i] = sum(temp2pos[i,]>0.95)
    cat(btrue[i], " & ",round(bmean500[i],digits=DIGITS)," & ", round(bmean1000[i],digits=DIGITS), " & ",
        round(bpos500[i],digits=DIGITS), " & ",round(bpos1000[i],digits=DIGITS),
        " & ",bpos500_num[i], " & ",bpos1000_num[i], "\\\\ \n", sep="")
  }
  
  
  rstr1 = paste0("True b")
  rstr2 = paste0("Mean b (n=500)")
  rstr3 = paste0("Mean b (n=1000)")
  rstr4 = paste0("P(b>0) (n=500)")
  rstr5 = paste0("P(b>0) (n=1000)")
  DIGITS=3
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
  
  
  library(robustbase)
  library(litteR)
  ggd1 = data.frame(b = factor(rep(seq(from=-0.8,to=0.8,by=0.1), each=nsims)), best = c(t(bests)), bpositive = c(t(bposs)))
  bmeanavg = rowMeans(bests)
  #ggd$b = as.factor(ggd$b)
  #write.table(ggd,"n500_bfits.txt")
  #temp = read.table("n500_bfits.txt")
  gg1 = ggplot(data=ggd1,aes(x=b,y=best)) + theme_bw() + xlab("True b") + ylab("Posterior marginal mean b") +
    #geom_boxplot(aes(x=b,y=best)) +
    stat_adj_boxplot() +
    stat_adj_boxplot_outlier()+
    geom_line(data=data.frame(bx=factor(seq(from=-0.8,to=0.8,by=0.1)),bgrid=seq(from=-0.8,to=0.8,by=0.1)),aes(x=factor(bx),y=bgrid, group=1), 
              linewidth=1.2, col="blue") +
    #   geom_line(data=data.frame(bx=factor(bgrid),bavg=bmeanavg),aes(x=factor(bx),y=bavg, group=1), 
    #            linewidth=1.2, col="red") +
    #geom_line(data=data.frame(bx=factor(bgrid),bavg=bmeanavg),aes(x=factor(bx),y=bavg, group=1), linewidth=1.2, col="blue") +
    #geom_boxplot(aes(x=b,y=best)) +
    ggtitle("(a) Posterior marginal mean", subtitle=paste0("n = ",500))
  
  gg2 = ggplot(data=ggd1,aes(x=b,y=bpositive)) + theme_bw() + xlab("True b") + ylab("P(b>0)") +
    #geom_boxplot() #+
    stat_adj_boxplot() +
    stat_adj_boxplot_outlier()+
    ggtitle("(c) Probability of positive trend", subtitle=paste0("n = ",500))
  #ggboth = ggarrange(gg1,gg2)
  #ggboth
  #ggsave(paste0("accuracytest-n",n,"-v3-28800x9600.eps"),plot=ggboth, device=cairo_ps, width=28800,
   #      height=9600, units="px", dpi=1800, limitsize=FALSE)
  
  
  #first_half = read.table()
  
  ggd2 = data.frame(b=factor(rep(seq(from=-0.8,to=0.8,by=0.1), each=nsims)), best=c(t(bestfull)),
                    bpositive=c(t(bposfull)))
  #ggd = data.frame(b = factor(rep(bgrid, each=nsims)), best = c(t(bests2)), bpositive = c(t(bposs2)))
  bmeanavg2 = rowMeans(bestfull)
  #ggd$b = as.factor(ggd$b)
  #write.table(ggd,"n500_bfits.txt")
  #temp = read.table("n500_bfits.txt")
  gg3 = ggplot(data=ggd2,aes(x=b,y=best)) + theme_bw() + xlab("True b") + ylab("Posterior marginal mean b") +
    #geom_boxplot() +
    stat_adj_boxplot() +
    stat_adj_boxplot_outlier()+
    geom_line(data=data.frame(bx=factor(bgrid),bgrid=bgrid),aes(x=factor(bx),y=bgrid, group=1), 
              linewidth=1.2, col="blue") +
    #   geom_line(data=data.frame(bx=factor(bgrid),bavg=bmeanavg),aes(x=factor(bx),y=bavg, group=1), 
    #            linewidth=1.2, col="red") +
    #geom_line(data=data.frame(bx=factor(bgrid),bavg=bmeanavg),aes(x=factor(bx),y=bavg, group=1), linewidth=1.2, col="blue") +
    #geom_boxplot(aes(x=b,y=best)) +
    ggtitle("(b) Posterior marginal mean", subtitle=paste0("n = ",1000))
  
  gg4 = ggplot(data=ggd2,aes(x=b,y=bpositive)) + theme_bw() + xlab("True b") + ylab("P(b>0)") +
    #geom_boxplot(aes(x=b,y=bpositive)) +
    stat_adj_boxplot() +
    stat_adj_boxplot_outlier()+
    ggtitle("(d) Probability of positive trend", subtitle=paste0("n = ",1000))
  ggboth = ggarrange(gg1,gg3,gg2,gg4, nrow=2,ncol=2)
  ggboth
  ggsave(paste0("accuracytest-full","-28800x19200-rev.eps"),plot=ggboth, device=cairo_ps, width=28800,
         height=19200, units="px", dpi=1800, limitsize=FALSE)
  
  
}