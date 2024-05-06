#source("rgeneric.ews.ar1.new.R")
if(FALSE){
  
  lp333 <- function(theta) {return(
    dnorm(theta[1], sd=10, log=TRUE) +
      dnorm(theta[2], sd=10, log=TRUE) +
      dnorm(theta[3], sd=10, log=TRUE))
  }
  
  library(bremla)
  
  data("events_rasmussen")
  #data("event_intervals")
  data("NGRIP_5cm")
  
  nevents = length(event_intervals[,1])  
  
  plot(NGRIP_5cm$age,NGRIP_5cm$d18O,type="l")
  abline(v=events_rasmussen$age)
  
  GS_onsets = which(grepl("Start of GS", events_rasmussen$event, fixed = TRUE))
  GI_onsets = which(grepl("Start of GI", events_rasmussen$event, fixed = TRUE))
  
  
  #Clear_GS_onsets0 = c(9,21,23,25,29,31,33,37,41,43,45,47,51,57,61,63,65,69,71,75,77)
  #Clear_GI_onsets0 = c(8,16,22,24,26,30,32,36,40,42,44,46,50,54,60,62,64,68,70,74,76)
  
  Clear_GS_onsets = c(9,21,23,25,29,31,33,37,41,43,45,47,51,55,61,65,71,77)
  Clear_GI_onsets = c(8,16,22,24,26,30,32,36,40,42,44,46,50,54,60,64,70,76)
  
  
  #time = seq(11728,59920,5)
  maxage = 59920.5
  minage= 11728
  whichind = NGRIP_5cm$age <= maxage & NGRIP_5cm$age >= minage
  time = rev(NGRIP_5cm$age[whichind])
  proxy = rev(NGRIP_5cm$d18O[whichind])
  
  end0 = events_rasmussen$age[Clear_GS_onsets[length(Clear_GS_onsets)]]
  start0 = events_rasmussen$age[Clear_GI_onsets[1]]
  int0 = which(time>start0 & time<end0)
  
  
  
  plot(time[int0],proxy[int0],type="l",xlim=rev(range(time[int0])),col="Red")
  
  bluecolor = rgb(0,115,150, maxColorValue=255)
  redcolor = rgb(203,51,59, maxColorValue=255)
  blackcolor = rgb(0,51,73, maxColorValue=255)
  yellowcolor = rgb(242,169,0, maxColorValue=255)

    ggy = ggplot() + theme_bw() + xlab("Time (yr b2k)") + 
    ylab(expression(paste(delta^18,"O (permil)"))) +
    theme(text=element_text(size=16), plot.title = element_text(size=22)) + 
    ggtitle("Stadial and interstadial periods",subtitle=expression(paste("NGRIP ",delta^18,"O record")))+
      #xlim(rev(range(time)))
      xlim(c(58560, min(time)))

  for (i in 1:17) {
    end = events_rasmussen$age[Clear_GS_onsets[i]]
    start = events_rasmussen$age[Clear_GI_onsets[i]]
    
    int = which(time>start & time<end)
    lines(x=time[int],y=proxy[int],col="Blue")
    ggd = data.frame(xx=time[int],yy=proxy[int])
    ggy = ggy + geom_line(data=ggd, aes(x=xx,y=yy),col=bluecolor)
    if(i <17){
      after = events_rasmussen$age[Clear_GS_onsets[i+1]]
      int2 = which(time>end & time<after)
      ggd = data.frame(xx=time[int2],yy=proxy[int2])
      ggy = ggy + geom_line(data=ggd, aes(x=xx,y=yy),col=redcolor)
    }
    ggdv = c(time[int[1]], time[int[length(int)]])
    print(ggdv)
    ggy = ggy + geom_vline(xintercept=ggdv[1], col=bluecolor)
    ggy = ggy + geom_vline(xintercept=ggdv[2], col=redcolor)
    #print(plot(x=time[int],y=proxy[int],type="l",main=i))
  }
    print(ggy)
    ggsave(paste0("events2.eps"),plot=ggy, device=cairo_ps, width=14,
           height=4.7, limitsize=FALSE)
    
  events = events_rasmussen$age
  abline(v=events,col="gray",lwd=0.8)
  abline(v=c(events[Clear_GS_onsets]),col="Blue",lwd=1.5)
  abline(v=c(events[Clear_GI_onsets]),col="Red",lwd=0.9)
  
  
  
  
  ###### start things #####
  
  data_proxy = list()
  data_time = list()
  for (i in 1:length(Clear_GI_onsets)) {
    end = events_rasmussen$age[Clear_GS_onsets[i]]
    start = events_rasmussen$age[Clear_GI_onsets[i]]
    int = which(time>start & time<end)
    data_proxy[[i]] = proxy[int]
    data_time[[i]] = time[int]
  }
  
  
  
  
  
  
  #### fit things
  library(INLA)
  
  
  
  library(ggplot2)
  library(ggpubr)
  ggps0 = c()
  ggps1 = c()
  ggps2 = c()
  ggpsrw2 = c()
  ggms1 = c()
  ggms2 = c()
  ggmsrw2 = c()
  rres0 = c()
  rres1 = c()
  rres2 = c()
  rresrw2 = c()
  
  bpos0 = numeric(length(data_time))
  bpos1 = numeric(length(data_time))
  bpos2 = numeric(length(data_time))
  bposrw2 = numeric(length(data_time))
  bmean0 = numeric(length(data_time))
  bmean1 = numeric(length(data_time))
  bmean2 = numeric(length(data_time))
  bmeanrw2 = numeric(length(data_time))
  amean0 = numeric(length(data_time))
  amean1 = numeric(length(data_time))
  amean2 = numeric(length(data_time))
  ameanrw2 = numeric(length(data_time))
  
  
  # m = get("inla.models", inla.get.inlaEnv())
  # m$latent$rw2$min.diff = NULL
  # assign("inla.models", m, inla.get.inlaEnv())
  
  #bypass = rep(FALSE, length(data_time))
  #bypass[2] = TRUE
  for( i in 1:length(data_time)){
    cat("running event #",i,", n = ",length(data_time[[i]]),"\n",sep="")
    fliptime=TRUE
    n = length(data_time[[i]])
    
    pproxy = data_proxy[[i]]
    ttime = data_time[[i]]
    time_norm = (ttime-min(ttime))/max(ttime-min(ttime))
    time_normalized=time_norm
    if(fliptime) time_normalized = 1-time_norm #make it increasing to avoid errors
    
    trend1 = time_normalized
    trend2 = time_normalized^2
    
    
    
    #plot(pproxy)
    cat("\nNo trend:\n",sep="")
    res0 = inla.ews(data=pproxy,print.progress=TRUE, timesteps = ttime, log.prior=lp333)
    #lines(res0$results$summary$trend$mean,col="blue",lwd=5)
    cat("\nLinear trend:\n",sep="")
    res1 = inla.ews(data=data.frame(y=pproxy,trend1=trend1), timesteps = ttime,log.prior=lp333,
                    formula=y~1+trend1,print.progress=TRUE)
    #lines(res1$results$summary$trend$mean,col="red",lwd=5)
    cat("\nSquare trend:\n",sep="")
    res2 = inla.ews(data=data.frame(y=pproxy,trend1=trend1,trend2=trend2), timesteps=ttime,log.prior=lp333,
                    formula=y~1+trend1+trend2,print.progress=TRUE)
    #lines(res2$results$summary$trend$mean,col="green",lwd=5)
    cat("\nRW2 trend:\n",sep="")
    # if(bypass[i]){
    #   #m = get("inla.models", inla.get.inlaEnv())
    #   #m$latent$rw2$min.diff = NULL
    #   #assign("inla.models", m, inla.get.inlaEnv())
    #   resrw2 = inla.ews(data=data.frame(y=pproxy,idx=time_normalized), timesteps = ttime,
    #                     formula=y~1+f(inla.group(idx),model="crw2"),print.progress=TRUE)
    # }else{
    ##doing some risky stuff to avoid error
    m = get("inla.models", inla.get.inlaEnv()) 
    m.old = m$latent$crw2$min.diff
    m$latent$crw2$min.diff = NULL
    assign("inla.models", m, inla.get.inlaEnv())
    # resrw2 = inla.ews(data=data.frame(y=pproxy,idx=time_normalized*(n-1)+1), timesteps = ttime,
    resrw2 = inla.ews(data=data.frame(y=pproxy,idx=ttime), timesteps = ttime,log.prior=lp333,
                      formula=y~1+f(idx,model="crw2"),print.progress=TRUE)
    ##setting things back the way they were
    m$latent$crw2$min.diff = m.old  
    assign("inla.models", m, inla.get.inlaEnv())
    # }
    
    #lines(resrw2$results$summary$trend$mean,col="orange",lwd=5)
    
    ggp0 = ggplot(data=as.data.frame(res0$results$summary$phi0), aes(ttime)) + theme_bw() +
      ylim(range(0,1,res0$results$summary$phi0$q0.025[2:n],res0$results$summary$phi0$q0.975[2:n])) + 
      xlim(c(ttime[1],tail(ttime,1))) +
      geom_line(data=data.frame(ttime=ttime[2:n],phimean=res0$results$summary$phi$mean[2:n]),
                aes(x=ttime,y=phimean),col="gray") +
      geom_ribbon(aes(ymin=q0.025,ymax=q0.975),col="red",fill="red",alpha=0.3,linewidth=1.2) +
      geom_line(aes(y=mean),col="blue", linewidth=1.2) +
      xlab("Time (yr b2k)") + ylab(expression(paste(phi,"(t)"))) +
      ggtitle(paste0("Event #",i,": No trend"), 
              subtitle = paste0("P(b>0) = ",round(res0$results$summary$b$prob_positive,digits=3)))
    
    ggp1 = ggplot(data=as.data.frame(res1$results$summary$phi0), aes(ttime)) + theme_bw() +
      ylim(range(0,1,res1$results$summary$phi0$q0.025[2:n],res1$results$summary$phi0$q0.975[2:n])) + 
      xlim(c(ttime[1],tail(ttime,1))) +
      geom_line(data=data.frame(ttime=ttime[2:n],phimean=res1$results$summary$phi$mean[2:n]),
                aes(x=ttime,y=phimean),col="gray") +
      geom_ribbon(aes(ymin=q0.025,ymax=q0.975),col="red",fill="red",alpha=0.3,linewidth=1.2) +
      geom_line(aes(y=mean),col="blue", linewidth=1.2) +
      xlab("Time (yr b2k)") + ylab(expression(paste(phi,"(t)"))) +
      ggtitle(paste0("Event #",i,": Linear trend"), 
              subtitle = paste0("P(b>0) = ",round(res1$results$summary$b$prob_positive,digits=3)))
    
    ggp2 = ggplot(data=as.data.frame(res2$results$summary$phi0), aes(ttime)) + theme_bw() +
      ylim(range(0,1,res2$results$summary$phi0$q0.025[2:n],res2$results$summary$phi0$q0.975[2:n])) + 
      xlim(c(ttime[1],tail(ttime,1))) +
      geom_line(data=data.frame(ttime=ttime[2:n],phimean=res2$results$summary$phi$mean[2:n]),
                aes(x=ttime,y=phimean),col="gray") +
      geom_ribbon(aes(ymin=q0.025,ymax=q0.975),col="red",fill="red",alpha=0.3,linewidth=1.2) +
      geom_line(aes(y=mean),col="blue", linewidth=1.2) +
      xlab("Time (yr b2k)") + ylab(expression(paste(phi,"(t)"))) +
      ggtitle(paste0("Event #",i,": Square trend"), 
              subtitle = paste0("P(b>0) = ",round(res2$results$summary$b$prob_positive,digits=3)))
    
    ggprw2 = ggplot(data=as.data.frame(resrw2$results$summary$phi0), aes(ttime)) + theme_bw() +
      ylim(range(0,1,resrw2$results$summary$phi0$q0.025[2:n],resrw2$results$summary$phi0$q0.975[2:n])) + 
      xlim(c(ttime[1],tail(ttime,1))) +
      geom_line(data=data.frame(ttime=ttime[2:n],phimean=resrw2$results$summary$phi$mean[2:n]),
                aes(x=ttime,y=phimean),col="gray") +
      geom_ribbon(aes(ymin=q0.025,ymax=q0.975),col="red",fill="red",alpha=0.3,linewidth=1.2) +
      geom_line(aes(y=mean),col="blue", linewidth=1.2) +
      xlab("Time (yr b2k)") + ylab(expression(paste(phi,"(t)"))) +
      ggtitle(paste0("Event #",i,": RW2 trend"), 
              subtitle = paste0("P(b>0) = ",round(resrw2$results$summary$b$prob_positive,digits=3)))
    
    print(ggarrange(ggp0,ggp1,ggp2,ggprw2))
    
    ggm0 = ggplot(data=as.data.frame(res0$results$summary$trend),aes(ttime)) + theme_bw() +
      ylim(range(pproxy,res0$results$summary$trend$quant0.025,res0$results$summary$trend$quant0.975)) + 
      xlim(c(ttime[1],tail(ttime,1))) +
      xlab("Time (yr b2k)") + ylab(expression(paste(delta^18,"O (permil)"))) +
      ggtitle(paste0("Event #",i,": No trend")) +
      geom_line(data=data.frame(y=pproxy,ttime=ttime),aes(x=ttime,y=y)) +
      geom_ribbon(aes(ymin=quant0.025,ymax=quant0.975), fill="red",col="red",alpha=0.3,linewidth=1.2)+
      geom_line(aes(y=mean),col="blue",linewidth=1.2)
    
    ggm1 = ggplot(data=as.data.frame(res1$results$summary$trend),aes(ttime)) + theme_bw() +
      ylim(range(pproxy,res1$results$summary$trend$quant0.025,res1$results$summary$trend$quant0.975)) + 
      xlim(c(ttime[1],tail(ttime,1))) +
      xlab("Time (yr b2k)") + ylab(expression(paste(delta^18,"O (permil)"))) +
      ggtitle(paste0("Event #",i,": Linear trend")) +
      geom_line(data=data.frame(y=pproxy,ttime=ttime),aes(x=ttime,y=y)) +
      geom_ribbon(aes(ymin=quant0.025,ymax=quant0.975), fill="red",col="red",alpha=0.3,linewidth=1.2)+
      geom_line(aes(y=mean),col="blue",linewidth=1.2)
    
    ggm2 = ggplot(data=as.data.frame(res2$results$summary$trend),aes(ttime)) + theme_bw() +
      ylim(range(pproxy,res2$results$summary$trend$quant0.025,res2$results$summary$trend$quant0.975)) + 
      xlim(c(ttime[1],tail(ttime,1))) +
      xlab("Time (yr b2k)") + ylab(expression(paste(delta^18,"O (permil)"))) +
      ggtitle(paste0("Event #",i,": Square trend")) +
      geom_line(data=data.frame(y=pproxy,ttime=ttime),aes(x=ttime,y=y)) +
      geom_ribbon(aes(ymin=quant0.025,ymax=quant0.975), fill="red",col="red",alpha=0.3,linewidth=1.2)+
      geom_line(aes(y=mean),col="blue",linewidth=1.2)
    
    ggmrw2 = ggplot(data=as.data.frame(resrw2$results$summary$trend),aes(ttime)) + theme_bw() +
      ylim(range(pproxy,resrw2$results$summary$trend$quant0.025,resrw2$results$summary$trend$quant0.975)) + 
      xlim(c(ttime[1],tail(ttime,1))) +
      xlab("Time (yr b2k)") + ylab(expression(paste(delta^18,"O (permil)"))) +
      ggtitle(paste0("Event #",i,": RW2 trend")) +
      geom_line(data=data.frame(y=pproxy,ttime=ttime),aes(x=ttime,y=y)) +
      geom_ribbon(aes(ymin=quant0.025,ymax=quant0.975), fill="red",col="red",alpha=0.3,linewidth=1.2)+
      geom_line(aes(y=mean),col="blue",linewidth=1.2)
    
    print(ggarrange(ggm0,ggm1,ggm2,ggmrw2))
    
    bpos0[i] = res0$results$summary$b$prob_positive
    bpos1[i] = res1$results$summary$b$prob_positive
    bpos2[i] = res2$results$summary$b$prob_positive
    bposrw2[i] = resrw2$results$summary$b$prob_positive
    
    bmean0[i] = res0$results$summary$b$mean
    bmean1 = res1$results$summary$b$mean
    bmean2 = res2$results$summary$b$mean
    bmeanrw2 = resrw2$results$summary$b$mean
    
    amean0 = res0$results$summary$a$mean
    amean1 = res1$results$summary$a$mean
    amean2 = res2$results$summary$a$mean
    ameanrw2 = resrw2$results$summary$a$mean
    
    cat("P(b>0): \n","\tNo trend:\t",
        round(bpos0[i],digits=4),"\n\tLinear trend:\t",
        round(bpos1[i],digits=4),"\n\tSquare trend:\t",
        round(bpos2[i],digits=4),"\n\tRW2 trend:\t",
        round(bposrw2[i],digits=4),"\n",sep="")
    
    
    ggps0 = c(ggps0, list(ggps0=ggp0))
    ggps1 = c(ggps1, list(ggps1=ggp1))
    ggps2 = c(ggps2, list(ggps2=ggp2))
    ggpsrw2 = c(ggpsrw2, list(ggpsrw2=ggprw2))
    ggms1 = c(ggms1, list(ggpm1=ggm1))
    ggms2 = c(ggms2, list(ggpm2=ggm2))
    ggmsrw2 = c(ggmsrw2, list(ggpmrw2=ggmrw2))
    rres0 = c(rres0, list(rres0=res0))
    rres1 = c(rres1, list(rres1=res1))
    rres2 = c(rres2, list(rres1=res2))
    rresrw2 = c(rresrw2, list(rres1=resrw2))
  }#end for
  
  
  rypdals = c("$p=0.2$", "$p=0.008$", "$-$","$-$","$p=0.13$", "$-$","$-$","$-$","$p=0.16$",
              "$-$","$-$","$-$","$0.39$","$-$","$-$","$-$","$-$")
  boers = c("$-$","$p<0.05$","$-$","$p<0.05$","$-$","$p<0.05$","$-$","$-$","$-$","$-$",
            "$p<0.05$","$-$","$p<0.05$","$p<0.05$","$p<0.05$","$-$","$-$")
  for(i in 1:17){
    cat(i," & ", round(bpos0[i],digits=4)," & ", round(bpos1[i],digits=4), " & ", round(bpos2[i],digits=4), " & ",
        round(bposrw2[i],digits=4), " & ", rypdals[i]," & ", boers[i]," \\\\ \n",sep="")
  }
  endat = 17
  cat("Ensemble & ", round(mean(bpos0[1:endat]),digits=4), " & ", round(mean(bpos1[1:endat]),digits=4), " & ",
      round(mean(bpos2[1:endat]),digits=4)," & ",round(mean(bposrw2[1:endat]),digits=4)," & ",
      "$-$"," & ","$-$","\\\\ \n")
  #
  
  fulltime = time
  fullproxy = proxy
  fulln = length(fulltime)
  
  pmean0 = pmean1 = pmean2 = pmeanrw2 = rep(NA,fulln)
  p0mean0 = p0mean1 = p0mean2 = p0meanrw2 = rep(NA,fulln)
  pL0 = pL1 = pL2 = pLrw2 = rep(NA,fulln)
  pU0 = pU1 = pU2 = pUrw2 = numeric(fulln)
  mmean0 = mmean1 = mmean2 = mmeanrw2 = rep(NA,fulln)
  mL0 = mL1 = mL2 = mLrw2 = rep(NA,fulln)
  mU0 = mU1 = mU2 = mUrw2 = rep(NA,fulln)
  
  for(i in 1:length(bpos0)){
    do.rev = TRUE
    ind0 = which(time==data_time[[i]][1])
    indn = which(time==data_time[[i]][length(data_time[[i]])])
    ind = ind0:indn
    if(do.rev){
      pmean0[ind] = rres0[[i]]$results$summary$phi$mean
      p0mean0[ind] = rres0[[i]]$results$summary$phi0$mean
      pL0[ind] = rres0[[i]]$results$summary$phi$q0.025
      pU0[ind] = rres0[[i]]$results$summary$phi$q0.975
      pmean1[ind] = rres1[[i]]$results$summary$phi$mean
      p0mean1[ind] = rres1[[i]]$results$summary$phi0$mean
      pL1[ind] = rres1[[i]]$results$summary$phi$q0.025
      pU1[ind] = rres1[[i]]$results$summary$phi$q0.975
      pmean2[ind] = rres2[[i]]$results$summary$phi$mean
      p0mean2[ind] = rres2[[i]]$results$summary$phi0$mean
      pL2[ind] = rres2[[i]]$results$summary$phi$q0.025
      pU2[ind] = rres2[[i]]$results$summary$phi$q0.975
      pmeanrw2[ind] = rresrw2[[i]]$results$summary$phi$mean
      p0meanrw2[ind] = rresrw2[[i]]$results$summary$phi0$mean
      pLrw2[ind] = rresrw2[[i]]$results$summary$phi$q0.025
      pUrw2[ind] = rresrw2[[i]]$results$summary$phi$q0.975
    }
  }
  
  
  
  
  
  ### FROM LUC
  
  
  
  df_plot_RW2 = list()
  df_plot_Lin = list()
  df_plot_Quad = list()
  df_plot_Nt = list()
  
  for (i in 1:length(data_time)) {
    
    
    df_plot_RW2[[i]]= data.frame(Time=data_time[[i]],
                                 Phi=rresrw2[[i]]$results$summary$phi0$mean,
                                 Phi1=rresrw2[[i]]$results$summary$phi$mean,
                                 phi_Low = rresrw2[[i]]$results$summary$phi0$q0.025,
                                 phi_Upp= rresrw2[[i]]$results$summary$phi0$q0.975,
                                 data=rresrw2[[i]]$.args$data$y,
                                 trend=rresrw2[[i]]$results$summary$trend$mean,
                                 trend_low=rresrw2[[i]]$results$summary$tren$quant0.025,
                                 trend_upp=rresrw2[[i]]$results$summary$tren$quant0.975
    )
    
    df_plot_Lin[[i]]= data.frame(Time=data_time[[i]],
                                 Phi=rres1[[i]]$results$summary$phi0$mean,
                                 Phi1=rresrw2[[i]]$results$summary$phi$mean,
                                 phi_Low = rres1[[i]]$results$summary$phi0$q0.025,
                                 phi_Upp= rres1[[i]]$results$summary$phi0$q0.975,
                                 data=rres1[[i]]$.args$data$y,
                                 trend=rres1[[i]]$results$summary$trend$mean,
                                 trend_low=rres1[[i]]$results$summary$tren$quant0.025,
                                 trend_upp=rres1[[i]]$results$summary$tren$quant0.975
    )
    
    df_plot_Quad[[i]]= data.frame(Time=data_time[[i]],
                                  Phi=rres2[[i]]$results$summary$phi0$mean,
                                  Phi1=rresrw2[[i]]$results$summary$phi$mean,
                                  phi_Low = rres2[[i]]$results$summary$phi0$q0.025,
                                  phi_Upp= rres2[[i]]$results$summary$phi0$q0.975,
                                  data=rres2[[i]]$.args$data$y,
                                  trend=rres2[[i]]$results$summary$trend$mean,
                                  trend_low=rres2[[i]]$results$summary$tren$quant0.025,
                                  trend_upp=rres2[[i]]$results$summary$tren$quant0.975
    )
    
    df_plot_Nt[[i]]= data.frame(Time=data_time[[i]],
                                Phi=rres0[[i]]$results$summary$phi0$mean,
                                Phi1=rresrw2[[i]]$results$summary$phi$mean,
                                phi_Low = rres0[[i]]$results$summary$phi0$q0.025,
                                phi_Upp= rres0[[i]]$results$summary$phi0$q0.975,
                                data=rres0[[i]]$.args$data,
                                trend=rres0[[i]]$results$summary$trend$mean,
                                trend_low=rres0[[i]]$results$summary$tren$quant0.025,
                                trend_upp=rres0[[i]]$results$summary$tren$quant0.975
    )
    
    
  }
  
  library(patchwork)
  library(ggbreak)
  
  trendselect = "sq"
  
  dfselect = list()
  for(i in 1:17){
    if(trendselect %in% c("no","nt","no trend")){
      dfselect = c(dfselect, list(df_plot_Nt[[i]]))
      selectstr = "No trend"
      trendstr = "(a) No trend"
      savestr = "-notrend"
      bposselect = bpos0
    }else if(trendselect %in% c("lin", "linear")){
      dfselect = c(dfselect, list(df_plot_Lin[[i]]))
      selectstr = "Linear trend"
      trendstr = "(b) Linear trend"
      savestr = "-lintrend"
      bposselect = bpos1
    }else if(trendselect %in% c("sq", "2")){
      dfselect = c(dfselect, list(df_plot_Quad[[i]]))
      selectstr = "Square trend"
      trendstr = "(c) Square trend"
      savestr = "-sqtrend"
      bposselect = bpos2
    }else if(trendselect %in% c("rw2")){
      dfselect = c(dfselect, list(df_plot_RW2[[i]]))
      selectstr = "RW2 trend"
      trendstr = "(d) RW2 trend"
      savestr = "-rw2trend"
      bposselect = bposrw2
    }
  }
  
  ## Plot 1
  
  gg_phi_Lin1 <- ggplot(data=dfselect[[1]],aes(Time)) + 
    geom_ribbon(data=dfselect[[1]],aes(ymin=phi_Low,ymax=phi_Upp),col="red",fill="red",alpha=0.3,linewidth=0.5) +
    #geom_line(data=dfselect[[1]],aes(y=Phi1),col="gray")+
    geom_line(data=dfselect[[1]],aes(y=Phi),col="blue")+
    theme_bw()+ labs(title="Lag-One Autocorrelation Evolution", subtitle= paste0(selectstr))+
    ylab("")+
    xlab("")+
    annotate(geom="text", x=mean(range(data=dfselect[[1]]$Time)), y=0.9, label=paste0("Event #",as.character(1)))+
    annotate(geom="text", x=mean(range(data=dfselect[[1]]$Time)), y=0.8, label=paste0("P(b>0) = ",as.character(round(bposselect[1],2))))+
    theme(text = element_text(size = 20))
  for (i in 2:3) {
    gg_phi_Lin1<- gg_phi_Lin1 + 
      geom_ribbon(data=dfselect[[i]],aes(ymin=phi_Low,ymax=phi_Upp),col="red",fill="red",alpha=0.3,linewidth=0.5)+
      #geom_line(data=dfselect[[i]],aes(y=Phi1),col="gray")+
      geom_line(data=dfselect[[i]],aes(x=Time,y=Phi),col="blue") +
      annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.9, label=paste0("Event #",as.character(i)))+
      annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.8, label=paste0("P(b>0) = ",as.character(round(bposselect[i],2))))
  }
  gg_phi_Lin1 <- gg_phi_Lin1 + 
    geom_vline(xintercept=c(events[Clear_GS_onsets]),color ="Blue") +
    geom_vline(xintercept=c(events[Clear_GI_onsets]),color ="red")+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
    coord_cartesian(ylim=c(0,1),xlim =c(max(dfselect[[3]]$Time)+100, min(dfselect[[1]]$Time)-50),expand = FALSE)#+
  # theme(panel.grid.major = element_blank(), 
  #        panel.grid.minor = element_blank())
  
  ## Plot 2
  
  gg_phi_Lin2 <- ggplot(data=dfselect[[4]],aes(Time)) + 
    geom_ribbon(data=dfselect[[4]],aes(ymin=phi_Low,ymax=phi_Upp),col="red",fill="red",alpha=0.3,linewidth=0.5) +
    #geom_line(data=dfselect[[4]],aes(y=Phi1),col="gray")+
    geom_line(data=dfselect[[4]],aes(y=Phi),col="blue")+
    theme_bw()+ #labs(title="Linratic Detrending ")+
    ylab(expression(paste(phi,"(t)")))+
    xlab("")+#Time (yr b2k)")+
    annotate(geom="text", x=mean(range(data=dfselect[[4]]$Time)), y=0.9, label=paste0("Event #",as.character(4)))+
    annotate(geom="text", x=mean(range(data=dfselect[[4]]$Time)), y=0.8, label=paste0("P(b>0) = ",as.character(round(bposselect[4],2))))+
    theme(text = element_text(size = 20))
  for (i in 5:11) {
    gg_phi_Lin2<- gg_phi_Lin2 + 
      geom_ribbon(data=dfselect[[i]],aes(ymin=phi_Low,ymax=phi_Upp),col="red",fill="red",alpha=0.3,linewidth=0.5)+
      #geom_line(data=dfselect[[i]],aes(x=Time,y=Phi1),col="gray") +
      geom_line(data=dfselect[[i]],aes(x=Time,y=Phi),col="blue") +
      annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.9, label=paste0("Event #",as.character(i)))+
      annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.8, label=paste0("P(b>0) = ",as.character(round(bposselect[i],2))))
  }
  
  gg_phi_Lin2 <- gg_phi_Lin2 + 
    geom_vline(xintercept=c(events[Clear_GS_onsets]),color ="Blue") +
    geom_vline(xintercept=c(events[Clear_GI_onsets]),color ="red")+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
    coord_cartesian(ylim=c(0,1),xlim =c(max(dfselect[[11]]$Time)+100, min(dfselect[[4]]$Time)-50),expand = FALSE)#+
  # theme(panel.grid.major = element_blank(), 
  #        panel.grid.minor = element_blank())
  
  ## Plot 3
  
  gg_phi_Lin3 <- ggplot(data=dfselect[[12]],aes(Time)) + 
    geom_ribbon(data=dfselect[[12]],aes(ymin=phi_Low,ymax=phi_Upp),col="red",fill="red",alpha=0.3,linewidth=0.5) +
    #geom_line(data=dfselect[[12]],aes(y=Phi1),col="gray")+
    geom_line(data=dfselect[[12]],aes(y=Phi),col="blue")+
    theme_bw()+ #labs(title="Linratic Detrending ")+
    ylab("")+#expression(paste(delta^18,"O (permil)")))+
    xlab("Time (yr b2k)")+
    annotate(geom="text", x=mean(range(data=dfselect[[12]]$Time)), y=0.9, label=paste("Event #",as.character(12)))+
    annotate(geom="text", x=mean(range(data=dfselect[[12]]$Time)), y=0.8, label=paste("P(b>0) = ",as.character(round(bposselect[12],2))))+
    theme(text = element_text(size = 20))
  for (i in 13:17) {
    gg_phi_Lin3<- gg_phi_Lin3 + 
      geom_ribbon(data=dfselect[[i]],aes(ymin=phi_Low,ymax=phi_Upp),col="red",fill="red",alpha=0.3,linewidth=0.5)+
      #geom_line(data=dfselect[[i]],aes(x=Time,y=Phi1),col="gray") +
      geom_line(data=dfselect[[i]],aes(x=Time,y=Phi),col="blue") +
      annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.9, label=paste0("Event #",as.character(i)))+
      annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.8, label=paste0("P(b>0) = ",as.character(round(bposselect[i],2))))
  }
  
  gg_phi_Lin3 <- gg_phi_Lin3 + 
    geom_vline(xintercept=c(events[Clear_GS_onsets[12:17]]),color ="Blue") +
    geom_vline(xintercept=c(events[Clear_GI_onsets[12:17]]),color ="red")+
    scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
    coord_cartesian(ylim=c(0,1),xlim =c(max(dfselect[[17]]$Time)+100, min(dfselect[[12]]$Time)-50),expand = FALSE)#+
  #  theme(panel.grid.major = element_blank(), 
  #       panel.grid.minor = element_blank())
  
  
  
  
  
  
  gg_final_Phi = gg_phi_Lin1/gg_phi_Lin2/gg_phi_Lin3
  
  
  print(gg_final_Phi) ########  To Print
  
#  ggsave(paste0("AC1_RW2_evolution.eps"),plot=gg_final_Phi, device=cairo_ps, width=28800,
 #        height=19000, units="px", dpi=1800, limitsize=FALSE)
  ggsave(paste0("AC1_RW2_evolution.eps"),plot=gg_final_Phi, device=cairo_ps, width=14,
         height=10.5, limitsize=FALSE)
#  ggsave(paste0("AC1","RW2","_evolution.eps"),plot=plotlist[[3]], device=cairo_ps, width=8,
#         height=6, units="cm")
  #ggsave("breakpoint-plot-8x4-vague.eps",plot=ggboth, device=cairo_ps, width=8,
  #       height=6)
  
  ## trend
  
  trendset = c("no","lin","sq","rw2")
  plotlist = list()
  for(ttt in 1:4){
    trendselect = trendset[ttt]
    dfselect = list()
    for(i in 1:17){
      if(trendselect %in% c("no","nt","no trend")){
        dfselect = c(dfselect, list(df_plot_Nt[[i]]))
        selectstr = "No trend"
        trendstr = "(a) No trend"
        savestr = "-notrend"
      }else if(trendselect %in% c("lin", "linear")){
        dfselect = c(dfselect, list(df_plot_Lin[[i]]))
        selectstr = "Linear trend"
        trendstr = "(b) Linear trend"
        savestr = "-lintrend"
      }else if(trendselect %in% c("sq", "2")){
        dfselect = c(dfselect, list(df_plot_Quad[[i]]))
        selectstr = "Square trend"
        trendstr = "(c) Square trend"
        savestr = "-sqtrend"
      }else if(trendselect %in% c("rw2")){
        dfselect = c(dfselect, list(df_plot_RW2[[i]]))
        selectstr = "RW2 trend"
        trendstr = "(d) RW2 trend"
        savestr = "-rw2trend"
      }
    }
    
    
    
    ## Plot 1
    
    gg_trend_Lin1 <- ggplot(data=dfselect[[1]],aes(Time)) + 
      geom_line(data=data.frame(time=time,y=proxy),aes(x=time,y=proxy),col="gray")+
      geom_line(data=dfselect[[1]],aes(y=trend),col="blue")+
      geom_ribbon(data=dfselect[[1]],aes(ymin=trend_low,ymax=trend_upp),col="red",fill="red",alpha=0.3,linewidth=0.5) +
      theme_bw()+ labs(title=paste0(selectstr))+
      theme(text=element_text(size=15), plot.title = element_text(size=20))+
      ylab("")+
      xlab("")#+
    #annotate(geom="text", x=mean(range(data=dfselect[[1]]$Time)), y=0.9, label=paste0("Event #",as.character(1)))+
    #annotate(geom="text", x=mean(range(data=dfselect[[1]]$Time)), y=0.8, label=paste0("P(b>0) = ",as.character(round(bposrw2[1],2))))+
    #theme(text = element_text(size = 20))
    for (i in 2:3) {
      gg_trend_Lin1<- gg_trend_Lin1 + 
        geom_line(data=dfselect[[i]],aes(x=Time,y=trend),col="blue") +
        geom_ribbon(data=dfselect[[i]],aes(ymin=trend_low,ymax=trend_upp),col="red",fill="red",alpha=0.3,linewidth=0.5)#+
      #annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.9, label=paste0("Event #",as.character(i)))+
      #annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.8, label=paste0("P(b>0) = ",as.character(round(bposrw2[i],2))))
    }
    gg_trend_Lin1 <- gg_trend_Lin1 + 
      geom_vline(xintercept=c(events[Clear_GS_onsets]),color ="Blue") +
      geom_vline(xintercept=c(events[Clear_GI_onsets]),color ="red")+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
      coord_cartesian(ylim=range(proxy),xlim =c(max(dfselect[[3]]$Time)+100, min(dfselect[[1]]$Time)-50),expand = FALSE)#+
    #theme(panel.grid.major = element_blank(), 
    #panel.grid.minor = element_blank())
    
    ## Plot 2
    
    gg_trend_Lin2 <- ggplot(data=dfselect[[4]],aes(Time)) + 
      geom_line(data=data.frame(time=time,y=proxy),aes(x=time,y=proxy),col="gray")+
      geom_line(data=dfselect[[4]],aes(y=trend),col="blue")+
      theme(text=element_text(size=15))+
      geom_ribbon(data=dfselect[[4]],aes(ymin=trend_low,ymax=trend_upp),col="red",fill="red",alpha=0.3,linewidth=0.5) +
      
      theme_bw()+ #labs(title="Linratic Detrending ")+
      ylab(expression(paste(delta^18,"O (permil)")))+
      xlab("")#+#Time (yr b2k)") #+
    #annotate(geom="text", x=mean(range(data=dfselect[[4]]$Time)), y=0.9, label=paste0("Event #",as.character(4)))+
    #annotate(geom="text", x=mean(range(data=dfselect[[4]]$Time)), y=0.8, label=paste0("P(b>0) = ",as.character(round(bposrw2[4],2))))+
    #theme(text = element_text(size = 20))
    for (i in 5:11) {
      gg_trend_Lin2<- gg_trend_Lin2 + 
        geom_line(data=dfselect[[i]],aes(x=Time,y=trend),col="blue") +
        geom_ribbon(data=dfselect[[i]],aes(ymin=trend_low,ymax=trend_upp),col="red",fill="red",alpha=0.3,linewidth=0.5)#+
      #annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.9, label=paste0("Event #",as.character(i)))+
      #annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.8, label=paste0("P(b>0) = ",as.character(round(bposrw2[i],2))))
    }
    
    gg_trend_Lin2 <- gg_trend_Lin2 + 
      theme(text=element_text(size=15))+
      geom_vline(xintercept=c(events[Clear_GS_onsets]),color ="Blue") +
      geom_vline(xintercept=c(events[Clear_GI_onsets]),color ="red")+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
      coord_cartesian(ylim=range(proxy),xlim =c(max(dfselect[[11]]$Time)+100, min(dfselect[[4]]$Time)-50),expand = FALSE)#+
    #  theme(panel.grid.major = element_blank(), 
    #        panel.grid.minor = element_blank())
    
    ## Plot 3
    
    gg_trend_Lin3 <- ggplot(data=dfselect[[12]],aes(Time)) + 
      theme(text=element_text(size=15))+
      geom_line(data=data.frame(time=time,y=proxy),aes(x=time,y=proxy),col="gray")+
      geom_line(data=dfselect[[12]],aes(y=trend),col="blue")+
      geom_ribbon(data=dfselect[[12]],aes(ymin=trend_low,ymax=trend_upp),col="red",fill="red",alpha=0.3,linewidth=0.5) +
      theme_bw()+ #labs(title="Linratic Detrending ")+
      ylab("")+#expression(paste(delta^18,"O (permil)")))+
      xlab("Time (yr b2k)")#+
    #annotate(geom="text", x=mean(range(data=dfselect[[12]]$Time)), y=0.9, label=paste("Event #",as.character(12)))+
    #annotate(geom="text", x=mean(range(data=dfselect[[12]]$Time)), y=0.8, label=paste("P(b>0) = ",as.character(round(bposrw2[12],2))))+
    #theme(text = element_text(size = 20))
    for (i in 13:17) {
      gg_trend_Lin3<- gg_trend_Lin3 + 
        geom_line(data=dfselect[[i]],aes(x=Time,y=trend),col="blue") +
        geom_ribbon(data=dfselect[[i]],aes(ymin=trend_low,ymax=trend_upp),col="red",fill="red",alpha=0.3,linewidth=0.5)#3+
      #annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.9, label=paste0("Event #",as.character(i)))+
      #annotate(geom="text", x=mean(range(data=dfselect[[i]]$Time)), y=0.8, label=paste0("P(b>0) = ",as.character(round(bposrw2[i],2))))
    }
    
    gg_trend_Lin3 <- gg_trend_Lin3 + 
      theme(text=element_text(size=15))+
      geom_vline(xintercept=c(events[Clear_GS_onsets[12:17]]),color ="Blue") +
      geom_vline(xintercept=c(events[Clear_GI_onsets[12:17]]),color ="red")+
      scale_x_continuous(breaks = scales::pretty_breaks(n = 15)) +
      coord_cartesian(ylim=range(proxy),xlim =c(max(dfselect[[17]]$Time)+100, 
                                                min(dfselect[[12]]$Time)-50),expand = FALSE)#+
    #  theme(panel.grid.major = element_blank(), 
    #        panel.grid.minor = element_blank())
    
    gg_final_trend = gg_trend_Lin1/gg_trend_Lin2/gg_trend_Lin3
    
    
    print(gg_final_trend) ########  To Print
    
    # ggsave(paste0("trend",selectstr,"_evolution.eps"),plot=gg_final_trend, device=cairo_ps, width=28800,
    #        height=14400, units="px", dpi=1800, limitsize=FALSE)
    # 
    plotlist = c(plotlist, list(gg_final_trend))  
    
  }
  ggar = ggarrange(plotlist[[2]],plotlist[[3]],plotlist[[4]],nrow=4,ncol=1)
  
  ggall = plotlist[[1]]/plotlist[[2]]/plotlist[[3]]/plotlist[[4]]
  
  
  #ggsave("alltrendplotsc.eps",plot=ggar, device=cairo_ps, width=43200,
  #       height=57600, units="px", dpi=3600, limitsize=FALSE) c(14,10.5)
  ggsave("alltrendplotsc2.eps",plot=ggar, device=cairo_ps, width=14,
         height=18.5, limitsize=FALSE)
  
  #ggsave(paste0("AC1",selectstr,"_evolution.eps"),plot=gg_final_Phi, device=cairo_ps, width=8,
  #       height=10.5, units="cm", dpi=1800, limitsize=FALSE)
  getwd()
  
  margplots = list()
  ggp0 = ggplot() + theme_bw() + xlab("b") + ylab("Density") +
    theme(text = element_text(size = 20)) +
    geom_vline(aes(xintercept=0), linewidth=1.2, linetype="dotted")
  for(i in 1:17 ){
    ggd0 = data.frame(bx = rres2[[i]]$results$marginals$b[,1], by=rres2[[i]]$results$marginals$b[,2])
    ggdp = data.frame(lower=rres2[[i]]$results$summary$b$quant0.025, 
                      upper=rres2[[i]]$results$summary$b$quant0.975, zero=0,
                      bpos=rres2[[i]]$results$summary$b$prob_positive)
    ggdribbon = data.frame(x=ggd0$bx[ggd0$bx>=ggdp$lower & ggd0$bx <=ggdp$upper], 
                           ymax=ggd0$by[ggd0$bx>=ggdp$lower & ggd0$bx <=ggdp$upper],
                           ymin=numeric(sum(ggd0$bx>=ggdp$lower & ggd0$bx <=ggdp$upper)))
    
    ggmarg = ggp0 + geom_ribbon(data=ggdribbon,mapping=aes(x=x,ymin=ymin,ymax=ymax),
                                fill="red", alpha=0.3)+
      geom_line(data=ggd0,mapping=aes(x=bx,y=by)) +
      ggtitle(paste0("Event # ",i), subtitle=paste0("P(b>0) = ",round(ggdp$bpos,digits=3)))
    margplots[[i]] = ggmarg
  }
  
  ggallb = ggarrange(plotlist=margplots,ncol=4,nrow=5)
  
  print(ggallb)
  
  ggsave("ggallb.eps",plot=ggallb, device=cairo_ps, width=15,
         height=18.5, limitsize=FALSE)
  
  ####
  plotlistlin1=list()
  plotlistsq1=list()
  plotlistrw21=list()
  start1=1
  end1=6
  for(i in start1:end1){
    ggdlin = data.frame(time=df_plot_Lin[[i]]$Time, data=df_plot_Lin[[i]]$data,
                     trendL=df_plot_Lin[[i]]$trend_low,
                     trendU=df_plot_Lin[[i]]$trend_upp,
                     trend=df_plot_Lin[[i]]$trend)
    ggdsq = data.frame(time=df_plot_Quad[[i]]$Time, data=df_plot_Quad[[i]]$data,
                        trendL=df_plot_Quad[[i]]$trend_low,
                        trendU=df_plot_Quad[[i]]$trend_upp,
                        trend=df_plot_Quad[[i]]$trend)
    ggdrw2 = data.frame(time=df_plot_RW2[[i]]$Time, data=df_plot_RW2[[i]]$data,
                       trendL=df_plot_RW2[[i]]$trend_low,
                       trendU=df_plot_RW2[[i]]$trend_upp,
                       trend=df_plot_RW2[[i]]$trend)

    ggplotlin = ggplot(ggdlin, aes(x=time)) + theme_bw() + xlab("")+ylab("") +
      geom_line(aes(y=data),col="gray") + xlim(rev(range(ggdlin$time)))+
      theme(text=element_text(size=16), plot.title = element_text(size=16, hjust=0.5))+
      geom_ribbon(aes(ymin=trendL,ymax=trendU), col="red", fill="red",alpha=0.3)+
      geom_line(aes(y=trend), col="blue") +ggtitle(paste0("Event #",i))
    ggplotsq = ggplot(ggdsq, aes(x=time)) + theme_bw() + xlab("")+ylab("") +
      geom_line(aes(y=data),col="gray") + xlim(rev(range(ggdlin$time)))+
      theme(text=element_text(size=16), plot.title = element_text(size=22)) +
      geom_ribbon(aes(ymin=trendL,ymax=trendU), col="red", fill="red",alpha=0.3)+
      geom_line(aes(y=trend), col="blue")
    ggplotrw2 = ggplot(ggdrw2, aes(x=time)) + theme_bw() + xlab("")+ylab("") +
      geom_line(aes(y=data),col="gray") + xlim(rev(range(ggdlin$time)))+
      theme(text=element_text(size=16), plot.title = element_text(size=22)) +
      geom_ribbon(aes(ymin=trendL,ymax=trendU), col="red", fill="red",alpha=0.3)+
      geom_line(aes(y=trend), col="blue")
    
    if(i!=end1){
      ggplotlin = ggplotlin + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      ggplotsq = ggplotsq + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      ggplotrw2 = ggplotrw2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    
    plotlistlin1[[i]] = ggplotlin
    plotlistsq1[[i]] = ggplotsq
    plotlistrw21[[i]] = ggplotrw2
    
  }
  library(grid)
  #tlin = text_grob("Linear trend", size=20)
  titlesize=20
  xsize=18
  ysize=18
  
  gglin1 = ggarrange(plotlist=rev(plotlistlin1),nrow=1)  #+ annotate("text",x=16000,y=-37, label="test")
    ggl1 = annotate_figure(gglin1,top=text_grob("Linear trend",hjust=0,x=0,size=titlesize))
    
  ggsq1 = ggarrange(plotlist=rev(plotlistsq1),nrow=1)  
  ggsq1 = annotate_figure(ggsq1,top=text_grob("Square trend",hjust=0,x=0,size=titlesize))
  
  ggrw21 = ggarrange(plotlist=rev(plotlistrw21),nrow=1)  
  ggrw21 = annotate_figure(ggrw21,top=text_grob("RW2 trend",hjust=0,x=0,size=titlesize))
  
  ggfirst = ggarrange(ggl1,ggsq1,ggrw21, nrow=3,ncol=1)
  ggf = annotate_figure(ggfirst,bottom=text_grob("Time (yr b2k)",size=xsize), 
                        left=text_grob(expression(paste(delta^18,"O (permil)")),size=ysize,rot=90))
  
  # line 2
  plotlistlin2=list()
  plotlistsq2=list()
  plotlistrw22=list()
  start2=7
  end2=12
  for(i in start2:end2){
    ggdlin = data.frame(time=df_plot_Lin[[i]]$Time, data=df_plot_Lin[[i]]$data,
                        trendL=df_plot_Lin[[i]]$trend_low,
                        trendU=df_plot_Lin[[i]]$trend_upp,
                        trend=df_plot_Lin[[i]]$trend)
    ggdsq = data.frame(time=df_plot_Quad[[i]]$Time, data=df_plot_Quad[[i]]$data,
                       trendL=df_plot_Quad[[i]]$trend_low,
                       trendU=df_plot_Quad[[i]]$trend_upp,
                       trend=df_plot_Quad[[i]]$trend)
    ggdrw2 = data.frame(time=df_plot_RW2[[i]]$Time, data=df_plot_RW2[[i]]$data,
                        trendL=df_plot_RW2[[i]]$trend_low,
                        trendU=df_plot_RW2[[i]]$trend_upp,
                        trend=df_plot_RW2[[i]]$trend)
    
    ggplotlin = ggplot(ggdlin, aes(x=time)) + theme_bw() + xlab("")+ylab("") +
      geom_line(aes(y=data),col="gray") +
      theme(text=element_text(size=16), plot.title = element_text(size=16, hjust=0.5))+
      geom_ribbon(aes(ymin=trendL,ymax=trendU), col="red", fill="red",alpha=0.3)+
      geom_line(aes(y=trend), col="blue")+ xlim(rev(range(ggdlin$time)))+ggtitle(paste0("Event #",i))
    ggplotsq = ggplot(ggdsq, aes(x=time)) + theme_bw() + xlab("")+ylab("") +
      geom_line(aes(y=data),col="gray") +
      theme(text=element_text(size=16), plot.title = element_text(size=22)) +
      geom_ribbon(aes(ymin=trendL,ymax=trendU), col="red", fill="red",alpha=0.3)+
      geom_line(aes(y=trend), col="blue")+ xlim(rev(range(ggdlin$time)))
    ggplotrw2 = ggplot(ggdrw2, aes(x=time)) + theme_bw() + xlab("")+ylab("") +
      geom_line(aes(y=data),col="gray") +
      theme(text=element_text(size=16), plot.title = element_text(size=22)) +
      geom_ribbon(aes(ymin=trendL,ymax=trendU), col="red", fill="red",alpha=0.3)+
      geom_line(aes(y=trend), col="blue")+ xlim(rev(range(ggdlin$time)))
    
    if(i!=end2){
      
      ggplotlin = ggplotlin + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      ggplotsq = ggplotsq + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      ggplotrw2 = ggplotrw2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    
    plotlistlin2[[i-start2+1]] = ggplotlin
    plotlistsq2[[i-start2+1]] = ggplotsq
    plotlistrw22[[i-start2+1]] = ggplotrw2
    
  }
  
  gglin2 = ggarrange(plotlist=rev(plotlistlin2),nrow=1)  #+ annotate("text",x=16000,y=-37, label="test")
  ggl2 = annotate_figure(gglin2,top=text_grob("Linear trend",hjust=0,x=0,size=titlesize))
  
  ggsq2 = ggarrange(plotlist=rev(plotlistsq2),nrow=1)  
  ggsq2 = annotate_figure(ggsq2,top=text_grob("Square trend",hjust=0,x=0,size=titlesize))
  
  ggrw22 = ggarrange(plotlist=rev(plotlistrw22),nrow=1)  
  ggrw22 = annotate_figure(ggrw22,top=text_grob("RW2 trend",hjust=0,x=0,size=titlesize))
  
  ggsec = ggarrange(ggl2,ggsq2,ggrw22, nrow=3,ncol=1)
  ggs = annotate_figure(ggsec,bottom=text_grob("Time (yr b2k)",size=xsize), 
                        left=text_grob(expression(paste(delta^18,"O (permil)")),size=ysize,rot=90))
  
  
  
  
  plotlistlin3=list()
  plotlistsq3=list()
  plotlistrw23=list()
  start3=13
  end3=17
  for(i in start3:end3){
    
    ggdlin = data.frame(time=df_plot_Lin[[i]]$Time, data=df_plot_Lin[[i]]$data,
                        trendL=df_plot_Lin[[i]]$trend_low,
                        trendU=df_plot_Lin[[i]]$trend_upp,
                        trend=df_plot_Lin[[i]]$trend)
    ggdsq = data.frame(time=df_plot_Quad[[i]]$Time, data=df_plot_Quad[[i]]$data,
                       trendL=df_plot_Quad[[i]]$trend_low,
                       trendU=df_plot_Quad[[i]]$trend_upp,
                       trend=df_plot_Quad[[i]]$trend)
    ggdrw2 = data.frame(time=df_plot_RW2[[i]]$Time, data=df_plot_RW2[[i]]$data,
                        trendL=df_plot_RW2[[i]]$trend_low,
                        trendU=df_plot_RW2[[i]]$trend_upp,
                        trend=df_plot_RW2[[i]]$trend)
    
    ggplotlin = ggplot(ggdlin, aes(x=time)) + theme_bw() + xlab("")+ylab("") +
      geom_line(aes(y=data),col="gray") +
      theme(text=element_text(size=16), plot.title = element_text(size=16)) +
      geom_ribbon(aes(ymin=trendL,ymax=trendU), col="red", fill="red",alpha=0.3)+
      geom_line(aes(y=trend), col="blue")+ xlim(rev(range(ggdlin$time)))+ggtitle(paste0("Event #",i))
    ggplotsq = ggplot(ggdsq, aes(x=time)) + theme_bw() + xlab("")+ylab("") +
      geom_line(aes(y=data),col="gray") +
      theme(text=element_text(size=16), plot.title = element_text(size=22)) +
      geom_ribbon(aes(ymin=trendL,ymax=trendU), col="red", fill="red",alpha=0.3)+
      geom_line(aes(y=trend), col="blue")+ xlim(rev(range(ggdlin$time)))
    ggplotrw2 = ggplot(ggdrw2, aes(x=time)) + theme_bw() + xlab("")+ylab("") +
      geom_line(aes(y=data),col="gray") +
      theme(text=element_text(size=16), plot.title = element_text(size=22)) +
      geom_ribbon(aes(ymin=trendL,ymax=trendU), col="red", fill="red",alpha=0.3)+
      geom_line(aes(y=trend), col="blue")+ xlim(rev(range(ggdlin$time)))
    
    if(i!=end3){
      
      ggplotlin = ggplotlin + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      ggplotsq = ggplotsq + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
      ggplotrw2 = ggplotrw2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
    }
    
    plotlistlin3[[i-start3+1]] = ggplotlin
    plotlistsq3[[i-start3+1]] = ggplotsq
    plotlistrw23[[i-start3+1]] = ggplotrw2
    
  }
  
  gglin3 = ggarrange(plotlist=rev(plotlistlin3),nrow=1)  #+ annotate("text",x=16000,y=-37, label="test")
  ggl3 = annotate_figure(gglin3,top=text_grob("Linear trend",hjust=0,x=0,size=titlesize))
  
  ggsq3 = ggarrange(plotlist=rev(plotlistsq3),nrow=1)  
  ggsq3 = annotate_figure(ggsq3,top=text_grob("Square trend",hjust=0,x=0,size=titlesize))
  
  ggrw23 = ggarrange(plotlist=rev(plotlistrw23),nrow=1)  
  ggrw23 = annotate_figure(ggrw23,top=text_grob("RW2 trend",hjust=0,x=0,size=titlesize))
  
  ggthird = ggarrange(ggl3,ggsq3,ggrw23, nrow=3,ncol=1)
  ggt = annotate_figure(ggthird,bottom=text_grob("Time (yr b2k)",size=xsize), 
                        left=text_grob(expression(paste(delta^18,"O (permil)")),size=ysize,rot=90))
  ggff = ggf + theme(plot.margin = unit(c(0,0,2,0),"lines"))
  ggss = ggs + theme(plot.margin = unit(c(0,0,2,0),"lines"))
  ggtt = ggt + theme(plot.margin = unit(c(0,0,2,0),"lines"))
  
  
  ggabsall = ggarrange(ggff,ggss,ggtt,nrow=3)
  
  ggsave("ggalltrends.eps",plot=ggabsall, device=cairo_ps, width=22,
         height=22, limitsize=FALSE)
  
  ####
  
  
  margplots = list()
  ggp0 = ggplot() + theme_bw() + xlab("b") + ylab("Density") +
    theme(text = element_text(size = 20)) +
    geom_vline(aes(xintercept=0), linewidth=1.2, linetype="dotted")
  for(i in 1:17 ){
    q0.05 = inla.qmarginal(0.05, rres2[[i]]$results$marginals$b)
    ggd0 = data.frame(bx = rres2[[i]]$results$marginals$b[,1], by=rres2[[i]]$results$marginals$b[,2])
    ggdp = data.frame(lower=rres2[[i]]$results$summary$b$quant0.025, 
                      upper=rres2[[i]]$results$summary$b$quant0.975, zero=0,
                      bpos=rres2[[i]]$results$summary$b$prob_positive,
                      q0.05=q0.05)
    ggdribbon = data.frame(x=ggd0$bx[ggd0$bx>=q0.05], 
                           ymax=ggd0$by[ggd0$bx>=q0.05],
                           ymin=numeric(sum(ggd0$bx>=q0.05)))
    
    ggmarg = ggp0 + geom_ribbon(data=ggdribbon,mapping=aes(x=x,ymin=ymin,ymax=ymax),
                                fill="red", alpha=0.3)+
      geom_line(data=ggd0,mapping=aes(x=bx,y=by)) +
      ggtitle(paste0("Event # ",i), subtitle=paste0("P(b>0) = ",round(ggdp$bpos,digits=3)))
    
    margplots[[i]] = ggmarg
  }
  
  ggallb = ggarrange(plotlist=margplots,ncol=5,nrow=4)
  
  print(ggallb)
  
  ggsave("ggallb2.eps",plot=ggallb, device=cairo_ps, width=22,
         height=18, limitsize=FALSE)

}
