#' Gather results from INLA fit
#' This function gathers and formats results from INLA fit, assuming a nested time dependent AR(1) model.
#' 
#' @param object S3 object of type \code{inla.ews} which includes result from 
#' \code{inla}-program.
#' @param nsims Integer stating the number of simulations used in Monte Carlo estimation of parameters
#' not obtained directly from INLA. Default 10000. 
#' @param print.progress boolean indicating if progress should be printed to screen.
#' @return Returns the same S3 object of class \code{inla.ews} as included in the 
#' input arguments, but appends summary statistics and Monte Carlo simulations in 
#' \code{object\$results}.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords INLA early warning signal summary 
#' @importFrom stats density dnorm sd
#' @importFrom matrixStats rowMedians rowSds rowQuantiles rowMeans2
resultgather_ar2 = function(object,nsims=10000,print.progress=FALSE){
  n = nrow(object$.args$inladata)
  #n = length(object$.args$data)
  r = object$inlafit
  df = object$.args$inladata
  time=object$.args$time_normalized
  time_normalized=object$.args$time_normalized
  
  #nrand = length(r$summary.random)-1
  if(length(object$.args$forcing)>0){
    nextrahyps = length(r$summary.hyperpar$mean)-6
  }else{
    nextrahyps = length(r$summary.hyperpar$mean)-5
  }
  
  
  sigma_est=INLA::inla.emarginal(function(x)1/sqrt(exp(x)),r$marginals.hyperpar$`Theta3 for idy`)
  sigma_marg=INLA::inla.tmarginal(function(x)1/sqrt(exp(x)),r$marginals.hyperpar$`Theta3 for idy`)
  
  rekke = diff(range(object$.args$time_normalized))
  b_phi_est=INLA::inla.emarginal(function(x) -1/rekke+2/rekke*1/(1+exp(-x)) ,r$marginals.hyperpar$`Theta1 for idy` )
  b_phi_marg=INLA::inla.tmarginal(function(x) -1/rekke+2/rekke*1/(1+exp(-x)) ,r$marginals.hyperpar$`Theta1 for idy`)
  b_phi_positive_prob = 1-INLA::inla.pmarginal(0,b_phi_marg)
  b_phi_zmarg = INLA::inla.zmarginal(b_phi_marg,silent=TRUE)
  b_phi_zmarg$mode = INLA::inla.mmarginal(b_phi_marg)
  b_phi_low = b_phi_zmarg$quant0.025
  b_phi_high= b_phi_zmarg$quant0.975
  
  
  b_rho_est=INLA::inla.emarginal(function(x) -1/rekke+2/rekke*1/(1+exp(-x)) ,r$marginals.hyperpar$`Theta4 for idy` )
  b_rho_marg=INLA::inla.tmarginal(function(x) -1/rekke+2/rekke*1/(1+exp(-x)) ,r$marginals.hyperpar$`Theta4 for idy`)
  b_rho_positive_prob = 1-INLA::inla.pmarginal(0,b_rho_marg)
  b_rho_zmarg = INLA::inla.zmarginal(b_rho_marg,silent=TRUE)
  b_rho_zmarg$mode = INLA::inla.mmarginal(b_rho_marg)
  b_rho_low = b_rho_zmarg$quant0.025
  b_rho_high= b_rho_zmarg$quant0.975
  
  
  
  
  
  hypersamples = INLA::inla.hyperpar.sample(nsims,r)
  
  if(object$.args$model %in% c("ar2","ar(2)","2")){
    if(length(object$.args$forcing)>0){
      sigmaf_marg = INLA::inla.tmarginal(function(x)1/sqrt(exp(x)),r$marginals.hyperpar$`Theta6 for idy`)
      #F0_marg = object$inlafit$marginals.hyperpar$`Theta5 for idy`
    }
    
    b_phi_sims = -1/rekke+2/rekke*1/(1+exp(-hypersamples[,1+nextrahyps]))
    b_rho_sims = -1/rekke+2/rekke*1/(1+exp(-hypersamples[,4+nextrahyps]))
    #b = -1/r+2/r*1/(1+exp(-theta[2]))
    a_phi_sims = numeric(nsims)
    a_rho_sims = numeric(nsims)
    
    if(print.progress){
      cat("Producing ",nsims," monte carlo samples for 'a' hyperparameters and 'phi,rho' vectors...",sep="")
    }
    
    phi_sims = matrix(NA,nrow=n,ncol=nsims)
    phi0_sims = matrix(NA,nrow=n,ncol=nsims)
    lambda_sims = matrix(NA,nrow=n,ncol=nsims)
    
    rho_sims = matrix(NA,nrow=n,ncol=nsims)
    
    
    for(i in 1:nsims){
      low_phi=min(b_phi_sims[i]*time[1],b_phi_sims[i]*time[n])
      high_phi=max(b_phi_sims[i]*time[1],b_phi_sims[i]*time[n])
      a_phi_sims[i] = -low_phi + (1-high_phi+low_phi)/(1+exp(-hypersamples[i,2+nextrahyps]))
      lambda_sims[,i] = -log(a_phi_sims[i]+b_phi_sims[i]*time)
      #phi_sims[,i] = a_sims[i]+b_sims[i]*df$time_normalized
      cc = 1/(length(time_normalized-1))
      
      phi_sims[,i] = exp(-lambda_sims[,i]*c(0,diff(time_normalized))/cc)
      phi0_sims[,i] = a_phi_sims[i]+b_phi_sims[i]*time_normalized
      
      
      low_rho=min(b_rho_sims[i]*time[1],b_rho_sims[i]*time[n])
      high_rho=max(b_rho_sims[i]*time[1],b_rho_sims[i]*time[n])
      a_rho_sims[i] = -low_rho + (1-high_rho+low_rho)/(1+exp(-hypersamples[i,5+nextrahyps]))
      cc = 1/(length(time_normalized-1))
      rho_sims[,i] = a_rho_sims[i]+b_rho_sims[i]*time_normalized
      
      
    }
    if(print.progress){
      cat(" completed!\n",sep="")
    }
    a_phi_marg = density(a_phi_sims); a_phi_marg=data.frame(x=a_phi_marg$x,y=a_phi_marg$y)
    a_phi_zmarg = INLA::inla.zmarginal(a_phi_marg,silent=TRUE)
    a_phi_zmarg$mode = INLA::inla.mmarginal(a_phi_marg)
    
    a_rho_marg = density(a_rho_sims); a_rho_marg=data.frame(x=a_rho_marg$x,y=a_rho_marg$y)
    a_rho_zmarg = INLA::inla.zmarginal(a_rho_marg,silent=TRUE)
    a_rho_zmarg$mode = INLA::inla.mmarginal(a_rho_marg)
    
    if(print.progress){
      cat("Computing summary statistics for each 'phi and rho'...",sep="")
    }
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
    rho_means = rowMeans(rho_sims)
    rho_median = numeric(n)
    rho_sd = numeric(n)
    rho_lower = numeric(n)
    rho_upper = numeric(n)
    rho_qlower = numeric(n)
    rho_qmid = numeric(n)
    rho_qupper = numeric(n)
    rho_mode = numeric(n)
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
      
      dens = density(rho_sims[i,])
      dens = data.frame(x=dens$x,y=dens$y)
      zmarg = INLA::inla.zmarginal(dens,silent=TRUE)
      rho_median[i]=zmarg$quant0.5; rho_sd[i] = zmarg$sd
      hpds = INLA::inla.hpdmarginal(0.95,dens)
      rho_lower[i]=hpds[1]; rho_upper[i]=hpds[2]
      rho_qlower[i]=zmarg$quant0.025
      rho_qmid[i]=zmarg$quant0.5
      rho_qupper[i]=zmarg$quant0.975
      rho_mode[i] = INLA::inla.mmarginal(dens)
      
    }
    phi_means[1] = NA
    phi_median[1] = NA
    phi_lower[1] = NA
    phi_upper[1] = NA
    phi_qmid[1] = NA
    phi_qlower[1] = NA
    phi_qupper[1] = NA
    phi_mode[1] = NA
    
    rho_means[1] = NA
    rho_median[1] = NA
    rho_lower[1] = NA
    rho_upper[1] = NA
    rho_qmid[1] = NA
    rho_qlower[1] = NA
    rho_qupper[1] = NA
    rho_mode[1] = NA
    if(print.progress){
      cat(" completed!\n",sep="")
    }
  }
  
  if(print.progress){
    cat("Sampling fixed and random effects...",sep="")
  }
  # compute trends
  
  
  trendmean=numeric(n)
  trendupper=numeric(n)
  trendlower=numeric(n)
  trendsamps = matrix(0,nrow=n,ncol=nsims)
  if(length(r$summary.fixed$mean)>0){
    fixnames = rownames(r$summary.fixed)
    
    nfix = length(fixnames)
    nrand = length(r$summary.random)
    postsamps = INLA::inla.posterior.sample(nsims,r)
    fixsamps = matrix(0,nrow=nfix,ncol=nsims)
    
    if(nrand>1){
      randsamps = matrix(NA,nrow=n,ncol=nsims)
    }
    for(i in 1:nsims){
      latents = postsamps[[i]]$latent
      ndpred = length(object$inlafit$summary.linear.predictor$mean)
      ndrand = 0
      for(ss in 1:length(object$inlafit$summary.random)){
        ndrand = ndrand + length(object$inlafit$summary.random[[ss]]$mean)
      }
      fixsamps[,i] = latents[ndpred+ndrand+1:nfix]
      trendsamps[,i] = numeric(n)
      for(k in 1:nfix){
        if(fixnames[k] == "(Intercept)"){
          trendsamps[,i] = trendsamps[,i] + rep(fixsamps[k,i],n)
        }else{
          trendsamps[,i] = trendsamps[,i] + df[[fixnames[k]]]*fixsamps[k,i]
        }
      }
      if(nrand>1){
        for(k in 1:(nrand-1)){
          trendsamps[,i] = trendsamps[,i] + latents[n*k+1:n]
        }
      }
      
    }
    
    # for(i in 1:n){
    #   dens0 = density(trendsamps[i,]);dens = data.frame(x=dens0$x,y=dens0$y)
    #   zm = INLA::inla.zmarginal(dens,silent=TRUE)
    #   trendmean[i] = zm$mean
    #   trendupper[i] = zm$quant0.025
    #   trendlower[i] = zm$quant0.975
    # }
    trendmean = rowMeans2(trendsamps)
    trendlower = rowQuantiles(trendsamps, probs=0.025)
    trendupper = rowQuantiles(trendsamps, probs=0.975)
    
  }
  if(print.progress){
    cat(" completed!\n",sep="")
  }
  
  
  
  
  if(length(object$.args$forcing)>0){
    if(print.progress){
      cat("Simulating forcing response...",sep="")
    }
    Fmean = numeric(n)
    Flower = numeric(n)
    Fupper = numeric(n)
    muveksamps = matrix(0,nrow=n,ncol=nsims)
    kappa_samp = exp(hypersamples[,6+nextrahyps])
    forcing=object$.args$forcing
    #F0samp = hypersamples[,5+nextrahyps]
    for(s in 1:nsims){
      
      sfsamp = 1/sqrt(kappa_samp[s])
      pars = c(sfsamp, b_phi_sims[s], a_phi_sims[s], kappa_samp[s], 0)#c(kappa_eps,b,a,kappa_f,F0)
      
      muveksamps[,s] = mucomputer(pars=pars,
                                forcing = forcing,#object$.args$forcing,
                               time_norm = time_normalized, as.theta=FALSE)
      #zz = (object$.args$forcing+F0samp[i])*sfsamp[i]
      
      #struktur = exp(-lambda_sims[s]*time_normalized)
      
      #for(i in 1:n){
      #  muveksamps[i,s] = rev(struktur[1:i])%*%zz[1:i]
      #}
      # lambda_est = -log(a_phi_sims[s] + b_phi_sims[s]*time_normalized)
      # zz_est = (forcing)/sqrt(kappa_samp)
      # for(k in 1:n){
      #   for(i in 1:k){
      #     muveksamps[k,s] = nu_est[k] + zz_est[i]*exp(-lambda_est[k]*(time_normalized[k]-time_normalized[i]))
      #   }
      # }
      # muveksamps[,s]= (1/(sqrt(2*lambda_est)))*muveksamps[,s]
      
      
      
      
    }
    # for(i in 1:n){
    #   dens0 = density(muveksamps[i,]); dens=data.frame(x=dens0$x,y=dens0$y)
    #   zm = INLA::inla.zmarginal(dens,silent=TRUE)
    #   Fmean[i] = zm$mean
    #   Flower[i] = zm$quant0.025
    #   Fupper[i] = zm$quant0.975
    # }
    Fmean = rowMeans2(muveksamps)
    Flower = rowQuantiles(muveksamps, probs=0.025)
    Fupper = rowQuantiles(muveksamps, probs=0.975)
    
    alltrendsamps = trendsamps + muveksamps
    # alltrendmean = numeric(n)
    # alltrendupper = numeric(n)
    # alltrendlower = numeric(n)
    if(diff(range(alltrendsamps))>0){
      # for(i in 1:n){
      #   dens0 = density(alltrendsamps[i,]); dens=data.frame(x=dens0$x,y=dens0$y)
      #   zm = INLA::inla.zmarginal(dens,silent=TRUE)
      #   alltrendmean[i] = zm$mean
      #   alltrendlower[i] = zm$quant0.025
      #   alltrendupper[i] = zm$quant0.975
      # }
      alltrendmean = rowMeans2(alltrendsamps)
      alltrendlower = rowQuantiles(alltrendsamps, probs=0.025)
      alltrendupper = rowQuantiles(alltrendsamps, probs=0.975)
    }
  }else{
    alltrendmean=trendmean
    alltrendupper=trendupper
    alltrendlower=trendlower
  }
  
  
  
  
  
  if(print.progress){
    cat("completed!\n",sep="")
  }
  
  
  if(print.progress){
    cat("Storing everything in object$results..\n",sep="")
  }
  object$results = list(marginals = list(a_phi=as.data.frame(a_phi_marg), b_phi=as.data.frame(b_phi_marg),
                                         a_rho=as.data.frame(a_rho_marg), b_rho=as.data.frame(b_rho_marg), 
                                         sigma=as.data.frame(sigma_marg)))
  object$results$summary = list(a_phi = a_phi_zmarg,#list(INLA::inla.zmarginal(a_marg,silent=TRUE))[[1]],
                                b_phi = b_phi_zmarg,#list(INLA::inla.zmarginal(b_marg,silent=TRUE))[[1]],
                                a_rho = a_rho_zmarg,#list(INLA::inla.zmarginal(a_marg,silent=TRUE))[[1]],
                                b_rho = b_rho_zmarg,#list(INLA::inla.zmarginal(b_marg,silent=TRUE))[[1]],
                                sigma = list(INLA::inla.zmarginal(sigma_marg,silent=TRUE))[[1]])
  object$results$summary$b_phi$prob_positive=b_phi_positive_prob
  object$results$summary$b_rho$prob_positive=b_rho_positive_prob
  
  object$results$summary$phi = list(mean=phi_means,median=phi_median,sd=phi_sd,
                                    q0.025=phi_qlower,q0.5=phi_qmid,q0.975=phi_qupper,
                                    hpd0.95lower=phi_lower,hpd0.95upper=phi_upper)
  object$results$summary$phi0 = list(mean=phi0_means,median=phi0_median,sd=phi0_sd,
                                     q0.025=phi0_qlower,q0.5=phi0_qmid,q0.975=phi0_qupper,
                                     hpd0.95lower=phi0_lower,hpd0.95upper=phi0_upper)
  object$results$summary$rho = list(mean=rho_means,median=rho_median,sd=rho_sd,
                                    q0.025=rho_qlower,q0.5=rho_qmid,q0.975=rho_qupper,
                                    hpd0.95lower=rho_lower,hpd0.95upper=rho_upper)
  object$results$simulations = list(phi_sims = phi_sims, phi0_sims=phi0_sims,rho_sims = rho_sims)
  
  object$results$summary$trend = list(mean=trendmean, quant0.025=trendlower, quant0.975=trendupper)
  object$results$summary$alltrend = list(mean=alltrendmean, quant0.025=alltrendlower, quant0.975=alltrendupper)
  if(length(object$.args$forcing)>0){
    object$results$marginals$sigmaf = as.data.frame(sigmaf_marg)
   # object$results$marginals$F0 = F0_marg
    
    object$results$summary$sigmaf = INLA::inla.zmarginal(sigmaf_marg,silent=TRUE)
    object$results$summary$sigmaf$mode = INLA::inla.mmarginal(sigmaf_marg)
    #object$results$summary$F0 = INLA::inla.zmarginal(F0_marg,silent=TRUE)
    #object$results$summary$F0$mode = INLA::inla.mmarginal(F0_marg)
    
    object$results$summary$Fresponse = list(mean=Fmean, quant0.025=Flower, quant0.975=Fupper)
  }
  
  # if(!is.null(object$inlafit$summary.fixed)){
  #   object$results$fixed = object$inlafit$summary.fixed     
  # }else{
  #   object$results$fixed = mean(object$.args$data[1:20])
  # }
  
  return(object)
}
