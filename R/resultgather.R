#' Gather results from INLA fit
#' This function gathers and formats results from INLA fit.
#' 
#' @param object S3 object of type \code{inla.ews} which includes result from 
#' \code{inla}-program.
#' @param print.progress boolean indicating if progress should be printed to screen.
#' @return Returns the same S3 object of class \code{inla.ews} as included in the 
#' input arguments, but appends summary statistics and Monte Carlo simulations in 
#' \code{object\$results}.
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords INLA early warning signal summary
#' @importFrom stats density dnorm
#' @importFrom matrixStats rowMedians rowSds
resultgather <- function(object,print.progress){
  n = length(object$.args$data)
  r = object$inlafit
  df = object$.args$inladata
  time=df$time_normalized
  sigma_est=INLA::inla.emarginal(function(x)1/sqrt(exp(x)),object$inlafit$marginals.hyperpar$`Theta1 for idy`)
  sigma_marg=INLA::inla.tmarginal(function(x)1/sqrt(exp(x)),object$inlafit$marginals.hyperpar$`Theta1 for idy`)
  
  rekke = diff(range(df$time_normalized))
  b_est=INLA::inla.emarginal(function(x) -1/rekke+2/rekke*1/(1+exp(-x)) ,r$marginals.hyperpar$`Theta2 for idy` )
  b_marg=INLA::inla.tmarginal(function(x) -1/rekke+2/rekke*1/(1+exp(-x)) ,r$marginals.hyperpar$`Theta2 for idy`)
  b_positive_prob = 1-INLA::inla.pmarginal(0,b_marg)
  b_zmarg = INLA::inla.zmarginal(b_marg,silent=TRUE)
  b_zmarg$mode = INLA::inla.mmarginal(b_marg)
  b_low = b_zmarg$quant0.025
  b_high= b_zmarg$quant0.975
  nsims = 30000
  hypersamples = INLA::inla.hyperpar.sample(nsims,r)
  
  if(object$.args$model %in% c("ar1","ar(1)","1", "ar1g")){
    if(length(object$.args$forcing)>0){
      sigmaf_marg = INLA::inla.tmarginal(function(x)1/sqrt(exp(x)),object$inlafit$marginals.hyperpar$`Theta4 for idy`)
      F0_marg = object$inlafit$marginals.hyperpar$`Theta5 for idy`
    }
    
    b_sims = -1/rekke+2/rekke*1/(1+exp(-hypersamples[,2]))
    #b = -1/r+2/r*1/(1+exp(-theta[2]))
    a_sims = numeric(nsims)
    if(print.progress){
      cat("Producing ",nsims," monte carlo samples for 'a' hyperparameter and 'phi' vector...",sep="")
    }
      
    phi_sims = matrix(NA,nrow=n,ncol=nsims)
    for(i in 1:nsims){
      low=min(b_sims[i]*df$time_normalized[1],b_sims[i]*df$time_normalized[n])
      high=max(b_sims[i]*df$time_normalized[1],b_sims[i]*df$time_normalized[n])
      a_sims[i] = -low + (1-high+low)/(1+exp(-hypersamples[i,3]))
      phi_sims[,i] = a_sims[i]+b_sims[i]*df$time_normalized
    }
    if(print.progress){
      cat(" completed!\n",sep="")
    }
    a_marg = density(a_sims); a_marg=data.frame(x=a_marg$x,y=a_marg$y)
    a_zmarg = INLA::inla.zmarginal(a_marg,silent=TRUE)
    a_zmarg$mode = INLA::inla.mmarginal(a_marg)
    
    if(print.progress){
      cat("Computing summary statistics for each 'phi'...",sep="")
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
    }
    if(print.progress){
      cat(" completed!\n",sep="")
    }
  }else if(object$.args$model %in% c("fgn","lrd")){
    if(length(object$.args$forcing)>0){
      sigmaf_marg = INLA::inla.tmarginal(function(x)1/sqrt(exp(x)),object$inlafit$marginals.hyperpar$`Theta4 for idy`)
      F0_marg = object$inlafit$marginals.hyperpar$`Theta5 for idy`
    }
    bmax = 0.5/(time[n]-time[1])
    bmin = -bmax
    b_marg = INLA::inla.tmarginal(function(x)bmin+(bmax-bmin)/(1+exp(-x)),r$marginals.hyperpar$`Theta2 for idy`)
    bsims = bmin+(bmax-bmin)/(1+exp(-hypersamples[,2]))
    a_sims = numeric(nsims)
    if(print.progress){
      cat("Producing ",nsims," monte carlo samples for 'a' hyperparameter and 'H' vector...",sep="")
    }
    H_sims = matrix(NA,nrow=n,ncol=nsims)
    for(i in 1:nsims){
      low = min(bsims[i]*time[1],bsims[i]*time[n])
      high = max(bsims[i]*time[1],bsims[i]*time[n])
      amin = 0.5-low
      amax = 1-high
      a_sims[i] = amin + (amax-amin)/(1+exp(-hypersamples[i,3]))
      H_sims[,i] = a_sims[i]+bsims[i]*time
    }
    if(print.progress){
      cat(" completed!\n",sep="")
    }
    a_marg = density(a_sims); a_marg=data.frame(x=a_marg$x,y=a_marg$y)
    a_zmarg = INLA::inla.zmarginal(a_marg,silent=TRUE)
    a_zmarg$mode = INLA::inla.mmarginal(a_marg)
    
    if(print.progress){
      cat("Computing summary statistics for each 'H'...",sep="")
    }
    Hmean = rowMeans(H_sims)
    Hmedian = matrixStats::rowMedians(H_sims)
    Hsds = matrixStats::rowSds(H_sims)
    Hlower = numeric(n)
    Hupper = numeric(n)
    Hqlower = numeric(n)
    Hqmid = numeric(n)
    Hqupper = numeric(n)
    Hmode = numeric(n)
    for(i in 1:n){
      dens =density(H_sims[i,]); dens = data.frame(x=dens$x,y=dens$y)
      zmarg = INLA::inla.zmarginal(dens,silent=TRUE)
      Hmode[i] = INLA::inla.mmarginal(dens)
      hpds = INLA::inla.hpdmarginal(0.95,dens)
      Hlower[i]=hpds[1]; Hupper[i]=hpds[2]
      Hqlower[i] = zmarg$quant0.025
      Hqmid[i] = zmarg$quant0.5
      Hqupper[i] = zmarg$quant0.975
    }
    if(print.progress){
      cat(" completed!\n",sep="")
    }
    
  }
  
  
  if(print.progress){
    cat("Storing everything in object$results..\n",sep="")
  }
  object$results = list(marginals = list(a=a_marg, b=b_marg, sigma=sigma_marg))
  object$results$summary = list(a = a_zmarg,#list(INLA::inla.zmarginal(a_marg,silent=TRUE))[[1]],
                                b = b_zmarg,#list(INLA::inla.zmarginal(b_marg,silent=TRUE))[[1]],
                                sigma = list(INLA::inla.zmarginal(sigma_marg,silent=TRUE))[[1]])
  object$results$summary$b$prob_positive=b_positive_prob
  if(tolower(object$.args$model) %in% c("ar1","ar(1)","1","ar1g")){
    object$results$summary$phi = list(mean=phi_means,median=phi_median,sd=phi_sd,
                                      q0.025=phi_qlower,q0.5=phi_qmid,q0.975=phi_qupper,
                                      hpd0.95lower=phi_lower,hpd0.95upper=phi_upper)
    object$results$simulations = list(phi_sims = phi_sims)
    
  }else if(tolower(object$.args$model) %in% c("fgn","lrd")){
    object$results$summary$H = list(mean=Hmean,median=Hmedian,sd=Hsds,
                                    q0.025=Hqlower,q0.5=Hqmid,q0.975=Hqupper,
                                    hpd0.95lower=Hlower,hpd0.95upper=Hupper)
    object$results$simulations = list(H_sims = H_sims)
    # if(length(object$.args$forcing)>0){
    #   object$results$marginals$sigmaf = sigmaf_marg
    #   object$results$marginals$F0 = F0_marg
    #   
    #   object$results$summary$sigmaf = INLA::inla.zmarginal(sigmaf_marg,silent=TRUE)
    #   object$results$summary$sigmaf = INLA::inla.zmarginal(sigmaf_marg,silent=TRUE)
    #   object$results$summary$F0 = INLA::inla.zmarginal(F0_marg,silent=TRUE)
    # }
  }
  if(length(object$.args$forcing)>0){
    object$results$marginals$sigmaf = sigmaf_marg
    object$results$marginals$F0 = F0_marg
    
    object$results$summary$sigmaf = INLA::inla.zmarginal(sigmaf_marg,silent=TRUE)
    object$results$summary$sigmaf$mode = INLA::inla.mmarginal(sigmaf_marg)
    object$results$summary$F0 = INLA::inla.zmarginal(F0_marg,silent=TRUE)
    object$results$summary$F0$mode = INLA::inla.mmarginal(F0_marg)
  }
  
  if(!is.null(object$inlafit$summary.fixed)){
    object$results$fixed = object$inlafit$summary.fixed     
  }else{
    object$results$fixed = mean(object$.args$data[1:20])
  }
  
  return(object)
}