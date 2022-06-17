#' Bayesian analysis of early warning signals
#'
#' Fits a Bayesian hierarchical model with time-dependent correlation using INLA.
#'
#' @param data Numeric describing the observations from which early warning signals is to be detected.
#' @param forcing Numeric describing potential forcing to the system. If no forcing
#' is included, set to \code{numeric(0)} (default).
#' @param formula Formula describing the linear predictor (NOT YET IMPLEMENTED)
#' @param model Character string describing which noise model should be used.
#' Currently supports \code{"ar1"} (default) and \code{"fgn"}.
#' @param compute.mu Should the forced response be computed? Set to \code{0} if not, 
#' \code{1} if only mean and standard deviation should be computed or \code{2} if quantiles
#' should be computed (this is slower).
#' @param inla.options List giving options (such as step length and the number of threads to be used),
#' to the \code{inla}-program. Any absent options are set to the default \code{INLA} options.
#' @param print.progress boolean indicating whether or not progress should be printed to screen.
#' @param memory.true numeric describing the true memory values. Can be used in 
#' \code{\link{plot.inla.ews}}.
#' @return Returns an S3 object of class \code{inla.ews}. This object contains 
#' all results from the \code{inla}-program and associated summary statistics.
#' 
#' @examples
#' \donttest{
#' ### AR(1) simulation example ###
#' n = 300
#' sigma = 1
#' a=0.2
#' b=0.7/n
#' F0 = 3
#' sigmaf = 0.3
#' time = 1:n
#' phis = a+b*time
#' noise=ar1_timedep_sim(n,phis=phis)
#' 
#' forcing = arima.sim(model=list(ar=c(0.95)),n=n,sd=sqrt(1-0.95^2))+1:n/n
#' zz = sigmaf*(F0+forcing)
#' 
#' lambdas = phis-1
#' muvek = mu.computer(forcing,sigmaf,F0,memory=phis,model="ar1")
#' 
#' data = noise + muvek
#' 
#' object = inla.ews(data,forcing,model="ar1",compute.mu=2, print.progress=TRUE,
#'                     memory.true=phis)
#' summary(object)
#' plot(object)
#' 
#' 
#' ### fGn simulation example ###
#' n=200
#' sigma = 1.2
#' time=1:n
#' a = 0.6
#' b = 0.35/n
#' Hs = a+b*time
#' 
#' data = fgn_timedep_sim(n,Hs=Hs)
#' 
#' object = inla.ews(data,model="fgn",memory.true=Hs)
#' summary(object)
#' plot(object)
#' 
#' }
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords INLA early warning signal
#' @export
inla.ews <- function(data, forcing=numeric(0), formula=NULL, model="ar1",compute.mu=0,
                     inla.options=NULL,print.progress=FALSE,
                     memory.true=NULL){

  catch = tryCatch(attachNamespace("INLA"),error=function(x){})
  if(length(find.package("INLA",quiet=TRUE))==0){
    stop("This function requires INLA. Please install at www.R-INLA.org or by calling 'install.packages(\"INLA\", repos=c(getOption(\"repos\"), INLA=\"https://inla.r-inla-download.org/R/testing\"), dep=TRUE)' from R.")
  }
  
  time.start = Sys.time()
  
  if(print.progress){
    cat("Initiating inla.ews..\n",sep="")
  }
  tid.start = proc.time()[[3]]
  
  inla.ews.call = sys.call(which=1)
  if(sum(is.na(forcing))>0) stop("Forcing contains NA values")
  
  old.digits=getOption("digits")
  if(old.digits < 7){
    options(digits = 7)
  }
  
  n=length(data)
  time=1:n
  intercept = mean(data[1:20])
  df = data.frame(y = data-intercept,data = data,idy=1:n,time=time)

  if(print.progress){
    cat("Setting up rgeneric model..\n",sep="")
  }
  if(tolower(model) %in% c("ar1","ar(1)","1")){
    if(length(forcing)>0){
      
      rgen_model = INLA::inla.rgeneric.define(rgeneric.ews.ar1.forcing,n=n,
                                              time=df$time,forcing=forcing)
    }else{
      rgen_model = INLA::inla.rgeneric.define(rgeneric.ews.ar1,n=n,time=df$time)
    }
  }else if(tolower(model) %in% c("fgn","lrd")){
    warning("This model is considerably slower than the AR(1) process. Expect longer computational time.")
    if(length(forcing)>0){
      
      rgen_model = INLA::inla.rgeneric.define(rgeneric.ews.fgn.forcing,n=n,
                                              time=df$time,forcing=forcing)
    }else{
      rgen_model = INLA::inla.rgeneric.define(rgeneric.ews.fgn,n=n,time=df$time)
    }
    
  }
  
  # find initial guesses for hyperparameters
  if(length(forcing)>0){
    sigma_ini = (data[n]-data[1])/(forcing[n]-forcing[1])
    F0_ini = data[1]/sigma_ini-forcing[1]
    
  }else{
    
  }
  formula = y ~ -1+ f(idy, model=rgen_model) 
  
  
  
  
  if(print.progress){
    cat("Initiating inla program..\n",sep="")
  }
   r = INLA::inla(formula,data=df,family="gaussian",control.family = list(initial=12,fixed=TRUE),
            verbose=FALSE,num.threads = 1)
  # r <- do.call(INLA::inla,c(list(formula=formula,data=df,inla.options)))
   if(print.progress){
     cat("Completed INLA optimization..\n",sep="")
   }
  object = list(formula=formula,inlafit=r,.args=list(data=data,
                                                     forcing=forcing,
                                                     call = inla.ews.call,
                                                     inladata=df,
                                                     model=model,
                                                     intercept=intercept,
                                                     compute.mu=compute.mu,
                                                     inla.options=inla.options,
                                                     memory.true=memory.true))
  class(object) = "inla.ews"
  if(print.progress){
    cat("Collecting results and perform transformation to user scaling..\n",sep="")
  }
  time.start.gather = Sys.time()
  object = resultgather(object,print.progress=print.progress)
  time.gather = difftime(Sys.time(), time.start.gather,units="secs")[[1]]
  
  if(compute.mu >0){
    #intercept = object$results$fixed
    if(compute.mu==2) object=forcingmaker(object,quick=FALSE,intercept=intercept,
                                          print.progress=print.progress)
    if(compute.mu==1) object=forcingmaker(object,quick=TRUE,intercept=intercept,
                                          print.progress=print.progress)
  }
  
  time.total = difftime(Sys.time(), time.start,units="secs")[[1]]
  
  object$cpu.used = list(inla=object$inlafit$cpu.used[4],gather=time.gather,total=time.total)
  if(print.progress){
    cat("inla.ews concluded in ",time.total," seconds!\n",sep="")
  }
  options(digits=old.digits)
  return(object)
}

if(FALSE){
  n=300
  sigma1 = 1
  time=1:n
  a = 0.6
  b = 0.35/n
  Hs = a+b*time
  F0 = -3
  sigmaf=0.1
  sigma = 1.2
  #Hs = rep(0.75,n)
  
  sigmat = sigmamaker(n,sigma,Hs)
  sigmachol = chol(sigmat)
  noise = sigmachol%*%rnorm(n)
  forcing = arima.sim(model=list(ar=c(0.9)),n)+1:n/n*10
  z = sigmaf*(forcing+F0)
  
  struct = (1:n-0.5)^(Hs-3/2)
  muvek=numeric(n)
  for(i in 1:n){
    muvek[i] = rev(struct[1:i])%*%z[1:i]
  }
  y=muvek+noise
  data = y
  forcing=z#numeric(0)
  
  object = inla.ews(data,forcing,model="ar1",inla.options=inla.options,
                    print.progress=TRUE)
  ####
  if(TRUE){
    plot(object$results$summary$H$median,type="l",ylim=c(0.5,1))
    lines(object$results$summary$H$hpd0.95lower,col="red")
    lines(object$results$summary$H$hpd0.95upper,col="red")
    lines(a+b*1:n,col="gray")
  }
}


if(FALSE){
  
  n=50
  
  sigma1=1; kappa1=1/sigma1^2
  a=0.2
  b=0.7/n
  F0 = 3
  sigmaf = 2
  time = 1:n
  phis = a+b*time
  s=numeric(n)
  s[1] = rnorm(1,mean=0,sd=sigma1)
  for(i in 2:n){
    s[i] = rnorm(1, mean=phis[i]*s[i-1],sd=sigma1)
  }
  
  z = arima.sim(model=list(ar=c(0.95)),n=n,sd=sqrt(1-0.95^2))+1:n/n*7
  zz = sigmaf*(F0+z)
  
  lambdas = phis-1
  
  struct = exp(lambdas*(1:n)-0.5)
  muvek=numeric(n)
  for(i in 1:n){
    muvek[i] = rev(struct[1:i])%*%zz[1:i]
  }
  y = s + muvek
  
  data=s
  forcing=numeric(0)
  model="ar1"
  
  # radius = 30
  # phiL=acf(y[1:radius],plot=FALSE)$acf[2]
  # phiH=acf(y[(n-radius):n],plot=FALSE)$acf[2]
  # bhat = (phiH-phiL)/n
  # ahat = phiL-bhat
  # varhat = (var(y[1:radius])*(1-phiL^2)+var(y[(n-radius):n])*(1-phiH^2))/2
  # 
  # theta_sigma = log(1/varhat)
  # theta_b = log((1+bhat)/(1-bhat))
  # theta_a = log( (min(b*time)+ahat)/(1-max(b*time)-ahat) )
  # initheta = c(theta_sigma,theta_b,theta_a)
  # if(length(forcing)>0) initheta=c(initheta,c(0,0))
  # #inla.options=list(control.mode=list(theta=initheta,restart=TRUE))
  # 
  # object = inla.ews(data,forcing,model=model,inla.options=inla.options,
  #                   modecontrol = list(theta=initheta,restart=TRUE))
  object = inla.ews(data,forcing,model=model,inla.options=inla.options)
  ####
  if(TRUE){
    plot(object$results$summary$phi$median,type="l",ylim=c(0,1))
    lines(object$results$summary$phi$hpd0.95lower,col="red")
    lines(object$results$summary$phi$hpd0.95upper,col="red")
    lines(a+b*1:n,col="gray")
    
  }
  stop("stopper her")
  
  
  par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
  plot(y,type="l",xlab="Time",ylab="Observation",main=paste0("varying-memory: a=",a,", b*n=",b*n))
  lines(muvek,col="red")
  
  df = data.frame(idy=1:n,y=y,time=time)
  library(INLA)
  
  rgen_model = inla.rgeneric.define(rgeneric.ar1.varphi.forcing,n=n,time=df$time,forcing=z)
  formula = y~-1+f(idy,model=rgen_model)
  r = inla(formula,data=df,control.family = list(initial=12,fixed=TRUE))
  summary(r)
  
  
  truevals = list(a=a,b=b,phis=phis,sigma=sigma,sigmaf=sigmaf,F0=F0)
  resultfunc(r,time,truevals=truevals)
  
}