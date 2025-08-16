#' Bayesian analysis of early warning signals
#'
#' Fits a Bayesian hierarchical model with time-dependent correlation using INLA.
#' @param data Numeric describing the observations from which early warning signals is to be detected.
#' @param forcing Numeric describing potential forcing to the system. If no forcing
#' is included, set to \code{numeric(0)} (default).
#' @param formula Formula describing the linear predictor
#' @param model Character string describing which noise model should be used.
#' Currently only supports \code{"ar1"} (default).
#' @param compute.mu Should the forced response be computed? Set to \code{0} if not, 
#' \code{1} if only mean and standard deviation should be computed or \code{2} if quantiles
#' should be computed (this is slower).
#' @param timesteps Numeric vector describing time steps of input variables (if steps are not equidistant).
#' @param log.prior Function expressing the logarithm of the joint prior for the 
#' hyperparameters in internal scaling. \code{theta[1]} corresponds to \code{log(kappa)}, 
#' \code{theta[2]} corresponds to the numerically stable parameterization of \code{b} 
#' and \code{theta[3]} corresponds to the parameterization of \code{a}. see examples 
#' for a demonstration. If \code{log.prior=NULL}, default values will be selected. These are
#' the gamma distribution on \code{kappa} with \code{shape=1} and rate=\code{0.01}, and uniform
#' priors on \code{a} and \code{b}. 
#' @param stepsize Numeric stating the step size used in the optimization procedure in 
#' the INLA program. If convergence cannot be achieved, then changing this will sometimes help. If 
#' an array of length>=2 is given, INLA will rerun, starting at the previous optima, with given stepsizes. 
#' Default value is \code{stepsize=0.01}.
#' @param num.threads Integer stating how many cores should be used by the INLA program. 
#' For rgeneric models (as is used here) stability is sometimes improved if \code{num.threads=1}, 
#' which is the default value here. Increasing this might improve runtime. 
#' @param nsims Integer stating the number of simulations used in Monte Carlo estimation of parameters
#' not obtained directly from INLA. Default 10000. 
#' @param do.cgeneric Boolean stating whether or not the model should be specified using
#' \code{cgeneric} rather than \code{rgeneric} in cases where this has been implemented. 
#' Gives improved speed, but bugs are more likely. PC priors are not implemented 
#' for this model, and any adjustments to the source code will be more difficult.
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
#' set.seed(123)
#' n = 1000
#' sigma = 1
#' a=0.3
#' b=0.5
#' time = 1:n
#' phis = a+b*seq(0,1,length.out=n)
#' data=ar1_timedep_sim(n,sigma=sigma,a=a,b=b)
#' 
#' 
#' object = inla.ews(data,model="ar1", print.progress=TRUE,
#'                     memory.true=phis)
#' summary(object)
#' plot(object)
#' 
#' ## customize prior ##
#' my.log.prior = function(theta){ #numeric of length 3 (log(kappa), theta_b, theta_a)
#'   lprior = dnorm(theta[1], sd=1, log=TRUE) + 
#'            dnorm(theta[2], sd=1, log=TRUE) + 
#'            dnorm(theta[3], sd=1, log=TRUE)
#'   return (lprior)
#' }
#' object2 = inla.ews(data, log.prior=my.log.prior,model="ar1",compute.mu=FALSE, print.progress=TRUE,
#'                     memory.true=phis)
#' 
#' summary(object2)
#' }
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @keywords INLA early warning signal
#' @export
#' @importFrom stats as.formula
inla.ews <- function(data, forcing=numeric(0), formula=NULL, model="ar1",compute.mu=0,
                     timesteps=NULL, log.prior=NULL, stepsize=0.005,num.threads=1, nsims=10000,
                     do.cgeneric=FALSE,
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
  
  inla.options = set.options(inla.options, inla.options.default()) #fill missing default settings
  
  
  inla.options$num.threads=num.threads
  
  if(is.null(inla.options$control.compute)){
    inla.options$control.compute = list(config=TRUE)
  }
  
  
  
  if(is.null(dim(data))){
    n = length(data)
  }else{
    n = nrow(data)
  }
  
  #intercept = mean(data[1:20])
  if(!is.null(timesteps)){
    time=timesteps
  }else{
    timesteps=1:n
    time=timesteps
  }
  
  #time_normalized = time-time[1]; time_normalized=time_normalized/time_normalized[n]
  time_normalized = (time-min(time))/max(time-min(time))
  if(time[1]>time[n]){
    flippedtime = TRUE
    time_normalized = 1-time_normalized
  }else{
    flippedtime=FALSE
  }
  
  if(is.null(dim(data))){
    df = data.frame(y = data,data = data,idy=1:n,time_normalized=time_normalized,
                    time=time)
  }else{
    df=data
  }
  
  
  
  if(print.progress){
    cat("Setting up rgeneric model..\n",sep="")
  }
  if(tolower(model) %in% c("ar1","ar(1)","1")){
    if(length(forcing)>0){
      require(INLA.ews, quietly=TRUE)
      
      if(!do.cgeneric){
        rgen_model = INLA::inla.rgeneric.define(rgeneric.ar1.forcing,n=n,
                                                time=time_normalized,
                                                my.log.prior=log.prior,
                                                forcing=forcing
        )
      }else{
        cscript.path=c()
        for(i in 1:length(.libPaths())){
          if("INLA.ews" %in% list.files(.libPaths()[i])){
            cscript.path=paste0(.libPaths()[i],"/INLA.ews/libs/INLA.ews.so")
          }
        }
        if(length(cscript.path)==0){
          stop("Could not find package directory, please make sure that INLA.climate is installed within one of the libraries displayed by '.libPaths()'.")
        }
        rgen_model <- INLA::inla.cgeneric.define(model = "inla_cgeneric_timedep_forcing",
                                                 #shlib = "src/cgeneric.so", 
                                                 #shlib = "libs/INLA.ews.so", 
                                                 shlib = cscript.path, 
                                                 n = n, debug=FALSE,
                                                 time=as.numeric(time_normalized), 
                                                 forcing=as.numeric(forcing)
        )
      }
      
    }else{
      rgen_model = INLA::inla.rgeneric.define(rgeneric.ar1,n=n,
                                              time=time_normalized,
                                              my.log.prior=log.prior)
    }
  }else if(tolower(model) %in% c("ar2","ar(2)","2")){
    if(length(forcing)>0){
      require(INLA.ews, quietly=TRUE)
      
      if(!do.cgeneric){
        # rgen_model = INLA::inla.rgeneric.define(rgeneric.ar2.forcing,n=n,
        #                                         time=time_normalized,
        #                                         my.log.prior=log.prior,
        #                                         forcing=forcing
        # )
        
        rgen_model = INLA::inla.rgeneric.define(rgeneric.ar2.forcing,n=n,
                                                #time=seq(from=0,to=1,length.out=n),
                                                time=time_normalized, #sometimes different results if the line above is uncommented instead
                                                my.log.prior=log.prior,
                                                forcing=forcing
        )
        
      }else{
        stop("cgeneric is not yet implemented for the AR(2) model.")
      }
      
    }else{
      rgen_model = INLA::inla.rgeneric.define(rgeneric.ar2,n=n,
                                              time=time_normalized,
                                              my.log.prior=log.prior)
    }
  }
  
  
  if(is.null(formula)){
    fstr0 = "y ~ 1"
    formula2 = y ~ 1+ f(idy, model=rgen_model)
    if(is.null(dim(data))){
      df = data.frame(y=data,idy=1:n)
    }else{
      df$idy = 1:n
    }
  }else{
    fstr0 = deparse(formula)
    fstr = paste0(fstr0,"+ f(idy, model=rgen_model) ")
    formula2 = as.formula(fstr)
    if(is.null(dim(data))){
      df = data.frame(y=data,idy=1:n)
    }else{
      df$idy = 1:n
    }
  }
  
  
  
  # r = INLA::inla(formula,data=df,family="gaussian",control.family = list(initial=12,fixed=TRUE),
  #          verbose=FALSE,num.threads = 1)
  
  nstepsizes = length(stepsize)
  for(i in 1:nstepsizes){
    if(print.progress && nstepsizes==1){
      cat("Initiating inla program..\n",sep="")
      
    }else if(print.progress && nstepsizes>=1){
      cat("Initiating inla program. Iteration ",i," of ",nstepsizes,": stepsize = ",stepsize[i],".\n",sep="")
      
    }
    if(i >= 2){
      inla.options$control.mode$theta = r$summary.hyperpar$mode
      inla.options$control.mode$result = r #use previous theta and x-mode
      inla.options$control.mode$restart = TRUE
      
    }
    inla.options$control.inla$h = stepsize[i]
    
    
    r <- tryCatch(
      do.call(INLA::inla, c(list(formula = formula2,data = df,family = "gaussian"),inla.options) ),
      error=warning
    )
    
    
    
    if(is.character(r)){
      feil = "\n Convergence can sometimes be improved by changing the step size."
      stop(paste0(r,feil))
    }
    if(print.progress){
      cat("Completed INLA optimization in ",format(r$cpu.used[4],digits=3)," seconds..\n",sep="")
    }
  }
  
  object = list(formula=formula2,inlafit=r,.args=list(data=data,
                                                      forcing=forcing,
                                                      call = inla.ews.call,
                                                      inladata=df,
                                                      model=model,
                                                      timesteps=timesteps,
                                                      time_normalized=time_normalized,
                                                      inputformula = formula,
                                                      #intercept=intercept,
                                                      compute.mu=compute.mu,
                                                      inla.options=inla.options,
                                                      memory.true=memory.true))
  
  
  
  class(object) = "inla.ews"
  if(print.progress){
    cat("Collecting results and perform transformation to user scaling..\n",sep="")
  }
  time.start.gather = Sys.time()
  if(tolower(model) %in% c("ar1", "ar(1)", "1")){
    object = resultgather(object, nsims=nsims,print.progress=print.progress)
  }else if(tolower(model) %in% c("ar2", "ar(2)", "2")){
    object = resultgather_ar2(object, nsims=nsims,print.progress=print.progress)
  }
  
  time.gather = difftime(Sys.time(), time.start.gather,units="secs")[[1]]
  
  
  if(print.progress){
    cat("Results collected in ", format(time.gather,digits=3)," seconds..\n",sep="")
  }
  
  # if(compute.mu >0){
  #   #intercept = object$results$fixed
  #   if(compute.mu==2) object=forcingmaker(object,quick=FALSE,intercept=intercept,
  #                                         print.progress=print.progress)
  #   if(compute.mu==1) object=forcingmaker(object,quick=TRUE,intercept=intercept,
  #                                         print.progress=print.progress)
  # }
  
  time.total = difftime(Sys.time(), time.start,units="secs")[[1]]
  
  object$cpu.used = list(inla=object$inlafit$cpu.used[4],gather=time.gather,total=time.total)
  if(print.progress){
    cat("inla.ews concluded in ",time.total," seconds!\n",sep="")
  }
  options(digits=old.digits)
  return(object)
}
