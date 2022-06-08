#' Summary-function for \code{inla.ews} object
#'
#' Displays a clear summary for the information stored in the \code{inla.ews} S3 object 
#' returned by the \code{\link{inla.ews}} function.
#'
#' @param object \code{inla.ews} S3 class. Output of \code{\link{inla.ews}} function.
#' @param digits Integer indicating the number of decimal places to be used.
#' @param ... Other arguments.
#'
#' @return  Returns an object of class \code{summary.inla.ews}.
#' @examples 
#' \donttest{
#' n = 500
#' sigma = 1
#' a=0.2
#' b=0.7/n
#' time = 1:n
#' phis = a+b*time
#' s=numeric(n)
#' s[1] = rnorm(1,mean=0,sd=sigma)
#' for(i in 2:n){
#'   s[i] = rnorm(1, mean=phis[i]*s[i-1],sd=sigma)
#' }
#' 
#' object = inla.ews(data=s,model="ar1", memory.true=phis)
#' summary(object)
#' }
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords inla.ews summary
#' @export
summary.inla.ews = function(object,digits=4L,...){
  out = list()
  out = c(out,model=object$.args$model)
  maxlength=2048L
  if(sum(nchar(object$.args$call)) > maxlength){
    out=c(out, list(call=paste0( substr(deparse(object$.args$call),1L,maxlength),"...") ) )
  }else{
    out=c(out, list(call=object$.args$call))
  }
  cpu = c(object$cpu.used$inla,object$cpu.used$gather,object$cpu.used$total)
  names(cpu) = c("Running INLA","Post processing","Total")
  out = c(out, list(cpu=cpu))
  hypers = matrix(round(
    c(object$results$summary$a$mean,object$results$summary$a$sd,
      object$results$summary$a$quant0.025,object$results$summary$a$quant0.5,
      object$results$summary$a$quant0.975,
      object$results$summary$b$mean,object$results$summary$b$sd,
      object$results$summary$b$quant0.025,object$results$summary$b$quant0.5,
      object$results$summary$b$quant0.975,
      object$results$summary$sigma$mean,object$results$summary$sigma$sd,
      object$results$summary$sigma$quant0.025,object$results$summary$sigma$quant0.5,
      object$results$summary$sigma$quant0.975
    ),digits=digits), nrow=3,byrow=TRUE)
  
  colnames(hypers) = c("mean","sd","0.025quant","0.5quant","0.975quant")
  rownames = c("a","b","sigma")
  if(length(object$.args$forcing)>0){
    is.forcing=TRUE
    rownames = c(rownames,c("sigma_f","F0"))
    hypers = rbind(hypers, round(c(
      object$results$summary$sigmaf$mean,object$results$summary$sigmaf$sd,
      object$results$summary$sigmaf$quant0.025,object$results$summary$sigmaf$quant0.5,
      object$results$summary$sigmaf$quant0.975),digits=digits
      ))
    hypers = rbind(hypers, round(c(
      object$results$summary$F0$mean,object$results$summary$F0$sd,
      object$results$summary$F0$quant0.025,object$results$summary$F0$quant0.5,
      object$results$summary$F0$quant0.975),digits=digits
    ))
  }else{
    is.forcing=FALSE
  }
  rownames(hypers)=rownames
  out = c(out, hypers=list(hypers))
  out = c(out, b_positive = object$results$summary$b$prob_positive)
  
  if(tolower(object$.args$model) %in% c("ar1","ar(1)","1")){
    memorystart = c(object$results$summary$phi$mean[1],
                    object$results$summary$phi$sd[1],
                    object$results$summary$phi$q0.025[1],
                    object$results$summary$phi$q0.5[1],
                    object$results$summary$phi$q0.975[1])
    memoryend = c(object$results$summary$phi$mean[n],
               object$results$summary$phi$sd[n],
               object$results$summary$phi$q0.025[n],
               object$results$summary$phi$q0.5[n],
               object$results$summary$phi$q0.975[n])
    rownames=c("phi_1","phi_n")
  }else if(tolower(object$.args$model) %in% c("fgn","lrd")){
    memorystart = c(object$results$summary$H$mean[1],
                    object$results$summary$H$sd[1],
                    object$results$summary$H$q0.025[1],
                    object$results$summary$H$q0.5[1],
                    object$results$summary$H$q0.975[1])
    memoryend = c(object$results$summary$H$mean[n],
                  object$results$summary$H$sd[n],
                  object$results$summary$H$q0.025[n],
                  object$results$summary$H$q0.5[n],
                  object$results$summary$H$q0.975[n])
    rownames=c("H_1","H_n")
  }
  
  memorymat = rbind(memorystart,memoryend)
  memorymat = round(memorymat,digits=digits)
  colnames(memorymat) = c("mean","sd","0.025quant","0.5quant","0.975quant")
  rownames(memorymat) = rownames
  out=c(out,list(memory=memorymat))
  n=length(object$.args$data)
  out = c(out,list(n=n))
  
  
  if(!is.null(object$inlafit$dic)){
    out=c(out,list(dic=object$inlafit$dic))
  }
  if(!is.null(object$inlafit$waic)){
    out=c(out,list(waic=object$inlafit$waic))
  }
  if(!is.null(object$inlafit$mlik)){
    out=c(out,list(mlik=object$inlafit$mlik))
  }
  
  out=c(out,list(is.forcing=is.forcing))
  #out=c(out, list(family=object$inla.result$family))
  class(out) = "summary.inla.ews"
  return(out)
}

#' Print function for the \code{inla.ews} object
#'
#' Displays the return object from \code{summary.inla.ews} to an easily readable format.
#'
#' @param x \code{summary.inla.ews} S3 class. Output of \code{\link{summary.inla.ews}} function.
#' @param digits Integer indicating the number of decimal places to be used.
#' @param ... Other arguments.
#'
#' @examples 
#' \donttest{
#' n = 500
#' sigma = 1
#' a=0.2
#' b=0.7/n
#' time = 1:n
#' phis = a+b*time
#' s=numeric(n)
#' s[1] = rnorm(1,mean=0,sd=sigma)
#' for(i in 2:n){
#'   s[i] = rnorm(1, mean=phis[i]*s[i-1],sd=sigma)
#' }
#' 
#' object = inla.ews(data=s,model="ar1", memory.true=phis)
#' sumob = summary(object)
#' print(sumob)
#' }
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews},\link{summary.inla.ews}}
#' @keywords inla.ews print.summary
#' @export
print.summary.inla.ews = function(x,digits=4L,...){
  cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
  cat("Time used:\n")
  print(round(x$cpu,digits=digits))
  #cat("\n",sep="")
  cat("\nPosterior marginal distributions for all hyperparameters have been computed.\n",sep="")
  if(x$is.forcing){
    cat("\nSummary statistics for using ",x$model," model (with forcing):\n",sep="")
  }else{
    cat("\nSummary statistics using ",x$model," model:\n",sep="")
  }
  print(x$hypers)
  
  cat("\nSummary for first and last memory variable:\n",sep="")
  print(x$memory)
  
  
  cat("\nprobability of positive trend is ",x$b_positive,"\n",sep="")
  if(!is.null(x$dic)){
    cat(paste0("Deviance Information Criterion (DIC) ...: ",
               format(x$dic$dic, digits=digits, nsmall=2), "\n",
               "Effective number of parameters .........: ",
               format(x$dic$p.eff, digits=digits, nsmall=2), "\n\n"))
  }
  
  if (!is.null(x$waic)){
    cat(paste0("Watanabe-Akaike information criterion (WAIC) ...: ",
               format(x$waic$waic, digits=digits, nsmall=2), "\n",
               "Effective number of parameters .................: ",
               format(x$waic$p.eff, digits=digits, nsmall=2), "\n\n"))
  }
  
  if(!is.null(x$mlik)){
    cat(paste("\nMarginal log-Likelihood: ", format(x$mlik[2], digits=digits, nsmall=2),"\n",sep=""))
  }
  
}