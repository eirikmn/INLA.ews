#' Print function for a \code{summary.inla.ews} object
#'
#' Displays the return object from \code{summary.inla.ews} to an easily 
#' readable format.
#'
#' @param x \code{inla.ews} S3 class. Output of \code{\link{inla.ews}} function.
#' @param digits Integer indicating the number of decimal places to print.
#' @param ... Other arguments.
#'
#' @examples 
#' \donttest{
#' n = 200
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
#' print(object)
#' }
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords inla.ews print
#' @export
print.inla.ews = function(x,digits=4L,...){
  cat("Call:\n")
  cat(deparse(x$.args$call),"\n\n",sep="")
  cat("Time used:\n",sep="")
  cpu = round(c(x$cpu.used$inla,x$cpu.used$gather,x$cpu.used$total),digits=digits)
  names(cpu) = c("Running INLA","Post processing","Total")
  print(cpu)
  
  hypers = matrix(round(
    c(x$results$summary$a$mean,x$results$summary$a$sd,
      x$results$summary$a$quant0.025,x$results$summary$a$quant0.5,
      x$results$summary$a$quant0.975,
      x$results$summary$b$mean,x$results$summary$b$sd,
      x$results$summary$b$quant0.025,x$results$summary$b$quant0.5,
      x$results$summary$b$quant0.975,
      x$results$summary$sigma$mean,x$results$summary$sigma$sd,
      x$results$summary$sigma$quant0.025,x$results$summary$sigma$quant0.5,
      x$results$summary$sigma$quant0.975
    ),digits=digits), nrow=3,byrow=TRUE)
  
  colnames(hypers) = c("mean","sd","0.025quant","0.5quant","0.975quant")
  rownames = c("a","b","sigma")
  if(length(x$.args$forcing)>0){
    is.forcing=TRUE
    rownames = c(rownames,c("sigma_f","F0"))
    hypers = rbind(hypers, round(c(
      x$results$summary$sigmaf$mean,x$results$summary$sigmaf$sd,
      x$results$summary$sigmaf$quant0.025,x$results$summary$sigmaf$quant0.5,
      x$results$summary$sigmaf$quant0.975),digits=digits
    ))
    hypers = rbind(hypers, round(c(
      x$results$summary$F0$mean,x$results$summary$F0$sd,
      x$results$summary$F0$quant0.025,x$results$summary$F0$quant0.5,
      x$results$summary$F0$quant0.975),digits=digits
    ))
  }else{
    is.forcing=FALSE
  }
  rownames(hypers)=rownames
  
  if(length(x$.args$forcing)>0){
    cat("\nSummary statistics for using ",x$model," model (with forcing):\n",sep="")
  }else{
    cat("\nSummary statistics using ",x$model," model:\n",sep="")
  }
  print(hypers)
  
  cat("\nProbability of positive trend is ",x$results$summary$b$prob_positive,"\n",sep="")
  
  if(!is.null(x$inlafit$mlik)){
    cat(paste("\nMarginal log-Likelihood: ", format(x$inlafit$mlik[2], digits=digits, nsmall=2),"\n",sep=""))
  }
}