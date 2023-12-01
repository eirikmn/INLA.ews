

#' Simulate time dependent AR(1) series
#' 
#' This function produces samples from a given time dependent AR(1) process.
#' 
#' @param n The length of the simulated time series.
#' @param sigma The standard deviation of the innovations.
#' @param a Intercept in the evolution of the lag-one correlation.
#' @param b Slope in the evolution of the lag-one correlation.
#' @param phis Numeric of length \code{n}. If evolution of lag-one correlation is 
#' to be given explicitly this is done here. Overrides \code{a} and \code{b}.
#' @return Returns the simulated time series as a \code{numeric} object.
#' 
#' @examples 
#' n = 200
#' sims = ar1_timedep_sim(n,sigma=1,a=0.2,b=0.7)
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords simulation ar1 timedep 
#' @export
#' @importFrom stats rnorm
ar1_timedep_sim <- function(n,sigma=1,a=0.2,b=0.7,phis=NULL){
  if(is.null(phis)){
    phis = a+b*seq(0,1,length.out=n)
  }
  noise=numeric(n)
  noise[1] = rnorm(1,mean=0,sd=sigma)
  for(i in 2:n){
    noise[i] = rnorm(1, mean=phis[i]*noise[i-1],sd=sigma)
  }
  return(as.numeric(noise))
}


#' Default variables in events
#'
#' Sets the default variables in the list \code{events} used to specify the climatic
#' periods and separating events used in the linear predictor. The list contains the following arguments:
#' \itemize{
#'   \item{\code{num.threads} }{Integer describing the number of cores used by the INLA program. 
#'   For rgeneric models, stabiltiy is sometimes improved by setting this equal to \code{1} (default).}
#'   \item{\code{control.inla} }{List containing \code{h=0.01} and \code{restart=1}. 
#'   \code{h} denotes the step size. Changing this might help reach convergence. 
#'   \code{restart} is an integer indicating how many times INLA should be restart to 
#'   at the found optimum. This might help improve the accuracy of the optimization procedure.}
#'   \item{\code{control.family} }{List with settings to tell INLA to set the variance of the 
#'   Gaussian likelihood to be fixed and equal to \code{exp(-12)}.}
#' }
#' @return Returns a list including default values for all variables in \code{inla.options}.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords inla.ews inla options default
#'
inla.options.default <- function(){
  return(list(
    num.threads=1,
    control.mode=list(restart=TRUE),
    control.inla=list(h=0.005,restart=1),
    control.family = list(hyper = list(prec = list(initial = 12, fixed=TRUE)))
  )
  )
}

#' Import default arguments
#'
#' Fills out missing arguments in list \code{opt} with default arguments in list
#' \code{default.opt}.
#' @param opt List object with different specifications.
#' @param default.opt List of default variables corresponding to \code{opt}.
#'
#' @return Returns the \code{opt} list, but with values from \code{default.opt} inserted
#' in missing values.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords inla.ews default
set.options <- function(opt,default.opt){
  temp = default.opt
  
  if(length(opt)>0){
    for(i in 1:length(opt)){
      if(names(opt)[i] %in% names(default.opt)){
        if(!is.list(opt[[i]])){
          temp[[ names(opt)[i] ]] <- opt[[i]]
        }else{
          for(j in 1:length(opt[[i]])){
            temp[[ names(opt)[i] ]][[names(opt[[i]])[j]]] <- opt[[i]][[j]]
          }
        }
      }else{
        temp[[names(opt)[i]]] <- opt[[i]]
      }
    }
  }
  return(temp)
}
