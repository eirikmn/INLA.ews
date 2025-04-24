

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
  n=length(phis)
  time_norm=seq(0,1,length.out=n)
  #lambdas = -log(a+b*time_norm)
  lambdas = -log(phis)
  
  cc=1/(n-1)
  #phis = exp(-lambdas*c(1,diff(time_norm)/cc))
  sims=numeric(n)
  sims[1] = rnorm(1,sd=sigma/sqrt(2*lambdas[1]))
  for(i in 2:n){
    sims[i] = phis[i]*sims[i-1] + rnorm(1,sd=sigma/(sqrt(2*lambdas[i])))
  }
  return(as.numeric(sims))
}


#' Default variables in events
#'
#' Sets the default variables in the list \code{events} used to specify the climatic
#' periods and separating events used in the linear predictor. The list contains the following arguments:
#' \describe{
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
    control.inla=list(h=0.005),
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

#' Transform from internal scaling
#'
#' Transforms parameters from internal scaling to user friendly parameters. 
#' Also includes other information transformed from parameters
#' @param theta \code{numeric} of parameters in internal scaling.
#' @param time_norm \code{numeric} that includes normalized time steps.
#' @param do.forcing \code{Boolean} that indicates whether or not forcing is used.
#'
#' @return Returns a list with parameters in user scaling (\code{pars}) as well as 
#' other useful information obtained from the hyperparameters.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords inla.ews default
#' @export
from.theta = function(theta, time_norm, do.forcing=FALSE){
  n=length(time_norm)
  if(do.forcing){
    pars = numeric(5)
    
    kappa_f = exp(theta[4])
    kappafs = kappa_f*2*lambdas
    F0 = theta[5]
    pars[4]=1/sqrt(kappa_f)
    pars[5]=F0
  }else{
    pars = numeric(3)
  }
  
  kappa_eps = exp(theta[1])
  r = diff(range(time_norm))
  b = -1/r+2/r*1/(1+exp(-theta[2]))
  
  low = min(b*time_norm[1],b*time_norm[n])
  high = max(b*time_norm[1],b*time_norm[n])
  
  a = -low + (1-high+low)/(1+exp(-theta[3]))
  
  pars[1]=1/sqrt(kappa_eps)
  pars[2]=b
  pars[3]=a
  lambdas = -log(a+b*time_norm)
  kappas = kappa_eps*2*lambdas
  
  cc = 1/(n-1)
  phis = c(exp(-lambdas*c(1,diff(time_norm)/cc)) ) #rescale
  if(do.forcing){
    return(list(pars=pars, phis = phis, kappa_eps = kappa_eps, a=a,b=b,#kappa_f=kappa_f,F0=F0
                lambdas=lambdas, kappas=kappas, sigmas=1/kappas,
                kappa_f=kappa_f, kappafs=kappafs, sigmafs=1/sqrt(kappafs), 
                F0=F0))
  }else{
    return(list(pars=pars, phis = phis, kappa_eps = kappa_eps, a=a,b=b,#kappa_f=kappa_f,F0=F0
                lambdas=lambdas, kappas=kappas, sigmas=1/kappas))
  }
  
}

#' Transform to internal scaling
#'
#' Transforms parameters to internal scaling from user friendly parameters. 
#' Also includes other information transformed from parameters
#' @param pars \code{numeric} of parameters in user scaling.
#' @param time_norm \code{numeric} that includes normalized time steps.
#' @param do.forcing \code{Boolean} that indicates whether or not forcing is used.
#'
#' @return Returns a list with parameters in internal scaling (\code{theta}) as well as 
#' other useful information obtained from the hyperparameters.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords inla.ews default
#' @export
to.theta = function(pars, time_norm, do.forcing=FALSE){
  n=length(time_norm)
  if(do.forcing){
    theta = numeric(5)
    
    theta[4] = log(1/pars[4]^2)
    theta[5] = pars[5]
  }else{
    theta = numeric(3)
  }
  kappa_eps = 1/pars[1]^2
  a=pars[3]
  b=pars[2]
  
  theta[1] = log(1/pars[1]^2)
  
  r = diff(range(time_norm))
  theta[2] = log((1+pars[2])/(1-pars[2]))
  low = min(pars[2]*time_norm[1],pars[2]*time_norm[n])
  high = max(pars[2]*time_norm[1],pars[2]*time_norm[n])
  theta[3] = log((pars[3]+low)/(1-high-pars[3]))
  
  lambdas = -log(pars[3]+pars[2]*time_norm)
  kappas = kappa_eps*2*lambdas
  
  cc = 1/(n-1)
  phis = c(exp(-lambdas*c(1,diff(time_norm)/cc)) ) #rescale
  
  kappa_f = exp(theta[4])
  kappafs = kappa_f*2*lambdas
  F0 = theta[5]
  
  if(do.forcing){
    return(list(theta=theta, phis = phis, kappa_eps = kappa_eps, a=a,b=b,#kappa_f=kappa_f,F0=F0
                lambdas=lambdas, kappas=kappas, sigmas=1/kappas,
                kappa_f=kappa_f, kappafs=kappafs, sigmafs=1/sqrt(kappafs), 
                F0=F0))
  }else{
    return(list(theta=theta, phis = phis, kappa_eps = kappa_eps, a=a,b=b,#kappa_f=kappa_f,F0=F0
                lambdas=lambdas, kappas=kappas, sigmas=1/kappas))
  }
}

#' Compute forcing response
#'
#' Computes forcing response given parameters in either internal or user scaling using Monte Carlo sampling.
#' @param pars \code{numeric} of parameters.
#' @param forcing \code{numeric} that includes known forcing component.
#' @param time_norm \code{numeric} that includes normalized time steps.
#' @param as.theta \code{Boolean} that indicates whether or not \code{pars} is given 
#' in internal or user scaling.
#'
#' @return Returns a \code{numeric} with the estimated forcing response.
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords inla.ews default
#' @export
mucomputer = function(pars,forcing,time_norm,as.theta=FALSE){
  if(as.theta){
    allinfo = from.theta(pars,time_norm,do.forcing=TRUE)
    parstheta = pars
    parsttemp = allinfo$pars
  }else{
    parstemp=pars
    allinfo = to.theta(pars,time_norm,do.forcing=TRUE)
    parstheta = allinfo$theta
  }
  n=length(time_norm)
  a=allinfo$a; b=allinfo$b
  sigmafs=allinfo$sigmafs
  F0=allinfo$F0
  zz = sigmafs*(F0+forcing)
  lambdas = -log(a+b*time_norm)
  muvek = numeric(n)
  #
  c_mu_ar1(muvek,z=zz, n=n,lambda=lambdas,time=time_norm)
  #
  
  # if(diff(range(diff(time_norm)))<10^(-12)){
  #   struktur = exp(-lambdas*time_norm)
  #   for(i in 1:n){
  #     muvek[i] = rev(struktur[1:i])%*%zz[1:i]
  #   }
  # }else{
  #   c_mu_ar1(muvek,z=zz, n=n,lambda=lambdas,time=time_norm)
  #   
  #   #for(k in 1:n){
  #   #  for(s in 1:k){
  #   #    muvek[k] = muvek[k] + zz[s]*exp(-lambdas[k]*(time_norm[k]-time_norm[s]))
  #   #  }
  #   
  # }
  
  return(muvek)
}