# Numerical optimization

#converts from theta values
interpret.theta = function(theta, timee, nn) {
  
  tau = exp(15) 
  r = diff(range(timee))
  b_Y = -1/r+2/r*1/(1+exp(-theta[1]))
  low_Y = min(b_Y*timee[1],b_Y*timee[nn])
  high_Y = max(b_Y*timee[1],b_Y*timee[nn])
  a_Y = -low_Y + (1-high_Y+low_Y)/(1+exp(-theta[2]))
  lambdas_Y = -log(a_Y+b_Y*timee)
  cc_Y = 1/(nn-1)
  phis = c(exp(-lambdas_Y*c(1,diff(timee)/cc_Y)) ) 
  
  kappa_eps_V = exp(theta[3])
  kappa2_eps_V=kappa_eps_V*2*lambdas_Y
  
  
  
  
  C=1
  b_V = -1/r+2/r*1/(1+exp(-theta[4]))
  low_V = min(b_V*timee[1],b_V*timee[nn])
  high_V = max(b_V*timee[1],b_V*timee[nn])
  a_V = -low_V + (1-high_V+low_V)/(1+exp(-theta[5]))
  lambdas = -log(a_V+b_V*timee)
  rhos=a_V+b_V*timee
  cc = 1/(nn-1)
  rhos = c(exp(-lambdas*c(1,diff(timee)/cc)) ) 
  
  
  kappa_f = exp(theta[6])
  kappa2fs = kappa_f*2*lambdas_Y
  #F0 = theta[7]
  F0 = 0
  
  # print(theta)
  return(list(phis = phis,a_Y=a_Y,b_Y=b_Y,
              tau=tau,
              C=C,
              lambdas_Y=lambdas_Y,
              rhos=rhos,a_V=a_V,b_V=b_V,
              kappa_eps_V = kappa_eps_V, 
              kappa2_eps_V=kappa2_eps_V,
              kappa_f=kappa_f, kappa2fs=kappa2fs,F0=F0))
}

## Computes mean vector
mu = function(theta, nn, timee, fforcing) {
  
  
  #require("INLA.ews",quietly=TRUE)
  hyperparam = interpret.theta(theta, timee, nn)
  
  sf = 1/sqrt(hyperparam$kappa_f)
  kappa_f = hyperparam$kappa_f
  kappa2fs = hyperparam$kappa2fs
  #F0 = fforcing[1]
  F0= hyperparam$F0
  
  zz = (fforcing+F0)/sqrt(kappa_f)
  
  innerstruct = timee
  llambdas = hyperparam$lambdas_Y
  
  
  muvek = numeric(nn)
  
  
  for(k in 1:nn){
    for(s in 1:k){
      muvek[k] = muvek[k] + zz[s]*exp(-llambdas[k]*(timee[k]-timee[s]))
    }
  }
  
  
  muvek = muvek*(1/(sqrt(2*hyperparam$lambdas_Y)))
  
  return(c(muvek))
  
  
}

## Create precision matrix

Q = function(theta, n, timee){
  
  
  nn=n*2
  
  params = interpret.theta(theta, timee, n)
  phis = params$phis
  kappa2s = params$kappa2_eps_V
  tau=params$tau
  C=params$C
  rhos=params$rhos
  
  
  
  
  ii = c(1,(nn/2),((nn/2)+1),nn,
         2:(nn/2-1), (nn/2+2):(nn-1),
         1:((nn/2)-1),((nn/2)+1):(nn-1),
         1:((nn/2)), 
         1:((nn/2)-1))
  
  jj = c(1,(nn/2),((nn/2)+1),nn,
         2:(nn/2-1), (nn/2+2):(nn-1),
         2:(nn/2),((nn/2)+2):nn, 
         ((nn/2)+1):(nn),  
         ((nn/2)+2):(nn))
  
  
  
  
  xx = c(tau+tau*phis[2]^2,tau,C^2*tau+kappa2s[2]*rhos[2]^2+kappa2s[1],kappa2s[nn/2]+tau*C^2,
         tau+tau*phis[3:(nn/2)]^2, C^2*tau+kappa2s[2:(nn/2-1)]+kappa2s[3:(nn/2)]*rhos[3:(nn/2)]^2,
         -tau*phis[2:(nn/2)],-kappa2s[2:(nn/2)]*rhos[2:(nn/2)],
         rep(-C*tau,(nn/2)),
         C*tau*phis[2:(nn/2)])
  
  # cat("range x: ", range(xx),"\n")
  Q = Matrix::sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
  return (Q)
}

# Likelihood function
likelihood <- function(theta, y, z, time=NULL, log=TRUE){
  n=length(y)
  if(is.null(time)){
    time = seq(from=0,to=1,length.out=n)
  }
  params = interpret.theta(theta, time, n)
  
  meanvek = mu(theta, n,time, z)
  x = y-meanvek
  Qmat = Q(theta, n,time)[1:n,1:n]
  # library(Matrix)
  cf <- Cholesky(Qmat, LDL = FALSE, perm = TRUE)   # CHOLMOD
  L  <- as(cf, "sparseMatrix")
  logdetQ <- 2 * sum(log(Diagonal(x = diag(L))@x))
  log.norm.const <- -0.5 * nrow(Qmat) * log(2*pi) + 0.5 * logdetQ
  #
  #   R = sparseMatrix(i=ii,j=jj,x=xx,symmetric=TRUE)
  # normconst = -n/2*log(2*pi) + 0.5*log(1-phi^2) -n*log(innosx) 
  
  exponent = as.numeric(-0.5*(t(x)%*%Qmat%*%(x)))
  
  if(log){
    return(log.norm.const + exponent)
  }else{
    return(exp(log.norm.const+exponent))
  }
}

## Likelihood optimizer
likelicost <- function(theta,y,z, time=NULL){
  
  cost = -likelihood(theta,y,z,time)
  return(cost)
}


