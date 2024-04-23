

set.seed(123)
n=1000
a = 0.3
b = 0.2

time = sort(1:n + rnorm(n,sd=0.1))
#time = 1:n
time_norm = (time-min(time))/(max(time)-min(time))
time_normalized=time_norm
lambdas = -log(a+b*time_norm)
cc=1/(n-1)
phis = exp(-lambdas*c(1,diff(time_norm)/cc))

plot(phis)
lines(exp(-lambdas),col="red",lwd=3)

forcing = arima.sim(list(ar=c(0.95)), n=n)

#phis=rep(phi,n)
F0=0
sigma=5
sigma_f = 0.1


sims=numeric(n)
sims[1] = rnorm(1,sd=sigma/sqrt(2*lambdas[1]))
for(i in 2:n){
  sims[i] = phis[i]*sims[i-1] + rnorm(1,sd=sigma/sqrt(lambdas[i]))
}
#sims = sigma * arima.sim(list(ar=c(phi)), n=n, sd=sqrt(1-phi^2))
muvek = mucomputer(pars=c(sigma,b,a,sigma_f,F0),forcing,time_norm,as.theta=FALSE)

intercept=0
y = sims + muvek + intercept
plot(y)
lines(muvek,col="red",lwd=5)

#zzz = (forcing+F0)*sigma_f
theta=numeric(5)
theta[1] = log(1/sigma^2)
theta[2] = log( (1+b)/(1-b) )
alower = -min(b,0); aupper = 1-max(b,0)
theta[3] = log( (a-alower)/(aupper-a) )
theta[4] = log(1/sigma_f^2)
lambdas = -log(a+b*time_norm)
kappa2fs = 2*lambdas/sigma_f^2

zzz = (forcing+F0)/sqrt(kappa2fs)
###### 
## EMN: HERE
rrc = inla.ews(data=y, forcing=forcing, timesteps = time, print.progress=TRUE,
               do.cgeneric=TRUE)
summary(rrc)


#cgeneric here


######
cscript.path=c()
for(i in 1:length(.libPaths())){
  if("INLA.ews" %in% list.files(.libPaths()[i])){
    cscript.path=paste0(.libPaths()[i],"/INLA.ews/libs/INLA.ews.so")
  }
}
if(length(cscript.path)==0){
  stop("Could not find package directory, please make sure that INLA.climate is installed within one of the libraries displayed by '.libPaths()'.")
}


cmodel <- INLA::inla.cgeneric.define(model = "inla_cgeneric_timedep_forcing",
                                     #shlib = "src/cgeneric.so", 
                                     shlib = cscript.path, 
                                     n = n, debug=FALSE,time=as.numeric(time), forcing=forcing
)
qq= inla.cgeneric.q(cmodel)

rc <- inla(
  y ~ -1 + f(idx, model = cmodel),
  data = data.frame(y=y, idx = 1:n),
  control.family = list(hyper = list(prec = list(initial = 12, fixed = TRUE))))

summary(rc)

2*lambdas[i]/sigma_f


rrr = inla.ews(y,forcing,timesteps = time,print.progress=TRUE,nsims=10000)
summary(rrr)
#### 




#### #### ###
#### #### ###
#### #### ###
rrm = inla.rgeneric.define(rgeneric.ar1.forcing.2, 
                           n=n,debug=FALSE,
                           time=time_normalized,forcing=forcing
)
qr = inla.rgeneric.q(rrm,cmd = c("graph", "Q", "mu", "initial", "log.norm.const", "log.prior", "quit"),theta=rep(0.1,5))
rr <- inla(
  y ~ -1 + f(idx, model = rrm),
  data = data.frame(y=y, idx = 1:n),verbose=TRUE,
  control.family = list(hyper = list(prec = list(initial = 12, fixed = TRUE))))
summary(rr)

#### #### ###
#### #### ###
#### #### ###
#
#rrr=inla.ews(y,forcing=forcing,timesteps=time,formula=y~-1, nsims=10)
#summary(rrr$inlafit)
testr = rc

testr = rr

nsims=1000
hyps = inla.hyperpar.sample(nsims,testr)

bsims = -1 + 2/(1+exp(-hyps[,2]))
sfsims = 1/sqrt(exp(hyps[,4]))
F0sims = hyps[,5]
asims= numeric(0)

phisims = matrix(NA,nrow=n,ncol=nsims)
phi0sims = matrix(NA,nrow=n,ncol=nsims)
lambdasims = matrix(NA,nrow=n,ncol=nsims)
muveksims = matrix(NA,nrow=n,ncol=nsims)

for(i in 1:nsims){
  alower = -min(bsims[i],0)
  aupper = 1-max(bsims[i],0)
  asims[i] = alower + (aupper-alower)/(1+exp(-hyps[i,3]))
  phi0sims[,i] = asims[i]+bsims[i]*time_norm
  lambdasims[,i] = -log(phi0sims[,i])
  
  muveksims[,i] = mucomputer(pars=c(0,bsims[i],asims[i],sfsims[i], F0sims[i]),
                              forcing = forcing,
                              time_norm = time_normalized, as.theta=FALSE)
}

mumean = numeric(n)
muupper = numeric(n)
mulower = numeric(n)
for(i in 1:nsims){
  dens0 = density(muveksims[i,]); dens=data.frame(x=dens0$x,y=dens0$y)
  zm = inla.zmarginal(dens,silent=TRUE)
  mumean[i] = zm$mean
  mulower[i] = zm$quant0.025
  muupper[i] = zm$quant0.975
}
plot(y,type="l", lwd=1)
lines(muvek,type="l", lwd=5,col="red")
lines(mumean,type="l", lwd=5,col="blue")

#

resmatr = data.frame(a=mean(asims),b=mean(bsims),sf=mean(sfsims),F0=mean(sfsims))
resmatc = data.frame(a=mean(asims),b=mean(bsims),sf=mean(sfsims),F0=mean(sfsims))
print(resmatr)
