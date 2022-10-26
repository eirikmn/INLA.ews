if(FALSE){
  


data = NGRIP_HR_f2$V1
proxy=rev(data[-1])
time = time_test2
do.detrend = TRUE
model="fgn"
ewsres = c()
probpos = numeric(length(Clear_GI_onsets))
Hlist = c()
ameans = numeric(length(Clear_GI_onsets))
bmeans = numeric(length(Clear_GI_onsets))
cmeans = numeric(length(Clear_GI_onsets))

plot(time,proxy,xlab="Time",ylab="d18O",type="l",col="red",xlim=rev(range(time)))
for(i in 1:length(Clear_GI_onsets)){
  cat("running model ",model," with detrending=",do.detrend,": iteration #",i,"\n",sep="")
  startage = events_rasmussen$age[Clear_GS_onsets[i]]
  endage = events_rasmussen$age[Clear_GI_onsets[i]]
  indices = rev(which(time<startage & time>endage))
  windowy = proxy[indices]
  windowx = time[indices]
  lines(windowx,windowy,col="blue")
  
  n = length(windowx)
  timenorm = seq(0,1,length.out=n)
  df = data.frame(y = windowy,idy=1:n)
  rgen_model = INLA::inla.rgeneric.define(rgeneric.ews.fgn2,n=n,time=1:n)
  formula = y ~ -1+ f(idy, model=rgen_model)
  r = inla(formula,data=df,family="gaussian", control.family=list(initial=12,fixed=TRUE),
           control.compute=list(config=TRUE),verbose=TRUE,
           control.inla=list(h=0.001))
  
  probpos[i] = ewsfit$results$summary$b$prob_positive
  ewsres = c(ewsres, list(ewsfit))
}



df = data.frame(y = noise,idy=1:n)

rgen_model = INLA::inla.rgeneric.define(rgeneric.ews.fgn2,n=n,time=time)

formula = y ~ -1+ f(idy, model=rgen_model)

r = inla(formula,data=df,family="gaussian", control.family=list(initial=12,fixed=TRUE),
control.compute=list(config=TRUE))

summary(r)

sig = inla.emarginal(function(x)1/sqrt(exp(x)),r$marginals.hyperpar$`Theta1 for idy`)
cmean = inla.emarginal(function(x)exp(x),r$marginals.hyperpar$`Theta4 for idy`)

nsims = 5000
hyperpars = inla.hyperpar.sample(nsims,r)
csamps = exp(hyperpars[,4])
sbar = 1-0.5
bsamps = -sbar+2*sbar/(1+exp(-hyperpars[,3]))
asamps = numeric(nsims)
memat = matrix(NA,n,nsims)
for(i in 1:nsims){
   alower = 0.5 -min(bsamps[i]*timenorm)
   aupper = 1 -max(bsamps[i]*timenorm)
   asamps[i] = alower+(aupper-alower)/(1+exp(-hyperpars[i,2]))
   memat[,i] = asamps[i] + bsamps[i]*timenorm^csamps[i]
}
mean(asamps)
mean(bsamps)
mean(csamps)
memean = numeric(n)
meupper = numeric(n)
melower=numeric(n)
for(i in 1:n){
   dens = density(memat[i,]); dens=data.frame(x=dens$x,y=dens$y)
   zme = inla.zmarginal(dens,silent=TRUE)
   memean[i] = zme$mean
   melower[i] = zme$quant0.025
   meupper[i] = zme$quant0.975
}
plot(timenorm,meupper,type="l",col="red",ylim=c(0.5,1))
lines(timenorm,melower,type="l",col="red")
lines(timenorm,memean,type="l",col="blue")
lines(timenorm,phis,col="black")



data = NGRIP_HR2$V1
proxy=rev(data[-1])
time = time_test2
do.detrend = FALSE
model="ar1"
ewsres = c()
probpos = numeric(length(Clear_GI_onsets))
{
  plot(time,proxy,xlab="Time",ylab="d18O",type="l",col="red",xlim=rev(range(time)))
  for(i in 1:length(Clear_GI_onsets)){
    startage = events_rasmussen$age[Clear_GS_onsets[i]]
    endage = events_rasmussen$age[Clear_GI_onsets[i]]
    indices = which(time<startage & time>endage)
    windowy = proxy[indices]
    windowx = time[indices]
    lines(windowx,windowy,col="blue")
    if(do.detrend){
      lmfit = lm(y~1+x,data=data.frame(y=windowy,x=windowx))
      ewsfit <- inla.ews(lmfit$residuals,model=model)
    }else{
      ewsfit <- inla.ews(windowy,model=model)
    }
    probpos[i] = ewsfit$results$summary$b$prob_positive
    ewsres = c(ewsres, list(ewsfit))
  }}
{plot(probpos,type="b");abline(h=0.5,lty=2);abline(h=0.95,lty=3)}
summary(ewsres[[1]])






}