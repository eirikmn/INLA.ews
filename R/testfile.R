if(FALSE){
n = 200
time=1:n
a = 0.6
b = 0.35/n
Hs = a+b*time
F0 = -3
sigmaf=0.1
sigma = 1.2
#Hs = rep(0.75,n)
library(INLA.ews)
sigmat = sigmaHmaker(sigma,Hs)
sigmachol = chol(sigmat)
noise = sigmachol%*%rnorm(n)
forcing = arima.sim(model=list(ar=c(0.9)),n)+1:n/n*10
z = sigmaf*(forcing+F0)

struct = (1:n-0.5)^(Hs-3/2)
muvek=numeric(n)
for(i in 1:n){
  muvek[i] = rev(struct[1:i])%*%z[1:i]
}

y=muvek+noise + 200
plot(y,type="l",col="grey",lwd=1.1)
#'
 object = INLA.ews::inla.ews(y,forcing,model="fgn",memory.true=Hs,compute.mu=2,
                             print.progress=TRUE)
 summary(object)
 plot(object)
 
 #############
 greenland = read.csv("/Users/emy016/Dropbox/Postdoc/timedep_LRD/working_files/Data_greenland_runoff.csv")
 
 y = greenland$Runoff_data
 forcing = greenland$Forcing
 object = INLA.ews::inla.ews(y,forcing,model="fgn",compute.mu=2,
                             print.progress=TRUE)
 summary(object)
 plot(object)
 
 
}

#######
######
if(FALSE){
  n=400
  Hs=seq(0.6,0.9,length.out=n)
  F0 = 1;sigmaf=0.2
  
  forcing=arima.sim(model=list(ar=0.9),sd=sqrt(1-0.9^2),n=n)+1:n/n/2
  
  noise = fgn_timedep_sim(n,Hs=Hs)
  muvek = mu.computer(forcing,sigmaf,F0,memory=Hs,model="fgn")
  
  T0 = 200
  y = noise + muvek+T0
  plot(y,type="l",lwd=2); lines(muvek+T0,col="red",lwd=2)
  
  objectfgn = inla.ews(y,forcing,model="fgn",compute.mu=2,memory.true = Hs,
                       print.progress=TRUE)
  summary(objectfgn)
  plot(objectfgn)
}
if(FALSE){
  n=400
  phis=seq(0.9,0.2,length.out=n)
  F0 = 1;sigmaf=0.2
  
  forcing=arima.sim(model=list(ar=0.9),sd=sqrt(1-0.9^2),n=n)+1:n/n/2
  
  noise = ar1_timedep_sim(n,phis=phis)
  muvek = mu.computer(forcing,sigmaf,F0,memory=phis,model="ar1")
  
  T0 = 200
  y = noise + muvek+T0
  plot(y,type="l",lwd=2); lines(muvek+T0,col="red",lwd=2)
  
  objectar = inla.ews(y,forcing,model="ar1",compute.mu=2,memory.true = phis,
                       print.progress=TRUE)
  summary(objectar)
  plot(objectar)
}
if(FALSE){
  cmuvek = numeric(n)
  compute_mu_fgn(cmuvek, n, zz,  Hs)
  muvek=numeric(n)
  #struct=numeric(nn)
  struct = (1:n-0.5)^(Hs-3/2)
  for(i in 1:n){
    for(j in 1:i){
      muvek[i] = muvek[i] + zz[j]*(i-j+0.5)^(Hs[i]-3/2) #struct[i-j+1]
    }
  }
  
}

