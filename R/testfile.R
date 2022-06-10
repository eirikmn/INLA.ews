if(FALSE){
n = 300
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

y=muvek+noise
plot(y,type="l",col="grey",lwd=1.1)
#'
 object = INLA.ews::inla.ews(y,forcing,model="fgn",memory.true=Hs)
 summary(object)
 plot(object)
 
 plot(object)
}