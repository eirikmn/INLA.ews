if(FALSE){
  ## cusp catastrophe 
  pot = function(x,mu,lam) x^4/4-lam*x^2/2-mu*x
  dpot = function(x,mu,lam)x^3-lam*x-mu
  dpot = function(x,mu,lam)-x^3+lam*x+mu
  
  disc = function(mu,lam) mu^2-4/27*lam^3
  
  lam=1
  bp1 = -2/3*sqrt(lam^3/3)
  bp2 = 2/3*sqrt(lam^3/3)
  
  
  
  #plot(xx[2:length(xx)],diff(pot(xx,mu,lam)),type="l")
  xx = seq(from=-2.5,to=2,length.out=300)
  
  par(mfrow=c(1,1))
  
  mu = bp1-0.1
  plot(xx,2*pot(xx,mu,lam),type="l",ylim=c(-1.6,3),xlim=c(-2,2),xlab="State variable x",
       ylab=expression(paste("2V(x, ",lambda,", ",mu[3],")")),
       main="a) Before transition")
  polyvec = c(-mu,-lam,0,1)
  res = polyroot(polyvec)
  r = sort(Re(res))
  #r = uniroot(dpot2,c(-2.5,2))
  points(r[c(1,3)],2*pot(r[c(1,3)],mu,lam),pch=19)
  points(r[2],2*pot(r[2],mu,lam))
  
  
  mu = bp1
  plot(xx,2*pot(xx,mu,lam),type="l",ylim=c(-1.6,3),xlim=c(-2,2),xlab="State variable x",
       ylab=expression(paste("2V(x, ",lambda,", ",mu[2],")")),
       main="b) At transition")
  polyvec = c(-mu,-lam,0,1)
  res = polyroot(polyvec)
  r = sort(Re(res))
  #r = uniroot(dpot2,c(-2.5,2))
  points(r[1],2*pot(r[1],mu,lam))
  points(r[3],2*pot(r[3],mu,lam),pch=19)
  
  
  
  mu = bp1+0.1
  plot(xx,2*pot(xx,mu,lam),type="l",ylim=c(-1.6,3),xlim=c(-2,2),xlab="State variable x",
       ylab=expression(paste("2V(x, ",lambda,", ",mu[1],")")),
       main="c) After transition")
  polyvec = c(-mu,-lam,0,1)
  res = polyroot(polyvec)
  r = sort(Re(res[abs(Im(res))<0.0001]))
  points(min(r),2*pot(min(r),mu,lam),pch=19)
  #points(r,2*pot(r,mu,lam))
  
  
  
  
  
  mu = bp2-0.1
  plot(xx,2*pot(xx,mu,lam),type="l",ylim=c(-1.6,3),xlim=c(-2,2),xlab="State variable x",
       ylab=expression(paste("2V(x, ",lambda,", ",mu[3],")")),
       main="a) Before transition")
  polyvec = c(-mu,-lam,0,1)
  res = polyroot(polyvec)
  r = sort(Re(res))
  #r = uniroot(dpot2,c(-2.5,2))
  points(r[c(1,3)],2*pot(r[c(1,3)],mu,lam),pch=19)
  points(r[2],2*pot(r[2],mu,lam))
  
  
  mu = bp2
  plot(xx,2*pot(xx,mu,lam),type="l",ylim=c(-1.6,3),xlim=c(-2,2),xlab="State variable x",
       ylab=expression(paste("2V(x, ",lambda,", ",mu[2],")")),
       main="b) At transition")
  polyvec = c(-mu,-lam,0,1)
  res = polyroot(polyvec)
  r = sort(Re(res))
  #r = uniroot(dpot2,c(-2.5,2))
  points(r[1],2*pot(r[1],mu,lam))
  points(r[3],2*pot(r[3],mu,lam),pch=19)
  
  
  
  mu = bp2+0.1
  plot(xx,2*pot(xx,mu,lam),type="l",ylim=c(-1.6,3),xlim=c(-2,2),xlab="State variable x",
       ylab=expression(paste("2V(x, ",lambda,", ",mu[1],")")),
       main="c) After transition")
  polyvec = c(-mu,-lam,0,1)
  res = polyroot(polyvec)
  r = sort(Re(res[abs(Im(res))<0.0001]))
  points(min(r),2*pot(min(r),mu,lam),pch=19)
  #points(r,2*pot(r,mu,lam))
  
  
  
  
  
  ####### 
  
  
  muu = seq(from=bp1-0.1,to=bp2+0.1,length.out=1000)
  #muu = seq(from=-0.6,to=0.6,length.out=500)
  muu = seq(from=-1,to=1,length.out=500)
  nt=length(muu)
  rootmat = matrix(NA,3,nt)
  for(i in 1:nt){
    polyvec = c(muu[i],lam,0,-1)
    res = sort(polyroot(polyvec))
    isreal = abs(Im(res))<0.00001
    res[!isreal] = NA
    rootmat[,i]=res
  }
  
  plot(muu,Re(rootmat)[3,],type="l",ylim=c(-1.5,1.5))
  lines(muu,Re(rootmat)[2,],lty=2)
  lines(muu,Re(rootmat)[1,])
  
  sim_cat = read.csv("/Users/emy016/Dropbox/Postdoc/timedep_LRD/working_files/sim_catastrophe.csv")
  lines(muu,sim_cat$X.1.,col="red")
  
  ggd = data.frame(muu=muu,lower=Re(rootmat)[1,],upper=Re(rootmat)[3,],mid=Re(rootmat)[2,],
                   sim=sim_cat$X.1.)
  
  ggp = ggplot(data=ggd,aes(x=muu))+theme_bw()+
    geom_line(aes(y=lower)) + geom_line(aes(y=mid),linetype="dashed") + geom_line(aes(y=upper))+
    geom_line(aes(y=sim),color="red") + 
    xlab(expression(paste("Control parameter ",mu,""))) +
    ylab(expression(paste("State variable ",x,"")))
  print(ggp)
  
  ###
  
  mu = 0.3
  polyvec = c(-mu,-lam,0,1)
  res = polyroot(polyvec)
  r1 = sort(Re(res))
  
  xx = seq(from=-1.9,to=1.9,length.out=300)
  
  ggp1 = ggplot(data=data.frame(xx=xx,pot2 = 2*pot(xx,0.3,lam)),aes(x=xx)) + 
    theme_bw()+
    geom_line(aes(y=pot2))+
    geom_point(data=data.frame(r=r1[c(1,3)],y=2*pot(r1[c(1,3)],0.3,lam)),
               aes(x=r,y=y),size=2.5) +
    geom_point(aes(x=r1[2],y=2*pot(r1[2],0.3,lam)),size=2.5,shape=1) +
    xlab("State variable x")+ylab(expression(paste("2V(x)")))+#,lambda,", ",mu[1],")")))+
    labs(title="(a) Before transition") + 
    ylim(c(-1.8,1.8))+xlim(c(-1.5,2))
  print(ggp1)
  
  ###
  
  mu = bp2
  polyvec = c(-mu,-lam,0,1)
  res = polyroot(polyvec)
  r2 = sort(Re(res))
  ggp2 = ggplot(data=data.frame(xx=xx,pot2 = 2*pot(xx,bp2,lam)),aes(x=xx)) + 
    theme_bw()+
    geom_line(aes(y=pot2))+
    geom_point(aes(x=r2[1],y=2*pot(r2[1],bp2,lam)),size=2.5,shape=1) +
    geom_point(aes(x=r2[3],y=2*pot(r2[3],bp2,lam)),size=2.5) +
    xlab("State variable x")+ylab(expression(paste("2V(x)")))+#, ",lambda,", ",mu[1],")")))+
    labs(title="(b) At transition")+ 
    ylim(c(-1.8,1.8))+xlim(c(-1.5,2))
  print(ggp2)
  ###
  
  mu = 0.5
  polyvec = c(-mu,-lam,0,1)
  res = polyroot(polyvec)
  r3 = sort(Re(res[abs(Im(res))<0.0001]))
  
  
  ggd = data.frame(xx=xx,pot2 = 2*pot(xx,mu,lam))
  
  ggp3 = ggplot(data=data.frame(xx=xx,pot2 = 2*pot(xx,mu,lam)),aes(x=xx)) + 
    theme_bw()+
    geom_line(aes(y=pot2))+
    geom_point(aes(x=min(r3),y=2*pot(min(r3),0.5,lam)),size=3) +
    xlab("State variable x")+ylab(expression(paste("2V(x)")))+#, ",lambda,", ",mu[1],")")))+
    labs(title="(c) After transition")+ 
    ylim(c(-1.8,1.8))+xlim(c(-1.5,2))
  print(ggp3)
  
  
  library(patchwork)
  ggp1+ggp2+ggp3
  
  ####### simulate ###
  
  
  muu = seq(from=bp2-0.1,to=bp1+0.1,length.out=1000)
  #muu=rev(muu)
  muu=-muu
  nt=length(muu)
  rootmat = matrix(NA,3,nt)
  xsims = numeric(nt);xsims[1]=-1.186831
  
  sigma=0.1;lam=1
  dpots=numeric(nt)
  #dpot = function(x,mu,lam)x^3-lam*x+mu
  dpot2 = function(x,mu,lam)x^3-lam*x+mu
  loweq = numeric(nt)
  eqdiff = numeric(nt)
  
  for(i in 1:nt){
    polyvec = c(muu[i],-lam,0,1)
    res = sort(polyroot(polyvec))
    isreal = abs(Im(res))<0.00001
    res[!isreal] = NA
    rootmat[,i]=res
    #
    loweq[i] = Re(rootmat[1,i])
    if(i>1){
      dx = -dpot2(xsims[i-1],muu[i-1],lam) + sigma*rnorm(1)
      xsims[i] = xsims[i-1] + dx
      dpots[i]=dpot2(xsims[i-1],muu[i],lam)
      
      eqdiff[i] = loweq[i]-xsims[i]
      #xsims[i] = xsims[i-1]*(1+lam)+muu[i-1]-xsims[i-1]^3 + sigma*rnorm(1)
      #xsims[i] = rcusp(1,beta=1,alpha=-muu[i])
    }
    
  }
  id=2
  eps = 0.01
  eq=loweq[id]; mu = muu[id]
  
  dpot2(eq+eps,mu,1)
  (eq^3+3*eq^2*eps+3*eq*eps^2+eps^3)-1*eq-1*eps+mu
  3*eq^2*eps+3*eq*eps^2+eps^3-1*eps
  #
  
  plot(Re(rootmat)[3,],type="l",ylim=c(-1.5,1.5),xlab="Time",ylab="State")
  lines(Re(rootmat)[2,],lty=2)
  lines(Re(rootmat)[1,])
  lines(xsims,col="red")
  #
  
  plot(eqdiff)
  idd = 520:550
  plot(loweq[idd]);lines(xsims[idd],col="red")
  plot(xsims-Re(rootmat[1,]))
  
  eq=Re(rootmat)[1,1]
  
  po = -dpot2(-1,muu[1],1)
  
  -dpot2(-1.1,muu[3],1)
  plot(dpots,type="l")
  
  
  
}
