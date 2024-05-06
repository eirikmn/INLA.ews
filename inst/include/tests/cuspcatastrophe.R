if(FALSE){
  library(readr)
  catsimpath = "/cloud/project/sim_catastrophe.csv"
  #data = read_csv(catsimpath)
  simdata = scan(catsimpath)
  
  xi=1
  
  Vfunc = function(x, xi=1, mu=1){
    x^4/4 -xi*x^2/2  -mu*x
  }
  Vdiff = function(x, xi=1, mu=1){
    -x^3 + xi*x + mu
  }
  
  
  mu1 = -2/3*sqrt(xi^3/3)
  mu2 = 2/3*sqrt(xi^3/3)
  
  
  #tipping points
  scale = 2
  
  mus = c(0.15, mu2, 1)
  xx = seq(-1.5,2,length.out=1000)
  
  x1s = sort(Re(polyroot(c(mus[1],xi,0,-1))))
  x2s = sort(Re(polyroot(c(mus[2],xi,0,-1))))
  x3s = max(Re(polyroot(c(mus[3],xi,0,-1))))
  
  y1s = scale*Vfunc(x1s, mu=mus[1])
  y2s = scale*Vfunc(x2s, mu=mus[2])
  y3s = scale*Vfunc(x3s, mu=mus[3])
  
  ggd = data.frame(x=xx,
                   v1 = scale*Vfunc(xx,mu=mus[1]),
                   v2 = scale*Vfunc(xx,mu=mus[2]),
                   v3 = scale*Vfunc(xx,mu=mus[3]))
  
  ggd2_1 = data.frame(x1=x1s[c(1,3)],x2=x1s[2],
                      y1=y1s[c(1,3)],y2=y1s[2])
  ggd2_2 = data.frame(x1=x2s[c(1,3)],x2=x2s[2],
                      y1=y2s[c(1,3)],y2=y2s[2])
  ggd2_3 = data.frame(x1=x3s[c(1,3)],x2=x3s[2],
                      y1=y3s[c(1,3)],y2=y3s[2])
  #ggd2 = data.frame(x1=x1s,x2=x2s,x3=x3s, y1=y1s,y2=y2s,y3=y3s)
  
  library(ggplot2)
  library(ggpubr)
  #plot a
  gg0 = ggplot(ggd,aes(x=xx)) + theme_bw() + xlab("State variable x") + 
    ylab("2V(x)") + theme(text=element_text(size=16), plot.title = element_text(size=22)) 
  gg1 = gg0 + geom_line(aes(y=v1)) + ggtitle("(a) Before transition") +
    geom_point(data=ggd2_1,aes(x=x1, y=y1), col="black",size=2.5) +
    geom_point(data=ggd2_1, aes(x=x2, y=y2),col="red",size=2.5)
  gg2 = gg0 + geom_line(aes(y=v2)) + ggtitle("(b) At transition") +
    geom_point(data=ggd2_2,aes(x=x1, y=y1), col="black",size=2.5) +
    geom_point(data=ggd2_2, aes(x=x2, y=y2),col="red",size=2.5)
  gg3 = gg0 + geom_line(aes(y=v3)) + ggtitle("(c) After transition") +
    geom_point(data=ggd2_3,aes(x=x1, y=y1), col="black",size=2.5) #+
  #geom_point(data=ggd2_3, aes(x=x2, y=y2),col="red",size=2.5)
  
  ggtip = ggarrange(gg1,gg2,gg3, nrow=1)
  
  ggsave(paste0("ggtipping-15x5.eps"),plot=ggtip, device=cairo_ps, width=15,
         height=5, limitsize=FALSE)
  
  sim = simdata[-1]
  muu = seq(-1,1,length.out=length(sim))
  lower = upper = middle = rep(NA,length(muu))
  for(i in 1:length(muu)){
    roots = sort(Re(polyroot(c(muu[i],xi,0,-1))))
    if(muu[i] < mu1){
      lower[i] = roots[1]
    }else{
      if(muu[i]<mu2){
        lower[i] = roots[1]
        middle[i] = roots[2]
        upper[i] = roots[3]
      }else{
        upper[i] = roots[3]
      }
    }
  }
  
  ggdb = data.frame(muu=muu,lower=lower,middle=middle,upper=upper,
                    sim=sim)
  ggb = ggplot(ggdb, aes(x=muu)) + theme_bw() + 
    ylab("State variable x") + 
    xlab(expression(paste("Control parameter ",mu))) + 
    ggtitle("(d) Bifurcation diagram")+
    theme(text=element_text(size=16), plot.title = element_text(size=22)) +
    geom_line(aes(y=lower))+
    geom_line(aes(y=upper))+
    geom_line(aes(y=middle), linetype="dashed") +
    geom_line(aes(y=sim),col="red")
  ggb
  ggsave(paste0("bifurcation.eps"),plot=ggb, device=cairo_ps, width=10,
         height=7, limitsize=FALSE)
  ggb2 = ggb + theme(plot.margin = unit(c(2,10,1,7),"lines")) #(top,right,bottom,left)
  ggall = ggarrange(ggtip,ggb2,nrow=2, heights=c(0.4,0.6))
  ggall
  ggsave(paste0("ggtipping_all.eps"),plot=ggall, device=cairo_ps, width=15,
         height=12, limitsize=FALSE)
  
  
}