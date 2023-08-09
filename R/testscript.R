if(FALSE){
  
  library(devtools)
  check(args="--no-examples")
  
  set.seed(1)
  ### testing ###
  
  plotoptions = list(plot.hyper=FALSE, plot.forced=FALSE, plot.memory=TRUE,use.median=TRUE,memory.true=TRUE)
  # sjekk om AR1 snudd også blir positiv
  
  n = 200
  a = 0.55
  b = 0.4
  Hs = seq(from = a,to = a+b,length.out=n)
  Hs2 = rev(Hs)
  y = fgn_timedep_sim(n,sigma=1,a=a,b=b)
  yy = rev(y)
  yyy = fgn_timedep_sim(n,sigma=1,Hs = Hs2)
  
  r = inla.ews(y,model="fgn",memory.true=Hs)
  rr = inla.ews(yy,model="fgn",memory.true=Hs2)
  rrr = inla.ews(yyy,model="fgn",memory.true=Hs2)
  summary(r)
  summary(rr)
  summary(rrr)
  plot(r,plotoptions)
  plot(rr,plotoptions)
  plot(rrr,plotoptions)
  
  t1 = r$inlafit$summary.hyperpar$mode
  t2 = rr$inlafit$summary.hyperpar$mode
  t3 = rrr$inlafit$summary.hyperpar$mode
  
  #sigma er feil
  
}

if(FALSE){
  time = seq(from=0,to=1,length.out=n)
  
  theta = t1
  envir = environment(sys.call()[[1]])
  q1 = Q()
  S1 = as.matrix(covmatter())
  
  theta = t2
  envir = environment(sys.call()[[1]])
  q2 = Q()
  S2 = as.matrix(covmatter())
  
  tt = c(3.2, -2,0.8)
  theta = tt
  envir = environment(sys.call()[[1]])
  qt = Q()
  St = as.matrix(covmatter())
  
  plot(diag(S1))
  plot(diag(St))
  plot(S1[1,])
  plot(St[1,])
}


if(FALSE){ 
  #
  if(TRUE){
    #tau = exp(15)
    
    
    hyperg = function(a,b,c,z,trunc){
      
      kk = 0:(trunc-1)
      
      cp = cumprod((a+kk)*(b+kk)/(c+kk)*z/(1:trunc)  ) 
      
      s = 1 + sum(cp)
      
      return(s)
      
    }
    
    
    incBeta = function(z,a,b,trim){
      
      hg = 1/a*z^a*hyperg(a,1-b,1+a,z,trunc=trim) #take this with a grain of salt, found by trial-and-error
      
      return(hg)
    }
    
    gammas = function(a,b,n,t,s){
      if((b *(s + t - 2 *max(s, t)))/n==0){
        0
      }
      else{
        ( gamma(-1 + 2 *a + (b *(s + t))/n)*gamma( 1 - 2 *a - (2*b *max(s, t))/n))/gamma((b *(s + t - 2 *max(s, t)))/n)
      }
    }
    
    
    
    R2test = function(a,b,n,t,s){
      
      (1/(2*a*n+b*(s+t)))*n*( sqrt( (a + (b* max(s,t)/n))*(-1 + 2 *a + (2* b* max(s,t)/n))/beta(
        2 - 2 *(a + (b *max(s,t)/n)), -(1/2) + a + (b *max(s,t)/n)))*sqrt( (a + (b* min(s,t)/n))*(-1 + 2 *a + (2* b* min(s,t)/n))/beta(
          2 - 2 *(a + (b *min(s,t)/n)), -(1/2) + a + (b *min(s,t)/n))))*(beta(-(1/2) + a + (b *max(s, t))/n, -((
            2* (-1 + a)* n + b* max(s, t) + b *min(s, t))/n))* beta((n + b* max(s, t) - b* min(s, t))/ n, ((-1 + 2* a)* n + b* max(s, t) + b* min(s, t))/n)* min(s, t)^( 2* a + (b* (s + t))/n) + 
              beta(-(1/2) + a + (b *min(s, t))/n, -((2*(-1 + a)*n +b*max(s, t) + b*min(s, t))/n))*
              ((-incBeta(min(s, t)/max(s, t),  1-2*a-(2*b*max(s, t))/n, -1 + 2 *a + (b *(s + t))/n,100) 
                + gammas(a,b,n,t,s)
                
              )*min(s, t)^(2 *a + (b* (s + t))/n) + (n *
                                                       hyperg(2 - 2* a - (b *(s + t))/n, (n - b* max(s, t) + b *min(s, t))/n, (2 *n - b* max(s, t) + b* min(s, t))/n, min(s, t)/max(s, t),100) 
                                                     *max(s, t)^(2*a+(b*(s + t))/n)*(min(s, t)/max(s, t))^(1+(b*(-max(s, t)+min(s, t)))/n))/(n-b*max(s, t)+b*min(s, t))))
    }
    
    Rfgn = function(a,b,n,t,s){
      R2test(a,b,n,t+1,s+1) - R2test(a,b,n,t+1,s) - R2test(a,b,n,t,s+1) + R2test(a,b,n,t,s)
    }
    
    interpret.theta = function() {
      if(!is.null(envir)){
        timee=get("time",envir)
        nn=get("n",envir)
      }
      kappa = exp(theta[1])
      bmax = 0.5/(timee[nn]-timee[1])
      bmin = -bmax
      b = bmin+(bmax-bmin)/(1+exp(-theta[2]))
      
      low = min(b*timee[1],b*timee[nn])
      high = max(b*timee[1],b*timee[nn])
      amin = 0.5-low
      amax = 1-high
      a = amin + (amax-amin)/(1+exp(-theta[3]))
      
      Hs = a+b*timee
      return(list(Hs = Hs, kappa = kappa, a=a,b=b,amin=amin,amax=amax,bmin=bmin,bmax=bmax))
    }
    
    covmatter = function(){
      if(!is.null(envir)){
        nn=get("n",envir)
      }
      
      hyperparam = interpret.theta()
      Hs = hyperparam$Hs
      kappa = hyperparam$kappa
      sx = 1/sqrt(hyperparam$kappa)
      a = hyperparam$a
      b=hyperparam$b
      
      covmat = data.frame(row.names = NULL,col.names=NULL)
      
      for (i in 1:nn) {
        for (j in 1:nn) {
          #covmat[i,j] =round( Rfgn(a,b,nn,i,j),3)
          covmat[i,j] = Rfgn(a,b,nn,i,j)
        }
      }   
      return(covmat)
    }
    
    Q = function()  {
      if(!is.null(envir)){
        nn=get("n",envir)
      }
      
      hyperparam = interpret.theta()
      Hs = hyperparam$Hs
      kappa = hyperparam$kappa
      sx = 1/sqrt(hyperparam$kappa)
      a = hyperparam$a
      b=hyperparam$b
      
      covmat = data.frame(row.names = NULL,col.names=NULL)
      
      for (i in 1:nn) {
        for (j in 1:nn) {
          #covmat[i,j] =round( Rfgn(a,b,nn,i,j),3)
          covmat[i,j] = Rfgn(a,b,nn,i,j)
        }
      }    
      
      # url = paste0("https://www.wolframcloud.com/obj/e965ae9a-187b-469b-9d5c-c9f5eaa14656?n1=",
      #              nn,"&a1=",a,"&b1=",b,"&sigma1=",sx)
      # 
      # df = scan(url)
      # 
      # covmat = matrix(df,ncol=n)
      
      return (solve(covmat))
    }
  }
  }
