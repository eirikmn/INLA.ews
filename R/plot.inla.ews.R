#' Plot inla.ews model
#'
#' Plots results from inla.ews S3 class.
#'
#' @param x \code{inla.ews} S3 class. Output of \code{\link{inla.ews}} function
#' @param plot.options list with settings for plot.
#' @param postscript Boolean variable indicating if postscript files should be produced instead.
#' @param pdf Boolean variable indicating if pdf files should be produced instead.
#' @param prefix The prefix for created files. Additional numbering is added.
#' @param ... Additional arguments to \code{postscripts()}, \code{pdf()} or \code{dev.new()}.
#'
#' @examples 
#' \donttest{
#' n = 200
#' sigma = 1
#' a=0.2
#' b=0.7/n
#' time = 1:n
#' phis = a+b*time
#' data=ar1_timedep_sim(n,phis=phis)
#' 
#' object = inla.ews(data,model="ar1", print.progress=TRUE,
#'                     memory.true=phis)
#' plot(object)
#' 
#' }
#'
#' @author Eirik Myrvoll-Nilsen, \email{eirikmn91@gmail.com}
#' @seealso \code{\link{inla.ews}}
#' @keywords inla.ews plot
#' @importFrom ggplot2 ggplot xlab ylab geom_line theme_bw geom_ribbon aes ylim 
#' geom_point labs ggtitle geom_segment
#' @importFrom graphics abline par
#' @importFrom rlang .data
#' @importFrom grDevices dev.cur dev.new dev.off
#' @export
plot.inla.ews <- function(x,
                          plot.options=list(plot.hyper=TRUE,
                                            plot.forced=TRUE,
                                            plot.memory=TRUE,
                                            use.median=TRUE,memory.true=TRUE
                                            ),
                          postscript=FALSE,
                          pdf=FALSE,
                          prefix = "INLA.ews.plots/figure-",
                          ...){
  if(!postscript && !pdf){
    dev=getOption("device")
    if(!is.character(dev) || dev != "RStudioGD"){
      dev.new(...)
    }
  }else{
    dir = dirname(prefix)
    if (!file.exists(dir) && nchar(dir) > 0L) {
      dir.create(dir, recursive=TRUE)
    } else {
      stopifnot(file.info(dir)$isdir)
    }
  }
  
  figure.count = 1L
  oldpar = par()
  if(!is.null(x$results) ){
    
    if(plot.options$plot.hyper){
      margnames = names(x$results$marginals)
      for(i in 1:length(x$results$marginals)){
        figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
        LOW = x$results$summary[[margnames[i]]]$quant0.025
        UPP = x$results$summary[[margnames[i]]]$quant0.975
        
        margi = x$results$marginals[[margnames[i]]]
        margix = margi[,1]
        margiy = margi[,2]
        intx = margix[ margix <= UPP & margix >= LOW]
        inty = margiy[ margix <= UPP & margix >= LOW]
        
        # mm = margic
        # bb=INLA::inla.zmarginal(mm,silent=TRUE)
        
        #intx = x$results$marginals[[i]]$x[ x$results$marginals[[i]]$x < UPP & x$results$marginals[[i]]$x > LOW]
        #inty = x$results$marginals[[i]]$y[ x$results$marginals[[i]]$x < UPP & x$results$marginals[[i]]$x > LOW]
        ggp = ggplot()  + theme_bw() +
          xlab(names(x$results$marginal)[i]) + ylab("Density") + 
          ggtitle(paste0("Posterior distribution: ",names(x$results$marginals)[i])) +
          geom_ribbon(data=data.frame(intx=intx,inty=inty,lower=numeric(length(intx))),
                      mapping=aes(x=.data$intx,ymin=.data$lower,ymax=.data$inty), color="red",fill="red",alpha=0.3,linewidth=0) +
          geom_segment(aes(x=LOW,xend=LOW,y=0,yend=inty[1]),col="red",linewidth=1.) +
          geom_segment(aes(x=UPP,xend=UPP,y=0,yend=inty[length(inty)]),col="red",linewidth=1.) +
          geom_line(data=x$results$marginals[[margnames[i]]],mapping=aes(x=.data$x,y=.data$y), linewidth=1.2)
        
        print(ggp)
        if(postscript || pdf){
          if (names(dev.cur()) != "null device") {
            dev.off()
          }
        }
      }
    }
    
  }

  
  if(plot.options$plot.forced && max(abs(x$results$summary$alltrend$mean)) > 0.0001){
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    
    ggd = data.frame(y=x$.args$inladata$y, time=x$.args$timesteps, mean=x$results$summary$alltrend$mean, 
                     lower=x$results$summary$alltrend$quant0.025, upper=x$results$summary$alltrend$quant0.975)
    
    ggp = ggplot(ggd, aes(x=.data$time)) + geom_line(aes(y=.data$y),col="gray") + theme_bw() +
      geom_ribbon(aes(ymin=.data$lower, ymax=.data$upper),col="red",fill="red",alpha=0.3,linewidth=1.2) +
      geom_line(aes(y=.data$mean),col="blue",linewidth=1.2)
    #
    ## EMN: KOMT HIT
    #
    # 
    # par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
    # ggd = data.frame(y=x$.args$data,time=x$.args$inladata$time,
    #                  mean=x$forced$mean)
    # if(!is.null(x$forced$quant0.025)){
    #   ggd$lower = x$forced$quant0.025
    #   ggd$upper = x$forced$quant0.975
    # }
    # gg = ggplot(ggd,aes(.data$time)) + theme_bw()+ 
    #   xlab("Time")+ylab("Observations")+
    #   labs(title="Forced response")+
    #   geom_line(aes(y=.data$y),color="gray",size=0.7,alpha=0.9)
    # if(!is.null(x$forced$quant0.025)){
    #   #gg = gg + geom_line(aes(y=.data$lower),col="")
    #   gg = gg + geom_ribbon(aes(ymin=.data$lower,ymax=.data$upper),color="red",
    #                         fill="red",alpha=0.3,size=0.9)
    # }
    # gg = gg + geom_line(aes(y=.data$mean),color="blue",size=1.1)
    # print(gg)
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
  }
  
  if(plot.options$plot.memory){
    if(tolower(x$.args$model) %in% c("ar1","ar(1)","1","ar1g")){
      plot.df = data.frame(time=x$.args$timesteps,
                           mean=x$results$summary$phi$mean,
                           median=x$results$summary$phi$median,
                           lower=x$results$summary$phi$hpd0.95lower,
                           upper=x$results$summary$phi$hpd0.95upper)
      ylim=c(0,1)
      ylab=expression(paste("lag-one correlation: ",phi))
    }else if(tolower(x$.args$model) %in% c("fgn","lrd")){
      plot.df = data.frame(time=x$.args$inladata$time,
                           mean=x$results$summary$H$mean,
                           median=x$results$summary$H$median,
                           lower=x$results$summary$H$hpd0.95lower,
                           upper=x$results$summary$H$hpd0.95upper)
      ylim=c(0.5,1)
      ylab="Hurst exponent: H"
    }else if(tolower(x$.args$model) %in% c("ar1","ar(1)","1","ar1g")){
      plot.df = data.frame(time=x$.args$timesteps,
                           mean=x$results$summary$phi$mean,
                           median=x$results$summary$phi$median,
                           lower=x$results$summary$phi$hpd0.95lower,
                           upper=x$results$summary$phi$hpd0.95upper)
      ylim=c(0,1)
      ylab=expression(paste("lag-one correlation: ",phi))
    }
      
    figure.count <- new.plot(postscript,pdf,prefix,figure.count,...) +1
    par(mfrow=c(1,1),mar=c(5,4,4,2)+0.1)
    gg2 <- ggplot(data=plot.df,aes(x=.data$time)) + 
      geom_ribbon(aes(ymin=.data$lower,ymax=.data$upper),color="red",fill="red",alpha=0.3)+
      theme_bw()+ylab(ylab)+xlab("Time")+ylim(ylim)+
      labs(title="Evolution of memory parameter")
    if(plot.options$use.median){
      gg2 <- gg2 + geom_line(aes(y=.data$median),color="blue")
    }else{
      gg2 <- gg2 + geom_line(aes(y=.data$mean),color="blue")
    }
    if(plot.options$memory.true && !is.null(x$.args$memory.true)){
      gg2 <- gg2 + 
        geom_line(data=data.frame(time=x$.args$inladata$time,
                                  true=x$.args$memory.true),
                  aes(y=.data$true),col="black",alpha=0.3)
    }
    print(gg2)
    if(postscript || pdf){
      if (names(dev.cur()) != "null device") {
        dev.off()
      }
    }
    
  }
  return(invisible(x))
}



new.plot = function(postscript,pdf,prefix,figure.count,...)
{
  
  #dev = getOption("device")
  if(postscript && pdf){
    stop("Multiple file types have been seleced.")
  }
  if(postscript) {
    ending=".eps"
  }else if(pdf){
    ending=".pdf"
  }
  if(!postscript && !pdf){
    dev=getOption("device")
    if(!is.character(dev) || dev != "RStudioGD"){
      dev.new(...)
    }
  }else{
    
    file.found=FALSE
    while(!file.found){
      filename=paste(prefix,figure.count,ending,sep="")
      
      if(file.exists(filename)){
        figure.count <- figure.count +1L
      }else{
        file.found=TRUE
      }
    }
    if(postscript){
      postscript(file=filename,...)
    }else if(pdf){
      pdf(file=filename,...)
    }
  }
  return (invisible(figure.count))
}
