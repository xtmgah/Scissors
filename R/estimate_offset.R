#'
#' @export
estimate_offset = function(data.centered,rawmat,exonset,
                           smoothness=0.7,draw.plot=FALSE,main="Gene") {

  medscale = data.centered$msf; exonset=process.obj$exonset;
  n = ncol(data.centered$outdata);

  case.sse = apply(data$outdata,2,FUN=bisqaure.sse)
  exonbase = c()
  for (j in 1:nrow(exonset)) {
    exonbase = c(exonbase,c(exonset[j,2]:exonset[j,3]))
  }
  goodcase = which(apply(rawmat[exonbase,],2,mean)>10)

  ## estimate g by lowess
  g = rep(1,n)
  if (length(goodcase)>floor(0.1*n)) {
    lowess1 = lowess(case.sse[goodcase] ~ medscale[goodcase],f=smoothness)

    tau = lowess1[[2]][which.min((lowess1[[1]]-1)^2)[1]]
    if (tau<0) {
      tau=min(lowess1[[2]][which(lowess1[[2]]>0)])
    }
    ## Rescale case.sse
    case.sse.adj = sqrt(case.sse/tau)

    lowess2 = as.list(NULL);
    lowess2$x = lowess1[[1]][which(lowess1[[2]]>0)]
    lowess2$y = sqrt(lowess1[[2]][which(lowess1[[2]]>0)]/tau)

    g[goodcase[order(medscale[goodcase])[which(lowess1[[2]]>0)]]] = lowess2$y
    g[which(g<1)] = 1
    g[which(medscale<1)] = 1
  }

  outdata = sweep(x=data.centered$outdata,2,g,FUN="/")

  if (draw.plot) {
    colmat=rep("black",length(goodcase))
    colmat[order(medscale[goodcase],decreasing=TRUE)] = rainbow(n=length(goodcase),start=0,end=0.8)

    plot(medscale[goodcase],case.sse.adj[goodcase],col=colmat,main=main)
    lines(lowess2,lwd=2)
    abline(v=1,h=1,lty=2)
    points(medscale[order(medscale)],g[order(medscale)],type="l",lwd=3,col="red")

    legend("topright",bty="n",
           legend=c(paste("min =",round(min(lowess2[[2]]),digits=3)),
                    paste("max =",round(max(lowess2[[2]]),digits=3))))
    legend("topleft",bty="n",
           legend=c(paste("# cases involved =",length(goodcase)),
                    paste("gene length =",length(exonbase)),
                    paste("smoothness  =",round(smoothness,digits=1))))
  }
  return(list(g=g,outdata=outdata))
}
