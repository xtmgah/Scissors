#'
#' @export
normalize_data = function(datalog,rawmat,exonset,loop=TRUE,
                          center.type=1,smoothness=0.7,draw.plot=FALSE,main="Gene", ...) {

  data.centered = center_data(data=datalog,type=center.type, ...)
  g1.offset = estimate_offset(data.centered=data.centered,rawmat=rawmat,exonset=exonset,
                                    smoothness=smoothness,draw.plot=F,main=main)
  msf = data.centered$msf
  datactr = sweep(x=data.centered$outdata,2,g1.offset$g,FUN="/")
  goodcase = g1.offset$goodcase

  ## loop until no different variations
  if (sum((g1.offset$g-1)^2)<1e-10) {
    cat("No adjustment applied. Data do not have enough expression.","\n")
    cat("........Plots are omitted.","\n")
    g2.offset=g1.offset
  } else {
    if (loop) {
      g2.offset = g1.offset
      k = 1
      while ((k<5) & (max(g2.offset$g)>1.05)) {
        g2.offset = estimate_offset(msf=msf,cenmat=datactr,rawmat=rawmat,exonset=exonset,
                                    smoothness=smoothness,draw.plot=F)
        datactr = sweep(x=datactr,2,g2.offset$g,FUN="/")
        k = k+1
      }
      if (draw.plot) {
        plot_offset(offset.obj=g1.offset,draw.legend=T,
                    main=paste(GeneName,": before normalization"))
        plot_offset(offset.obj=g2.offset,draw.legend=T,
                    main=paste(GeneName,": after normalization"))
      }
    } else {
      g2.offset = g1.offset
      if (draw.plot) {
        plot_offset(offset.obj=g1.offset,draw.legend=T,main=GeneName)
      }
    }

  }
  return(list(outdata=datactr,msf=msf,g1.offset=g1.offset,g2.offset=g2.offset,
              data.center=data.centered$data.center,goodcase=goodcase))
}
