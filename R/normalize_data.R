#'
#' @export
normalize_data = function(datalog,rawmat,exonset,
                          center.type=1,smoothness=0.7,draw.plot=FALSE,main="Gene", ...) {

  data.centered = center_data(data=datalog,type=center.type, ...)
  g1.offset = estimate_offset(data.centered=data.centered,rawmat=rawmat,exonset=exonset,
                                    smoothness=smoothness,draw.plot=F,main=main)
  msf = data.centered$msf
  datactr = sweep(x=data.centered$outdata,2,g1.offset$g,FUN="/")
  
  ## loop until no different variations 
  g2 = g1.offset$g
  while (max(g2)>1.05) {
    g2.offset = estimate_offset(msf=msf,cenmat=datactr,rawmat=rawmat,exonset=exonset,
                                smoothness=smoothness,draw.plot=F)
    g2 = g2.offset$g
    datactr = sweep(x=datactr,2,g2,FUN="/")
  }
  if (draw.plot) {
    plot_offset(offset.obj=g1.offset,draw.legend=T,main=GeneName)
    plot_offset(offset.obj=g2.offset,draw.legend=T,main=GeneName)
  }
  return(list(outdata=datactr,msf=msf,g1.offset=g1.offset,g2.offset=g2.offset,
              data.center=data.centered$data.center))
}
