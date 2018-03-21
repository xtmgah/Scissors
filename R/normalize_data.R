#'
#' @export
normalize_data = function(datalog,rawmat,exonset,
                          center.type=1,smoothness=0.7,draw.plot=FALSE,main="Gene", ...) {

  data.centered = center_data(data=datalog,type=center.type, ...)
  g = estimate_offset(data.centered=data.centered,rawmat=rawmat,exonset=exonset,
                                    smoothness=smoothness,draw.plot=draw.plot,main=main)

  outdata = sweep(x=data.centered$outdata,2,g,FUN="/")
  return(list(outdata=outdata,msf=data.centered$msf,g=g,
              data.center=data.centered$data.center))
}
