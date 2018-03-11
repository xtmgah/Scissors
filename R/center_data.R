#' Centering data
#' @export
center_data = function(data,type=1, ...) {
  msfres = scale.factor(X=data, ...)
  data.center = msfres$mean.vec
  msf=msfres$msf

  if (type==0) {
    outdata = data - data.center
  } else if (type==1) {
    outdata = msfres$NewX
  }
  return(list(outdata=outdata,msf=msf,data.center=data.center))
}
