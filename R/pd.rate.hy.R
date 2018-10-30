#' Projection depth
#' @export
pd.rate.hy = function(x,qrsc=FALSE) {
  # projection depth
  m = mad(x)
  if (m<1e-5) {
    return(rep(0,length(x)))
  } else {
    if (qrsc) {
      rsc=compScales(x)
      return((x-median(x))/rsc$sa)
    } else {
      return((x-median(x))/m)
    }
  }
}
