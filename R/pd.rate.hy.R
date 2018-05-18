#' Projection depth
#' @export
pd.rate.hy = function(x) {
  # projection depth
  m = mad(x)
  if (m<1e-5) {
    return(rep(0,length(x)))
  } else {
    return((x-median(x))/m)
  }
}
