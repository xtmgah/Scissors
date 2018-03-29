#' Log-transforming data
#'
#' This function is used to get
#' @param data data to be transformed
#' @param logshift.val
#' @keywords log shift parameter
#' @export
#' @examples

logtransform_data = function(data,logshift.val=NULL, param.grid=NULL, ...) {
  # Find exon bases out of intron-contained coverage
  #     & Take log-transformation.
  if (is.null(logshift.val)) {
    logshift.val = getShiftParam(X=data, plot=FALSE, param.grid=param.grid,...)$optim.param
  }
  outdata = log10(data + logshift.val) - log10(logshift.val) ;
  return(list(outdata=outdata,logshift.val=logshift.val))
}
