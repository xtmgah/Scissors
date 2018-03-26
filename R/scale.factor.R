#'
#' @export
scale.factor = function(X,average="median",adjval=NULL) {
  ## % Obtain scale factors based on a adjusted center
  mean.vec = apply(X,1,FUN=function(x){adj.center(x,average=average,adjval=adjval)})
  if (sum(mean.vec^2)>0) {
    msf = as.vector(t(mean.vec)%*%X/(sum(mean.vec^2)));
  } else {
    msf = rep(0,ncol(X))
  }
  NewX = X - mean.vec%*%t(msf)
  return(list(msf=msf,mean.vec=mean.vec,NewX=NewX))
}
